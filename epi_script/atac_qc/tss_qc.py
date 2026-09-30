"""ATAC cut-site TSS enrichment with ENCODE-style edge normalization.

Set ``ATAC_TSS_ARTIFACTS`` to ``profile`` (default), ``matrix``, or
``heatmap`` to select how many optional diagnostic artifacts are retained.
All modes calculate the same enrichment score and write its profile.

The offset histogram is accumulated per contig with ``numpy`` instead of
handing a ``403,699 x 400`` interval BED to ``bedtools coverage``.  Both
formulations count exactly the same (TSS, cut-site) pairs: for every TSS the
same ``range(-window, window, bin_size)`` bins are summed with the same
``genomic_end <= 0`` skip and ``genomic_start = max(genomic_start, 0)``
clamping.  The numpy path avoids writing the ~5 GB bins file and the
multi-hour ``bedtools coverage`` run on large pooled tagAligns.
"""

from __future__ import annotations

import gzip
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterator, Optional

import numpy as np


TSS_METHOD = "tn5_cutsite_strand_aware_10bp_edge100bp"
TSS_ARTIFACT_LEVELS = ("profile", "matrix", "heatmap")
MATRIX_TSS_LIMIT = 5000
EDGE_BP = 100

# awk 输出定宽记录 "contig 序号 + 坐标 + 换行"，便于 numpy 按块批量解析。
_CUTSITE_ID_WIDTH = 4
_CUTSITE_POS_WIDTH = 10
_CUTSITE_RECORD_WIDTH = _CUTSITE_ID_WIDTH + _CUTSITE_POS_WIDTH + 1
_CUTSITE_MAX_POS = 10 ** _CUTSITE_POS_WIDTH - 1
_CUTSITE_READ_CHUNK = 1 << 26
_MAX_CONTIGS = 10 ** _CUTSITE_ID_WIDTH - 1

_CUTSITE_AWK_TEMPLATE = """BEGIN {{ FS = "{sep}" }}
NR == FNR {{ id[$1] = $2; next }}
{{
    idx = id[$1]
    if (idx == "") next
    if ($6 == "+") pos = $2
    else if ($6 == "-") pos = $3 - 1
    else next
    if (pos < 0 || pos > {max_pos}) next
    printf "{fmt}\\n", idx, pos
}}
"""


def resolve_tss_artifact_level(artifacts: Optional[str] = None) -> str:
    """Return the requested diagnostic-artifact level.

    The explicit argument is useful to callers with their own CLI.  When it is
    omitted, use the process-level setting so existing callers remain backward
    compatible without needing another positional argument.
    """
    requested = artifacts if artifacts is not None else os.environ.get(
        "ATAC_TSS_ARTIFACTS", "profile"
    )
    level = requested.strip().lower()
    if level not in TSS_ARTIFACT_LEVELS:
        raise ValueError(
            "ATAC TSS artifact level must be one of "
            f"{', '.join(TSS_ARTIFACT_LEVELS)}; got {requested!r}"
        )
    return level


def score_profile(
    offset_counts: dict[int, int],
    tss_count: int,
    window: int,
    bin_size: int,
    edge_bp: int = EDGE_BP,
) -> tuple[Optional[float], dict[int, float]]:
    raw_profile = {
        offset: count / tss_count
        for offset, count in sorted(offset_counts.items())
    }
    edge_values = [
        value for offset, value in raw_profile.items()
        if offset < -window + edge_bp or offset >= window - edge_bp
    ]
    edge_mean = sum(edge_values) / len(edge_values) if edge_values else 0.0
    if edge_mean <= 0:
        return None, {}
    normalized = {offset: value / edge_mean for offset, value in raw_profile.items()}
    return max(normalized.values(), default=None), normalized


def _load_tss(
    tss_bed: Path,
) -> tuple[int, dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]]:
    """读入 TSS BED，返回 TSS 总数和按 contig/链分组的中心坐标。

    每个分组的第二个返回值是该 TSS 在 BED 中的 0 基行号，用于回溯 matrix
    模式选取的行；链不是 ``+``/``-`` 的一律按旧实现当作 ``+``。
    """
    contig_order: list[str] = []
    centers: dict[str, list[int]] = {}
    row_ids: dict[str, list[int]] = {}
    strands: dict[str, list[str]] = {}
    tss_count = 0
    with tss_bed.open("r", encoding="utf-8", errors="replace") as source:
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            chrom = fields[0]
            strand = fields[5] if len(fields) >= 6 and fields[5] in {"+", "-"} else "+"
            if chrom not in centers:
                contig_order.append(chrom)
                centers[chrom] = []
                row_ids[chrom] = []
                strands[chrom] = []
            centers[chrom].append((start + end) // 2)
            row_ids[chrom].append(tss_count)
            strands[chrom].append(strand)
            tss_count += 1
    grouped: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {}
    for chrom in contig_order:
        by_strand: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        chrom_strands = np.array(strands[chrom], dtype=object)
        chrom_centers = np.array(centers[chrom], dtype=np.int64)
        chrom_rows = np.array(row_ids[chrom], dtype=np.int64)
        for strand in ("+", "-"):
            mask = chrom_strands == strand
            by_strand[strand] = (chrom_centers[mask], chrom_rows[mask])
        grouped[chrom] = by_strand
    return tss_count, grouped


def _parse_cut_site_block(block: bytes) -> tuple[np.ndarray, np.ndarray]:
    rows = np.frombuffer(block, dtype=np.uint8).reshape(-1, _CUTSITE_RECORD_WIDTH)
    contig = (
        rows[:, :_CUTSITE_ID_WIDTH]
        .copy()
        .view(f"S{_CUTSITE_ID_WIDTH}")
        .ravel()
        .astype(np.int32)
    )
    position = (
        rows[:, _CUTSITE_ID_WIDTH:_CUTSITE_ID_WIDTH + _CUTSITE_POS_WIDTH]
        .copy()
        .view(f"S{_CUTSITE_POS_WIDTH}")
        .ravel()
        .astype(np.int64)
    )
    return contig, position


def _iter_cut_site_blocks(
    contig_ids: dict[str, int], tagalign_file: Path, workdir: Path
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """流式产出 (contig 序号, cut-site 坐标) 的定宽记录块。

    cut site 定义与旧实现一致：``$6 == "+"`` 取 ``$2``，``$6 == "-"`` 取
    ``$3 - 1``，其余 strand 值直接丢弃（旧 awk 也只输出这两种）。
    """
    map_file = workdir / "cut_site_contigs.tsv"
    with map_file.open("w", encoding="utf-8", newline="") as handle:
        for name, idx in contig_ids.items():
            handle.write(f"{name}\t{idx}\n")
    program_file = workdir / "cut_sites.awk"
    program_file.write_text(
        _CUTSITE_AWK_TEMPLATE.format(
            sep="\\t",
            max_pos=_CUTSITE_MAX_POS,
            fmt=f"%0{_CUTSITE_ID_WIDTH}d%0{_CUTSITE_POS_WIDTH}d",
        ),
        encoding="utf-8",
    )
    command = (
        "set -o pipefail; "
        f"zcat -f {shlex.quote(str(tagalign_file))} | "
        f"awk -f {shlex.quote(str(program_file))} {shlex.quote(str(map_file))} -"
    )
    process = subprocess.Popen(
        ["bash", "-c", command],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout = process.stdout
    if stdout is None:
        raise RuntimeError("cut-site 提取失败：无法读取 awk 输出")
    buffer = b""
    try:
        while True:
            chunk = stdout.read(_CUTSITE_READ_CHUNK)
            if not chunk:
                break
            buffer += chunk
            usable = len(buffer) - (len(buffer) % _CUTSITE_RECORD_WIDTH)
            if usable == 0:
                continue
            block, buffer = buffer[:usable], buffer[usable:]
            yield _parse_cut_site_block(block)
    finally:
        stdout.close()
    stderr = process.stderr.read().decode("utf-8", errors="replace") if process.stderr else ""
    if process.stderr is not None:
        process.stderr.close()
    returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(
            f"cut-site 提取失败 (exit {returncode}): {stderr.strip() or 'awk/zcat failed'}"
        )


def _collect_cut_site_counts(
    contig_ids: dict[str, int], tagalign_file: Path, workdir: Path
) -> dict[str, np.ndarray]:
    """按 contig 汇总每个碱基上的 cut-site 数（等价于旧实现的 ``-b`` 输入）。"""
    per_contig: dict[int, list[np.ndarray]] = {}
    for contig, position in _iter_cut_site_blocks(contig_ids, tagalign_file, workdir):
        order = np.argsort(contig, kind="stable")
        sorted_contig = contig[order]
        boundaries = np.flatnonzero(np.diff(sorted_contig)) + 1
        keys = sorted_contig[np.concatenate(([0], boundaries))]
        for key, group in zip(keys, np.split(order, boundaries)):
            per_contig.setdefault(int(key), []).append(position[group])
    counts: dict[str, np.ndarray] = {}
    for name, idx in contig_ids.items():
        parts = per_contig.pop(idx, None)
        if not parts:
            continue
        position = np.concatenate(parts) if len(parts) > 1 else parts[0]
        counts[name] = np.bincount(position, minlength=int(position.max()) + 1)
    return counts


def _tss_inputs(
    tagalign_file: Path, tss_bed: Path
) -> tuple[int, dict[str, dict[str, tuple[np.ndarray, np.ndarray]]], dict[str, np.ndarray]]:
    """一次读入 TSS 与 tagAlign，产出 profile 与 matrix 共用的中间结果。"""
    tss_count, grouped = _load_tss(tss_bed)
    if tss_count == 0:
        return 0, grouped, {}
    if len(grouped) > _MAX_CONTIGS:
        raise ValueError(
            f"TSS BED 含 {len(grouped)} 个 contig，超过定宽编码上限 {_MAX_CONTIGS}"
        )
    contig_ids = {chrom: index + 1 for index, chrom in enumerate(grouped)}
    with tempfile.TemporaryDirectory(prefix="atac_tss_") as tmpdir:
        counts_by_contig = _collect_cut_site_counts(
            contig_ids, tagalign_file, Path(tmpdir)
        )
    return tss_count, grouped, counts_by_contig


def _cumulative_counts(counts: np.ndarray) -> np.ndarray:
    cumulative = np.zeros(counts.size + 1, dtype=np.int64)
    np.cumsum(counts, dtype=np.int64, out=cumulative[1:])
    return cumulative


def _bin_counts(
    centers: np.ndarray,
    strand: str,
    offset: int,
    bin_size: int,
    cumulative: np.ndarray,
    upper: int,
) -> np.ndarray:
    """按旧 ``_write_tss_bins`` 的 bin 定义，统计每个 TSS 该 offset 上的 cut-site 数。"""
    if strand == "+":
        low = centers + offset
        high = low + bin_size
    else:
        high = centers - offset
        low = high - bin_size
    np.clip(low, 0, upper, out=low)
    np.clip(high, 0, upper, out=high)
    return cumulative[high] - cumulative[low]


def _accumulate_bins(
    counts_by_contig: dict[str, np.ndarray],
    grouped: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
    window: int,
    bin_size: int,
    tss_count: int,
    want_totals: bool,
) -> tuple[dict[int, int], Optional[np.ndarray]]:
    offsets = np.arange(-window, window, bin_size, dtype=np.int64)
    offset_counts: dict[int, int] = {int(offset): 0 for offset in offsets}
    totals = np.zeros(tss_count, dtype=np.int64) if want_totals else None
    for chrom, by_strand in grouped.items():
        counts = counts_by_contig.get(chrom)
        if counts is None or counts.size == 0:
            continue
        cumulative = _cumulative_counts(counts)
        upper = counts.size
        for strand in ("+", "-"):
            centers, row_ids = by_strand[strand]
            if centers.size == 0:
                continue
            for offset in offsets:
                values = _bin_counts(
                    centers, strand, int(offset), bin_size, cumulative, upper
                )
                offset_counts[int(offset)] += int(values.sum())
                if totals is not None:
                    totals[row_ids] += values
    return offset_counts, totals


def _matrix_rows(totals: Optional[np.ndarray], limit: int = MATRIX_TSS_LIMIT) -> np.ndarray:
    """按 ``(该 TSS 的 bin 计数之和, TSS 行号)`` 从大到小选出前 ``limit`` 行。"""
    if totals is None or totals.size == 0:
        return np.zeros(0, dtype=np.int64)
    order = np.lexsort((np.arange(totals.size, dtype=np.int64), totals))
    return order[-limit:][::-1]


def _matrix_values(
    counts_by_contig: dict[str, np.ndarray],
    grouped: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]],
    selected: np.ndarray,
    tss_count: int,
    window: int,
    bin_size: int,
    edge_mean: float,
) -> tuple[list[int], np.ndarray]:
    offsets = np.arange(-window, window, bin_size, dtype=np.int64)
    positions = np.full(tss_count, -1, dtype=np.int64)
    positions[selected] = np.arange(selected.size, dtype=np.int64)
    matrix = np.zeros((selected.size, offsets.size), dtype=np.float64)
    for chrom, by_strand in grouped.items():
        counts = counts_by_contig.get(chrom)
        if counts is None or counts.size == 0:
            continue
        cumulative = _cumulative_counts(counts)
        upper = counts.size
        for strand in ("+", "-"):
            centers, row_ids = by_strand[strand]
            if centers.size == 0:
                continue
            columns = positions[row_ids]
            keep = columns >= 0
            if not keep.any():
                continue
            centers = centers[keep]
            columns = columns[keep]
            for column, offset in enumerate(offsets):
                values = _bin_counts(
                    centers, strand, int(offset), bin_size, cumulative, upper
                )
                matrix[columns, column] = values / edge_mean
    return [int(offset) for offset in offsets], matrix


def compute_offset_counts(
    tagalign_file: Path,
    tss_bed: Path,
    window: int = 2000,
    bin_size: int = 10,
) -> tuple[int, dict[int, int]]:
    """按 TSS offset 统计 cut-site 总数，返回 ``(tss_count, offset_counts)``。"""
    tss_count, grouped, counts_by_contig = _tss_inputs(tagalign_file, tss_bed)
    if tss_count == 0:
        return 0, {}
    offset_counts, _ = _accumulate_bins(
        counts_by_contig, grouped, window, bin_size, tss_count, want_totals=False
    )
    return tss_count, offset_counts


def compute_tss_enrichment(
    tagalign_file: Path,
    tss_bed: Optional[Path],
    output_file: Path,
    window: int = 2000,
    bin_size: int = 10,
    artifacts: Optional[str] = None,
) -> dict[str, object]:
    artifact_level = resolve_tss_artifact_level(artifacts)
    base = {
        "tss_enrichment_score": None, "tss_profile_file": "", "tss_method": TSS_METHOD,
        "tss_matrix_file": "", "tss_heatmap_file": "", "tss_artifacts": artifact_level,
    }
    if tss_bed is None:
        return {**base, "warnings": ["tss_not_configured"]}
    if not tss_bed.exists():
        return {**base, "warnings": [f"tss_bed_missing:{tss_bed}"]}
    if not tagalign_file.exists():
        return {**base, "warnings": [f"tss_missing_tagalign:{tagalign_file}"]}
    missing_tools = [tool for tool in ("awk", "zcat") if shutil.which(tool) is None]
    if missing_tools:
        return {**base, "warnings": ["tss_tools_not_found:" + ",".join(missing_tools)]}

    output_file.parent.mkdir(parents=True, exist_ok=True)
    want_matrix = artifact_level in {"matrix", "heatmap"}
    try:
        # ATAC_align/finalize already applies +4/-5 to tagAlign intervals. Use
        # the strand-specific shifted 5' coordinate as a one-base cut site.
        tss_count, grouped, counts_by_contig = _tss_inputs(tagalign_file, tss_bed)
    except (RuntimeError, ValueError) as error:
        return {**base, "warnings": [f"tss_failed:{error}"]}
    if tss_count == 0:
        return {**base, "warnings": ["tss_bed_empty"]}
    offset_counts, totals = _accumulate_bins(
        counts_by_contig, grouped, window, bin_size, tss_count, want_totals=want_matrix
    )
    score, profile = score_profile(offset_counts, tss_count, window, bin_size)
    if score is None:
        return {**base, "warnings": ["tss_edge_background_zero"]}
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        handle.write("offset_start\toffset_end\tnormalized_cutsite_enrichment\n")
        for offset, value in sorted(profile.items()):
            handle.write(f"{offset}\t{offset + bin_size}\t{value:.8f}\n")
    matrix_file: Path | None = None
    heatmap_file: Path | None = None
    plot_warning: list[str] = []
    if want_matrix:
        edge_offsets = [
            offset for offset in sorted(profile)
            if offset < -window + EDGE_BP or offset >= window - EDGE_BP
        ]
        edge_mean = sum(
            offset_counts[offset] / tss_count for offset in edge_offsets
        ) / len(edge_offsets)
        # Keep only the strongest MATRIX_TSS_LIMIT TSS rows for a bounded
        # diagnostic matrix; profile-only QC avoids this second pass entirely.
        selected = _matrix_rows(totals)
        offsets, matrix = _matrix_values(
            counts_by_contig, grouped, selected, tss_count, window, bin_size, edge_mean
        )
        matrix_file = output_file.with_name(
            output_file.name.replace(".tss_profile.tsv", ".tss_matrix.tsv.gz")
        )
        with gzip.open(matrix_file, "wt", encoding="utf-8") as handle:
            handle.write("tss_id\t" + "\t".join(map(str, offsets)) + "\n")
            for row_numbers, values in zip(selected, matrix):
                handle.write(
                    str(int(row_numbers) + 1) + "\t" +
                    "\t".join(f"{value:.6f}" for value in values) + "\n"
                )
        if artifact_level == "heatmap":
            heatmap_file = output_file.with_name(
                output_file.name.replace(".tss_profile.tsv", ".tss_heatmap.png")
            )
            try:
                import matplotlib
                matplotlib.use("Agg")
                from matplotlib import pyplot as plt
                image_rows = matrix.tolist()
                if image_rows:
                    fig, axis = plt.subplots(figsize=(6, 8))
                    axis.imshow(
                        image_rows,
                        aspect="auto",
                        interpolation="nearest",
                        cmap="Reds",
                        vmin=0,
                        vmax=max(5, score),
                    )
                    axis.set_xlabel("Distance from TSS (10 bp bins)")
                    axis.set_ylabel("TSS ordered by cut-site signal")
                    fig.tight_layout()
                    fig.savefig(heatmap_file, dpi=150)
                    plt.close(fig)
            except ImportError:
                plot_warning.append("matplotlib_not_found_for_tss_heatmap")
    return {
        "tss_enrichment_score": score,
        "tss_profile_file": str(output_file),
        "tss_method": TSS_METHOD,
        "tss_matrix_file": str(matrix_file) if matrix_file is not None else "",
        "tss_heatmap_file": (
            str(heatmap_file)
            if heatmap_file is not None and heatmap_file.exists()
            else ""
        ),
        "tss_artifacts": artifact_level,
        "warnings": plot_warning,
    }
