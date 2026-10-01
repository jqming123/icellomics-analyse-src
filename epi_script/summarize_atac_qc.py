"""Summarize ATAC-seq QC metrics for a pipeline project.

FRiP is computed from one-record-per-fragment BED when available. Legacy
tagAlign inputs remain readable but are explicitly labelled as record-level:

  reads_in_peaks = bedtools intersect -a pooled.tn5.tagAlign -b peaks -wa -u
  total_reads = number of records in pooled.tn5.tagAlign
  frip = reads_in_peaks / total_reads

Default project layout follows the ATAC pipeline:

  <project>/0_data/sample_run_map.tsv
  <project>/1_result/1_alignment/<RUN>/fastp/<RUN>_fastp.json
  <project>/1_result/1_alignment/<RUN>/bowtie2/<RUN>.nodup.bam.idxstats
  <project>/1_result/2_tagalign/<BIOSAMPLE>.tn5.tagAlign.gz
  <project>/1_result/3_peak_calling/<BIOSAMPLE>_peaks.narrowPeak
  <project>/1_result/5_qc/atac_qc_summary.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Mapping, Optional, Sequence, cast

from atac_qc.alignment_qc import (
    read_idxstats,
    read_picard_dup_metrics,
    summarize_raw_bam_flagstat,
    read_picard_dup_details,
)
from atac_qc.blacklist_qc import compute_blacklist_qc
from atac_qc.common import (
    DEFAULT_BASE_DIR,
    DEFAULT_MITO_CONTIGS,
    QC_COLUMNS,
    add_optional,
    format_float,
    format_int,
    split_values,
)
from atac_qc.fastp_qc import read_fastp_reads
from atac_qc.fragment_qc import summarize_fragment_sizes
from atac_qc.library_complexity_qc import compute_library_complexity
from atac_qc.frip_qc import (
    FRIP_METHOD,
    LEGACY_FRIP_METHOD,
    compute_records_in_peaks,
    count_all_records,
    ensure_frip_tools,
    validate_frip_contigs,
)
from atac_qc.peak_qc import summarize_peak_file
from atac_qc.sample_map import (
    SAFE_ID,
    VALID_LAYOUTS,
    ReplicateRecord,
    detect_layout,
    group_replicates,
    infer_samples_from_tagalign,
    read_replicate_records,
    summarize_layout,
)
from atac_qc.status_rules import QC_MODES, evaluate_qc_statuses
from atac_qc.tss_qc import compute_tss_enrichment
from validate_atac_blacklist_resource import sha256_file


SCRIPT_DIR = Path(__file__).resolve().parent
CONFIG_PATH = SCRIPT_DIR / "epi_qc_config.sh"
CONFIG_KEYS = (
    "BASE_DIR",
    "PROJECT_DIR",
    "REF_DIR",
    "REF_GENOME",
    "TSS_BED",
    "BLACKLIST_BED",
    "BLACKLIST_STATUS",
    "MITO_CONTIGS",
    "QC_THRESHOLD_PROFILE",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize ATAC-seq QC and compute bedtools tagAlign FRiP.",
        epilog=(
            "Common usage: python summarize_atac_qc.py <PROJECT_NAME> <REF_NAME>. "
            "Explicit --project-dir, --tss-bed, and --blacklist-bed are still supported "
            "and override config-derived paths."
        ),
    )
    parser.add_argument(
        "project_name_arg",
        nargs="?",
        metavar="PROJECT_NAME",
        help="Project name under the configured base dir, e.g. PRJNA728969_HEK293.",
    )
    parser.add_argument(
        "ref_name_arg",
        nargs="?",
        metavar="REF_NAME",
        help="Reference name configured in epi_config.sh, e.g. hg38_Ensembl.",
    )
    project_group = parser.add_mutually_exclusive_group()
    project_group.add_argument(
        "--project-name",
        help="Project name under --base-dir, e.g. PRJNA728969_HEK293.",
    )
    project_group.add_argument(
        "--project-dir",
        type=Path,
        help="Full path to an ATAC pipeline project directory.",
    )
    parser.add_argument(
        "--ref-name",
        help="Reference name configured in epi_config.sh. Used to resolve TSS/blacklist paths.",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=None,
        help=f"Base directory used with --project-name. Default: {DEFAULT_BASE_DIR}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output TSV path. Default: <project>/1_result/5_qc/atac_qc_summary.tsv",
    )
    parser.add_argument(
        "--mito-contigs",
        help="Comma/space separated mitochondrial contig names. Default: config MITO_CONTIGS or MT.",
    )
    parser.add_argument(
        "--skip-frip",
        action="store_true",
        help="Skip bedtools FRiP calculation and only summarize existing QC files.",
    )
    parser.add_argument(
        "--blacklist-bed",
        type=Path,
        help="Optional blacklist BED for blacklist fraction calculation.",
    )
    parser.add_argument(
        "--tss-bed",
        type=Path,
        help="Optional 1-bp TSS BED for TSS enrichment calculation.",
    )
    parser.add_argument(
        "--threshold-profile",
        help="QC profile in atac_qc_thresholds.json; overrides REF_NAME config.",
    )
    parser.add_argument(
        "--mode",
        choices=sorted(QC_MODES),
        default="release",
        help=(
            "QC scope. release preserves the complete current-pipeline gate; "
            "core evaluates principal ATAC evidence without release-only artefacts; "
            "existing-results is for historical outputs and never claims that old "
            "peaks were reprocessed by the current release workflow. Default: release."
        ),
    )
    parser.add_argument(
        "--run-manifest",
        type=Path,
        help=(
            "Optional TSV with historical per-run file paths. Paths are resolved "
            "relative to the manifest. Requires biosample_id; see script help/doc."
        ),
    )
    parser.add_argument(
        "--biosample-manifest",
        type=Path,
        help=(
            "Optional TSV with historical pooled/peak/fragment paths. Paths are "
            "resolved relative to the manifest. Requires biosample_id."
        ),
    )
    args = parser.parse_args()

    if args.project_name_arg and (args.project_name or args.project_dir):
        parser.error("Do not combine positional PROJECT_NAME with --project-name or --project-dir.")
    if args.ref_name_arg and args.ref_name:
        parser.error("Do not combine positional REF_NAME with --ref-name.")
    if args.project_name_arg and not (args.ref_name_arg or args.ref_name):
        parser.error("Positional usage requires both PROJECT_NAME and REF_NAME.")
    if not (args.project_name_arg or args.project_name or args.project_dir):
        parser.error("Provide PROJECT_NAME REF_NAME, --project-name, or --project-dir.")

    if args.project_name_arg:
        args.project_name = args.project_name_arg
    if args.ref_name_arg:
        args.ref_name = args.ref_name_arg

    if args.mode == "existing-results" and args.biosample_manifest is None:
        parser.error("--mode existing-results requires --biosample-manifest; do not infer legacy inputs")
    if args.mode == "release":
        if not args.ref_name:
            parser.error("--mode release requires --ref-name so release resources come from epi_config.sh")
        release_overrides = {
            "--blacklist-bed": args.blacklist_bed,
            "--tss-bed": args.tss_bed,
            "--threshold-profile": args.threshold_profile,
            "--mito-contigs": args.mito_contigs,
            "--skip-frip": args.skip_frip,
            "--run-manifest": args.run_manifest,
            "--biosample-manifest": args.biosample_manifest,
        }
        disallowed = [flag for flag, value in release_overrides.items() if value]
        if disallowed:
            parser.error(
                "release QC does not accept configuration/input overrides: "
                + ", ".join(disallowed)
            )

    return args


def load_ref_config(ref_name: Optional[str], project_name: Optional[str]) -> dict[str, str]:
    if not ref_name:
        return {}
    if not CONFIG_PATH.exists():
        raise RuntimeError(f"config file does not exist: {CONFIG_PATH}")

    commands = [
        "set -euo pipefail",
        f"export REF_NAME={shlex.quote(ref_name)}",
    ]
    if project_name:
        commands.append(f"export PROJECT_NAME={shlex.quote(project_name)}")
    else:
        commands.append("unset PROJECT_NAME")
    commands.extend(
        [
            f"source {shlex.quote(str(CONFIG_PATH))}",
            f"for key in {' '.join(CONFIG_KEYS)}; do",
            "  printf '%s\\t%s\\n' \"$key\" \"${!key-}\"",
            "done",
        ]
    )

    result = subprocess.run(
        ["bash", "-c", "\n".join(commands)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or "failed to load epi_config.sh"
        raise RuntimeError(message)

    config: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if "\t" not in line:
            continue
        key, value = line.split("\t", 1)
        config[key] = value
    return config


def _path_or_none(value: Optional[str]) -> Optional[Path]:
    if not value:
        return None
    return Path(value)


RUN_MANIFEST_PATH_COLUMNS = frozenset(
    {
        "fastp_json",
        "raw_bam",
        "raw_mito_qc",
        "dup_metrics",
        "pre_mito_idxstats",
        "final_idxstats",
        "final_bam",
    }
)
BIOSAMPLE_MANIFEST_PATH_COLUMNS = frozenset(
    {
        "tagalign_file",
        "fragment_file",
        "peak_file",
        "fragment_bams",
        "reproducibility_json",
        "fingerprint_json",
    }
)


def _manifest_path(value: str, manifest_dir: Path) -> str:
    """Return a normalized path, resolving relative values beside a manifest."""

    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest_dir / path
    return str(path.resolve())


def _manifest_evidence(value: str, manifest_dir: Path) -> str:
    """Keep a literal evidence note or normalize a likely evidence-file path.

    ``literal:...`` is an escape hatch for text that happens to contain a
    slash.  JSON/TXT/TSV and path-like values are resolved relative to the
    manifest, which lets a run manifest point at a future
    ``filter_provenance.json`` sidecar without introducing another CLI option.
    """

    value = value.strip()
    if value.startswith("literal:"):
        return value[len("literal:"):].strip()
    suffix = Path(value).suffix.lower()
    if "/" in value or "\\" in value or suffix in {".json", ".txt", ".tsv", ".log"}:
        return _manifest_path(value, manifest_dir)
    return value


def _safe_manifest_id(value: str, label: str, path: Path, line_number: int) -> str:
    if not value or not SAFE_ID.fullmatch(value):
        raise ValueError(f"{path}:{line_number}: unsafe or empty {label}: {value!r}")
    return value


def _manifest_values(value: str) -> list[str]:
    return split_values(value.replace(";", ",")) if value else []


def _read_manifest_rows(path: Path, kind: str) -> list[dict[str, str]]:
    """Read a tabular manifest and normalize declared path values.

    Both manifests deliberately use TSV so that they remain reviewable and do
    not require a YAML/JSON dependency on the cluster.  Empty/comment rows are
    ignored, while duplicate biosample rows are rejected by the specific
    readers below.
    """

    if not path.exists():
        raise ValueError(f"{kind} manifest does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{kind} manifest has no header: {path}")
        headers = [header.strip() if header else "" for header in reader.fieldnames]
        if "biosample_id" not in headers:
            raise ValueError(f"{kind} manifest requires a biosample_id column: {path}")
        rows: list[dict[str, str]] = []
        path_columns = (
            RUN_MANIFEST_PATH_COLUMNS if kind == "run" else BIOSAMPLE_MANIFEST_PATH_COLUMNS
        )
        for line_number, raw_row in enumerate(reader, start=2):
            row = {
                (key or "").strip(): (value or "").strip()
                for key, value in raw_row.items()
            }
            if not any(row.values()) or row.get("biosample_id", "").startswith("#"):
                continue
            row["biosample_id"] = _safe_manifest_id(
                row.get("biosample_id", ""), "biosample_id", path, line_number
            )
            for key in path_columns:
                if row.get(key):
                    if key == "fragment_bams":
                        row[key] = ";".join(
                            _manifest_path(item, path.parent)
                            for item in _manifest_values(row[key])
                        )
                    else:
                        row[key] = _manifest_path(row[key], path.parent)
            if kind == "run" and row.get("mapq_filter_evidence"):
                row["mapq_filter_evidence"] = _manifest_evidence(
                    row["mapq_filter_evidence"], path.parent
                )
            row["_manifest_path"] = str(path.resolve())
            row["_manifest_line"] = str(line_number)
            rows.append(row)
    return rows


def read_run_manifest(path: Path) -> "OrderedDict[str, list[dict[str, str]]]":
    """Read optional per-run paths used to QC historical result directories.

    Accepted columns are ``biosample_id`` plus any of ``run_id``,
    ``biological_replicate``, ``layout``, ``fastp_json``, ``raw_bam``,
    ``raw_mito_qc``, ``dup_metrics``, ``pre_mito_idxstats``,
    ``final_idxstats``, ``final_bam``, ``mapq_filter_threshold`` and
    ``mapq_filter_evidence``.  ``run_id`` is optional only to support old
    aggregated outputs; a deterministic manifest-local identifier is assigned
    when it is absent.
    """

    grouped: "OrderedDict[str, list[dict[str, str]]]" = OrderedDict()
    seen_runs: dict[str, str] = {}
    for index, row in enumerate(_read_manifest_rows(path, "run"), start=1):
        line_number = int(row["_manifest_line"])
        run_id = row.get("run_id", "") or f"manifest_run_{index}"
        row["run_id"] = _safe_manifest_id(run_id, "run_id", path, line_number)
        previous_biosample = seen_runs.get(run_id)
        if previous_biosample and previous_biosample != row["biosample_id"]:
            raise ValueError(
                f"{path}:{line_number}: run {run_id!r} is already assigned to "
                f"{previous_biosample}"
            )
        seen_runs[run_id] = row["biosample_id"]
        layout = row.get("layout", "").upper() or "AUTO"
        if layout not in VALID_LAYOUTS:
            raise ValueError(f"{path}:{line_number}: layout must be PE, SE, or AUTO; got {layout!r}")
        row["layout"] = layout
        replicate = row.get("biological_replicate", "").strip()
        if replicate:
            replicate = replicate if replicate.lower().startswith("rep") else f"rep{replicate}"
            row["biological_replicate"] = _safe_manifest_id(
                replicate, "biological_replicate", path, line_number
            )
            row["_legacy_replicate"] = "false"
        else:
            row["biological_replicate"] = "rep1"
            row["_legacy_replicate"] = "true"
        grouped.setdefault(row["biosample_id"], []).append(row)
    return grouped


def read_biosample_manifest(path: Path) -> "OrderedDict[str, dict[str, str]]":
    """Read optional pooled/peak paths used to QC historical results.

    Accepted columns are ``biosample_id`` plus any of ``run_ids``,
    ``biological_replicates``, ``layout``, ``tagalign_file``,
    ``fragment_file``, ``peak_file``, ``fragment_bams``,
    ``reproducibility_json``, ``fingerprint_json``, ``peak_provenance``,
    ``blacklist_peak_filter_status`` and ``tagalign_tn5_shifted``. The last
    value is a provenance declaration (for example ``true``), not a path.
    """

    result: "OrderedDict[str, dict[str, str]]" = OrderedDict()
    for row in _read_manifest_rows(path, "biosample"):
        biosample = row["biosample_id"]
        if biosample in result:
            raise ValueError(
                f"{path}:{row['_manifest_line']}: duplicate biosample_id {biosample!r}"
            )
        layout = row.get("layout", "").upper()
        if layout and layout not in VALID_LAYOUTS:
            raise ValueError(
                f"{path}:{row['_manifest_line']}: layout must be PE, SE, or AUTO; got {layout!r}"
            )
        row["layout"] = layout or "AUTO"
        result[biosample] = row
    return result


def project_dir_from_args(args: argparse.Namespace, ref_config: dict[str, str]) -> Path:
    if args.project_dir:
        return args.project_dir.resolve()
    if args.project_name:
        if args.base_dir:
            return (args.base_dir / args.project_name).resolve()
        if ref_config.get("PROJECT_DIR"):
            return Path(ref_config["PROJECT_DIR"]).resolve()
        if ref_config.get("BASE_DIR"):
            return (Path(ref_config["BASE_DIR"]) / args.project_name).resolve()
        return (DEFAULT_BASE_DIR / args.project_name).resolve()
    raise ValueError("project directory could not be resolved")


def resolve_mito_contigs(args: argparse.Namespace, ref_config: dict[str, str]) -> list[str]:
    if args.mito_contigs:
        return split_values(args.mito_contigs)
    if ref_config.get("MITO_CONTIGS"):
        return split_values(ref_config["MITO_CONTIGS"])
    return list(DEFAULT_MITO_CONTIGS)


def read_raw_mito_qc(path: Path) -> tuple[Optional[int], Optional[int]]:
    if not path.exists():
        return None, None
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return None, None
    try:
        return int(rows[0]["mapped_reads"]), int(rows[0]["mitochondrial_reads"])
    except (KeyError, TypeError, ValueError):
        return None, None


def _manifest_file(
    record: Optional[Mapping[str, str]], key: str, default: Path
) -> Path:
    """Use an explicitly declared legacy path or retain the standard layout."""

    value = record.get(key, "") if record else ""
    return Path(value) if value else default


def _mapq_value_from_json(value: object) -> Optional[str]:
    """Find common MAPQ-threshold keys in a small provenance JSON payload."""

    keys = {
        "mapq_filter_threshold",
        "mapq_threshold",
        "minimum_mapping_quality",
        "minimum_mapq",
        "mapq_min",
    }
    if isinstance(value, dict):
        for key, nested in value.items():
            if str(key).lower() in keys and nested not in (None, ""):
                return str(nested)
        for nested in value.values():
            found = _mapq_value_from_json(nested)
            if found is not None:
                return found
    if isinstance(value, list):
        for nested in value:
            found = _mapq_value_from_json(nested)
            if found is not None:
                return found
    return None


MAPQ_PROVENANCE_FIELDS = (
    "mapq_filter_operator",
    "samtools_view_exclude_flag",
    "samtools_view_required_flag",
    "raw_bam_alignment_record_count",
    "post_mapq_filter_alignment_record_count",
    "post_final_filter_alignment_record_count",
)


def _mapq_details_from_json(value: object) -> dict[str, str]:
    """Extract the exact run-level filtering fields written by ATAC_align.sh."""

    if not isinstance(value, dict):
        return {}
    return {
        key: str(value[key])
        for key in MAPQ_PROVENANCE_FIELDS
        if value.get(key) not in (None, "")
    }


def _mapq_details_from_evidence_tsv(path: Path, run_id: str) -> tuple[dict[str, str], str]:
    """Extract one run's filtering counts from a legacy evidence TSV.

    A current run writes a per-run ``*.mapping_filter.provenance.json``.
    Legacy projects instead carry a single project-level evidence TSV whose
    rows are keyed by ``run_id``; the same path is therefore declared for every
    run in the legacy manifest, so the counts must be looked up per run rather
    than summed over the whole file.

    Returns the provenance details plus a warning string when the run has no
    usable counts, so a partially readable file cannot silently bias retention.
    """

    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return {}, f"mapq_evidence_unreadable:{path}"

    exclude_flag = ""
    header: Optional[list[str]] = None
    row: Optional[dict[str, str]] = None
    for line in lines:
        if line.startswith("#"):
            match = re.search(r"(?:^|\s)-F\s+(\d+)", line)
            if match and not exclude_flag:
                exclude_flag = match.group(1)
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        if fields and fields[0] == run_id:
            row = dict(zip(header, fields))
            break
    if row is None:
        return {}, f"mapq_evidence_run_missing:{run_id}:{path}"

    if row.get("layout") == "PE":
        post_mapq_records = row.get("tmp_records", "")
        required_flag = "2"
    else:
        post_mapq_records = row.get("filt_records", "")
        required_flag = "null"
    counts = {
        "raw_bam_alignment_record_count": row.get("raw_total", ""),
        "post_mapq_filter_alignment_record_count": post_mapq_records,
        "post_final_filter_alignment_record_count": row.get("filt_records", ""),
    }
    if any(not value.isdigit() for value in counts.values()):
        return {}, f"mapq_evidence_count_unavailable:{run_id}:{path}"

    details = {
        "mapq_filter_operator": ">=",
        "samtools_view_exclude_flag": exclude_flag,
        "samtools_view_required_flag": required_flag,
        **counts,
    }
    return details, ""


def _read_mapq_evidence(value: str) -> tuple[str, Optional[str], str]:
    """Return displayed evidence, optional threshold extracted from JSON, and status."""

    if not value:
        return "", None, ""
    path = Path(value)
    if path.exists() and path.is_file():
        extracted: Optional[str] = None
        if path.suffix.lower() == ".json":
            try:
                extracted = _mapq_value_from_json(json.loads(path.read_text(encoding="utf-8")))
            except (OSError, ValueError, json.JSONDecodeError):
                # A declared file is still useful provenance even when its
                # schema is not one this QC script knows how to parse.
                extracted = None
        return str(path), extracted, "EVIDENCE_FILE"
    looks_like_path = (
        "/" in value
        or "\\" in value
        or path.suffix.lower() in {".json", ".txt", ".tsv", ".log"}
    )
    if looks_like_path:
        return value, None, "EVIDENCE_PATH_MISSING"
    return value, None, "DECLARED"


def _summarize_mapq_provenance(
    records: Sequence[Optional[Mapping[str, str]]],
    qc_mode: str,
) -> dict[str, object]:
    """Summarize optional MAPQ-filter evidence without inferring it from BAMs."""

    thresholds: list[str] = []
    evidence_values: list[str] = []
    statuses: list[str] = []
    warnings: list[str] = []
    provenance_values: dict[str, list[str]] = {
        key: [] for key in MAPQ_PROVENANCE_FIELDS
    }
    for record in records:
        if not record:
            continue
        declared = record.get("mapq_filter_threshold", "").strip()
        evidence, extracted, evidence_status = _read_mapq_evidence(
            record.get("mapq_filter_evidence", "").strip()
        )
        if declared:
            thresholds.append(declared)
        elif extracted:
            thresholds.append(extracted)
        if evidence:
            evidence_values.append(evidence)
            statuses.append(evidence_status)
            if evidence_status == "EVIDENCE_FILE":
                suffix = Path(evidence).suffix.lower()
                if suffix == ".json":
                    try:
                        payload = json.loads(Path(evidence).read_text(encoding="utf-8"))
                    except (OSError, ValueError, json.JSONDecodeError):
                        payload = {}
                    for key, detail in _mapq_details_from_json(payload).items():
                        provenance_values[key].append(detail)
                elif suffix == ".tsv":
                    details, detail_warning = _mapq_details_from_evidence_tsv(
                        Path(evidence), str(record.get("run_id", ""))
                    )
                    if detail_warning:
                        warnings.append(detail_warning)
                    for key, detail in details.items():
                        provenance_values[key].append(detail)

    unique_thresholds = list(dict.fromkeys(thresholds))
    unique_evidence = list(dict.fromkeys(evidence_values))
    if "EVIDENCE_PATH_MISSING" in statuses:
        warnings.append("mapq_filter_evidence_path_missing")
    if unique_thresholds:
        for value in unique_thresholds:
            try:
                float(value)
            except ValueError:
                warnings.append(f"mapq_filter_threshold_not_numeric:{value}")
        if len(unique_thresholds) > 1:
            warnings.append("mapq_filter_threshold_inconsistent_across_runs")
        if "EVIDENCE_FILE" in statuses:
            status = "EVIDENCE_FILE"
        elif "EVIDENCE_PATH_MISSING" in statuses:
            status = "EVIDENCE_PATH_MISSING"
        elif unique_evidence:
            status = "DECLARED"
        else:
            status = "DECLARED_THRESHOLD_NO_EVIDENCE"
            warnings.append("mapq_filter_threshold_declared_without_evidence")
    elif unique_evidence:
        if "EVIDENCE_PATH_MISSING" in statuses:
            status = "EVIDENCE_PATH_MISSING"
        else:
            status = "EVIDENCE_DECLARED_NO_THRESHOLD"
            warnings.append("mapq_filter_evidence_declared_without_threshold")
    elif qc_mode == "existing-results":
        status = "NOT_AVAILABLE_LEGACY"
        warnings.append("mapq_filter_not_available_legacy")
    else:
        status = "NOT_DECLARED"
        warnings.append("mapq_filter_not_declared")
    def joined_values(key: str) -> str:
        return ";".join(dict.fromkeys(provenance_values[key]))

    def summed_record_counts(key: str) -> Optional[int]:
        values = provenance_values[key]
        if not values:
            return None
        try:
            parsed = [int(value) for value in values]
        except ValueError:
            warnings.append(f"mapq_provenance_count_not_integer:{key}")
            return None
        return sum(parsed)

    raw_records = summed_record_counts("raw_bam_alignment_record_count")
    post_mapq_records = summed_record_counts("post_mapq_filter_alignment_record_count")
    post_final_records = summed_record_counts("post_final_filter_alignment_record_count")
    return {
        "mapq_filter_threshold": ";".join(unique_thresholds),
        "mapq_filter_evidence": ";".join(unique_evidence),
        "mapq_filter_status": status,
        "mapq_filter_operator": joined_values("mapq_filter_operator"),
        "mapq_excluded_sam_flags": joined_values("samtools_view_exclude_flag"),
        "mapq_required_sam_flags": joined_values("samtools_view_required_flag"),
        "mapq_raw_alignment_records": raw_records,
        "mapq_post_filter_alignment_records": post_mapq_records,
        "mapq_post_final_filter_alignment_records": post_final_records,
        "mapq_filter_retention": (
            post_mapq_records / raw_records
            if raw_records is not None and raw_records > 0 and post_mapq_records is not None
            else None
        ),
        "warnings": warnings,
    }


def collect_run_qc(
    project_dir: Path,
    align_dir: Path,
    run_ids: Sequence[str],
    mito_contigs: Sequence[str],
    replicate_ids: Sequence[str] = (),
    biosample: str = "",
    declared_layout: str = "",
    run_manifest_records: Optional[Mapping[str, Mapping[str, str]]] = None,
    qc_mode: str = "release",
    artifact_dir: Optional[Path] = None,
) -> dict[str, object]:
    layouts: list[str] = []
    missing: list[str] = []
    fastp_raw_reads: Optional[int] = None
    fastp_clean_reads: Optional[int] = None
    examined_reads: Optional[int] = None
    duplicate_reads: Optional[int] = None
    optical_duplicate_reads: Optional[int] = None
    raw_bam_total_reads: Optional[int] = None
    raw_bam_mapped_reads: Optional[int] = None
    raw_bam_paired_reads: Optional[int] = None
    raw_bam_proper_pair_reads: Optional[int] = None
    nodup_pre_mito_reads: Optional[int] = None
    final_nodup_reads: Optional[int] = None
    mitochondrial_contig_reads: Optional[int] = None
    warnings: list[str] = []
    raw_mito_mapped_reads: Optional[int] = None
    raw_mito_reads: Optional[int] = None
    gc_bias_statuses: list[str] = []
    manifest_records: list[Optional[Mapping[str, str]]] = []

    for run_id in run_ids:
        manifest_record = run_manifest_records.get(run_id) if run_manifest_records else None
        bowtie_dir = align_dir / run_id / "bowtie2"
        # Current/reanalysis jobs may emit a per-run provenance sidecar without
        # requiring a historical manifest.  Discover only explicit, narrowly
        # named candidates; never infer a MAPQ threshold from a filtered BAM.
        if not manifest_record:
            for candidate in (
                bowtie_dir / f"{run_id}.mapping_filter.provenance.json",
                bowtie_dir / f"{run_id}.filter_provenance.json",
                bowtie_dir / "filter_provenance.json",
            ):
                if candidate.exists():
                    manifest_record = {"mapq_filter_evidence": str(candidate)}
                    break
        manifest_records.append(manifest_record)
        manifest_layout = manifest_record.get("layout", "AUTO") if manifest_record else "AUTO"
        layouts.append(manifest_layout if manifest_layout != "AUTO" else detect_layout(project_dir, run_id))

        fastp_json = _manifest_file(
            manifest_record, "fastp_json", align_dir / run_id / "fastp" / f"{run_id}_fastp.json"
        )
        raw_reads, clean_reads = read_fastp_reads(fastp_json)
        if not fastp_json.exists():
            missing.append(str(fastp_json))
        fastp_raw_reads = add_optional(fastp_raw_reads, raw_reads)
        fastp_clean_reads = add_optional(fastp_clean_reads, clean_reads)

        raw_bam = _manifest_file(manifest_record, "raw_bam", bowtie_dir / f"{run_id}.raw.bam")
        raw_alignment_qc = summarize_raw_bam_flagstat(raw_bam)
        raw_bam_total_reads_value = cast(Optional[int], raw_alignment_qc["raw_bam_total_reads"])
        raw_bam_mapped_reads_value = cast(Optional[int], raw_alignment_qc["raw_bam_mapped_reads"])
        raw_bam_paired_reads_value = cast(Optional[int], raw_alignment_qc["raw_bam_paired_reads"])
        raw_bam_proper_pair_reads_value = cast(
            Optional[int], raw_alignment_qc["raw_bam_proper_pair_reads"]
        )
        raw_alignment_warnings = cast(list[str], raw_alignment_qc["warnings"])
        raw_bam_total_reads = add_optional(
            raw_bam_total_reads,
            raw_bam_total_reads_value,
        )
        raw_bam_mapped_reads = add_optional(
            raw_bam_mapped_reads,
            raw_bam_mapped_reads_value,
        )
        raw_bam_paired_reads = add_optional(
            raw_bam_paired_reads,
            raw_bam_paired_reads_value,
        )
        raw_bam_proper_pair_reads = add_optional(
            raw_bam_proper_pair_reads,
            raw_bam_proper_pair_reads_value,
        )
        warnings.extend(raw_alignment_warnings)

        raw_mito_file = _manifest_file(
            manifest_record, "raw_mito_qc", bowtie_dir / f"{run_id}.raw.mito.qc.tsv"
        )
        raw_mapped, raw_mito = read_raw_mito_qc(raw_mito_file)
        if not raw_mito_file.exists():
            missing.append(str(raw_mito_file))
        raw_mito_mapped_reads = add_optional(raw_mito_mapped_reads, raw_mapped)
        raw_mito_reads = add_optional(raw_mito_reads, raw_mito)

        # Legacy fallback. New phase-1 jobs overwrite these values below with
        # replicate-level metrics produced after technical-run merging.
        metrics_file = _manifest_file(manifest_record, "dup_metrics", bowtie_dir / f"{run_id}.dup.qc")
        run_examined, run_duplicates = read_picard_dup_metrics(metrics_file)
        if not metrics_file.exists() and not replicate_ids:
            missing.append(str(metrics_file))
        examined_reads = add_optional(examined_reads, run_examined)
        duplicate_reads = add_optional(duplicate_reads, run_duplicates)

        pre_mito_idxstats = _manifest_file(
            manifest_record,
            "pre_mito_idxstats",
            bowtie_dir / f"{run_id}.nodup.with_mito.bam.idxstats",
        )
        pre_total, pre_mito = read_idxstats(pre_mito_idxstats, mito_contigs)
        if pre_mito_idxstats.exists():
            nodup_pre_mito_reads = add_optional(nodup_pre_mito_reads, pre_total)
            mitochondrial_contig_reads = add_optional(mitochondrial_contig_reads, pre_mito)

        final_idxstats = _manifest_file(
            manifest_record, "final_idxstats", bowtie_dir / f"{run_id}.nodup.bam.idxstats"
        )
        final_total, _ = read_idxstats(final_idxstats, mito_contigs)
        if not final_idxstats.exists() and not replicate_ids:
            missing.append(str(final_idxstats))
        final_nodup_reads = add_optional(final_nodup_reads, final_total)

    layout = declared_layout or summarize_layout(layouts)
    if layout == "unknown" and raw_bam_paired_reads and raw_bam_paired_reads > 0:
        layout = "PE"

    if replicate_ids and biosample:
        examined_reads = duplicate_reads = nodup_pre_mito_reads = final_nodup_reads = mitochondrial_contig_reads = None
        optical_duplicate_reads = None
        complexity_rows: list[dict[str, object]] = []
        for replicate in replicate_ids:
            rep_dir = align_dir / "replicates" / biosample / replicate
            prefix = f"{biosample}.{replicate}"
            metrics_file = rep_dir / f"{prefix}.dup.qc"
            dup_details = read_picard_dup_details(metrics_file)
            rep_examined = cast(Optional[int], dup_details["examined_reads"])
            rep_duplicates = cast(Optional[int], dup_details["duplicate_reads"])
            examined_reads = add_optional(examined_reads, rep_examined)
            duplicate_reads = add_optional(duplicate_reads, rep_duplicates)
            optical_duplicate_reads = add_optional(
                optical_duplicate_reads, cast(Optional[int], dup_details["optical_duplicate_reads"])
            )
            pre_file = rep_dir / f"{prefix}.nodup.with_mito.bam.idxstats"
            final_file = rep_dir / f"{prefix}.nodup.bam.idxstats"
            pre_total, pre_mito = read_idxstats(pre_file, mito_contigs)
            final_total, _ = read_idxstats(final_file, mito_contigs)
            nodup_pre_mito_reads = add_optional(nodup_pre_mito_reads, pre_total)
            mitochondrial_contig_reads = add_optional(mitochondrial_contig_reads, pre_mito)
            final_nodup_reads = add_optional(final_nodup_reads, final_total)
            for required_file in (metrics_file, pre_file, final_file):
                if not required_file.exists():
                    missing.append(str(required_file))
            # A QC rerun must never create or replace a sidecar below the
            # historical alignment directory.  In particular,
            # ``existing-results`` is permitted to read old BAMs but writes
            # every newly derived artifact below the mode-specific QC root.
            complexity_file = (
                artifact_dir / "library_complexity" / biosample / f"{prefix}.lib_complexity.qc"
                if artifact_dir is not None
                else rep_dir / f"{prefix}.lib_complexity.qc"
            )
            complexity = compute_library_complexity(
                rep_dir / f"{prefix}.filt.bam",
                layout,
                complexity_file,
            )
            complexity_rows.append(complexity)
            warnings.extend(cast(list[str], complexity["warnings"]))
            gc_status_file = rep_dir / f"{prefix}.gc_bias.status.txt"
            if gc_status_file.exists():
                gc_bias_statuses.append(gc_status_file.read_text(encoding="utf-8").strip().split(":", 1)[0])
            else:
                missing.append(str(gc_status_file))
    else:
        warnings.append("legacy_run_level_duplicate_metrics")
        complexity_rows = []

    if gc_bias_statuses and all(value == "PASS" for value in gc_bias_statuses):
        gc_bias_status = "PASS"
    elif gc_bias_statuses:
        gc_bias_status = "WARN"
    else:
        gc_bias_status = None

    duplication_fraction = None
    if examined_reads and duplicate_reads is not None:
        duplication_fraction = duplicate_reads / examined_reads
    optical_duplicate_fraction = None
    if examined_reads and optical_duplicate_reads is not None:
        optical_duplicate_fraction = optical_duplicate_reads / examined_reads

    def sum_complexity(key: str) -> Optional[int]:
        values = [cast(Optional[int], row.get(key)) for row in complexity_rows]
        clean = [value for value in values if value is not None]
        return sum(clean) if clean and len(clean) == len(values) else None

    def worst_complexity(key: str) -> Optional[float]:
        values = [cast(Optional[float], row.get(key)) for row in complexity_rows]
        clean = [value for value in values if value is not None]
        return min(clean) if clean and len(clean) == len(values) else None

    mapping_rate = None
    if raw_bam_total_reads and raw_bam_mapped_reads is not None:
        mapping_rate = raw_bam_mapped_reads / raw_bam_total_reads

    proper_pair_rate = None
    if raw_bam_paired_reads and raw_bam_proper_pair_reads is not None:
        proper_pair_rate = raw_bam_proper_pair_reads / raw_bam_paired_reads

    raw_mito_fraction = None
    if raw_mito_mapped_reads and raw_mito_reads is not None:
        raw_mito_fraction = raw_mito_reads / raw_mito_mapped_reads

    final_usable_fraction = None
    if fastp_clean_reads and final_nodup_reads is not None:
        final_usable_fraction = final_nodup_reads / fastp_clean_reads

    mito_removed_reads = None
    mito_filter_fraction = None
    if nodup_pre_mito_reads is not None and final_nodup_reads is not None:
        mito_removed_reads = max(nodup_pre_mito_reads - final_nodup_reads, 0)
        if nodup_pre_mito_reads > 0:
            mito_filter_fraction = mito_removed_reads / nodup_pre_mito_reads

    mapq_qc = _summarize_mapq_provenance(manifest_records, qc_mode)
    warnings.extend(cast(list[str], mapq_qc["warnings"]))

    return {
        "layout": layout,
        "fastp_raw_reads": fastp_raw_reads,
        "fastp_clean_reads": fastp_clean_reads,
        "picard_examined_reads": examined_reads,
        "picard_duplicate_reads": duplicate_reads,
        "picard_optical_duplicate_reads": optical_duplicate_reads,
        "optical_duplicate_fraction": optical_duplicate_fraction,
        "duplication_fraction": duplication_fraction,
        "library_total_fragments": sum_complexity("library_total_fragments"),
        "library_distinct_fragments": sum_complexity("library_distinct_fragments"),
        "library_one_read_fragments": sum_complexity("library_one_read_fragments"),
        "library_two_read_fragments": sum_complexity("library_two_read_fragments"),
        "nrf": worst_complexity("nrf"),
        "pbc1": worst_complexity("pbc1"),
        "pbc2": worst_complexity("pbc2"),
        "library_complexity_file": ";".join(
            str(row["library_complexity_file"]) for row in complexity_rows
            if row.get("library_complexity_file")
        ),
        "gc_bias_status": gc_bias_status,
        "raw_bam_total_reads": raw_bam_total_reads,
        "raw_bam_mapped_reads": raw_bam_mapped_reads,
        "mapping_rate": mapping_rate,
        "raw_bam_paired_reads": raw_bam_paired_reads,
        "raw_bam_proper_pair_reads": raw_bam_proper_pair_reads,
        "proper_pair_rate": proper_pair_rate,
        "mapq_filter_threshold": mapq_qc["mapq_filter_threshold"],
        "mapq_filter_evidence": mapq_qc["mapq_filter_evidence"],
        "mapq_filter_status": mapq_qc["mapq_filter_status"],
        "mapq_filter_operator": mapq_qc["mapq_filter_operator"],
        "mapq_excluded_sam_flags": mapq_qc["mapq_excluded_sam_flags"],
        "mapq_required_sam_flags": mapq_qc["mapq_required_sam_flags"],
        "mapq_raw_alignment_records": mapq_qc["mapq_raw_alignment_records"],
        "mapq_post_filter_alignment_records": mapq_qc["mapq_post_filter_alignment_records"],
        "mapq_post_final_filter_alignment_records": mapq_qc["mapq_post_final_filter_alignment_records"],
        "mapq_filter_retention": mapq_qc["mapq_filter_retention"],
        "raw_mito_mapped_reads": raw_mito_mapped_reads,
        "raw_mito_reads": raw_mito_reads,
        "raw_mito_fraction": raw_mito_fraction,
        "nodup_pre_mito_reads": nodup_pre_mito_reads,
        "final_nodup_reads": final_nodup_reads,
        "final_usable_fraction": final_usable_fraction,
        "mitochondrial_contig_reads": mitochondrial_contig_reads,
        "mitochondrial_filter_removed_reads": mito_removed_reads,
        "mitochondrial_filter_fraction": mito_filter_fraction,
        "warnings": warnings,
        "missing_files": missing,
    }


def collect_peak_and_frip_qc(
    biosample: str,
    pooled_dir: Path,
    peaks_dir: Path,
    artifact_dir: Path,
    skip_frip: bool,
    blacklist_bed: Optional[Path],
    tss_bed: Optional[Path],
    tagalign_override: Optional[Path] = None,
    fragment_override: Optional[Path] = None,
    peak_override: Optional[Path] = None,
    qc_mode: str = "release",
    peak_provenance: str = "",
    blacklist_peak_filter_status: str = "",
    tagalign_tn5_shifted: str = "",
    frip_tool_error: str = "",
) -> dict[str, object]:
    tagalign_file = tagalign_override or (pooled_dir / f"{biosample}.tn5.tagAlign.gz")
    fragment_file = fragment_override or (pooled_dir / f"{biosample}.fragments.bed.gz")
    peak_file = peak_override or (peaks_dir / f"{biosample}_peaks.narrowPeak")
    missing: list[str] = []
    warnings: list[str] = []
    tagalign_total_reads = None
    reads_in_peaks = None
    frip_fraction = None
    frip_percent = None
    fragment_total = None
    fragments_in_peaks = None
    fragment_frip_fraction = None
    frip_input_unit = "fragment" if fragment_file.exists() else "tagAlign_record_legacy"
    frip_input = fragment_file if fragment_file.exists() else tagalign_file
    frip_qc_status = "SKIPPED" if skip_frip else "MISSING_INPUT"
    frip_contig_intersection = ""

    try:
        peak_qc = summarize_peak_file(peak_file)
    except (OSError, EOFError, UnicodeError, ValueError) as exc:
        peak_qc = {"peak_count": None, "peak_total_bp": None, "peak_median_width": None}
        warnings.append(f"peak_qc_failed:{exc}")
    if not tagalign_file.exists():
        missing.append(str(tagalign_file))
    if not peak_file.exists():
        missing.append(str(peak_file))

    if not skip_frip and frip_tool_error:
        frip_qc_status = "TOOLS_UNAVAILABLE"
        warnings.append(f"frip_tools_unavailable:{frip_tool_error}")
    elif not skip_frip and frip_input.exists() and peak_file.exists():
        try:
            shared_contigs = validate_frip_contigs(frip_input, peak_file)
        except ValueError as exc:
            frip_qc_status = "CONTIG_MISMATCH"
            warnings.append(f"frip_contig_mismatch:{exc}")
        except (OSError, EOFError, UnicodeError) as exc:
            frip_qc_status = "FAILED"
            warnings.append(f"frip_contig_check_failed:{exc}")
        else:
            frip_contig_intersection = ",".join(sorted(shared_contigs))
            try:
                total_records = count_all_records(frip_input)
                if total_records and total_records > 0:
                    overlap_records = compute_records_in_peaks(frip_input, peak_file)
                    frip_value = overlap_records / total_records
                    if fragment_file.exists():
                        fragment_total = total_records
                        fragments_in_peaks = overlap_records
                        fragment_frip_fraction = frip_value
                    else:
                        tagalign_total_reads = total_records
                        reads_in_peaks = overlap_records
                        frip_fraction = frip_value
                        warnings.append("legacy_tagalign_record_frip_not_fragment_frip")
                    frip_percent = frip_value * 100
                    frip_qc_status = "COMPUTED"
                else:
                    if fragment_file.exists():
                        fragment_total, fragments_in_peaks, fragment_frip_fraction = 0, 0, 0.0
                    else:
                        tagalign_total_reads, reads_in_peaks, frip_fraction = 0, 0, 0.0
                    frip_percent = 0.0
                    frip_qc_status = "NO_INPUT_RECORDS"
            except (OSError, EOFError, UnicodeError, ValueError, RuntimeError) as exc:
                frip_qc_status = "FAILED"
                warnings.append(f"frip_calculation_failed:{exc}")

    if blacklist_bed is not None and tagalign_total_reads is None and tagalign_file.exists():
        try:
            tagalign_total_reads = count_all_records(tagalign_file)
        except (OSError, EOFError, UnicodeError, ValueError) as exc:
            warnings.append(f"tagalign_record_count_failed:{exc}")

    blacklist_qc = compute_blacklist_qc(tagalign_file, blacklist_bed, tagalign_total_reads)
    warnings.extend(cast(list[str], blacklist_qc["warnings"]))

    # A current-pipeline tagAlign is named and generated as Tn5-shifted.  For
    # historical inputs the manifest must make that provenance explicit; the
    # score is still useful diagnostically when unknown, but cannot be treated
    # as fully method-comparable in the core QC status.
    declared_shift = tagalign_tn5_shifted.strip().lower()
    if qc_mode != "existing-results":
        tss_shift_status = "CURRENT_PIPELINE_TN5_SHIFTED"
    elif declared_shift in {"1", "true", "yes", "shifted", "tn5_shifted"}:
        tss_shift_status = "DECLARED_TN5_SHIFTED"
    elif declared_shift in {"0", "false", "no", "unshifted", "not_shifted"}:
        tss_shift_status = "DECLARED_NOT_TN5_SHIFTED"
        warnings.append("tss_tagalign_not_tn5_shifted_diagnostic_only")
    else:
        tss_shift_status = "UNKNOWN_LEGACY"
        warnings.append("tss_tagalign_shift_unknown_diagnostic_only")

    tss_profile_file = artifact_dir / "tss_enrichment" / f"{biosample}.tss_profile.tsv"
    tss_qc = compute_tss_enrichment(tagalign_file, tss_bed, tss_profile_file)
    warnings.extend(cast(list[str], tss_qc["warnings"]))

    # A post-hoc overlap calculation cannot retroactively turn a historical
    # peak file into a blacklist-filtered peak set.  Record this explicitly so
    # downstream Markdown/JSON reports cannot make that claim by accident.
    if not peak_provenance:
        peak_provenance = "legacy_existing_peak" if qc_mode == "existing-results" else "not_recorded"
    if not blacklist_peak_filter_status:
        blacklist_peak_filter_status = "NOT_RECALLED" if qc_mode == "existing-results" else "NOT_RECORDED"

    return {
        "tagalign_total_reads": tagalign_total_reads,
        "reads_in_peaks": reads_in_peaks,
        "frip_fraction": frip_fraction,
        "frip_percent": frip_percent,
        "peak_count": peak_qc["peak_count"],
        "peak_total_bp": peak_qc["peak_total_bp"],
        "peak_median_width": peak_qc["peak_median_width"],
        "fragment_total": fragment_total,
        "fragments_in_peaks": fragments_in_peaks,
        "fragment_frip_fraction": fragment_frip_fraction,
        "frip_input_unit": frip_input_unit,
        "frip_method": "" if skip_frip else (FRIP_METHOD if fragment_file.exists() else LEGACY_FRIP_METHOD),
        "frip_qc_status": frip_qc_status,
        "frip_contig_intersection": frip_contig_intersection,
        "blacklist_reads": blacklist_qc["blacklist_reads"],
        "blacklist_fraction": blacklist_qc["blacklist_fraction"],
        "blacklist_qc_status": blacklist_qc["blacklist_qc_status"],
        "blacklist_qc_scope": blacklist_qc["blacklist_qc_scope"],
        "blacklist_contig_intersection": blacklist_qc["blacklist_contig_intersection"],
        "peak_provenance": peak_provenance,
        "blacklist_peak_filter_status": blacklist_peak_filter_status,
        "tss_enrichment_score": tss_qc["tss_enrichment_score"],
        "tss_method": tss_qc["tss_method"],
        "tss_tagalign_shift_status": tss_shift_status,
        "tss_artifacts": tss_qc["tss_artifacts"],
        "tss_profile_file": tss_qc["tss_profile_file"],
        "tss_matrix_file": tss_qc["tss_matrix_file"],
        "tss_heatmap_file": tss_qc["tss_heatmap_file"],
        "warnings": warnings,
        "missing_files": missing,
        "tagalign_file": tagalign_file,
        "peak_file": peak_file,
    }


def collect_reproducibility_qc(
    project_dir: Path,
    biosample: str,
    reproducibility_override: Optional[Path] = None,
    fingerprint_override: Optional[Path] = None,
) -> dict[str, object]:
    path = reproducibility_override or (
        project_dir / "1_result" / "4_reproducibility" / biosample /
        f"{biosample}.reproducibility.json"
    )
    fingerprint_path = fingerprint_override or (
        project_dir / "1_result" / "5_qc" / "fingerprint" / biosample /
        f"{biosample}.fingerprint.json"
    )
    warnings: list[str] = []
    missing_files: list[str] = []
    fingerprint_status = None
    if fingerprint_path.exists():
        try:
            with fingerprint_path.open("r", encoding="utf-8") as handle:
                fingerprint_data = json.load(handle)
            if not isinstance(fingerprint_data, dict):
                raise ValueError("fingerprint JSON top level must be an object")
            fingerprint_status = fingerprint_data.get("status")
        except (OSError, json.JSONDecodeError, ValueError) as exc:
            fingerprint_status = "INVALID_JSON"
            warnings.append(f"fingerprint_json_invalid:{fingerprint_path}:{exc}")
    else:
        missing_files.append(str(fingerprint_path))
    if not path.exists():
        return {
            "self_consistency_ratio": None,
            "rescue_ratio": None,
            "reproducibility_status": None,
            "reproducibility_json": str(path),
            "fingerprint_status": fingerprint_status,
            "fingerprint_json": str(fingerprint_path),
            "peak_blacklist_applied": None,
            "peak_blacklist_status": "",
            "peak_blacklist_bed": "",
            "peak_blacklist_runtime_bed": "",
            "peak_blacklist_sha256": "",
            "peak_blacklist_source": "",
            "peak_pipeline_mode": "",
            "peak_blacklist_validation_status": "",
            "peak_blacklist_validation_record": "",
            "peak_blacklist_validation_bed": "",
            "peak_blacklist_validation_sha256": "",
            "peak_blacklist_validation_reference_fai": "",
            "peak_blacklist_validation_reference_fai_sha256": "",
            "peak_blacklist_validation_interval_count": None,
            "warnings": warnings + [f"reproducibility_json_missing:{path}"],
            "missing_files": [str(path)] + missing_files,
        }
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = json.load(handle)
        if not isinstance(data, dict):
            raise ValueError("reproducibility JSON top level must be an object")
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        return {
            "self_consistency_ratio": None,
            "rescue_ratio": None,
            "reproducibility_status": "INVALID_JSON",
            "reproducibility_json": str(path),
            "fingerprint_status": fingerprint_status,
            "fingerprint_json": str(fingerprint_path),
            "peak_blacklist_applied": None,
            "peak_blacklist_status": "",
            "peak_blacklist_bed": "",
            "peak_blacklist_runtime_bed": "",
            "peak_blacklist_sha256": "",
            "peak_blacklist_source": "",
            "peak_pipeline_mode": "",
            "peak_blacklist_validation_status": "",
            "peak_blacklist_validation_record": "",
            "peak_blacklist_validation_bed": "",
            "peak_blacklist_validation_sha256": "",
            "peak_blacklist_validation_reference_fai": "",
            "peak_blacklist_validation_reference_fai_sha256": "",
            "peak_blacklist_validation_interval_count": None,
            "warnings": warnings + [f"reproducibility_json_invalid:{path}:{exc}"],
            "missing_files": missing_files,
        }
    return {
        "self_consistency_ratio": data.get("self_consistency_ratio"),
        "rescue_ratio": data.get("rescue_ratio"),
        "reproducibility_status": data.get("reproducibility_status"),
        "reproducibility_json": str(path),
        "fingerprint_status": fingerprint_status,
        "fingerprint_json": str(fingerprint_path),
        "peak_blacklist_applied": data.get("blacklist_applied"),
        "peak_blacklist_status": data.get("blacklist_status") or "",
        "peak_blacklist_bed": data.get("blacklist_bed") or "",
        "peak_blacklist_runtime_bed": data.get("blacklist_runtime_bed") or "",
        "peak_blacklist_sha256": data.get("blacklist_sha256") or "",
        "peak_blacklist_source": data.get("blacklist_source") or "",
        "peak_pipeline_mode": data.get("peak_pipeline_mode") or "",
        "peak_blacklist_validation_status": data.get("blacklist_validation_status") or "",
        "peak_blacklist_validation_record": data.get("blacklist_validation_record") or "",
        "peak_blacklist_validation_bed": data.get("blacklist_validation_bed") or "",
        "peak_blacklist_validation_sha256": data.get("blacklist_validation_sha256") or "",
        "peak_blacklist_validation_reference_fai": (
            data.get("blacklist_validation_reference_fai") or ""
        ),
        "peak_blacklist_validation_reference_fai_sha256": (
            data.get("blacklist_validation_reference_fai_sha256") or ""
        ),
        "peak_blacklist_validation_interval_count": data.get("blacklist_validation_interval_count"),
        "warnings": warnings,
        "missing_files": missing_files,
    }


def format_summary_row(
    project_dir: Path,
    biosample: str,
    run_ids: Sequence[str],
    run_qc: dict[str, object],
    peak_frip_qc: dict[str, object],
    fragment_qc: dict[str, object],
    replicate_ids: Sequence[str],
    legacy_sample_map: bool,
    threshold_profile: str,
    reproducibility_qc: dict[str, object],
    qc_mode: str = "release",
    blacklist_config_status: str = "",
    blacklist_config_bed: str = "",
    blacklist_config_sha256: str = "",
    reference_config_fai: str = "",
    reference_config_fai_sha256: str = "",
) -> dict[str, object]:
    missing_files = list(cast(Sequence[str], run_qc["missing_files"])) + list(
        cast(Sequence[str], peak_frip_qc["missing_files"])
    ) + list(cast(Sequence[str], reproducibility_qc["missing_files"]))
    warnings = (
        list(cast(Sequence[str], run_qc["warnings"]))
        + list(cast(Sequence[str], peak_frip_qc["warnings"]))
        + list(cast(Sequence[str], fragment_qc["warnings"]))
        + list(cast(Sequence[str], reproducibility_qc["warnings"]))
    )
    blacklist_applied_value = reproducibility_qc["peak_blacklist_applied"]
    if isinstance(blacklist_applied_value, bool):
        peak_blacklist_applied = str(blacklist_applied_value).lower()
    elif blacklist_applied_value is None:
        peak_blacklist_applied = ""
    else:
        peak_blacklist_applied = str(blacklist_applied_value).lower()
    metrics: dict[str, object] = {
        "project_name": project_dir.name,
        "biosample_id": biosample,
        "biological_replicates": ",".join(replicate_ids),
        "biological_replicate_count": str(len(replicate_ids)),
        "legacy_sample_map": str(legacy_sample_map).lower(),
        "qc_mode": qc_mode,
        "threshold_profile": threshold_profile,
        "run_ids": ",".join(run_ids),
        "run_count": str(len(run_ids)),
        "layout": str(cast(object, run_qc["layout"])),
        "fastp_raw_reads": format_int(cast(Optional[int], run_qc["fastp_raw_reads"])),
        "fastp_clean_reads": format_int(cast(Optional[int], run_qc["fastp_clean_reads"])),
        "picard_examined_reads": format_int(cast(Optional[int], run_qc["picard_examined_reads"])),
        "picard_duplicate_reads": format_int(cast(Optional[int], run_qc["picard_duplicate_reads"])),
        "picard_optical_duplicate_reads": format_int(cast(Optional[int], run_qc["picard_optical_duplicate_reads"])),
        "optical_duplicate_fraction": format_float(cast(Optional[float], run_qc["optical_duplicate_fraction"])),
        "duplication_fraction": format_float(cast(Optional[float], run_qc["duplication_fraction"])),
        "library_total_fragments": format_int(cast(Optional[int], run_qc["library_total_fragments"])),
        "library_distinct_fragments": format_int(cast(Optional[int], run_qc["library_distinct_fragments"])),
        "library_one_read_fragments": format_int(cast(Optional[int], run_qc["library_one_read_fragments"])),
        "library_two_read_fragments": format_int(cast(Optional[int], run_qc["library_two_read_fragments"])),
        "nrf": format_float(cast(Optional[float], run_qc["nrf"])),
        "pbc1": format_float(cast(Optional[float], run_qc["pbc1"])),
        "pbc2": format_float(cast(Optional[float], run_qc["pbc2"])),
        "library_complexity_file": str(run_qc["library_complexity_file"]),
        "gc_bias_status": str(run_qc["gc_bias_status"] or ""),
        "raw_bam_total_reads": format_int(cast(Optional[int], run_qc["raw_bam_total_reads"])),
        "raw_bam_mapped_reads": format_int(cast(Optional[int], run_qc["raw_bam_mapped_reads"])),
        "mapping_rate": format_float(cast(Optional[float], run_qc["mapping_rate"])),
        "raw_bam_paired_reads": format_int(cast(Optional[int], run_qc["raw_bam_paired_reads"])),
        "raw_bam_proper_pair_reads": format_int(cast(Optional[int], run_qc["raw_bam_proper_pair_reads"])),
        "proper_pair_rate": format_float(cast(Optional[float], run_qc["proper_pair_rate"])),
        "mapq_filter_threshold": str(run_qc["mapq_filter_threshold"]),
        "mapq_filter_evidence": str(run_qc["mapq_filter_evidence"]),
        "mapq_filter_status": str(run_qc["mapq_filter_status"]),
        "mapq_filter_operator": str(run_qc["mapq_filter_operator"]),
        "mapq_excluded_sam_flags": str(run_qc["mapq_excluded_sam_flags"]),
        "mapq_required_sam_flags": str(run_qc["mapq_required_sam_flags"]),
        "mapq_raw_alignment_records": format_int(cast(Optional[int], run_qc["mapq_raw_alignment_records"])),
        "mapq_post_filter_alignment_records": format_int(cast(Optional[int], run_qc["mapq_post_filter_alignment_records"])),
        "mapq_post_final_filter_alignment_records": format_int(cast(Optional[int], run_qc["mapq_post_final_filter_alignment_records"])),
        "mapq_filter_retention": format_float(cast(Optional[float], run_qc["mapq_filter_retention"])),
        "raw_mito_mapped_reads": format_int(cast(Optional[int], run_qc["raw_mito_mapped_reads"])),
        "raw_mito_reads": format_int(cast(Optional[int], run_qc["raw_mito_reads"])),
        "raw_mito_fraction": format_float(cast(Optional[float], run_qc["raw_mito_fraction"])),
        "nodup_pre_mito_reads": format_int(cast(Optional[int], run_qc["nodup_pre_mito_reads"])),
        "final_nodup_reads": format_int(cast(Optional[int], run_qc["final_nodup_reads"])),
        "final_usable_fraction": format_float(cast(Optional[float], run_qc["final_usable_fraction"])),
        "mitochondrial_contig_reads": format_int(cast(Optional[int], run_qc["mitochondrial_contig_reads"])),
        "mitochondrial_filter_removed_reads": format_int(
            cast(Optional[int], run_qc["mitochondrial_filter_removed_reads"])
        ),
        "mitochondrial_filter_fraction": format_float(
            cast(Optional[float], run_qc["mitochondrial_filter_fraction"])
        ),
        "tagalign_total_reads": format_int(cast(Optional[int], peak_frip_qc["tagalign_total_reads"])),
        "reads_in_peaks": format_int(cast(Optional[int], peak_frip_qc["reads_in_peaks"])),
        "frip_fraction": format_float(cast(Optional[float], peak_frip_qc["frip_fraction"])),
        "frip_percent": format_float(cast(Optional[float], peak_frip_qc["frip_percent"]), digits=4),
        "fragment_total": format_int(cast(Optional[int], peak_frip_qc["fragment_total"])),
        "fragments_in_peaks": format_int(cast(Optional[int], peak_frip_qc["fragments_in_peaks"])),
        "fragment_frip_fraction": format_float(cast(Optional[float], peak_frip_qc["fragment_frip_fraction"])),
        "frip_input_unit": str(peak_frip_qc["frip_input_unit"]),
        "peak_count": format_int(cast(Optional[int], peak_frip_qc["peak_count"])),
        "peak_total_bp": format_int(cast(Optional[int], peak_frip_qc["peak_total_bp"])),
        "peak_median_width": format_float(cast(Optional[float], peak_frip_qc["peak_median_width"]), digits=1),
        "frip_method": str(peak_frip_qc["frip_method"]),
        "frip_qc_status": str(peak_frip_qc["frip_qc_status"]),
        "frip_contig_intersection": str(peak_frip_qc["frip_contig_intersection"]),
        "fragment_count": format_int(cast(Optional[int], fragment_qc["fragment_count"])),
        "fragment_mean_size": format_float(cast(Optional[float], fragment_qc["fragment_mean_size"]), digits=2),
        "fragment_median_size": format_int(cast(Optional[int], fragment_qc["fragment_median_size"])),
        "fragment_p10_size": format_int(cast(Optional[int], fragment_qc["fragment_p10_size"])),
        "fragment_p90_size": format_int(cast(Optional[int], fragment_qc["fragment_p90_size"])),
        "nfr_fraction": format_float(cast(Optional[float], fragment_qc["nfr_fraction"])),
        "mono_nucleosome_fraction": format_float(cast(Optional[float], fragment_qc["mono_nucleosome_fraction"])),
        "di_nucleosome_fraction": format_float(cast(Optional[float], fragment_qc["di_nucleosome_fraction"])),
        "fragment_periodicity_score": format_float(cast(Optional[float], fragment_qc["fragment_periodicity_score"])),
        "fragment_histogram_file": str(fragment_qc["fragment_histogram_file"]),
        "fragment_input_status": str(fragment_qc["fragment_input_status"]),
        "blacklist_reads": format_int(cast(Optional[int], peak_frip_qc["blacklist_reads"])),
        "blacklist_fraction": format_float(cast(Optional[float], peak_frip_qc["blacklist_fraction"])),
        "blacklist_qc_status": str(peak_frip_qc["blacklist_qc_status"]),
        "blacklist_qc_scope": str(peak_frip_qc["blacklist_qc_scope"]),
        "blacklist_contig_intersection": str(peak_frip_qc["blacklist_contig_intersection"]),
        "peak_provenance": str(peak_frip_qc["peak_provenance"]),
        "blacklist_peak_filter_status": str(peak_frip_qc["blacklist_peak_filter_status"]),
        "blacklist_config_status": blacklist_config_status,
        "blacklist_config_bed": blacklist_config_bed,
        "blacklist_config_sha256": blacklist_config_sha256,
        "reference_config_fai": reference_config_fai,
        "reference_config_fai_sha256": reference_config_fai_sha256,
        "peak_pipeline_mode": str(reproducibility_qc["peak_pipeline_mode"]),
        "peak_blacklist_applied": peak_blacklist_applied,
        "peak_blacklist_status": str(reproducibility_qc["peak_blacklist_status"]),
        "peak_blacklist_bed": str(reproducibility_qc["peak_blacklist_bed"]),
        "peak_blacklist_runtime_bed": str(reproducibility_qc["peak_blacklist_runtime_bed"]),
        "peak_blacklist_sha256": str(reproducibility_qc["peak_blacklist_sha256"]),
        "peak_blacklist_source": str(reproducibility_qc["peak_blacklist_source"]),
        "peak_blacklist_validation_status": str(
            reproducibility_qc["peak_blacklist_validation_status"]
        ),
        "peak_blacklist_validation_record": str(
            reproducibility_qc["peak_blacklist_validation_record"]
        ),
        "peak_blacklist_validation_bed": str(
            reproducibility_qc["peak_blacklist_validation_bed"]
        ),
        "peak_blacklist_validation_sha256": str(
            reproducibility_qc["peak_blacklist_validation_sha256"]
        ),
        "peak_blacklist_validation_reference_fai": str(
            reproducibility_qc["peak_blacklist_validation_reference_fai"]
        ),
        "peak_blacklist_validation_reference_fai_sha256": str(
            reproducibility_qc["peak_blacklist_validation_reference_fai_sha256"]
        ),
        "peak_blacklist_validation_interval_count": format_int(
            cast(Optional[int], reproducibility_qc["peak_blacklist_validation_interval_count"])
        ),
        "tss_enrichment_score": format_float(
            cast(Optional[float], peak_frip_qc["tss_enrichment_score"]),
            digits=4,
        ),
        "tss_method": str(peak_frip_qc["tss_method"]),
        "tss_tagalign_shift_status": str(peak_frip_qc["tss_tagalign_shift_status"]),
        "tss_artifacts": str(peak_frip_qc["tss_artifacts"]),
        "tss_profile_file": str(peak_frip_qc["tss_profile_file"]),
        "tss_matrix_file": str(peak_frip_qc["tss_matrix_file"]),
        "tss_heatmap_file": str(peak_frip_qc["tss_heatmap_file"]),
        "self_consistency_ratio": format_float(cast(Optional[float], reproducibility_qc["self_consistency_ratio"])),
        "rescue_ratio": format_float(cast(Optional[float], reproducibility_qc["rescue_ratio"])),
        "reproducibility_status": str(reproducibility_qc["reproducibility_status"] or ""),
        "reproducibility_json": str(reproducibility_qc["reproducibility_json"]),
        "fingerprint_status": str(reproducibility_qc["fingerprint_status"] or ""),
        "fingerprint_json": str(reproducibility_qc["fingerprint_json"]),
        "core_qc_completeness": "",
        "missing_files": ";".join(missing_files),
        "tagalign_file": str(peak_frip_qc["tagalign_file"]),
        "peak_file": str(peak_frip_qc["peak_file"]),
    }
    if legacy_sample_map:
        warnings.append("legacy_sample_map_biological_replicates_unknown")
    metrics.update(
        evaluate_qc_statuses(
            metrics,
            warnings,
            profile_name=threshold_profile,
            mode=qc_mode,
        )
    )
    return metrics


def summarize_project(
    project_dir: Path,
    output_file: Path,
    mito_contigs: Sequence[str],
    skip_frip: bool,
    blacklist_bed: Optional[Path],
    tss_bed: Optional[Path],
    threshold_profile: str,
    mode: str = "release",
    run_manifest: Optional[Path] = None,
    biosample_manifest: Optional[Path] = None,
    blacklist_config_status: str = "",
    blacklist_config_bed: str = "",
    blacklist_config_sha256: str = "",
    reference_config_fai: str = "",
    reference_config_fai_sha256: str = "",
) -> None:
    if mode not in QC_MODES:
        raise ValueError(f"Unknown QC mode: {mode}; expected one of {sorted(QC_MODES)}")
    if mode == "existing-results" and biosample_manifest is None:
        raise ValueError(
            "existing-results mode requires a biosample manifest; do not infer legacy inputs"
        )
    if mode != "existing-results" and (run_manifest is not None or biosample_manifest is not None):
        raise ValueError(
            "run/biosample manifests are supported only in existing-results mode; "
            "core and release QC must use the current sample_run_map.tsv and canonical outputs"
        )
    result_dir = project_dir / "1_result"
    align_dir = result_dir / "1_alignment"
    pooled_dir = result_dir / "2_tagalign"
    peaks_dir = result_dir / "3_peak_calling"
    qc_dir = result_dir / "5_qc"
    # Never share generated QC artifacts between historical, core and release
    # runs.  A baseline summary must keep pointing at the exact profile and
    # histogram from which it was interpreted.
    artifact_dir = qc_dir / "artifacts" / mode
    map_file = project_dir / "0_data" / "sample_run_map.tsv"

    run_manifest_by_sample = read_run_manifest(run_manifest) if run_manifest else OrderedDict()
    biosample_manifest_by_sample = (
        read_biosample_manifest(biosample_manifest) if biosample_manifest else OrderedDict()
    )

    # A supplied run manifest is authoritative for historical layouts.  An
    # existing-results baseline must never fall back to the current project
    # sample map: it may describe a later rerun in the same project directory.
    # Without a legacy run manifest, only the biosample-manifest ``run_ids``
    # are eligible for historical run-level QC.
    if run_manifest:
        records: list[ReplicateRecord] = []
        for biosample, manifest_rows in run_manifest_by_sample.items():
            for row in manifest_rows:
                records.append(
                    ReplicateRecord(
                        biosample,
                        row["biological_replicate"],
                        row["run_id"],
                        row["layout"],
                        row["_legacy_replicate"] == "true",
                    )
                )
    elif mode == "existing-results":
        records = []
    elif map_file.exists():
        records = read_replicate_records(map_file)
    else:
        raise RuntimeError(f"Current {mode} QC requires sample map: {map_file}")

    if mode == "existing-results":
        undeclared = set(run_manifest_by_sample) - set(biosample_manifest_by_sample)
        if undeclared:
            raise ValueError(
                "Each legacy run-manifest biosample must have one biosample-manifest row; "
                f"missing: {', '.join(sorted(undeclared))}"
            )
    grouped_replicates = group_replicates(records)
    sample_runs: "OrderedDict[str, list[str]]" = OrderedDict()
    for record in records:
        sample_runs.setdefault(record.biosample_id, [])
        if record.run_id not in sample_runs[record.biosample_id]:
            sample_runs[record.biosample_id].append(record.run_id)

    # A biosample manifest may describe pooled historical data that has no
    # recoverable run-level files.  It still deserves a row with its measured
    # TSS/FRiP/fragment evidence rather than being discarded by discovery.
    for biosample, settings in biosample_manifest_by_sample.items():
        sample_runs.setdefault(biosample, [])
        for run_id in _manifest_values(settings.get("run_ids", "")):
            run_id = _safe_manifest_id(
                run_id, "run_id", biosample_manifest or Path("<biosample-manifest>"),
                int(settings.get("_manifest_line", "0")),
            )
            if run_id not in sample_runs[biosample]:
                sample_runs[biosample].append(run_id)
    if mode != "existing-results":
        inferred_samples = infer_samples_from_tagalign(pooled_dir)
        for biosample, runs in inferred_samples.items():
            sample_runs.setdefault(biosample, runs)

    if not sample_runs:
        raise RuntimeError(
            f"No biosamples found from {map_file} or {pooled_dir}/*.tagAlign.gz"
        )

    frip_tool_error = ""
    if not skip_frip:
        try:
            ensure_frip_tools()
        except RuntimeError as exc:
            # A historical baseline should still report its recoverable MAPQ,
            # fragment-size and TSS evidence when bedtools is unavailable.
            # Core/release rows carry FRiP as unavailable and the release gate
            # will remain incomplete rather than aborting the entire project.
            frip_tool_error = str(exc)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=QC_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        for biosample, run_ids in sample_runs.items():
            replicate_groups = grouped_replicates.get(biosample, {})
            manifest_settings = biosample_manifest_by_sample.get(biosample, {})
            manifest_replicates: list[str] = []
            for value in _manifest_values(manifest_settings.get("biological_replicates", "")):
                normalized = value if value.lower().startswith("rep") else f"rep{value}"
                manifest_replicates.append(
                    _safe_manifest_id(
                        normalized,
                        "biological_replicates",
                        biosample_manifest or Path("<biosample-manifest>"),
                        int(manifest_settings.get("_manifest_line", "0")),
                    )
                )
            replicate_ids = list(replicate_groups) or list(dict.fromkeys(manifest_replicates))
            legacy_map = bool(replicate_groups) and any(
                record.legacy_two_column
                for group in replicate_groups.values()
                for record in group
            )
            if not replicate_groups and not manifest_replicates:
                legacy_map = bool(run_ids) or mode == "existing-results"
            replicate_dir_presence = [
                (align_dir / "replicates" / biosample / replicate).exists()
                for replicate in replicate_ids
            ]
            # Historical baselines must not auto-discover replicate directories:
            # those paths can belong to a later rerun.  They may still provide
            # explicitly declared run/fragment paths through their manifests.
            # For core/release, if any new-style output exists, require every
            # declared replicate; this prevents a partial set from passing.
            use_replicate_ids = (
                replicate_ids
                if mode != "existing-results" and any(replicate_dir_presence)
                else []
            )
            run_manifest_records = {
                row["run_id"]: row
                for row in run_manifest_by_sample.get(biosample, [])
            }
            explicit_layout = manifest_settings.get("layout", "AUTO")
            declared_layout = explicit_layout if explicit_layout != "AUTO" else (
                next(
                    iter(
                        {
                            record.layout
                            for group in replicate_groups.values()
                            for record in group
                            if record.layout != "AUTO"
                        }
                    )
                )
                if len(
                    {
                        record.layout
                        for group in replicate_groups.values()
                        for record in group
                        if record.layout != "AUTO"
                    }
                ) == 1
                else ""
            )
            run_qc = collect_run_qc(
                project_dir, align_dir, run_ids, mito_contigs,
                replicate_ids=use_replicate_ids, biosample=biosample,
                declared_layout=declared_layout,
                run_manifest_records=run_manifest_records,
                qc_mode=mode,
                artifact_dir=artifact_dir,
            )
            peak_frip_qc = collect_peak_and_frip_qc(
                biosample,
                pooled_dir,
                peaks_dir,
                artifact_dir,
                skip_frip,
                blacklist_bed,
                tss_bed,
                tagalign_override=_path_or_none(manifest_settings.get("tagalign_file")),
                fragment_override=_path_or_none(manifest_settings.get("fragment_file")),
                peak_override=_path_or_none(manifest_settings.get("peak_file")),
                qc_mode=mode,
                peak_provenance=manifest_settings.get("peak_provenance", ""),
                blacklist_peak_filter_status=manifest_settings.get(
                    "blacklist_peak_filter_status", ""
                ),
                tagalign_tn5_shifted=manifest_settings.get("tagalign_tn5_shifted", ""),
                frip_tool_error=frip_tool_error,
            )
            manifest_fragment_bams = [
                Path(value) for value in _manifest_values(manifest_settings.get("fragment_bams", ""))
            ]
            if manifest_fragment_bams:
                fragment_bams = manifest_fragment_bams
            elif use_replicate_ids:
                fragment_bams = [
                    align_dir / "replicates" / biosample / replicate /
                    f"{biosample}.{replicate}.nodup.bam"
                    for replicate in use_replicate_ids
                ]
            else:
                fragment_bams = [
                    _manifest_file(
                        run_manifest_records.get(run_id),
                        "final_bam",
                        align_dir / run_id / "bowtie2" / f"{run_id}.nodup.bam",
                    )
                    for run_id in run_ids
                ]
            fragment_qc = summarize_fragment_sizes(
                fragment_bams,
                artifact_dir / "fragment_size" / f"{biosample}.fragment_size.tsv",
                str(run_qc["layout"]),
            )
            reproducibility_qc = collect_reproducibility_qc(
                project_dir,
                biosample,
                reproducibility_override=_path_or_none(manifest_settings.get("reproducibility_json")),
                fingerprint_override=_path_or_none(manifest_settings.get("fingerprint_json")),
            )
            writer.writerow(
                format_summary_row(
                    project_dir,
                    biosample,
                    run_ids,
                    run_qc,
                    peak_frip_qc,
                    fragment_qc,
                    replicate_ids,
                    legacy_map,
                    threshold_profile,
                    reproducibility_qc,
                    mode,
                    blacklist_config_status,
                    blacklist_config_bed,
                    blacklist_config_sha256,
                    reference_config_fai,
                    reference_config_fai_sha256,
                )
            )


def main() -> int:
    args = parse_args()
    try:
        ref_config = load_ref_config(args.ref_name, args.project_name)
        project_dir = project_dir_from_args(args, ref_config)
        if args.mode == "release" and not ref_config.get("BLACKLIST_STATUS"):
            raise RuntimeError(
                "release QC config did not provide BLACKLIST_STATUS; verify epi_config.sh"
            )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    if not project_dir.exists():
        print(f"Error: project directory does not exist: {project_dir}", file=sys.stderr)
        return 1

    default_summary_name = (
        "atac_qc_summary.tsv"
        if args.mode == "release"
        else f"atac_qc_summary_{args.mode}.tsv"
    )
    output_file = args.output or (project_dir / "1_result" / "5_qc" / default_summary_name)
    mito_contigs = resolve_mito_contigs(args, ref_config)
    blacklist_bed = args.blacklist_bed or _path_or_none(ref_config.get("BLACKLIST_BED"))
    tss_bed = args.tss_bed or _path_or_none(ref_config.get("TSS_BED"))
    threshold_profile = args.threshold_profile or ref_config.get(
        "QC_THRESHOLD_PROFILE", "custom_nonhuman_phase23"
    )
    blacklist_config_bed = ""
    blacklist_config_sha256 = ""
    if ref_config.get("BLACKLIST_STATUS", "").upper() == "CONFIGURED" and blacklist_bed:
        blacklist_config_bed = str(blacklist_bed.expanduser().resolve())
        if blacklist_bed.is_file():
            blacklist_config_sha256 = sha256_file(blacklist_bed)
    reference_config_fai = ""
    reference_config_fai_sha256 = ""
    if ref_config.get("REF_GENOME"):
        configured_fai = Path(ref_config["REF_GENOME"] + ".fai")
        reference_config_fai = str(configured_fai.expanduser().resolve())
        if configured_fai.is_file():
            reference_config_fai_sha256 = sha256_file(configured_fai)

    try:
        summarize_project(
            project_dir=project_dir,
            output_file=output_file,
            mito_contigs=mito_contigs,
            skip_frip=args.skip_frip,
            blacklist_bed=blacklist_bed,
            tss_bed=tss_bed,
            threshold_profile=threshold_profile,
            mode=args.mode,
            run_manifest=args.run_manifest,
            biosample_manifest=args.biosample_manifest,
            blacklist_config_status=ref_config.get("BLACKLIST_STATUS", ""),
            blacklist_config_bed=blacklist_config_bed,
            blacklist_config_sha256=blacklist_config_sha256,
            reference_config_fai=reference_config_fai,
            reference_config_fai_sha256=reference_config_fai_sha256,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote QC summary: {output_file}")
    if args.ref_name:
        print(f"Used REF_NAME config: {args.ref_name}")
    if blacklist_bed:
        print(f"Blacklist BED: {blacklist_bed}")
    if tss_bed:
        print(f"TSS BED: {tss_bed}")
    print(f"QC mode: {args.mode}")
    if args.run_manifest:
        print(f"Run manifest: {args.run_manifest}")
    if args.biosample_manifest:
        print(f"Biosample manifest: {args.biosample_manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
