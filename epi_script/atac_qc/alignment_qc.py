"""Alignment, duplicate, and mitochondrial QC parsers."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path
from typing import Mapping, Optional, Sequence

from atac_qc.common import as_int


def flagstat_count(counts: Mapping[str, int], label: str) -> Optional[int]:
    """Return a flagstat count, tolerating annotated ``-O tsv`` labels.

    ``samtools flagstat -O tsv`` labels the leading row
    ``total (QC-passed reads + QC-failed reads)``, while the default text output
    labels it ``in total``.  Matching the ``total`` row by prefix keeps both
    layouts readable instead of silently dropping the denominator.
    """

    if label in counts:
        return counts[label]
    if label == "total":
        if "in total" in counts:
            return counts["in total"]
        for key, value in counts.items():
            if key.startswith("total"):
                return value
    return None


def parse_flagstat_output(text: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        tsv_fields = line.split("\t")
        if len(tsv_fields) >= 3:
            pass_count = as_int(tsv_fields[0])
            if pass_count is not None:
                counts[tsv_fields[2].strip()] = pass_count
                continue

        match = re.match(r"^(\d+)\s+\+\s+\d+\s+(.+?)(?:\s+\(|$)", line)
        if match:
            counts[match.group(2).strip()] = int(match.group(1))
    return counts


def read_or_run_flagstat(bam_file: Path) -> tuple[dict[str, int], list[str]]:
    warnings: list[str] = []
    sidecar = Path(f"{bam_file}.flagstat")
    if sidecar.exists():
        return parse_flagstat_output(sidecar.read_text(encoding="utf-8", errors="replace")), warnings
    if not bam_file.exists():
        return {}, [f"missing_raw_bam_and_flagstat:{bam_file}"]
    if shutil.which("samtools") is None:
        return {}, ["samtools_not_found_for_flagstat"]

    result = subprocess.run(
        ["samtools", "flagstat", "-O", "tsv", str(bam_file)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        message = result.stderr.strip() or "samtools flagstat failed"
        warnings.append(f"flagstat_failed:{bam_file}:{message}")
        return {}, warnings
    return parse_flagstat_output(result.stdout), warnings


def summarize_raw_bam_flagstat(bam_file: Path) -> dict[str, object]:
    counts, warnings = read_or_run_flagstat(bam_file)
    total_reads = flagstat_count(counts, "total")
    mapped_reads = flagstat_count(counts, "mapped")
    paired_reads = counts.get("paired in sequencing")
    proper_pair_reads = counts.get("properly paired")
    return {
        "raw_bam_total_reads": total_reads,
        "raw_bam_mapped_reads": mapped_reads,
        "raw_bam_paired_reads": paired_reads,
        "raw_bam_proper_pair_reads": proper_pair_reads,
        "warnings": warnings,
    }


def read_idxstats(
    idxstats_file: Path,
    mito_contigs: Sequence[str],
) -> tuple[Optional[int], Optional[int]]:
    if not idxstats_file.exists():
        return None, None

    total = 0
    mito = 0
    mito_set = set(mito_contigs)
    with idxstats_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            contig = fields[0]
            mapped = as_int(fields[2]) or 0
            total += mapped
            if contig in mito_set:
                mito += mapped
    return total, mito


def read_picard_dup_metrics(metrics_file: Path) -> tuple[Optional[int], Optional[int]]:
    metrics = read_picard_dup_details(metrics_file)
    return metrics["examined_reads"], metrics["duplicate_reads"]


def read_picard_dup_details(metrics_file: Path) -> dict[str, Optional[int]]:
    empty = {"examined_reads": None, "duplicate_reads": None, "optical_duplicate_reads": None}
    if not metrics_file.exists():
        return empty

    header: Optional[list[str]] = None
    with metrics_file.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("##"):
                continue
            fields = line.split("\t")
            if "UNPAIRED_READS_EXAMINED" in fields and "READ_PAIRS_EXAMINED" in fields:
                header = fields
                continue
            if header is None:
                continue
            values = dict(zip(header, fields))
            unpaired_examined = as_int(values.get("UNPAIRED_READS_EXAMINED")) or 0
            read_pairs_examined = as_int(values.get("READ_PAIRS_EXAMINED")) or 0
            unpaired_duplicates = as_int(values.get("UNPAIRED_READ_DUPLICATES")) or 0
            read_pair_duplicates = as_int(values.get("READ_PAIR_DUPLICATES")) or 0
            examined_reads = unpaired_examined + 2 * read_pairs_examined
            duplicate_reads = unpaired_duplicates + 2 * read_pair_duplicates
            optical_pairs = as_int(values.get("READ_PAIR_OPTICAL_DUPLICATES"))
            return {
                "examined_reads": examined_reads,
                "duplicate_reads": duplicate_reads,
                "optical_duplicate_reads": None if optical_pairs is None else 2 * optical_pairs,
            }
    return empty
