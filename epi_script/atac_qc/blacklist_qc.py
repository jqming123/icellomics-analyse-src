"""Blacklist-overlap QC helpers."""

from __future__ import annotations

import gzip
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Optional


def _contigs_from_intervals(path: Path, max_records: Optional[int] = None) -> set[str]:
    """Read contig names from a BED/tagAlign-like text or gzip file.

    Blacklists are small enough to read completely.  TagAlign can be very
    large, so a bounded prefix is sufficient to catch the dangerous global
    naming mismatch (for example ``chr1`` versus ``1``) before an overlap of
    zero is reported as a biological result.
    """

    opener = gzip.open if path.suffix == ".gz" else open
    contigs: set[str] = set()
    records = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if not fields or not fields[0]:
                continue
            contigs.add(fields[0])
            records += 1
            if max_records is not None and records >= max_records:
                break
    return contigs


def _format_contigs(contigs: set[str], limit: int = 10) -> str:
    shown = sorted(contigs)[:limit]
    suffix = ",..." if len(contigs) > limit else ""
    return ",".join(shown) + suffix


def compute_blacklist_qc(
    tagalign_file: Path,
    blacklist_bed: Optional[Path],
    tagalign_total_reads: Optional[int],
) -> dict[str, object]:
    """Measure tagAlign overlap with a blacklist without claiming peak filtering.

    This is intentionally a *post-hoc tagAlign overlap* measurement.  It does
    not establish that an existing peak file was blacklist-filtered; callers
    must record peak provenance separately.
    """

    scope = "posthoc_tagalign_overlap"
    if blacklist_bed is None:
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "NOT_CONFIGURED",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": "",
            "warnings": ["blacklist_not_configured"],
        }
    if not blacklist_bed.exists():
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "MISSING",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": "",
            "warnings": [f"blacklist_missing:{blacklist_bed}"],
        }
    if not tagalign_file.exists():
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "MISSING_TAGALIGN",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": "",
            "warnings": [f"blacklist_missing_tagalign:{tagalign_file}"],
        }
    if not tagalign_total_reads:
        return {
            "blacklist_reads": 0,
            "blacklist_fraction": None,
            "blacklist_qc_status": "NO_TAGALIGN_RECORDS",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": "",
            "warnings": ["blacklist_no_tagalign_records"],
        }
    try:
        blacklist_contigs = _contigs_from_intervals(blacklist_bed)
        tagalign_contigs = _contigs_from_intervals(tagalign_file, max_records=100000)
    except (OSError, EOFError, ValueError) as exc:
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "FAILED",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": "",
            "warnings": [f"blacklist_contig_check_failed:{exc}"],
        }
    shared_contigs = blacklist_contigs & tagalign_contigs
    if not shared_contigs:
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "CONTIG_MISMATCH",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": "",
            "warnings": [
                "blacklist_contig_mismatch:"
                f"blacklist={_format_contigs(blacklist_contigs)};"
                f"tagalign_prefix={_format_contigs(tagalign_contigs)}"
            ],
        }
    if shutil.which("bedtools") is None or shutil.which("zcat") is None:
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "TOOLS_NOT_FOUND",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": ",".join(sorted(shared_contigs)),
            "warnings": ["blacklist_tools_not_found"],
        }

    cmd = (
        "set -o pipefail; "
        f"zcat -f {shlex.quote(str(tagalign_file))} | "
        f"bedtools intersect -a stdin -b {shlex.quote(str(blacklist_bed))} -wa -u | "
        "wc -l"
    )
    result = subprocess.run(
        cmd,
        shell=True,
        executable="/bin/bash",
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        message = result.stderr.strip() or "bedtools blacklist intersect failed"
        return {
            "blacklist_reads": None,
            "blacklist_fraction": None,
            "blacklist_qc_status": "FAILED",
            "blacklist_qc_scope": scope,
            "blacklist_contig_intersection": ",".join(sorted(shared_contigs)),
            "warnings": [f"blacklist_failed:{message}"],
        }

    blacklist_reads = int(result.stdout.strip())
    return {
        "blacklist_reads": blacklist_reads,
        "blacklist_fraction": blacklist_reads / tagalign_total_reads,
        "blacklist_qc_status": "COMPUTED",
        "blacklist_qc_scope": scope,
        "blacklist_contig_intersection": ",".join(sorted(shared_contigs)),
        "warnings": [],
    }
