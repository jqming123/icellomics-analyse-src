"""bedtools tagAlign FRiP helpers."""

from __future__ import annotations

import gzip
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Optional


FRIP_METHOD = "bedtools_fragment_any_overlap"
LEGACY_FRIP_METHOD = "bedtools_tagalign_record_any_overlap_legacy"


def interval_contigs(path: Path, max_records: int | None = None) -> set[str]:
    """Return observed BED/tagAlign contigs, using a bounded input prefix.

    A global ``chr1`` versus ``1`` mismatch must not be interpreted as a FRiP
    of zero.  The peak file is small enough to inspect in full; a prefix of a
    large interval input keeps this guard inexpensive before the normal
    full-file bedtools calculation.
    """

    opener = gzip.open if path.suffix == ".gz" else open
    contigs: set[str] = set()
    records = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t", 1)
            if not fields or not fields[0]:
                continue
            contigs.add(fields[0])
            records += 1
            if max_records is not None and records >= max_records:
                break
    return contigs


def validate_frip_contigs(interval_file: Path, peak_file: Path) -> set[str]:
    """Return shared contigs or raise when a naming mismatch is apparent."""

    input_contigs = interval_contigs(interval_file, max_records=100000)
    peak_contigs = interval_contigs(peak_file)
    if not input_contigs:
        raise ValueError(f"FRiP interval input has no BED records: {interval_file}")
    if not peak_contigs:
        raise ValueError(f"FRiP peak input has no BED records: {peak_file}")
    shared = input_contigs & peak_contigs
    if not shared:
        raise ValueError(
            "FRiP interval/peak contigs do not overlap; verify assembly and "
            f"chromosome naming. interval_prefix={sorted(input_contigs)[:5]}, "
            f"peaks={sorted(peak_contigs)[:5]}"
        )
    return shared


def ensure_frip_tools() -> None:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools is required for FRiP calculation but was not found in PATH")
    if shutil.which("zcat") is None:
        raise RuntimeError("zcat is required for FRiP calculation but was not found in PATH")


def count_all_records(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    count = 0
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as handle:
        for _ in handle:
            count += 1
    return count


def compute_reads_in_peaks(tagalign_file: Path, peak_file: Path) -> int:
    cmd = (
        "set -o pipefail; "
        f"zcat -f {shlex.quote(str(tagalign_file))} | "
        f"bedtools intersect -a stdin -b {shlex.quote(str(peak_file))} -wa -u | "
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
        message = result.stderr.strip() or "bedtools intersect failed"
        raise RuntimeError(message)
    return int(result.stdout.strip())


def compute_records_in_peaks(interval_file: Path, peak_file: Path) -> int:
    """Count input interval records overlapping at least one peak."""
    return compute_reads_in_peaks(interval_file, peak_file)
