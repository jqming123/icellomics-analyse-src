"""Peak-file QC helpers."""

from __future__ import annotations

import gzip
import statistics
from pathlib import Path
from typing import Optional


def count_bed_records(path: Path) -> Optional[int]:
    if not path.exists():
        return None
    count = 0
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip() and line[:1] != "#":
                count += 1
    return count


def summarize_peak_file(path: Path) -> dict[str, Optional[float] | Optional[int]]:
    if not path.exists():
        return {
            "peak_count": None,
            "peak_total_bp": None,
            "peak_median_width": None,
        }

    widths: list[int] = []
    opener = gzip.open if path.suffix == ".gz" else open
    mode = "rt" if path.suffix == ".gz" else "r"
    with opener(path, mode, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line[:1] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            width = max(end - start, 0)
            widths.append(width)

    if not widths:
        return {
            "peak_count": 0,
            "peak_total_bp": 0,
            "peak_median_width": None,
        }

    return {
        "peak_count": len(widths),
        "peak_total_bp": sum(widths),
        "peak_median_width": float(statistics.median(widths)),
    }
