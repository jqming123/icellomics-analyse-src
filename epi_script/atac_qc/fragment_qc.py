"""Fragment-size QC for paired-end ATAC-seq BAM files."""

from __future__ import annotations

import math
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import Sequence


def _unavailable_fragment_result(
    warnings: list[str], input_status: str = "UNAVAILABLE"
) -> dict[str, object]:
    """Return an unavailable (not failed) fragment-size result.

    An absent historical BAM inventory means there was no measurement.  It is
    materially different from reading a valid PE BAM and observing zero valid
    positive TLEN records, which is a real QC failure.
    """

    return {
        "fragment_count": None,
        "fragment_mean_size": None,
        "fragment_median_size": None,
        "fragment_p10_size": None,
        "fragment_p90_size": None,
        "nfr_fraction": None,
        "mono_nucleosome_fraction": None,
        "di_nucleosome_fraction": None,
        "fragment_periodicity_score": None,
        "fragment_histogram_file": "",
        "fragment_input_status": input_status,
        "warnings": warnings,
    }


def fragment_periodicity_metrics(histogram: Counter[int], total: int) -> dict[str, float | None]:
    if total <= 0:
        return {"nfr_fraction": None, "mono_nucleosome_fraction": None,
                "di_nucleosome_fraction": None, "fragment_periodicity_score": None}
    nfr = sum(count for size, count in histogram.items() if size <= 100)
    mono = sum(count for size, count in histogram.items() if 180 <= size <= 247)
    di = sum(count for size, count in histogram.items() if 315 <= size <= 473)
    trough = sum(count for size, count in histogram.items() if 101 <= size <= 179 or 248 <= size <= 314)
    return {
        "nfr_fraction": nfr / total,
        "mono_nucleosome_fraction": mono / total,
        "di_nucleosome_fraction": di / total,
        "fragment_periodicity_score": (mono + di) / trough if trough else None,
    }


def _percentile_from_histogram(histogram: Counter[int], total: int, percentile: float) -> int | None:
    if total <= 0:
        return None
    target = max(1, math.ceil(total * percentile))
    seen = 0
    for size in sorted(histogram):
        seen += histogram[size]
        if seen >= target:
            return size
    return None


def summarize_fragment_sizes(
    bam_files: Sequence[Path],
    output_file: Path,
    layout: str,
    max_fragment_size: int = 2000,
) -> dict[str, object]:
    if "PE" not in layout:
        return _unavailable_fragment_result(
            ["fragment_size_requires_pe"], "NOT_APPLICABLE"
        )
    if not bam_files:
        return _unavailable_fragment_result(["fragment_size_bam_not_available"])
    if shutil.which("samtools") is None:
        return _unavailable_fragment_result(["samtools_not_found_for_fragment_size"])

    warnings: list[str] = []
    histogram: Counter[int] = Counter()
    successfully_read_bams = 0
    for bam_file in bam_files:
        if not bam_file.exists():
            warnings.append(f"missing_fragment_bam:{bam_file}")
            continue

        proc = subprocess.Popen(
            ["samtools", "view", str(bam_file)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if proc.stdout is None or proc.stderr is None:
            warnings.append(f"fragment_size_failed:{bam_file}:failed_to_open_samtools_stream")
            continue

        for line in proc.stdout:
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            try:
                template_len = int(fields[8])
            except ValueError:
                continue
            if template_len <= 0:
                continue
            if template_len > max_fragment_size:
                continue
            histogram[template_len] += 1

        stderr = proc.stderr.read().strip()
        return_code = proc.wait()
        if return_code != 0:
            message = stderr or "samtools view failed"
            warnings.append(f"fragment_size_failed:{bam_file}:{message}")
        else:
            successfully_read_bams += 1

    total = sum(histogram.values())
    input_status = "COMPLETE" if successfully_read_bams == len(bam_files) else "PARTIAL"
    if total == 0:
        if successfully_read_bams == 0:
            warnings.append("fragment_size_bam_not_available")
            return _unavailable_fragment_result(warnings)
        warnings.append("fragment_size_no_records")
        return {
            "fragment_count": 0,
            "fragment_mean_size": None,
            "fragment_median_size": None,
            "fragment_p10_size": None,
            "fragment_p90_size": None,
            "nfr_fraction": None,
            "mono_nucleosome_fraction": None,
            "di_nucleosome_fraction": None,
            "fragment_periodicity_score": None,
            "fragment_histogram_file": "",
            "fragment_input_status": input_status,
            "warnings": warnings,
        }

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        handle.write("fragment_size\tcount\n")
        for size in sorted(histogram):
            handle.write(f"{size}\t{histogram[size]}\n")

    mean_size = sum(size * count for size, count in histogram.items()) / total
    periodicity = fragment_periodicity_metrics(histogram, total)
    return {
        "fragment_count": total,
        "fragment_mean_size": mean_size,
        "fragment_median_size": _percentile_from_histogram(histogram, total, 0.5),
        "fragment_p10_size": _percentile_from_histogram(histogram, total, 0.1),
        "fragment_p90_size": _percentile_from_histogram(histogram, total, 0.9),
        **periodicity,
        "fragment_histogram_file": str(output_file),
        "fragment_input_status": input_status,
        "warnings": warnings,
    }
