"""fastp JSON parsing helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from atac_qc.common import as_int


def read_fastp_metrics(fastp_json: Path) -> dict[str, Optional[float | int]]:
    """Return core pre/post-filter metrics without discarding fastp detail."""
    keys = (
        "raw_reads", "clean_reads", "raw_bases", "clean_bases",
        "raw_q20_rate", "clean_q20_rate", "raw_q30_rate", "clean_q30_rate",
        "raw_gc_content", "clean_gc_content",
    )
    empty: dict[str, Optional[float | int]] = {key: None for key in keys}
    if not fastp_json.exists():
        return empty
    with fastp_json.open("r", encoding="utf-8") as handle:
        summary = json.load(handle).get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})
    return {
        "raw_reads": as_int(before.get("total_reads")),
        "clean_reads": as_int(after.get("total_reads")),
        "raw_bases": as_int(before.get("total_bases")),
        "clean_bases": as_int(after.get("total_bases")),
        "raw_q20_rate": before.get("q20_rate"),
        "clean_q20_rate": after.get("q20_rate"),
        "raw_q30_rate": before.get("q30_rate"),
        "clean_q30_rate": after.get("q30_rate"),
        "raw_gc_content": before.get("gc_content"),
        "clean_gc_content": after.get("gc_content"),
    }


def read_fastp_reads(fastp_json: Path) -> tuple[Optional[int], Optional[int]]:
    metrics = read_fastp_metrics(fastp_json)
    return as_int(metrics["raw_reads"]), as_int(metrics["clean_reads"])
