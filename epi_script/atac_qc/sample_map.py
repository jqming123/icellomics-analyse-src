"""Sample/run discovery and sequencing-layout helpers.

The preferred sample map has three or four tab-separated columns::

    biosample_id    biological_replicate    run_id    layout

``layout`` is optional and may be PE, SE, or AUTO.  The legacy two-column
format remains readable, but all runs are assigned to biological replicate
``rep1`` and the record is marked as legacy so QC cannot silently claim that
replicate-aware processing was performed.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SAFE_ID = re.compile(r"^[A-Za-z0-9._-]+$")
VALID_LAYOUTS = {"PE", "SE", "AUTO"}
HEADER_NAMES = {"biosample", "biosample_id", "sample", "sample_id"}


@dataclass(frozen=True)
class ReplicateRecord:
    biosample_id: str
    biological_replicate: str
    run_id: str
    layout: str = "AUTO"
    legacy_two_column: bool = False


def _normalise_replicate(value: str) -> str:
    value = value.strip()
    if not value:
        raise ValueError("biological_replicate cannot be empty")
    return value if value.lower().startswith("rep") else f"rep{value}"


def read_replicate_records(map_file: Path) -> list[ReplicateRecord]:
    """Read and validate a sample map while preserving input order."""
    records: list[ReplicateRecord] = []
    if not map_file.exists():
        return records

    seen_runs: dict[str, tuple[str, str]] = {}
    with map_file.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            row = [value.strip() for value in row]
            if not row or not row[0] or row[0].startswith("#"):
                continue
            if row[0].lower() in HEADER_NAMES:
                continue
            if len(row) < 2:
                raise ValueError(f"{map_file}:{line_number}: expected at least 2 columns")

            legacy = len(row) == 2
            if legacy:
                biosample, run_id = row[:2]
                replicate, layout = "rep1", "AUTO"
            else:
                biosample, replicate, run_id = row[:3]
                layout = row[3].upper() if len(row) >= 4 and row[3] else "AUTO"
                replicate = _normalise_replicate(replicate)

            for label, value in (("biosample_id", biosample), ("biological_replicate", replicate), ("run_id", run_id)):
                if not value or not SAFE_ID.fullmatch(value):
                    raise ValueError(f"{map_file}:{line_number}: unsafe or empty {label}: {value!r}")
            if layout not in VALID_LAYOUTS:
                raise ValueError(f"{map_file}:{line_number}: layout must be PE, SE, or AUTO; got {layout!r}")

            group = (biosample, replicate)
            if run_id in seen_runs and seen_runs[run_id] != group:
                previous = "/".join(seen_runs[run_id])
                raise ValueError(f"{map_file}:{line_number}: run {run_id!r} is already assigned to {previous}")
            seen_runs[run_id] = group
            record = ReplicateRecord(biosample, replicate, run_id, layout, legacy)
            if record not in records:
                records.append(record)
    return records


def group_replicates(records: Iterable[ReplicateRecord]) -> "OrderedDict[str, OrderedDict[str, list[ReplicateRecord]]]":
    grouped: "OrderedDict[str, OrderedDict[str, list[ReplicateRecord]]]" = OrderedDict()
    for record in records:
        grouped.setdefault(record.biosample_id, OrderedDict())
        grouped[record.biosample_id].setdefault(record.biological_replicate, [])
        grouped[record.biosample_id][record.biological_replicate].append(record)
    return grouped


def read_sample_run_map(map_file: Path) -> "OrderedDict[str, list[str]]":
    """Compatibility view used by older callers: biosample -> run IDs."""
    sample_runs: "OrderedDict[str, list[str]]" = OrderedDict()
    for record in read_replicate_records(map_file):
        sample_runs.setdefault(record.biosample_id, [])
        if record.run_id not in sample_runs[record.biosample_id]:
            sample_runs[record.biosample_id].append(record.run_id)
    return sample_runs


def infer_samples_from_tagalign(pooled_dir: Path) -> "OrderedDict[str, list[str]]":
    samples: "OrderedDict[str, list[str]]" = OrderedDict()
    for tagalign in sorted(pooled_dir.glob("*.tagAlign.gz")):
        biosample = tagalign.name.split(".", 1)[0]
        samples.setdefault(biosample, [])
    return samples


def detect_layout(project_dir: Path, run_id: str) -> str:
    reads_dir = project_dir / "1_result" / "0_fastq" / run_id / "reads"
    raw_r1 = reads_dir / f"{run_id}_1.fastq.gz"
    raw_r2 = reads_dir / f"{run_id}_2.fastq.gz"
    raw_se = reads_dir / f"{run_id}.fastq.gz"
    if raw_r1.exists() and raw_r2.exists():
        return "PE"
    if raw_r1.exists() or raw_r2.exists():
        return "incomplete_PE"
    if raw_se.exists():
        return "SE"
    return "unknown"


def summarize_layout(layouts: Sequence[str]) -> str:
    clean = [item for item in layouts if item and item != "unknown"]
    if not clean:
        return "unknown"
    unique = sorted(set(clean))
    if len(unique) == 1:
        return unique[0]
    return "mixed:" + ",".join(unique)


def _resolved_layout(records: Sequence[ReplicateRecord], project_dir: Path | None) -> str:
    explicit = {record.layout for record in records if record.layout != "AUTO"}
    detected = {detect_layout(project_dir, record.run_id) for record in records} if project_dir else set()
    layouts = explicit or {item for item in detected if item != "unknown"}
    if len(layouts) != 1:
        return summarize_layout(sorted(layouts))
    return next(iter(layouts))


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate/query an ATAC sample-run map")
    parser.add_argument("--map", required=True, type=Path)
    parser.add_argument("--project-dir", type=Path)
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--list-biosamples", action="store_true")
    parser.add_argument("--biosample")
    parser.add_argument("--emit-replicates", action="store_true")
    args = parser.parse_args()

    records = read_replicate_records(args.map)
    if not records:
        raise SystemExit(f"No records found in {args.map}")
    grouped = group_replicates(records)
    if args.list_biosamples:
        print("\n".join(grouped))
    if args.emit_replicates:
        if not args.biosample or args.biosample not in grouped:
            raise SystemExit("--emit-replicates requires a biosample present in the map")
        for replicate, rep_records in grouped[args.biosample].items():
            layout = _resolved_layout(rep_records, args.project_dir)
            legacy = "true" if any(record.legacy_two_column for record in rep_records) else "false"
            print("\t".join((replicate, layout, legacy, ",".join(record.run_id for record in rep_records))))
    if args.validate and any(record.legacy_two_column for record in records):
        print("WARNING: legacy two-column map detected; all runs are treated as biological replicate rep1", file=__import__("sys").stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
