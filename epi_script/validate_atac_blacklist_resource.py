#!/usr/bin/env python3
"""Validate an assembly-matched ATAC blacklist/exclusion BED once per resource.

This tool is intentionally independent of peak calling.  It verifies that a
resource configured for release has valid BED intervals, matches every contig
and coordinate bound in the exact FASTA index used by the pipeline, and has a
stable content checksum.  It does not liftover, rename, filter, or otherwise
alter its input.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO


class ResourceValidationError(ValueError):
    """Raised when a blacklist resource cannot be safely used with a FASTA."""


def _open_text(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8", errors="strict") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8", errors="strict"
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fai(path: Path) -> dict[str, int]:
    contigs: dict[str, int] = {}
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 2 or not fields[0]:
                raise ResourceValidationError(f"{path}:{line_number}: malformed FASTA index record")
            try:
                length = int(fields[1])
            except ValueError as exc:
                raise ResourceValidationError(
                    f"{path}:{line_number}: FASTA contig length is not an integer"
                ) from exc
            if length <= 0 or fields[0] in contigs:
                raise ResourceValidationError(f"{path}:{line_number}: invalid/duplicate FASTA contig {fields[0]!r}")
            contigs[fields[0]] = length
    if not contigs:
        raise ResourceValidationError(f"FASTA index is empty: {path}")
    return contigs


def validate_blacklist_resource(
    bed: Path,
    reference_fai: Path,
    *,
    expected_intervals: int | None = None,
    require_ensembl_grcm39_contigs: bool = False,
) -> dict[str, object]:
    """Validate BED syntax, sorting, bounds, names, and optional record count."""

    bed = bed.expanduser().resolve()
    reference_fai = reference_fai.expanduser().resolve()
    if not bed.is_file():
        raise ResourceValidationError(f"Blacklist BED does not exist: {bed}")
    if not reference_fai.is_file():
        raise ResourceValidationError(f"Reference FASTA index does not exist: {reference_fai}")

    reference_contigs = read_fai(reference_fai)
    record_count = 0
    intervals_by_contig: dict[str, int] = {}
    seen_contigs: set[str] = set()
    closed_contigs: set[str] = set()
    current_contig = ""
    previous_start = -1
    allowed_mouse = {str(number) for number in range(1, 20)} | {"X", "Y", "MT"}

    try:
        with _open_text(bed) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.rstrip("\r\n")
                if not stripped or stripped.startswith(("#", "track", "browser")):
                    continue
                fields = stripped.split()
                if len(fields) < 3:
                    raise ResourceValidationError(f"{bed}:{line_number}: BED record has fewer than 3 columns")
                contig = fields[0]
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise ResourceValidationError(
                        f"{bed}:{line_number}: BED start/end are not integers"
                    ) from exc
                if start < 0 or end <= start:
                    raise ResourceValidationError(f"{bed}:{line_number}: invalid BED interval {contig}:{start}-{end}")
                if contig not in reference_contigs:
                    raise ResourceValidationError(
                        f"{bed}:{line_number}: contig {contig!r} is absent from {reference_fai}"
                    )
                if end > reference_contigs[contig]:
                    raise ResourceValidationError(
                        f"{bed}:{line_number}: interval {contig}:{start}-{end} exceeds reference length "
                        f"{reference_contigs[contig]}"
                    )
                if require_ensembl_grcm39_contigs and contig not in allowed_mouse:
                    raise ResourceValidationError(
                        f"{bed}:{line_number}: {contig!r} is not an Ensembl GRCm39 primary contig"
                    )

                if contig != current_contig:
                    if current_contig:
                        closed_contigs.add(current_contig)
                    if contig in closed_contigs:
                        raise ResourceValidationError(
                            f"{bed}:{line_number}: contig {contig!r} reappears after another contig; BED is not grouped"
                        )
                    current_contig = contig
                    previous_start = -1
                if start < previous_start:
                    raise ResourceValidationError(
                        f"{bed}:{line_number}: start coordinate decreases within contig {contig!r}"
                    )
                previous_start = start
                seen_contigs.add(contig)
                intervals_by_contig[contig] = intervals_by_contig.get(contig, 0) + 1
                record_count += 1
    except (OSError, UnicodeError, gzip.BadGzipFile) as exc:
        raise ResourceValidationError(f"Unable to read blacklist BED {bed}: {exc}") from exc

    if record_count == 0:
        raise ResourceValidationError(f"Blacklist BED has no intervals: {bed}")
    if expected_intervals is not None and record_count != expected_intervals:
        raise ResourceValidationError(
            f"Expected {expected_intervals} blacklist intervals, observed {record_count}: {bed}"
        )

    return {
        "status": "VERIFIED",
        "blacklist_bed": str(bed),
        "blacklist_sha256": sha256_file(bed),
        "reference_fai": str(reference_fai),
        "reference_fai_sha256": sha256_file(reference_fai),
        "interval_count": record_count,
        "contigs": sorted(seen_contigs),
        "intervals_by_contig": intervals_by_contig,
        "expected_intervals": expected_intervals,
        "ensembl_grcm39_primary_contigs_required": require_ensembl_grcm39_contigs,
        "validated_utc": datetime.now(timezone.utc).isoformat(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bed", required=True, type=Path, help="Plain or gzip-compressed BED resource")
    parser.add_argument("--reference-fai", required=True, type=Path, help="Exact FASTA .fai used for alignment")
    parser.add_argument("--expected-intervals", type=int, help="Optional exact record count")
    parser.add_argument(
        "--require-ensembl-grcm39-contigs",
        action="store_true",
        help="Reject chr-prefixed/non-primary contigs; allow only 1..19, X, Y, MT",
    )
    parser.add_argument("--output", type=Path, help="Optional JSON validation record")
    args = parser.parse_args()

    try:
        result = validate_blacklist_resource(
            args.bed,
            args.reference_fai,
            expected_intervals=args.expected_intervals,
            require_ensembl_grcm39_contigs=args.require_ensembl_grcm39_contigs,
        )
    except ResourceValidationError as exc:
        parser.error(str(exc))

    payload = json.dumps(result, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
        print(f"Wrote validation record: {args.output}")
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
