#!/usr/bin/env python3
"""Merge per-project JBrowse ATAC track manifests into one cell-line manifest.

The pipeline step 5 writes one manifest per project, for example:

  PRJNA915236_HEK293/1_result/4_jbrowse_manifest/
    atac_jbrowse_track_manifest.PRJNA915236_HEK293.tsv

This script scans those project manifests and creates one merged TSV per cell
line. The merged file keeps the same columns as the original manifests so it
can be used directly by a later DS0005 browser-track import step.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


EXPECTED_COLUMNS = [
    "dataset",
    "project_id",
    "cell_line",
    "sample_id",
    "species_id",
    "assembly_name",
    "ref_name",
    "signal_bw_path",
    "peaks_path",
    "peaks_index_path",
    "summits_path",
    "summits_index_path",
    "track_meta_path",
    "assembly_meta_path",
]

REQUIRED_PATH_COLUMNS = [
    "signal_bw_path",
    "peaks_path",
    "peaks_index_path",
    "track_meta_path",
    "assembly_meta_path",
]

OPTIONAL_PATH_COLUMNS = [
    "summits_path",
    "summits_index_path",
]


@dataclass(frozen=True)
class ManifestRow:
    values: dict[str, str]
    source_manifest: Path

    @property
    def dataset(self) -> str:
        return self.values["dataset"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge per-project JBrowse ATAC manifests by cell line."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Root containing project directories such as PRJNA915236_HEK293.",
    )
    parser.add_argument(
        "--cell-line",
        required=True,
        help="Cell line name to merge, for example HEK293.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output merged manifest TSV path.",
    )
    parser.add_argument(
        "--jbrowse-data-root",
        type=Path,
        help=(
            "Optional JBrowse data root used to validate relative paths, "
            "for example /hpcdisk1/.../jbrowse2_data."
        ),
    )
    parser.add_argument(
        "--manifest-glob",
        default="*/1_result/4_jbrowse_manifest/atac_jbrowse_track_manifest.*.tsv",
        help="Glob relative to input-root. Default scans pipeline step 5 output.",
    )
    parser.add_argument(
        "--allow-conflicting-duplicates",
        action="store_true",
        help=(
            "Keep the first row when the same dataset appears with different "
            "metadata. By default conflicting duplicates stop the merge."
        ),
    )
    parser.add_argument(
        "--allow-missing-files",
        action="store_true",
        help=(
            "Do not fail if --jbrowse-data-root validation finds missing files. "
            "Missing paths are still reported to stderr."
        ),
    )
    return parser.parse_args()


def discover_manifests(input_root: Path, manifest_glob: str, cell_line: str) -> list[Path]:
    candidates = sorted(input_root.glob(manifest_glob))
    suffix = f"_{cell_line}"
    manifests = []
    for path in candidates:
        project_name = path.stem.replace("atac_jbrowse_track_manifest.", "", 1)
        if project_name.endswith(suffix):
            manifests.append(path)
    return manifests


def read_manifest(path: Path, cell_line: str) -> list[ManifestRow]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise SystemExit(
                f"Unexpected header in {path}\n"
                f"Found:    {reader.fieldnames}\n"
                f"Expected: {EXPECTED_COLUMNS}"
            )

        rows = []
        for line_number, row in enumerate(reader, start=2):
            cleaned = {column: (row.get(column) or "").strip() for column in EXPECTED_COLUMNS}
            if not cleaned["dataset"]:
                raise SystemExit(f"Missing dataset in {path}:{line_number}")
            if cleaned["cell_line"] != cell_line:
                raise SystemExit(
                    f"Unexpected cell_line in {path}:{line_number}: "
                    f"{cleaned['cell_line']} != {cell_line}"
                )
            rows.append(ManifestRow(values=cleaned, source_manifest=path))
    return rows


def same_manifest_values(left: ManifestRow, right: ManifestRow) -> bool:
    return all(left.values[column] == right.values[column] for column in EXPECTED_COLUMNS)


def merge_rows(
    rows: Iterable[ManifestRow],
    *,
    allow_conflicting_duplicates: bool,
) -> list[ManifestRow]:
    merged: dict[str, ManifestRow] = {}
    duplicate_count = 0

    for row in rows:
        existing = merged.get(row.dataset)
        if existing is None:
            merged[row.dataset] = row
            continue

        duplicate_count += 1
        if same_manifest_values(existing, row):
            continue

        message = (
            f"Conflicting duplicate dataset: {row.dataset}\n"
            f"First source: {existing.source_manifest}\n"
            f"Second source: {row.source_manifest}"
        )
        if allow_conflicting_duplicates:
            print(f"Warning: {message}", file=sys.stderr)
            continue
        raise SystemExit(message)

    sorted_rows = sorted(
        merged.values(),
        key=lambda item: (item.values["project_id"], item.values["sample_id"], item.dataset),
    )
    print(f"Duplicate dataset rows skipped: {duplicate_count}")
    return sorted_rows


def validate_paths(
    rows: Sequence[ManifestRow],
    jbrowse_data_root: Path,
    *,
    allow_missing_files: bool,
) -> None:
    missing: list[str] = []

    for row in rows:
        for column in REQUIRED_PATH_COLUMNS:
            relative_path = row.values[column]
            if relative_path and not (jbrowse_data_root / relative_path).is_file():
                missing.append(f"{row.dataset}\t{column}\t{relative_path}")

        for column in OPTIONAL_PATH_COLUMNS:
            relative_path = row.values[column]
            if relative_path and not (jbrowse_data_root / relative_path).is_file():
                missing.append(f"{row.dataset}\t{column}\t{relative_path}")

    if not missing:
        print("Path validation passed.")
        return

    print("Missing files found during path validation:", file=sys.stderr)
    for item in missing[:50]:
        print(f"  {item}", file=sys.stderr)
    if len(missing) > 50:
        print(f"  ... {len(missing) - 50} more missing paths", file=sys.stderr)

    if not allow_missing_files:
        raise SystemExit(
            "Path validation failed. Re-run with --allow-missing-files to write "
            "the merged manifest anyway."
        )


def write_manifest(rows: Sequence[ManifestRow], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=EXPECTED_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row.values)


def summarize(rows: Sequence[ManifestRow], manifests: Sequence[Path], output: Path) -> None:
    projects = sorted({row.values["project_id"] for row in rows})
    assemblies = sorted({f"{row.values['species_id']}:{row.values['assembly_name']}" for row in rows})

    print(f"Input manifests: {len(manifests)}")
    print(f"Merged datasets: {len(rows)}")
    print(f"Projects:        {len(projects)}")
    print(f"Assemblies:      {', '.join(assemblies)}")
    print(f"Output written:  {output}")


def main() -> None:
    args = parse_args()
    input_root = args.input_root.resolve()
    if not input_root.is_dir():
        raise SystemExit(f"Input root does not exist or is not a directory: {input_root}")

    manifests = discover_manifests(input_root, args.manifest_glob, args.cell_line)
    if not manifests:
        raise SystemExit(
            f"No manifests found for cell line {args.cell_line} under {input_root}"
        )

    all_rows: list[ManifestRow] = []
    for manifest in manifests:
        all_rows.extend(read_manifest(manifest, args.cell_line))

    merged_rows = merge_rows(
        all_rows,
        allow_conflicting_duplicates=args.allow_conflicting_duplicates,
    )

    if args.jbrowse_data_root:
        validate_paths(
            merged_rows,
            args.jbrowse_data_root.resolve(),
            allow_missing_files=args.allow_missing_files,
        )

    write_manifest(merged_rows, args.output.resolve())
    summarize(merged_rows, manifests, args.output.resolve())


if __name__ == "__main__":
    main()
