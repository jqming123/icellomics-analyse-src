#!/usr/bin/env python3
"""Prepare DS0005 import assets for ATAC peak analysis results.

This wrapper reuses merge_atac_peak_results.py to generate:
1. A merged TSV/CSV that matches tb_atac_peak
2. A manifest summary file
3. A ready-to-edit LOAD DATA SQL template
4. Validation SQL snippets

It is designed so DS0002 can remain unchanged while DS0005 is populated.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence


XLS_SUFFIX = "_peaks.xls"
NARROWPEAK_SUFFIX = "_peaks.narrowPeak"
DEFAULT_PROJECT_REGEX = r"^(PRJ[^_]+)_(.+)$"
OUTPUT_COLUMNS = [
    "dataset",
    "species_id",
    "cell_line",
    "sample_id",
    "project_id",
    "chr",
    "start",
    "end",
    "pileup",
    "p_value",
    "q_value",
    "fold_enrichment",
    "peak_name",
    "peak_region",
    "peak_source_file",
]
MANIFEST_COLUMNS = [
    "summary_level",
    "cell_line",
    "project_id",
    "sample_id",
    "dataset_count",
    "sample_count",
    "source_file_count",
    "peak_count",
]
CELL_LINE_META_COLUMNS = [
    "cell_line",
    "dataset_count",
    "peak_count",
    "sort_order",
    "updated_at",
]
DATASET_META_COLUMNS = [
    "cell_line",
    "dataset",
    "project_id",
    "sample_count",
    "peak_count",
    "updated_at",
]
EXPECTED_XLS_HEADER = [
    "chr",
    "start",
    "end",
    "length",
    "abs_summit",
    "pileup",
    "-log10(pvalue)",
    "fold_enrichment",
    "-log10(qvalue)",
    "name",
]
SUMMARY_GROUP_SPECS = [
    ("cell_line", lambda item: (item.cell_line,), ("cell_line",)),
    ("project", lambda item: (item.cell_line, item.project_id), ("cell_line", "project_id")),
    (
        "sample",
        lambda item: (item.cell_line, item.project_id, item.sample_id),
        ("cell_line", "project_id", "sample_id"),
    ),
]


@dataclass
class InputFileMeta:
    path: Path
    file_type: str
    project_id: str
    cell_line: str
    sample_id: str
    dataset: str
    relative_source: str


@dataclass
class PeakRecord:
    dataset: str
    species_id: str
    cell_line: str
    sample_id: str
    project_id: str
    chrom: str
    start: str
    end: str
    pileup: str
    p_value: str
    q_value: str
    fold_enrichment: str
    peak_name: str
    peak_region: str
    peak_source_file: str

    def as_dict(self) -> dict[str, str]:
        return {
            "dataset": self.dataset,
            "species_id": self.species_id,
            "cell_line": self.cell_line,
            "sample_id": self.sample_id,
            "project_id": self.project_id,
            "chr": self.chrom,
            "start": self.start,
            "end": self.end,
            "pileup": self.pileup,
            "p_value": self.p_value,
            "q_value": self.q_value,
            "fold_enrichment": self.fold_enrichment,
            "peak_name": self.peak_name,
            "peak_region": self.peak_region,
            "peak_source_file": self.peak_source_file,
        }


@dataclass
class FileSummary:
    dataset: str
    cell_line: str
    project_id: str
    sample_id: str
    file_type: str
    peak_count: int
    source_file: str

    def as_dict(self) -> dict[str, str]:
        return {
            "dataset": self.dataset,
            "cell_line": self.cell_line,
            "project_id": self.project_id,
            "sample_id": self.sample_id,
            "file_type": self.file_type,
            "peak_count": str(self.peak_count),
            "source_file": self.source_file,
        }


@dataclass
class AggregateSummaryRow:
    summary_level: str
    cell_line: str
    project_id: str
    sample_id: str
    dataset_count: int
    sample_count: int
    source_file_count: int
    peak_count: int

    def as_dict(self) -> dict[str, str]:
        return {
            "summary_level": self.summary_level,
            "cell_line": self.cell_line,
            "project_id": self.project_id,
            "sample_id": self.sample_id,
            "dataset_count": str(self.dataset_count),
            "sample_count": str(self.sample_count),
            "source_file_count": str(self.source_file_count),
            "peak_count": str(self.peak_count),
        }


@dataclass
class CellLineMetaRow:
    cell_line: str
    dataset_count: int
    peak_count: int
    sort_order: int
    updated_at: str

    def as_dict(self) -> dict[str, str]:
        return {
            "cell_line": self.cell_line,
            "dataset_count": str(self.dataset_count),
            "peak_count": str(self.peak_count),
            "sort_order": str(self.sort_order),
            "updated_at": self.updated_at,
        }


@dataclass
class DatasetMetaRow:
    cell_line: str
    dataset: str
    project_id: str
    sample_count: int
    peak_count: int
    updated_at: str

    def as_dict(self) -> dict[str, str]:
        return {
            "cell_line": self.cell_line,
            "dataset": self.dataset,
            "project_id": self.project_id,
            "sample_count": str(self.sample_count),
            "peak_count": str(self.peak_count),
            "updated_at": self.updated_at,
        }


DEFAULT_DATABASE_NAME = "big_industry_cell_atac_results_ds0005"
DEFAULT_CELL_LINE_META_TABLE_NAME = "tb_atac_cell_line_meta"
DEFAULT_DATASET_META_TABLE_NAME = "tb_atac_dataset_meta"
DEFAULT_TABLE_NAME_PREFIX = "tb_atac_peak_"
CELL_LINE_METADATA = {
    "hek293": ("HEK293", "hek293"),
    "cho": ("CHO", "cho"),
    "hela": ("HeLa", "hela"),
    "h9": ("H9", "h9"),
    "cap": ("CAP", "cap"),
    "mrc-5": ("MRC-5", "mrc-5"),
    "ht1080": ("HT1080", "ht1080"),
    "wi-38": ("WI-38", "wi-38"),
    "pk-15": ("PK-15", "pk-15"),
    "vero": ("Vero", "vero"),
    "bt": ("BT", "bt"),
    "mdck": ("MDCK", "mdck"),
    "sp2/0": ("SP2/0", "sp2-0"),
    "sp2-0": ("SP2/0", "sp2-0"),
    "cef": ("CEF", "cef"),
    "df-1": ("DF-1", "df-1"),
    "huh-7": ("Huh-7", "huh-7"),
    "jurkat": ("Jurkat", "jurkat"),
    "k562": ("K562", "k562"),
    "sk-br-3": ("SK-BR-3", "sk-br-3"),
    "2bs": ("2BS", "2bs"),
}


@dataclass(frozen=True)
class TargetCellLine:
    cell_line: str
    slug: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare merged ATAC import files and SQL templates for DS0005."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help="Root directory containing epigenome project result folders.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Root directory where the cell-line-specific output directory will be created.",
    )
    parser.add_argument(
        "--input-format",
        choices=("xls", "narrowpeak", "auto"),
        default="auto",
        help="Input file format. auto prefers peaks.xls and falls back to narrowPeak.",
    )
    parser.add_argument(
        "--format",
        choices=("tsv", "csv"),
        default="tsv",
        help="Merged output delimiter format. Default: tsv.",
    )
    parser.add_argument(
        "--species-id",
        default="3",
        help="species_id written to the merged output. Default: 3.",
    )
    parser.add_argument(
        "--cell-lines",
        nargs="+",
        required=True,
        help="Exactly one cell line to process, for example --cell-lines HEK293.",
    )
    parser.add_argument(
        "--dataset-template",
        default="{project_id}_{sample_id}",
        help="Python format string for dataset names.",
    )
    parser.add_argument(
        "--peak-region-default",
        default="Unannotated",
        help="Default peak_region value written for every row.",
    )
    parser.add_argument(
        "--chrom-prefix",
        default="",
        help="Optional chromosome prefix to prepend when missing, for example chr.",
    )
    parser.add_argument(
        "--project-regex",
        default=DEFAULT_PROJECT_REGEX,
        help="Regex used to parse project directory names.",
    )
    parser.add_argument(
        "--exclude-projects",
        nargs="+",
        default=(),
        help=(
            "Optional project IDs to skip during the merge, "
            "e.g. --exclude-projects PRJNA682619."
        ),
    )
    parser.add_argument(
        "--encoding",
        default="utf-8",
        help="Input file encoding. Default: utf-8.",
    )
    parser.add_argument(
        "--database-name",
        default=DEFAULT_DATABASE_NAME,
        help=f"Target schema name. Default: {DEFAULT_DATABASE_NAME}.",
    )
    parser.add_argument(
        "--cell-line-meta-table-name",
        default=DEFAULT_CELL_LINE_META_TABLE_NAME,
        help=(
            "Target cell line metadata table name. "
            f"Default: {DEFAULT_CELL_LINE_META_TABLE_NAME}."
        ),
    )
    parser.add_argument(
        "--dataset-meta-table-name",
        default=DEFAULT_DATASET_META_TABLE_NAME,
        help=(
            "Target dataset metadata table name. "
            f"Default: {DEFAULT_DATASET_META_TABLE_NAME}."
        ),
    )
    parser.add_argument(
        "--merged-filename",
        help="Optional override for the merged import filename.",
    )
    parser.add_argument(
        "--manifest-filename",
        help="Optional override for the manifest summary filename.",
    )
    parser.add_argument(
        "--cell-line-meta-filename",
        help="Optional override for the generated cell line metadata filename.",
    )
    parser.add_argument(
        "--dataset-meta-filename",
        help="Optional override for the generated dataset metadata filename.",
    )
    parser.add_argument(
        "--load-sql-filename",
        help="Optional override for the generated LOAD DATA SQL filename.",
    )
    parser.add_argument(
        "--check-sql-filename",
        help="Optional override for the generated validation SQL filename.",
    )
    parser.add_argument(
        "--sql-data-dir",
        default="/data/industry_cellline_data/ATAC_peak_data",
        help=(
            "When set, replaces the output-root prefix in LOAD DATA LOCAL INFILE paths "
            "with this directory. Default: /data/industry_cellline_data/ATAC_peak_data."
        ),
    )
    return parser.parse_args()


def resolve_target_cell_line(cell_lines: Sequence[str]) -> TargetCellLine:
    normalized = [
        item.strip()
        for item in cell_lines
        if item is not None and item.strip()
    ]
    if len(normalized) != 1:
        raise SystemExit(
            "--cell-lines must contain exactly one cell line in sharded mode, "
            "for example --cell-lines HEK293."
        )

    key = normalized[0].lower()
    metadata = CELL_LINE_METADATA.get(key)
    if metadata is None:
        raise SystemExit(
            f"Unsupported cell line for sharded ATAC import: {normalized[0]}\n"
            "Please add its cell_line -> slug mapping to this script first."
        )
    return TargetCellLine(cell_line=metadata[0], slug=metadata[1])


def normalize_cell_line_filters(target: TargetCellLine) -> set[str]:
    return {
        target.cell_line.strip().lower(),
        target.slug.strip().lower(),
    }


def should_keep_cell_line(cell_line: str, cell_line_filters: set[str]) -> bool:
    if not cell_line_filters:
        return True
    return cell_line.strip().lower() in cell_line_filters


def output_directory_for(output_root: Path, target: TargetCellLine) -> Path:
    safe_name = target.cell_line.replace("/", "_").replace("\\", "_")
    return output_root.resolve() / safe_name


def default_filename(filename: str | None, generated_name: str) -> str:
    if filename and filename.strip():
        return filename
    return generated_name


def extract_project_meta(file_path: Path, project_pattern: re.Pattern[str]) -> tuple[str, str]:
    for parent in file_path.parents:
        match = project_pattern.match(parent.name)
        if match:
            return match.group(1), match.group(2)
    raise SystemExit(
        f"Could not parse project_id and cell_line from path: {file_path}\n"
        f"Expected a directory name matching regex: {project_pattern.pattern}"
    )


def derive_sample_id(file_name: str) -> str:
    for suffix in (XLS_SUFFIX, NARROWPEAK_SUFFIX):
        if file_name.endswith(suffix):
            return file_name[: -len(suffix)]
    return Path(file_name).stem


def build_input_meta(
    file_path: Path,
    file_type: str,
    input_root: Path,
    project_pattern: re.Pattern[str],
    dataset_template: str,
) -> InputFileMeta:
    project_id, cell_line = extract_project_meta(file_path, project_pattern)
    sample_id = derive_sample_id(file_path.name)

    try:
        relative_source = str(file_path.resolve().relative_to(input_root))
    except ValueError:
        relative_source = file_path.name

    dataset = dataset_template.format(
        cell_line=cell_line,
        project_id=project_id,
        sample_id=sample_id,
    )

    return InputFileMeta(
        path=file_path,
        file_type=file_type,
        project_id=project_id,
        cell_line=cell_line,
        sample_id=sample_id,
        dataset=dataset,
        relative_source=relative_source,
    )


def discover_input_files(args: argparse.Namespace) -> List[InputFileMeta]:
    input_root = args.input_root.resolve()
    if not input_root.exists():
        raise SystemExit(f"Input root does not exist: {input_root}")

    project_pattern = re.compile(args.project_regex)
    target_cell_line = resolve_target_cell_line(args.cell_lines)
    cell_line_filters = normalize_cell_line_filters(target_cell_line)
    exclude_projects = set(args.exclude_projects or ())
    files_by_sample_key: dict[tuple[str, str, str], InputFileMeta] = {}
    scan_specs = []
    if args.input_format in {"xls", "auto"}:
        scan_specs.append(("xls", XLS_SUFFIX, True))
    if args.input_format in {"narrowpeak", "auto"}:
        scan_specs.append(("narrowpeak", NARROWPEAK_SUFFIX, args.input_format == "narrowpeak"))

    # Only <project_dir>/1_result/3_peak_calling is treated as the analysis
    # result location; recursive scans may pick up stale or backup files.
    for project_dir in sorted(
        child
        for child in input_root.iterdir()
        if child.is_dir() and project_pattern.match(child.name)
    ):
        peak_dir = project_dir / "1_result" / "3_peak_calling"
        if not peak_dir.is_dir():
            continue
        for file_type, suffix, should_override in scan_specs:
            for file_path in sorted(peak_dir.glob(f"*{suffix}")):
                meta = build_input_meta(
                    file_path=file_path,
                    file_type=file_type,
                    input_root=input_root,
                    project_pattern=project_pattern,
                    dataset_template=args.dataset_template,
                )
                if meta.project_id in exclude_projects:
                    continue
                if not should_keep_cell_line(meta.cell_line, cell_line_filters):
                    continue
                key = (meta.project_id, meta.cell_line, meta.sample_id)
                if should_override or key not in files_by_sample_key:
                    files_by_sample_key[key] = meta

    files = sorted(files_by_sample_key.values(), key=lambda item: (item.cell_line, item.project_id, item.sample_id))
    if not files:
        raise SystemExit(
            f"No peak files found under {input_root} for input format {args.input_format}."
        )
    return files


def normalize_chromosome(raw_value: str, chrom_prefix: str) -> str:
    chrom = raw_value.strip()
    if not chrom_prefix:
        return chrom
    if chrom.startswith(chrom_prefix):
        return chrom
    return f"{chrom_prefix}{chrom}"


def safe_int_text(value: str) -> str:
    return str(int(float(value)))


def iter_data_lines(path: Path, encoding: str) -> Iterable[str]:
    with path.open("r", encoding=encoding, newline="") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if stripped and not stripped.startswith("#"):
                yield stripped


def build_peak_record(
    meta: InputFileMeta,
    species_id: str,
    peak_region_default: str,
    chrom_prefix: str,
    *,
    chrom: str,
    start: str,
    end: str,
    pileup: str,
    p_value: str,
    q_value: str,
    fold_enrichment: str,
    peak_name: str,
) -> PeakRecord:
    return PeakRecord(
        dataset=meta.dataset,
        species_id=str(species_id),
        cell_line=meta.cell_line,
        sample_id=meta.sample_id,
        project_id=meta.project_id,
        chrom=normalize_chromosome(chrom, chrom_prefix),
        start=safe_int_text(start),
        end=safe_int_text(end),
        pileup=pileup,
        p_value=p_value,
        q_value=q_value,
        fold_enrichment=fold_enrichment,
        peak_name=peak_name,
        peak_region=peak_region_default,
        peak_source_file=meta.relative_source,
    )


def read_macs_xls(
    meta: InputFileMeta,
    species_id: str,
    peak_region_default: str,
    chrom_prefix: str,
    encoding: str,
) -> List[PeakRecord]:
    lines = list(iter_data_lines(meta.path, encoding))

    if not lines:
        raise SystemExit(f"{meta.path} contains no data rows.")

    header = lines[0].split()
    if header != EXPECTED_XLS_HEADER:
        raise SystemExit(
            f"{meta.path} has an unexpected peaks.xls header:\n"
            f"Found:    {header}\n"
            f"Expected: {EXPECTED_XLS_HEADER}"
        )

    records: List[PeakRecord] = []
    for line in lines[1:]:
        fields = line.split()
        if len(fields) != len(EXPECTED_XLS_HEADER):
            raise SystemExit(
                f"{meta.path} has a malformed data row with {len(fields)} columns:\n{line}"
            )

        chrom, start, end, _length, _summit, pileup, p_value, fold_enrichment, q_value, peak_name = fields
        records.append(
            build_peak_record(
                meta,
                species_id,
                peak_region_default,
                chrom_prefix,
                chrom=chrom,
                start=start,
                end=end,
                pileup=pileup,
                p_value=p_value,
                q_value=q_value,
                fold_enrichment=fold_enrichment,
                peak_name=peak_name,
            )
        )

    return records


def read_narrowpeak(
    meta: InputFileMeta,
    species_id: str,
    peak_region_default: str,
    chrom_prefix: str,
    encoding: str,
) -> List[PeakRecord]:
    records: List[PeakRecord] = []
    for line in iter_data_lines(meta.path, encoding):
        fields = line.split()
        if len(fields) < 10:
            raise SystemExit(
                f"{meta.path} has a malformed narrowPeak row with {len(fields)} columns:\n{line}"
            )

        chrom, chrom_start, chrom_end, peak_name, _score, _strand, signal_value, p_value, q_value, _peak = fields[:10]
        records.append(
            build_peak_record(
                meta,
                species_id,
                peak_region_default,
                chrom_prefix,
                chrom=chrom,
                start=str(int(chrom_start) + 1),
                end=chrom_end,
                pileup="\\N",
                p_value=p_value,
                q_value=q_value,
                fold_enrichment=signal_value,
                peak_name=peak_name,
            )
        )

    return records


def read_peak_records(
    meta: InputFileMeta,
    species_id: str,
    peak_region_default: str,
    chrom_prefix: str,
    encoding: str,
) -> tuple[List[PeakRecord], FileSummary]:
    if meta.file_type == "xls":
        records = read_macs_xls(meta, species_id, peak_region_default, chrom_prefix, encoding)
    elif meta.file_type == "narrowpeak":
        records = read_narrowpeak(meta, species_id, peak_region_default, chrom_prefix, encoding)
    else:
        raise SystemExit(f"Unsupported file type: {meta.file_type}")

    summary = FileSummary(
        dataset=meta.dataset,
        cell_line=meta.cell_line,
        project_id=meta.project_id,
        sample_id=meta.sample_id,
        file_type=meta.file_type,
        peak_count=len(records),
        source_file=meta.relative_source,
    )
    return records, summary


def write_output(records: Iterable[PeakRecord], output_path: Path, output_format: str) -> None:
    write_dict_rows(
        output_path,
        OUTPUT_COLUMNS,
        (record.as_dict() for record in records),
        output_format,
    )


def group_file_summaries(
    manifest_rows: Sequence[FileSummary],
    key_fn,
) -> List[tuple[tuple[str, ...], List[FileSummary]]]:
    grouped: dict[tuple[str, ...], List[FileSummary]] = {}
    for item in manifest_rows:
        key = key_fn(item)
        grouped.setdefault(key, []).append(item)
    return sorted(grouped.items(), key=lambda item: item[0])


def build_aggregate_summary_rows(
    manifest_rows: Sequence[FileSummary],
) -> List[AggregateSummaryRow]:
    rows: List[AggregateSummaryRow] = []
    for summary_level, key_fn, field_names in SUMMARY_GROUP_SPECS:
        for key, items in group_file_summaries(manifest_rows, key_fn=key_fn):
            key_values = dict(zip(field_names, key))
            rows.append(
                AggregateSummaryRow(
                    summary_level=summary_level,
                    cell_line=key_values.get("cell_line", ""),
                    project_id=key_values.get("project_id", ""),
                    sample_id=key_values.get("sample_id", ""),
                    dataset_count=len({item.dataset for item in items}),
                    sample_count=len({(item.cell_line, item.project_id, item.sample_id) for item in items}),
                    source_file_count=len(items),
                    peak_count=sum(item.peak_count for item in items),
                )
            )
    return rows


def delimiter_for(output_format: str) -> str:
    return "\t" if output_format == "tsv" else ","


def write_dict_rows(
    output_path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[dict[str, str]],
    output_format: str,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter_for(output_format),
            quoting=csv.QUOTE_MINIMAL,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_manifest(
    manifest_rows: Sequence[FileSummary],
    output_path: Path,
    output_format: str,
) -> None:
    write_dict_rows(
        output_path,
        MANIFEST_COLUMNS,
        (row.as_dict() for row in build_aggregate_summary_rows(manifest_rows)),
        output_format,
    )


def build_cell_line_meta_rows(file_summaries: Sequence[FileSummary]) -> List[CellLineMetaRow]:
    rows: List[CellLineMetaRow] = []
    updated_at = "CURRENT_TIMESTAMP"
    grouped = group_file_summaries(file_summaries, key_fn=lambda item: (item.cell_line,))
    for sort_order, ((cell_line,), items) in enumerate(grouped, start=1):
        rows.append(
            CellLineMetaRow(
                cell_line=cell_line,
                dataset_count=len({item.dataset for item in items}),
                peak_count=sum(item.peak_count for item in items),
                sort_order=sort_order,
                updated_at=updated_at,
            )
        )
    return rows


def build_dataset_meta_rows(file_summaries: Sequence[FileSummary]) -> List[DatasetMetaRow]:
    rows: List[DatasetMetaRow] = []
    updated_at = "CURRENT_TIMESTAMP"
    grouped = group_file_summaries(
        file_summaries,
        key_fn=lambda item: (item.cell_line, item.dataset, item.project_id),
    )
    for (cell_line, dataset, project_id), items in grouped:
        rows.append(
            DatasetMetaRow(
                cell_line=cell_line,
                dataset=dataset,
                project_id=project_id,
                sample_count=len({item.sample_id for item in items}),
                peak_count=sum(item.peak_count for item in items),
                updated_at=updated_at,
            )
        )
    return rows


def write_cell_line_meta(
    file_summaries: Sequence[FileSummary],
    output_path: Path,
    output_format: str,
) -> None:
    write_dict_rows(
        output_path,
        CELL_LINE_META_COLUMNS,
        (row.as_dict() for row in build_cell_line_meta_rows(file_summaries)),
        output_format,
    )


def write_dataset_meta(
    file_summaries: Sequence[FileSummary],
    output_path: Path,
    output_format: str,
) -> None:
    write_dict_rows(
        output_path,
        DATASET_META_COLUMNS,
        (row.as_dict() for row in build_dataset_meta_rows(file_summaries)),
        output_format,
    )


def print_summary(records: Sequence[PeakRecord], file_summaries: Sequence[FileSummary]) -> None:
    file_count = len(file_summaries)
    peak_count = len(records)
    project_counter = Counter(item.project_id for item in file_summaries)
    cell_line_counter = Counter(item.cell_line for item in file_summaries)

    print(f"Discovered peak result files: {file_count}")
    print(f"Total peak rows written:      {peak_count}")
    print("")
    print("Files by cell line:")
    for cell_line, count in sorted(cell_line_counter.items()):
        print(f"  {cell_line}: {count}")
    print("")
    print("Files by project:")
    for project_id, count in sorted(project_counter.items()):
        print(f"  {project_id}: {count}")


def to_sql_path_text(path: Path) -> str:
    return str(path.resolve()).replace("\\", "/")


def _sql_data_file_path(data_path: Path, output_root: Path, sql_data_dir: str | None) -> str:
    if sql_data_dir:
        return str(Path(sql_data_dir) / data_path.resolve().relative_to(output_root.resolve()))
    return to_sql_path_text(data_path)


def write_load_sql(
    sql_path: Path,
    data_path: Path,
    cell_line_meta_path: Path,
    dataset_meta_path: Path,
    database_name: str,
    table_name: str,
    cell_line_meta_table_name: str,
    dataset_meta_table_name: str,
    sql_data_dir: str | None = None,
    output_root: Path | None = None,
) -> None:
    sql_path.parent.mkdir(parents=True, exist_ok=True)
    sql_text = f"""USE `{database_name}`;

LOAD DATA LOCAL INFILE '{_sql_data_file_path(data_path, output_root, sql_data_dir)}'
INTO TABLE `{table_name}`
CHARACTER SET utf8mb4
FIELDS TERMINATED BY '\\t'
OPTIONALLY ENCLOSED BY '"'
LINES TERMINATED BY '\\n'
IGNORE 1 LINES
(
  dataset,
  species_id,
  cell_line,
  sample_id,
  project_id,
  chr,
  start,
  end,
  pileup,
  p_value,
  q_value,
  fold_enrichment,
  peak_name,
  peak_region,
  peak_source_file
);

CREATE TEMPORARY TABLE `tmp_atac_cell_line_meta` LIKE `{cell_line_meta_table_name}`;

LOAD DATA LOCAL INFILE '{_sql_data_file_path(cell_line_meta_path, output_root, sql_data_dir)}'
INTO TABLE `tmp_atac_cell_line_meta`
CHARACTER SET utf8mb4
FIELDS TERMINATED BY '\\t'
OPTIONALLY ENCLOSED BY '"'
LINES TERMINATED BY '\\n'
IGNORE 1 LINES
(
  cell_line,
  dataset_count,
  peak_count,
  sort_order,
  @updated_at
)
SET updated_at = NOW();

INSERT INTO `{cell_line_meta_table_name}` (
  cell_line,
  dataset_count,
  peak_count,
  sort_order,
  updated_at
)
SELECT
  cell_line,
  dataset_count,
  peak_count,
  sort_order,
  updated_at
FROM `tmp_atac_cell_line_meta`
ON DUPLICATE KEY UPDATE
  dataset_count = `{cell_line_meta_table_name}`.dataset_count + VALUES(dataset_count),
  peak_count = `{cell_line_meta_table_name}`.peak_count + VALUES(peak_count),
  sort_order = LEAST(`{cell_line_meta_table_name}`.sort_order, VALUES(sort_order)),
  updated_at = NOW();

DROP TEMPORARY TABLE `tmp_atac_cell_line_meta`;

CREATE TEMPORARY TABLE `tmp_atac_dataset_meta` LIKE `{dataset_meta_table_name}`;

LOAD DATA LOCAL INFILE '{_sql_data_file_path(dataset_meta_path, output_root, sql_data_dir)}'
INTO TABLE `tmp_atac_dataset_meta`
CHARACTER SET utf8mb4
FIELDS TERMINATED BY '\\t'
OPTIONALLY ENCLOSED BY '"'
LINES TERMINATED BY '\\n'
IGNORE 1 LINES
(
  cell_line,
  dataset,
  project_id,
  sample_count,
  peak_count,
  @updated_at
)
SET updated_at = NOW();

INSERT INTO `{dataset_meta_table_name}` (
  cell_line,
  dataset,
  project_id,
  sample_count,
  peak_count,
  updated_at
)
SELECT
  cell_line,
  dataset,
  project_id,
  sample_count,
  peak_count,
  updated_at
FROM `tmp_atac_dataset_meta`
ON DUPLICATE KEY UPDATE
  sample_count = `{dataset_meta_table_name}`.sample_count + VALUES(sample_count),
  peak_count = `{dataset_meta_table_name}`.peak_count + VALUES(peak_count),
  updated_at = NOW();

DROP TEMPORARY TABLE `tmp_atac_dataset_meta`;
"""
    sql_path.write_text(sql_text, encoding="utf-8")


def write_check_sql(
    sql_path: Path,
    database_name: str,
    table_name: str,
    cell_line_meta_table_name: str,
    dataset_meta_table_name: str,
) -> None:
    sql_path.parent.mkdir(parents=True, exist_ok=True)
    sql_text = f"""USE `{database_name}`;

SELECT COUNT(*) AS total_rows
FROM `{table_name}`;

SELECT dataset, COUNT(*) AS peak_rows
FROM `{table_name}`
GROUP BY dataset
ORDER BY dataset;

SELECT cell_line, COUNT(*) AS peak_rows
FROM `{table_name}`
GROUP BY cell_line
ORDER BY cell_line;

SELECT *
FROM `{cell_line_meta_table_name}`
ORDER BY sort_order, cell_line;

SELECT *
FROM `{dataset_meta_table_name}`
ORDER BY cell_line, dataset;

SELECT project_id, COUNT(*) AS peak_rows
FROM `{table_name}`
GROUP BY project_id
ORDER BY project_id;

SELECT chr, COUNT(*) AS peak_rows
FROM `{table_name}`
GROUP BY chr
ORDER BY chr;

SELECT *
FROM `{table_name}`
ORDER BY peak_id DESC
LIMIT 20;
"""
    sql_path.write_text(sql_text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    target_cell_line = resolve_target_cell_line(args.cell_lines)
    output_dir = output_directory_for(args.output_root, target_cell_line)
    table_name = f"{DEFAULT_TABLE_NAME_PREFIX}{target_cell_line.slug}"
    data_ext = "tsv" if args.format == "tsv" else "csv"

    merged_output = output_dir / default_filename(
        args.merged_filename,
        f"merged_atac_peaks_ds0005_{target_cell_line.slug}.{data_ext}",
    )
    manifest_output = output_dir / default_filename(
        args.manifest_filename,
        f"merged_atac_peak_manifest_ds0005_{target_cell_line.slug}.{data_ext}",
    )
    cell_line_meta_output = output_dir / default_filename(
        args.cell_line_meta_filename,
        f"atac_cell_line_meta_ds0005_{target_cell_line.slug}.{data_ext}",
    )
    dataset_meta_output = output_dir / default_filename(
        args.dataset_meta_filename,
        f"atac_dataset_meta_ds0005_{target_cell_line.slug}.{data_ext}",
    )
    load_sql_output = output_dir / default_filename(
        args.load_sql_filename,
        f"load_{table_name}_ds0005.sql",
    )
    check_sql_output = output_dir / default_filename(
        args.check_sql_filename,
        f"validate_{table_name}_ds0005.sql",
    )

    files = discover_input_files(args)

    merged_records: List[PeakRecord] = []
    file_summaries: List[FileSummary] = []
    for meta in files:
        records, summary = read_peak_records(
            meta=meta,
            species_id=str(args.species_id),
            peak_region_default=args.peak_region_default,
            chrom_prefix=args.chrom_prefix,
            encoding=args.encoding,
        )
        merged_records.extend(records)
        file_summaries.append(summary)

    write_output(merged_records, merged_output, args.format)
    write_manifest(file_summaries, manifest_output, args.format)
    write_cell_line_meta(file_summaries, cell_line_meta_output, args.format)
    write_dataset_meta(file_summaries, dataset_meta_output, args.format)
    write_load_sql(
        sql_path=load_sql_output,
        data_path=merged_output,
        cell_line_meta_path=cell_line_meta_output,
        dataset_meta_path=dataset_meta_output,
        database_name=args.database_name,
        table_name=table_name,
        cell_line_meta_table_name=args.cell_line_meta_table_name,
        dataset_meta_table_name=args.dataset_meta_table_name,
        sql_data_dir=args.sql_data_dir,
        output_root=args.output_root,
    )
    write_check_sql(
        sql_path=check_sql_output,
        database_name=args.database_name,
        table_name=table_name,
        cell_line_meta_table_name=args.cell_line_meta_table_name,
        dataset_meta_table_name=args.dataset_meta_table_name,
    )

    print_summary(merged_records, file_summaries)
    print("")
    print(f"Target cell line:       {target_cell_line.cell_line}")
    print(f"Target shard table:     {table_name}")
    print(f"Output directory:       {output_dir}")
    print("")
    print(f"Merged output written: {merged_output}")
    print(f"Manifest written:      {manifest_output}")
    print(f"Cell line meta written:{cell_line_meta_output}")
    print(f"Dataset meta written:  {dataset_meta_output}")
    print(f"Load SQL written:      {load_sql_output}")
    print(f"Check SQL written:     {check_sql_output}")


if __name__ == "__main__":
    main()
