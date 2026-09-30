#!/usr/bin/env python3
"""Prepare DS0002 transcriptome comparison import assets for one cell line.

The script scans one cell-line DEG result directory with this structure:

    HEK293/
      PRJEB46065_ADAR_KO_vs_WT/
        diff_genes.csv
        count_transformation_normTransform.csv
        count_transformation_vst.csv
        count_transformation_rlog.csv
        PRJEB46065_ADAR_KO_vs_WT.tsv
        ...

It generates:
1. One import-ready shard TSV for the selected cell line.
2. Metadata TSVs for tb_comparison_cell_line_meta and tb_comparison_study_meta.
3. LOAD DATA SQL to import all generated TSVs.
4. Validation SQL and a manifest TSV for quick inspection.

The generated shard TSV is aligned with the current DS0002 shard schema and
keeps normTransform, VST, and rlog expression arrays in the same result table.

example usage:
python prepare_ds0002_transcriptome_import.py \
  --groups-root /path/to/groups/HEK293 \
  --output-dir ./ds0002_hek293_output \
  --cell-line-metadata cell_line_metadata.tsv \
  --cell-line HEK293 \
  --sql-data-dir /data/industry_cellline_data/trs_data/hek293
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence


DEFAULT_DATABASE_NAME = "big_industry_cell_analysis_results_ds0002"
DEFAULT_RESULT_TABLE_PREFIX = "tb_comparison_gene_exp_"
DEFAULT_CELL_LINE_META_TABLE = "tb_comparison_cell_line_meta"
DEFAULT_STUDY_META_TABLE = "tb_comparison_study_meta"
DEFAULT_COMPARISON_DIR_GLOB = "*"
DEFAULT_REQUIRED_FILES = (
    "diff_genes.csv",
    "count_transformation_normTransform.csv",
    "count_transformation_vst.csv",
    "count_transformation_rlog.csv",
)
EXCLUDED_CELL_LINES = {"cap", "2bs"}
SHARD_OUTPUT_COLUMNS = [
    "gene_name",
    "gene_alias",
    "cellline_name",
    "bioproject_accession",
    "sample_list_control",
    "sample_list_experimental",
    "avg_norm_transform_control",
    "avg_norm_transform_experimental",
    "median_norm_transform_control",
    "median_norm_transform_experimental",
    "norm_transform_list_control",
    "norm_transform_list_experimental",
    "avg_vst_control",
    "avg_vst_experimental",
    "median_vst_control",
    "median_vst_experimental",
    "vst_list_control",
    "vst_list_experimental",
    "avg_rlog_control",
    "avg_rlog_experimental",
    "median_rlog_control",
    "median_rlog_experimental",
    "rlog_list_control",
    "rlog_list_experimental",
    "fold_change_value",
    "p_value",
    "comparison",
    "significant",
    "q_value",
    "ontology_id",
    "tissue",
    "primary_experimental",
]
CELL_LINE_META_COLUMNS = [
    "cell_line",
    "cell_line_slug",
    "study_count",
    "comparison_count",
    "gene_row_count",
    "sort_order",
]
STUDY_META_COLUMNS = [
    "cell_line",
    "cell_line_slug",
    "bioproject_accession",
    "comparison",
    "gene_row_count",
    "sort_order",
]
MANIFEST_COLUMNS = [
    "cell_line",
    "cell_line_slug",
    "bioproject_accession",
    "comparison",
    "comparison_dir",
    "control_sample_count",
    "experimental_sample_count",
    "gene_row_count",
]


@dataclass(frozen=True)
class CellLineInfo:
    display_name: str
    slug: str
    sort_order: int


@dataclass
class ComparisonBundle:
    cell_line: CellLineInfo
    bioproject_accession: str
    comparison: str
    comparison_dir: Path
    control_samples: List[str]
    experimental_samples: List[str]
    gene_rows: List[dict[str, str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate DS0002 import TSVs and LOAD DATA SQL for one cell line."
    )
    parser.add_argument(
        "--groups-root",
        type=Path,
        required=True,
        help="Directory of the target cell line under the DEG/groups tree, e.g. .../groups/HEK293",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for generated TSV and SQL files.",
    )
    parser.add_argument(
        "--sql-data-dir",
        type=Path,
        help=(
            "Root path used in generated LOAD DATA LOCAL INFILE statements after "
            "--output-dir is copied to the database/import server. If omitted, "
            "local absolute paths under --output-dir are used."
        ),
    )
    parser.add_argument(
        "--cell-line-metadata",
        type=Path,
        required=True,
        help="Path to cell_line_metadata.tsv.",
    )
    parser.add_argument(
        "--cell-line",
        required=True,
        help="Cell line name or slug for the selected --groups-root directory, e.g. HEK293 or hek293.",
    )
    parser.add_argument(
        "--comparison-dir-glob",
        default=DEFAULT_COMPARISON_DIR_GLOB,
        help="Glob used to find comparison directories under each cell-line directory. Default: *",
    )
    parser.add_argument(
        "--padj-cutoff",
        type=float,
        default=0.05,
        help="Significance cutoff for q_value / padj. Default: 0.05",
    )
    parser.add_argument(
        "--database",
        default=DEFAULT_DATABASE_NAME,
        help=f"Target MySQL database name. Default: {DEFAULT_DATABASE_NAME}",
    )
    parser.add_argument(
        "--result-table-prefix",
        default=DEFAULT_RESULT_TABLE_PREFIX,
        help=f"Result shard table prefix. Default: {DEFAULT_RESULT_TABLE_PREFIX}",
    )
    parser.add_argument(
        "--cell-line-meta-table",
        default=DEFAULT_CELL_LINE_META_TABLE,
        help=f"Cell line metadata table name. Default: {DEFAULT_CELL_LINE_META_TABLE}",
    )
    parser.add_argument(
        "--study-meta-table",
        default=DEFAULT_STUDY_META_TABLE,
        help=f"Study metadata table name. Default: {DEFAULT_STUDY_META_TABLE}",
    )
    parser.add_argument(
        "--gene-alias-map",
        type=Path,
        help="Optional TSV/CSV containing gene_name and gene_alias columns.",
    )
    return parser.parse_args()


def read_cell_line_metadata(path: Path) -> Dict[str, CellLineInfo]:
    if not path.exists():
        raise SystemExit(f"cell_line_metadata.tsv was not found: {path}")
    metadata: Dict[str, CellLineInfo] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            display_name = (row.get("cell_line_name") or "").strip()
            slug = (row.get("cell_line_slug") or "").strip()
            dataset_id = (row.get("dataset_id") or "").strip()
            if not display_name or not slug:
                continue
            if slug.lower() in EXCLUDED_CELL_LINES:
                continue
            sort_order = int(dataset_id) if dataset_id.isdigit() else 0
            info = CellLineInfo(display_name=display_name, slug=slug, sort_order=sort_order)
            metadata[display_name.lower()] = info
            metadata[slug.lower()] = info
    return metadata


def read_gene_alias_map(path: Path | None) -> Dict[str, str]:
    if path is None:
        return {}
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    alias_map: Dict[str, str] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {"gene_name", "gene_alias"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise SystemExit(f"{path} must contain columns: gene_name, gene_alias")
        for row in reader:
            gene_name = (row.get("gene_name") or "").strip()
            gene_alias = (row.get("gene_alias") or "").strip()
            if gene_name:
                alias_map[gene_name] = gene_alias
    return alias_map


def resolve_cell_line_target(args: argparse.Namespace, metadata: Dict[str, CellLineInfo]) -> tuple[CellLineInfo, Path]:
    info = metadata.get(args.cell_line.strip().lower())
    if info is None:
        raise SystemExit(f"Unsupported or excluded cell line: {args.cell_line}")
    if not args.groups_root.is_dir():
        raise SystemExit(f"--groups-root is not a directory: {args.groups_root}")
    return info, args.groups_root


def find_comparison_dirs(cell_line_dir: Path, pattern: str) -> List[Path]:
    dirs: List[Path] = []
    for child in sorted(cell_line_dir.glob(pattern)):
        if not child.is_dir():
            continue
        missing = [name for name in DEFAULT_REQUIRED_FILES if not (child / name).exists()]
        if missing:
            print(f"[WARNING] 跳过不完整的差异分析目录: {child.name}")
            print(f"          缺少文件: {', '.join(missing)}")
            continue
        dirs.append(child)
    return dirs


def normalize_sample_id(sample_id: str) -> str:
    normalized = sample_id.replace("-", ".").strip()
    if normalized and normalized[0].isdigit():
        normalized = f"X{normalized}"
    return normalized


def parse_bioproject_and_comparison(name: str) -> tuple[str, str]:
    if "_" not in name:
        return "", name
    prefix, remainder = name.split("_", 1)
    if prefix.upper().startswith("PRJ"):
        return prefix, name
    return "", name


def read_condition_groups(path: Path) -> tuple[List[str], List[str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "sample" not in reader.fieldnames or "condition" not in reader.fieldnames:
            raise SystemExit(f"{path} must contain sample and condition columns.")
        control_samples: List[str] = []
        experimental_samples: List[str] = []
        for row in reader:
            raw_sample = (row.get("sample") or "").strip()
            condition = (row.get("condition") or "").strip().lower()
            if not raw_sample:
                continue
            sample = normalize_sample_id(raw_sample)
            if condition == "control":
                control_samples.append(sample)
            else:
                experimental_samples.append(sample)
    if not control_samples or not experimental_samples:
        raise SystemExit(f"{path} does not contain both control and experimental samples.")
    return control_samples, experimental_samples


def read_csv_matrix(path: Path) -> tuple[List[str], Dict[str, Dict[str, float]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise SystemExit(f"{path} is empty.") from exc
        if len(header) < 2:
            raise SystemExit(f"{path} does not look like a matrix file.")
        samples = [item.strip().strip('"') for item in header[1:]]
        matrix: Dict[str, Dict[str, float]] = {}
        for row in reader:
            if not row:
                continue
            gene_name = row[0].strip().strip('"')
            if not gene_name:
                continue
            values: Dict[str, float] = {}
            for sample, text in zip(samples, row[1:]):
                value_text = text.strip().strip('"')
                if value_text == "":
                    continue
                values[sample] = float(value_text)
            matrix[gene_name] = values
    return samples, matrix


def read_diff_genes(path: Path) -> List[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise SystemExit(f"{path} is empty.") from exc
        required = {"baseMean", "log2FoldChange", "pvalue", "padj"}
        header_names = [item.strip() for item in header]
        if not required.issubset(set(header_names)):
            raise SystemExit(f"{path} is missing one or more required columns: {sorted(required)}")
        rows: List[dict[str, str]] = []
        for row in reader:
            if not row:
                continue
            record = {header_names[idx]: row[idx].strip() if idx < len(row) else "" for idx in range(len(header_names))}
            gene_name = row[0].strip().strip('"')
            if not gene_name:
                continue
            record["gene_name"] = gene_name
            rows.append(record)
    return rows


def json_list(values: Sequence[float]) -> str:
    return json.dumps([float(f"{value:.10g}") for value in values], separators=(",", ":"))


def json_text_list(values: Sequence[str]) -> str:
    return json.dumps(list(values), ensure_ascii=False, separators=(",", ":"))


def format_float(value: float | None) -> str:
    if value is None:
        return ""
    if math.isnan(value) or math.isinf(value):
        return ""
    return f"{value:.10g}"


def safe_float_text(value: str) -> float:
    text = value.strip()
    if text == "" or text.lower() in {"na", "nan", "null"}:
        return math.inf
    try:
        return float(text)
    except ValueError:
        return math.inf


def padj_to_significant(padj_text: str, cutoff: float) -> str:
    try:
        return "1" if float(padj_text) <= cutoff else "0"
    except ValueError:
        return "0"


def build_metric_fields(values_by_sample: Dict[str, float], control_samples: Sequence[str], experimental_samples: Sequence[str], prefix: str) -> dict[str, str]:
    control_values = [values_by_sample[sample] for sample in control_samples if sample in values_by_sample]
    experimental_values = [values_by_sample[sample] for sample in experimental_samples if sample in values_by_sample]
    if len(control_values) != len(control_samples) or len(experimental_values) != len(experimental_samples):
        missing = [
            sample
            for sample in list(control_samples) + list(experimental_samples)
            if sample not in values_by_sample
        ]
        raise SystemExit(f"Missing matrix values for samples: {', '.join(missing)}")
    return {
        f"avg_{prefix}_control": format_float(statistics.mean(control_values) if control_values else None),
        f"avg_{prefix}_experimental": format_float(statistics.mean(experimental_values) if experimental_values else None),
        f"median_{prefix}_control": format_float(statistics.median(control_values) if control_values else None),
        f"median_{prefix}_experimental": format_float(statistics.median(experimental_values) if experimental_values else None),
        f"{prefix}_list_control": json_list(control_values),
        f"{prefix}_list_experimental": json_list(experimental_values),
    }


def build_comparison_rows(
    info: CellLineInfo,
    comparison_dir: Path,
    gene_alias_map: Dict[str, str],
    padj_cutoff: float,
) -> ComparisonBundle:
    comparison_name = comparison_dir.name
    bioproject_accession, comparison_value = parse_bioproject_and_comparison(comparison_name)
    condition_file = comparison_dir / f"{comparison_name}.tsv"
    if not condition_file.exists():
        raise SystemExit(f"Missing condition file for comparison: {condition_file}")

    control_samples, experimental_samples = read_condition_groups(condition_file)
    _, norm_transform_matrix = read_csv_matrix(comparison_dir / "count_transformation_normTransform.csv")
    _, vst_matrix = read_csv_matrix(comparison_dir / "count_transformation_vst.csv")
    _, rlog_matrix = read_csv_matrix(comparison_dir / "count_transformation_rlog.csv")
    diff_rows = read_diff_genes(comparison_dir / "diff_genes.csv")

    shard_rows: List[dict[str, str]] = []
    for diff_row in diff_rows:
        gene_name = diff_row["gene_name"]
        norm_values = norm_transform_matrix.get(gene_name)
        vst_values = vst_matrix.get(gene_name)
        rlog_values = rlog_matrix.get(gene_name)
        if norm_values is None or vst_values is None or rlog_values is None:
            raise SystemExit(
                f"Gene {gene_name} is missing from one or more transformed count matrices in {comparison_dir}"
            )

        record: dict[str, str] = {
            "gene_name": gene_name,
            "gene_alias": gene_alias_map.get(gene_name, ""),
            "cellline_name": info.display_name,
            "bioproject_accession": bioproject_accession,
            "sample_list_control": json_text_list(control_samples),
            "sample_list_experimental": json_text_list(experimental_samples),
            "fold_change_value": diff_row.get("log2FoldChange", ""),
            "p_value": diff_row.get("pvalue", ""),
            "comparison": comparison_value,
            "significant": padj_to_significant(diff_row.get("padj", ""), padj_cutoff),
            "q_value": diff_row.get("padj", ""),
            "ontology_id": "",
            "tissue": "",
            "primary_experimental": "",
        }
        record.update(build_metric_fields(norm_values, control_samples, experimental_samples, "norm_transform"))
        record.update(build_metric_fields(vst_values, control_samples, experimental_samples, "vst"))
        record.update(build_metric_fields(rlog_values, control_samples, experimental_samples, "rlog"))
        shard_rows.append(record)

    shard_rows.sort(
        key=lambda item: (
            item["comparison"],
            safe_float_text(item["q_value"]),
            safe_float_text(item["p_value"]),
            item["gene_name"],
        )
    )
    return ComparisonBundle(
        cell_line=info,
        bioproject_accession=bioproject_accession,
        comparison=comparison_value,
        comparison_dir=comparison_dir,
        control_samples=control_samples,
        experimental_samples=experimental_samples,
        gene_rows=shard_rows,
    )


def ensure_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, str]]) -> None:
    ensure_output_dir(path.parent)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_load_sql(
    output_dir: Path,
    sql_data_dir: Path | None,
    database: str,
    result_table_prefix: str,
    cell_line_meta_table: str,
    study_meta_table: str,
    shard_files: Dict[str, Path],
    cell_line_rows: List[dict[str, str]],
) -> str:
    lines = [
        f"USE `{database}`;",
        "SET NAMES utf8mb4;",
        "",
    ]
    target_cell_lines = [row["cell_line"] for row in cell_line_rows]
    if target_cell_lines:
        quoted = ", ".join(sql_quote(value) for value in target_cell_lines)
        lines.append(f"DELETE FROM `{study_meta_table}` WHERE `cell_line` IN ({quoted});")
        lines.append(f"DELETE FROM `{cell_line_meta_table}` WHERE `cell_line` IN ({quoted});")
        lines.append("")

    for slug, shard_path in shard_files.items():
        lines.append(f"TRUNCATE TABLE `{result_table_prefix}{slug}`;")
        lines.append(
            build_load_data_sql(
                shard_path,
                f"{result_table_prefix}{slug}",
                SHARD_OUTPUT_COLUMNS,
                output_dir,
                sql_data_dir,
            )
        )
        lines.append("")

    meta_dir = output_dir / "meta"
    lines.append(
        build_load_data_sql(
            meta_dir / "comparison_cell_line_meta.tsv",
            cell_line_meta_table,
            CELL_LINE_META_COLUMNS,
            output_dir,
            sql_data_dir,
        )
    )
    lines.append("")
    lines.append(
        build_load_data_sql(
            meta_dir / "comparison_study_meta.tsv",
            study_meta_table,
            STUDY_META_COLUMNS,
            output_dir,
            sql_data_dir,
        )
    )
    lines.append("")
    return "\n".join(lines).strip() + "\n"


def build_load_data_sql(
    path: Path,
    table_name: str,
    columns: Sequence[str],
    output_dir: Path,
    sql_data_dir: Path | None,
) -> str:
    column_sql = ",\n  ".join(f"`{column}`" for column in columns)
    file_path = resolve_sql_data_path(path, output_dir, sql_data_dir).replace("'", "''")
    return (
        f"LOAD DATA LOCAL INFILE '{file_path}'\n"
        f"INTO TABLE `{table_name}`\n"
        "CHARACTER SET utf8mb4\n"
        "FIELDS TERMINATED BY '\\t'\n"
        "OPTIONALLY ENCLOSED BY '\"'\n"
        "LINES TERMINATED BY '\\n'\n"
        "IGNORE 1 LINES\n"
        "(\n"
        f"  {column_sql}\n"
        ");"
    )


def resolve_sql_data_path(path: Path, output_dir: Path, sql_data_dir: Path | None) -> str:
    if sql_data_dir is None:
        return path.resolve().as_posix()
    try:
        relative_path = path.resolve().relative_to(output_dir.resolve())
    except ValueError as exc:
        raise SystemExit(f"Generated TSV path is not under --output-dir: {path}") from exc
    return (sql_data_dir / relative_path).as_posix()


def build_validate_sql(
    database: str,
    result_table_prefix: str,
    cell_line_meta_table: str,
    study_meta_table: str,
    cell_line_rows: List[dict[str, str]],
) -> str:
    lines = [f"USE `{database}`;", ""]
    for row in cell_line_rows:
        lines.append(
            f"SELECT '{row['cell_line']}' AS cell_line, COUNT(*) AS row_count "
            f"FROM `{result_table_prefix}{row['cell_line_slug']}`;"
        )
    if cell_line_rows:
        quoted = ", ".join(sql_quote(row["cell_line"]) for row in cell_line_rows)
        lines.append("")
        lines.append(
            f"SELECT `cell_line`, `study_count`, `comparison_count`, `gene_row_count`\n"
            f"FROM `{cell_line_meta_table}`\n"
            f"WHERE `cell_line` IN ({quoted})\n"
            "ORDER BY `sort_order`, `cell_line`;"
        )
        lines.append("")
        lines.append(
            f"SELECT `cell_line`, `bioproject_accession`, `comparison`, `gene_row_count`\n"
            f"FROM `{study_meta_table}`\n"
            f"WHERE `cell_line` IN ({quoted})\n"
            "ORDER BY `cell_line`, `sort_order`, `bioproject_accession`, `comparison`;"
        )
    return "\n".join(lines).strip() + "\n"


def sql_quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def main() -> None:
    args = parse_args()
    metadata = read_cell_line_metadata(args.cell_line_metadata)
    gene_alias_map = read_gene_alias_map(args.gene_alias_map)
    target_info, target_dir = resolve_cell_line_target(args, metadata)

    bundles_by_slug: Dict[str, List[ComparisonBundle]] = {}
    manifest_rows: List[dict[str, str]] = []

    comparison_dirs = find_comparison_dirs(target_dir, args.comparison_dir_glob)
    if not comparison_dirs:
        raise SystemExit(f"No valid comparison directories found under: {target_dir}")

    bundles: List[ComparisonBundle] = []
    for comparison_dir in comparison_dirs:
        bundle = build_comparison_rows(target_info, comparison_dir, gene_alias_map, args.padj_cutoff)
        bundles.append(bundle)
        manifest_rows.append(
            {
                "cell_line": target_info.display_name,
                "cell_line_slug": target_info.slug,
                "bioproject_accession": bundle.bioproject_accession,
                "comparison": bundle.comparison,
                "comparison_dir": comparison_dir.resolve().as_posix(),
                "control_sample_count": str(len(bundle.control_samples)),
                "experimental_sample_count": str(len(bundle.experimental_samples)),
                "gene_row_count": str(len(bundle.gene_rows)),
            }
        )
    bundles_by_slug[target_info.slug] = bundles

    if not bundles_by_slug:
        raise SystemExit("No valid comparison directories were found for the selected cell line.")

    shard_dir = args.output_dir / "shards"
    meta_dir = args.output_dir / "meta"
    sql_dir = args.output_dir / "sql"
    ensure_output_dir(shard_dir)
    ensure_output_dir(meta_dir)
    ensure_output_dir(sql_dir)

    cell_line_rows: List[dict[str, str]] = []
    study_rows: List[dict[str, str]] = []
    shard_files: Dict[str, Path] = {}

    for slug, bundles in sorted(bundles_by_slug.items()):
        info = bundles[0].cell_line
        shard_path = shard_dir / f"comparison_gene_exp.{slug}.tsv"
        shard_rows: List[dict[str, str]] = []
        for bundle in bundles:
            shard_rows.extend(bundle.gene_rows)
        write_tsv(shard_path, SHARD_OUTPUT_COLUMNS, shard_rows)
        shard_files[slug] = shard_path

        cell_line_rows.append(
            {
                "cell_line": info.display_name,
                "cell_line_slug": info.slug,
                "study_count": str(len({bundle.bioproject_accession for bundle in bundles})),
                "comparison_count": str(len(bundles)),
                "gene_row_count": str(sum(len(bundle.gene_rows) for bundle in bundles)),
                "sort_order": str(info.sort_order),
            }
        )
        sorted_bundles = sorted(bundles, key=lambda item: (item.bioproject_accession, item.comparison))
        for index, bundle in enumerate(sorted_bundles):
            study_rows.append(
                {
                    "cell_line": info.display_name,
                    "cell_line_slug": info.slug,
                    "bioproject_accession": bundle.bioproject_accession,
                    "comparison": bundle.comparison,
                    "gene_row_count": str(len(bundle.gene_rows)),
                    "sort_order": str(index),
                }
            )

    cell_line_rows.sort(key=lambda item: (int(item["sort_order"]), item["cell_line"]))
    study_rows.sort(key=lambda item: (item["cell_line"], int(item["sort_order"]), item["bioproject_accession"], item["comparison"]))
    manifest_rows.sort(key=lambda item: (item["cell_line"], item["bioproject_accession"], item["comparison"]))

    write_tsv(meta_dir / "comparison_cell_line_meta.tsv", CELL_LINE_META_COLUMNS, cell_line_rows)
    write_tsv(meta_dir / "comparison_study_meta.tsv", STUDY_META_COLUMNS, study_rows)
    write_tsv(args.output_dir / "comparison_import_manifest.tsv", MANIFEST_COLUMNS, manifest_rows)

    load_sql = build_load_sql(
        output_dir=args.output_dir,
        sql_data_dir=args.sql_data_dir,
        database=args.database,
        result_table_prefix=args.result_table_prefix,
        cell_line_meta_table=args.cell_line_meta_table,
        study_meta_table=args.study_meta_table,
        shard_files=shard_files,
        cell_line_rows=cell_line_rows,
    )
    (sql_dir / "load_ds0002_transcriptome_import.sql").write_text(load_sql, encoding="utf-8")

    validate_sql = build_validate_sql(
        database=args.database,
        result_table_prefix=args.result_table_prefix,
        cell_line_meta_table=args.cell_line_meta_table,
        study_meta_table=args.study_meta_table,
        cell_line_rows=cell_line_rows,
    )
    (sql_dir / "validate_ds0002_transcriptome_import.sql").write_text(validate_sql, encoding="utf-8")

    print(f"Prepared {len(shard_files)} shard TSV file(s) in {shard_dir}")
    print(f"Prepared metadata TSV files in {meta_dir}")
    print(f"Prepared SQL files in {sql_dir}")


if __name__ == "__main__":
    main()
