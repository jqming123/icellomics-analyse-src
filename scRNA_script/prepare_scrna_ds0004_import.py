#!/usr/bin/env python3
"""Prepare DS0004 import assets from a cluster-centric scRNA-seq h5ad file.

This script targets the standalone single-cell database:
  big_industry_cell_singlecell_ds0004

It writes import-ready TSV files for the unified tables:
  - scrna_dataset
  - scrna_cell
  - scrna_umap
  - scrna_gene
  - scrna_gene_exp
  - scrna_cluster_top_gene
  - scrna_sample_cluster_stat
  - scrna_cell_signature_score  (optional; only when obs holds UCell__/AUCell__ score columns)

It also emits a LOAD DATA SQL template and a JSON manifest.
The schema assumes the dataset is primarily organized by clusters rather than
stable cell-type annotations.

The generated LOAD DATA shell script expects the TSV files to be present on the
database server. Use --sql-import-root to set that server-side directory. It is
only written into the generated shell script as DATA_DIR; it does not control
the local --output-dir. By default, DATA_DIR is:
  /data/industry_cellline_data/scrna_data/<cell_line>/<project_accession>
If project_accession is unavailable, dataset_key is used for the final path
segment so the generated path does not contain a null value.

Usage example (Windows PowerShell)::

  python `
    CellLine-back-master/scripts/prepare_scrna_ds0004_import.py `
    --input-h5ad ref_data/seurat_obj_annotated.h5ad `
    --output-dir test_data/scrna_ds0004 `
    --singlecell-dataset-id 400001 `
    --cell-line HEK293

Usage example (Linux)::

  python3 prepare_scrna_ds0004_import.py \
    --input-h5ad ref_data/seurat_obj_annotated.h5ad \
    --output-dir test_data/scrna_ds0004 \
    --singlecell-dataset-id 400001 \
    --cell-line HEK293
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


DEFAULT_DATABASE_NAME = "big_industry_cell_singlecell_ds0004"
DEFAULT_UMAP_KEY = "X_umap"
DEFAULT_DELIMITER_NAME = "tsv"
DEFAULT_DATASET_TABLE = "scrna_dataset"
DEFAULT_CELL_TABLE = "scrna_cell"
DEFAULT_UMAP_TABLE = "scrna_umap"
DEFAULT_GENE_TABLE = "scrna_gene"
DEFAULT_GENE_EXP_TABLE = "scrna_gene_exp"
DEFAULT_STAT_TABLE = "scrna_sample_cluster_stat"
DEFAULT_CLUSTER_TOP_GENE_TABLE = "scrna_cluster_top_gene"
DEFAULT_CELL_SIGNATURE_SCORE_TABLE = "scrna_cell_signature_score"
DEFAULT_TABLE_PREFIX = "scrna"
DEFAULT_SOURCE_CODE = "DS0004"
DEFAULT_DATASET_FIELD = "dataset"
DEFAULT_SAMPLE_FIELD = "sample_id"
DEFAULT_CLUSTER_FIELD = "seurat_clusters"
DEFAULT_ORIG_IDENT_FIELD = "orig.ident"
DEFAULT_CELL_ANNOTATION_FIELD = "cell_type"
DEFAULT_CONFIDENCE_FIELD = "annotation_confidence"
DEFAULT_NCOUNT_FIELD = "nCount_RNA"
DEFAULT_NFEATURE_FIELD = "nFeature_RNA"
DEFAULT_PERCENT_MT_FIELD = "percent.mt"
DEFAULT_PERCENT_RIBO_FIELD = "percent.ribo"
DEFAULT_DOUBLET_SCORE_FIELD = "doublet_score"
DEFAULT_PREDICTED_DOUBLET_FIELD = "predicted_doublet"
DEFAULT_DOUBLET_THRESHOLD_FIELD = "doublet_threshold"
DEFAULT_CELL_CYCLE_S_FIELD = "S.Score"
DEFAULT_CELL_CYCLE_G2M_FIELD = "G2M.Score"
DEFAULT_CELL_CYCLE_PHASE_FIELD = "Phase"
DEFAULT_SIGNATURE_SCORE_PREFIXES = "UCell__,AUCell__"
DEFAULT_CELL_LINE_TSV = "docs/cell_line_metadata.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare DS0004 scRNA-seq import files from an h5ad file."
    )
    parser.add_argument("--input-h5ad", type=Path, required=True, help="Input h5ad file.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for generated TSV/SQL files.")
    parser.add_argument(
        "--singlecell-dataset-id",
        type=int,
        required=True,
        help="Stable numeric singlecell_dataset_id to be used across all six tables.",
    )
    parser.add_argument(
        "--dataset-field",
        default=DEFAULT_DATASET_FIELD,
        help="obs column holding the dataset label. Default: dataset.",
    )
    parser.add_argument("--dataset-key", default=None, help="Override dataset_key. Defaults to the single unique dataset field value or h5ad stem.")
    parser.add_argument("--dataset-name", default=None, help="Override dataset_name. Defaults to dataset_key.")
    parser.add_argument("--dataset-label", default=None, help="Optional human-readable dataset label.")
    parser.add_argument("--project-accession", default=None, help="Optional project accession override, for example PRJNA1256660.")
    parser.add_argument("--cell-line", type=str, required=True, help="Cell line name, e.g. HEK293. Must match a cell_line_name in the cell-line metadata TSV.")
    parser.add_argument(
        "--cell-line-tsv",
        type=Path,
        default=DEFAULT_CELL_LINE_TSV,
        help="Path to cell line metadata TSV (cell_line_name, cell_line_slug, species_name, reference_genome). Default: docs/cell_line_metadata.tsv.",
    )
    parser.add_argument("--umap-key", default=DEFAULT_UMAP_KEY, help="obsm key for UMAP coordinates. Default: X_umap.")
    parser.add_argument("--sample-field", default=DEFAULT_SAMPLE_FIELD, help="obs column for sample identifiers. Default: sample_id.")
    parser.add_argument("--cluster-field", default=DEFAULT_CLUSTER_FIELD, help="obs column for clusters. Default: seurat_clusters.")
    parser.add_argument("--condition-field", default=None, help="Optional obs column for experimental condition or group.")
    parser.add_argument("--orig-ident-field", default=DEFAULT_ORIG_IDENT_FIELD, help="Optional obs column for original sample identity. Default: orig.ident.")
    parser.add_argument(
        "--cell-annotation-field",
        default=DEFAULT_CELL_ANNOTATION_FIELD,
        help="Optional obs column for supplemental cell annotation. Default: cell_type.",
    )
    parser.add_argument(
        "--confidence-field",
        default=DEFAULT_CONFIDENCE_FIELD,
        help="Optional obs column for annotation confidence. Default: annotation_confidence.",
    )
    parser.add_argument("--ncount-field", default=DEFAULT_NCOUNT_FIELD, help="Optional obs column for total RNA count.")
    parser.add_argument("--nfeature-field", default=DEFAULT_NFEATURE_FIELD, help="Optional obs column for feature count.")
    parser.add_argument("--percent-mt-field", default=DEFAULT_PERCENT_MT_FIELD, help="Optional obs column for mitochondrial percentage.")
    parser.add_argument("--percent-ribo-field", default=DEFAULT_PERCENT_RIBO_FIELD, help="Optional obs column for ribosomal percentage.")
    parser.add_argument("--doublet-score-field", default=DEFAULT_DOUBLET_SCORE_FIELD, help="Optional obs column for doublet score.")
    parser.add_argument(
        "--predicted-doublet-field",
        default=DEFAULT_PREDICTED_DOUBLET_FIELD,
        help="Optional obs column for predicted doublet flag.",
    )
    parser.add_argument(
        "--doublet-threshold-field",
        default=DEFAULT_DOUBLET_THRESHOLD_FIELD,
        help="Optional obs column for doublet threshold.",
    )
    parser.add_argument(
        "--cell-cycle-s-field",
        default=DEFAULT_CELL_CYCLE_S_FIELD,
        help="Optional obs column for the S-phase score written to scrna_cell. Default: S.Score.",
    )
    parser.add_argument(
        "--cell-cycle-g2m-field",
        default=DEFAULT_CELL_CYCLE_G2M_FIELD,
        help="Optional obs column for the G2M-phase score written to scrna_cell. Default: G2M.Score.",
    )
    parser.add_argument(
        "--cell-cycle-phase-field",
        default=DEFAULT_CELL_CYCLE_PHASE_FIELD,
        help="Optional obs column for the cell-cycle phase written to scrna_cell. Default: Phase.",
    )
    parser.add_argument(
        "--signature-score-prefixes",
        default=DEFAULT_SIGNATURE_SCORE_PREFIXES,
        help=(
            "Comma-separated obs column prefixes exported to the cell signature score long table. "
            "score_type is the prefix without a trailing '__'. Default: UCell__,AUCell__."
        ),
    )
    parser.add_argument(
        "--cluster-label-prefix",
        default="Cluster ",
        help="Prefix used when constructing cluster_label from cluster_id. Default: 'Cluster '.",
    )
    parser.add_argument(
        "--gene-symbol-source",
        choices=("auto", "index", "gene_name"),
        default="auto",
        help="How to derive gene_symbol for scrna_gene. Default: auto.",
    )
    parser.add_argument(
        "--gene-name-field",
        default="gene_name",
        help="var column used as gene_name and as gene_symbol when --gene-symbol-source gene_name. Default: gene_name.",
    )
    parser.add_argument(
        "--raw-mode",
        choices=("auto", "include", "skip"),
        default="auto",
        help="Whether to export raw.X into scrna_gene_exp as value_type=raw_count. Default: auto.",
    )
    parser.add_argument(
        "--delimiter",
        choices=("tsv", "csv"),
        default=DEFAULT_DELIMITER_NAME,
        help="Output delimiter format. Default: tsv.",
    )
    parser.add_argument(
        "--gene-exp-progress-interval",
        type=int,
        default=1000,
        help="Progress interval while writing scrna_gene_exp rows. Default: 1000 genes.",
    )
    parser.add_argument(
        "--cluster-top-gene-limit",
        type=int,
        default=2000,
        help="Number of top expressed genes to export per cluster. Default: 2000.",
    )
    parser.add_argument(
        "--mysql-user",
        default="root",
        help="MySQL username for the generated import bash script.",
    )
    parser.add_argument(
        "--mysql-password",
        default="RootMySQL@Song123!",
        help="MySQL password for the generated import bash script. If not set, the script uses -p without password (expects ~/.my.cnf or MYSQL_PWD).",
    )
    parser.add_argument(
        "--mysql-socket",
        default="/data/mysql3306/data/mysql.sock",
        help="MySQL socket path for the generated import bash script (--socket=...).",
    )
    parser.add_argument(
        "--sql-import-root",
        default=None,
        help=(
            "Server-side directory where TSV files will be placed for the generated LOAD DATA shell script. "
            "Defaults to '/data/industry_cellline_data/scrna_data/<cell_line>/<project_accession>'. "
            "If project_accession is unavailable, dataset_key is used as the final path segment."
        ),
    )
    parser.add_argument(
        "--table-prefix",
        default=DEFAULT_TABLE_PREFIX,
        help=(
            "Target table prefix for generated TSV/import files. "
            "Default: scrna, which targets the public scrna_* tables. "
            "Example: scrna_starCL2 targets scrna_starCL2_* tables."
        ),
    )
    parser.add_argument(
        "--datasource-code",
        default=DEFAULT_SOURCE_CODE,
        help="datasource_code value written to the dataset row. Default: DS0004.",
    )
    parser.add_argument(
        "--import-script-stem",
        default=None,
        help="Optional generated import shell script stem. Defaults to load_scrna_ds0004 for scrna, otherwise load_<table-prefix>.",
    )
    parser.add_argument(
        "--manifest-stem",
        default=None,
        help="Optional generated manifest JSON stem. Defaults to scrna_manifest for scrna, otherwise <table-prefix>_manifest.",
    )
    return parser.parse_args()


def slugify(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "-", value.strip())
    text = text.strip("-").lower()
    return text or "dataset"


def configure_table_names(table_prefix: str) -> None:
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", table_prefix):
        raise ValueError(
            "--table-prefix must start with a letter and contain only letters, numbers, and underscores."
        )

    global DEFAULT_DATASET_TABLE
    global DEFAULT_CELL_TABLE
    global DEFAULT_UMAP_TABLE
    global DEFAULT_GENE_TABLE
    global DEFAULT_GENE_EXP_TABLE
    global DEFAULT_STAT_TABLE
    global DEFAULT_CLUSTER_TOP_GENE_TABLE
    global DEFAULT_CELL_SIGNATURE_SCORE_TABLE

    DEFAULT_DATASET_TABLE = f"{table_prefix}_dataset"
    DEFAULT_CELL_TABLE = f"{table_prefix}_cell"
    DEFAULT_UMAP_TABLE = f"{table_prefix}_umap"
    DEFAULT_GENE_TABLE = f"{table_prefix}_gene"
    DEFAULT_GENE_EXP_TABLE = f"{table_prefix}_gene_exp"
    DEFAULT_STAT_TABLE = f"{table_prefix}_sample_cluster_stat"
    DEFAULT_CLUSTER_TOP_GENE_TABLE = f"{table_prefix}_cluster_top_gene"
    DEFAULT_CELL_SIGNATURE_SCORE_TABLE = f"{table_prefix}_cell_signature_score"


def default_import_script_stem(table_prefix: str) -> str:
    if table_prefix == DEFAULT_TABLE_PREFIX:
        return "load_scrna_ds0004"
    return f"load_{table_prefix}"


def default_manifest_stem(table_prefix: str) -> str:
    if table_prefix == DEFAULT_TABLE_PREFIX:
        return "scrna_manifest"
    return f"{table_prefix}_manifest"


def default_sql_import_root(cell_line: str, project_accession: str | None, dataset_key: str) -> str:
    project_segment = project_accession or dataset_key
    return f"/data/industry_cellline_data/scrna_data/{cell_line}/{project_segment}"


def compute_sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def parse_project_and_cell_line(dataset_value: str) -> tuple[str | None, str | None]:
    match = re.match(r"^((?:PRJ|SRP)[^_]+)_(.+)$", dataset_value)
    if not match:
        return None, None
    return match.group(1), match.group(2)


def load_cell_line_metadata(tsv_path: Path) -> dict[str, dict[str, str | None]]:
    """Load cell line metadata TSV, keyed by cell_line_name."""
    if not tsv_path.exists():
        raise FileNotFoundError(f"Cell line metadata TSV not found: {tsv_path}")
    df = pd.read_csv(str(tsv_path), sep="\t", dtype=str, keep_default_na=False)
    required_cols = {"cell_line_name", "cell_line_slug", "species_name", "reference_genome"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Cell line TSV missing columns: {', '.join(sorted(missing))}")
    metadata: dict[str, dict[str, str | None]] = {}
    for _, row in df.iterrows():
        name = row["cell_line_name"].strip()
        if not name:
            continue
        metadata[name] = {
            "cell_line_slug": to_nullable_string(row.get("cell_line_slug")),
            "species_name": to_nullable_string(row.get("species_name")),
            "reference_genome": to_nullable_string(row.get("reference_genome")),
        }
    return metadata


def require_column(df: pd.DataFrame, column: str) -> str:
    if column in df.columns:
        return column
    raise KeyError(f"Required column not found in obs: {column}")


def find_column(df: pd.DataFrame, column: str | None) -> str | None:
    if column is None:
        return None
    if column in df.columns:
        return column
    return None


def get_single_unique_value(df: pd.DataFrame, column: str | None) -> str | None:
    if column is None:
        return None
    values = df[column].dropna().astype(str).unique().tolist()
    if not values:
        return None
    if len(values) > 1:
        raise ValueError(
            f"Expected one unique value in obs['{column}'], but found {len(values)}: {values[:10]}"
        )
    return values[0]


def to_nullable_string(value: Any) -> str | None:
    if value is None:
        return None
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "":
        return None
    return text


def to_nullable_bool_int(value: Any) -> int | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, bool):
        return int(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return 1
    if text in {"false", "f", "0", "no", "n"}:
        return 0
    return None


def normalize_optional_string_series(series: pd.Series) -> pd.Series:
    return series.map(to_nullable_string)


def normalize_required_string_series(series: pd.Series, column_name: str) -> pd.Series:
    normalized = normalize_optional_string_series(series)
    missing_count = int(normalized.isna().sum())
    if missing_count > 0:
        raise ValueError(
            f"Required column contains {missing_count} empty/null values: {column_name}"
        )
    return normalized.astype(str)


def build_cluster_order_map(cluster_values: Iterable[str]) -> dict[str, int]:
    unique_values = sorted(set(cluster_values), key=cluster_sort_key)
    order_map: dict[str, int] = {}
    for index, value in enumerate(unique_values):
        if re.fullmatch(r"\d+", value):
            order_map[value] = int(value)
        else:
            order_map[value] = index
    return order_map


def cluster_sort_key(value: str) -> tuple[int, Any]:
    if re.fullmatch(r"\d+", value):
        return (0, int(value))
    return (1, value)


def materialize_sparse_matrix(matrix: Any) -> sparse.csr_matrix:
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if sparse.isspmatrix_csr(matrix):
        return matrix
    if sparse.isspmatrix(matrix):
        return matrix.tocsr()
    return sparse.csr_matrix(matrix)


def resolve_gene_symbol_series(var_df: pd.DataFrame, gene_symbol_source: str, gene_name_field: str) -> pd.Series:
    if gene_symbol_source == "index":
        return pd.Series(var_df.index.astype(str), index=var_df.index)
    if gene_symbol_source == "gene_name":
        if gene_name_field not in var_df.columns:
            raise KeyError(f"gene_name field not found in var: {gene_name_field}")
        return normalize_required_string_series(var_df[gene_name_field], gene_name_field)
    if gene_name_field in var_df.columns:
        gene_name_series = normalize_optional_string_series(var_df[gene_name_field])
        if gene_name_series.notna().all():
            return gene_name_series.astype(str)
    return pd.Series(var_df.index.astype(str), index=var_df.index)


def sanitize_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value)


def write_dict_rows(path: Path, columns: list[str], rows: Iterable[dict[str, Any]], delimiter: str) -> int:
    row_count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=delimiter,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_export_value(row.get(column)) for column in columns})
            row_count += 1
    return row_count


def format_export_value(value: Any) -> Any:
    if value is None:
        return r"\N"
    try:
        if pd.isna(value):
            return r"\N"
    except TypeError:
        pass
    if isinstance(value, bool):
        return int(value)
    return value


def truncate_nullable(value: Any, max_length: int) -> str | None:
    text = to_nullable_string(value)
    if text is None:
        return None
    return text[:max_length]


def dense_1d(value: Any) -> np.ndarray:
    array = np.asarray(value).ravel()
    return array.astype(float, copy=False)


def generate_cluster_top_gene_rows(
    matrix: sparse.csr_matrix,
    gene_rows: pd.DataFrame,
    obs_export: pd.DataFrame,
    singlecell_dataset_id: int,
    value_type: str,
    limit: int,
) -> Iterable[dict[str, Any]]:
    safe_limit = max(1, int(limit))
    clusters = (
        obs_export[["cluster_id", "cluster_label", "cluster_order"]]
        .drop_duplicates()
        .sort_values(["cluster_order", "cluster_id"], kind="stable")
    )

    for _, cluster in clusters.iterrows():
        cluster_id = str(cluster["cluster_id"])
        cell_indices = obs_export.loc[obs_export["cluster_id"] == cluster_id, "cell_index"].astype(int).to_numpy()
        cluster_cell_count = int(len(cell_indices))
        if cluster_cell_count == 0:
            continue

        print(f"[cluster-top-gene] {value_type}: cluster {cluster_id}, cells={cluster_cell_count}")
        cluster_matrix = matrix[cell_indices, :]
        sum_values = dense_1d(cluster_matrix.sum(axis=0))
        nonzero_counts = np.asarray(cluster_matrix.getnnz(axis=0)).ravel().astype(int)
        mean_values = sum_values / cluster_cell_count

        top_df = pd.DataFrame(
            {
                "gene_index": np.arange(matrix.shape[1], dtype=int),
                "mean_expression": mean_values,
                "sum_expression": sum_values,
                "nonzero_cell_count": nonzero_counts,
            }
        )
        top_df = top_df.merge(
            gene_rows[["gene_index", "gene_symbol", "gene_name"]],
            on="gene_index",
            how="left",
        )
        top_df["gene_symbol_sort"] = top_df["gene_symbol"].fillna("").astype(str)
        top_df = top_df.sort_values(
            ["mean_expression", "nonzero_cell_count", "gene_symbol_sort"],
            ascending=[False, False, True],
            kind="stable",
        ).head(safe_limit)

        for rank_no, (_, row) in enumerate(top_df.iterrows(), start=1):
            nonzero_cell_count = int(row["nonzero_cell_count"])
            yield {
                "singlecell_dataset_id": singlecell_dataset_id,
                "cluster_id": truncate_nullable(cluster_id, 50),
                "cluster_label": truncate_nullable(cluster["cluster_label"], 50),
                "cluster_order": cluster["cluster_order"],
                "gene_index": int(row["gene_index"]),
                "gene_symbol": truncate_nullable(row["gene_symbol"], 50) or "",
                "gene_name": truncate_nullable(row["gene_name"], 50),
                "value_type": value_type,
                "rank_no": rank_no,
                "mean_expression": float(row["mean_expression"]),
                "sum_expression": float(row["sum_expression"]),
                "nonzero_cell_count": nonzero_cell_count,
                "cluster_cell_count": cluster_cell_count,
                "expressing_cell_ratio": nonzero_cell_count / cluster_cell_count,
            }


def write_cluster_top_gene_file(
    path: Path,
    singlecell_dataset_id: int,
    matrix: sparse.csr_matrix,
    gene_rows: pd.DataFrame,
    obs_export: pd.DataFrame,
    value_type: str,
    limit: int,
    delimiter: str,
) -> int:
    columns = [
        "singlecell_dataset_id",
        "cluster_id",
        "cluster_label",
        "cluster_order",
        "gene_index",
        "gene_symbol",
        "gene_name",
        "value_type",
        "rank_no",
        "mean_expression",
        "sum_expression",
        "nonzero_cell_count",
        "cluster_cell_count",
        "expressing_cell_ratio",
    ]
    return write_dict_rows(
        path,
        columns,
        generate_cluster_top_gene_rows(
            matrix,
            gene_rows,
            obs_export,
            singlecell_dataset_id,
            value_type,
            limit,
        ),
        delimiter,
    )


def generate_gene_exp_rows(
    matrix: sparse.csr_matrix,
    singlecell_dataset_id: int,
    value_type: str,
    progress_interval: int,
) -> Iterable[dict[str, Any]]:
    csc = matrix.tocsc()
    total_genes = csc.shape[1]
    for gene_index in range(total_genes):
        start = csc.indptr[gene_index]
        end = csc.indptr[gene_index + 1]
        cell_indices = csc.indices[start:end]
        values = csc.data[start:end]
        if gene_index % progress_interval == 0:
            print(f"[gene-exp] {value_type}: processed {gene_index}/{total_genes} genes")
        if len(cell_indices) == 0:
            yield {
                "singlecell_dataset_id": singlecell_dataset_id,
                "gene_index": gene_index,
                "value_type": value_type,
                "nonzero_cell_count": 0,
                "cell_index_list": "",
                "exp_value_list": "",
                "min_value": 0.0,
                "max_value": 0.0,
            }
            continue
        cell_index_list = "|".join(str(int(idx)) for idx in cell_indices.tolist())
        exp_value_list = "|".join(format(float(value), ".10g") for value in values.tolist())
        yield {
            "singlecell_dataset_id": singlecell_dataset_id,
            "gene_index": gene_index,
            "value_type": value_type,
            "nonzero_cell_count": int(len(cell_indices)),
            "cell_index_list": cell_index_list,
            "exp_value_list": exp_value_list,
            "min_value": float(values.min()),
            "max_value": float(values.max()),
        }


def write_gene_exp_file(
    path: Path,
    singlecell_dataset_id: int,
    matrix: sparse.csr_matrix,
    value_type: str,
    delimiter: str,
    progress_interval: int,
) -> int:
    columns = [
        "singlecell_dataset_id",
        "gene_index",
        "value_type",
        "nonzero_cell_count",
        "cell_index_list",
        "exp_value_list",
        "min_value",
        "max_value",
    ]
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter=delimiter,
            extrasaction="ignore",
            lineterminator="\n",
        )
        row_count = 0
        if handle.tell() == 0:
            writer.writeheader()
        for row in generate_gene_exp_rows(matrix, singlecell_dataset_id, value_type, progress_interval):
            writer.writerow({column: format_export_value(row.get(column)) for column in columns})
            row_count += 1
    return row_count


def write_load_data_bash(
    output_path: Path,
    database_name: str,
    delimiter_name: str,
    file_map: dict[str, Path],
    singlecell_dataset_id: int,
    mysql_user: str | None,
    mysql_password: str | None,
    mysql_socket: str | None,
    sql_import_root: str | None,
) -> None:
    delimiter = r"\t" if delimiter_name == "tsv" else ","
    table_columns = {
        DEFAULT_DATASET_TABLE: [
            "singlecell_dataset_id",
            "datasource_code",
            "dataset_key",
            "dataset_name",
            "dataset_label",
            "project_accession",
            "cell_line",
            "cell_line_slug",
            "species_name",
            "reference_genome",
            "sample_count",
            "cell_count",
            "gene_count",
            "cluster_count",
            "umap_key",
            "default_cluster_field",
            "default_sample_field",
            "expression_value_type",
            "raw_value_type",
            "has_raw",
            "source_h5ad_path",
            "source_file_hash",
            "notes",
        ],
        DEFAULT_CELL_TABLE: [
            "singlecell_dataset_id",
            "cell_index",
            "barcode",
            "sample_id",
            "orig_ident",
            "condition_label",
            "cluster_id",
            "cluster_label",
            "cluster_order",
            "cell_annotation",
            "annotation_confidence",
            "predicted_doublet",
            "doublet_score",
            "doublet_threshold",
            "ncount_rna",
            "nfeature_rna",
            "percent_mt",
            "percent_ribo",
            "s_score",
            "g2m_score",
            "cell_cycle_phase",
        ],
        DEFAULT_UMAP_TABLE: [
            "singlecell_dataset_id",
            "cell_index",
            "umap_key",
            "umap_1",
            "umap_2",
        ],
        DEFAULT_CELL_SIGNATURE_SCORE_TABLE: [
            "singlecell_dataset_id",
            "cell_index",
            "score_type",
            "signature_name",
            "score_value",
        ],
        DEFAULT_GENE_TABLE: [
            "singlecell_dataset_id",
            "gene_index",
            "gene_symbol",
            "gene_name",
            "highly_variable",
        ],
        DEFAULT_GENE_EXP_TABLE: [
            "singlecell_dataset_id",
            "gene_index",
            "value_type",
            "nonzero_cell_count",
            "cell_index_list",
            "exp_value_list",
            "min_value",
            "max_value",
        ],
        DEFAULT_STAT_TABLE: [
            "singlecell_dataset_id",
            "sample_id",
            "condition_label",
            "cluster_id",
            "cluster_label",
            "cluster_order",
            "cell_count",
            "sample_cell_count",
            "proportion",
        ],
        DEFAULT_CLUSTER_TOP_GENE_TABLE: [
            "singlecell_dataset_id",
            "cluster_id",
            "cluster_label",
            "cluster_order",
            "gene_index",
            "gene_symbol",
            "gene_name",
            "value_type",
            "rank_no",
            "mean_expression",
            "sum_expression",
            "nonzero_cell_count",
            "cluster_cell_count",
            "expressing_cell_ratio",
        ],
    }

    lines: list[str] = [
        "#!/bin/bash",
        "set -euo pipefail",
        "",
        f"MYSQL_USER='{mysql_user or '<MYSQL_USER>'}';",
    ]
    if mysql_password:
        lines.append(f"MYSQL_PASSWORD='{mysql_password}';")
    lines.extend([
        f"MYSQL_SOCKET='{mysql_socket or '<MYSQL_SOCKET>'}';",
        f"DB_NAME='{database_name}';",
        f"DATA_DIR='{sql_import_root or '<DATA_DIR>'}';",
        "",
        f"echo '=== Starting scRNA DS0004 import (dataset_id={singlecell_dataset_id}) ==='",
        "",
    ])

    import_order = [
        DEFAULT_DATASET_TABLE,
        DEFAULT_GENE_TABLE,
        DEFAULT_CELL_TABLE,
        DEFAULT_UMAP_TABLE,
        DEFAULT_CELL_SIGNATURE_SCORE_TABLE,
        DEFAULT_GENE_EXP_TABLE,
        DEFAULT_CLUSTER_TOP_GENE_TABLE,
        DEFAULT_STAT_TABLE,
    ]
    import_order = [table_name for table_name in import_order if table_name in file_map]

    password_arg = '-p"${MYSQL_PASSWORD}"' if mysql_password else "-p"
    socket_arg = "--socket=${MYSQL_SOCKET}" if mysql_socket else ""

    for table_name in import_order:
        file_path = file_map[table_name]
        base_name = file_path.name
        columns = table_columns[table_name]
        column_list = ", ".join(columns)

        lines.extend([
            f"echo 'Importing {table_name}...'",
            f"nohup mysql -u\"${{MYSQL_USER}}\" {password_arg} \\",
        ])
        if mysql_socket:
            lines.append(f"  {socket_arg} \\")
        lines.extend([
            f"  --local-infile=1 \\",
            f"  -D \"${{DB_NAME}}\" \\",
            f"  -e \"LOAD DATA LOCAL INFILE '${{DATA_DIR}}/{base_name}'",
            f"      INTO TABLE {table_name}",
            f"      CHARACTER SET utf8mb4",
            f"      FIELDS TERMINATED BY '{delimiter}'",
            f"      ENCLOSED BY '\\\"'",
            f"      LINES TERMINATED BY '\\\\n'",
            f"      IGNORE 1 LINES",
            f"      ({column_list});\" \\",
            f"  > {table_name}.import.out 2>&1 &",
            "",
        ])

    lines.extend([
        "wait",
        "echo '=== All imports completed. ==='",
        "",
    ])

    output_path.write_text("\n".join(lines), encoding="utf-8")
    output_path.chmod(0o755)


def main() -> None:
    args = parse_args()
    configure_table_names(args.table_prefix)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.input_h5ad.exists():
        raise FileNotFoundError(f"h5ad file not found: {args.input_h5ad}")

    print(f"[load] reading {args.input_h5ad}")
    adata = ad.read_h5ad(args.input_h5ad)
    obs = adata.obs.copy()
    var_df = adata.var.copy()

    sample_field = require_column(obs, args.sample_field)
    cluster_field = require_column(obs, args.cluster_field)
    dataset_field = find_column(obs, args.dataset_field)
    orig_ident_field = find_column(obs, args.orig_ident_field)
    condition_field = find_column(obs, args.condition_field)
    cell_annotation_field = find_column(obs, args.cell_annotation_field)
    confidence_field = find_column(obs, args.confidence_field)
    ncount_field = find_column(obs, args.ncount_field)
    nfeature_field = find_column(obs, args.nfeature_field)
    percent_mt_field = find_column(obs, args.percent_mt_field)
    percent_ribo_field = find_column(obs, args.percent_ribo_field)
    doublet_score_field = find_column(obs, args.doublet_score_field)
    predicted_doublet_field = find_column(obs, args.predicted_doublet_field)
    doublet_threshold_field = find_column(obs, args.doublet_threshold_field)
    cell_cycle_s_field = find_column(obs, args.cell_cycle_s_field)
    cell_cycle_g2m_field = find_column(obs, args.cell_cycle_g2m_field)
    cell_cycle_phase_field = find_column(obs, args.cell_cycle_phase_field)

    if args.umap_key not in adata.obsm:
        raise KeyError(f"UMAP key not found in obsm: {args.umap_key}")

    dataset_value = get_single_unique_value(obs, dataset_field) if dataset_field else None
    dataset_key = args.dataset_key or dataset_value or args.input_h5ad.stem
    dataset_name = args.dataset_name or dataset_key
    dataset_label = args.dataset_label

    inferred_project_accession, _ = parse_project_and_cell_line(dataset_key)
    project_accession = args.project_accession or inferred_project_accession

    cell_line = args.cell_line
    cell_line_tsv_path = args.cell_line_tsv
    if not cell_line_tsv_path.is_absolute():
        cell_line_tsv_path = (Path(__file__).resolve().parent.parent / cell_line_tsv_path).resolve()
    cell_line_meta = load_cell_line_metadata(cell_line_tsv_path)
    if cell_line not in cell_line_meta:
        available = ", ".join(sorted(cell_line_meta.keys()))
        raise ValueError(
            f"Cell line '{cell_line}' not found in {cell_line_tsv_path}. Available: {available}"
        )
    meta = cell_line_meta[cell_line]
    cell_line_slug = meta["cell_line_slug"]
    species_name = meta["species_name"]
    reference_genome = meta["reference_genome"]

    sql_import_root = args.sql_import_root or default_sql_import_root(cell_line, project_accession, dataset_key)

    sample_series = normalize_required_string_series(obs[sample_field], sample_field)
    cluster_series = normalize_required_string_series(obs[cluster_field], cluster_field)
    cluster_order_map = build_cluster_order_map(cluster_series.tolist())

    obs_export = pd.DataFrame(
        {
            "singlecell_dataset_id": args.singlecell_dataset_id,
            "cell_index": list(range(adata.n_obs)),
            "barcode": obs.index.astype(str),
            "sample_id": sample_series,
            "orig_ident": normalize_optional_string_series(obs[orig_ident_field]) if orig_ident_field else None,
            "condition_label": normalize_optional_string_series(obs[condition_field]) if condition_field else None,
            "cluster_id": cluster_series,
            "cluster_label": cluster_series.map(lambda value: f"{args.cluster_label_prefix}{value}"),
            "cluster_order": cluster_series.map(cluster_order_map),
            "cell_annotation": normalize_optional_string_series(obs[cell_annotation_field]) if cell_annotation_field else None,
            "annotation_confidence": normalize_optional_string_series(obs[confidence_field]) if confidence_field else None,
            "predicted_doublet": obs[predicted_doublet_field].map(to_nullable_bool_int) if predicted_doublet_field else None,
            "doublet_score": pd.to_numeric(obs[doublet_score_field], errors="coerce") if doublet_score_field else None,
            "doublet_threshold": pd.to_numeric(obs[doublet_threshold_field], errors="coerce") if doublet_threshold_field else None,
            "ncount_rna": pd.to_numeric(obs[ncount_field], errors="coerce") if ncount_field else None,
            "nfeature_rna": pd.to_numeric(obs[nfeature_field], errors="coerce") if nfeature_field else None,
            "percent_mt": pd.to_numeric(obs[percent_mt_field], errors="coerce") if percent_mt_field else None,
            "percent_ribo": pd.to_numeric(obs[percent_ribo_field], errors="coerce") if percent_ribo_field else None,
            "s_score": pd.to_numeric(obs[cell_cycle_s_field], errors="coerce") if cell_cycle_s_field else None,
            "g2m_score": pd.to_numeric(obs[cell_cycle_g2m_field], errors="coerce") if cell_cycle_g2m_field else None,
            "cell_cycle_phase": normalize_optional_string_series(obs[cell_cycle_phase_field]) if cell_cycle_phase_field else None,
        }
    )

    umap_values = adata.obsm[args.umap_key]
    if umap_values.shape[1] < 2:
        raise ValueError(f"UMAP embedding {args.umap_key} does not have two columns.")
    umap_rows = pd.DataFrame(
        {
            "singlecell_dataset_id": args.singlecell_dataset_id,
            "cell_index": list(range(adata.n_obs)),
            "umap_key": args.umap_key,
            "umap_1": umap_values[:, 0],
            "umap_2": umap_values[:, 1],
        }
    )

    score_prefixes = [prefix.strip() for prefix in args.signature_score_prefixes.split(",") if prefix.strip()]
    signature_columns: list[tuple[str, str, str]] = []
    for prefix in score_prefixes:
        score_type = prefix[:-2] if prefix.endswith("__") else prefix
        for column in obs.columns:
            if column.startswith(prefix):
                signature_columns.append((score_type, column[len(prefix):], column))
    signature_columns.sort(key=lambda item: (item[0], item[1]))

    if signature_columns:
        cell_index = np.arange(adata.n_obs, dtype=np.int64)
        signature_rows = pd.concat(
            [
                pd.DataFrame(
                    {
                        "singlecell_dataset_id": args.singlecell_dataset_id,
                        "cell_index": cell_index,
                        "score_type": score_type,
                        "signature_name": signature_name,
                        "score_value": pd.to_numeric(obs[column], errors="coerce").to_numpy(),
                    }
                )
                for score_type, signature_name, column in signature_columns
            ],
            ignore_index=True,
        )
        print(f"[signature] exporting {len(signature_columns)} score columns -> {len(signature_rows)} rows")
    else:
        signature_rows = None
        print(
            "[signature] no obs columns matched prefixes "
            f"'{args.signature_score_prefixes}'; cell signature score table will be skipped"
        )

    gene_symbol_series = resolve_gene_symbol_series(var_df, args.gene_symbol_source, args.gene_name_field)
    gene_name_series = (
        normalize_optional_string_series(var_df[args.gene_name_field])
        if args.gene_name_field in var_df.columns
        else pd.Series([None] * adata.n_vars, index=var_df.index, dtype=object)
    )
    highly_variable_series = (
        var_df["highly_variable"].map(to_nullable_bool_int)
        if "highly_variable" in var_df.columns
        else pd.Series([None] * adata.n_vars, index=var_df.index)
    )
    gene_rows = pd.DataFrame(
        {
            "singlecell_dataset_id": args.singlecell_dataset_id,
            "gene_index": list(range(adata.n_vars)),
            "gene_symbol": gene_symbol_series.astype(str).tolist(),
            "gene_name": gene_name_series.tolist(),
            "highly_variable": highly_variable_series.tolist(),
        }
    )

    stat_rows = (
        obs_export.groupby(
            ["singlecell_dataset_id", "sample_id", "condition_label", "cluster_id", "cluster_label", "cluster_order"],
            dropna=False,
        )
        .size()
        .reset_index(name="cell_count")
    )
    sample_totals = stat_rows.groupby(["singlecell_dataset_id", "sample_id"], dropna=False)["cell_count"].sum().rename("sample_cell_count")
    stat_rows = stat_rows.merge(sample_totals, on=["singlecell_dataset_id", "sample_id"], how="left")
    stat_rows["proportion"] = stat_rows["cell_count"] / stat_rows["sample_cell_count"]

    export_raw = False
    if args.raw_mode == "include":
        export_raw = True
    elif args.raw_mode == "auto" and adata.raw is not None:
        export_raw = True

    dataset_row = {
        "singlecell_dataset_id": args.singlecell_dataset_id,
        "datasource_code": args.datasource_code,
        "dataset_key": dataset_key,
        "dataset_name": dataset_name,
        "dataset_label": dataset_label,
        "project_accession": project_accession,
        "cell_line": cell_line,
        "cell_line_slug": cell_line_slug,
        "species_name": species_name,
        "reference_genome": reference_genome,
        "sample_count": int(obs_export["sample_id"].nunique()),
        "cell_count": int(adata.n_obs),
        "gene_count": int(adata.n_vars),
        "cluster_count": int(obs_export["cluster_id"].nunique()),
        "umap_key": args.umap_key,
        "default_cluster_field": args.cluster_field,
        "default_sample_field": args.sample_field,
        "expression_value_type": "normalized_log1p",
        "raw_value_type": "raw_count" if export_raw else None,
        "has_raw": 1 if export_raw else 0,
        "source_h5ad_path": str(args.input_h5ad.resolve()),
        "source_file_hash": compute_sha256(args.input_h5ad),
        "notes": (
            f"Imported from h5ad; cluster_field={args.cluster_field}; "
            f"sample_field={args.sample_field}; dataset_field={args.dataset_field or ''}"
        ),
    }

    base_name = sanitize_filename(f"{args.singlecell_dataset_id}_{dataset_key}")
    delimiter = "\t" if args.delimiter == "tsv" else ","
    dataset_file = args.output_dir / f"{DEFAULT_DATASET_TABLE}.{base_name}.{args.delimiter}"
    cell_file = args.output_dir / f"{DEFAULT_CELL_TABLE}.{base_name}.{args.delimiter}"
    umap_file = args.output_dir / f"{DEFAULT_UMAP_TABLE}.{base_name}.{args.delimiter}"
    gene_file = args.output_dir / f"{DEFAULT_GENE_TABLE}.{base_name}.{args.delimiter}"
    gene_exp_file = args.output_dir / f"{DEFAULT_GENE_EXP_TABLE}.{base_name}.{args.delimiter}"
    cluster_top_gene_file = args.output_dir / f"{DEFAULT_CLUSTER_TOP_GENE_TABLE}.{base_name}.{args.delimiter}"
    stat_file = args.output_dir / f"{DEFAULT_STAT_TABLE}.{base_name}.{args.delimiter}"
    signature_file = (
        args.output_dir / f"{DEFAULT_CELL_SIGNATURE_SCORE_TABLE}.{base_name}.{args.delimiter}"
        if signature_rows is not None
        else None
    )
    import_script_stem = args.import_script_stem or default_import_script_stem(args.table_prefix)
    manifest_stem = args.manifest_stem or default_manifest_stem(args.table_prefix)
    bash_file = args.output_dir / f"{sanitize_filename(import_script_stem)}.{base_name}.sh"
    manifest_file = args.output_dir / f"{sanitize_filename(manifest_stem)}.{base_name}.json"

    print("[write] dataset/cell/umap/gene/stat files")
    write_dict_rows(dataset_file, list(dataset_row.keys()), [dataset_row], delimiter)
    write_dict_rows(cell_file, list(obs_export.columns), obs_export.to_dict("records"), delimiter)
    write_dict_rows(umap_file, list(umap_rows.columns), umap_rows.to_dict("records"), delimiter)
    write_dict_rows(gene_file, list(gene_rows.columns), gene_rows.to_dict("records"), delimiter)
    write_dict_rows(stat_file, list(stat_rows.columns), stat_rows.to_dict("records"), delimiter)

    if signature_rows is not None:
        print("[write] cell signature score file")
        signature_rows.to_csv(
            signature_file,
            sep=delimiter,
            index=False,
            na_rep=r"\N",
            lineterminator="\n",
        )

    print("[write] gene expression file")
    matrix = materialize_sparse_matrix(adata.X)
    gene_exp_row_count = write_gene_exp_file(
        gene_exp_file,
        args.singlecell_dataset_id,
        matrix,
        "normalized_log1p",
        delimiter,
        args.gene_exp_progress_interval,
    )
    print("[write] cluster top gene file")
    cluster_top_gene_row_count = write_cluster_top_gene_file(
        cluster_top_gene_file,
        args.singlecell_dataset_id,
        matrix,
        gene_rows,
        obs_export,
        "normalized_log1p",
        args.cluster_top_gene_limit,
        delimiter,
    )
    del matrix
    gc.collect()

    if export_raw:
        if adata.raw is None:
            raise ValueError("raw-mode requested raw export, but adata.raw is missing.")
        raw_var_names = list(adata.raw.var_names.astype(str))
        current_var_names = list(adata.var_names.astype(str))
        if raw_var_names != current_var_names:
            raise ValueError(
                "adata.raw.var_names do not match adata.var_names. "
                "Disable raw export or harmonize the h5ad first."
            )
        raw_matrix = materialize_sparse_matrix(adata.raw.X)
        gene_exp_row_count += write_gene_exp_file(
            gene_exp_file,
            args.singlecell_dataset_id,
            raw_matrix,
            "raw_count",
            delimiter,
            args.gene_exp_progress_interval,
        )
        del raw_matrix
        gc.collect()

    file_map = {
        DEFAULT_DATASET_TABLE: dataset_file,
        DEFAULT_CELL_TABLE: cell_file,
        DEFAULT_UMAP_TABLE: umap_file,
        DEFAULT_GENE_TABLE: gene_file,
        DEFAULT_GENE_EXP_TABLE: gene_exp_file,
        DEFAULT_CLUSTER_TOP_GENE_TABLE: cluster_top_gene_file,
        DEFAULT_STAT_TABLE: stat_file,
    }
    if signature_file is not None:
        file_map[DEFAULT_CELL_SIGNATURE_SCORE_TABLE] = signature_file

    row_counts = {
        DEFAULT_DATASET_TABLE: 1,
        DEFAULT_CELL_TABLE: int(len(obs_export)),
        DEFAULT_UMAP_TABLE: int(len(umap_rows)),
        DEFAULT_GENE_TABLE: int(len(gene_rows)),
        DEFAULT_GENE_EXP_TABLE: int(gene_exp_row_count),
        DEFAULT_CLUSTER_TOP_GENE_TABLE: int(cluster_top_gene_row_count),
        DEFAULT_STAT_TABLE: int(len(stat_rows)),
    }
    if signature_rows is not None:
        row_counts[DEFAULT_CELL_SIGNATURE_SCORE_TABLE] = int(len(signature_rows))

    write_load_data_bash(
        bash_file,
        DEFAULT_DATABASE_NAME,
        args.delimiter,
        file_map,
        args.singlecell_dataset_id,
        args.mysql_user,
        args.mysql_password,
        args.mysql_socket,
        sql_import_root,
    )

    manifest = {
        "database_name": DEFAULT_DATABASE_NAME,
        "table_prefix": args.table_prefix,
        "datasource_code": args.datasource_code,
        "sql_import_root": sql_import_root,
        "singlecell_dataset_id": args.singlecell_dataset_id,
        "dataset_key": dataset_key,
        "dataset_name": dataset_name,
        "project_accession": project_accession,
        "cell_line": cell_line,
        "cell_line_slug": cell_line_slug,
        "input_h5ad": str(args.input_h5ad.resolve()),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "umap_key": args.umap_key,
        "cluster_field": args.cluster_field,
        "sample_field": args.sample_field,
        "condition_field": args.condition_field,
        "raw_exported": export_raw,
        "signature_score_prefixes": args.signature_score_prefixes,
        "signature_score_column_count": len(signature_columns),
        "files": {table_name: str(path.resolve()) for table_name, path in file_map.items()},
        "import_script": str(bash_file.resolve()),
        "row_counts": row_counts,
    }
    manifest_file.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")

    print("[done]")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
