#!/usr/bin/env python
"""Export an AnnData H5AD object to sparse Matrix Market input for the R scorer.

The R scoring workflow uses raw UMI counts for AUCell and UCell.  This helper
therefore exports ``adata.raw.X`` when it is available; otherwise it exports
``adata.X``.  It also carries forward all cell metadata and an existing UMAP.

Requirements: Python packages anndata, scipy, pandas and numpy.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import io as spio
from scipy import sparse


def make_unique(values: pd.Series) -> list[str]:
    """Make feature names unique while keeping the first occurrence unchanged."""
    result: list[str] = []
    seen: dict[str, int] = {}
    for value in values.fillna("").astype(str):
        base = value.strip() or "unknown_feature"
        if base not in seen:
            seen[base] = 0
            result.append(base)
        else:
            seen[base] += 1
            result.append(f"{base}.{seen[base]}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export raw H5AD counts, metadata and UMAP for R functional-state scoring."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input .h5ad file")
    parser.add_argument("--output-dir", required=True, type=Path, help="Output directory")
    args = parser.parse_args()

    if not args.input.is_file():
        raise FileNotFoundError(f"H5AD file does not exist: {args.input}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.input, backed="r")
    source = adata.raw if adata.raw is not None else adata
    source_label = "raw" if adata.raw is not None else "X"

    matrix = source.X.to_memory() if hasattr(source.X, "to_memory") else source.X
    if not sparse.issparse(matrix):
        matrix = sparse.csr_matrix(matrix)
    matrix = matrix.tocsc()  # genes x cells is required by Matrix Market / AUCell

    var = source.var.copy()
    gene_ids = pd.Series(source.var_names.astype(str), index=var.index)
    if "gene_name" in var.columns:
        gene_symbols = pd.Series(var["gene_name"].to_numpy(), index=var.index)
        gene_symbols = gene_symbols.where(gene_symbols.notna(), gene_ids)
        gene_symbols = gene_symbols.mask(gene_symbols.astype(str).str.strip() == "", gene_ids)
    else:
        gene_symbols = gene_ids
    gene_symbols = make_unique(gene_symbols)

    spio.mmwrite(args.output_dir / "matrix.mtx", matrix.transpose().tocoo())
    pd.Series(adata.obs_names.astype(str)).to_csv(
        args.output_dir / "barcodes.tsv", index=False, header=False
    )
    features = pd.DataFrame(
        {
            "gene_id": gene_ids.to_numpy(),
            "gene_symbol": gene_symbols,
            "feature_type": "Gene Expression",
        }
    )
    features.to_csv(args.output_dir / "features.tsv", sep="\t", index=False, header=False)

    metadata = adata.obs.copy()
    for column in metadata.columns:
        if pd.api.types.is_categorical_dtype(metadata[column]):
            metadata[column] = metadata[column].astype(str)
    metadata.index = adata.obs_names.astype(str)
    metadata.index.name = "cell"
    metadata.to_csv(args.output_dir / "metadata.csv")

    if "X_umap" in adata.obsm:
        umap = np.asarray(adata.obsm["X_umap"])
        if umap.ndim == 2 and umap.shape[0] == adata.n_obs and umap.shape[1] >= 2:
            pd.DataFrame(
                {"cell": adata.obs_names.astype(str), "UMAP_1": umap[:, 0], "UMAP_2": umap[:, 1]}
            ).to_csv(args.output_dir / "umap.csv", index=False)

    pd.DataFrame(
        {
            "item": ["input_h5ad", "expression_source", "n_cells", "n_genes"],
            "value": [str(args.input.resolve()), source_label, adata.n_obs, source.n_vars],
        }
    ).to_csv(args.output_dir / "export_manifest.csv", index=False)
    adata.file.close()
    print(f"Exported {source.n_vars} genes x {adata.n_obs} cells from {source_label} to {args.output_dir}")


if __name__ == "__main__":
    main()
