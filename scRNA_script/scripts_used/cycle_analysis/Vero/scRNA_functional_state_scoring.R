#!/usr/bin/env Rscript

# Functional-state scoring for one scRNA-seq dataset.
#
# Method:
#   * MSigDB Hallmark collection H is a fixed functional-state library.
#   * UCell is the primary per-cell score; AUCell is calculated from the same
#     gene sets as a continuous-score comparison.
#   * Seurat CellCycleScoring adds S.Score, G2M.Score and Phase to every cell.
#
# Inputs accepted by --input:
#   1. A Seurat .rds object;
#   2. A 10X Matrix Market directory (matrix.mtx[.gz], features.tsv[.gz],
#      barcodes.tsv[.gz]);
#   3. An .h5ad file.  For this case the companion Python exporter is called
#      automatically.  It uses adata.raw.X when available.
#
# Vero uses local Chlorocebus sabaeus Hallmark and cell-cycle ortholog tables.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(UCell)
  library(AUCell)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

option_list <- list(
  make_option(c("--input"), type = "character", default = NULL,
              help = "Input .rds, .h5ad, or 10X Matrix Market directory [required]"),
  make_option(c("--output-dir"), dest = "output_dir", type = "character", default = NULL,
              help = "Directory for result tables, figures and Seurat object [required]"),
  make_option(c("--species"), type = "character", default = "Chlorocebus sabaeus",
              help = "Species label recorded in output metadata [default %default]"),
  make_option(c("--db-species"), dest = "db_species", type = "character", default = "CS",
              help = "MSigDB database label recorded in output metadata [default %default]"),
  make_option(c("--hallmark-file"), dest = "hallmark_file", type = "character",
              default = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/scRNA_script/LTC/cycle_analysis/Vero/MSigDB_Hallmark/msigdb_hallmark_CS.csv",
              help = "Local MSigDB Hallmark CSV [default %default]"),
  make_option(c("--cell-cycle-file"), dest = "cell_cycle_file", type = "character",
              default = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/scRNA_script/LTC/cycle_analysis/Vero/MSigDB_Hallmark/cell_cycle_genes_CS.csv",
              help = "Local Chlorocebus sabaeus cell-cycle gene CSV [default %default]"),
  make_option(c("--assay"), type = "character", default = "RNA",
              help = "Seurat assay holding RNA counts [default %default]"),
  make_option(c("--cluster-col"), dest = "cluster_col", type = "character", default = "seurat_clusters",
              help = "Existing cell-level cluster column [default %default]"),
  make_option(c("--python"), type = "character", default = "python",
              help = "Python executable with anndata/scipy; used only for .h5ad input [default %default]"),
  make_option(c("--h5ad-export-dir"), dest = "h5ad_export_dir", type = "character", default = NULL,
              help = "Reusable Matrix Market export directory for .h5ad input [default: output-dir/h5ad_export]"),
  make_option(c("--reuse-h5ad-export"), dest = "reuse_h5ad_export", action = "store_true", default = FALSE,
              help = "Reuse a complete prior H5AD export instead of exporting again [default %default]"),
  make_option(c("--ucell-max-rank"), dest = "ucell_max_rank", type = "integer", default = NULL,
              help = "UCell maxRank. Default is the median detected-gene count, capped at gene count"),
  make_option(c("--aucell-max-rank"), dest = "aucell_max_rank", type = "integer", default = NULL,
              help = "AUCell aucMaxRank. Default is min(5%% of genes, 5th percentile of detected genes)"),
  make_option(c("--min-signature-genes"), dest = "min_signature_genes", type = "integer", default = 5,
              help = "Minimum matched genes required to score a Hallmark set [default %default]"),
  make_option(c("--min-cycle-genes"), dest = "min_cycle_genes", type = "integer", default = 10,
              help = "Minimum matched genes required in each cell-cycle set [default %default]"),
  make_option(c("--ncores"), type = "integer", default = 1,
              help = "CPU cores for UCell/AUCell; use 1 on Windows [default %default]"),
  make_option(c("--seed"), type = "integer", default = 1234,
              help = "Seed used if the script must create a UMAP [default %default]")
)

parser <- OptionParser(
  option_list = option_list,
  description = "Score Hallmark functional states with UCell and AUCell, then score cell cycle."
)
opt <- parse_args(parser)
set.seed(opt$seed)

if (is.null(opt$input) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("Both --input and --output-dir are required.", call. = FALSE)
}
if (!file.exists(opt$input) && !dir.exists(opt$input)) {
  stop("Input does not exist: ", opt$input, call. = FALSE)
}
if (!file.exists(opt$hallmark_file)) stop("Hallmark CSV does not exist: ", opt$hallmark_file, call. = FALSE)
if (!file.exists(opt$cell_cycle_file)) stop("Cell-cycle CSV does not exist: ", opt$cell_cycle_file, call. = FALSE)
if (opt$ncores < 1) stop("--ncores must be >= 1.", call. = FALSE)

input_path <- normalizePath(opt$input, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
object_dir <- file.path(output_dir, "seurat")
for (directory in c(output_dir, table_dir, figure_dir, object_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

script_path <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_path) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_path), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

`%||%` <- function(left, right) {
  if (is.null(left) || length(left) == 0 || is.na(left[[1]]) || identical(left[[1]], "")) right else left
}

read_tabular <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    read.delim(gzfile(path), header = FALSE, stringsAsFactors = FALSE, ...)
  } else {
    read.delim(path, header = FALSE, stringsAsFactors = FALSE, ...)
  }
}

read_lines <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    readLines(gzfile(path))
  } else {
    readLines(path)
  }
}

find_first_existing <- function(paths, label) {
  found <- paths[file.exists(paths)]
  if (length(found) == 0) stop("Missing ", label, ". Checked: ", paste(paths, collapse = "; "), call. = FALSE)
  found[[1]]
}

read_10x_directory <- function(data_dir, project_name) {
  matrix_file <- find_first_existing(file.path(data_dir, c("matrix.mtx", "matrix.mtx.gz")), "matrix.mtx")
  feature_file <- find_first_existing(file.path(data_dir, c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz")), "features.tsv/genes.tsv")
  barcode_file <- find_first_existing(file.path(data_dir, c("barcodes.tsv", "barcodes.tsv.gz")), "barcodes.tsv")

  counts <- if (grepl("\\.gz$", matrix_file, ignore.case = TRUE)) {
    Matrix::readMM(gzfile(matrix_file))
  } else {
    Matrix::readMM(matrix_file)
  }
  counts <- methods::as(counts, "dgCMatrix")
  features <- read_tabular(feature_file)
  barcodes <- read_lines(barcode_file)
  if (nrow(features) != nrow(counts) || length(barcodes) != ncol(counts)) {
    stop("Matrix dimensions do not agree with features or barcodes in ", data_dir, call. = FALSE)
  }
  genes <- if (ncol(features) >= 2) features[[2]] else features[[1]]
  genes[is.na(genes) | genes == ""] <- features[[1]][is.na(genes) | genes == ""]
  rownames(counts) <- make.unique(as.character(genes))
  colnames(counts) <- make.unique(as.character(barcodes))

  metadata_path <- file.path(data_dir, "metadata.csv")
  metadata <- data.frame(row.names = colnames(counts))
  if (file.exists(metadata_path)) {
    input_metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    shared_cells <- intersect(rownames(metadata), rownames(input_metadata))
    for (column in colnames(input_metadata)) {
      metadata[[column]] <- NA
      metadata[shared_cells, column] <- input_metadata[shared_cells, column]
    }
  }
  Seurat::CreateSeuratObject(counts = counts, project = project_name, meta.data = metadata)
}

restore_umap <- function(object, data_dir, assay) {
  umap_path <- file.path(data_dir, "umap.csv")
  if (!file.exists(umap_path)) return(object)
  umap <- read.csv(umap_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("cell", "UMAP_1", "UMAP_2")
  if (!all(required %in% colnames(umap))) {
    warning("umap.csv does not contain cell, UMAP_1 and UMAP_2; existing UMAP will not be restored.")
    return(object)
  }
  umap <- umap[match(colnames(object), umap$cell), required, drop = FALSE]
  if (anyNA(umap$UMAP_1) || anyNA(umap$UMAP_2)) {
    warning("Some cells in umap.csv do not match the count matrix; existing UMAP will not be restored.")
    return(object)
  }
  embedding <- as.matrix(umap[, c("UMAP_1", "UMAP_2")])
  rownames(embedding) <- colnames(object)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  object[["umap"]] <- Seurat::CreateDimReducObject(embeddings = embedding, key = "UMAP_", assay = assay)
  object
}

get_assay_matrix <- function(object, assay, layer) {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = layer)
  )
}

has_data_layer <- function(object, assay) {
  layers <- tryCatch(SeuratObject::Layers(object[[assay]]), error = function(e) character())
  if ("data" %in% layers) return(TRUE)
  data_matrix <- tryCatch(get_assay_matrix(object, assay, "data"), error = function(e) NULL)
  !is.null(data_matrix) && length(data_matrix@x) > 0
}

ensure_normalized <- function(object, assay) {
  if (!has_data_layer(object, assay)) {
    message("Normalizing RNA counts for CellCycleScoring.")
    object <- Seurat::NormalizeData(object, assay = assay, verbose = FALSE)
  }
  object
}

ensure_umap <- function(object, assay, cluster_col, seed) {
  if ("umap" %in% names(object@reductions)) return(object)
  message("No UMAP was supplied; creating a standard RNA PCA/UMAP embedding for visualization.")
  object <- Seurat::FindVariableFeatures(object, assay = assay, selection.method = "vst", nfeatures = min(2000, nrow(object)), verbose = FALSE)
  object <- Seurat::ScaleData(object, assay = assay, features = Seurat::VariableFeatures(object), verbose = FALSE)
  n_pcs <- min(30, length(Seurat::VariableFeatures(object)), ncol(object) - 1)
  if (n_pcs < 2) stop("At least three cells are required to calculate a UMAP.", call. = FALSE)
  object <- Seurat::RunPCA(object, assay = assay, features = Seurat::VariableFeatures(object), npcs = n_pcs, seed.use = seed, verbose = FALSE)
  dimensions <- seq_len(n_pcs)
  if (!(cluster_col %in% colnames(object[[]]))) {
    object <- Seurat::FindNeighbors(object, reduction = "pca", dims = dimensions, verbose = FALSE)
    object <- Seurat::FindClusters(object, resolution = 0.8, random.seed = seed, verbose = FALSE)
    object[[cluster_col]] <- object$seurat_clusters
  }
  Seurat::RunUMAP(object, reduction = "pca", dims = dimensions, seed.use = seed, verbose = FALSE)
}

match_symbols <- function(query_symbols, feature_symbols) {
  feature_keys <- toupper(as.character(feature_symbols))
  query_keys <- toupper(as.character(query_symbols))
  matched_index <- match(query_keys, feature_keys)
  unique(feature_symbols[matched_index[!is.na(matched_index)]])
}

write_csv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(data, connection, row.names = FALSE, quote = TRUE)
}

cluster_median_matrix <- function(score_matrix, clusters) {
  cluster_levels <- sort(unique(as.character(clusters)))
  result <- matrix(NA_real_, nrow = length(cluster_levels), ncol = ncol(score_matrix),
                   dimnames = list(cluster_levels, colnames(score_matrix)))
  for (index in seq_along(cluster_levels)) {
    cells <- as.character(clusters) == cluster_levels[[index]]
    result[index, ] <- apply(score_matrix[cells, , drop = FALSE], 2, median, na.rm = TRUE)
  }
  result
}

zscore_columns <- function(score_matrix) {
  z_scores <- scale(score_matrix)
  constant_columns <- apply(score_matrix, 2, function(x) {
    standard_deviation <- sd(x, na.rm = TRUE)
    standard_deviation == 0 || is.na(standard_deviation)
  })
  z_scores[, constant_columns] <- 0
  z_scores[is.na(z_scores)] <- 0
  z_scores
}

write_heatmap <- function(score_matrix, output_file, title) {
  grDevices::pdf(output_file, width = 15, height = max(5, nrow(score_matrix) * 0.45 + 2))
  pheatmap::pheatmap(
    score_matrix,
    cluster_rows = nrow(score_matrix) > 1,
    cluster_cols = ncol(score_matrix) > 1,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    main = title,
    fontsize_col = 6,
    fontsize_row = 9,
    border_color = NA
  )
  grDevices::dev.off()
}

save_plot_pages <- function(plots, output_file, width, height) {
  grDevices::pdf(output_file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (plot_item in plots) print(plot_item)
}

message("Input: ", input_path)
message("Output directory: ", output_dir)

if (grepl("\\.rds$", input_path, ignore.case = TRUE)) {
  object <- readRDS(input_path)
  if (!inherits(object, "Seurat")) stop("The .rds input is not a Seurat object.", call. = FALSE)
} else if (grepl("\\.h5ad$", input_path, ignore.case = TRUE)) {
  export_dir <- opt$h5ad_export_dir %||% file.path(output_dir, "h5ad_export")
  export_dir <- normalizePath(export_dir, winslash = "/", mustWork = FALSE)
  required_export <- file.path(export_dir, c("matrix.mtx", "features.tsv", "barcodes.tsv", "metadata.csv"))
  if (!opt$reuse_h5ad_export || !all(file.exists(required_export))) {
    exporter <- file.path(script_dir, "export_h5ad_for_functional_scoring.py")
    if (!file.exists(exporter)) stop("H5AD exporter was not found beside this script: ", exporter, call. = FALSE)
    dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
    message("Exporting raw H5AD counts to Matrix Market format.")
    exporter_args <- c(shQuote(exporter), "--input", shQuote(input_path), "--output-dir", shQuote(export_dir))
    export_output <- system2(opt$python, args = exporter_args, stdout = TRUE, stderr = TRUE)
    export_status <- attr(export_output, "status")
    if (!is.null(export_status) && export_status != 0) {
      stop("H5AD export failed:\n", paste(export_output, collapse = "\n"), call. = FALSE)
    }
    message(paste(export_output, collapse = "\n"))
  } else {
    message("Reusing existing H5AD export: ", export_dir)
  }
  object <- read_10x_directory(export_dir, basename(tools::file_path_sans_ext(input_path)))
  object <- restore_umap(object, export_dir, opt$assay)
} else if (dir.exists(input_path)) {
  object <- read_10x_directory(input_path, basename(input_path))
  object <- restore_umap(object, input_path, opt$assay)
} else {
  stop("Supported inputs are .rds, .h5ad, or a 10X Matrix Market directory.", call. = FALSE)
}

if (!(opt$assay %in% SeuratObject::Assays(object))) {
  stop("Assay '", opt$assay, "' is not present. Available assays: ", paste(SeuratObject::Assays(object), collapse = ", "), call. = FALSE)
}
SeuratObject::DefaultAssay(object) <- opt$assay
object <- ensure_normalized(object, opt$assay)
object <- ensure_umap(object, opt$assay, opt$cluster_col, opt$seed)
if (!(opt$cluster_col %in% colnames(object[[]]))) {
  stop("Cluster column '", opt$cluster_col, "' is absent after input processing.", call. = FALSE)
}

counts <- get_assay_matrix(object, opt$assay, "counts")
if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "dgCMatrix")
feature_symbols <- rownames(counts)
detected_genes <- Matrix::colSums(counts > 0)
median_detected <- as.integer(round(stats::median(detected_genes)))
q05_detected <- as.integer(floor(stats::quantile(detected_genes, probs = 0.05, names = FALSE)))

ucell_max_rank <- opt$ucell_max_rank %||% min(median_detected, nrow(counts))
aucell_default_rank <- min(as.integer(ceiling(0.05 * nrow(counts))), q05_detected)
aucell_max_rank <- opt$aucell_max_rank %||% aucell_default_rank
if (ucell_max_rank < 1 || aucell_max_rank < 1) stop("Calculated rank cut-off is < 1; inspect the count matrix.", call. = FALSE)

message("Loading local MSigDB Hallmark collection: ", opt$hallmark_file)
hallmark_table <- read.csv(opt$hallmark_file, stringsAsFactors = FALSE, check.names = FALSE)
required_hallmark_columns <- c("gene_symbol", "gs_name")
if (!all(required_hallmark_columns %in% colnames(hallmark_table))) {
  stop("Hallmark CSV must contain columns: ", paste(required_hallmark_columns, collapse = ", "), call. = FALSE)
}
if ("gs_collection" %in% colnames(hallmark_table)) hallmark_table <- hallmark_table[hallmark_table$gs_collection == "H", , drop = FALSE]
hallmark_table <- hallmark_table[!is.na(hallmark_table$gene_symbol) & nzchar(hallmark_table$gene_symbol) & !is.na(hallmark_table$gs_name) & nzchar(hallmark_table$gs_name), , drop = FALSE]
hallmark_sets <- split(hallmark_table$gene_symbol, hallmark_table$gs_name)
hallmark_sets <- lapply(hallmark_sets, function(x) unique(as.character(x)))
if (length(hallmark_sets) != 50) warning("Expected 50 Hallmark sets, but the local CSV contains ", length(hallmark_sets), ".")
matched_hallmark_sets <- lapply(hallmark_sets, match_symbols, feature_symbols = feature_symbols)
coverage <- data.frame(
  signature = names(hallmark_sets),
  genes_in_msigdb = lengths(hallmark_sets),
  genes_matched = lengths(matched_hallmark_sets),
  coverage_fraction = round(lengths(matched_hallmark_sets) / pmax(1, lengths(hallmark_sets)), 4),
  stringsAsFactors = FALSE
)
coverage$scored <- coverage$genes_matched >= opt$min_signature_genes
write.csv(coverage, file.path(table_dir, "hallmark_gene_set_coverage.csv"), row.names = FALSE, quote = TRUE)
active_hallmark_sets <- matched_hallmark_sets[coverage$scored]
if (length(active_hallmark_sets) == 0) {
  stop("No Hallmark gene set passed --min-signature-genes. Check species and gene identifiers.", call. = FALSE)
}
if (any(!coverage$scored)) {
  warning(sum(!coverage$scored), " Hallmark sets have fewer than ", opt$min_signature_genes,
          " matched genes and were not scored. See hallmark_gene_set_coverage.csv.")
}

message("Scoring ", length(active_hallmark_sets), " Hallmark sets with UCell (maxRank = ", ucell_max_rank, ").")
ucell_input <- active_hallmark_sets
names(ucell_input) <- paste0("UCell__", names(active_hallmark_sets))
object <- UCell::AddModuleScore_UCell(
  obj = object,
  features = ucell_input,
  maxRank = ucell_max_rank,
  ncores = opt$ncores,
  assay = opt$assay,
  slot = "counts",
  name = NULL
)
ucell_output_columns <- vapply(names(active_hallmark_sets), function(signature) {
  candidates <- c(paste0("UCell__", signature), paste0("UCell__", signature, "_UCell"))
  present <- candidates[candidates %in% colnames(object[[]])]
  if (length(present) != 1) stop("Could not identify UCell metadata column for ", signature, call. = FALSE)
  present[[1]]
}, character(1))
names(ucell_output_columns) <- names(active_hallmark_sets)

message("Scoring the same Hallmark sets with AUCell (aucMaxRank = ", aucell_max_rank, ").")
set.seed(opt$seed)
aucell_rankings <- AUCell::AUCell_buildRankings(
  exprMat = counts,
  plotStats = FALSE,
  nCores = opt$ncores,
  verbose = FALSE
)
aucell_scores <- AUCell::AUCell_calcAUC(
  geneSets = active_hallmark_sets,
  rankings = aucell_rankings,
  nCores = opt$ncores,
  aucMaxRank = aucell_max_rank,
  verbose = FALSE
)
aucell_matrix <- t(as.matrix(AUCell::getAUC(aucell_scores)))
aucell_matrix <- aucell_matrix[colnames(object), names(active_hallmark_sets), drop = FALSE]
aucell_output_columns <- paste0("AUCell__", colnames(aucell_matrix))
names(aucell_output_columns) <- colnames(aucell_matrix)
for (signature in names(aucell_output_columns)) {
  object[[aucell_output_columns[[signature]]]] <- aucell_matrix[, signature]
}
rm(aucell_rankings, aucell_scores)
gc(verbose = FALSE)

message("Calculating Seurat S and G2/M cell-cycle scores for every cell.")
cycle_table <- read.csv(opt$cell_cycle_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("phase", "gene_symbol") %in% colnames(cycle_table))) stop("Cell-cycle CSV must contain phase and gene_symbol columns.", call. = FALSE)
s_reference <- unique(cycle_table$gene_symbol[toupper(cycle_table$phase) == "S"])
g2m_reference <- unique(cycle_table$gene_symbol[toupper(cycle_table$phase) == "G2M"])
s_features <- match_symbols(s_reference, feature_symbols)
g2m_features <- match_symbols(g2m_reference, feature_symbols)
cycle_coverage <- data.frame(
  gene_set = c("S", "G2M"),
  genes_in_reference = c(length(s_reference), length(g2m_reference)),
  genes_matched = c(length(s_features), length(g2m_features)),
  stringsAsFactors = FALSE
)
write.csv(cycle_coverage, file.path(table_dir, "cell_cycle_gene_set_coverage.csv"), row.names = FALSE)
if (any(cycle_coverage$genes_matched < opt$min_cycle_genes)) {
  stop("Too few cell-cycle genes matched the matrix. Check the gene identifiers or provide a species-appropriate input matrix.", call. = FALSE)
}
set.seed(opt$seed)
object <- Seurat::CellCycleScoring(
  object = object,
  s.features = s_features,
  g2m.features = g2m_features,
  set.ident = FALSE
)
# Preserve completed per-cell scores before generating figures. This lets a
# plotting issue be inspected without losing the expensive AUCell calculation.
scored_rds <- file.path(object_dir, "functional_state_scored_seurat.rds")
saveRDS(object, scored_rds, compress = "gzip")

# Reporting is executed in a fresh R process.  This avoids a known incompatibility
# between Seurat 5.2 FeaturePlot objects and patchwork in this environment.
report_script <- file.path(script_dir, "scRNA_functional_state_report.R")
if (!file.exists(report_script)) stop("Reporting script was not found: ", report_script, call. = FALSE)
rscript_executable <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript_executable)) rscript_executable <- file.path(R.home("bin"), "Rscript")
report_args <- c(
  shQuote(report_script),
  paste0("--input-rds=", shQuote(normalizePath(scored_rds, winslash = "/", mustWork = TRUE))),
  paste0("--output-dir=", shQuote(output_dir)),
  paste0("--cluster-col=", shQuote(opt$cluster_col)),
  paste0("--species=", shQuote(opt$species)),
  paste0("--db-species=", shQuote(opt$db_species)),
  paste0("--ucell-max-rank=", ucell_max_rank),
  paste0("--aucell-max-rank=", aucell_max_rank),
  paste0("--ncores=", opt$ncores),
  paste0("--min-signature-genes=", opt$min_signature_genes),
  paste0("--min-cycle-genes=", opt$min_cycle_genes),
  paste0("--seed=", opt$seed)
)
message("Writing tables and figures from the scored Seurat checkpoint.")
report_output <- system2(rscript_executable, args = report_args, stdout = TRUE, stderr = TRUE)
report_status <- attr(report_output, "status")
if (!is.null(report_status) && report_status != 0) {
  stop("Functional-state reporting failed:\n", paste(report_output, collapse = "\n"), call. = FALSE)
}
message(paste(report_output, collapse = "\n"))
message("Completed successfully.")
message("Scored Seurat object: ", scored_rds)

if (FALSE) {
metadata <- object[[]]
clusters <- as.character(metadata[[opt$cluster_col]])
ucell_matrix <- as.matrix(metadata[, unname(ucell_output_columns), drop = FALSE])
colnames(ucell_matrix) <- names(ucell_output_columns)
aucell_matrix <- as.matrix(metadata[, unname(aucell_output_columns), drop = FALSE])
colnames(aucell_matrix) <- names(aucell_output_columns)

ucell_cluster_medians <- cluster_median_matrix(ucell_matrix, clusters)
aucell_cluster_medians <- cluster_median_matrix(aucell_matrix, clusters)
write.csv(ucell_cluster_medians, file.path(table_dir, "hallmark_ucell_cluster_medians.csv"), row.names = TRUE, quote = TRUE)
write.csv(aucell_cluster_medians, file.path(table_dir, "hallmark_aucell_cluster_medians.csv"), row.names = TRUE, quote = TRUE)
write_heatmap(zscore_columns(ucell_cluster_medians),
              file.path(figure_dir, "hallmark_ucell_cluster_medians_zscore.pdf"),
              "Hallmark cluster medians: UCell (z-score by signature)")
write_heatmap(zscore_columns(aucell_cluster_medians),
              file.path(figure_dir, "hallmark_aucell_cluster_medians_zscore.pdf"),
              "Hallmark cluster medians: AUCell (z-score by signature)")

agreement <- data.frame(
  signature = names(active_hallmark_sets),
  ucell_column = unname(ucell_output_columns[names(active_hallmark_sets)]),
  aucell_column = unname(aucell_output_columns[names(active_hallmark_sets)]),
  spearman_rho = vapply(names(active_hallmark_sets), function(signature) {
    suppressWarnings(stats::cor(ucell_matrix[, signature], aucell_matrix[, signature], method = "spearman"))
  }, numeric(1)),
  stringsAsFactors = FALSE
)
write.csv(agreement, file.path(table_dir, "ucell_aucell_spearman_agreement.csv"), row.names = FALSE, quote = TRUE)

phase_levels <- c("G1", "S", "G2M")
phase_table <- as.data.frame(table(
  cluster = factor(clusters, levels = sort(unique(clusters))),
  phase = factor(as.character(metadata$Phase), levels = phase_levels)
), stringsAsFactors = FALSE)
colnames(phase_table)[3] <- "cell_count"
phase_table$proportion <- phase_table$cell_count / ave(phase_table$cell_count, phase_table$cluster, FUN = sum)
write.csv(phase_table, file.path(table_dir, "cell_cycle_phase_by_cluster.csv"), row.names = FALSE, quote = TRUE)

SeuratObject::Idents(object) <- opt$cluster_col
cycle_feature_plot <- Seurat::FeaturePlot(
  object, features = c("S.Score", "G2M.Score"), reduction = "umap", order = TRUE, ncol = 2
)
cycle_phase_plot <- Seurat::DimPlot(object, reduction = "umap", group.by = "Phase", label = TRUE, repel = TRUE)
cycle_violin_plot <- Seurat::VlnPlot(
  object, features = c("S.Score", "G2M.Score"), group.by = opt$cluster_col, pt.size = 0, ncol = 2
)
phase_bar_plot <- ggplot(phase_table, aes(x = cluster, y = proportion, fill = phase)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cluster", y = "Cells", fill = "Cell-cycle phase", title = "Cell-cycle phase composition by cluster") +
  theme_classic(base_size = 11)
save_plot_pages(
  list(cycle_feature_plot, cycle_phase_plot),
  file.path(figure_dir, "cell_cycle_umap.pdf"),
  width = 12, height = 8
)
save_plot_pages(
  list(cycle_violin_plot, phase_bar_plot),
  file.path(figure_dir, "cell_cycle_by_cluster.pdf"),
  width = 14, height = 8
)

selected_hallmarks <- intersect(
  c(
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_APOPTOSIS"
  ),
  names(active_hallmark_sets)
)
if (length(selected_hallmarks) > 0) {
  hallmark_umap_plots <- list()
  for (signature in selected_hallmarks) {
    hallmark_umap_plots[[length(hallmark_umap_plots) + 1]] <-
      Seurat::FeaturePlot(object, features = ucell_output_columns[[signature]], reduction = "umap", order = TRUE,
                          min.cutoff = "q05", max.cutoff = "q95")
    hallmark_umap_plots[[length(hallmark_umap_plots) + 1]] <-
      Seurat::FeaturePlot(object, features = aucell_output_columns[[signature]], reduction = "umap", order = TRUE,
                          min.cutoff = "q05", max.cutoff = "q95")
  }
  save_plot_pages(
    hallmark_umap_plots,
    file.path(figure_dir, "selected_hallmark_ucell_aucell_umap.pdf"),
    width = 7, height = 6
  )
}

scatter_hallmarks <- head(agreement$signature[order(agreement$spearman_rho, decreasing = TRUE, na.last = TRUE)], 8)
if (length(scatter_hallmarks) > 0) {
  agreement_plots <- lapply(scatter_hallmarks, function(signature) {
    plot_data <- data.frame(UCell = ucell_matrix[, signature], AUCell = aucell_matrix[, signature])
    ggplot(plot_data, aes(x = UCell, y = AUCell)) +
      geom_point(size = 0.25, alpha = 0.12, color = "#2166AC") +
      labs(title = signature, subtitle = paste0("Spearman rho = ", round(agreement$spearman_rho[agreement$signature == signature], 3))) +
      theme_classic(base_size = 10)
  })
  save_plot_pages(
    agreement_plots,
    file.path(figure_dir, "ucell_aucell_agreement_scatter.pdf"),
    width = 6, height = 5
  )
}

cell_level_columns <- unique(c(
  opt$cluster_col, "S.Score", "G2M.Score", "Phase",
  unname(ucell_output_columns), unname(aucell_output_columns)
))
cell_level_scores <- data.frame(cell = rownames(metadata), metadata[, cell_level_columns, drop = FALSE], check.names = FALSE)
write_csv_gz(cell_level_scores, file.path(table_dir, "cell_level_functional_and_cell_cycle_scores.csv.gz"))

parameters <- data.frame(
  parameter = c("input", "hallmark_file", "cell_cycle_file", "species", "db_species", "assay", "cluster_col", "n_cells", "n_genes",
                "median_detected_genes", "q05_detected_genes", "ucell_max_rank", "aucell_max_rank",
                "hallmark_sets_returned", "hallmark_sets_scored", "min_signature_genes", "min_cycle_genes", "ncores"),
  value = c(input_path, normalizePath(opt$hallmark_file, winslash = "/", mustWork = TRUE), normalizePath(opt$cell_cycle_file, winslash = "/", mustWork = TRUE), opt$species, opt$db_species, opt$assay, opt$cluster_col, ncol(counts), nrow(counts),
            median_detected, q05_detected, ucell_max_rank, aucell_max_rank, length(hallmark_sets),
            length(active_hallmark_sets), opt$min_signature_genes, opt$min_cycle_genes, opt$ncores),
  stringsAsFactors = FALSE
)
write.csv(parameters, file.path(table_dir, "scoring_parameters.csv"), row.names = FALSE, quote = TRUE)
package_versions <- data.frame(
  package = c("Seurat", "UCell", "AUCell", "pheatmap", "optparse"),
  version = vapply(c("Seurat", "UCell", "AUCell", "pheatmap", "optparse"),
                   function(package) as.character(utils::packageVersion(package)), character(1)),
  stringsAsFactors = FALSE
)
write.csv(package_versions, file.path(table_dir, "package_versions.csv"), row.names = FALSE)

saveRDS(object, file.path(object_dir, "functional_state_scored_seurat.rds"), compress = "gzip")
message("Completed successfully.")
message("Scored Seurat object: ", file.path(object_dir, "functional_state_scored_seurat.rds"))
}
