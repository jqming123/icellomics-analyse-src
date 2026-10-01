#!/usr/bin/env Rscript

# Recalculate UCell at the same rank window used by AUCell, while preserving
# the existing UCell (default window) and AUCell scores in a scored Seurat RDS.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(UCell)
  library(Matrix)
  library(pheatmap)
})

option_list <- list(
  make_option(c("--input-rds"), dest = "input_rds", type = "character", default = NULL,
              help = "Existing functional_state_scored_seurat.rds [required]"),
  make_option(c("--output-dir"), dest = "output_dir", type = "character", default = NULL,
              help = "Directory for the rank-window comparison [required]"),
  make_option(c("--species"), type = "character", default = "Mus musculus"),
  make_option(c("--db-species"), dest = "db_species", type = "character", default = "MM"),
  make_option(c("--hallmark-file"), dest = "hallmark_file", type = "character",
              default = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/scRNA_script/LTC/cycle_analysis/SP20/MSigDB_Hallmark/msigdb_hallmark_MM.csv"),
  make_option(c("--cluster-col"), dest = "cluster_col", type = "character", default = "seurat_clusters"),
  make_option(c("--old-ucell-max-rank"), dest = "old_ucell_max_rank", type = "integer", default = NULL,
              help = "Original UCell maxRank stored in the input object [default: read from --manifest]"),
  make_option(c("--ucell-max-rank"), dest = "ucell_max_rank", type = "integer", default = NULL,
              help = "New UCell maxRank matched to AUCell aucMaxRank [default: read from --manifest]"),
  make_option(c("--manifest"), dest = "manifest", type = "character", default = NULL,
              help = "scoring_parameters.csv from functional scoring; supplies --old-ucell-max-rank (ucell_max_rank) and --ucell-max-rank (aucell_max_rank) when the explicit flags are omitted [default %default]"),
  make_option(c("--min-signature-genes"), dest = "min_signature_genes", type = "integer", default = 5),
  make_option(c("--ncores"), type = "integer", default = 1),
  make_option(c("--save-rds"), dest = "save_rds", action = "store_true", default = FALSE,
              help = "Save a duplicate Seurat RDS containing matched-rank UCell columns [default %default]")
)
parser <- OptionParser(option_list = option_list, description = "Compare UCell at a matched rank window against existing AUCell scores.")
opt <- parse_args(parser)
if (is.null(opt$input_rds) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("--input-rds and --output-dir are required.", call. = FALSE)
}
if (!file.exists(opt$hallmark_file)) stop("Hallmark CSV does not exist: ", opt$hallmark_file, call. = FALSE)
if (!file.exists(opt$input_rds)) stop("Input RDS does not exist: ", opt$input_rds, call. = FALSE)
if (is.null(opt$old_ucell_max_rank) || is.null(opt$ucell_max_rank)) {
  if (is.null(opt$manifest)) {
    stop("--old-ucell-max-rank and --ucell-max-rank are required (or supply --manifest).", call. = FALSE)
  }
  if (!file.exists(opt$manifest)) stop("Manifest does not exist: ", opt$manifest, call. = FALSE)
  manifest_table <- read.csv(opt$manifest, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("parameter", "value") %in% colnames(manifest_table))) {
    stop("Manifest must contain 'parameter' and 'value' columns.", call. = FALSE)
  }
  manifest_value <- function(name) {
    row <- manifest_table[manifest_table$parameter == name, "value"]
    if (length(row) != 1 || is.na(row)) stop("Manifest missing parameter: ", name, call. = FALSE)
    as.integer(row)
  }
  if (is.null(opt$old_ucell_max_rank)) opt$old_ucell_max_rank <- manifest_value("ucell_max_rank")
  if (is.null(opt$ucell_max_rank)) opt$ucell_max_rank <- manifest_value("aucell_max_rank")
}
if (opt$old_ucell_max_rank < 1 || opt$ucell_max_rank < 1 || opt$ncores < 1) stop("Rank and core count must be positive.", call. = FALSE)
old_ucell_rank <- opt$old_ucell_max_rank
matched_rank <- opt$ucell_max_rank

output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
for (directory in c(output_dir, table_dir, figure_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

get_counts <- function(object, assay) {
  tryCatch(
    SeuratObject::GetAssayData(object, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = "counts")
  )
}

match_symbols <- function(query_symbols, feature_symbols) {
  feature_keys <- toupper(as.character(feature_symbols))
  matched <- match(toupper(as.character(query_symbols)), feature_keys)
  unique(feature_symbols[matched[!is.na(matched)]])
}

cluster_median_matrix <- function(score_matrix, clusters) {
  cluster_levels <- sort(unique(as.character(clusters)))
  result <- matrix(NA_real_, nrow = length(cluster_levels), ncol = ncol(score_matrix),
                   dimnames = list(cluster_levels, colnames(score_matrix)))
  for (index in seq_along(cluster_levels)) {
    in_cluster <- as.character(clusters) == cluster_levels[[index]]
    result[index, ] <- apply(score_matrix[in_cluster, , drop = FALSE], 2, median, na.rm = TRUE)
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

write_heatmap <- function(score_matrix, output_file, title, limits = NULL) {
  grDevices::pdf(output_file, width = 15, height = max(5, nrow(score_matrix) * 0.45 + 2))
  arguments <- list(
    mat = score_matrix,
    cluster_rows = nrow(score_matrix) > 1,
    cluster_cols = ncol(score_matrix) > 1,
    main = title,
    fontsize_col = 6,
    fontsize_row = 9,
    border_color = NA
  )
  if (is.null(limits)) {
    arguments$color <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
  } else {
    arguments$color <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
    arguments$breaks <- seq(limits[[1]], limits[[2]], length.out = 102)
  }
  do.call(pheatmap::pheatmap, arguments)
  grDevices::dev.off()
}

write_csv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(data, connection, row.names = FALSE, quote = TRUE)
}

agreement_metrics <- function(left_scores, right_scores, left_cluster, right_cluster, signatures) {
  data.frame(
    signature = signatures,
    cell_spearman = vapply(signatures, function(signature) {
      suppressWarnings(stats::cor(left_scores[, signature], right_scores[, signature], method = "spearman"))
    }, numeric(1)),
    cluster_spearman = vapply(signatures, function(signature) {
      suppressWarnings(stats::cor(left_cluster[, signature], right_cluster[, signature], method = "spearman"))
    }, numeric(1)),
    cluster_z_pearson = vapply(signatures, function(signature) {
      suppressWarnings(stats::cor(zscore_columns(left_cluster)[, signature], zscore_columns(right_cluster)[, signature]))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

message("Loading scored Seurat object.")
object <- readRDS(opt$input_rds)
if (!inherits(object, "Seurat")) stop("Input RDS is not a Seurat object.", call. = FALSE)
assay <- SeuratObject::DefaultAssay(object)
metadata <- object[[]]
if (!(opt$cluster_col %in% colnames(metadata))) stop("Cluster column is absent: ", opt$cluster_col, call. = FALSE)

existing_ucell_columns <- grep("^UCell__HALLMARK_", colnames(metadata), value = TRUE)
existing_aucell_columns <- grep("^AUCell__HALLMARK_", colnames(metadata), value = TRUE)
signatures <- intersect(sub("^UCell__", "", existing_ucell_columns), sub("^AUCell__", "", existing_aucell_columns))
if (length(signatures) == 0) stop("No shared existing UCell/AUCell Hallmark scores were found.", call. = FALSE)

counts <- get_counts(object, assay)
if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "dgCMatrix")
message("Retrieving the same Hallmark collection and matching genes to the count matrix.")
hallmark_table <- read.csv(opt$hallmark_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("gene_symbol", "gs_name") %in% colnames(hallmark_table))) stop("Hallmark CSV must contain gene_symbol and gs_name columns.", call. = FALSE)
if ("gs_collection" %in% colnames(hallmark_table)) hallmark_table <- hallmark_table[hallmark_table$gs_collection == "H", , drop = FALSE]
hallmark_sets <- split(hallmark_table$gene_symbol, hallmark_table$gs_name)
matched_sets <- lapply(hallmark_sets, match_symbols, feature_symbols = rownames(counts))
matched_sets <- matched_sets[names(matched_sets) %in% signatures]
matched_sets <- matched_sets[lengths(matched_sets) >= opt$min_signature_genes]
signatures <- intersect(signatures, names(matched_sets))
matched_sets <- matched_sets[signatures]
if (length(signatures) == 0) stop("No Hallmark signature passed the matched-gene threshold.", call. = FALSE)

message("Adding UCell scores with maxRank = ", opt$ucell_max_rank, ".")
new_features <- matched_sets
names(new_features) <- paste0("UCell", opt$ucell_max_rank, "__", names(matched_sets))
object <- UCell::AddModuleScore_UCell(
  obj = object,
  features = new_features,
  maxRank = opt$ucell_max_rank,
  ncores = opt$ncores,
  assay = assay,
  slot = "counts",
  name = NULL
)
metadata <- object[[]]
new_columns <- vapply(signatures, function(signature) {
  base_name <- paste0("UCell", opt$ucell_max_rank, "__", signature)
  candidates <- c(base_name, paste0(base_name, "_UCell"))
  present <- candidates[candidates %in% colnames(metadata)]
  if (length(present) != 1) stop("Could not identify new UCell column for ", signature, call. = FALSE)
  present[[1]]
}, character(1))
names(new_columns) <- signatures

old_ucell_columns <- setNames(paste0("UCell__", signatures), signatures)
aucell_columns <- setNames(paste0("AUCell__", signatures), signatures)
old_ucell <- as.matrix(metadata[, unname(old_ucell_columns), drop = FALSE])
new_ucell <- as.matrix(metadata[, unname(new_columns), drop = FALSE])
aucell <- as.matrix(metadata[, unname(aucell_columns), drop = FALSE])
colnames(old_ucell) <- signatures
colnames(new_ucell) <- signatures
colnames(aucell) <- signatures
clusters <- as.character(metadata[[opt$cluster_col]])

old_cluster <- cluster_median_matrix(old_ucell, clusters)
new_cluster <- cluster_median_matrix(new_ucell, clusters)
aucell_cluster <- cluster_median_matrix(aucell, clusters)
write.csv(new_cluster, file.path(table_dir, paste0("hallmark_ucell", opt$ucell_max_rank, "_cluster_medians.csv")), row.names = TRUE, quote = TRUE)
write_heatmap(zscore_columns(new_cluster), file.path(figure_dir, paste0("hallmark_ucell", opt$ucell_max_rank, "_cluster_medians_zscore.pdf")),
              paste0("Hallmark cluster medians: UCell maxRank = ", opt$ucell_max_rank))

agreement_old_aucell <- agreement_metrics(old_ucell, aucell, old_cluster, aucell_cluster, signatures)
agreement_new_aucell <- agreement_metrics(new_ucell, aucell, new_cluster, aucell_cluster, signatures)
agreement_old_new <- agreement_metrics(old_ucell, new_ucell, old_cluster, new_cluster, signatures)
comparison <- setNames(data.frame(
  signature = signatures,
  agreement_old_aucell$cell_spearman,
  agreement_new_aucell$cell_spearman,
  agreement_new_aucell$cell_spearman - agreement_old_aucell$cell_spearman,
  agreement_old_aucell$cluster_spearman,
  agreement_new_aucell$cluster_spearman,
  agreement_new_aucell$cluster_spearman - agreement_old_aucell$cluster_spearman,
  agreement_old_aucell$cluster_z_pearson,
  agreement_new_aucell$cluster_z_pearson,
  agreement_new_aucell$cluster_z_pearson - agreement_old_aucell$cluster_z_pearson,
  stringsAsFactors = FALSE
), c(
  "signature",
  paste0("cell_spearman_ucell", old_ucell_rank, "_aucell"),
  paste0("cell_spearman_ucell", matched_rank, "_aucell"),
  "cell_spearman_change",
  paste0("cluster_spearman_ucell", old_ucell_rank, "_aucell"),
  paste0("cluster_spearman_ucell", matched_rank, "_aucell"),
  "cluster_spearman_change",
  paste0("cluster_z_pearson_ucell", old_ucell_rank, "_aucell"),
  paste0("cluster_z_pearson_ucell", matched_rank, "_aucell"),
  "cluster_z_pearson_change"
))
write.csv(comparison, file.path(table_dir, "ucell_rank_window_aucell_agreement.csv"), row.names = FALSE, quote = TRUE)
write.csv(agreement_old_new, file.path(table_dir, paste0("ucell", old_ucell_rank, "_ucell", matched_rank, "_agreement.csv")), row.names = FALSE, quote = TRUE)

old_difference <- zscore_columns(old_cluster) - zscore_columns(aucell_cluster)
new_difference <- zscore_columns(new_cluster) - zscore_columns(aucell_cluster)
shared_limit <- max(abs(c(old_difference, new_difference)), na.rm = TRUE)
write_heatmap(old_difference, file.path(figure_dir, paste0("ucell", old_ucell_rank, "_minus_aucell_cluster_zscore.pdf")),
              paste0("Cluster z-score difference: UCell ", old_ucell_rank, " minus AUCell ", matched_rank), c(-shared_limit, shared_limit))
write_heatmap(new_difference, file.path(figure_dir, paste0("ucell", matched_rank, "_minus_aucell_cluster_zscore.pdf")),
              paste0("Cluster z-score difference: UCell ", matched_rank, " minus AUCell ", matched_rank), c(-shared_limit, shared_limit))

summary_table <- data.frame(
  comparison = c(paste0("UCell ", old_ucell_rank, " vs AUCell ", matched_rank),
                 paste0("UCell ", matched_rank, " vs AUCell ", matched_rank),
                 paste0("UCell ", old_ucell_rank, " vs UCell ", matched_rank)),
  median_cell_spearman = c(
    median(agreement_old_aucell$cell_spearman, na.rm = TRUE),
    median(agreement_new_aucell$cell_spearman, na.rm = TRUE),
    median(agreement_old_new$cell_spearman, na.rm = TRUE)
  ),
  median_cluster_spearman = c(
    median(agreement_old_aucell$cluster_spearman, na.rm = TRUE),
    median(agreement_new_aucell$cluster_spearman, na.rm = TRUE),
    median(agreement_old_new$cluster_spearman, na.rm = TRUE)
  ),
  median_cluster_z_pearson = c(
    median(agreement_old_aucell$cluster_z_pearson, na.rm = TRUE),
    median(agreement_new_aucell$cluster_z_pearson, na.rm = TRUE),
    median(agreement_old_new$cluster_z_pearson, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_table, file.path(table_dir, "ucell_rank_window_comparison_summary.csv"), row.names = FALSE, quote = TRUE)

cell_level <- data.frame(cell = rownames(metadata), new_ucell, check.names = FALSE)
write_csv_gz(cell_level, file.path(table_dir, paste0("cell_level_ucell", opt$ucell_max_rank, "_scores.csv.gz")))
if (opt$save_rds) saveRDS(object, file.path(output_dir, paste0("functional_state_scored_ucell", opt$ucell_max_rank, ".rds")), compress = "gzip")

message("Completed UCell rank-window comparison for ", length(signatures), " Hallmark signatures.")
