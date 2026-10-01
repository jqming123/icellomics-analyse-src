#!/usr/bin/env Rscript

# Produce tables and figures from a Seurat object already scored by
# scRNA_functional_state_scoring.R.  This script deliberately plots extracted
# UMAP coordinates with ggplot2 rather than composing Seurat plot objects, so
# it is compatible with Seurat 5.2 and the installed patchwork version.

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(pheatmap)
})

options_list <- list(
  make_option(c("--input-rds"), dest = "input_rds", type = "character", default = NULL,
              help = "Scored Seurat RDS [required]"),
  make_option(c("--output-dir"), dest = "output_dir", type = "character", default = NULL,
              help = "Result directory [required]"),
  make_option(c("--cluster-col"), dest = "cluster_col", type = "character", default = "seurat_clusters",
              help = "Cluster field in Seurat metadata [default %default]"),
  make_option(c("--species"), type = "character", default = "Canis lupus familiaris"),
  make_option(c("--db-species"), dest = "db_species", type = "character", default = "CLF"),
  make_option(c("--ucell-max-rank"), dest = "ucell_max_rank", type = "integer", default = NA_integer_),
  make_option(c("--aucell-max-rank"), dest = "aucell_max_rank", type = "integer", default = NA_integer_),
  make_option(c("--ncores"), type = "integer", default = 1),
  make_option(c("--min-signature-genes"), dest = "min_signature_genes", type = "integer", default = NA_integer_),
  make_option(c("--min-cycle-genes"), dest = "min_cycle_genes", type = "integer", default = NA_integer_),
  make_option(c("--seed"), type = "integer", default = 1234,
              help = "Global random seed [default %default]")
)
parser <- OptionParser(option_list = options_list, description = "Write functional-state score tables and figures from a scored Seurat object.")
opt <- parse_args(parser)
if (is.null(opt$input_rds) || is.null(opt$output_dir)) {
  print_help(parser)
  stop("--input-rds and --output-dir are required.", call. = FALSE)
}
if (!file.exists(opt$input_rds)) stop("RDS file does not exist: ", opt$input_rds, call. = FALSE)

output_dir <- normalizePath(opt$output_dir, winslash = "/", mustWork = FALSE)
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
for (directory in c(output_dir, table_dir, figure_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)

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

write_heatmap <- function(score_matrix, output_file, title) {
  grDevices::pdf(output_file, width = 15, height = max(5, nrow(score_matrix) * 0.45 + 2))
  pheatmap::pheatmap(
    score_matrix,
    cluster_rows = nrow(score_matrix) > 1,
    cluster_cols = ncol(score_matrix) > 1,
    color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
    main = title, fontsize_col = 6, fontsize_row = 9, border_color = NA
  )
  grDevices::dev.off()
}

save_plot_pages <- function(plots, output_file, width, height) {
  grDevices::pdf(output_file, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (plot_item in plots) print(plot_item)
}

write_csv_gz <- function(data, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(data, connection, row.names = FALSE, quote = TRUE)
}

make_score_umap <- function(embedding, values, title) {
  plot_data <- data.frame(UMAP_1 = embedding[, 1], UMAP_2 = embedding[, 2], score = values)
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = score)) +
    geom_point(size = 0.22, alpha = 0.7) +
    scale_color_gradientn(colors = c("#440154", "#21908C", "#FDE725"), na.value = "grey85") +
    coord_equal() +
    labs(title = title, color = "Score") +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))
}

make_group_umap <- function(embedding, groups, title, legend_title) {
  plot_data <- data.frame(UMAP_1 = embedding[, 1], UMAP_2 = embedding[, 2], group = as.factor(groups))
  ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = group)) +
    geom_point(size = 0.22, alpha = 0.75) +
    coord_equal() +
    labs(title = title, color = legend_title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5))
}

message("Loading scored Seurat object: ", opt$input_rds)
object <- readRDS(opt$input_rds)
if (!inherits(object, "Seurat")) stop("--input-rds is not a Seurat object.", call. = FALSE)
metadata <- object[[]]
if (!(opt$cluster_col %in% colnames(metadata))) stop("Cluster column is absent: ", opt$cluster_col, call. = FALSE)
if (!all(c("S.Score", "G2M.Score", "Phase") %in% colnames(metadata))) stop("Cell-cycle columns are absent from the RDS.", call. = FALSE)
if (!("umap" %in% names(object@reductions))) stop("The scored object does not contain a UMAP reduction.", call. = FALSE)

# Seurat::CellCycleScoring assigns G1 when both S and G2/M scores are below
# zero; it does not generate a separate G1 gene-set score.  This derived index
# is positive exactly under that G1 rule and supports a continuous cluster
# violin plot. Higher values indicate lower S and G2/M programme activity. It
# should not be interpreted as an independent G1 gene-set enrichment score.
metadata$G1.RelativeScore <- -pmax(metadata$S.Score, metadata$G2M.Score)

ucell_columns <- grep("^UCell__HALLMARK_", colnames(metadata), value = TRUE)
aucell_columns <- grep("^AUCell__HALLMARK_", colnames(metadata), value = TRUE)
ucell_signatures <- sub("^UCell__", "", ucell_columns)
aucell_signatures <- sub("^AUCell__", "", aucell_columns)
common_signatures <- intersect(ucell_signatures, aucell_signatures)
if (length(common_signatures) == 0) stop("No shared UCell/AUCell Hallmark metadata columns were found.", call. = FALSE)
ucell_columns <- setNames(paste0("UCell__", common_signatures), common_signatures)
aucell_columns <- setNames(paste0("AUCell__", common_signatures), common_signatures)
ucell_matrix <- as.matrix(metadata[, unname(ucell_columns), drop = FALSE])
colnames(ucell_matrix) <- common_signatures
aucell_matrix <- as.matrix(metadata[, unname(aucell_columns), drop = FALSE])
colnames(aucell_matrix) <- common_signatures
clusters <- as.character(metadata[[opt$cluster_col]])

message("Writing Hallmark cluster summaries and UCell/AUCell agreement.")
ucell_cluster_medians <- cluster_median_matrix(ucell_matrix, clusters)
aucell_cluster_medians <- cluster_median_matrix(aucell_matrix, clusters)
write.csv(ucell_cluster_medians, file.path(table_dir, "hallmark_ucell_cluster_medians.csv"), row.names = TRUE, quote = TRUE)
write.csv(aucell_cluster_medians, file.path(table_dir, "hallmark_aucell_cluster_medians.csv"), row.names = TRUE, quote = TRUE)
write_heatmap(zscore_columns(ucell_cluster_medians), file.path(figure_dir, "hallmark_ucell_cluster_medians_zscore.pdf"),
              "Hallmark cluster medians: UCell (z-score by signature)")
write_heatmap(zscore_columns(aucell_cluster_medians), file.path(figure_dir, "hallmark_aucell_cluster_medians_zscore.pdf"),
              "Hallmark cluster medians: AUCell (z-score by signature)")

agreement <- data.frame(
  signature = common_signatures,
  ucell_column = unname(ucell_columns[common_signatures]),
  aucell_column = unname(aucell_columns[common_signatures]),
  spearman_rho = vapply(common_signatures, function(signature) {
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

# A cluster label is categorical, so Pearson/Spearman correlation with a
# numeric cell-cycle score is not meaningful.  Quantify score--cluster
# association with a distribution-free Kruskal-Wallis test and its rank-based
# effect sizes.  Quantify Phase--cluster association with chi-square and
# Cramer's V; standardized residuals identify phase enrichment/depletion in
# individual clusters.  Cell-level p values are descriptive because cells from
# a single library are not independent biological replicates.
kruskal_association <- function(score, cluster, score_type) {
  valid <- is.finite(score) & !is.na(cluster)
  score <- score[valid]
  cluster <- factor(cluster[valid], levels = sort(unique(as.character(cluster[valid]))))
  test <- stats::kruskal.test(score ~ cluster)
  n <- length(score)
  k <- nlevels(cluster)
  statistic <- unname(test$statistic)
  eta_squared_h <- if (n > k) max(0, min(1, (statistic - k + 1) / (n - k))) else NA_real_
  rank_epsilon_squared <- if (n > 1) max(0, min(1, statistic / (n - 1))) else NA_real_
  data.frame(
    score_type = score_type,
    n_cells = n,
    n_clusters = k,
    kruskal_wallis_chi_sq = statistic,
    degrees_freedom = unname(test$parameter),
    p_value = test$p.value,
    rank_eta_squared_H = eta_squared_h,
    rank_epsilon_squared = rank_epsilon_squared,
    stringsAsFactors = FALSE
  )
}

cycle_score_association <- do.call(rbind, list(
  kruskal_association(metadata$S.Score, clusters, "S.Score"),
  kruskal_association(metadata$G2M.Score, clusters, "G2M.Score"),
  kruskal_association(metadata$G1.RelativeScore, clusters, "G1.RelativeScore")
))
cycle_score_association$p_value_holm <- stats::p.adjust(cycle_score_association$p_value, method = "holm")
cycle_score_association$p_value_display <- ifelse(
  cycle_score_association$p_value == 0,
  "< 2.2e-16",
  format.pval(cycle_score_association$p_value, eps = 2.2e-16, digits = 3)
)
write.csv(cycle_score_association, file.path(table_dir, "cell_cycle_score_cluster_association.csv"), row.names = FALSE, quote = TRUE)

phase_cluster_table <- table(
  cluster = factor(clusters, levels = sort(unique(clusters))),
  phase = factor(as.character(metadata$Phase), levels = phase_levels)
)
phase_cluster_test <- stats::chisq.test(phase_cluster_table, correct = FALSE)
phase_cluster_n <- sum(phase_cluster_table)
phase_cramers_v <- sqrt(unname(phase_cluster_test$statistic) /
  (phase_cluster_n * min(nrow(phase_cluster_table) - 1, ncol(phase_cluster_table) - 1)))
phase_cluster_association <- data.frame(
  n_cells = phase_cluster_n,
  n_clusters = nrow(phase_cluster_table),
  n_phases = ncol(phase_cluster_table),
  chi_square = unname(phase_cluster_test$statistic),
  degrees_freedom = unname(phase_cluster_test$parameter),
  p_value = phase_cluster_test$p.value,
  cramers_v = phase_cramers_v,
  stringsAsFactors = FALSE
)
phase_cluster_association$p_value_display <- ifelse(
  phase_cluster_association$p_value == 0,
  "< 2.2e-16",
  format.pval(phase_cluster_association$p_value, eps = 2.2e-16, digits = 3)
)
write.csv(phase_cluster_association, file.path(table_dir, "cell_cycle_phase_cluster_association.csv"), row.names = FALSE, quote = TRUE)
phase_cluster_residuals <- as.data.frame(as.table(phase_cluster_test$stdres), stringsAsFactors = FALSE)
colnames(phase_cluster_residuals) <- c("cluster", "phase", "standardized_residual")
phase_cluster_residuals$observed_cell_count <- as.vector(phase_cluster_test$observed)
phase_cluster_residuals$expected_cell_count <- as.vector(phase_cluster_test$expected)
write.csv(phase_cluster_residuals, file.path(table_dir, "cell_cycle_phase_cluster_standardized_residuals.csv"), row.names = FALSE, quote = TRUE)

embedding <- Seurat::Embeddings(object, reduction = "umap")
embedding <- embedding[rownames(metadata), 1:2, drop = FALSE]

# =============================================================================
# 追加：functional_state_annotated_umap.pdf 绘制
# 规则：
# 1. 复用前面已经计算完成的每个 cluster 的 Hallmark UCell 中位数矩阵
#    ucell_cluster_medians，不在本段重复计算 cluster 中位数。
# 2. 对每个 Hallmark 在不同 cluster 间进行 z-score 标准化。
# 3. 每个 cluster 选择 z-score 最高的 Hallmark。
# 4. 将其作为该 cluster 的相对占优功能状态显示在 UMAP 上。
# =============================================================================
ucell_cluster_z <- zscore_columns(ucell_cluster_medians)

# 每个 cluster 选择 z-score 最高的 Hallmark，作为相对占优功能状态。
dominant_hallmark <- vapply(
  rownames(ucell_cluster_z),
  function(cluster_id) {
    scores <- ucell_cluster_z[cluster_id, ]
    colnames(ucell_cluster_z)[which.max(scores)]
  },
  character(1)
)
names(dominant_hallmark) <- rownames(ucell_cluster_z)

# 将 cluster 级功能状态映射回每个细胞，并把名称转换为仅首字母大写。
format_hallmark <- function(x) {
  label <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x)))
  paste0(toupper(substr(label, 1, 1)), substr(label, 2, nchar(label)))
}
cell_functional_state <- unname(dominant_hallmark[clusters])
short_functional_state <- format_hallmark(cell_functional_state)
functional_plot_data <- data.frame(
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  cluster = clusters,
  functional_state = short_functional_state,
  stringsAsFactors = FALSE
)
functional_centroids <- aggregate(
  cbind(UMAP_1, UMAP_2) ~ cluster,
  data = functional_plot_data,
  FUN = median
)
functional_centroids$functional_state <- format_hallmark(
  unname(dominant_hallmark[functional_centroids$cluster])
)

# 在 cluster 中位坐标附近放置标签；长名称自动换行，标签自动避让。
functional_centroids$label <- vapply(
  seq_len(nrow(functional_centroids)),
  function(index) {
    wrapped_state <- paste(strwrap(functional_centroids$functional_state[index], width = 24), collapse = "\n")
    paste0("C", functional_centroids$cluster[index], "\n", wrapped_state)
  },
  character(1)
)

# 生成传统最终注释图形式的功能状态 UMAP。
functional_state_umap <- ggplot(
  functional_plot_data,
  aes(x = UMAP_1, y = UMAP_2, color = functional_state)
) +
  geom_point(size = 0.28, alpha = 0.75) +
  ggrepel::geom_text_repel(
    data = functional_centroids,
    aes(x = UMAP_1, y = UMAP_2, label = label),
    inherit.aes = FALSE,
    seed = opt$seed,
    size = 2.7,
    fontface = "bold",
    box.padding = 0.65,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.color = "grey55",
    segment.size = 0.3,
    max.overlaps = Inf,
    max.iter = 20000,
    force = 1.5
  ) +
  scale_x_continuous(expand = expansion(mult = 0.14)) +
  scale_y_continuous(expand = expansion(mult = 0.14)) +
  coord_equal(clip = "off") +
  labs(
    title = "Functional-state annotated UMAP",
    subtitle = "Cluster labels show the highest relative Hallmark UCell z-score",
    color = "Dominant hallmark"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    plot.margin = margin(12, 18, 12, 18)
  )
ggsave(
  file.path(figure_dir, "functional_state_annotated_umap.pdf"),
  functional_state_umap,
  width = 12.5,
  height = 9
)
# =============================================================================

message("Writing cell-cycle figures.")
cycle_umap_plots <- list(
  make_score_umap(embedding, metadata$S.Score, "S phase score"),
  make_score_umap(embedding, metadata$G2M.Score, "G2/M phase score"),
  make_group_umap(embedding, metadata$Phase, "Cell-cycle phase", "Phase")
)
save_plot_pages(cycle_umap_plots, file.path(figure_dir, "cell_cycle_umap.pdf"), width = 7, height = 6)

cycle_long <- rbind(
  data.frame(cluster = clusters, score_type = "S.Score", score = metadata$S.Score),
  data.frame(cluster = clusters, score_type = "G2M.Score", score = metadata$G2M.Score),
  data.frame(cluster = clusters, score_type = "G1.RelativeScore", score = metadata$G1.RelativeScore)
)
cycle_long$score_type <- factor(
  cycle_long$score_type,
  levels = c("S.Score", "G2M.Score", "G1.RelativeScore"),
  labels = c("S phase score", "G2/M phase score", "G1-relative score*")
)
cycle_violin_plot <- ggplot(cycle_long, aes(x = cluster, y = score, fill = cluster)) +
  geom_violin(scale = "width", trim = TRUE, color = NA) +
  facet_wrap(~score_type, scales = "free_y") +
  labs(
    x = "Cluster", y = "Score", title = "Cell-cycle scores by cluster",
    subtitle = "* G1-relative score = -max(S phase score, G2/M phase score); higher values indicate lower S/G2M activity"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none", plot.subtitle = element_text(size = 8))
phase_bar_plot <- ggplot(phase_table, aes(x = cluster, y = proportion, fill = phase)) +
  geom_col(width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cluster", y = "Proportion of cells", fill = "Cell-cycle phase", title = "Cell-cycle phase composition by cluster") +
  theme_classic(base_size = 11)
save_plot_pages(list(cycle_violin_plot, phase_bar_plot), file.path(figure_dir, "cell_cycle_by_cluster.pdf"), width = 10, height = 6)

selected_hallmarks <- intersect(
  c(
    "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS", "HALLMARK_HYPOXIA",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_APOPTOSIS"
  ), common_signatures
)
if (length(selected_hallmarks) > 0) {
  hallmark_umap_plots <- list()
  for (signature in selected_hallmarks) {
    hallmark_umap_plots[[length(hallmark_umap_plots) + 1]] <-
      make_score_umap(embedding, ucell_matrix[, signature], paste(signature, "UCell"))
    hallmark_umap_plots[[length(hallmark_umap_plots) + 1]] <-
      make_score_umap(embedding, aucell_matrix[, signature], paste(signature, "AUCell"))
  }
  save_plot_pages(hallmark_umap_plots, file.path(figure_dir, "selected_hallmark_ucell_aucell_umap.pdf"), width = 7, height = 6)
}

scatter_hallmarks <- head(agreement$signature[order(agreement$spearman_rho, decreasing = TRUE, na.last = TRUE)], 8)
agreement_plots <- lapply(scatter_hallmarks, function(signature) {
  plot_data <- data.frame(UCell = ucell_matrix[, signature], AUCell = aucell_matrix[, signature])
  ggplot(plot_data, aes(x = UCell, y = AUCell)) +
    geom_point(size = 0.25, alpha = 0.12, color = "#2166AC") +
    labs(title = signature, subtitle = paste0("Spearman rho = ", round(agreement$spearman_rho[agreement$signature == signature], 3))) +
    theme_classic(base_size = 10)
})
if (length(agreement_plots) > 0) {
  save_plot_pages(agreement_plots, file.path(figure_dir, "ucell_aucell_agreement_scatter.pdf"), width = 6, height = 5)
}

cell_level_columns <- unique(c(opt$cluster_col, "S.Score", "G2M.Score", "G1.RelativeScore", "Phase", unname(ucell_columns), unname(aucell_columns)))
cell_level_scores <- data.frame(cell = rownames(metadata), metadata[, cell_level_columns, drop = FALSE], check.names = FALSE)
write_csv_gz(cell_level_scores, file.path(table_dir, "cell_level_functional_and_cell_cycle_scores.csv.gz"))

counts <- tryCatch(
  SeuratObject::GetAssayData(object, assay = SeuratObject::DefaultAssay(object), layer = "counts"),
  error = function(e) SeuratObject::GetAssayData(object, assay = SeuratObject::DefaultAssay(object), slot = "counts")
)
detected <- Matrix::colSums(counts > 0)
parameters <- data.frame(
  parameter = c("input_scored_rds", "species", "db_species", "cluster_col", "n_cells", "n_genes",
                "median_detected_genes", "q05_detected_genes", "ucell_max_rank", "aucell_max_rank", "hallmark_sets_scored", "ncores", "min_signature_genes", "min_cycle_genes", "seed"),
  value = c(normalizePath(opt$input_rds, winslash = "/", mustWork = TRUE), opt$species, opt$db_species, opt$cluster_col,
            ncol(object), nrow(object), round(stats::median(detected)), floor(stats::quantile(detected, 0.05)),
            opt$ucell_max_rank, opt$aucell_max_rank, length(common_signatures), opt$ncores, opt$min_signature_genes, opt$min_cycle_genes, opt$seed),
  stringsAsFactors = FALSE
)
write.csv(parameters, file.path(table_dir, "scoring_parameters.csv"), row.names = FALSE, quote = TRUE)
packages <- c("Seurat", "UCell", "AUCell", "pheatmap", "optparse")
package_versions <- data.frame(
  package = packages,
  version = vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) as.character(utils::packageVersion(package)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
write.csv(package_versions, file.path(table_dir, "package_versions.csv"), row.names = FALSE)

message("Report completed successfully.")
