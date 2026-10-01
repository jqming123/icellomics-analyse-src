#!/usr/bin/env Rscript

# =============================================================================
# scRNA_downstream_cycle_analysis.R
# 功能：承接 Cell Ranger matrix，完成 QC、Scrublet、Harmony、
#       降维聚类和差异分析；随后调用配套脚本完成 Hallmark UCell/AUCell
#       功能状态与细胞周期分析（functional-state profiling）。
# =============================================================================

# ------------------------------ 参数解析 -------------------------------------
.libPaths(c("/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/single_cell_1/lib/R/library", .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(future)
  library(reticulate)
  library(pheatmap)
  library(reshape2)
  library(gridExtra)
  library(RColorBrewer)
  library(viridis)
})

option_list <- list(
  make_option(c("--project_path"), type = "character", default = NULL,
              help = "项目根目录 [必需]"),
  make_option(c("--sample_list"), type = "character", default = NULL,
              help = "指定样本ID，逗号分隔 [可选]"),
  make_option(c("--out_dir"), type = "character", default = "annotation/preprocessing",
              help = "输出目录（相对于项目根目录） [默认 %default]"),
  make_option(c("--resolution"), type = "numeric", default = 0.8,
              help = "聚类分辨率 [默认 %default]"),
  make_option(c("--nfeatures"), type = "integer", default = 2000,
              help = "高变基因数目 [默认 %default]"),
  make_option(c("--seed"), type = "integer", default = 1234,
              help = "Global random seed [default %default]"),
  make_option(c("--nFeature_min"), type = "integer", default = 200,
              help = "最小基因数 [默认 %default]"),
  make_option(c("--nFeature_max"), type = "integer", default = 8000,
              help = "最大基因数 [默认 %default]"),
  make_option(c("--nCount_min"), type = "integer", default = 500,
              help = "最小 UMI 数 [默认 %default]"),
  make_option(c("--nCount_max"), type = "integer", default = 50000,
              help = "最大 UMI 数 [默认 %default]"),
  make_option(c("--mt_threshold"), type = "integer", default = 15,
              help = "线粒体基因百分比阈值 [默认 %default]"),
  make_option(c("--mt_pattern"), type = "character", default = "^(MT-|mt-)",
              help = "线粒体基因识别正则表达式 [默认 %default]"),
  make_option(c("--run_scrublet"), type = "logical", default = TRUE,
              help = "是否运行 Scrublet 双细胞检测 [默认 %default]"),
  make_option(c("--convert_to_h5ad"), type = "logical", default = FALSE,
              help = "是否转换为 h5ad 格式 [默认 %default]"),
  make_option(c("--annotation_script_dir"), type = "character", default = NULL,
              help = "功能状态分析脚本所在目录 [默认：本脚本所在目录]"),
  make_option(c("--annotation_output_dir"), type = "character", default = "annotation/functional_state_results",
              help = "功能状态分析输出目录（相对于项目根目录） [默认 %default]"),
  make_option(c("--annotation_ncores"), type = "integer", default = 12,
              help = "UCell/AUCell 使用的 CPU 核数 [默认 %default]"),
  make_option(c("--run_functional_annotation"), type = "logical", default = TRUE,
              help = "是否运行功能状态与细胞周期分析 [默认 %default]"),
  make_option(c("--species"), type = "character", default = "Mus musculus",
              help = "传给 msigdbr 的物种名称 [默认 %default]"),
  make_option(c("--db_species"), type = "character", default = "MM",
              help = "MSigDB 数据库物种（小鼠为 MM） [默认 %default]")
)

opt_parser <- OptionParser(option_list = option_list, description = "通用 scRNA-seq QC/聚类 + functional-state profiling")
opt <- parse_args(opt_parser)
set.seed(opt$seed)
cat("Global random seed:", opt$seed, "\n")

if (is.null(opt$project_path)) {
  print_help(opt_parser)
  stop("错误：必须提供 --project_path", call. = FALSE)
}
if (opt$annotation_ncores < 1) {
  stop("错误：--annotation_ncores 必须 >= 1", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
current_script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
annotation_script_dir <- if (is.null(opt$annotation_script_dir) || !nzchar(opt$annotation_script_dir)) {
  current_script_dir
} else {
  normalizePath(opt$annotation_script_dir, winslash = "/", mustWork = TRUE)
}

# ------------------------------ 全局设置 -------------------------------------

options(future.globals.maxSize = 16000 * 1024^2)

out_dir <- file.path(opt$project_path, opt$out_dir)
fig_dir <- file.path(out_dir, "figures")
rds_dir <- file.path(out_dir, "seurat_objects")
tab_dir <- file.path(out_dir, "tables")
doublet_dir <- file.path(out_dir, "doublet_detection")
for (d in c(out_dir, fig_dir, rds_dir, tab_dir, doublet_dir,
            file.path(tab_dir, "pca"), file.path(tab_dir, "umap"),
            file.path(tab_dir, "tsne"), file.path(tab_dir, "clustering"),
            file.path(tab_dir, "diffexp"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

qc_stats_file <- file.path(tab_dir, "qc_stats.csv")
if (file.exists(qc_stats_file)) {
  file.remove(qc_stats_file)
}

if (opt$run_scrublet || opt$convert_to_h5ad) {
  use_python("/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/single_cell_1/bin/python", required = TRUE)
}

start_time <- Sys.time()
cat("\n", strrep("*", 70), "\n")
cat("Start time: ", format(start_time), "\n")
cat("Project path: ", opt$project_path, "\n")
cat("Output directory: ", out_dir, "\n")
cat(strrep("*", 70), "\n\n")

# =============================================================================
# 辅助函数
# =============================================================================

# 读取10X矩阵（支持.gz）
read_10x_gz <- function(data_dir) {
  if (!dir.exists(data_dir)) {
    warning("目录不存在，跳过: ", data_dir)
    return(NULL)
  }
  matrix_file <- file.path(data_dir, "matrix.mtx.gz")
  barcodes_file <- file.path(data_dir, "barcodes.tsv.gz")
  features_file <- file.path(data_dir, "features.tsv.gz")
  if (!all(file.exists(matrix_file, barcodes_file, features_file))) {
    warning("文件不完整: ", data_dir)
    return(NULL)
  }
  counts <- readMM(gzfile(matrix_file))
  barcodes <- read.table(gzfile(barcodes_file), header = FALSE, stringsAsFactors = FALSE)[,1]
  features <- read.table(gzfile(features_file), header = FALSE, stringsAsFactors = FALSE)
  if (ncol(features) >= 2) rownames(counts) <- features[,2] else rownames(counts) <- features[,1]
  colnames(counts) <- barcodes
  rownames(counts) <- make.unique(rownames(counts))
  return(counts)
}

# 获取样本列表
get_sample_ids <- function(project_path) {
  expr_dir <- file.path(project_path, "3_expression_result")
  if (!dir.exists(expr_dir)) stop("未找到 3_expression_result 目录: ", expr_dir)
  if (!is.null(opt$sample_list)) {
    samples <- trimws(unlist(strsplit(opt$sample_list, ",")))
  } else {
    samples <- list.dirs(expr_dir, full.names = FALSE, recursive = FALSE)
    samples <- samples[samples != ""]
  }
  if (length(samples) == 0) stop("未找到任何样本")
  cat("检测到样本: ", paste(samples, collapse = ", "), "\n")
  return(samples)
}

# 基础质量控制与过滤
qc_filter <- function(seurat_obj, sample_id) {
  cat("  计算基础质控指标（下限 + 线粒体）...\n")
  seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = opt$mt_pattern)
  seurat_obj$percent.ribo <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA > opt$nFeature_min &
                         nCount_RNA > opt$nCount_min &
                         percent.mt < opt$mt_threshold)
  after <- ncol(seurat_obj)
  cat("    过滤前细胞数: ", before, ", 过滤后: ", after,
      " (保留率: ", round(after/before*100, 1), "%)\n", sep = "")

  thresholds <- list(
    nFeature_min = opt$nFeature_min, nFeature_max = NA_real_,
    nCount_min = opt$nCount_min, nCount_max = NA_real_,
    mt_threshold = opt$mt_threshold
  )
  save_qc_stats(sample_id, before, after, thresholds, stage = "basic_qc")

  return(seurat_obj)
}

upper_qc_filter <- function(seurat_obj, sample_id) {
  before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA <= opt$nFeature_max &
                         nCount_RNA <= opt$nCount_max)
  after <- ncol(seurat_obj)
  cat("  上界 QC（Scrublet 之后）：过滤前细胞数: ", before, ", 过滤后: ", after,
      " (保留率: ", round(after/before*100, 1), "%)\n", sep = "")

  thresholds <- list(
    nFeature_min = NA_real_, nFeature_max = opt$nFeature_max,
    nCount_min = NA_real_, nCount_max = opt$nCount_max,
    mt_threshold = NA_real_
  )
  save_qc_stats(sample_id, before, after, thresholds, stage = "upper_qc")

  return(seurat_obj)
}

# Scrublet 双细胞检测
run_scrublet <- function(seurat_obj, sample_id, doublet_dir) {

  if (!opt$run_scrublet) {
    cat("  双细胞检测未启用（--run_scrublet FALSE）\n")
    return(seurat_obj)
  }

  cat("  运行 Scrublet 双细胞检测...\n")
  suppressPackageStartupMessages(library(reticulate))

  # ---------------------------------------------------------------------------
  # 1. 检查 Python Scrublet / scipy
  # ---------------------------------------------------------------------------
  scrublet <- tryCatch(
    import("scrublet"),
    error = function(e) {
      stop("Scrublet 模块加载失败，分析终止: ", e$message)
    }
  )

  scipy_io <- tryCatch(
    import("scipy.io"),
    error = function(e) {
      stop("scipy.io 加载失败，分析终止: ", e$message)
    }
  )

  # ---------------------------------------------------------------------------
  # 2. 提取原始 counts
  #    Seurat: gene × cell
  #    Scrublet: cell × gene
  # ---------------------------------------------------------------------------
  counts <- GetAssayData(
    seurat_obj,
    assay = "RNA",
    layer = "counts"
  )

  cat(
    "    Scrublet 输入：",
    ncol(counts), " cells × ",
    nrow(counts), " genes\n",
    sep = ""
  )

  # 用 Matrix Market 在 R -> Python 间传递 sparse matrix（保持稀疏，不做 as.matrix）
  tmp_mtx <- tempfile(fileext = ".mtx")

  Matrix::writeMM(counts, file = tmp_mtx)

  # ---------------------------------------------------------------------------
  # 3. 读入 / 转置 / 运行 Scrublet —— 全部在 Python 侧完成
  # ---------------------------------------------------------------------------
  py_run_string(sprintf("
import scipy.io
import scrublet as scr

_mtx = scipy.io.mmread(r'%s').T.tocsr()

_scrub = scr.Scrublet(
    counts_matrix=_mtx,
    random_state=%d,
    # n_neighbors 不指定：由 Scrublet 按细胞数自动取 round(0.5*sqrt(n_cells))
    # expected_doublet_rate 使用 Scrublet 官方默认 0.1
)

# 其余参数（sim_doublet_ratio=2、n_neighbors 自动选择，以及 scrub_doublets 的
# min_counts=3 / min_cells=3 / min_gene_variability_pctl=85 / n_prin_comps=30）
# 均为 Scrublet 官方默认值，此处不显式传参，直接复用默认行为
scrublet_scores, scrublet_pred = _scrub.scrub_doublets()

if scrublet_pred is None:
    raise RuntimeError(
        'Scrublet failed to determine an automatic doublet threshold.'
    )

scrublet_threshold = float(_scrub.threshold_)
", tmp_mtx, opt$seed))

  unlink(tmp_mtx)

  # ---------------------------------------------------------------------------
  # 4. 把 Scrublet 结果放回 Seurat（只回传数值向量/标量，矩阵不跨边界）
  # ---------------------------------------------------------------------------
  doublet_scores <- as.numeric(py$scrublet_scores)
  predicted_doublets <- as.logical(py$scrublet_pred)

  if (length(doublet_scores) != ncol(seurat_obj)) {
    stop(
      "Scrublet 返回细胞数与 Seurat 对象不一致：",
      length(doublet_scores), " vs ", ncol(seurat_obj)
    )
  }

  seurat_obj$doublet_score <- doublet_scores
  seurat_obj$predicted_doublet <- predicted_doublets

  threshold <- as.numeric(py$scrublet_threshold)
  seurat_obj$doublet_threshold <- threshold

  # ---------------------------------------------------------------------------
  # 5. 先保存“所有细胞”的 Scrublet 结果
  # ---------------------------------------------------------------------------
  doublet_df <- data.frame(
    cell = colnames(seurat_obj),
    doublet_score = seurat_obj$doublet_score,
    predicted_doublet = seurat_obj$predicted_doublet,
    nCount_RNA = seurat_obj$nCount_RNA,
    nFeature_RNA = seurat_obj$nFeature_RNA,
    doublet_threshold = threshold,
    stringsAsFactors = FALSE
  )

  write.csv(
    doublet_df,
    file = file.path(
      doublet_dir,
      paste0(sample_id, "_doublet_info.csv")
    ),
    row.names = FALSE,
    quote = FALSE
  )

  # ---------------------------------------------------------------------------
  # 6. 先用全部细胞画图
  # ---------------------------------------------------------------------------
  meta_all <- seurat_obj@meta.data

  p1 <- ggplot(meta_all, aes(x = doublet_score)) +
    geom_histogram(
      bins = 50,
      fill = "steelblue",
      alpha = 0.7
    ) +
    geom_vline(
      xintercept = threshold,
      color = "red",
      linetype = "dashed"
    ) +
    labs(
      title = paste(sample_id, "Scrublet doublet score"),
      x = "Doublet score",
      y = "Cell count"
    ) +
    theme_bw()

  p2 <- ggplot(
    meta_all,
    aes(
      x = nCount_RNA,
      y = doublet_score,
      color = predicted_doublet
    )
  ) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(
      values = c(
        "FALSE" = "gray70",
        "TRUE" = "red"
      )
    ) +
    geom_hline(
      yintercept = threshold,
      color = "red",
      linetype = "dashed"
    ) +
    theme_bw()

  p3 <- ggplot(
    meta_all,
    aes(
      x = nFeature_RNA,
      y = doublet_score,
      color = predicted_doublet
    )
  ) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(
      values = c(
        "FALSE" = "gray70",
        "TRUE" = "red"
      )
    ) +
    geom_hline(
      yintercept = threshold,
      color = "red",
      linetype = "dashed"
    ) +
    theme_bw()

  combined <- (p1 + p2 + p3) +
    plot_annotation(
      title = paste(sample_id, "Scrublet results")
    )

  out_file <- file.path(
    doublet_dir,
    paste0("doublet_plot_", sample_id, ".pdf")
  )

  ggsave(
    out_file,
    combined,
    width = 15,
    height = 5
  )

  # ---------------------------------------------------------------------------
  # 7. 保存汇总统计
  # ---------------------------------------------------------------------------
  before <- ncol(seurat_obj)
  n_doublet <- sum(
    seurat_obj$predicted_doublet,
    na.rm = TRUE
  )

  cat(
    "    检测到 doublets: ",
    n_doublet,
    " / ",
    before,
    " (",
    round(n_doublet / before * 100, 2),
    "%)\n",
    sep = ""
  )

  # ---------------------------------------------------------------------------
  # 8. 删除 doublets
  # ---------------------------------------------------------------------------
  seurat_obj <- subset(
    seurat_obj,
    subset = predicted_doublet == FALSE
  )

  after <- ncol(seurat_obj)

  cat(
    "    双细胞过滤前: ",
    before,
    ", 过滤后: ",
    after,
    " (保留率: ",
    round(after / before * 100, 1),
    "%)\n",
    sep = ""
  )

  cat(
    "  ✓ Scrublet 结果表: ",
    file.path(
      doublet_dir,
      paste0(sample_id, "_doublet_info.csv")
    ),
    "\n"
  )

  cat(
    "  ✓ Scrublet 图: ",
    out_file,
    "\n"
  )

  return(seurat_obj)
}

# 预处理：标准化 + 高变基因
preprocess <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                     nfeatures = opt$nfeatures, verbose = FALSE)
  return(seurat_obj)
}

# 绘制质控图
qc_plots <- function(seurat_obj, sample_id) {
  meta <- seurat_obj@meta.data
  p1 <- ggplot(meta, aes(x = sample_id, y = nFeature_RNA)) +
    geom_violin(fill = "#E69F00", alpha = 0.7) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "nFeature_RNA", x = NULL, y = "Gene count")
  p2 <- ggplot(meta, aes(x = sample_id, y = nCount_RNA)) +
    geom_violin(fill = "#56B4E9", alpha = 0.7) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "nCount_RNA", x = NULL, y = "UMI count")
  p3 <- ggplot(meta, aes(x = sample_id, y = percent.mt)) +
    geom_violin(fill = "#009E73", alpha = 0.7) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "percent.mt", x = NULL, y = "Mitochondrial %")
  
  p4 <- ggplot(meta, aes(x = sample_id, y = percent.ribo)) +
    geom_violin(fill = "#CC79A7", alpha = 0.7) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "percent.ribo", x = NULL, y = "Ribosomal %")
  
  p5 <- ggplot(meta, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_viridis_c() +
    geom_hline(yintercept = c(opt$nFeature_min, opt$nFeature_max), color = "red", linetype = "dashed") +
    geom_vline(xintercept = c(opt$nCount_min, opt$nCount_max), color = "blue", linetype = "dashed") +
    theme_bw() +
    labs(title = "nCount vs nFeature (thresholds: red=genes, blue=UMI)")
  
  p6 <- ggplot(meta, aes(x = nCount_RNA, y = percent.mt, color = nFeature_RNA)) +
    geom_point(alpha = 0.6, size = 0.5) +
    scale_color_viridis_c() +
    geom_hline(yintercept = opt$mt_threshold, color = "green", linetype = "dashed") +
    theme_bw() +
    labs(title = "nCount vs percent.mt (green line: mt threshold)")
  
  top_row <- (p1 | p2 | p3 | p4) + plot_annotation(title = paste("QC -", sample_id))
  bottom_row <- (p5 | p6)
  p <- top_row / bottom_row + plot_layout(heights = c(1, 1))
  
  out_file <- file.path(fig_dir, paste0("qc_", sample_id, ".pdf"))
  ggsave(out_file, p, width = 16, height = 10)
  cat("  ✓ 质控组合图已保存: ", out_file, "\n")
}

add_ribo_percent <- function(seurat_obj, pattern = "^RP[SL]") {
  seurat_obj$percent.ribo <- PercentageFeatureSet(seurat_obj, pattern = pattern)
  return(seurat_obj)
}

save_qc_stats <- function(sample_id, before, after, thresholds, stage) {
  # 容错：基础 QC / 上界 QC 各自只填用到的阈值，未用到的记为 NA
  get_th <- function(k) {
    v <- thresholds[[k]]
    if (is.null(v) || (length(v) == 1L && is.na(v))) NA else v
  }
  stats_df <- data.frame(
    sample = sample_id,
    stage = stage,
    cells_before = before,
    cells_after = after,
    retention_rate = round(after/before*100, 2),
    nFeature_min = get_th("nFeature_min"),
    nFeature_max = get_th("nFeature_max"),
    nCount_min = get_th("nCount_min"),
    nCount_max = get_th("nCount_max"),
    mt_threshold = get_th("mt_threshold")
  )
  stats_file <- file.path(tab_dir, "qc_stats.csv")
  if (!file.exists(stats_file)) {
    write.csv(stats_df, stats_file, row.names = FALSE, quote = FALSE)
  } else {
    existing <- read.csv(stats_file)
    write.csv(rbind(existing, stats_df), stats_file, row.names = FALSE, quote = FALSE)
  }
  cat("  ✓ 质控统计表已保存: ", stats_file, "\n")
}

# =============================================================================
# 主流程
# =============================================================================
cat(">>> 1. 读取样本并创建 Seurat 对象\n")
samples <- get_sample_ids(opt$project_path)
seurat_list <- list()

for (sid in samples) {
  cat("  处理样本: ", sid, "\n")
  matrix_dir <- file.path(opt$project_path, "3_expression_result", sid, "outs", "filtered_feature_bc_matrix")
  counts <- read_10x_gz(matrix_dir)
  if (is.null(counts)) next
  seurat <- CreateSeuratObject(counts = counts, project = sid,
                               min.cells = 3, min.features = 200)
  seurat$sample_id <- sid
  seurat$dataset <- basename(opt$project_path)

  saveRDS(seurat, file = file.path(rds_dir, paste0("seurat_obj_raw_", sid, ".rds")))
  cat("    原始对象已保存: seurat_obj_raw_", sid, ".rds\n", sep = "")

  seurat <- qc_filter(seurat, sid)
  if (ncol(seurat) == 0) {
    cat("    样本 ", sid, " 基础 QC 后无细胞，跳过\n")
    next
  }
  seurat <- run_scrublet(seurat, sid, doublet_dir)
  if (ncol(seurat) == 0) {
    cat("    样本 ", sid, " Scrublet 去除 doublet 后无细胞，跳过\n")
    next
  }
  seurat <- upper_qc_filter(seurat, sid)
  if (ncol(seurat) == 0) {
    cat("    样本 ", sid, " 上界 QC 后无细胞，跳过\n")
    next
  }
  qc_plots(seurat, sid)

  saveRDS(seurat, file = file.path(rds_dir, paste0("seurat_obj_filtered_", sid, ".rds")))
  cat("    过滤后对象已保存: seurat_obj_filtered_", sid, ".rds\n", sep = "")

  seurat <- preprocess(seurat)
  seurat_list[[sid]] <- seurat
  cat("    最终保留细胞数: ", ncol(seurat), "\n")
}

if (length(seurat_list) == 0) stop("没有有效样本，退出")
cat("成功读取 ", length(seurat_list), " 个样本\n")

# ------------------------------ 整合 -----------------------------------------
cat("\n>>> 2. 数据整合\n")
library(harmony)

if (length(seurat_list) == 1) {
  integrated_obj <- seurat_list[[1]]
  cat("  仅有一个样本，跳过 Harmony 批次校正，直接进行降维聚类\n")
  DefaultAssay(integrated_obj) <- "RNA"
  integrated_obj <- NormalizeData(integrated_obj, verbose = FALSE)
  integrated_obj <- FindVariableFeatures(integrated_obj, selection.method = "vst",
                                         nfeatures = opt$nfeatures, verbose = FALSE)
  integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
  integrated_obj <- RunPCA(integrated_obj, npcs = 50, seed.use = opt$seed, verbose = FALSE)
  pca_stdev <- integrated_obj[["pca"]]@stdev
  n_elbow <- min(30, length(pca_stdev))
  x <- seq_len(n_elbow)
  y <- pca_stdev[seq_len(n_elbow)]
  x1 <- x[1]; y1 <- y[1]
  x2 <- x[n_elbow]; y2 <- y[n_elbow]
  dist_to_line <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
  selected_npcs <- which.max(dist_to_line)
  selected_npcs <- max(5, selected_npcs)
  cat("Automatically selected PCs:", selected_npcs, "\n")
  integrated_obj <- RunUMAP(integrated_obj, dims = 1:selected_npcs, seed.use = opt$seed, verbose = FALSE)
  integrated_obj <- RunTSNE(integrated_obj, dims = 1:selected_npcs, seed.use = opt$seed, verbose = FALSE)
  integrated_obj <- FindNeighbors(integrated_obj, dims = 1:selected_npcs, verbose = FALSE)
  integrated_obj <- FindClusters(integrated_obj, resolution = opt$resolution, random.seed = opt$seed, verbose = FALSE)
} else {
  integrated_obj <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list)
  )
  DefaultAssay(integrated_obj) <- "RNA"
  integrated_obj <- NormalizeData(integrated_obj, verbose = FALSE)
  integrated_obj <- FindVariableFeatures(integrated_obj, selection.method = "vst",
                                         nfeatures = opt$nfeatures, verbose = FALSE)
  integrated_obj <- ScaleData(integrated_obj, verbose = FALSE)
  integrated_obj <- RunPCA(integrated_obj, npcs = 50, seed.use = opt$seed, verbose = FALSE)
  pca_stdev <- integrated_obj[["pca"]]@stdev
  n_elbow <- min(30, length(pca_stdev))
  x <- seq_len(n_elbow)
  y <- pca_stdev[seq_len(n_elbow)]
  x1 <- x[1]; y1 <- y[1]
  x2 <- x[n_elbow]; y2 <- y[n_elbow]
  dist_to_line <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / sqrt((y2 - y1)^2 + (x2 - x1)^2)
  selected_npcs <- which.max(dist_to_line)
  selected_npcs <- max(5, selected_npcs)
  cat("Automatically selected PCs:", selected_npcs, "\n")
  set.seed(opt$seed)
  integrated_obj <- RunHarmony(integrated_obj, group.by.vars = "sample_id", theta = 2,
                               dims.use = 1:selected_npcs, verbose = FALSE)
  DefaultAssay(integrated_obj) <- "RNA"
  integrated_obj <- RunUMAP(integrated_obj, reduction = "harmony", dims = 1:selected_npcs, seed.use = opt$seed, verbose = FALSE)
  integrated_obj <- RunTSNE(integrated_obj, reduction = "harmony", dims = 1:selected_npcs, seed.use = opt$seed, verbose = FALSE)
  integrated_obj <- FindNeighbors(integrated_obj, reduction = "harmony", dims = 1:selected_npcs, verbose = FALSE)
  integrated_obj <- FindClusters(integrated_obj, resolution = opt$resolution, random.seed = opt$seed, verbose = FALSE)
}

# ------------------------------ 自动保存选中的 PC 数 -----------------------
cat("\n>>> 自动选中的 PC 数\n")
write.csv(data.frame(selected_npcs = selected_npcs, elbow_search_range = n_elbow),
          file.path(tab_dir, "pca", "selected_npcs.csv"), row.names = FALSE)

# ------------------------------ 保存 PCA 详细信息 ---------------------------
cat("\n>>> 保存 PCA 详细信息...\n")

pca_loadings <- Loadings(integrated_obj, reduction = "pca")
if (!is.null(pca_loadings) && ncol(pca_loadings) > 0) {
  pca_loadings_df <- data.frame(
    gene = rownames(pca_loadings),
    pca_loadings,
    row.names = NULL
  )
  write.csv(pca_loadings_df,
            file = file.path(tab_dir, "pca", "components.csv"),
            row.names = FALSE, quote = FALSE)
  cat("  ✓ PCA 载荷矩阵已保存: ", file.path(tab_dir, "pca", "components.csv"), "\n")
}

pca_stdev <- integrated_obj[["pca"]]@stdev
if (!is.null(pca_stdev)) {
  variance <- pca_stdev^2
  pct_variance <- variance / sum(variance) * 100
  cumulative_variance <- cumsum(pct_variance)
  pca_var_df <- data.frame(
    PC = 1:length(pca_stdev),
    standard_deviation = pca_stdev,
    variance = variance,
    percent_variance = pct_variance,
    cumulative_variance = cumulative_variance
  )
  write.csv(pca_var_df,
            file = file.path(tab_dir, "pca", "variance.csv"),
            row.names = FALSE, quote = FALSE)
  cat("  ✓ PCA 方差解释已保存: ", file.path(tab_dir, "pca", "variance.csv"), "\n")
}

elbow_plot <- ElbowPlot(integrated_obj, ndims = 50) +
  ggtitle("PCA Elbow Plot") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(fig_dir, "pca_elbow.pdf"), elbow_plot, width = 10, height = 6)
cat("  ✓ PCA 肘部图已保存: ", file.path(fig_dir, "pca_elbow.pdf"), "\n")

default_assay <- DefaultAssay(integrated_obj)

if (default_assay == "integrated") {
  hvg_assay <- "integrated"
} else {
  hvg_assay <- "RNA"
}

hvg_info <- HVFInfo(integrated_obj, assay = hvg_assay)
if (!is.null(hvg_info) && nrow(hvg_info) > 0) {
  all_genes_info <- hvg_info
  all_genes_info$gene <- rownames(all_genes_info)
  all_genes_info <- all_genes_info[, c("gene", setdiff(colnames(all_genes_info), "gene"))]
  write.csv(all_genes_info,
            file = file.path(tab_dir, "pca", "dispersion.csv"),
            row.names = FALSE, quote = FALSE)
  cat("  ✓ 所有基因离散度已保存: ", file.path(tab_dir, "pca", "dispersion.csv"), "\n")
  
  top_hvgs <- hvg_info %>%
    arrange(desc(variance.standardized)) %>%
    head(opt$nfeatures)
  top_hvgs$gene <- rownames(top_hvgs)
  top_hvgs <- top_hvgs[, c("gene", setdiff(colnames(top_hvgs), "gene"))]
  write.csv(top_hvgs,
            file = file.path(tab_dir, "pca", "features_selected.csv"),
            row.names = FALSE, quote = FALSE)
  cat("  ✓ 高变基因列表已保存: ", file.path(tab_dir, "pca", "features_selected.csv"), "\n")
}

cat("  整合后细胞数: ", ncol(integrated_obj), "\n")
cat("  聚类数: ", length(unique(integrated_obj$seurat_clusters)), "\n")

# ------------------------------ 保存降维坐标与聚类 ---------------------------
cat("\n>>> 3. 保存降维坐标和聚类结果\n")

pca_emb <- Embeddings(integrated_obj, "pca")
write.csv(data.frame(cell = rownames(pca_emb), 
                     sample_id = integrated_obj$sample_id,
                     dataset = integrated_obj$dataset,
                     pca_emb),
          file = file.path(tab_dir, "pca", "pca_coords.csv"), row.names = FALSE, quote = FALSE)
umap_emb <- Embeddings(integrated_obj, "umap")
write.csv(data.frame(cell = rownames(umap_emb), 
                     sample_id = integrated_obj$sample_id,
                     dataset = integrated_obj$dataset,
                     umap_emb),
          file = file.path(tab_dir, "umap", "umap_coords.csv"), row.names = FALSE, quote = FALSE)
tsne_emb <- Embeddings(integrated_obj, "tsne")
write.csv(data.frame(cell = rownames(tsne_emb), 
                     sample_id = integrated_obj$sample_id,
                     dataset = integrated_obj$dataset,
                     tsne_emb),
          file = file.path(tab_dir, "tsne", "tsne_coords.csv"), row.names = FALSE, quote = FALSE)

cluster_df <- data.frame(cell = colnames(integrated_obj),
                         sample = integrated_obj$sample_id,
                         seurat_cluster = integrated_obj$seurat_clusters)
write.csv(cluster_df, file = file.path(tab_dir, "clustering", "seurat_clusters.csv"),
          row.names = FALSE, quote = FALSE)

p_umap_dataset <- DimPlot(integrated_obj, reduction = "umap", group.by = "sample_id", pt.size = 0.2)
p_umap_cluster <- DimPlot(integrated_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
p_tsne <- DimPlot(integrated_obj, reduction = "tsne", group.by = "seurat_clusters", label = TRUE)
combined <- (p_umap_dataset + p_umap_cluster) / p_tsne + plot_annotation(title = "Integration Results")
ggsave(file.path(fig_dir, "integration_plots.pdf"), combined, width = 16, height = 12)

# ------------------------------ 差异表达分析 ---------------------------------
cat("\n>>> 4. 差异表达分析（FindAllMarkers）\n")
DefaultAssay(integrated_obj) <- "RNA"

if (packageVersion("Seurat") >= package_version("5.0.0")) {
  cat("  合并 RNA assay 的 layers...\n")
  integrated_obj <- JoinLayers(integrated_obj)
}

integrated_obj <- NormalizeData(integrated_obj, verbose = FALSE)

all_markers <- FindAllMarkers(integrated_obj, 
                              only.pos = TRUE, 
                              min.pct = 0.01,
                              logfc.threshold = 0.1, 
                              test.use = "wilcox", 
                              verbose = FALSE,
                              slot = "data",
                              return.thresh = 1)

if (nrow(all_markers) > 0) {
  write.csv(all_markers, 
            file = file.path(tab_dir, "diffexp", "all_markers.csv"),
            row.names = FALSE, quote = FALSE)
  
  sig_markers <- all_markers %>% filter(p_val_adj < 0.05)
  if (nrow(sig_markers) > 0) {
    write.csv(sig_markers,
              file = file.path(tab_dir, "diffexp", "significant_markers.csv"),
              row.names = FALSE, quote = FALSE)
    cat("  显著差异表达基因已保存: ", file.path(tab_dir, "diffexp", "significant_markers.csv"), "\n")
  }
  
  top5 <- sig_markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>% ungroup()
  write.csv(top5, 
            file = file.path(tab_dir, "diffexp", "top5_markers_per_cluster.csv"),
            row.names = FALSE, quote = FALSE)
  
  clusters <- unique(sig_markers$cluster)
  for (cl in clusters) {
    cl_markers <- sig_markers %>%
      filter(cluster == cl) %>%
      arrange(p_val_adj, desc(avg_log2FC)) %>%
      head(10)
    if (nrow(cl_markers) > 0) {
      write.csv(cl_markers,
                file = file.path(tab_dir, "diffexp", paste0("cluster_", cl, "_top10_markers.csv")),
                row.names = FALSE, quote = FALSE)
    }
  }
  cat("  各聚类 top10 marker 基因已保存\n")
  
  cat("  差异表达分析完成，找到 ", nrow(all_markers), " 个标记基因\n")
  
} else {
  cat("  警告：未找到差异表达基因（可能由于细胞数太少或聚类间差异小）\n")
  write.csv(data.frame(), 
            file = file.path(tab_dir, "diffexp", "all_markers.csv"), 
            row.names = FALSE)
  write.csv(data.frame(), 
            file = file.path(tab_dir, "diffexp", "top5_markers_per_cluster.csv"), 
            row.names = FALSE)
}

# =============================================================================
# 分析阶段说明
# =============================================================================
# 下方保存聚类后对象，并将其交给功能状态与细胞周期分析脚本。

# ------------------------------ 保存聚类后对象 -------------------------------
cat("\n>>> 5. 保存聚类后、功能状态分析前的 Seurat 对象\n")
pre_profiling_rds <- file.path(rds_dir, "seurat_obj_integrated_pre_profiling.rds")
saveRDS(integrated_obj, file = pre_profiling_rds, compress = "gzip")
cat("  ✓ 已保存: ", pre_profiling_rds, "\n")

# ------------------------------ 功能状态与细胞周期分析 -----------------------
if (opt$run_functional_annotation) {
  cat("\n>>> 6. 运行功能状态与细胞周期分析（functional-state profiling）\n")
  scoring_script <- file.path(annotation_script_dir, "scRNA_functional_state_scoring.R")
  report_script <- file.path(annotation_script_dir, "scRNA_functional_state_report.R")
  if (!file.exists(scoring_script)) stop("未找到功能状态评分脚本: ", scoring_script, call. = FALSE)
  if (!file.exists(report_script)) stop("未找到功能状态报告脚本: ", report_script, call. = FALSE)

  annotation_out <- if (grepl("^/", opt$annotation_output_dir)) {
    opt$annotation_output_dir
  } else {
    file.path(opt$project_path, opt$annotation_output_dir)
  }
  dir.create(annotation_out, recursive = TRUE, showWarnings = FALSE)

  rscript_executable <- file.path(R.home("bin"), "Rscript")
  if (.Platform$OS.type == "windows") {
    windows_rscript <- file.path(R.home("bin"), "Rscript.exe")
    if (file.exists(windows_rscript)) rscript_executable <- windows_rscript
  }
  scoring_args <- c(
    shQuote(scoring_script),
    "--input", shQuote(pre_profiling_rds),
    "--output-dir", shQuote(annotation_out),
    "--species", shQuote(opt$species),
    "--db-species", shQuote(opt$db_species),
    "--assay", "RNA",
    "--cluster-col", "seurat_clusters",
    "--ncores", as.character(opt$annotation_ncores),
    "--seed", as.character(opt$seed)
  )
  scoring_output <- system2(rscript_executable, args = scoring_args, stdout = TRUE, stderr = TRUE)
  scoring_status <- attr(scoring_output, "status")
  cat(paste(scoring_output, collapse = "\n"), "\n")
  if (!is.null(scoring_status) && scoring_status != 0) {
    stop("功能状态与细胞周期分析失败（退出码 ", scoring_status, "）", call. = FALSE)
  }

  scored_rds <- file.path(annotation_out, "seurat", "functional_state_scored_seurat.rds")
  if (!file.exists(scored_rds)) stop("评分脚本结束，但未找到最终 RDS: ", scored_rds, call. = FALSE)
  integrated_obj <- readRDS(scored_rds)
  cat("  ✓ 功能状态与细胞周期分析完成: ", scored_rds, "\n")
} else {
  cat("\n>>> 6. 已按参数跳过功能状态与细胞周期分析\n")
}

# ------------------------------ 转换为 h5ad（可选） -------------------------
if (opt$convert_to_h5ad) {
  cat("\n>>> 7. 转换为 h5ad 格式\n")
  suppressPackageStartupMessages(library(reticulate))
  
  h5ad_dir <- file.path(out_dir, "h5ad")
  dir.create(h5ad_dir, recursive = TRUE, showWarnings = FALSE)
  out_h5ad <- file.path(h5ad_dir, "seurat_obj_functional_state_profiled.h5ad")
  
  cat("  提取数据...\n")
  counts <- GetAssayData(integrated_obj, assay = "RNA", layer = "counts")
  
  metadata <- integrated_obj@meta.data
  for (col in names(metadata)) {
    if (is.factor(metadata[[col]])) metadata[[col]] <- as.character(metadata[[col]])
  }
  
  var_df <- data.frame(gene_name = rownames(integrated_obj), row.names = rownames(integrated_obj))
  hvgs <- VariableFeatures(integrated_obj)
  var_df$highly_variable <- rownames(var_df) %in% hvgs
  
  obsm_list <- list()
  if ("pca" %in% Reductions(integrated_obj)) {
    obsm_list[["X_pca"]] <- Embeddings(integrated_obj, "pca")
  }
  if ("umap" %in% Reductions(integrated_obj)) {
    obsm_list[["X_umap"]] <- Embeddings(integrated_obj, "umap")
  }
  if ("tsne" %in% Reductions(integrated_obj)) {
    obsm_list[["X_tsne"]] <- Embeddings(integrated_obj, "tsne")
  }
  
  tmp_mtx <- tempfile(fileext = ".mtx")
  tmp_obs <- tempfile(fileext = ".csv")
  tmp_var <- tempfile(fileext = ".csv")
  tmp_obsm <- list()
  
  Matrix::writeMM(counts, file = tmp_mtx)
  write.csv(metadata, file = tmp_obs, row.names = TRUE)
  write.csv(var_df,   file = tmp_var, row.names = TRUE)
  
  for (nm in names(obsm_list)) {
    f <- tempfile(fileext = ".csv")
    df <- as.data.frame(obsm_list[[nm]])
    df$cell_id <- rownames(df)
    df <- df[, c("cell_id", setdiff(colnames(df), "cell_id"))]
    write.csv(df, file = f, row.names = FALSE)
    tmp_obsm[[nm]] <- f
  }
  
  cat("  调用 Python 构建 h5ad...\n")
  py_run_string(sprintf("
import anndata as ad
ad.settings.allow_write_nullable_strings = True
import pandas as pd
import scipy.io
import scanpy as sc

counts = scipy.io.mmread(r'%s').T.tocsr()
obs = pd.read_csv(r'%s', index_col=0)
var = pd.read_csv(r'%s', index_col=0)

adata = ad.AnnData(X=counts, obs=obs, var=var)
adata.raw = adata.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

for name, file in r.tmp_obsm.items():
    coords = pd.read_csv(file, index_col='cell_id')
    adata.obsm[name] = coords.loc[adata.obs_names, :].values.astype(float)

adata.write(r'%s')
print('h5ad 构建完成')
", tmp_mtx, tmp_obs, tmp_var, out_h5ad))
  
  unlink(c(tmp_mtx, tmp_obs, tmp_var))
  for (f in tmp_obsm) unlink(f)
  
  if (file.exists(out_h5ad)) {
    cat("  ✓ h5ad 文件已生成：", out_h5ad, "\n")
    cat("    文件大小：", round(file.info(out_h5ad)$size / 1024^2, 2), "MB\n")
  } else {
    cat("  ⚠️ h5ad 文件生成失败\n")
  }
}

# ------------------------------ 完成 -----------------------------------------
end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")
cat("\n", strrep("*", 70), "\n")
cat("End time: ", format(end_time), "\n")
cat("Total time: ", round(total_time, 2), " minutes\n")
cat(strrep("*", 70), "\n")
cat("✅ scRNA-seq QC/聚类 + functional-state profiling 全部完成！\n")
cat("预处理输出目录: ", out_dir, "\n")
if (opt$run_functional_annotation) cat("功能状态分析输出目录: ", annotation_out, "\n")
