#!/usr/bin/env Rscript

# @File       :advance_img.R
# @Time       :2025/10/28
# @Author     :Gemini
# @Product    :RStudio
# @Version    :R 4.1.1
# @agroups_dir    :encode
# @Description:
#   本脚本用于对DEGs.R脚本生成的差异分析结果进行高级可视化。
#   它会读取差异分析的CSV文件和标准化的表达矩阵，然后生成火山图和热图。
# @Usage      : Rscript <this_script.R> <group_alias> <cl_name>
# @Example    : Rscript advance_img.R PRJNA974014_a_stressed_Day5 CHO

rm(list = ls())

######################################
# 1. 初始化和参数设置
######################################

cat("==================================================================\n")
cat("[INFO] 步骤 1: 初始化环境和设置参数 (高级可视化)...\n")
cat("==================================================================\n")

# --- 从命令行获取参数 ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("请在运行脚本时提供组别名 (group_name) 和细胞系名 (cl_name)。\n用法: Rscript advance_img.R <group_name> <cl_name>", call. = FALSE)
}
group_name <- args[1]
cl_name <- args[2]
cat(paste0("  - 接收到的组别名 (group_name): ", group_name, "\n"))
cat(paste0("  - 接收到的细胞系名 (cl_name): ", cl_name, "\n"))

# --- 动态构建变量 (路径必须与 DEGs.R 脚本中的定义一致) ---
groups_dir <- '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/groups'
outFolder <- file.path(groups_dir, cl_name, group_name)

# --- 定义输入文件的路径 (这些文件由 DEGs.R 脚本生成) ---
diff_results_file <- file.path(outFolder, "diff_genes.csv")
normalized_counts_file <- file.path(outFolder, "count_transformation_vst.csv")
condition_file <- file.path(outFolder, paste0(group_name, ".tsv")) # 需要原始分组信息用于热图标注

cat(paste0("  - 项目根目录 (groups_dir): ", groups_dir, "\n"))
cat(paste0("  - 输出/输入目录 (outFolder): ", outFolder, "\n"))
cat(paste0("  - 差异分析结果文件: ", diff_results_file, "\n"))
cat(paste0("  - 标准化表达矩阵文件: ", normalized_counts_file, "\n"))
cat(paste0("  - 分组信息文件: ", condition_file, "\n"))

# --- 检查输入文件是否存在 ---
if (!file.exists(diff_results_file) || !file.exists(normalized_counts_file) || !file.exists(condition_file)) {
  stop("一个或多个必需的输入文件不存在！请先成功运行 DEGs.R 脚本。", call. = FALSE)
}
# --- 绘图参数 (应与 DEGs.R 脚本中的设置保持一致) ---
pValue <- 0.05
fcValue <- 2
log2FC_threshold <- log2(fcValue)
cat(paste0("  - P值阈值 (pValue): ", pValue, "\n"))
cat(paste0("  - Fold Change阈值 (fcValue): ", fcValue, " (log2FC: ", round(log2FC_threshold, 2), ")\n\n"))


#####################################
# 2. 加载包和数据
#####################################
cat("==================================================================\n")
cat("[INFO] 步骤 2: 加载 R 包和已处理的数据...\n")
cat("==================================================================\n")

# 加载包
pacman::p_load(ggplot2, pheatmap, ggrepel, RColorBrewer, dplyr)
cat("  - 成功加载所有必需的 R 包。\n")

# 设置工作目录
setwd(groups_dir)
cat(paste0("  - 工作目录已设置为: ", getwd(), "\n"))

# 加载差异分析结果
res_df <- read.csv(diff_results_file, row.names = 1)
cat(paste0("  - 已加载差异分析结果。共包含 ", nrow(res_df), " 个基因。\n"))

# 加载标准化后的counts
normalized_counts <- read.csv(normalized_counts_file, row.names = 1)
cat(paste0("  - 已加载VST标准化的表达矩阵。维度: ", nrow(normalized_counts), " x ", ncol(normalized_counts), "。\n"))

# 加载并处理分组信息 (用于热图标注)
# 复制DEGs.R中的函数以确保样本名匹配
correct_IDs_func <-function(df,colname){
  if(!colname %in% colnames(df)) {
    stop(paste0("列名 ", colname, " 不存在"))
  }
  df[[colname]] <- gsub("-", ".", df[[colname]])
  df[[colname]] <- ifelse(grepl("^\\d", df[[colname]]), paste0("X", df[[colname]]), df[[colname]])
  df[df == ""] <- NA
  df <- na.omit(df) %>% distinct()
  return(df)
}
condition_df <- read.table(condition_file, header = T, fill=T, na.strings = "", sep="\t")
condition_df <- correct_IDs_func(condition_df, 'sample')
cat("  - 已加载并处理了样本分组信息。\n\n")


#####################################
# 3. 绘制火山图 (Volcano Plot)
#####################################
cat("==================================================================\n")
cat("[INFO] 步骤 3: 生成并保存火山图...\n")
cat("==================================================================\n")

# 为基因添加分类标签 (Up, Down, Not Sig)
res_df <- res_df %>%
  mutate(change = case_when(
    padj < pValue & log2FoldChange >= log2FC_threshold ~ "UP",
    padj < pValue & log2FoldChange <= -log2FC_threshold ~ "DOWN",
    TRUE ~ "NOT SIG"
  ))
cat(paste0("  - 差异基因分类统计:\n"))
print(table(res_df$change))

# 挑选出最显著的基因用于在图上标记
top_genes <- res_df %>%
  filter(padj < pValue, abs(log2FoldChange) >= log2FC_threshold) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  head(20) # 标记前20个最显著的基因

# 绘制火山图
volcano_plot <- ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = change), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("UP" = "#E64B35", "DOWN" = "#3C5488", "NOT SIG" = "grey")) +
  # 添加阈值线
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pValue), linetype = "dashed", color = "black") +
  # 添加基因标签
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  max.overlaps = Inf, size = 3, box.padding = 0.5) +
  labs(title = paste("Volcano Plot for", group_name, "in", cl_name), 
       x = "log2(Fold Change)",
       y = "-log10(Adjusted P-value)") +
  theme_bw(base_size = 14) +
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

# 保存火山图
volcano_plot_path <- file.path(outFolder, "volcano_plot.pdf")
ggsave(volcano_plot_path, volcano_plot, width = 10, height = 8)
cat(paste0("[OUTPUT] 火山图已成功保存至: ", volcano_plot_path, "\n\n"))


#####################################
# 4. 绘制热图 (Heatmap)
#####################################
cat("==================================================================\n")
cat("[INFO] 步骤 4: 生成并保存差异基因热图...\n")
cat("==================================================================\n")

# 筛选出前50个最显著的差异基因
top_degs <- res_df %>%
  filter(padj < pValue) %>%
  arrange(padj) %>%
  head(50)

if (nrow(top_degs) < 2) {
    cat(paste0("  - [WARNING] 显著差异基因数量不足 (padj < 0.05, n = ", nrow(top_degs), ")，需要至少 2 个基因才能生成热图。\n"))
} else {
    cat(paste0("  - 已筛选出前 ", nrow(top_degs), " 个最显著的差异基因用于绘制热图。\n"))

    # 从标准化的表达矩阵中提取这些基因的数据
    heatmap_data <- normalized_counts[rownames(top_degs), condition_df$sample]

    # 创建样本注释信息
    annotation_col <- data.frame(
      condition = factor(condition_df$condition),
      row.names = condition_df$sample
    )
    
    # 定义注释颜色
    ann_colors <- list(
        condition = c(control = "#8491B4", case = "#F39B7F")
    )

    # 绘制热图
    heatmap_path <- file.path(outFolder, "heatmap_top50_DEGs.pdf")
    pdf(heatmap_path, width = 10, height = 12)
    pheatmap(heatmap_data,
             scale = "row", # 对行进行Z-score标准化，突出相对变化
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             show_rownames = TRUE, # 如果基因名不乱，可以显示
             show_colnames = TRUE,
             annotation_col = annotation_col,
             annotation_colors = ann_colors,
             border_color = "white",
             fontsize_row = 8,
             main = paste("Heatmap of Top 50 DEGs for", group_name, "in", cl_name))
    dev.off() # 关闭PDF设备

    cat(paste0("[OUTPUT] 差异基因热图已成功保存至: ", heatmap_path, "\n\n"))
}


#####################################
# 5. 脚本结束
#####################################
cat("==================================================================\n")
cat("[SUCCESS] 高级可视化脚本已成功完成！\n")
cat("==================================================================\n")