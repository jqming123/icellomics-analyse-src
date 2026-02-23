#!/usr/bin/env Rscript

# @File       :DEGs
# @Time       :2022/8/27 00:35
# @Author     :ZhouBowen
# @Product    :DataSpell
# @Version    :R 4.1.1
# @Project    :encode
# @Description:
# @Usage      : Rscript <DEG_and_basic_img.R> <group_alias> <cell_line_name>
# @Example    : Rscript DEG_and_basic_img.R PRJNA974014_a_stressed_Day5 CHO

rm(list = ls())

######################################
# 1. 初始化和参数设置
######################################

cat("========================================================\n")
cat("[INFO] 步骤 1: 初始化环境和设置参数...\n")
cat("========================================================\n")

# --- 从命令行获取参数 ---
args <- commandArgs(trailingOnly = TRUE)

# 检查是否提供了至少两个参数
if (length(args) < 2) {
  stop("错误: 请提供组别名 (groupname) 和 细胞系名称 (cl_name).\n用法: Rscript DEG_and_basic_img.R <groupname> <cl_name>\n示例: Rscript this_script.R PRJNA974014_a_stressed_Day5 HEK293", call. = FALSE)
}

groupname <- args[1]
cl_name <- args[2]

cat(paste0("  - 接收到的组别名 (groupname): ", groupname, "\n"))
cat(paste0("  - 接收到的细胞系 (cl_name): ", cl_name, "\n"))

groups_dir <- "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/groups"

# --- 动态构建路径变量 ---
project <- file.path(groups_dir, cl_name)
outFolder <- file.path(project, groupname)
condition <- file.path(outFolder, paste0(groupname, ".tsv"))
expMatrix <- file.path(outFolder, paste0(groupname, ".count.tsv"))

cat(paste0("  - 项目根目录 (project): ", project, "\n"))
cat(paste0("  - 输出目录 (outFolder): ", outFolder, "\n"))
cat(paste0("  - 分组信息文件 (condition): ", condition, "\n"))
cat(paste0("  - 表达矩阵文件 (expMatrix): ", expMatrix, "\n"))


# --- 其他参数 ---
pValue <- 0.05
fcValue <- 2
minCounts <- 10
cat(paste0("  - P值阈值 (pValue): ", pValue, "\n"))
cat(paste0("  - Fold Change阈值 (fcValue): ", fcValue, "\n"))
cat(paste0("  - 最低计数过滤阈值 (minCounts): ", minCounts, "\n\n"))

#####################################
# 2. 加载包和设置环境
#####################################
cat("========================================================\n")
cat("[INFO] 步骤 2: 加载 R 包和设置工作环境...\n")
cat("========================================================\n")

# 定义需要加载的包
libs <- c("DESeq2", "pheatmap", "ggrepel", "RColorBrewer", "vsn", "ggplot2",
          "cowplot", "hexbin", "stringr", "reshape2", "ggsignif", "dplyr",
          "tidyr", "ggsci")

# 循环加载包，如果缺少包则直接报错停止，而不是尝试联网安装
for (lib in libs) {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    stop(paste0("错误: 未找到 R 包 '", lib, "'。请先在联网环境下安装该包。"), call. = FALSE)
  }
}

cat("  - 成功加载所有必需的 R 包。\n")

# 设置工作目录
setwd(project)
cat(paste0("  - 工作目录已设置为: ", getwd(), "\n"))

# 检查或创建输出目录
if (file.exists(outFolder)){
    cat(paste0("  - 输出目录 '", outFolder,"' 已存在。\n"))
}else{
    dir.create(outFolder, recursive = TRUE) # 使用 recursive = TRUE 更安全
    cat(paste0("  - 成功创建输出目录: '",outFolder,"'.\n"))
}

# --- 设置日志文件 ---
# 注意: 这个设置必须在 outFolder 创建之后
log_file_path <- file.path(outFolder, paste0(groupname, "_run.log"))

# --- 修正后的日志重定向 ---
# 使用一个 sink() 调用同时捕获标准输出和消息/错误
# split = TRUE 意味着输出也会同时显示在控制台，方便监控
# append = FALSE 确保每次运行都创建一个新的日志文件，覆盖旧的
sink(log_file_path, append = FALSE, type = c("output", "message"), split = TRUE)


# 在日志文件和控制台打印一些初始信息
cat("--- 日志开始 ---\n")
cat(paste("脚本启动时间:", Sys.time(), "\n"))
cat(paste("接收到的组别名为:", groupname, "\n"))
cat(paste("所有输出将被记录到:", log_file_path, "\n\n"))

#####################################
# 3. 数据加载和预处理
#####################################
cat("========================================================\n")
cat("[INFO] 步骤 3: 加载并预处理输入数据...\n")
cat("========================================================\n")

in_df <- read.table(expMatrix, sep="\t", row.names=1, header=T)
cat(paste0("  - 已加载表达矩阵。维度: ", nrow(in_df), " 行, ", ncol(in_df), " 列。\n"))

#去除都是0的行
rows_before <- nrow(in_df)
in_df <- in_df[which(rowSums(in_df) > 0),]
rows_after <- nrow(in_df)
cat(paste0("  - 移除了 ", rows_before - rows_after, " 个在所有样本中表达量均为0的基因。\n"))

#去除ensembl ID的版本号
rownames(in_df)<-factor(unlist(lapply(as.character(rownames(in_df)),function(x){strsplit(x, "\\.")[[1]][1]})))
cat("  - 已移除基因 Ensembl ID 的版本号 (例如, 从 ENSG00000223972.5 到 ENSG00000223972)。\n")

#数据长宽转换
condition_df <- read.table(condition,header = T,fill=T,na.strings = "",sep="\t")
cat(paste0("  - 已加载分组信息文件。包含 ", nrow(condition_df), " 个样本条目。\n"))

# 判断并修正df的sample列名
correct_IDs_func <-function(df,colname){
  if(!colname %in% colnames(df)) {
    stop(paste0("列名 ", colname, " 不存在"))
  }
  cat("  - 正在修正样本ID (将'-'替换为'.', 为数字开头的ID加'X'前缀)...\n")
  df[[colname]] <- gsub("-", ".", df[[colname]])
  df[[colname]] <- ifelse(grepl("^\\d", df[[colname]]), paste0("X", df[[colname]]), df[[colname]])
  df[df == ""] <- NA
  df <- na.omit(df) %>% distinct()
  cat("  - 样本ID修正完成。修正后的前3行分组信息如下:\n")
  print(head(df,3))
  return(df)
}
condition_df <- correct_IDs_func(condition_df,'sample')

#转化为因子
condition_df$condition <- factor(condition_df$condition)
# 确保 'control' 是参考水平
condition_df$condition <- relevel(condition_df$condition, ref = "control")
cat("  - 已将 'condition' 列转换为因子，并设置 'control' 为参考水平。\n")

#提取目标分析样本生成新矩阵
cts <- in_df[, condition_df$sample, drop = FALSE]
cat(paste0("  - 已根据分组信息文件筛选表达矩阵。最终分析矩阵维度: ", nrow(cts), " 基因 x ", ncol(cts), " 样本。\n"))
cat("  - 最终表达矩阵的前3行如下:\n")
print(head(cts,3))

# 原始数据count的分布可视化
cts_long <- log2(cts+1) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "count")

rawCount_boxplot <- ggplot(cts_long, aes(x = sample, y = count)) +
    geom_boxplot(aes(fill = sample), color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "sample", y = "log2(count+1)") +
    ggtitle("rawCount") +
    theme(legend.position = "none")
ggsave(file.path(outFolder, "rawCounts.pdf"), rawCount_boxplot, width=8, height=4)
cat(paste0("[OUTPUT] 原始 read counts 的箱线图已保存至: ", file.path(outFolder, "rawCounts.pdf"), "\n\n"))

#####################################
# 4. DESeq2 差异表达分析
#####################################
cat("========================================================\n")
cat("[INFO] 步骤 4: 执行 DESeq2 差异表达分析...\n")
cat("========================================================\n")

#构建DESeq2对象dds
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                             colData = condition_df,
                             design = ~ condition)
cat("  - DESeqDataSet 对象已成功创建。\n")
print(dds)

#数据过滤（默认>=10）
cat(paste0("  - 正在过滤低表达基因 (在所有样本中合计counts < ", minCounts, ")...\n"))
cat("  - 过滤前的维度: ")
print(dim(dds))
keep <- rowSums(counts(dds)) >= minCounts
dds <- dds[keep,]
cat("  - 过滤后的维度: ")
print(dim(dds))

cat("  - 确认实验设计中的条件水平: \n")
print(dds$condition)

#使用DESeq()函数进行差异分析流程；
cat("\n  - 正在运行 DESeq() 函数进行核心分析 (估算大小因子，基因离散度，拟合模型等)...\n")
dds <- DESeq(dds)
cat("  - DESeq() 分析完成。\n")

#使用results()函数提取分析结果；
cat("\n  - 正在提取 case vs control 的差异表达结果...\n")
res <- results(dds, contrast=c("condition","case","control"))
cat("  - 结果摘要如下:\n")
print(res)

#根据p value，对结果进行升序排列:
resOrdered <- res[order(res$pvalue),]
cat("  - 已按 p-value 对结果进行升序排序:\n")
print(resOrdered)

#统计adjusted p-values < 0.05 （默认）的基因数；
num_sig_genes <- sum(res$padj < pValue, na.rm=TRUE)
cat(paste0("  - 在 padj < ", pValue, " 的阈值下，共发现 ", num_sig_genes, " 个显著差异表达基因。\n\n"))

#####################################
# 5. 保存差异分析结果
#####################################
cat("========================================================\n")
cat("[INFO] 步骤 5: 导出差异分析结果文件...\n")
cat("========================================================\n")

#导出差异分析结果数据表格
resOrdered <- resOrdered[complete.cases(resOrdered), ]
resOrdered <- resOrdered[resOrdered$pvalue != 0, ]
write.csv(as.data.frame(resOrdered),
          file=file.path(outFolder, "diff_genes.csv"), quote = FALSE)
cat(paste0("[OUTPUT] 完整的差异基因列表已保存至: ", file.path(outFolder, "diff_genes.csv"), "\n"))

#如果只导出差异基因，可使用subset函数；
resSig <- subset(resOrdered, padj < pValue)
write.csv(as.data.frame(resSig),
          file=file.path(outFolder, paste0("diff_genes_padj", pValue, ".csv")), quote = FALSE)
cat(paste0("[OUTPUT] 显著差异基因列表 (padj < ", pValue, ") 已保存至: ", file.path(outFolder, paste0("diff_genes_padj", pValue, ".csv")), "\n\n"))


#####################################
# 6. 质控与可视化
#####################################
cat("========================================================\n")
cat("[INFO] 步骤 6: 质量控制和数据可视化...\n")
cat("========================================================\n")

# 计算异常值
pdf(file.path(outFolder, "outliers.pdf"), width=6,height=6)
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()
cat(paste0("[OUTPUT] Cook's distance 异常值检测图已保存至: ", file.path(outFolder, "outliers.pdf"), "\n"))


# 独立过滤结果
cat("  - 独立过滤的 alpha 值为: ", metadata(res)$alpha, "\n")
cat("  - 用于独立过滤的表达量阈值为:\n")
print(metadata(res)$filterThreshold)

# 根据样本数量选择不同的标准化方法
if (nrow(condition_df) <= 50) {
  cat("\n  - 样本数 (<=50)，将执行 VST, RLD 和 NTD 三种标准化方法。\n")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  ntd <- normTransform(dds)
  cat("  - 三种标准化方法均已完成。\n")

  #导出标准化后的数据；
  write.csv(as.data.frame(assay(vsd)),
            file=file.path(outFolder, "count_transformation_vst.csv"))
  cat(paste0("[OUTPUT] VST 标准化后的表达矩阵已保存。\n"))
  write.csv(as.data.frame(assay(rld)),
            file=file.path(outFolder, "count_transformation_rlog.csv"))
  cat(paste0("[OUTPUT] RLD 标准化后的表达矩阵已保存。\n"))
  write.csv(as.data.frame(assay(ntd)),
            file=file.path(outFolder, "count_transformation_normTransform.csv"))
  cat(paste0("[OUTPUT] NTD 标准化后的表达矩阵已保存。\n"))

  #标准化数据count的分布可视化
  vsd_norm_long <- as.data.frame(assay(vsd)) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "count")
  rld_norm_long <- as.data.frame(assay(rld)) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "count")
  ntd_norm_long <- as.data.frame(assay(ntd)) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "count")

  vsdCount_boxplot <- ggplot(vsd_norm_long, aes(x = sample, y = count)) +
      geom_boxplot(aes(fill = sample), color = "black") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "sample", y = "normalized count") +
      ggtitle("vsdCount") +
      theme(legend.position = "none")

  rldCount_boxplot <- ggplot(rld_norm_long, aes(x = sample, y = count)) +
      geom_boxplot(aes(fill = sample), color = "black") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "sample", y = "normalized count") +
      ggtitle("rldCount") +
      theme(legend.position = "none")

  ntdCount_boxplot <- ggplot(ntd_norm_long, aes(x = sample, y = count)) +
      geom_boxplot(aes(fill = sample), color = "black") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "sample", y = "normalized count") +
      ggtitle("ntdCount") +
      theme(legend.position = "none")

  pdf(file.path(outFolder, "normalizationCounts.pdf"), width=13,height=4)
  print(plot_grid(vsdCount_boxplot,
                      rldCount_boxplot,
                      ntdCount_boxplot,
                      ncol = 3))
  dev.off()
  cat(paste0("[OUTPUT] 三种标准化方法的效果箱线图已保存至: ", file.path(outFolder, "normalizationCounts.pdf"), "\n"))

} else {
  cat("\n  - 样本数 (>50)，为提高效率，仅执行 VST 标准化方法。\n")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  cat("  - VST 标准化已完成。\n")

  write.csv(as.data.frame(assay(vsd)),
            file=file.path(outFolder, "count_transformation_vst.csv"))
  cat(paste0("[OUTPUT] VST 标准化后的表达矩阵已保存。\n"))

  vsd_norm_long <- as.data.frame(assay(vsd)) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "count")

  vsdCount_boxplot <- ggplot(vsd_norm_long, aes(x = sample, y = count)) +
    geom_boxplot(aes(fill = sample), color = "black") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "sample", y = "normalized count") +
    ggtitle("vsdCount") +
    theme(legend.position = "none")

  pdf(file.path(outFolder, "normalizationCounts.pdf"), width=13,height=4)
  print(plot_grid(vsdCount_boxplot))
  dev.off()
  cat(paste0("[OUTPUT] VST 标准化方法的效果箱线图已保存至: ", file.path(outFolder, "normalizationCounts.pdf"), "\n"))
}

#绘制样本距离聚类热图；
cat("\n  - 正在计算样本间距离并绘制聚类热图...\n")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$sample

# 保存样本距离矩阵
write.csv(as.data.frame(sampleDistMatrix),
          file=file.path(outFolder, "sample_clusters.csv"))
cat(paste0("[OUTPUT] 样本间距离矩阵已保存至: ", file.path(outFolder, "sample_clusters.csv"), "\n"))


#####################################
# 7. 脚本结束
#####################################
cat("\n========================================================\n")
cat("[SUCCESS] 所有分析步骤已成功完成！\n")
cat("========================================================\n")

# --- 关闭日志重定向 ---
cat(paste("\n脚本结束运行时间:", Sys.time(), "\n"))
cat("--- 日志结束 ---\n")
# 依序关闭 sink 连接
sink(type = "message")
sink()