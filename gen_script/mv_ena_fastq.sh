#!/bin/bash

# --- 配置根目录 ---
PROJECT_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects"

# --- 检查输入参数 ---
if [ -z "$1" ]; then
    echo "错误: 未提供项目目录名。"
    echo "用法: bash mv_ena_fastq.sh <PROJECT_NAME>"
    echo "示例: bash mv_ena_fastq.sh PRJEB44115_H9"
    exit 1
fi

PROJECT_NAME=$1
PROJECT_DIR="${PROJECT_ROOT}/${PROJECT_NAME}"
SRA_DIR="${PROJECT_DIR}/00_data/dna_sra"
RAW_DIR="${PROJECT_DIR}/00_data/raw_fastq"

# 检查项目目录是否存在
if [ ! -d "$PROJECT_DIR" ]; then
    echo "错误: 找不到项目目录 $PROJECT_DIR"
    exit 1
fi

# 确保目标目录存在
mkdir -p "$RAW_DIR"

echo ">>> 开始处理项目: $PROJECT_NAME"

# 移动并重命名 fastq.gz 文件
# 查找所有 _1.fastq.gz 文件并循环处理
shopt -s nullglob # 防止找不到文件时循环报错
for f1 in "${SRA_DIR}"/*_1.fastq.gz; do
    # 提取样本基本名 (如 ERR5654582)
    sample=$(basename "$f1" "_1.fastq.gz")
    f2="${SRA_DIR}/${sample}_2.fastq.gz"

    echo "正在处理样本: $sample"

    # 处理 R1
    mv "$f1" "${RAW_DIR}/${sample}_r1.fq.gz"

    # 处理 R2 (如果存在)
    if [ -f "$f2" ]; then
        mv "$f2" "${RAW_DIR}/${sample}_r2.fq.gz"
    else
        echo "警告: 找不到 $sample 的 R2 文件"
    fi
done

echo ">>> 处理完成！"
echo "文件已移动至: $RAW_DIR"
ls -lh "$RAW_DIR"
