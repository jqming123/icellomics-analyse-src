#!/bin/bash
# 功能: 自动生成 GenomicsDBImport 所需的 sample_map.txt 文件
# 使用方式: 在 01_scripts/ 目录下运行 `bash generate_sample_map.sh <PROJECT_NAME>`
set -eo pipefail

# --- 获取参数 ---
if [ -z "$1" ] ; then
    echo "用法: $0 <PROJECT_NAME>" >&2
    echo "例子: $0 MyProject" >&2
    exit 1
fi

PROJECT_NAME="$1"
export PROJECT_NAME
echo "当前项目名称 (PROJECT_NAME): ${PROJECT_NAME}"

export REF_NAME="dont_need_ref"
# 引入配置
source ./config.sh

# gVCF 目录
GVCF_DIR="${RESULTS_DIR}/03_gvcf"
OUTPUT_FILE="${GVCF_DIR}/sample_map.txt"

echo "Generating sample_map.txt in ${GVCF_DIR} ..."

# 清空旧文件
> "${OUTPUT_FILE}"

# 遍历 gVCF 文件 (*.g.vcf.gz)，根据文件名提取样本名
for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
    # 样本名：去掉路径和扩展名
    sample=$(basename "$gvcf" .g.vcf.gz)
    echo -e "${sample}\t${gvcf}" >> "${OUTPUT_FILE}"
done

echo "Done. sample_map.txt generated at: ${OUTPUT_FILE}"

