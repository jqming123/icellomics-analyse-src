#!/bin/bash

# ==============================================================================
# 脚本功能: 智能地将 scaffolds 拆分成 N 个分区文件 (intervals)。
#           拆分原则是确保每个分区文件包含的基因组碱基总数大致相等，
#           从而实现下游并行任务的负载均衡。
# 运行方式: 在 01_scripts/ 目录下执行 `bash 01a_prepare_intervals_by_size.sh`
# ==============================================================================

# 引入项目配置文件
source ./config.sh

# --- 可调参数 ---
# 您希望将整个基因组拆分成多少个并行的任务包
SCATTER_COUNT=50

# --- 路径定义 ---
# 定义存放最终分区文件的目录
# 这个路径需要与 02_generate_genomicsdb_jobs.sh 脚本中的路径保持一致
INTERVAL_DIR="${RESULTS_DIR}/04_genomicsdb/intervals_by_size"

# --- 主逻辑 ---
echo "=== Starting Balanced Interval Preparation ==="

# 1. 创建输出目录
mkdir -p "$INTERVAL_DIR"
echo "Interval files will be created in: ${INTERVAL_DIR}"

# 2. 确保参考基因组索引 (.fai) 文件存在
FAI_FILE="${REF_GENOME}.fai"
if [ ! -f "${FAI_FILE}" ]; then
    echo "Reference index .fai not found. Creating it with 'samtools faidx'..."
    samtools faidx "${REF_GENOME}"
    echo "Index created at: ${FAI_FILE}"
fi

# 3. 计算基因组总大小和每个分区的目标大小
echo "Calculating target size for each scatter job..."
TOTAL_SIZE=$(awk '{s+=$2} END {print s}' "${FAI_FILE}")
TARGET_CHUNK_SIZE=$((TOTAL_SIZE / SCATTER_COUNT))

echo "Total genome size: ${TOTAL_SIZE} bp"
echo "Target size per job: ~${TARGET_CHUNK_SIZE} bp for ${SCATTER_COUNT} jobs."

# 4. 使用 awk 脚本进行智能拆分
echo "Splitting scaffolds into balanced groups..."
awk -v target_size="${TARGET_CHUNK_SIZE}" -v prefix="${INTERVAL_DIR}/scatter_" '
BEGIN {
    batch_num = 0;
    current_size = 0;
    # 使用三位数字格式化，例如 scatter_000, scatter_001
    outfile = sprintf("%s%03d", prefix, batch_num);
}
{
    if (current_size > 0 && (current_size + $2) > target_size) {
        close(outfile);
        batch_num++;
        outfile = sprintf("%s%03d", prefix, batch_num);
        current_size = 0;
    }
    print $1 > outfile;
    current_size += $2;
}
' "${FAI_FILE}"

# 5. 报告结果
num_files=$(ls -1 "${INTERVAL_DIR}/scatter_"* | wc -l)
echo "Successfully created ${num_files} interval files."
echo "=== Interval Preparation Complete ==="