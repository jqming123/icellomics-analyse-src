#!/usr/bin/env bash

# @File       :generate_featurecounts_job.sh
# @Description:Generate SLURM scripts for featureCounts (FRiP calculation) for all samples.
# @Usage      :bash generate_featurecounts_job.sh <PROJECT_NAME>

if [ "$#" -ne 1 ]; then
    echo "错误: 参数数量不正确！"
    echo "用法: bash generate_featurecounts_job.sh <PROJECT_NAME>"
    exit 1
fi

export PROJECT_NAME="$1"
# 强制使用不需要参考基因组的配置
export REF_NAME="dont_need_ref" 

# --- 加载项目配置 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"

if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

source "${CONFIG_PATH}"

echo "正在为项目 ${PROJECT_NAME} 生成 featureCounts 任务..."
echo "当前指定的队列：$QUEUE_NAME"

# --- 定义输入输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs/featureCounts"
LOG_DIR="${PROJECT_DIR}/3_logs/featureCounts"
RESULTS_DIR="${PROJECT_DIR}/1_result"
ALIGN_DIR="${RESULTS_DIR}/1_alignment"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"
FC_DIR="${RESULTS_DIR}/4_featureCounts"
MAP_FILE="${PROJECT_DIR}/0_data/sample_run_map.tsv"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${FC_DIR}"

if [ ! -f "${MAP_FILE}" ]; then
    echo "错误: 映射文件未找到: ${MAP_FILE}"
    exit 1
fi

# 声明关联数组用于存储 biosample 和对应的 BAM 列表
declare -A sample_bams

# 读取映射文件并组装 BAM 路径
while IFS=$'\t' read -r biosample run; do
    if [[ -n "$biosample" && -n "$run" ]]; then
        bam_path="${ALIGN_DIR}/${run}/bowtie2/${run}.nodup.bam"
        if [[ -z "${sample_bams[$biosample]}" ]]; then
            sample_bams[$biosample]="${bam_path}"
        else
            sample_bams[$biosample]="${sample_bams[$biosample]} ${bam_path}"
        fi
    fi
done < "${MAP_FILE}"

# 遍历每个 biosample 生成对应的 SLURM 脚本
for BIOSAMPLE_NAME in "${!sample_bams[@]}"; do
    BAM_FILES="${sample_bams[$BIOSAMPLE_NAME]}"
    PEAK_FILE="${PEAKS_DIR}/${BIOSAMPLE_NAME}_peaks.narrowPeak"
    SAF_FILE="${FC_DIR}/${BIOSAMPLE_NAME}_peaks.saf"
    OUT_COUNTS="${FC_DIR}/${BIOSAMPLE_NAME}.peaks.counts"
    
    JOB_SCRIPT_PATH="${JOB_DIR}/fc_${BIOSAMPLE_NAME}.sh"
    echo "正在为样本 ${BIOSAMPLE_NAME} 生成SLURM作业脚本: ${JOB_SCRIPT_PATH}"

    cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${BIOSAMPLE_NAME}_fc
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_SMALL}
#SBATCH --time=1-00:00:00
#SBATCH --output=${LOG_DIR}/fc_${BIOSAMPLE_NAME}_%j.log

echo "=========================================================="
echo "Job started on \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "=========================================================="

set -e

# --- 准备环境 ---
source "${CONFIG_PATH}" 
source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}" 

# --- 文件检查与准备 ---
if [ ! -f "${PEAK_FILE}" ]; then
    echo "错误: 未找到 peak 文件 ${PEAK_FILE}"
    exit 1
fi

echo "Generating SAF file..."
awk -F \$'\t' 'BEGIN {OFS=FS} {print \$4, \$1, \$2+1, \$3, "."}' "${PEAK_FILE}" > "${SAF_FILE}"

# --- 运行 featureCounts ---
echo "--- 开始运行 featureCounts ---"
FC_LOG="${FC_DIR}/${BIOSAMPLE_NAME}_fc.log"

# 运行并将标准错误(包含统计信息)重定向到日志文件
featureCounts -p -F SAF -a "${SAF_FILE}" --fracOverlap 0.2 -o "${OUT_COUNTS}" ${BAM_FILES} 2> "\${FC_LOG}"

# 将 featureCounts 的日志输出到 slurm 的 log 中方便查看
cat "\${FC_LOG}"

# --- 自动提取并计算 FRiP ---
# 提取形如 "Successfully assigned alignments : 3043655 (11.0%)" 中的百分比
FRIP_PERCENT=\$(grep "Successfully assigned alignments" "\${FC_LOG}" | grep -o '([0-9.]*%)' | tr -d '()')

echo "=========================================================="
echo "FRiP SUMMARY FOR BIOSAMPLE: ${BIOSAMPLE_NAME}"
if [ -n "\${FRIP_PERCENT}" ]; then
    echo "Fraction of Reads in Peaks (FRiP): \${FRIP_PERCENT}"
else
    echo "Fraction of Reads in Peaks (FRiP): Failed to parse from log."
fi
echo "=========================================================="

conda deactivate

echo "Job finished on \$(date)"
EOF

    chmod +x "${JOB_SCRIPT_PATH}"
done

echo "所有任务脚本生成完毕！请进入 ${JOB_DIR} 目录批量提交作业。"
