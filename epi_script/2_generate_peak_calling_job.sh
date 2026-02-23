#!/usr/bin/env bash

# @File       :2_generate_peak_calling_job.sh (Updated)
# @Description:Generate a SLURM script for ATAC-seq peak calling using MACS2.
# @Usage      :bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME> <FINAL_NAME>
# @Example    :bash 2_generate_peak_calling_job.sh PRJNA667472_CHO CH_Ensemble SRR12774932

if [ "$#" -ne 3 ]; then
    echo "错误: 参数数量不正确！"
    echo "用法: bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME> <FINAL_NAME>"
    exit 1
fi

# 1. 提取参数并设置为环境变量（供 epi_config.sh 使用）
export PROJECT_NAME="$1"
export REF_NAME="$2" 
FINAL_NAME="$3"

# 2. 加载项目配置
# 注意：必须在设置完 PROJECT_NAME 和 REF_NAME 后加载
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"

if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

# 加载配置，这会自动根据 REF_NAME 设置 GSIZE, REF_GENOME, BOWTIE2_INDEX 等
source "${CONFIG_PATH}"

# 3. 验证关键变量是否已成功加载 (可选，建议保留以增加脚本健壮性)
if [ -z "${GSIZE}" ]; then
    echo "错误: GSIZE 未设置，请检查 epi_config.sh 是否正确支持 '${REF_NAME}'" >&2
    exit 1
fi

# --- 接下来你可以直接使用变量，例如 $GSIZE, $PROJECT_DIR, $EPI_CONDA_ENV_NAME 等 ---
echo "正在为项目 ${PROJECT_NAME} (参考基因组: ${REF_NAME}, GSIZE: ${GSIZE}) 生成 Peak Calling 任务..."

# --- 定义输入输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs"
LOG_DIR="${PROJECT_DIR}/3_logs"
RESULTS_DIR="${PROJECT_DIR}/1_result"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"
INPUT_TAGALIGN="${POOLED_DIR}/${FINAL_NAME}.fixed.tagAlign.gz"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${PEAKS_DIR}"

# --- 生成SLURM作业脚本 ---
JOB_SCRIPT_PATH="${JOB_DIR}/callpeaks_${FINAL_NAME}.sh"
echo "正在生成SLURM作业脚本: ${JOB_SCRIPT_PATH}"

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${FINAL_NAME}_callpeaks
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/callpeaks_${FINAL_NAME}_%j.log

echo "=========================================================="
echo "Job started on \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "=========================================================="

set -e

# --- 准备环境和配置 ---
export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="${REF_NAME}" 
source "${CONFIG_PATH}" # 在SLURM脚本内部再次加载配置，以确保所有变量可用

source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}" # 激活ATAC_E4环境

# --- 运行 Peak Calling ---
echo "--- 开始运行 Peak Calling ---"
if [ ! -f "${INPUT_TAGALIGN}" ]; then
    echo "错误: 输入文件 ${INPUT_TAGALIGN} 未找到！"
    exit 1
fi

# 直接整合 callPeaks_ATAC_E4.sh 的 MACS2 命令
# 原始变量映射:
# in_file -> "${INPUT_TAGALIGN}"
# out_dir -> "${PEAKS_DIR}"
# sample_id -> "${FINAL_NAME}"

macs2 callpeak \\
-t "${INPUT_TAGALIGN}" \\
-n "${FINAL_NAME}" \\
--outdir "${PEAKS_DIR}" \\
-f BED -g ${GSIZE} -p 0.01 --shift -100 --extsize 200 \\
--nomodel -B --SPMR --keep-dup all --call-summits --buffer-size 1000

conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

chmod +x "${JOB_SCRIPT_PATH}"
echo "成功！请检查生成的脚本: ${JOB_SCRIPT_PATH}"
echo "使用以下命令提交作业: sbatch ${JOB_SCRIPT_PATH}"