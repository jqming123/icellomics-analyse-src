#!/usr/bin/env bash

# @File       :2_generate_peak_calling_job.sh (Updated)
# @Description:Generate SLURM scripts for ATAC-seq peak calling using MACS2 for all samples in 2_tagalign.
# @Usage      :bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME>
# @Example    :bash 2_generate_peak_calling_job.sh PRJNA667472_CHO CH_Ensembl

if [ "$#" -ne 2 ]; then
    echo "错误: 参数数量不正确！"
    echo "用法: bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME>"
    exit 1
fi

# 1. 提取参数并设置为环境变量（供 epi_config.sh 使用）
export PROJECT_NAME="$1"
export REF_NAME="$2" 

# 2. 加载项目配置
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"

if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

source "${CONFIG_PATH}"

if [ -z "${GSIZE}" ]; then
    echo "错误: GSIZE 未设置，请检查 epi_config.sh 是否正确支持 '${REF_NAME}'" >&2
    exit 1
fi

echo "正在为项目 ${PROJECT_NAME} (参考基因组: ${REF_NAME}, GSIZE: ${GSIZE}) 生成 Peak Calling 任务..."
echo "当前指定的队列：$QUEUE_NAME"

# --- 定义输入输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs/peak_calling"
LOG_DIR="${PROJECT_DIR}/3_logs/peak_calling"
RESULTS_DIR="${PROJECT_DIR}/1_result"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${PEAKS_DIR}"

# --- 遍历目录提取 BIOSAMPLE_NAME 并生成SLURM作业脚本 ---
# 匹配目录下的所有 tagAlign.gz 文件
for INPUT_TAGALIGN in "${POOLED_DIR}"/*.tagAlign.gz; do
    # 检查是否存在匹配的文件
    if [ ! -f "${INPUT_TAGALIGN}" ]; then
        echo "在 ${POOLED_DIR} 中未找到 *.tagAlign.gz 文件。"
        exit 0
    fi

    # 从文件名中提取 BIOSAMPLE_NAME (截取第一个点之前的部分)
    FILENAME=$(basename "${INPUT_TAGALIGN}")
    BIOSAMPLE_NAME="${FILENAME%%.*}"

    JOB_SCRIPT_PATH="${JOB_DIR}/callpeaks_${BIOSAMPLE_NAME}.sh"
    echo "正在为样本 ${BIOSAMPLE_NAME} 生成SLURM作业脚本: ${JOB_SCRIPT_PATH}"

    cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${BIOSAMPLE_NAME}_callpeaks
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/callpeaks_${BIOSAMPLE_NAME}_%j.log

echo "=========================================================="
echo "Job started on \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "=========================================================="

set -e

# --- 准备环境和配置 ---
export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="${REF_NAME}" 
source "${CONFIG_PATH}" 

source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}" 

# --- 运行 Peak Calling ---
echo "--- 开始运行 Peak Calling ---"
if [ ! -f "${INPUT_TAGALIGN}" ]; then
    echo "错误: 输入文件 ${INPUT_TAGALIGN} 未找到！"
    exit 1
fi

macs3 callpeak \\
-t "${INPUT_TAGALIGN}" \\
-n "${BIOSAMPLE_NAME}" \\
--outdir "${PEAKS_DIR}" \\
-f BED -g ${GSIZE} -p 0.01 --shift -100 --extsize 200 \\
--nomodel -B --SPMR --keep-dup all --call-summits --buffer-size 1000

conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

    chmod +x "${JOB_SCRIPT_PATH}"
done

echo "所有任务脚本生成完毕！请进入 ${JOB_DIR} 目录批量提交作业。"
