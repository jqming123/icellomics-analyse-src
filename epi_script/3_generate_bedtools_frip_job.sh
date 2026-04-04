#!/usr/bin/env bash

# @File       :3_generate_bedtools_frip_job.sh
# @Description:Generate SLURM scripts to calculate FRiP (Fraction of Reads in Peaks) for all samples.
# @Usage      :bash 3_generate_bedtools_frip_job.sh <PROJECT_NAME>
# @Example    :bash 3_generate_bedtools_frip_job.sh PRJDB10440_HEK293

if [ "$#" -ne 1 ]; then
    echo "错误: 参数数量不正确！"
    echo "用法: bash 3_generate_frip_job.sh <PROJECT_NAME>"
    exit 1
fi

# 1. 设置环境变量供 epi_config.sh 使用
export PROJECT_NAME="$1"
export REF_NAME="dont_need_ref" # 不需要参考基因组

# 2. 加载项目配置
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"

if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

source "${CONFIG_PATH}"

echo "正在为项目 ${PROJECT_NAME} 生成 FRiP 计算任务..."
echo "当前指定的队列：$QUEUE_NAME"

# --- 定义输入输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs/frip"
LOG_DIR="${PROJECT_DIR}/3_logs/frip"
RESULTS_DIR="${PROJECT_DIR}/1_result"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"
QC_DIR="${RESULTS_DIR}/4_qc/FRiP" # 用于存放 FRiP 计算结果

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${QC_DIR}"

# --- 遍历目录提取 BIOSAMPLE_NAME 并生成SLURM作业脚本 ---
# 匹配目录下的所有 tagAlign.gz 文件 (你的文件可能是 .tn5.tagAlign.gz)
for INPUT_TAGALIGN in "${POOLED_DIR}"/*.tagAlign.gz; do
    if [ ! -f "${INPUT_TAGALIGN}" ]; then
        echo "在 ${POOLED_DIR} 中未找到 *.tagAlign.gz 文件。"
        exit 0
    fi

    # 提取 BIOSAMPLE_NAME (截取第一个点之前的部分)
    FILENAME=$(basename "${INPUT_TAGALIGN}")
    BIOSAMPLE_NAME="${FILENAME%%.*}"
    
    # 对应的 Peak 文件路径
    PEAK_FILE="${PEAKS_DIR}/${BIOSAMPLE_NAME}_peaks.narrowPeak"
    OUTPUT_FRIP_FILE="${QC_DIR}/${BIOSAMPLE_NAME}_FRiP.tsv"

    JOB_SCRIPT_PATH="${JOB_DIR}/frip_${BIOSAMPLE_NAME}.sh"
    echo "正在为样本 ${BIOSAMPLE_NAME} 生成SLURM作业脚本: ${JOB_SCRIPT_PATH}"

    cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${BIOSAMPLE_NAME}_frip
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_SMALL}
#SBATCH --time=1-00:00:00
#SBATCH --output=${LOG_DIR}/frip_${BIOSAMPLE_NAME}_%j.log

echo "=========================================================="
echo "Job started on \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "=========================================================="

set -e

# --- 准备环境和配置 ---
source "${CONDA_PROFILE_PATH}"
conda activate "${EPI_CONDA_ENV_NAME}" 

# --- 检查输入文件 ---
if [ ! -f "${INPUT_TAGALIGN}" ]; then
    echo "错误: TagAlign 文件 ${INPUT_TAGALIGN} 未找到！"
    exit 1
fi

if [ ! -f "${PEAK_FILE}" ]; then
    echo "错误: Peak 文件 ${PEAK_FILE} 未找到！请确认是否已经完成 Peak Calling。"
    exit 1
fi

echo "--- 开始计算 FRiP ---"
echo "Biosample: ${BIOSAMPLE_NAME}"
echo "Reads: ${INPUT_TAGALIGN}"
echo "Peaks: ${PEAK_FILE}"

# 1. 统计落在 peak 内的 reads 数 (结合 zcat 和 bedtools 进程替换)
reads_in_peaks=\$(bedtools intersect -a <(zcat -f "${INPUT_TAGALIGN}") -b "${PEAK_FILE}" -wa -u | wc -l)

# 2. 统计 tagAlign 的总 reads 数
total_reads=\$(zcat -f "${INPUT_TAGALIGN}" | wc -l)

# 3. 计算 FRiP 并格式化输出
if [[ \$total_reads -gt 0 ]]; then
    frip=\$(bc <<< "scale=2; 100 * \$reads_in_peaks / \$total_reads")
    
    # 将结果打印到日志并保存到统一的结果文件中
    echo -e "Biosample\tTotal_Reads\tReads_in_Peaks\tFRiP" > "${OUTPUT_FRIP_FILE}"
    echo -e "${BIOSAMPLE_NAME}\t\${total_reads}\t\${reads_in_peaks}\t\${frip}%" >> "${OUTPUT_FRIP_FILE}"
    
    echo "结果: FRiP = \${frip}%"
else
    echo "警告: 总 Reads 数为 0，无法计算 FRiP。"
    echo -e "Biosample\tTotal_Reads\tReads_in_Peaks\tFRiP" > "${OUTPUT_FRIP_FILE}"
    echo -e "${BIOSAMPLE_NAME}\t0\t0\t0.00%" >> "${OUTPUT_FRIP_FILE}"
fi

conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

    chmod +x "${JOB_SCRIPT_PATH}"
done

echo "所有任务脚本生成完毕！请进入 ${JOB_DIR} 目录批量提交作业。"
