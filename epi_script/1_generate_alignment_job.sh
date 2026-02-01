#!/usr/bin/env bash

# @File       :1_generate_alignment_job.sh
# @Description:Generate a SLURM script for ATAC-seq alignment and pooling.
# @Usage      :bash 1_generate_alignment_job.sh <PROJECT_NAME> <REF_NAME> <FINAL_NAME> <SampleID1> [SampleID2] ...
# @Example    :bash 1_generate_alignment_job.sh PRJNA667472_CHO CriGri-PICRH-1.0 ATAC_pooled CriGri-PICRH-1.0 SRR12774931 SRR12774932  # 多个样本
# @Example    :bash 1_generate_alignment_job.sh PRJNA667472_CHO CriGri-PICRH-1.0 SRR12774931 SRR12774931             # 单个样本

if [ "$#" -lt 4 ]; then
    echo "错误: 参数不足！"
    echo "用法: bash 1_generate_alignment_job.sh <PROJECT_NAME> <REF_NAME> <FINAL_NAME> <SampleID1> [SampleID2] ..."
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="$2"
FINAL_NAME="$3" 
shift 3
SAMPLES=("$@")

# --- 加载项目配置 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi
source "${CONFIG_PATH}"

# --- 定义输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs"
LOG_DIR="${PROJECT_DIR}/3_logs"
SRA_DIR="${PROJECT_DIR}/0_data"
RESULTS_DIR="${PROJECT_DIR}/1_result"
FASTQ_DIR="${RESULTS_DIR}/0_fastq"
ALIGN_DIR="${RESULTS_DIR}/1_alignment"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${FASTQ_DIR}" "${ALIGN_DIR}" "${POOLED_DIR}" "${TMP_DIR}"

# --- 生成SLURM作业脚本 ---
JOB_SCRIPT_PATH="${JOB_DIR}/align_pool_${FINAL_NAME}.sh"
echo "正在生成SLURM作业脚本: ${JOB_SCRIPT_PATH}"

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${FINAL_NAME}_align
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/align_pool_${FINAL_NAME}_%j.log

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

# --- 循环处理每个样本 ---
SAMPLES=(${SAMPLES[@]})
TAGALIGN_FILES=() 

for sample_id in "\${SAMPLES[@]}"; do
    echo "--- 开始处理样本: \${sample_id} ---"

    # 1. SRA to FASTQ
    echo "步骤 1: SRA to FASTQ for \${sample_id}"
    SRA_FILE_PATH="${SRA_DIR}/\${sample_id}/\${sample_id}.sra"
    READS_DIR="${FASTQ_DIR}/\${sample_id}/reads"
    mkdir -p "\${READS_DIR}"

    if [ -f "\${SRA_FILE_PATH}" ]; then
        if [ ! -f "\${READS_DIR}/\${sample_id}.fastq.gz" ] && [ ! -f "\${READS_DIR}/\${sample_id}_1.fastq.gz" ]; then
            fasterq-dump -e \${SLURM_CPUS_PER_TASK} --split-3 --outdir "\${READS_DIR}" "\${SRA_FILE_PATH}"
            find "\${READS_DIR}" -name "*.fastq" -print0 | xargs -0 -P \${SLURM_CPUS_PER_TASK} pigz -f
        else
            echo "FASTQ文件已存在，跳过fasterq-dump。"
        fi
    else
        echo "错误: SRA文件 \${SRA_FILE_PATH} 不存在！"
        exit 1
    fi
    
    # 2. Alignment using ATAC_align.sh
    echo "步骤 2: Alignment for \${sample_id}"
    THREAD_MEM=\$(( ${MEM_MEDIUM%G} / ${THREADS} ))
    
    bash "\${EPI_SCRIPT_DIR}/ATAC_align.sh" \\
        "${FASTQ_DIR}" \\
        "${ALIGN_DIR}" \\
        "\${sample_id}" \\
        "\${SLURM_CPUS_PER_TASK}" \\
        "\${THREAD_MEM}"

    TAGALIGN_FILES+=("${ALIGN_DIR}/\${sample_id}/bowtie2/\${sample_id}.tn5.tagAlign.gz")
done

# --- 步骤 3: 根据样本数量决定合并或移动 TagAlign 文件 ---
NUM_SAMPLES=\${#SAMPLES[@]}

if [ "\$NUM_SAMPLES" -gt 1 ]; then
    echo "--- 检测到 \${NUM_SAMPLES} 个样本，开始合并 TagAlign 文件 ---"
    python "\${EPI_SCRIPT_DIR}/poolTagAligns.py" \\
        "${FINAL_NAME}" \\
        "\${TAGALIGN_FILES[@]}"
    
    # 将合并后的文件移动到指定的pooled目录
    mv "./${FINAL_NAME}.pooled.tn5.tagAlign.gz" "${POOLED_DIR}/${FINAL_NAME}.tn5.tagAlign.gz"
    echo "合并完成，文件: ${POOLED_DIR}/${FINAL_NAME}.pooled.tn5.tagAlign.gz"

elif [ "\$NUM_SAMPLES" -eq 1 ]; then
    echo "--- 检测到单个样本，跳过合并，直接重命名并移动文件 ---"
    # TAGALIGN_FILES数组中只有一个元素
    SINGLE_TAGALIGN_FILE="\${TAGALIGN_FILES[0]}"
    FINAL_OUTPUT_PATH="${POOLED_DIR}/${FINAL_NAME}.tn5.tagAlign.gz"
    
    mv "\${SINGLE_TAGALIGN_FILE}" "\${FINAL_OUTPUT_PATH}"
    echo "处理完成，文件: \${FINAL_OUTPUT_PATH}"
else
    echo "警告: 没有有效的 TagAlign 文件生成，未执行任何操作。"
fi


conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

chmod +x "${JOB_SCRIPT_PATH}"
echo "成功！请检查生成的脚本: ${JOB_SCRIPT_PATH}"
echo "使用以下命令提交作业: sbatch ${JOB_SCRIPT_PATH}"