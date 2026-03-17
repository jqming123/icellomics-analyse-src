#!/usr/bin/env bash

# @File       :0_generate_sra2fastq_job.sh
# @Description:Generate separate SLURM scripts for each SRA file found in 0_data.
# @Usage      :bash 0_generate_sra2fastq_job.sh <PROJECT_NAME>
# @Example    :bash 0_generate_sra2fastq_job.sh PRJNA667472_CHO

if [ "$#" -ne 1 ]; then
    echo "错误: 参数数量不正确！"
    echo "用法: bash 0_generate_sra2fastq_job.sh <PROJECT_NAME>"
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="dont_need_ref"  # 对应 config 中的新增项

# --- 加载项目配置以获取路径 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi
source "${CONFIG_PATH}"

# --- 定义并创建目录 ---
SRA_DIR="${PROJECT_DIR}/0_data"
# 新的脚本存放位置
SUB_JOB_DIR="${PROJECT_DIR}/2_jobs/sra2fastq"
LOG_DIR="${PROJECT_DIR}/3_logs/sra2fastq"
RESULTS_DIR="${PROJECT_DIR}/1_result"
FASTQ_DIR="${RESULTS_DIR}/0_fastq"

if [ ! -d "${SRA_DIR}" ]; then
    echo "错误: SRA目录不存在: ${SRA_DIR}"
    exit 1
fi

mkdir -p "${SUB_JOB_DIR}" "${LOG_DIR}" "${FASTQ_DIR}" "${TMP_DIR}"

# --- 遍历 SRA 目录并生成脚本 ---
echo "正在为各样本生成独立的 SLURM 脚本..."

# 计数器
COUNT=0

for sample_dir in "${SRA_DIR}"/*; do
    # 确保是目录
    if [ -d "${sample_dir}" ]; then
        sample_id=$(basename "${sample_dir}")
        SRA_FILE="${sample_dir}/${sample_id}.sra"

        # 检查 SRA 文件是否存在
        if [ ! -f "${SRA_FILE}" ]; then
            echo "警告: 跳过 ${sample_id}，未在文件夹中找到 .sra 文件。"
            continue
        fi

        JOB_SCRIPT_PATH="${SUB_JOB_DIR}/sra2fastq_${sample_id}.sh"
        
        # 生成单个样本的 SLURM 脚本
        cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=s2f_${sample_id}
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/sra2fastq_${sample_id}_%j.log

echo "=========================================================="
echo "Job started on \$(date)"
echo "Sample ID: ${sample_id}"
echo "=========================================================="

set -e

# --- 准备环境 ---
export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="dont_need_ref"
source "${CONFIG_PATH}"

source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}"

# --- 路径设置 ---
READS_DIR="${FASTQ_DIR}/${sample_id}/reads"
mkdir -p "\${READS_DIR}"

# --- 执行转换 ---
if [ ! -f "\${READS_DIR}/${sample_id}.fastq.gz" ] && [ ! -f "\${READS_DIR}/${sample_id}_1.fastq.gz" ]; then
    echo "正在转换: ${sample_id}"
    fasterq-dump -e \${SLURM_CPUS_PER_TASK} --split-3 --outdir "\${READS_DIR}" "${SRA_FILE}"
    
    echo "正在压缩: ${sample_id}"
    find "\${READS_DIR}" -name "*.fastq" -print0 | xargs -0 -P \${SLURM_CPUS_PER_TASK} pigz -f
else
    echo "FASTQ 文件已存在，跳过转换。"
fi

conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

        chmod +x "${JOB_SCRIPT_PATH}"
        ((COUNT++))
    fi
done

echo "完成！共生成了 ${COUNT} 个脚本。"
echo "脚本目录: ${SUB_JOB_DIR}"
echo "日志目录: ${LOG_DIR}"
echo ""
echo "你可以使用以下命令批量提交这些作业:"
echo "for f in ${SUB_JOB_DIR}/*.sh; do sbatch \$f; done"