#!/usr/bin/env bash

# @File       :generate_peak_calling_job_PESE.sh
# @Description:Generate a SLURM script for ATAC-seq peak calling using MACS3.
# @Usage      :bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME> <BIOSAMPLE_ID> <SEQ_TYPE>
# @Example    :bash 2_generate_peak_calling_job.sh PRJNA667472_CHO CH_Ensemble SRR12774932 PE

if [ "$#" -ne 4 ]; then
    echo "错误: 参数数量不正确！"
    echo "用法: bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME> <BIOSAMPLE_ID> <PE|SE>"
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="$2" 
BIOSAMPLE_ID="$3"
SEQ_TYPE=$(echo "$4" | tr '[:lower:]' '[:upper:]')

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

MAP_FILE="${PROJECT_DIR}/0_data/sample_run_map.tsv"
if [ ! -f "${MAP_FILE}" ]; then
    echo "错误: 未找到 sample_run_map.tsv: ${MAP_FILE}"
    exit 1
fi
dos2unix "${MAP_FILE}" 2>/dev/null

# 提取当前 BioSample ID 对应的所有 Run ID
RUN_IDS=($(awk -v bio="$BIOSAMPLE_ID" -F'\t' '$1 == bio {print $2}' "${MAP_FILE}"))

if [ ${#RUN_IDS[@]} -eq 0 ]; then
    echo "错误: 在 sample_run_map.tsv 中未找到 BioSample ID: ${BIOSAMPLE_ID}"
    exit 1
fi

echo "正在为项目 ${PROJECT_NAME} (参考基因组: ${REF_NAME}, GSIZE: ${GSIZE}, 测序类型: ${SEQ_TYPE}) 生成 Peak Calling 任务..."
echo "BioSample ID: ${BIOSAMPLE_ID} 包含的 Run IDs: ${RUN_IDS[*]}"

# --- 定义输入输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs"
LOG_DIR="${PROJECT_DIR}/3_logs"
RESULTS_DIR="${PROJECT_DIR}/1_result"
ALIGN_DIR="${RESULTS_DIR}/1_alignment"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"
MERGED_BAM_DIR="${RESULTS_DIR}/${BIOSAMPLE_ID}/merged_bam"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${PEAKS_DIR}"

# --- 生成SLURM作业脚本 ---
JOB_SCRIPT_PATH="${JOB_DIR}/callpeaks_${BIOSAMPLE_ID}.sh"
echo "正在生成SLURM作业脚本: ${JOB_SCRIPT_PATH}"

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${BIOSAMPLE_ID}_callpeaks
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/callpeaks_${BIOSAMPLE_ID}_%j.log

echo "=========================================================="
echo "Job started on \$(date)"
echo "Job ID: \${SLURM_JOB_ID}"
echo "=========================================================="

set -e

export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="${REF_NAME}" 
source "${CONFIG_PATH}"

source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}"

SEQ_TYPE="${SEQ_TYPE}"
BIOSAMPLE_ID="${BIOSAMPLE_ID}"

echo "--- 开始处理数据准备 ---"

if [ "\${SEQ_TYPE}" == "PE" ]; then
    RUN_ARRAY=(${RUN_IDS[*]})
    NUM_RUNS=\${#RUN_ARRAY[@]}
    
    if [ "\${NUM_RUNS}" -gt 1 ]; then
        echo "检测到多个 Run，开始合并 BAM 文件..."
        mkdir -p "${MERGED_BAM_DIR}"
        INPUT_FILE="${MERGED_BAM_DIR}/\${BIOSAMPLE_ID}.merged.nodup.bam"
        
        BAM_FILES=()
        for run_id in "\${RUN_ARRAY[@]}"; do
            BAM_FILES+=("${ALIGN_DIR}/\${run_id}/bowtie2/\${run_id}.nodup.bam")
        done
        
        samtools merge -@ \${SLURM_CPUS_PER_TASK} "\${INPUT_FILE}" "\${BAM_FILES[@]}"
        samtools index -@ \${SLURM_CPUS_PER_TASK} "\${INPUT_FILE}"
        echo "合并完成: \${INPUT_FILE}"
    else
        echo "单 Run，无需合并 BAM 文件。"
        INPUT_FILE="${ALIGN_DIR}/\${RUN_ARRAY[0]}/bowtie2/\${RUN_ARRAY[0]}.nodup.bam"
    fi
    
    MACS3_FORMAT="BAMPE"
    MACS3_EXTRA_PARAMS=""
    
elif [ "\${SEQ_TYPE}" == "SE" ]; then
    # 单端使用经过 Tn5 偏移校正并可能已在 pool 步骤合并的 tagAlign 文件
    INPUT_FILE="${POOLED_DIR}/${BIOSAMPLE_ID}.tn5.tagAlign.gz"
    MACS3_FORMAT="BED"
    MACS3_EXTRA_PARAMS="--nomodel --shift -100 --extsize 200"
else
    echo "错误的测序类型!"
    exit 1
fi

if [ ! -f "\${INPUT_FILE}" ]; then
    echo "错误: 输入文件 \${INPUT_FILE} 未找到！"
    exit 1
fi

echo "--- 开始运行 MACS3 Peak Calling ---"

macs3 callpeak \\
-t "\${INPUT_FILE}" \\
-n "\${BIOSAMPLE_ID}" \\
--outdir "${PEAKS_DIR}" \\
-f \${MACS3_FORMAT} -g ${GSIZE} -q 0.01 \${MACS3_EXTRA_PARAMS} \\
-B --SPMR --keep-dup all --call-summits --buffer-size 1000

conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

chmod +x "${JOB_SCRIPT_PATH}"
echo "成功！请检查生成的脚本: ${JOB_SCRIPT_PATH}"
echo "使用以下命令提交作业: sbatch ${JOB_SCRIPT_PATH}"
