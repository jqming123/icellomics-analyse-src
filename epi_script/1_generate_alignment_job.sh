#!/usr/bin/env bash

# @File       :1_generate_alignment_job.sh
# @Description:Generate SLURM scripts automatically from sample_run_map.tsv
# @Usage      :bash 1_generate_alignment_job.sh <PROJECT_NAME> <REF_NAME> 

if [ "$#" -lt 2 ]; then
    echo "错误: 参数不足！"
    echo "用法: bash 1_generate_alignment_job.sh <PROJECT_NAME> <REF_NAME>"
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="$2"

# --- 加载项目配置 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi
source "${CONFIG_PATH}"

echo "当前项目：$PROJECT_NAME"
echo "当前指定的参考基因组：$REF_NAME"
echo "当前指定的队列：$QUEUE_NAME"

# --- sample_run_map.tsv 路径 ---
MAP_FILE="${PROJECT_DIR}/0_data/sample_run_map.tsv"

if [ ! -f "${MAP_FILE}" ]; then
    echo "error: 未找到 sample_run_map.tsv: ${MAP_FILE}"
    exit 1
fi

dos2unix $MAP_FILE
TOTAL_JOBS=$(grep -v '^#' "${MAP_FILE}" | awk -F'\t' '{print $1}' | sort -u | wc -l)
echo "检测到 $TOTAL_JOBS 个唯一的 Biosample，准备生成脚本..."

# --- 定义输出目录 ---
JOB_DIR="${PROJECT_DIR}/2_jobs/align_pool"
LOG_DIR="${PROJECT_DIR}/3_logs/align_pool"
RESULTS_DIR="${PROJECT_DIR}/1_result"
FASTQ_DIR="${RESULTS_DIR}/0_fastq"
ALIGN_DIR="${RESULTS_DIR}/1_alignment"
RESULT_DIR="${RESULTS_DIR}/2_tagalign"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${ALIGN_DIR}" "${RESULT_DIR}" "${TMP_DIR}"

echo "正在解析 sample_run_map.tsv ..."
# 表格无表头。按第一列分组（第一列是BIOSAMPLE_ID，第二列是RunID）
# 跳过 # 注释行
grep -v '^#' "${MAP_FILE}" | awk -F'\t' '{print $1"\t"$2}' | \
sort | \
awk -F'\t' '
{
    biosample=$1
    run=$2
    if (biosample in runs) {
        runs[biosample]=runs[biosample]" "run
    } else {
        runs[biosample]=run
    }
}
END {
    for (b in runs) {
        print b"\t"runs[b]
    }
}' | while IFS=$'\t' read -r BIOSAMPLE_ID RUN_ID_STR
do

    echo "--------------------------------------------------"
    echo "处理 Biosample: ${BIOSAMPLE_ID}"
    echo "RunIDs: ${RUN_ID_STR}"

    # 转成数组
    RUN_IDS=(${RUN_ID_STR})

    JOB_SCRIPT_PATH="${JOB_DIR}/align_pool_${BIOSAMPLE_ID}.sh"

    cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${BIOSAMPLE_ID}_align
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/align_pool_${BIOSAMPLE_ID}_%j.log

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

RUN_IDS=(${RUN_ID_STR})
TAGALIGN_FILES=()

for run_id in "\${RUN_IDS[@]}"; do
    echo "--- 开始处理 Run: \${run_id} ---"
    
    THREAD_MEM=\$(( ${MEM_MEDIUM%G} / ${THREADS} ))
    
    bash "\${EPI_SCRIPT_DIR}/ATAC_align.sh" \\
        "${FASTQ_DIR}" \\
        "${ALIGN_DIR}" \\
        "\${run_id}" \\
        "${BIOSAMPLE_ID}" \\
        "\${SLURM_CPUS_PER_TASK}" \\
        "\${THREAD_MEM}"
        
    TAGALIGN_FILES+=("${ALIGN_DIR}/\${run_id}/bowtie2/\${run_id}.tn5.tagAlign.gz")
done

NUM_RUNS=\${#RUN_IDS[@]}

if [ "\$NUM_RUNS" -gt 1 ]; then
    echo "--- 合并 \${NUM_RUNS} 个 Run 数据到 Biosample: ${BIOSAMPLE_ID} ---"
    python "\${EPI_SCRIPT_DIR}/poolTagAligns.py" \\
        "${BIOSAMPLE_ID}" \\
        "\${TAGALIGN_FILES[@]}"
    mv "./${BIOSAMPLE_ID}.pooled.tn5.tagAlign.gz" "${RESULT_DIR}/${BIOSAMPLE_ID}.tn5.tagAlign.gz"
elif [ "\$NUM_RUNS" -eq 1 ]; then
    echo "--- 单 Run 处理 ---"
    mv "\${TAGALIGN_FILES[0]}" "${RESULT_DIR}/${BIOSAMPLE_ID}.tn5.tagAlign.gz"
fi

conda deactivate

echo "=========================================================="
echo "Job finished on \$(date)"
echo "=========================================================="
EOF
    echo "已生成: align_pool_${BIOSAMPLE_ID}.sh"
    ((COUNT++))
    chmod +x "${JOB_SCRIPT_PATH}"


done

echo "所有作业脚本生成完成！总计生成脚本数量：$TOTAL_JOBS"
echo "脚本所在路径：$JOB_DIR"
echo ""
echo "你可以使用以下命令批量提交这些作业:"
echo "for f in ${JOB_DIR}/*.sh; do sbatch \$f; done"
