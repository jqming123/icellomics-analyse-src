#!/usr/bin/env bash

# Generate one SLURM job per biosample. Technical runs are merged within each
# biological replicate before duplicate marking/removal.
# Preferred sample_run_map.tsv columns:
# biosample_id <TAB> biological_replicate <TAB> run_id [<TAB> PE|SE|AUTO]

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "用法: bash 1_generate_alignment_job.sh <PROJECT_NAME> <REF_NAME>" >&2
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="$2"
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}" >&2
    exit 1
fi
source "${CONFIG_PATH}"

MAP_FILE="${PROJECT_DIR}/0_data/sample_run_map.tsv"
MAP_PARSER="${EPI_SCRIPT_DIR}/atac_qc/sample_map.py"
if [ ! -f "${MAP_FILE}" ]; then
    echo "error: 未找到 sample_run_map.tsv: ${MAP_FILE}" >&2
    exit 1
fi

dos2unix "${MAP_FILE}"
python "${MAP_PARSER}" --map "${MAP_FILE}" --project-dir "${PROJECT_DIR}" --validate
mapfile -t BIOSAMPLES < <(python "${MAP_PARSER}" --map "${MAP_FILE}" --list-biosamples)
TOTAL_JOBS=${#BIOSAMPLES[@]}

JOB_DIR="${PROJECT_DIR}/2_jobs/align_pool"
LOG_DIR="${PROJECT_DIR}/3_logs/align_pool"
RESULTS_DIR="${PROJECT_DIR}/1_result"
FASTQ_DIR="${RESULTS_DIR}/0_fastq"
ALIGN_DIR="${RESULTS_DIR}/1_alignment"
TAGALIGN_DIR="${RESULTS_DIR}/2_tagalign"
mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${ALIGN_DIR}" "${TAGALIGN_DIR}" "${TMP_DIR}"

for BIOSAMPLE_ID in "${BIOSAMPLES[@]}"; do
    JOB_SCRIPT_PATH="${JOB_DIR}/align_pool_${BIOSAMPLE_ID}.sh"
    cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${BIOSAMPLE_ID}_align
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS}
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/align_pool_${BIOSAMPLE_ID}_%j.log

# nounset (-u) is deliberately not enabled here: the cluster conda
# activate/deactivate hooks reference optional variables (ZSH_VERSION,
# JAVA_HOME) without a default and would abort the job during activation.
set -eo pipefail
export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="${REF_NAME}"
source "${CONFIG_PATH}"
source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}"

mapfile -t REPLICATE_ROWS < <(python "\${EPI_SCRIPT_DIR}/atac_qc/sample_map.py" \
    --map "${MAP_FILE}" --project-dir "${PROJECT_DIR}" \
    --biosample "${BIOSAMPLE_ID}" --emit-replicates)

REP_TAGALIGNS=()
REP_FRAGMENTS=()
for row in "\${REPLICATE_ROWS[@]}"; do
    IFS=$'\t' read -r replicate layout legacy run_csv <<< "\${row}"
    IFS=',' read -r -a run_ids <<< "\${run_csv}"
    if [ "\${legacy}" = "true" ]; then
        echo "WARNING: legacy two-column map; all runs are processed as rep1" >&2
    fi
    for run_id in "\${run_ids[@]}"; do
        THREAD_MEM=$(( ${MEM_MEDIUM%G} / ${THREADS} ))
        bash "\${EPI_SCRIPT_DIR}/ATAC_align.sh" \
            "${FASTQ_DIR}" "${ALIGN_DIR}" "\${run_id}" "${BIOSAMPLE_ID}" \
            "\${SLURM_CPUS_PER_TASK}" "\${THREAD_MEM}" true
    done
    THREAD_MEM=$(( ${MEM_MEDIUM%G} / ${THREADS} ))
    bash "\${EPI_SCRIPT_DIR}/ATAC_finalize_replicate.sh" \
        "${BIOSAMPLE_ID}" "\${replicate}" "\${layout}" "${ALIGN_DIR}" "${TAGALIGN_DIR}" \
        "\${SLURM_CPUS_PER_TASK}" "\${THREAD_MEM}" "\${run_csv}"
    REP_TAGALIGNS+=("${TAGALIGN_DIR}/replicates/${BIOSAMPLE_ID}/${BIOSAMPLE_ID}.\${replicate}.tn5.tagAlign.gz")
    REP_FRAGMENTS+=("${TAGALIGN_DIR}/replicates/${BIOSAMPLE_ID}/${BIOSAMPLE_ID}.\${replicate}.fragments.bed.gz")
done

zcat -f "\${REP_TAGALIGNS[@]}" | gzip -nc > "${TAGALIGN_DIR}/${BIOSAMPLE_ID}.tn5.tagAlign.gz"
zcat -f "\${REP_FRAGMENTS[@]}" | gzip -nc > "${TAGALIGN_DIR}/${BIOSAMPLE_ID}.fragments.bed.gz"
printf "biosample_id\tbiological_replicates\n%s\t%s\n" "${BIOSAMPLE_ID}" "\${#REPLICATE_ROWS[@]}" \
    > "${TAGALIGN_DIR}/${BIOSAMPLE_ID}.replicate_manifest.tsv"
conda deactivate
EOF
    chmod +x "${JOB_SCRIPT_PATH}"
    echo "已生成: ${JOB_SCRIPT_PATH}"
done

echo "所有作业脚本生成完成：${TOTAL_JOBS} 个 Biosample。"
echo "提交命令: for f in ${JOB_DIR}/*.sh; do sbatch \$f; done"
