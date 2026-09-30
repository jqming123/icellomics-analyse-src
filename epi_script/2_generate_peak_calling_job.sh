#!/usr/bin/env bash

# Generate pooled/replicate MACS3 and optional IDR jobs, one per biosample.
# Usage: bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME>
#        [--mode core|release]
#        [--reproducibility-mode auto|pooled_only|true_rep_only|full_idr]

set -euo pipefail
if [ "$#" -lt 2 ]; then
    echo "用法: bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME> [--mode core|release] [--reproducibility-mode auto|pooled_only|true_rep_only|full_idr]" >&2
    exit 1
fi
export PROJECT_NAME="$1"
export REF_NAME="$2"
shift 2

PEAK_MODE="release"
REPRO_MODE_OVERRIDE=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --mode)
            [ "$#" -ge 2 ] || { echo "错误: --mode 需要一个值。" >&2; exit 1; }
            PEAK_MODE="$2"
            shift 2
            ;;
        --reproducibility-mode)
            [ "$#" -ge 2 ] || { echo "错误: --reproducibility-mode 需要一个值。" >&2; exit 1; }
            REPRO_MODE_OVERRIDE="$2"
            shift 2
            ;;
        -h|--help)
            echo "用法: bash 2_generate_peak_calling_job.sh <PROJECT_NAME> <REF_NAME> [--mode core|release] [--reproducibility-mode auto|pooled_only|true_rep_only|full_idr]"
            exit 0
            ;;
        *)
            echo "错误: 未识别的参数: $1" >&2
            exit 1
            ;;
    esac
done

case "${PEAK_MODE}" in
    core|release) ;;
    *) echo "错误: --mode 只能是 core 或 release。" >&2; exit 1 ;;
esac
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}" >&2
    exit 1
fi
source "${CONFIG_PATH}"

REPRO_MODE="${REPRO_MODE_OVERRIDE:-${REPRODUCIBILITY_MODE:-auto}}"
case "${REPRO_MODE}" in
    auto|pooled_only|true_rep_only|full_idr) ;;
    *) echo "错误: 无效的 reproducibility mode: ${REPRO_MODE}" >&2; exit 1 ;;
esac

MAP_FILE="${PROJECT_DIR}/0_data/sample_run_map.tsv"
MAP_PARSER="${EPI_SCRIPT_DIR}/atac_qc/sample_map.py"
python "${MAP_PARSER}" --map "${MAP_FILE}" --project-dir "${PROJECT_DIR}" --validate
mapfile -t BIOSAMPLES < <(python "${MAP_PARSER}" --map "${MAP_FILE}" --list-biosamples)

JOB_DIR="${PROJECT_DIR}/2_jobs/peak_reproducibility"
LOG_DIR="${PROJECT_DIR}/3_logs/peak_reproducibility"
mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${PROJECT_DIR}/1_result/3_peak_calling" \
    "${PROJECT_DIR}/1_result/4_reproducibility"

EFFECTIVE_BLACKLIST_STATUS="${BLACKLIST_STATUS:-NOT_CONFIGURED}"
if [ "${EFFECTIVE_BLACKLIST_STATUS}" = "CONFIGURED" ] && [ ! -s "${BLACKLIST_BED:-}" ]; then
    if [ "${PEAK_MODE}" = "release" ] && [ "${BLACKLIST_POLICY:-release-required}" = "release-required" ]; then
        echo "错误: ${REF_NAME} 的 BLACKLIST_STATUS=CONFIGURED，但 BLACKLIST_BED 为空、缺失或为空文件。" >&2
        echo "请修复 epi_config.sh 并安装/验证 assembly-matched blacklist 后再生成 release peak 作业。" >&2
        exit 1
    fi
    echo "警告: 配置声明 blacklist 但资源不可用；core peak 作业将保留 unfiltered peaks 并记录 MISSING。" >&2
    EFFECTIVE_BLACKLIST_STATUS="MISSING"
elif [ -n "${BLACKLIST_BED:-}" ] && [ -s "${BLACKLIST_BED}" ]; then
    # Validate the configured resource once at release-job generation, not
    # once per biosample.  Core QC intentionally remains recoverable.
    if [ "${PEAK_MODE}" = "release" ] && [ "${EFFECTIVE_BLACKLIST_STATUS}" = "CONFIGURED" ]; then
        VALIDATOR="${EPI_SCRIPT_DIR}/validate_atac_blacklist_resource.py"
        if [ ! -f "${VALIDATOR}" ] || [ ! -s "${REF_GENOME}.fai" ]; then
            echo "错误: release blacklist 验证需要 ${VALIDATOR} 和 ${REF_GENOME}.fai。" >&2
            exit 1
        fi
        VALIDATOR_ARGS=(python "${VALIDATOR}" --bed "${BLACKLIST_BED}" --reference-fai "${REF_GENOME}.fai")
        if [ "${REF_NAME}" = "Mouse_E_GRCm39" ]; then
            VALIDATOR_ARGS+=(--expected-intervals 3147 --require-ensembl-grcm39-contigs)
        fi
        "${VALIDATOR_ARGS[@]}"
    fi
elif [ -n "${BLACKLIST_BED:-}" ]; then
    # An unconfigured optional path may be useful in core work, but must not
    # silently claim that a release resource was applied.
    echo "警告: blacklist 文件不存在或为空；将不传递给 peak caller。" >&2
    EFFECTIVE_BLACKLIST_STATUS="MISSING"
fi

for BIOSAMPLE in "${BIOSAMPLES[@]}"; do
    JOB_SCRIPT="${JOB_DIR}/peak_reproducibility_${BIOSAMPLE}.sh"
    cat > "${JOB_SCRIPT}" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=${BIOSAMPLE}_peak_${PEAK_MODE}
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=7-00:00:00
#SBATCH --output=${LOG_DIR}/peak_reproducibility_${BIOSAMPLE}_%j.log

# nounset (-u) is deliberately not enabled here: the cluster conda
# activate/deactivate hooks reference optional variables (ZSH_VERSION,
# JAVA_HOME) without a default and would abort the job during activation.
set -eo pipefail
export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="${REF_NAME}"
source "${CONFIG_PATH}"
source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}"

command -v macs3 >/dev/null 2>&1 || { echo "ERROR: missing macs3" >&2; exit 1; }
BLACKLIST_RUNTIME_ARGS=()
if [ -n "\${BLACKLIST_BED:-}" ] && [ -s "\${BLACKLIST_BED}" ]; then
    BLACKLIST_RUNTIME_ARGS+=(--blacklist-bed "\${BLACKLIST_BED}")
fi
BLACKLIST_VALIDATION_ARGS=()
if [ "${PEAK_MODE}" = "release" ] && [ "\${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ]; then
    if [ -z "\${BLACKLIST_BED:-}" ] || [ ! -s "\${BLACKLIST_BED}" ]; then
        echo "ERROR: configured release blacklist is missing at runtime" >&2
        exit 1
    fi
    VALIDATOR="\${EPI_SCRIPT_DIR}/validate_atac_blacklist_resource.py"
    if [ ! -f "\${VALIDATOR}" ] || [ ! -s "\${REF_GENOME}.fai" ]; then
        echo "ERROR: release blacklist validation requires \${VALIDATOR} and \${REF_GENOME}.fai" >&2
        exit 1
    fi
    VALIDATION_RECORD="${PROJECT_DIR}/1_result/4_reproducibility/${BIOSAMPLE}/${BIOSAMPLE}.blacklist.validation.json"
    mkdir -p "${PROJECT_DIR}/1_result/4_reproducibility/${BIOSAMPLE}"
    VALIDATOR_ARGS=(python "\${VALIDATOR}" --bed "\${BLACKLIST_BED}" --reference-fai "\${REF_GENOME}.fai")
    if [ "${REF_NAME}" = "Mouse_E_GRCm39" ]; then
        VALIDATOR_ARGS+=(--expected-intervals 3147 --require-ensembl-grcm39-contigs)
    fi
    "\${VALIDATOR_ARGS[@]}" --output "\${VALIDATION_RECORD}"
    BLACKLIST_VALIDATION_ARGS=(--reference-fai "\${REF_GENOME}.fai" --blacklist-validation-json "\${VALIDATION_RECORD}")
fi
python "\${EPI_SCRIPT_DIR}/run_atac_reproducibility.py" \
    --project-dir "${PROJECT_DIR}" \
    --biosample "${BIOSAMPLE}" \
    --gsize "${GSIZE}" \
    --cap-peaks "${PEAK_CAP}" \
    --idr-threshold 0.05 \
    --reproducibility-mode "${REPRO_MODE}" \
    --peak-mode "${PEAK_MODE}" \
    --blacklist-status "\${BLACKLIST_STATUS:-NOT_CONFIGURED}" \
    --blacklist-source "\${BLACKLIST_SOURCE:-}" \
    "\${BLACKLIST_RUNTIME_ARGS[@]}" "\${BLACKLIST_VALIDATION_ARGS[@]}"
conda deactivate
EOF
    chmod +x "${JOB_SCRIPT}"
    echo "已生成: ${JOB_SCRIPT}"
done

echo "已生成 ${#BIOSAMPLES[@]} 个 peak/IDR 作业。"
echo "提交命令: for f in ${JOB_DIR}/*.sh; do sbatch \$f; done"
