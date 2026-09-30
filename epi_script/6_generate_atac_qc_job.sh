#!/usr/bin/env bash

# @File       :6_generate_atac_qc_job.sh
# @Description:Generate an ATAC-seq QC summary job. The same entry point
#              supports legacy existing-result QC, early core QC for a new
#              analysis, and release QC.
# @Usage      :bash 6_generate_atac_qc_job.sh <PROJECT_NAME> <REF_NAME> \
#              [--mode existing-results|core|release] \
#              [--run-manifest FILE] [--biosample-manifest FILE]

set -euo pipefail

usage() {
    cat >&2 <<'EOF'
Usage:
  bash 6_generate_atac_qc_job.sh <PROJECT_NAME> <REF_NAME> \
    [--mode existing-results|core|release] \
    [--run-manifest FILE] [--biosample-manifest FILE]

Modes:
  existing-results  Read-only legacy baseline QC. Requires --biosample-manifest;
                    does not require IDR, fingerprint, GC bias, a current
                    sample map, or a configured blacklist.
  core              QC of the four reviewer-facing core metrics for a new run.
                    Release-only artifacts do not block the core result.
  release           Full release QC (default). A configured assembly-matched
                    blacklist and fingerprint are required.

Manifest options are accepted only with existing-results. Core/release always
use the current sample_run_map.tsv and canonical pipeline outputs.
EOF
}

if [ "$#" -lt 2 ]; then
    usage
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="$2"
shift 2

QC_MODE="release"
RUN_MANIFEST=""
BIOSAMPLE_MANIFEST=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        --mode)
            [ "$#" -ge 2 ] || { echo "错误: --mode 需要一个值。" >&2; exit 1; }
            QC_MODE="$2"
            shift 2
            ;;
        --run-manifest)
            [ "$#" -ge 2 ] || { echo "错误: --run-manifest 需要一个文件。" >&2; exit 1; }
            RUN_MANIFEST="$2"
            shift 2
            ;;
        --biosample-manifest)
            [ "$#" -ge 2 ] || { echo "错误: --biosample-manifest 需要一个文件。" >&2; exit 1; }
            BIOSAMPLE_MANIFEST="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "错误: 未识别的参数: $1" >&2
            usage
            exit 1
            ;;
    esac
done

case "${QC_MODE}" in
    existing-results|core|release) ;;
    *)
        echo "错误: --mode 必须是 existing-results、core 或 release。" >&2
        exit 1
        ;;
esac

abspath() {
    local target="$1"
    local parent
    parent="$(dirname "${target}")"
    parent="$(cd "${parent}" && pwd -P)"
    printf '%s/%s\n' "${parent}" "$(basename "${target}")"
}

if [ -n "${RUN_MANIFEST}" ]; then
    [ -f "${RUN_MANIFEST}" ] || { echo "错误: run manifest 不存在: ${RUN_MANIFEST}" >&2; exit 1; }
    RUN_MANIFEST="$(abspath "${RUN_MANIFEST}")"
fi
if [ -n "${BIOSAMPLE_MANIFEST}" ]; then
    [ -f "${BIOSAMPLE_MANIFEST}" ] || { echo "错误: biosample manifest 不存在: ${BIOSAMPLE_MANIFEST}" >&2; exit 1; }
    BIOSAMPLE_MANIFEST="$(abspath "${BIOSAMPLE_MANIFEST}")"
fi
if [ "${QC_MODE}" = "existing-results" ] && [ -z "${BIOSAMPLE_MANIFEST}" ]; then
    echo "错误: existing-results 模式要求 --biosample-manifest，以避免猜测旧结果路径或文件名。" >&2
    exit 1
fi
if [ "${QC_MODE}" != "existing-results" ] && { [ -n "${RUN_MANIFEST}" ] || [ -n "${BIOSAMPLE_MANIFEST}" ]; }; then
    echo "错误: --run-manifest/--biosample-manifest 仅可用于 existing-results；core/release 必须使用当前 sample_run_map.tsv 和标准结果路径。" >&2
    exit 1
fi

# The server job sources this configuration file, not the local checkout copy.
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}" >&2
    exit 1
fi
source "${CONFIG_PATH}"

if [ -z "${PROJECT_DIR:-}" ] || [ -z "${EPI_SCRIPT_DIR:-}" ]; then
    echo "错误: PROJECT_DIR 或 EPI_SCRIPT_DIR 未设置，请检查 epi_config.sh。" >&2
    exit 1
fi
if [ ! -f "${EPI_SCRIPT_DIR}/summarize_atac_qc.py" ]; then
    echo "错误: 未找到 QC 脚本 ${EPI_SCRIPT_DIR}/summarize_atac_qc.py" >&2
    exit 1
fi

echo "正在为项目 ${PROJECT_NAME} (参考基因组: ${REF_NAME}) 生成 ATAC-seq ${QC_MODE} QC 任务..."
echo "当前指定的队列：${QUEUE_NAME}"
echo "内存配置：${MEM_QC}"
echo "时间限制：${QC_WALLTIME}"

JOB_DIR="${PROJECT_DIR}/2_jobs/atac_qc"
LOG_DIR="${PROJECT_DIR}/3_logs/atac_qc"
RESULTS_DIR="${PROJECT_DIR}/1_result"
QC_DIR="${RESULTS_DIR}/5_qc"
mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${QC_DIR}" "${QC_DIR}/artifacts/${QC_MODE}/fragment_size" \
    "${QC_DIR}/artifacts/${QC_MODE}/tss_enrichment" \
    "${QC_DIR}/artifacts/${QC_MODE}/library_complexity" "${QC_DIR}/json/${QC_MODE}"
SUMMARY_PATH="${QC_DIR}/atac_qc_summary_${QC_MODE}.tsv"
TSS_ARTIFACT_LEVEL="${ATAC_TSS_ARTIFACTS:-profile}"

PROJECT_ID="${PROJECT_NAME%%_*}"
JOB_SCRIPT_PATH="${JOB_DIR}/atac_qc_${QC_MODE}_${PROJECT_NAME}.sh"

MANIFEST_ARGUMENTS=""
REPORT_MANIFEST_ARGUMENTS=""
if [ -n "${RUN_MANIFEST}" ]; then
    printf -v quoted '%q' "${RUN_MANIFEST}"
    MANIFEST_ARGUMENTS+=" --run-manifest ${quoted}"
    REPORT_MANIFEST_ARGUMENTS+=" --run-manifest ${quoted}"
fi
if [ -n "${BIOSAMPLE_MANIFEST}" ]; then
    printf -v quoted '%q' "${BIOSAMPLE_MANIFEST}"
    MANIFEST_ARGUMENTS+=" --biosample-manifest ${quoted}"
    REPORT_MANIFEST_ARGUMENTS+=" --biosample-manifest ${quoted}"
fi

# A missing configured blacklist prevents release QC only. Core and legacy
# modes record the absence rather than hiding recoverable reviewer metrics.
SUMMARY_BLACKLIST_ARGUMENT=""
FINGERPRINT_BLACKLIST_ARGUMENT=""
if [ "${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ] && [ ! -s "${BLACKLIST_BED:-}" ]; then
    if [ "${QC_MODE}" = "release" ] && [ "${BLACKLIST_POLICY:-release-required}" = "release-required" ]; then
        echo "错误: ${REF_NAME} 的 BLACKLIST_STATUS=CONFIGURED，但 BLACKLIST_BED 为空、缺失或为空文件。" >&2
        echo "请修复 epi_config.sh 并安装/验证 assembly-matched blacklist 后再生成 release QC 作业。" >&2
        exit 1
    fi
    echo "警告: 配置声明 blacklist 但资源不可用；${QC_MODE} QC 将记录该项不可用。" >&2
elif [ -n "${BLACKLIST_BED:-}" ] && [ -s "${BLACKLIST_BED}" ]; then
    printf -v quoted_blacklist '%q' "${BLACKLIST_BED}"
    # Release summary must resolve the configured path itself; otherwise a
    # direct CLI override could bypass its configuration/validation gate.
    if [ "${QC_MODE}" != "release" ]; then
        SUMMARY_BLACKLIST_ARGUMENT=" --blacklist-bed ${quoted_blacklist}"
    fi
    FINGERPRINT_BLACKLIST_ARGUMENT="--blacklist-bed ${quoted_blacklist}"
elif [ -n "${BLACKLIST_BED:-}" ]; then
    echo "警告: blacklist 不存在或为空；${QC_MODE} QC 将继续，并记录该项不可用。" >&2
fi

# Do the costly-but-small resource integrity check once while generating a
# release job.  It is deliberately not imposed on core or existing-results
# QC, whose purpose is to retain recoverable evidence when release artefacts
# are absent.
if [ "${QC_MODE}" = "release" ] && [ "${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ]; then
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

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/usr/bin/env bash
# Generated by 6_generate_atac_qc_job.sh; do not edit by hand.
# QC mode: ${QC_MODE}
# Existing-result QC is read-only and does not modify BAM, tagAlign, peak, or track files.
#SBATCH --job-name=atacqc_${PROJECT_ID}
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_QC}
#SBATCH --time=${QC_WALLTIME}
#SBATCH --output=${LOG_DIR}/atac_qc_${QC_MODE}_${PROJECT_NAME}_%j.log

# The cluster conda activate/deactivate hooks reference optional variables
# (ZSH_VERSION, JAVA_HOME) without a default, so nounset (-u) aborts the job at
# activation time, before any QC code runs.  Keep errexit and pipefail only.
set -eo pipefail

QC_MODE="${QC_MODE}"
export ATAC_TSS_ARTIFACTS="${TSS_ARTIFACT_LEVEL}"
source "${CONDA_PROFILE_PATH}"
conda activate "${EPI_CONDA_ENV_NAME}"

ALIGN_DIR="${RESULTS_DIR}/1_alignment"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"
MAP_FILE="${PROJECT_DIR}/0_data/sample_run_map.tsv"

if [ "\${QC_MODE}" != "existing-results" ]; then
    python "${EPI_SCRIPT_DIR}/atac_qc/sample_map.py" \\
        --map "\${MAP_FILE}" --project-dir "${PROJECT_DIR}" --validate

    if [ ! -d "\${ALIGN_DIR}" ] || [ ! -d "\${POOLED_DIR}" ]; then
        echo "错误: core/release QC 需要当前 alignment 和 pooled tagAlign 目录。" >&2
        exit 1
    fi
    if ! ls "\${POOLED_DIR}"/*.fragments.bed.gz >/dev/null 2>&1; then
        echo "警告: 未找到 pooled fragment BED；核心 FRiP 将标为 PARTIAL。" >&2
    fi
    if ! ls "\${PEAKS_DIR}"/*_peaks.narrowPeak >/dev/null 2>&1; then
        echo "警告: 未找到 pooled peak；FRiP 与 peak 统计将标为 PARTIAL。" >&2
    fi
fi

if [ "\${QC_MODE}" = "release" ]; then
    mapfile -t BIOSAMPLES < <(python "${EPI_SCRIPT_DIR}/atac_qc/sample_map.py" \\
        --map "\${MAP_FILE}" --list-biosamples)
    for biosample in "\${BIOSAMPLES[@]}"; do
        python "${EPI_SCRIPT_DIR}/compute_atac_fingerprint_qc.py" \\
            --project-dir "${PROJECT_DIR}" --biosample "\${biosample}" \\
            --threads "\${SLURM_CPUS_PER_TASK}" ${FINGERPRINT_BLACKLIST_ARGUMENT}
    done
fi

echo "--- Running ATAC-seq \${QC_MODE} QC summary ---"
python "${EPI_SCRIPT_DIR}/summarize_atac_qc.py" \\
    "${PROJECT_NAME}" "${REF_NAME}" --mode "\${QC_MODE}" \\
    --output "${SUMMARY_PATH}"${MANIFEST_ARGUMENTS}${SUMMARY_BLACKLIST_ARGUMENT}
python "${EPI_SCRIPT_DIR}/generate_atac_qc_report.py" \\
    --project-dir "${PROJECT_DIR}" \\
    --summary "${SUMMARY_PATH}" \\
    --output "${QC_DIR}/atac_qc_report_\${QC_MODE}.md" \\
    --mode "\${QC_MODE}"${REPORT_MANIFEST_ARGUMENTS}

# Keep the historic release summary path for downstream callers. Core and
# existing-result summaries deliberately remain mode-qualified so a baseline
# cannot be overwritten by a later new-analysis run.
if [ "\${QC_MODE}" = "release" ]; then
    cp "${SUMMARY_PATH}" "${QC_DIR}/atac_qc_summary.tsv"
fi

conda deactivate
echo "QC summary: ${SUMMARY_PATH}"
echo "QC report: ${QC_DIR}/atac_qc_report_\${QC_MODE}.md"
EOF

chmod +x "${JOB_SCRIPT_PATH}"
echo "任务脚本生成完毕: ${JOB_SCRIPT_PATH}"
echo "提交命令: sbatch ${JOB_SCRIPT_PATH}"
