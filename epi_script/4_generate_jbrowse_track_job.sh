#!/usr/bin/env bash

# @File       :4_generate_jbrowse_track_job.sh
# @Description:Generate one SLURM job script per Biosample to convert verified
#              release ATAC-seq peak outputs into public JBrowse 2 track files.
# @Usage      :bash 4_generate_jbrowse_track_job.sh <PROJECT_NAME> <REF_NAME>
# @Example    :bash 4_generate_jbrowse_track_job.sh PRJNA728969_HEK293 hg38_Ensembl
# Generated conversion jobs are I/O-bound and use a single signal stream, so
# they deliberately request a small two-core allocation rather than peak-call resources.
# Core/diagnostic peak outputs are intentionally not published to this shared
# JBrowse namespace.

if [ "$#" -ne 2 ]; then
    echo "Error: invalid argument count."
    echo "Usage: bash 4_generate_jbrowse_track_job.sh <PROJECT_NAME> <REF_NAME>"
    exit 1
fi

export PROJECT_NAME="$1"
export REF_NAME="$2"

CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "Error: config file not found at ${CONFIG_PATH}"
    exit 1
fi
source "${CONFIG_PATH}"

if [ -z "${JBROWSE_DATA_ROOT:-}" ]; then
    echo "Error: JBROWSE_DATA_ROOT is not defined in ${CONFIG_PATH}"
    exit 1
fi

if [ -z "${SPECIES_ID:-}" ] || [ -z "${ASSEMBLY_NAME:-}" ] || [ "${ASSEMBLY_NAME}" = "NONE" ]; then
    echo "Error: SPECIES_ID/ASSEMBLY_NAME are not available for REF_NAME='${REF_NAME}' in ${CONFIG_PATH}"
    exit 1
fi

PROJECT_ID="${PROJECT_NAME%%_*}"
CELL_LINE="${PROJECT_NAME#*_}"

JOB_DIR="${PROJECT_DIR}/2_jobs/jbrowse_tracks"
LOG_DIR="${PROJECT_DIR}/3_logs/jbrowse_tracks"
RESULTS_DIR="${PROJECT_DIR}/1_result"
POOLED_DIR="${RESULTS_DIR}/2_tagalign"
PEAKS_DIR="${RESULTS_DIR}/3_peak_calling"

mkdir -p "${JOB_DIR}" "${LOG_DIR}"

FOUND_ANY=0

for PEAK_FILE in "${PEAKS_DIR}"/*_peaks.narrowPeak; do
    if [ ! -f "${PEAK_FILE}" ]; then
        continue
    fi

    FOUND_ANY=1
    FILE_NAME="$(basename "${PEAK_FILE}")"
    BIOSAMPLE_NAME="${FILE_NAME%%_peaks.narrowPeak}"
    BEDGRAPH_FILE="${PEAKS_DIR}/${BIOSAMPLE_NAME}_treat_pileup.bdg"
    TAGALIGN_FILE="${POOLED_DIR}/${BIOSAMPLE_NAME}.tn5.tagAlign.gz"
    REPRO_PEAK_FILE="${RESULTS_DIR}/4_reproducibility/${BIOSAMPLE_NAME}/${BIOSAMPLE_NAME}.reproducible_peaks.narrowPeak"
    REPRO_JSON="${RESULTS_DIR}/4_reproducibility/${BIOSAMPLE_NAME}/${BIOSAMPLE_NAME}.reproducibility.json"
    DATASET_ID="${PROJECT_ID}_${BIOSAMPLE_NAME}"
    JOB_SCRIPT_PATH="${JOB_DIR}/jbrowse_track_${BIOSAMPLE_NAME}.sh"

    cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=jb2trk_${BIOSAMPLE_NAME}
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_SMALL}
#SBATCH --time=1-00:00:00
#SBATCH --output=${LOG_DIR}/jbrowse_track_${BIOSAMPLE_NAME}_%j.log

set -euo pipefail

echo "=========================================================="
echo "Job started on \$(date)"
echo "Project name: ${PROJECT_NAME}"
echo "Biosample: ${BIOSAMPLE_NAME}"
echo "Dataset ID: ${DATASET_ID}"
echo "=========================================================="

export PROJECT_NAME="${PROJECT_NAME}"
export REF_NAME="${REF_NAME}"
source "${CONFIG_PATH}"

set +u
source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}"
set -u

JBROWSE_DATA_ROOT="${JBROWSE_DATA_ROOT}"
ASSEMBLY_DIR="\${JBROWSE_DATA_ROOT}/assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}"
TRACK_DIR="\${JBROWSE_DATA_ROOT}/tracks/atac/${CELL_LINE}/${DATASET_ID}"
HTTP_BASE="\${JBROWSE_HTTP_BASE:-}"
CHROM_SIZES="\${ASSEMBLY_DIR}/chrom.sizes"
SORT_THREADS="\${SLURM_CPUS_PER_TASK:-1}"

require_cmd() {
    command -v "\$1" >/dev/null 2>&1 || {
        echo "Error: required command not found: \$1" >&2
        exit 1
    }
}

require_cmd bgzip
require_cmd tabix
require_cmd bedGraphToBigWig
require_cmd sha256sum

if [ ! -f "\${CHROM_SIZES}" ]; then
    echo "Error: chrom.sizes not found: \${CHROM_SIZES}" >&2
    echo "Run 3_generate_jbrowse_assembly_job.sh first." >&2
    exit 1
fi

TMP_DIR=\$(mktemp -d)
trap 'rm -rf "\${TMP_DIR}"' EXIT

SIGNAL_MODE=""
SIGNAL_INPUT=""
SIGNAL_SOURCE_MODE=""
SIGNAL_BLACKLIST_APPLIED_JSON="false"
SIGNAL_BLACKLIST_STATUS="NOT_APPLICABLE"
SIGNAL_BLACKLIST_METHOD="not_applicable"
SIGNAL_BLACKLIST_SCOPE="not_applicable"
BLACKLIST_SCOPE="peaks_only_no_blacklist_resource"
mapfile -t PEAK_META < <(python "\${EPI_SCRIPT_DIR}/atac_qc/peak_provenance.py" \
    --pooled-peak "${PEAK_FILE}" \
    --reproducible-peak "${REPRO_PEAK_FILE}" \
    --reproducibility-json "${REPRO_JSON}" \
    --emit-lines)
if [ "\${#PEAK_META[@]}" -ne 10 ]; then
    echo "Error: failed to resolve peak/blacklist provenance for ${BIOSAMPLE_NAME}" >&2
    exit 1
fi
ACTIVE_PEAK_FILE="\${PEAK_META[0]}"
PEAK_MODE="\${PEAK_META[1]}"
BLACKLIST_APPLIED="\${PEAK_META[2]}"
BLACKLIST_STATUS_USED="\${PEAK_META[3]}"
BLACKLIST_BED_USED="\${PEAK_META[4]}"
BLACKLIST_SHA256_USED="\${PEAK_META[5]}"
BLACKLIST_SOURCE_USED="\${PEAK_META[6]}"
PEAK_PIPELINE_MODE="\${PEAK_META[7]}"
BLACKLIST_VALIDATION_STATUS="\${PEAK_META[8]}"
BLACKLIST_VALIDATION_SHA256="\${PEAK_META[9]}"
case "\${BLACKLIST_APPLIED}" in
    true|false) BLACKLIST_APPLIED_JSON="\${BLACKLIST_APPLIED}" ;;
    *) BLACKLIST_APPLIED_JSON="null" ;;
esac
if [ ! -s "\${ACTIVE_PEAK_FILE}" ]; then
    echo "Error: selected peak file is missing or empty: \${ACTIVE_PEAK_FILE}" >&2
    exit 1
fi
if [ "\${PEAK_PIPELINE_MODE}" != "release" ]; then
    echo "Error: public JBrowse tracks require a release peak sidecar; observed mode '\${PEAK_PIPELINE_MODE:-missing}'." >&2
    exit 1
fi
if [ "\${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ]; then
    if [ "\${BLACKLIST_APPLIED}" != "true" ] || [ "\${BLACKLIST_STATUS_USED}" != "APPLIED" ] || [ -z "\${BLACKLIST_SHA256_USED}" ]; then
        echo "Error: configured blacklist was not proven applied to the public peak track." >&2
        exit 1
    fi
    if [ "\${BLACKLIST_VALIDATION_STATUS}" != "VERIFIED" ] || [ "\${BLACKLIST_VALIDATION_SHA256}" != "\${BLACKLIST_SHA256_USED}" ]; then
        echo "Error: public peak blacklist lacks a matching VERIFIED validation record." >&2
        exit 1
    fi
    if [ -z "\${BLACKLIST_BED:-}" ] || [ ! -s "\${BLACKLIST_BED}" ]; then
        echo "Error: configured blacklist is unavailable while preparing a public track." >&2
        exit 1
    fi
    CURRENT_BLACKLIST_SHA256=\$(sha256sum "\${BLACKLIST_BED}" | awk '{print \$1}')
    if [ "\${CURRENT_BLACKLIST_SHA256}" != "\${BLACKLIST_SHA256_USED}" ]; then
        echo "Error: configured blacklist changed after release peak calling; regenerate release peaks." >&2
        exit 1
    fi
    require_cmd bedtools
    SIGNAL_BLACKLIST_APPLIED_JSON="true"
    SIGNAL_BLACKLIST_STATUS="MASKED_INTERVALS"
    SIGNAL_BLACKLIST_METHOD="bedtools_subtract_coordinate_mask"
    SIGNAL_BLACKLIST_SCOPE="masked_intervals"
    BLACKLIST_SCOPE="peaks_and_signal_mask"
else
    SIGNAL_BLACKLIST_STATUS="NOT_APPLICABLE_\${BLACKLIST_STATUS:-NOT_CONFIGURED}"
fi

mkdir -p "\${TRACK_DIR}"

if [ -f "${BEDGRAPH_FILE}" ]; then
    SIGNAL_SOURCE_MODE="macs3_treat_pileup_bdg"
    SIGNAL_INPUT="${BEDGRAPH_FILE}"
    sort --parallel="\${SORT_THREADS}" -k1,1 -k2,2n "${BEDGRAPH_FILE}" > "\${TMP_DIR}/signal.sorted.bdg"
elif [ -f "${TAGALIGN_FILE}" ]; then
    SIGNAL_SOURCE_MODE="pooled_tagalign_fallback"
    SIGNAL_INPUT="${TAGALIGN_FILE}"
    require_cmd bedtools
    zcat "${TAGALIGN_FILE}" | sort --parallel="\${SORT_THREADS}" -k1,1 -k2,2n > "\${TMP_DIR}/tagalign.sorted.bed"
    bedtools genomecov -bg -i "\${TMP_DIR}/tagalign.sorted.bed" -g "\${CHROM_SIZES}" > "\${TMP_DIR}/signal.sorted.bdg"
else
    echo "Error: neither bedGraph nor pooled tagAlign exists for ${BIOSAMPLE_NAME}" >&2
    exit 1
fi

if command -v bedClip >/dev/null 2>&1; then
    bedClip "\${TMP_DIR}/signal.sorted.bdg" "\${CHROM_SIZES}" "\${TMP_DIR}/signal.clipped.bdg"
    SIGNAL_BDG="\${TMP_DIR}/signal.clipped.bdg"
else
    SIGNAL_BDG="\${TMP_DIR}/signal.sorted.bdg"
fi

if [ "\${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ]; then
    # Mask the display coordinate only.  Reads are not removed and the source
    # pooled MACS3 bedGraph remains in the analysis result directory.
    bedtools subtract -a "\${SIGNAL_BDG}" -b "\${BLACKLIST_BED}" > "\${TMP_DIR}/signal.blacklist_masked.bdg"
    if [ ! -s "\${TMP_DIR}/signal.blacklist_masked.bdg" ]; then
        echo "Error: blacklist masking removed all signal intervals for ${BIOSAMPLE_NAME}." >&2
        exit 1
    fi
    SIGNAL_BDG="\${TMP_DIR}/signal.blacklist_masked.bdg"
fi

if [ "\${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ]; then
    SIGNAL_MODE="blacklist_masked_\${SIGNAL_SOURCE_MODE}"
else
    SIGNAL_MODE="\${SIGNAL_SOURCE_MODE}_no_blacklist_resource"
fi
bedGraphToBigWig "\${SIGNAL_BDG}" "\${CHROM_SIZES}" "\${TMP_DIR}/signal.bw"

sort --parallel="\${SORT_THREADS}" -k1,1 -k2,2n "\${ACTIVE_PEAK_FILE}" > "\${TMP_DIR}/peaks.sorted.narrowPeak"
bgzip -f -c "\${TMP_DIR}/peaks.sorted.narrowPeak" > "\${TMP_DIR}/peaks.narrowPeak.gz"
tabix -f -p bed "\${TMP_DIR}/peaks.narrowPeak.gz"

SUMMIT_ARGS=()
if [ "\${BLACKLIST_STATUS:-NOT_CONFIGURED}" = "CONFIGURED" ]; then
    SUMMIT_ARGS=(--blacklist-bed "\${BLACKLIST_BED}")
fi
python "\${EPI_SCRIPT_DIR}/create_atac_release_summits.py" \
    --narrowpeak "\${ACTIVE_PEAK_FILE}" \
    --output "\${TMP_DIR}/summits.release.bed" \
    "\${SUMMIT_ARGS[@]}"
sort --parallel="\${SORT_THREADS}" -k1,1 -k2,2n "\${TMP_DIR}/summits.release.bed" > "\${TMP_DIR}/summits.sorted.bed"
SUMMIT_COUNT=\$(wc -l < "\${TMP_DIR}/summits.sorted.bed")
bgzip -f -c "\${TMP_DIR}/summits.sorted.bed" > "\${TMP_DIR}/summits.bed.gz"
tabix -f -p bed "\${TMP_DIR}/summits.bed.gz"

cat > "\${TMP_DIR}/track.meta.json" <<JSON
{
  "dataset": "${DATASET_ID}",
  "project_id": "${PROJECT_ID}",
  "cell_line": "${CELL_LINE}",
  "sample_id": "${BIOSAMPLE_NAME}",
  "species_id": ${SPECIES_ID},
  "assembly_name": "${ASSEMBLY_NAME}",
  "ref_name": "${REF_NAME}",
  "signal_mode": "\${SIGNAL_MODE}",
  "signal_source_mode": "\${SIGNAL_SOURCE_MODE}",
  "signal_input": "\${SIGNAL_INPUT}",
  "signal_blacklist_applied": \${SIGNAL_BLACKLIST_APPLIED_JSON},
  "signal_blacklist_status": "\${SIGNAL_BLACKLIST_STATUS}",
  "signal_blacklist_method": "\${SIGNAL_BLACKLIST_METHOD}",
  "signal_blacklist_scope": "\${SIGNAL_BLACKLIST_SCOPE}",
  "peak_mode": "\${PEAK_MODE}",
  "peak_pipeline_mode": "\${PEAK_PIPELINE_MODE}",
  "peak_input": "\${ACTIVE_PEAK_FILE}",
  "peak_blacklist_applied": \${BLACKLIST_APPLIED_JSON},
  "peak_blacklist_status": "\${BLACKLIST_STATUS_USED}",
  "peak_blacklist_validation_status": "\${BLACKLIST_VALIDATION_STATUS}",
  "peak_blacklist_validation_sha256": "\${BLACKLIST_VALIDATION_SHA256}",
  "blacklist_scope": "\${BLACKLIST_SCOPE}",
  "blacklist_applied": \${BLACKLIST_APPLIED_JSON},
  "blacklist_status": "\${BLACKLIST_STATUS_USED}",
  "blacklist_bed": "\${BLACKLIST_BED_USED}",
  "blacklist_sha256": "\${BLACKLIST_SHA256_USED}",
  "blacklist_source": "\${BLACKLIST_SOURCE_USED}",
  "signal_bw_path": "tracks/atac/${CELL_LINE}/${DATASET_ID}/signal.bw",
  "peaks_path": "tracks/atac/${CELL_LINE}/${DATASET_ID}/peaks.narrowPeak.gz",
  "peaks_index_path": "tracks/atac/${CELL_LINE}/${DATASET_ID}/peaks.narrowPeak.gz.tbi",
  "summits_path": "tracks/atac/${CELL_LINE}/${DATASET_ID}/summits.bed.gz",
  "summits_index_path": "tracks/atac/${CELL_LINE}/${DATASET_ID}/summits.bed.gz.tbi",
  "summits_status": "DERIVED_FROM_FINAL_RELEASE_NARROWPEAK",
  "summits_source": "final_release_narrowPeak_column10_peak_offset",
  "summit_count": \${SUMMIT_COUNT},
  "publication_status": "RELEASE_READY",
  "http_base": "\${HTTP_BASE}"
}
JSON

# Do not alter an existing public track until every replacement payload has
# been built and indexed successfully.
mv -f "\${TMP_DIR}/signal.bw" "\${TRACK_DIR}/signal.bw"
mv -f "\${TMP_DIR}/peaks.narrowPeak.gz" "\${TRACK_DIR}/peaks.narrowPeak.gz"
mv -f "\${TMP_DIR}/peaks.narrowPeak.gz.tbi" "\${TRACK_DIR}/peaks.narrowPeak.gz.tbi"
mv -f "\${TMP_DIR}/summits.bed.gz" "\${TRACK_DIR}/summits.bed.gz"
mv -f "\${TMP_DIR}/summits.bed.gz.tbi" "\${TRACK_DIR}/summits.bed.gz.tbi"
mv -f "\${TMP_DIR}/track.meta.json" "\${TRACK_DIR}/track.meta.json"

set +u
conda deactivate

echo "=========================================================="
echo "Track directory prepared at: \${TRACK_DIR}"
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

    chmod +x "${JOB_SCRIPT_PATH}"
done

if [ "${FOUND_ANY}" -eq 0 ]; then
    echo "No *_peaks.narrowPeak files were found in ${PEAKS_DIR}"
    exit 0
fi

echo "Done. Generated JBrowse 2 track job scripts in:"
echo "  ${JOB_DIR}"
echo ""
echo "Submit them in batch with:"
echo "  for f in ${JOB_DIR}/*.sh; do sbatch \$f; done"
