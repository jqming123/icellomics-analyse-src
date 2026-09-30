#!/usr/bin/env bash

# @File       :5_generate_jbrowse_manifest_job.sh
# @Description:Generate one SLURM job script that scans prepared JBrowse 2
#              track directories and writes a manifest TSV for later backend
#              registration or static-config generation.
# @Usage      :bash 5_generate_jbrowse_manifest_job.sh <PROJECT_NAME> <REF_NAME>
# @Example    :bash 5_generate_jbrowse_manifest_job.sh PRJNA728969_HEK293 hg38_Ensembl

if [ "$#" -ne 2 ]; then
    echo "Error: invalid argument count."
    echo "Usage: bash 5_generate_jbrowse_manifest_job.sh <PROJECT_NAME> <REF_NAME>"
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

JOB_DIR="${PROJECT_DIR}/2_jobs/jbrowse_manifest"
LOG_DIR="${PROJECT_DIR}/3_logs/jbrowse_manifest"
MANIFEST_DIR="${PROJECT_DIR}/1_result/4_jbrowse_manifest"

mkdir -p "${JOB_DIR}" "${LOG_DIR}" "${MANIFEST_DIR}"

JOB_SCRIPT_PATH="${JOB_DIR}/build_jbrowse_manifest_${PROJECT_NAME}.sh"

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=jb2mft_${PROJECT_ID}
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MEM_MEDIUM}
#SBATCH --time=1-00:00:00
#SBATCH --output=${LOG_DIR}/build_jbrowse_manifest_${PROJECT_NAME}_%j.log

set -euo pipefail

echo "=========================================================="
echo "Job started on \$(date)"
echo "Project name: ${PROJECT_NAME}"
echo "=========================================================="

TRACK_BASE="${JBROWSE_DATA_ROOT}/tracks/atac/${CELL_LINE}"
MANIFEST_TSV="${MANIFEST_DIR}/atac_jbrowse_track_manifest.${PROJECT_NAME}.tsv"
ASSEMBLY_META="${JBROWSE_DATA_ROOT}/assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/assembly.meta.json"

mkdir -p "${MANIFEST_DIR}"

if [ ! -f "\${ASSEMBLY_META}" ]; then
    echo "Error: assembly metadata not found: \${ASSEMBLY_META}" >&2
    echo "Run 3_generate_jbrowse_assembly_job.sh first." >&2
    exit 1
fi

if [ ! -d "\${TRACK_BASE}" ]; then
    echo "Error: track base directory not found: \${TRACK_BASE}" >&2
    echo "Run 4_generate_jbrowse_track_job.sh and submit the generated jobs first." >&2
    exit 1
fi

{
    printf 'dataset\tproject_id\tcell_line\tsample_id\tspecies_id\tassembly_name\tref_name\tsignal_bw_path\tpeaks_path\tpeaks_index_path\tsummits_path\tsummits_index_path\ttrack_meta_path\tassembly_meta_path\n'

    shopt -s nullglob
    for TRACK_DIR in "\${TRACK_BASE}"/${PROJECT_ID}_*; do
        [ -d "\${TRACK_DIR}" ] || continue

        DATASET_ID=\$(basename "\${TRACK_DIR}")
        SAMPLE_ID="\${DATASET_ID#${PROJECT_ID}_}"

        SIGNAL_PATH="tracks/atac/${CELL_LINE}/\${DATASET_ID}/signal.bw"
        PEAKS_PATH="tracks/atac/${CELL_LINE}/\${DATASET_ID}/peaks.narrowPeak.gz"
        PEAKS_INDEX_PATH="tracks/atac/${CELL_LINE}/\${DATASET_ID}/peaks.narrowPeak.gz.tbi"
        SUMMITS_PATH="tracks/atac/${CELL_LINE}/\${DATASET_ID}/summits.bed.gz"
        SUMMITS_INDEX_PATH="tracks/atac/${CELL_LINE}/\${DATASET_ID}/summits.bed.gz.tbi"
        TRACK_META_PATH="tracks/atac/${CELL_LINE}/\${DATASET_ID}/track.meta.json"
        ASSEMBLY_META_PATH="assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/assembly.meta.json"

        if [ ! -f "\${TRACK_DIR}/signal.bw" ] || [ ! -f "\${TRACK_DIR}/peaks.narrowPeak.gz" ] || [ ! -f "\${TRACK_DIR}/peaks.narrowPeak.gz.tbi" ] || [ ! -f "\${TRACK_DIR}/summits.bed.gz" ] || [ ! -f "\${TRACK_DIR}/summits.bed.gz.tbi" ] || [ ! -f "\${TRACK_DIR}/track.meta.json" ]; then
            echo "Skip incomplete track directory: \${TRACK_DIR}" >&2
            continue
        fi

        # Public manifests accept only the complete release contract from
        # step 4: final peaks, blacklist-masked signal where a blacklist is
        # configured, and summits derived from the final release narrowPeak.
        if ! python -c '
import json, sys
data = json.load(open(sys.argv[1], encoding="utf-8"))
release_ok = data.get("publication_status") == "RELEASE_READY" and data.get("peak_pipeline_mode") == "release"
summit_ok = (
    data.get("summits_status") == "DERIVED_FROM_FINAL_RELEASE_NARROWPEAK"
    and data.get("summits_source") == "final_release_narrowPeak_column10_peak_offset"
    and bool(data.get("summits_path"))
    and bool(data.get("summits_index_path"))
)
peak_blacklist = data.get("peak_blacklist_applied")
if peak_blacklist is True:
    signal_ok = (
        data.get("signal_blacklist_applied") is True
        and data.get("signal_blacklist_status") == "MASKED_INTERVALS"
        and data.get("signal_blacklist_method") == "bedtools_subtract_coordinate_mask"
        and data.get("signal_blacklist_scope") == "masked_intervals"
    )
else:
    signal_ok = (
        peak_blacklist is False
        and data.get("signal_blacklist_applied") is False
        and str(data.get("signal_blacklist_status", "")).startswith("NOT_APPLICABLE_")
    )
sys.exit(0 if release_ok and summit_ok and signal_ok else 1)
' "\${TRACK_DIR}/track.meta.json"; then
            echo "Skip track lacking the current release peak/signal/summit contract: \${TRACK_DIR}" >&2
            continue
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "\${DATASET_ID}" \
            "${PROJECT_ID}" \
            "${CELL_LINE}" \
            "\${SAMPLE_ID}" \
            "${SPECIES_ID}" \
            "${ASSEMBLY_NAME}" \
            "${REF_NAME}" \
            "\${SIGNAL_PATH}" \
            "\${PEAKS_PATH}" \
            "\${PEAKS_INDEX_PATH}" \
            "\${SUMMITS_PATH}" \
            "\${SUMMITS_INDEX_PATH}" \
            "\${TRACK_META_PATH}" \
            "\${ASSEMBLY_META_PATH}"
    done
} > "\${MANIFEST_TSV}"

echo "Manifest written to: \${MANIFEST_TSV}"
echo "Assembly meta used:  \${ASSEMBLY_META}"
echo "Job finished on \$(date)"
EOF

chmod +x "${JOB_SCRIPT_PATH}"

echo "Done. Generated JBrowse 2 manifest job script:"
echo "  ${JOB_SCRIPT_PATH}"
echo ""
echo "Submit it with:"
echo "  sbatch ${JOB_SCRIPT_PATH}"
