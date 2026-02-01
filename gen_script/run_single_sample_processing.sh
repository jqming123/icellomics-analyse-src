#!/bin/bash
# ==================================
# Script: run_single_sample_processing.sh
# Description: Process a single sample BAM file to generate gVCF
# Author: gaoxiaojing project
# ==================================

set -eo pipefail

# 假设此脚本在 resources/src/ 目录下运行
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Initialize variables to capture arguments
PROJECT_NAME_ARG=""
REF_GENOME_ARG=""
# --- 参数解析 ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --sampleName) SAMPLE_NAME="$2"; shift 2;;
        --inputBam) INPUT_BAM="$2"; shift 2;;
        --ref) REF_GENOME_ARG="$2"; shift 2;; # Capture reference genome from argument
        --outDir) OUT_DIR="$2"; shift 2;;
        --logDir) LOG_DIR="$2"; shift 2;;
        --projectName) PROJECT_NAME_ARG="$2"; shift 2;; # Add this line to parse PROJECT_NAME
        *) echo "Unknown option $1"; exit 1;;
    esac
done
# Set PROJECT_NAME for config.sh. Prioritize argument, then environment.
if [ -n "${PROJECT_NAME_ARG}" ]; then
    export PROJECT_NAME="${PROJECT_NAME_ARG}"
# else, PROJECT_NAME is expected to be in the environment (exported by parent script)
# config.sh will validate if it's set.
fi

# Source config.sh AFTER PROJECT_NAME is potentially set.
source "${SCRIPT_DIR}/config.sh"

# If REF_GENOME_ARG was provided, it overrides the value from config.sh
if [ -n "${REF_GENOME_ARG}" ]; then
    REF_GENOME="${REF_GENOME_ARG}"
fi

if [[ -z "$SAMPLE_NAME" || -z "$INPUT_BAM" || -z "$REF_GENOME" || -z "$OUT_DIR" || -z "$PROJECT_NAME" ]]; then # Modified validation
    echo "Usage: $0 --sampleName <name> --inputBam <bam> --ref <reference.fna> --outDir <dir> [--logDir <dir>] --projectName <project_name>" # Modified usage message
   exit 1
fi

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

echo "[`date`] Project name: ${PROJECT_NAME}" # Added for clarity

# --- 激活 miniforge conda 环境 ---
# source ${CONDA_PROFILE_PATH}
# conda activate ${MAMBA_ENV_NAME}

echo "[`date`] Processing sample: ${SAMPLE_NAME}"
echo "[`date`] Input BAM: ${INPUT_BAM}"
echo "[`date`] Reference: ${REF_GENOME}"
echo "[`date`] Output dir: ${OUT_DIR}"

# --- 文件路径 ---
DEDUP_BAM="${OUT_DIR}/${SAMPLE_NAME}.dedup.bam"
METRICS_FILE="${OUT_DIR}/${SAMPLE_NAME}.metrics.txt"
GVCF_FILE="${OUT_DIR}/${SAMPLE_NAME}.g.vcf.gz"

# 检查参考基因组对应的.fai 和 .dict文件是否存在
if [ ! -f "${REF_GENOME}.fai" ]; then
  echo "未发现.fai文件...正在生成.fai文件"
  samtools faidx "${REF_GENOME}"
fi
if [ ! -f "${REF_GENOME%.fna}.dict" ] && [ ! -f "${REF_GENOME%.fa}.dict" ]; then
  echo "未发现.dict文件...正在生成.dict文件"
  gatk CreateSequenceDictionary -R "${REF_GENOME}" -O "${REF_GENOME%.*}.dict"
fi

# --- Step 1: MarkDuplicates ---
echo "[`date`] MarkDuplicates..."
gatk MarkDuplicates \
    -I "${INPUT_BAM}" \
    -O "${DEDUP_BAM}" \
    -M "${METRICS_FILE}" \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT

# --- Step 2: HaplotypeCaller to gVCF ---
echo "[`date`] Running HaplotypeCaller..."
gatk HaplotypeCaller \
    -R "${REF_GENOME}" \
    -I "${DEDUP_BAM}" \
    -O "${GVCF_FILE}" \
    -ERC GVCF

echo "[`date`] Finished processing sample: ${SAMPLE_NAME}"

