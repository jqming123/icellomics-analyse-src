#!/usr/bin/env bash
#SBATCH -p corexd192
#SBATCH -J build_Dog_E_UUGSD_hisat2_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 10
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_hisat2_index.log

set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
  mamba activate scRNA_env
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD"
FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
OUT_PREFIX="${GENOME_DIR}/hisat2_index/genome"

if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: FASTA not found: ${FASTA}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUT_PREFIX}")"
hisat2-build -p 10 "${FASTA}" "${OUT_PREFIX}"
echo "HISAT2 index built at: ${OUT_PREFIX}.*"
