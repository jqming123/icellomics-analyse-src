#!/usr/bin/env bash
#SBATCH -p corexd192
#SBATCH -J ref2bed_Dog_E_UUGSD_scRNA
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --time=200:00:00
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_ref2bed.log

set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
  mamba activate scRNA_env
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD"
REF_GTF="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.115.gtf"
REF_BED="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.115.for_scRNA.ref.bed"

if [[ ! -f "${REF_GTF}" ]]; then
  echo "ERROR: GTF not found: ${REF_GTF}" >&2
  exit 1
fi

echo "================================================="
echo "REF2BED START: $(date)"
echo "Input GTF: ${REF_GTF}"
echo "Output BED: ${REF_BED}"
echo "================================================="

gtf2bed < "${REF_GTF}" > "${REF_BED}"

echo "### REF2BED COMPLETED: $(date) ###"
