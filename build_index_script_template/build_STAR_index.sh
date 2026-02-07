#!/bin/bash
#SBATCH --partition=corexd192
#SBATCH --job-name=STAR_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --output=/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/CriGri-PICRH-1.0_Ensemble/logs/star_index_%j.log

# 1. 设置路径
REF_ROOT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
# 注意修改以下参数
GENOME_NAME=""
REF_GENOME_DIR="${REF_ROOT_PATH}/${GENOME_NAME}"

GENOME_FA="${REF_GENOME_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
ANNOTATION_GTF="${REF_GENOME_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.115.gtf"
LOG_DIR="${REF_GENOME_DIR}/logs"

# 2. 参数设置
NCPUS=20
SJDB_OVERHANG=100
STAR_GENOME_RAM=190000000000 

set -e

# 激活环境
eval "$(mamba shell hook --shell bash)"
mamba activate /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4

cd ${REF_GENOME_DIR}

echo "================================================="
echo "STAR INDEX SCRIPT START: $(date)"
echo "================================================="

mkdir -p star.index

STAR \
    --runMode genomeGenerate \
    --genomeDir star.index \
    --genomeFastaFiles "${GENOME_FA}" \
    --sjdbGTFfile "${ANNOTATION_GTF}" \
    --sjdbOverhang "${SJDB_OVERHANG}" \
    --runThreadN "${NCPUS}" \
    --limitGenomeGenerateRAM "${STAR_GENOME_RAM}"

echo "### STAR INDEX COMPLETED: $(date) ###"