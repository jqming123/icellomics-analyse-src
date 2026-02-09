#!/bin/bash
#SBATCH --partition=corexd192
#SBATCH --job-name=Dog_E_UUGSD_star_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --output=/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_star_index_%j.log

REF_GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD"
GENOME_FA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
ANNOTATION_GTF="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.115.gtf"
LOG_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs"
NCPUS=20
SJDB_OVERHANG=100
STAR_GENOME_RAM=190000000000

set -e

eval "$(mamba shell hook --shell bash)"
mamba activate /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4

cd ${REF_GENOME_DIR}
mkdir -p ${LOG_DIR}
mkdir -p star.index

echo "================================================="
echo "STAR INDEX SCRIPT START: $(date)"
echo "================================================="

STAR     --runMode genomeGenerate     --genomeDir star.index     --genomeFastaFiles "${GENOME_FA}"     --sjdbGTFfile "${ANNOTATION_GTF}"     --sjdbOverhang "${SJDB_OVERHANG}"     --runThreadN "${NCPUS}"     --limitGenomeGenerateRAM "${STAR_GENOME_RAM}"

echo "### STAR INDEX COMPLETED: $(date) ###"
