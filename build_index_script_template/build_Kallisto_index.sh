#!/bin/bash
#SBATCH --partition=corexd192
#SBATCH --job-name=Dog_E_UUGSD_kallisto_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --output=/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_kallisto_index_%j.log

REF_GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD"
GENOME_FA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
ANNOTATION_GTF="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.115.gtf"
LOG_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs"
PREFIX="UU_Cfam_GSD_1.0.115"

set -e

eval "$(mamba shell hook --shell bash)"
mamba activate /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4

cd ${REF_GENOME_DIR}
mkdir -p ${LOG_DIR}

echo "================================================="
echo "KALLISTO PREP & INDEX START: $(date)"
echo "================================================="

TRANSCRIPT_FA="${PREFIX}.transcript.fa"
GENE_BED="${PREFIX}.gene.bed"
GENE_FA="${PREFIX}.gene.fa"

echo "--> Extracting transcript sequences with gffread..."
gffread "${ANNOTATION_GTF}" -g "${GENOME_FA}" -w "${TRANSCRIPT_FA}"

echo "--> Extracting gene sequences..."
awk 'BEGIN{OFS="\t"} =="gene" { if (match(, /gene_id "([^"]+)"/, arr) && arr[1] != "") { print , -1, , arr[1], ".",  } }' "${ANNOTATION_GTF}" | sort -k1,1 -k2,2n -u > "${GENE_BED}"
bedtools getfasta -name -s -fi "${GENOME_FA}" -bed "${GENE_BED}" > "${GENE_FA}"

mkdir -p kallisto.index
echo "--> Indexing transcript sequences..."
kallisto index -i "kallisto.index/${PREFIX}.transcript.idx" "${TRANSCRIPT_FA}"

echo "--> Indexing gene sequences..."
kallisto index --make-unique -i "kallisto.index/${PREFIX}.gene.idx" "${GENE_FA}"
echo "### KALLISTO INDEX COMPLETED: $(date) ###"
