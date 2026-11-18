#!/bin/bash
#SBATCH -p corexd192
#SBATCH -J build_Dog_E_UUGSD_rsem_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 20
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH -o /path/to/user/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_rsem_index.log

set -e

echo "Job started at: $(date)"
echo "----------------------------------------------------"

REF_GENOME_DIR="/path/to/user/CellLine/resources/ref_genome/Dog_E_UUGSD"
GENOME_FNA_ORIG="/path/to/user/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
ANNOTATION_GTF="/path/to/user/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.115.gtf"
NCPUS=20

eval "$(mamba shell hook --shell bash)"
mamba activate /path/to/user/softwares/miniforge3/envs/RNAseq_E4
export PERL5LIB="/path/to/user/softwares/miniforge3/envs/RNAseq_E4/lib/perl5/5.32"

cd "${REF_GENOME_DIR}"
mkdir -p rsem.index

echo "### Building RSEM index... ###"
rsem-prepare-reference --gtf "${ANNOTATION_GTF}" -p "${NCPUS}" "${GENOME_FNA_ORIG}" rsem.index/reference

echo "### COMPLETED ###"
echo
echo "----------------------------------------------------"
echo "Job finished at: $(date)"
