#!/usr/bin/env bash
#SBATCH -J PRJNA593571_HEK293_upstream
#SBATCH -p corexd192
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/scRNA_projects/HEK293/PRJNA593571_HEK293/5_logs/PRJNA593571_HEK293_upstream.%j.log

set -eo pipefail

if [[ -f "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh" ]]; then
  source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
  conda activate "scRNA_env"
  echo "Successfully activate scRNA_env"
fi

cd "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/scRNA_script/GENtoolkit"

python /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/scRNA_script/GENtoolkit/GENToolkit_alter.py \
  -blt 10X \
  -rgf /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  -rgg /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensembl/Homo_sapiens.GRCh38.109.gtf \
  -ci /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensembl/cellranger_index \
  -bf /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensembl/hg38_Ensembl.109.for_scRNA.ref.bed \
  -hi /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensembl/hisat2_index/genome \
  -ri /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensembl/rsem.index/reference \
  -rd /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/scRNA_projects/HEK293/PRJNA593571_HEK293/1_raw \
  -rt sra \
  -ib index_exist \
  -sp /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/scRNA_env/bin \
  -da All_samples
