#!/bin/bash
#SBATCH -p corexd192
#SBATCH --job-name=Dog_E_UUGSD_bowtie_index
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 10
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --output=/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_bowtie_index.log

echo "脚本启动时间: $(date)"
source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate ATAC_E4

REF_GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD"
FASTA_FILE="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
INDEX_BASENAME="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel"

if [ -f "${FASTA_FILE}" ]; then
    echo "正在为 ${FASTA_FILE} 构建 Bowtie2 索引..."
    bowtie2-build "${FASTA_FILE}" "${INDEX_BASENAME}"
    echo "Bowtie2 索引构建完成。"
else
    echo "错误: FASTA 文件 ${FASTA_FILE} 不存在，无法构建索引。"
    exit 1
fi

echo "脚本结束时间: $(date)"
echo "作业成功完成。"
