#!/bin/bash
#SBATCH -p corexd192
#SBATCH -J build_Dog_E_UUGSD_bwa_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 10
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs/Dog_E_UUGSD_bwa_index.log

set -euo pipefail

REF_GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD"
GENOME_FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
LOG_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Dog_E_UUGSD/logs"
MAMBA_ENV_NAME="genome_env"

mkdir -p ${LOG_DIR}

echo "=========================================================="
echo "脚本启动时间: $(date)"
echo "参考基因组路径: ${GENOME_FASTA}"
echo "=========================================================="

echo "正在激活环境: ${MAMBA_ENV_NAME}..."
eval "$(mamba shell hook --shell bash)"
mamba activate ${MAMBA_ENV_NAME}
echo "${MAMBA_ENV_NAME} 环境激活成功。"
echo "bwa 的路径是: $(which bwa)"

echo "正在检查参考基因组文件是否存在..."
if [ ! -f "${GENOME_FASTA}" ]; then
    echo "错误: 参考基因组文件不存在: ${GENOME_FASTA}"
    exit 1
fi
echo "文件检查通过。"

echo "开始建立 BWA 索引..."
cd ${REF_GENOME_DIR}
bwa index ${GENOME_FASTA}
echo "BWA 索引成功创建。"

echo "验证生成的索引文件:"
ls -lh ${REF_GENOME_DIR} | grep "$(basename ${GENOME_FASTA})"

echo "gatk 的路径是: $(which gatk)"
cd ${REF_GENOME_DIR}
gatk CreateSequenceDictionary -R ${GENOME_FASTA}

echo "=========================================================="
echo "脚本结束时间: $(date)"
echo "作业成功完成。"
echo "=========================================================="
