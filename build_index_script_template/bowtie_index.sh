#!/bin/bash
#SBATCH -p corexd192          # 分区/队列名称，请根据您的集群实际情况修改
#SBATCH --job-name=CriGri-PICRH-1.0_Ensemble_bowtie_index
#SBATCH -N 1                # 节点数
#SBATCH --ntasks-per-node=1 # 每个节点的任务数
#SBATCH -c 10               # 每个任务的CPU核心数
#SBATCH --mem=200G          # 内存申请
#SBATCH --time=200:00:00    # 运行时间
#SBATCH --output=/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/CriGri-PICRH-1.0_Ensemble/bowtie_index.log

echo "脚本启动时间: $(date)"
source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate ATAC_E4

# 确保参考基因组FASTA文件存在
FASTA_FILE="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/CriGri-PICRH-1.0_Ensemble/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
INDEX_BASENAME="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/CriGri-PICRH-1.0_Ensemble/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel"

if [ -f "$FASTA_FILE" ]; then
    echo "正在为 $FASTA_FILE 构建 Bowtie2 索引..."
    bowtie2-build "$FASTA_FILE" "$INDEX_BASENAME"
    echo "Bowtie2 索引构建完成。"
else
    echo "错误: FASTA 文件 $FASTA_FILE 不存在，无法构建索引。"
    exit 1
fi

echo "脚本结束时间: $(date)"
echo "作业成功完成。"