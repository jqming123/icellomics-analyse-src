#!/bin/bash
#SBATCH -p corexd192          # 分区/队列名称，请根据您的集群实际情况修改
#SBATCH -J build_hg38_resm_idx   # 作业名称
#SBATCH -N 1                # 节点数
#SBATCH --ntasks-per-node=1 # 每个节点的任务数
#SBATCH -c 20               # 每个任务的CPU核心数
#SBATCH --mem=200G          # 内存申请
#SBATCH --time=200:00:00    # 运行时间
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/logs/rsem_index.log # 标准输出


set -e

# 打印作业开始时间
echo "Job started at: $(date)"
echo "----------------------------------------------------"

# 1. 设置参考基因组工作目录
REF_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble"

# 2. 设置原始参考文件
GENOME_FNA_ORIG="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
ANNOTATION_GTF="${REF_DIR}/Homo_sapiens.GRCh38.115.gtf" 

# 3. 设置线程数
NCPUS=20


###private environment###
eval "$(mamba shell hook --shell bash)"
mamba activate /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4
export PERL5LIB="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4/lib/perl5/5.32"


cd "${REF_DIR}"

echo "### Building RSEM index... ###"
mkdir -p rsem.index

rsem-prepare-reference --gtf "${ANNOTATION_GTF}" -p "${NCPUS}" "${GENOME_FNA_ORIG}" rsem.index/reference

echo "### COMPLETED ###"
echo

# 打印作业结束时间
echo "----------------------------------------------------"
echo "Job finished at: $(date)"