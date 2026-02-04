#!/usr/bin/env bash
#SBATCH -p corexd192          # 分区/队列名称，请根据您的集群实际情况修改
#SBATCH -J build_hg38_hisat2_idx   # 作业名称
#SBATCH -N 1                # 节点数
#SBATCH --ntasks-per-node=1 # 每个节点的任务数
#SBATCH -c 10               # 每个任务的CPU核心数
#SBATCH --mem=200G          # 内存申请
#SBATCH --time=200:00:00    # 运行时间
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/logs/hisat2_index.log

set -euo pipefail

# Activate mamba environment
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
  mamba activate scRNA_env
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

# Build HISAT2 index for hg38_Ensemble
GENOME_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble"
FASTA="${GENOME_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
OUT_PREFIX="${GENOME_DIR}/hisat2_index/genome"

if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: FASTA not found: ${FASTA}" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUT_PREFIX}")"

hisat2-build -p "${SLURM_CPUS_PER_TASK:-1}" "${FASTA}" "${OUT_PREFIX}"

echo "HISAT2 index built at: ${OUT_PREFIX}.*"
