#!/usr/bin/env bash
#SBATCH -p corexd192          # 分区/队列名称，请根据您的集群实际情况修改
#SBATCH -J build_Cattle_E_UOA_Wagyu_1_cellranger_idx   # 作业名称
#SBATCH -N 1                # 节点数
#SBATCH --ntasks-per-node=1 # 每个节点的任务数
#SBATCH -c 10               # 每个任务的CPU核心数
#SBATCH --mem=200G          # 内存申请
#SBATCH --time=200:00:00    # 运行时间
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/Cattle_E_UOA_Wagyu_1_Ensemble/logs/cellranger_index.log

set -euo pipefail

# Activate mamba environment
if command -v mamba >/dev/null 2>&1; then
  eval "$(mamba shell hook --shell bash)"
  mamba activate scRNA_env
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

# Build Cell Ranger index for Cattle_E_UOA_Wagyu_1_Ensemble
REF_ROOT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
# 注意修改以下参数
GENOME_NAME="CriGri-PICRH-1.0_Ensemble"
GENOME_DIR="${REF_ROOT_PATH}/${GENOME_NAME}"
FASTA="${GENOME_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
GTF="${GENOME_DIR}/Homo_sapiens.GRCh38.115.gtf"


if [[ ! -f "${FASTA}" ]]; then
  echo "ERROR: FASTA not found: ${FASTA}" >&2
  exit 1
fi
if [[ ! -f "${GTF}" ]]; then
  echo "ERROR: GTF not found: ${GTF}" >&2
  exit 1
fi

# 1. 定义索引所在目录
INDEX_NAME="cellranger_index"

# 2. 切换到目标存放目录
cd "${GENOME_DIR}"

# 3. 运行命令（--genome 仅使用名称）
cellranger mkref --genome="${INDEX_NAME}" --fasta="${FASTA}" --genes="${GTF}"

echo "Cell Ranger genome built at: ${GENOME_DIR}/${INDEX_NAME}"