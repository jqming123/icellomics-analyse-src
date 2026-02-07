#!/bin/bash
#SBATCH -p corexd192          # 分区/队列名称，请根据您的集群实际情况修改
#SBATCH -J build_hg38_bwa_idx   # 作业名称
#SBATCH -N 1                # 节点数
#SBATCH --ntasks-per-node=1 # 每个节点的任务数
#SBATCH -c 10               # 每个任务的CPU核心数
#SBATCH --mem=200G          # 内存申请
#SBATCH --time=200:00:00    # 运行时间
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/logs/bwa_index.log

set -e pipefail

# --- 1. 定义变量和路径 ---
REF_ROOT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
# 注意修改以下参数
GENOME_NAME=""
REF_GENOME_DIR="${REF_ROOT_PATH}/${GENOME_NAME}"
GENOME_FASTA="${REF_GENOME_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
LOG_DIR="${REF_GENOME_DIR}/logs"
MAMBA_ENV_NAME="genome_env"

# 创建日志目录 (如果不存在)
mkdir -p ${LOG_DIR}




# --- 2. 记录脚本启动信息 ---
echo "=========================================================="
echo "脚本启动时间: $(date)"
echo "参考基因组路径: ${GENOME_FASTA}"
echo "=========================================================="



# --- 3. 设置软件环境 ---
echo "正在激活环境: ${MAMBA_ENV_NAME}..."
# 如果 mamba 没有直接挂载 profile.d，仍然可以通过 conda.sh 激活

eval "$(mamba shell hook --shell bash)"
mamba activate ${MAMBA_ENV_NAME}

echo "${MAMBA_ENV_NAME}环境激活成功。"

echo "bwa的路径是: $(which bwa)"

# --- 4. 前置检查 ---
echo "正在检查参考基因组文件是否存在..."
if [ ! -f "${GENOME_FASTA}" ]; then
    echo "错误: 参考基因组文件不存在: ${GENOME_FASTA}"
    exit 1
fi
echo "文件检查通过。"

# --- 5. 执行核心命令：建立索引 ---
echo "开始建立 BWA索引..."
cd ${REF_GENOME_DIR}

bwa index ${GENOME_FASTA}

if [ $? -eq 0 ]; then
    echo "BWA索引成功创建。"
else
    echo "错误: BWA索引建立失败。请检查错误日志。"
    exit 1
fi
echo

# --- 6. 验证输出并结束 ---
echo "验证生成的索引文件:"
ls -lh ${REF_GENOME_DIR} | grep "$(basename ${GENOME_FASTA})"


# --- 7. 创建序列字典 (.dict文件)，GATK需要它来处理染色体信息
echo "gatk的路径是: $(which gatk)"
cd ${REF_GENOME_DIR}
gatk CreateSequenceDictionary -R ${GENOME_FASTA}

echo "=========================================================="
echo "脚本结束时间: $(date)"
echo "作业成功完成。"
echo "=========================================================="

