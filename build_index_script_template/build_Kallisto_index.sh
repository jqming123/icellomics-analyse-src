#!/bin/bash
#SBATCH --partition=corexd192
#SBATCH --job-name=Kallisto_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --output=/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/CriGri-PICRH-1.0_Ensemble/logs/kallisto_index_%j.log

# 1. 设置路径
REF_ROOT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
# 注意修改以下参数
GENOME_NAME=""
REF_GENOME_DIR="${REF_ROOT_PATH}/${GENOME_NAME}"
GENOME_FA="${REF_GENOME_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
ANNOTATION_GTF="${REF_GENOME_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.115.gtf"
LOG_DIR="${REF_GENOME_DIR}/logs"

# 设置前缀
PREFIX="CriGri-PICRH-1.0.115"

set -e

# 激活环境
eval "$(mamba shell hook --shell bash)"
mamba activate /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4

cd ${REF_GENOME_DIR}

echo "================================================="
echo "KALLISTO PREP & INDEX START: $(date)"
echo "================================================="

# 定义中间文件
TRANSCRIPT_FA="${PREFIX}.transcript.fa"
GENE_BED="${PREFIX}.gene.bed"
GENE_FA="${PREFIX}.gene.fa"

#--------------------------------------------------------------------------
# 步骤 1: 提取序列
#--------------------------------------------------------------------------
echo "--> Extracting transcript sequences with gffread..."
gffread "${ANNOTATION_GTF}" -g "${GENOME_FA}" -w "${TRANSCRIPT_FA}"

echo "--> Extracting gene sequences..."
# 【修改点 1】: 优化 awk 逻辑。
# 1. 增加 if 判断，只有成功匹配到 gene_id 且 ID 不为空时才打印，防止出现 "::coords" 格式的空 ID。
# 2. 增加 sort -u 对 BED 文件进行物理去重，防止 GTF 中同一基因存在多个重叠坐标导致序列重复。
awk 'BEGIN{OFS="\t"} $3=="gene" {
    if (match($9, /gene_id "([^"]+)"/, arr) && arr[1] != "") {
        print $1, $4-1, $5, arr[1], ".", $7
    }
}' "${ANNOTATION_GTF}" | sort -k1,1 -k2,2n -u > "${GENE_BED}"

# 【说明】: 使用 -name 参数会将 BED 第四列作为 FASTA 的 Header
bedtools getfasta -name -s -fi "${GENOME_FA}" -bed "${GENE_BED}" > "${GENE_FA}"

#--------------------------------------------------------------------------
# 步骤 2: 构建 Kallisto 索引
#--------------------------------------------------------------------------
mkdir -p kallisto.index

echo "--> Indexing transcript sequences..."
# 转录本索引通常 ID 较为规范，保持原样
kallisto index -i "kallisto.index/${PREFIX}.transcript.idx" "${TRANSCRIPT_FA}"

echo "--> Indexing gene sequences..."
# 【修改点 2】: 增加 --make-unique 参数。
# 即使 FASTA 中仍存在极个别重复命名的序列（例如同名基因位于不同染色体），
# 该参数会通过增加后缀（如 _1, _2）强制使其唯一，从而避免报错退出。
# 同时将输出后缀统一为 .idx 以适配你的定量脚本参数。
kallisto index --make-unique -i "kallisto.index/${PREFIX}.gene.idx" "${GENE_FA}"

echo "### KALLISTO INDEX COMPLETED: $(date) ###"