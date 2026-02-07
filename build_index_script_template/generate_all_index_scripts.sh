#!/usr/bin/env bash
set -euo pipefail

# 填写一次参数，自动生成各个索引构建的 SLURM 脚本。
# 运行环境：Linux bash。

########################################
# 可配置参数
########################################
REF_ROOT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"  # 参考基因组根目录
GENOME_NAME="CriGri-PICRH-1.0_Ensemble"                                          # 参考基因组目录名
FASTA_FILENAME="Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"        # 基因组 FASTA 文件名
GTF_FILENAME="Cricetulus_griseus_picr.CriGri-PICRH-1.0.115.gtf"                  # 基因注释 GTF 文件名
# Kallisto 输出前缀
KALLISTO_PREFIX="CriGri-PICRH-1.0.115"                                          

# 环境名称 / 路径
ENV_ATAC="ATAC_E4"                                                   # bowtie2 脚本使用的环境
ENV_BWA="genome_env"                                                 # bwa 脚本使用的 mamba 环境
ENV_SCRNA="scRNA_env"                                                # hisat2 / cellranger 脚本使用的 mamba 环境
ENV_RNASEQ="RNAseq_E4"                                               # STAR / Kallisto / RSEM 使用

# 其他参数
SLURM_PARTITION="corexd192"
SLURM_MEM="200G"
SLURM_TIME="200:00:00"
HISAT2_CPUS=10                                         # hisat2-build 线程数
STAR_SJDB_OVERHANG=100
STAR_GENOME_RAM=190000000000
RSEM_CPUS=20
STAR_CPUS=20
KALLISTO_CPUS=20
BOWTIE_CPUS=10
BWA_CPUS=10
CELLRANGER_CPUS=10

########################################
# 派生变量（通常无需修改）
########################################
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_DIR="${REF_ROOT_PATH}/${GENOME_NAME}"
FASTA_PATH="${GENOME_DIR}/${FASTA_FILENAME}"
GTF_PATH="${GENOME_DIR}/${GTF_FILENAME}"
LOG_DIR="${GENOME_DIR}/logs"
mkdir -p "${LOG_DIR}"

LOG_BOWTIE="${LOG_DIR}/${GENOME_NAME}_bowtie_index.log"
LOG_BWA="${LOG_DIR}/${GENOME_NAME}_bwa_index.log"
LOG_CELLRANGER="${LOG_DIR}/${GENOME_NAME}_cellranger_index.log"
LOG_HISAT2="${LOG_DIR}/${GENOME_NAME}_hisat2_index.log"
LOG_KALLISTO="${LOG_DIR}/${GENOME_NAME}_kallisto_index_%j.log"
LOG_RSEM="${LOG_DIR}/${GENOME_NAME}_rsem_index.log"
LOG_STAR="${LOG_DIR}/${GENOME_NAME}_star_index_%j.log"
LOG_REF2BED="${LOG_DIR}/${GENOME_NAME}_ref2bed.log"

########################################
# 生成脚本
########################################

cat >"${SCRIPT_DIR}/build_bowtie_index.sh" <<EOF
#!/bin/bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH --job-name=${GENOME_NAME}_bowtie_index
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c ${BOWTIE_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=${LOG_BOWTIE}

echo "脚本启动时间: \\$(date)"
source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate ${ENV_ATAC}

REF_GENOME_DIR="${GENOME_DIR}"
FASTA_FILE="${FASTA_PATH}"
INDEX_BASENAME="${GENOME_DIR}/${FASTA_FILENAME%.*}"

if [ -f "\\${FASTA_FILE}" ]; then
    echo "正在为 \\${FASTA_FILE} 构建 Bowtie2 索引..."
    bowtie2-build "\\${FASTA_FILE}" "\\${INDEX_BASENAME}"
    echo "Bowtie2 索引构建完成。"
else
    echo "错误: FASTA 文件 \\${FASTA_FILE} 不存在，无法构建索引。"
    exit 1
fi

echo "脚本结束时间: \\$(date)"
echo "作业成功完成。"
EOF

cat >"${SCRIPT_DIR}/build_BWA_Index.sh" <<EOF
#!/bin/bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -J build_${GENOME_NAME}_bwa_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c ${BWA_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH -o ${LOG_BWA}

set -euo pipefail

REF_GENOME_DIR="${GENOME_DIR}"
GENOME_FASTA="${FASTA_PATH}"
LOG_DIR="${LOG_DIR}"
MAMBA_ENV_NAME="${ENV_BWA}"

mkdir -p \\${LOG_DIR}

echo "=========================================================="
echo "脚本启动时间: \\$(date)"
echo "参考基因组路径: \\${GENOME_FASTA}"
echo "=========================================================="

echo "正在激活环境: \\${MAMBA_ENV_NAME}..."
eval "\\$(mamba shell hook --shell bash)"
mamba activate \\${MAMBA_ENV_NAME}
echo "\\${MAMBA_ENV_NAME} 环境激活成功。"
echo "bwa 的路径是: \\$(which bwa)"

echo "正在检查参考基因组文件是否存在..."
if [ ! -f "\\${GENOME_FASTA}" ]; then
    echo "错误: 参考基因组文件不存在: \\${GENOME_FASTA}"
    exit 1
fi
echo "文件检查通过。"

echo "开始建立 BWA 索引..."
cd \\${REF_GENOME_DIR}
bwa index \\${GENOME_FASTA}
echo "BWA 索引成功创建。"

echo "验证生成的索引文件:"
ls -lh \\${REF_GENOME_DIR} | grep "\\$(basename \\${GENOME_FASTA})"

echo "gatk 的路径是: \\$(which gatk)"
cd \\${REF_GENOME_DIR}
gatk CreateSequenceDictionary -R \\${GENOME_FASTA}

echo "=========================================================="
echo "脚本结束时间: \\$(date)"
echo "作业成功完成。"
echo "=========================================================="
EOF

cat >"${SCRIPT_DIR}/build_cellranger_index.sh" <<EOF
#!/usr/bin/env bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -J build_${GENOME_NAME}_cellranger_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c ${CELLRANGER_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH -o ${LOG_CELLRANGER}

set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
  eval "\\$(mamba shell hook --shell bash)"
  mamba activate ${ENV_SCRNA}
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

GENOME_DIR="${GENOME_DIR}"
FASTA="${FASTA_PATH}"
GTF="${GTF_PATH}"
INDEX_NAME="cellranger_index"

if [[ ! -f "\\${FASTA}" ]]; then
  echo "ERROR: FASTA not found: \\${FASTA}" >&2
  exit 1
fi
if [[ ! -f "\\${GTF}" ]]; then
  echo "ERROR: GTF not found: \\${GTF}" >&2
  exit 1
fi

cd "\\${GENOME_DIR}"
cellranger mkref --genome="\\${INDEX_NAME}" --fasta="\\${FASTA}" --genes="\\${GTF}"
echo "Cell Ranger genome built at: \\${GENOME_DIR}/\\${INDEX_NAME}"
EOF

cat >"${SCRIPT_DIR}/build_hisat2_index.sh" <<EOF
#!/usr/bin/env bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -J build_${GENOME_NAME}_hisat2_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c ${HISAT2_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH -o ${LOG_HISAT2}

set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
  eval "\\$(mamba shell hook --shell bash)"
  mamba activate ${ENV_SCRNA}
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

GENOME_DIR="${GENOME_DIR}"
FASTA="${FASTA_PATH}"
OUT_PREFIX="\\${GENOME_DIR}/hisat2_index/genome"

if [[ ! -f "\\${FASTA}" ]]; then
  echo "ERROR: FASTA not found: \\${FASTA}" >&2
  exit 1
fi

mkdir -p "\\$(dirname "\\${OUT_PREFIX}")"
hisat2-build -p "${HISAT2_CPUS}" "\\${FASTA}" "\\${OUT_PREFIX}"
echo "HISAT2 index built at: \\${OUT_PREFIX}.*"
EOF

cat >"${SCRIPT_DIR}/build_Kallisto_index.sh" <<EOF
#!/bin/bash
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --job-name=${GENOME_NAME}_kallisto_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${KALLISTO_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=${LOG_KALLISTO}

REF_GENOME_DIR="${GENOME_DIR}"
GENOME_FA="${FASTA_PATH}"
ANNOTATION_GTF="${GTF_PATH}"
LOG_DIR="${LOG_DIR}"
PREFIX="${KALLISTO_PREFIX}"

set -e

eval "\\$(mamba shell hook --shell bash)"
mamba activate ${ENV_RNASEQ}

cd \\${REF_GENOME_DIR}
mkdir -p \\${LOG_DIR}

echo "================================================="
echo "KALLISTO PREP & INDEX START: \\$(date)"
echo "================================================="

TRANSCRIPT_FA="\\${PREFIX}.transcript.fa"
GENE_BED="\\${PREFIX}.gene.bed"
GENE_FA="\\${PREFIX}.gene.fa"

echo "--> Extracting transcript sequences with gffread..."
gffread "\\${ANNOTATION_GTF}" -g "\\${GENOME_FA}" -w "\\${TRANSCRIPT_FA}"

echo "--> Extracting gene sequences..."
awk 'BEGIN{OFS="\t"} $3=="gene" { if (match($9, /gene_id "([^"]+)"/, arr) && arr[1] != "") { print $1, $4-1, $5, arr[1], ".", $7 } }' "\\${ANNOTATION_GTF}" | sort -k1,1 -k2,2n -u > "\\${GENE_BED}"
bedtools getfasta -name -s -fi "\\${GENOME_FA}" -bed "\\${GENE_BED}" > "\\${GENE_FA}"

mkdir -p kallisto.index
echo "--> Indexing transcript sequences..."
kallisto index -i "kallisto.index/\\${PREFIX}.transcript.idx" "\\${TRANSCRIPT_FA}"

echo "--> Indexing gene sequences..."
kallisto index --make-unique -i "kallisto.index/\\${PREFIX}.gene.idx" "\\${GENE_FA}"
echo "### KALLISTO INDEX COMPLETED: \\$(date) ###"
EOF

cat >"${SCRIPT_DIR}/build_RSEM_index.sh" <<EOF
#!/bin/bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -J build_${GENOME_NAME}_rsem_idx
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c ${RSEM_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH -o ${LOG_RSEM}

set -e

echo "Job started at: \\$(date)"
echo "----------------------------------------------------"

REF_GENOME_DIR="${GENOME_DIR}"
GENOME_FNA_ORIG="${FASTA_PATH}"
ANNOTATION_GTF="${GTF_PATH}"
NCPUS=${RSEM_CPUS}

eval "\\$(mamba shell hook --shell bash)"
mamba activate ${ENV_RNASEQ}
export PERL5LIB="${ENV_RNASEQ}/lib/perl5/5.32"

cd "\\${REF_GENOME_DIR}"
mkdir -p rsem.index

echo "### Building RSEM index... ###"
rsem-prepare-reference --gtf "\\${ANNOTATION_GTF}" -p "\\${NCPUS}" "\\${GENOME_FNA_ORIG}" rsem.index/reference

echo "### COMPLETED ###"
echo
echo "----------------------------------------------------"
echo "Job finished at: \\$(date)"
EOF

cat >"${SCRIPT_DIR}/build_STAR_index.sh" <<EOF
#!/bin/bash
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --job-name=${GENOME_NAME}_star_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${STAR_CPUS}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=${LOG_STAR}

REF_GENOME_DIR="${GENOME_DIR}"
GENOME_FA="${FASTA_PATH}"
ANNOTATION_GTF="${GTF_PATH}"
LOG_DIR="${LOG_DIR}"
NCPUS=${STAR_CPUS}
SJDB_OVERHANG=${STAR_SJDB_OVERHANG}
STAR_GENOME_RAM=${STAR_GENOME_RAM}

set -e

eval "\\$(mamba shell hook --shell bash)"
mamba activate ${ENV_RNASEQ}

cd \\${REF_GENOME_DIR}
mkdir -p \\${LOG_DIR}
mkdir -p star.index

echo "================================================="
echo "STAR INDEX SCRIPT START: \\$(date)"
echo "================================================="

STAR \
    --runMode genomeGenerate \
    --genomeDir star.index \
    --genomeFastaFiles "\\${GENOME_FA}" \
    --sjdbGTFfile "\\${ANNOTATION_GTF}" \
    --sjdbOverhang "\\${SJDB_OVERHANG}" \
    --runThreadN "\\${NCPUS}" \
    --limitGenomeGenerateRAM "\\${STAR_GENOME_RAM}"

echo "### STAR INDEX COMPLETED: \\$(date) ###"
EOF

cat >"${SCRIPT_DIR}/ref2bed_for_scRNA.sh" <<EOF
#!/usr/bin/env bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -J ref2bed_${GENOME_NAME}_scRNA
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --time=${SLURM_TIME}
#SBATCH -o ${LOG_REF2BED}

set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
  eval "\$(mamba shell hook --shell bash)"
  mamba activate ${ENV_SCRNA}
else
  echo "ERROR: mamba not found in PATH." >&2
  exit 1
fi

GENOME_DIR="${GENOME_DIR}"
REF_GTF="${GTF_PATH}"
REF_BED_BASE="$(basename "${GTF_PATH%.*}")"
REF_BED="${GENOME_DIR}/${REF_BED_BASE}.for_scRNA.ref.bed"

if [[ ! -f "\${REF_GTF}" ]]; then
  echo "ERROR: GTF not found: \${REF_GTF}" >&2
  exit 1
fi

echo "================================================="
echo "REF2BED START: \$(date)"
echo "Input GTF: \${REF_GTF}"
echo "Output BED: \${REF_BED}"
echo "================================================="

gtf2bed < "\${REF_GTF}" > "\${REF_BED}"

echo "### REF2BED COMPLETED: \$(date) ###"
EOF

chmod +x "${SCRIPT_DIR}"/build_*_index.sh 2>/dev/null || true
chmod +x "${SCRIPT_DIR}/ref2bed_for_scRNA.sh" 2>/dev/null || true
echo "All index scripts have been regenerated under ${SCRIPT_DIR}."
