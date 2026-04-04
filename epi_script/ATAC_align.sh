#!/usr/bin/env bash

# @File        :ATAC_align.sh 
# @Time        :2024/8/23 10:44
# @Author      :zhoubw (Updated for Picard compatibility)
# @Version     :2.1.0
# @Description :Main script for ATAC_pipeline. Added Read Group support for Picard.
# @Usage       :bash ATAC_align.sh <fastq_dir> <out_dir> <run_id> <biosample_id> <ncpus> <ramGB>

###parameter###
fastq_dir=$1
out_dir=$2
run_id=$3      # 对应 RG ID (如 ERR1951098)
sample_id=$4   # 对应 RG SM (如 BIOSAMPLE ID)
ncpus=$5       # 参数位置顺延
ramGB=$6       # 参数位置顺延
#############

# --- 加载项目配置 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

source "${CONFIG_PATH}"

# 检查所需变量
if [ -z "${BOWTIE2_INDEX}" ] || [ -z "${TMP_DIR}" ]; then
    echo "错误: 必要的环境变量 (BOWTIE2_INDEX, TMP_DIR) 未设置。" >&2
    exit 1
fi

# 定义 Read Group 字符串 (解决 Picard NullPointerException 的核心) 
RG_STR="--rg-id ${run_id} --rg SM:${sample_id} --rg LB:lib1 --rg PL:ILLUMINA"

### fastp ###
mkdir -p ${out_dir}/${run_id}/fastp
PE_raw=$(find "${fastq_dir}/${run_id}/reads/" -type f \( -name "${run_id}_1.fastq.gz" -o -name "${run_id}_2.fastq.gz" \))
SE_raw=$(find "${fastq_dir}/${run_id}/reads/" -type f -name "${run_id}.fastq.gz")

if [[ -n "$PE_raw" && -z "$SE_raw" ]]; then
    echo "...run fastp PE pipeline..."
    paired_end=0
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${run_id}/reads/${run_id}_1.fastq.gz \
    -I ${fastq_dir}/${run_id}/reads/${run_id}_2.fastq.gz \
    -o ${out_dir}/${run_id}/fastp/${run_id}.clean.R1.fastq.gz \
    -O ${out_dir}/${run_id}/fastp/${run_id}.clean.R2.fastq.gz \
    -j ${out_dir}/${run_id}/fastp/${run_id}_fastp.json \
    -h ${out_dir}/${run_id}/fastp/${run_id}_fastp.html \
    -R "${run_id}_fastp_report"
elif [[ -z "$PE_raw" && -n "$SE_raw" ]]; then
    echo "...run fastp SE pipeline..."
    paired_end=1
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${run_id}/reads/${run_id}.fastq.gz \
    -o ${out_dir}/${run_id}/fastp/${run_id}.clean.fastq.gz \
    -j ${out_dir}/${run_id}/fastp/${run_id}_fastp.json \
    -h ${out_dir}/${run_id}/fastp/${run_id}_fastp.html \
    -R "${run_id}_fastp_report"
else
    echo "ERROR: Raw sequencing file dose not exist!"
    exit 1
fi

### bowtie2 ###
mkdir -p ${out_dir}/${run_id}/bowtie2
cd ${out_dir}/${run_id}/bowtie2

if [ "$paired_end" -eq 0 ]; then
  echo "...run bowtie2 PE pipeline..."
  bowtie2 --mm --threads ${ncpus} -X2000 -q \
  ${RG_STR} \
  -x ${BOWTIE2_INDEX} \
  -1 ${out_dir}/${run_id}/fastp/${run_id}.clean.R1.fastq.gz \
  -2 ${out_dir}/${run_id}/fastp/${run_id}.clean.R2.fastq.gz | \
  samtools view -@ ${ncpus} -1 -S -b > ${run_id}.raw.bam

  samtools view -@ ${ncpus} -F 1804 -q 30 -f 2 -u ${run_id}.raw.bam | \
  samtools sort -@ ${ncpus} -n -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.tmp.bam
  samtools fixmate -@ ${ncpus} -r -O bam ${run_id}.tmp.bam ${run_id}.fixmate.bam
  samtools view -@ ${ncpus} -F 1804 -f 2 -u ${run_id}.fixmate.bam | \
  samtools sort -@ ${ncpus} -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.filt.bam

elif [ "$paired_end" -eq 1 ]; then
  echo "...run bowtie2 SE pipeline..."
  bowtie2 --mm --threads ${ncpus} \
  ${RG_STR} \
  -x ${BOWTIE2_INDEX} \
  -U ${out_dir}/${run_id}/fastp/${run_id}.clean.fastq.gz | \
  samtools view -@ ${ncpus} -1 -S -b > ${run_id}.raw.bam

  samtools view -@ ${ncpus} -F 1804 -q 30 -u ${run_id}.raw.bam | \
  samtools sort -@ ${ncpus} -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.filt.bam
fi

samtools index -@ ${ncpus} -b ${run_id}.filt.bam

# 使用 Picard 标记重复
picard -Xmx"${ramGB}G" MarkDuplicates \
INPUT=${run_id}.filt.bam \
OUTPUT=${run_id}.dupmark.bam \
METRICS_FILE=${run_id}.dup.qc \
VALIDATION_STRINGENCY=LENIENT \
USE_JDK_DEFLATER=TRUE \
USE_JDK_INFLATER=TRUE \
ASSUME_SORTED=TRUE \
REMOVE_DUPLICATES=FALSE

### 移除重复并生成最终 BAM ###
if [ "$paired_end" -eq 0 ]; then
  samtools view -@ ${ncpus} -F 1804 -f 2 -b ${run_id}.dupmark.bam > ${run_id}.nodup.bam
else
  samtools view -@ ${ncpus} -F 1804 -b ${run_id}.dupmark.bam > ${run_id}.nodup.bam
fi

samtools index -@ ${ncpus} -b ${run_id}.nodup.bam
samtools stats -@ ${ncpus} ${run_id}.nodup.bam > ${run_id}.nodup.bam.stats
samtools flagstat -@ ${ncpus} -O tsv ${run_id}.nodup.bam > ${run_id}.nodup.bam.flagstat

### TAG-Align 生成 ###
if [ "$paired_end" -eq 0 ]; then
  samtools sort -@ ${ncpus} -n ${run_id}.nodup.bam | \
  bedtools bamtobed -bedpe -mate1 -i ${run_id}.nodup.bam | \
  awk 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n", $1,$2,$3,$9,$4,$5,$6,$10}' | \
  gzip -nc > ${run_id}.nodup.tagAlign.gz
else
  bedtools bamtobed -i ${run_id}.nodup.bam | \
  awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
  gzip -nc > ${run_id}.nodup.tagAlign.gz
fi

# TN5 偏移校正
zcat -f ${run_id}.nodup.tagAlign.gz | awk 'BEGIN {OFS = "\t"} {
    if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5}
    if ($2 >= $3) {if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}}
    print $0
}' | gzip -nc > ${run_id}.tn5.tagAlign.gz