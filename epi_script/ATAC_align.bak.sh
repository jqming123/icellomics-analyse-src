#!/usr/bin/env bash

# @File       :ATAC_align.sh 
# @Time       :2024/8/23 10:44
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :ATAC_pipeline
# @Version    :2.0.0
# @Description:main script for ATAC_pipeline. Relies on environment variables from epi_config.sh
# @Usage      :bash ATAC_align.sh <fastq_dir> <out_dir> <run_id> <ncpus> <ramGB>

###software vision###
#fastp:0.23.2
#bowte2:2.5.3
#samtools:1.15.1
#picard:2.22.8
#bedtools:2.31.1
#macs2:2.2.4

# 注意: 此脚本依赖以下环境变量，请确保在运行前已通过 sourcing epi_config.sh 设置:
# - BOWTIE2_INDEX
# - TMP_DIR

###parameter###
fastq_dir=$1
out_dir=$2
run_id=$3
ncpus=$4
ramGB=$5
#############

# --- 加载项目配置 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

source "${CONFIG_PATH}"

# 检查所需变量是否已设置
if [ -z "${BOWTIE2_INDEX}" ] || [ -z "${TMP_DIR}" ]; then
    echo "错误: 必要的环境变量 (BOWTIE2_INDEX, TMP_DIR) 未设置。" >&2
    echo "请确保在调用此脚本前 sourcing 了正确的配置文件 (epi_config.sh)。" >&2
    exit 1
fi



###fastp###
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

    cd ${out_dir}/${run_id}/fastp
    json_file="${run_id}_fastp.json"
    read1_length=$(jq -r '.summary.after_filtering.read1_mean_length' "$json_file")
    read2_length=$(jq -r '.summary.after_filtering.read2_mean_length' "$json_file")
    echo "read1_mean_length: $read1_length"
    echo "read2_mean_length: $read2_length"

elif [[ -z "$PE_raw" && -n "$SE_raw" ]]; then
    echo "...run fastp SE pipeline..."
    paired_end=1
    fastp -g -q 5 -u 50 -n 5 \
    -i ${fastq_dir}/${run_id}/reads/${run_id}.fastq.gz \
    -o ${out_dir}/${run_id}/fastp/${run_id}.clean.fastq.gz \
    -j ${out_dir}/${run_id}/fastp/${run_id}_fastp.json \
    -h ${out_dir}/${run_id}/fastp/${run_id}_fastp.html \
    -R "${run_id}_fastp_report"

    cd ${out_dir}/${run_id}/fastp
    json_file="${run_id}_fastp.json"
    read1_length=$(jq -r '.summary.after_filtering.read1_mean_length' "$json_file")
    echo "read1_mean_length: $read1_length"
else
    echo "ERROR: Raw sequencing file dose not exist or have a wrong name!"
    exit 1
fi


###bowtie2###
mkdir -p ${out_dir}/${run_id}/bowtie2
cd ${out_dir}/${run_id}/bowtie2

if [ "$paired_end" -eq 0 ]; then
  echo "...run bowtie2 PE pipeline..."
  bowtie2 --mm --threads ${ncpus} -X2000 -q \
  -x ${BOWTIE2_INDEX} \
  -1 ${out_dir}/${run_id}/fastp/${run_id}.clean.R1.fastq.gz \
  -2 ${out_dir}/${run_id}/fastp/${run_id}.clean.R2.fastq.gz | \
  samtools view -@ ${ncpus} -1 -S -b > ${run_id}.raw.bam

  # 使用 ramGB 参数设置 samtools sort 的内存
  samtools view -@ ${ncpus} -F 1804 -q 30 -f 2 -u ${run_id}.raw.bam | \
  samtools sort -@ ${ncpus} -n -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.tmp.bam
  samtools fixmate -@ ${ncpus} -r -O bam ${run_id}.tmp.bam ${run_id}.fixmate.bam

  samtools view -@ ${ncpus} -F 1804 -f 2 -u ${run_id}.fixmate.bam |
  samtools sort -@ ${ncpus} -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.filt.bam

elif [ "$paired_end" -eq 1 ]; then
  echo "...run bowtie2 SE pipeline..."
  bowtie2 --mm --threads ${ncpus} \
  -x ${BOWTIE2_INDEX} \
  -U ${out_dir}/${run_id}/fastp/${run_id}.clean.fastq.gz | \
  samtools view -@ ${ncpus} -1 -S -b > ${run_id}.raw.bam

  # 使用 ramGB 参数设置 samtools sort 的内存
  samtools view -@ ${ncpus} -F 1804 -q 30 -u ${run_id}.raw.bam | \
  samtools sort -@ ${ncpus} -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.filt.bam
else
    echo "ERROR: ${run_id}.clean.R1|R2.fastq.gz does not exist!"
    exit 1
fi

samtools index -@ ${ncpus} -b ${run_id}.filt.bam

# 使用 ramGB 参数设置 Picard 的最大堆内存
picard -Xmx"${ramGB}G" MarkDuplicates \
INPUT=${run_id}.filt.bam \
OUTPUT=${run_id}.dupmark.bam \
METRICS_FILE=${run_id}.dup.qc \
VALIDATION_STRINGENCY=LENIENT \
USE_JDK_DEFLATER=TRUE \
USE_JDK_INFLATER=TRUE \
ASSUME_SORTED=TRUE \
REMOVE_DUPLICATES=FALSE

if [ "$paired_end" -eq 0 ]; then
  echo "...run remove duplication PE pipeline..."
  samtools view -@ ${ncpus} -F 1804 -f 2 -b ${run_id}.dupmark.bam > ${run_id}.nodup.bam
elif [ "$paired_end" -eq 1 ]; then
  echo "...run remove duplication SE pipeline..."
  samtools view -@ ${ncpus} -F 1804 -b ${run_id}.dupmark.bam > ${run_id}.nodup.bam
else
  echo "ERROR: ${run_id}.dupmark.bam does not exist!"
  exit 1
fi

samtools index -@ ${ncpus} -b ${run_id}.nodup.bam
samtools stats -@ ${ncpus} ${run_id}.nodup.bam > ${run_id}.nodup.bam.stats
samtools flagstat -@ ${ncpus} -O tsv ${run_id}.nodup.bam > ${run_id}.nodup.bam.flagstat

###TAG-Align###
if [ "$paired_end" -eq 0 ]; then
  echo "...run bam2ta PE pipeline..."
  bedtools bamtobed -bedpe -mate1 -i ${run_id}.nodup.bam | \
  awk 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n", $1,$2,$3,$9,$4,$5,$6,$10}' | \
  gzip -nc > ${run_id}.nodup.tagAlign.gz

elif [ "$paired_end" -eq 1 ]; then
  echo "...run bam2ta SE pipeline..."
  bedtools bamtobed -i ${run_id}.nodup.bam | \
  awk 'BEGIN{{OFS="\\t"}}{{$4="N";$5="1000";print $0}}\' | \
  gzip -nc > ${run_id}.nodup.tagAlign.gz
else
    echo "ERROR: ${run_id}.nodup.bam does not exist!"
    exit 1
fi

zcat -f ${run_id}.nodup.tagAlign.gz | awk 'BEGIN {OFS = "\t"} {
    if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5}
    if ($2 >= $3) {if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}}
    print $0
}' | gzip -nc > ${run_id}.tn5.tagAlign.gz