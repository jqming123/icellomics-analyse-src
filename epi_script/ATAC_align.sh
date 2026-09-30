#!/usr/bin/env bash

set -euo pipefail

# @File        :ATAC_align.sh 
# @Time        :2024/8/23 10:44
# @Author      :zhoubw (Updated for Picard compatibility)
# @Version     :2.2.0
# @Description :Run-level preprocessing/alignment. Replicate-aware jobs use alignment_only=true.
# @Usage       :bash ATAC_align.sh <fastq_dir> <out_dir> <run_id> <biosample_id> <ncpus> <ramGB> [alignment_only]

if [ "$#" -lt 6 ] || [ "$#" -gt 7 ]; then
    echo "Usage: bash ATAC_align.sh <fastq_dir> <out_dir> <run_id> <biosample_id> <ncpus> <ramGB> [alignment_only]" >&2
    exit 1
fi

###parameter###
fastq_dir=$1
out_dir=$2
run_id=$3      # 对应 RG ID (如 ERR1951098)
sample_id=$4   # 对应 RG SM (如 BIOSAMPLE ID)
ncpus=$5       # 参数位置顺延
ramGB=$6       # 参数位置顺延
alignment_only=${7:-false}
#############

# --- 加载项目配置 ---
CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "错误: 配置文件未找到于 ${CONFIG_PATH}"
    exit 1
fi

source "${CONFIG_PATH}"

# 检查所需变量
if [ -z "${BOWTIE2_INDEX:-}" ] || [ -z "${TMP_DIR:-}" ]; then
    echo "错误: 必要的环境变量 (BOWTIE2_INDEX, TMP_DIR) 未设置。" >&2
    exit 1
fi

# 定义 Read Group 字符串 (解决 Picard NullPointerException 的核心) 
RG_STR="--rg-id ${run_id} --rg SM:${sample_id} --rg LB:lib1 --rg PL:ILLUMINA"

### fastp ###
mkdir -p "${out_dir}/${run_id}/fastp"
READS_DIR="${fastq_dir}/${run_id}/reads"
RAW_R1="${READS_DIR}/${run_id}_1.fastq.gz"
RAW_R2="${READS_DIR}/${run_id}_2.fastq.gz"
RAW_SE="${READS_DIR}/${run_id}.fastq.gz"

if [[ -f "${RAW_R1}" && -f "${RAW_R2}" ]]; then
    echo "...run fastp PE pipeline (found complete PE reads$( [[ -f "${RAW_SE}" ]] && echo ', ignoring co-existing SE file'))..."
    paired_end=0
    fastp -g -q 5 -u 50 -n 5 \
    -i "${RAW_R1}" \
    -I "${RAW_R2}" \
    -o "${out_dir}/${run_id}/fastp/${run_id}.clean.R1.fastq.gz" \
    -O "${out_dir}/${run_id}/fastp/${run_id}.clean.R2.fastq.gz" \
    -j "${out_dir}/${run_id}/fastp/${run_id}_fastp.json" \
    -h "${out_dir}/${run_id}/fastp/${run_id}_fastp.html" \
    -R "${run_id}_fastp_report"
elif [[ -f "${RAW_R1}" || -f "${RAW_R2}" ]]; then
    echo "ERROR: Incomplete PE FASTQ files. Both files must exist: ${RAW_R1} and ${RAW_R2}" >&2
    exit 1
elif [[ -f "${RAW_SE}" ]]; then
    echo "...run fastp SE pipeline..."
    paired_end=1
    fastp -g -q 5 -u 50 -n 5 \
    -i "${RAW_SE}" \
    -o "${out_dir}/${run_id}/fastp/${run_id}.clean.fastq.gz" \
    -j "${out_dir}/${run_id}/fastp/${run_id}_fastp.json" \
    -h "${out_dir}/${run_id}/fastp/${run_id}_fastp.html" \
    -R "${run_id}_fastp_report"
else
    echo "ERROR: Raw sequencing file does not exist! Checked: ${READS_DIR}"
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
  -2 ${out_dir}/${run_id}/fastp/${run_id}.clean.R2.fastq.gz \
  2> ${run_id}.bowtie2.log | \
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
  -U ${out_dir}/${run_id}/fastp/${run_id}.clean.fastq.gz \
  2> ${run_id}.bowtie2.log | \
  samtools view -@ ${ncpus} -1 -S -b > ${run_id}.raw.bam

  samtools view -@ ${ncpus} -F 1804 -q 30 -u ${run_id}.raw.bam | \
  samtools sort -@ ${ncpus} -m "${ramGB}G" -O bam -T ${TMP_DIR} -o ${run_id}.filt.bam
fi

samtools index -@ ${ncpus} -b ${run_id}.filt.bam

# Raw mitochondrial fraction is measured before duplicate and mitochondrial
# removal. It is diagnostic and intentionally uses mapped read records.
samtools flagstat -@ ${ncpus} -O tsv "${run_id}.raw.bam" > "${run_id}.raw.bam.flagstat"
RAW_MAPPED=$(samtools view -@ ${ncpus} -c -F 4 "${run_id}.raw.bam")
RAW_MITO=$(samtools view -@ ${ncpus} -F 4 "${run_id}.raw.bam" | \
  awk -v mito_contigs="${MITO_CONTIGS:-}" 'BEGIN {
    n=split(mito_contigs, contigs, /[[:space:],]+/)
    for(i=1;i<=n;i++) if(contigs[i] != "") is_mito[contigs[i]]=1
  }
  ($3 in is_mito) {count++}
  END {print count+0}')
printf "mapped_reads\tmitochondrial_reads\n%s\t%s\n" "${RAW_MAPPED}" "${RAW_MITO}" > "${run_id}.raw.mito.qc.tsv"

# Record filtering evidence while the raw and filtered run-level BAMs are both
# present.  This sidecar is deliberately written before duplicate marking and
# mitochondrial filtering, so it documents the mapping-quality decision rather
# than attempting to infer it from the final nodup BAM alone.
if [ "$paired_end" -eq 0 ]; then
  FILTER_LAYOUT="PE"
  REQUIRED_FLAG=2
  POST_MAPQ_BAM="${run_id}.tmp.bam"
else
  FILTER_LAYOUT="SE"
  REQUIRED_FLAG=null
  POST_MAPQ_BAM="${run_id}.filt.bam"
fi
if [ "${alignment_only}" = "true" ]; then
  DUPLICATE_STAGE="deferred_to_ATAC_finalize_replicate_after_technical_run_merge"
else
  DUPLICATE_STAGE="pending_run_level_MarkDuplicates"
fi
RAW_ALIGNMENT_RECORDS=$(samtools view -@ "${ncpus}" -c "${run_id}.raw.bam")
POST_MAPQ_ALIGNMENT_RECORDS=$(samtools view -@ "${ncpus}" -c "${POST_MAPQ_BAM}")
POST_FINAL_FILTER_ALIGNMENT_RECORDS=$(samtools view -@ "${ncpus}" -c "${run_id}.filt.bam")
FILTER_PROVENANCE="${run_id}.mapping_filter.provenance.json"
{
  printf '{\n'
  printf '  "schema_version": 1,\n'
  printf '  "run_id": "%s",\n' "${run_id}"
  printf '  "layout": "%s",\n' "${FILTER_LAYOUT}"
  printf '  "mapq_filter_threshold": 30,\n'
  printf '  "mapq_filter_operator": ">=",\n'
  printf '  "samtools_view_exclude_flag": 1804,\n'
  printf '  "samtools_view_required_flag": %s,\n' "${REQUIRED_FLAG}"
  printf '  "raw_bam_alignment_record_count": %s,\n' "${RAW_ALIGNMENT_RECORDS}"
  printf '  "post_mapq_filter_alignment_record_count": %s,\n' "${POST_MAPQ_ALIGNMENT_RECORDS}"
  printf '  "post_final_filter_alignment_record_count": %s,\n' "${POST_FINAL_FILTER_ALIGNMENT_RECORDS}"
  printf '  "post_mapq_filter_bam": "%s",\n' "${POST_MAPQ_BAM}"
  printf '  "final_filter_bam": "%s",\n' "${run_id}.filt.bam"
  printf '  "duplicate_handling_stage": "%s",\n' "${DUPLICATE_STAGE}"
  printf '  "mitochondrial_filter_stage": "pending_after_duplicate_handling"\n'
  printf '}\n'
} > "${FILTER_PROVENANCE}"
echo "Wrote mapping-filter provenance: ${FILTER_PROVENANCE}"

if [ "${alignment_only}" = "true" ]; then
  echo "Run-level alignment complete; duplicate removal will occur after technical-run merging."
  exit 0
fi

# 使用 Picard 标记重复
picard -Xmx"${ramGB}G" MarkDuplicates \
INPUT=${run_id}.filt.bam \
OUTPUT=${run_id}.dupmark.bam \
METRICS_FILE=${run_id}.dup.qc \
VALIDATION_STRINGENCY=LENIENT \
USE_JDK_DEFLATER=TRUE \
USE_JDK_INFLATER=TRUE \
VERBOSITY=ERROR \
ASSUME_SORTED=TRUE \
REMOVE_DUPLICATES=FALSE

### 移除重复并生成最终 BAM ###
NODUP_WITH_MITO_BAM="${run_id}.nodup.with_mito.bam"
FINAL_NODUP_BAM="${run_id}.nodup.bam"

if [ "$paired_end" -eq 0 ]; then
  samtools view -@ ${ncpus} -F 1804 -f 2 -b ${run_id}.dupmark.bam > "${NODUP_WITH_MITO_BAM}"
else
  samtools view -@ ${ncpus} -F 1804 -b ${run_id}.dupmark.bam > "${NODUP_WITH_MITO_BAM}"
fi

if [ -n "${MITO_CONTIGS:-}" ]; then
  echo "...remove mitochondrial reads from final nodup BAM: ${MITO_CONTIGS}..."
  samtools index -@ ${ncpus} -b "${NODUP_WITH_MITO_BAM}"
  samtools idxstats "${NODUP_WITH_MITO_BAM}" > "${NODUP_WITH_MITO_BAM}.idxstats"
  samtools view -@ ${ncpus} -h "${NODUP_WITH_MITO_BAM}" | \
  awk -v mito_contigs="${MITO_CONTIGS}" 'BEGIN {
      n = split(mito_contigs, contigs, /[[:space:],]+/)
      for (i = 1; i <= n; i++) {
          if (contigs[i] != "") {
              is_mito[contigs[i]] = 1
          }
      }
  }
  /^@/ {print; next}
  {
      mate = $7
      if (!(($3 in is_mito) || (mate in is_mito))) {
          print
      }
  }' | \
  samtools view -@ ${ncpus} -b -o "${FINAL_NODUP_BAM}" -
  rm -f "${NODUP_WITH_MITO_BAM}" "${NODUP_WITH_MITO_BAM}.bai"
else
  mv "${NODUP_WITH_MITO_BAM}" "${FINAL_NODUP_BAM}"
fi

samtools index -@ ${ncpus} -b "${FINAL_NODUP_BAM}"
samtools stats -@ ${ncpus} "${FINAL_NODUP_BAM}" > "${FINAL_NODUP_BAM}.stats"
samtools flagstat -@ ${ncpus} -O tsv "${FINAL_NODUP_BAM}" > "${FINAL_NODUP_BAM}.flagstat"
samtools idxstats "${FINAL_NODUP_BAM}" > "${FINAL_NODUP_BAM}.idxstats"

### TAG-Align 生成 ###
if [ "$paired_end" -eq 0 ]; then
  samtools sort -@ ${ncpus} -n "${FINAL_NODUP_BAM}" -o ${run_id}.nodup.namesort.bam
  bedtools bamtobed -bedpe -mate1 -i ${run_id}.nodup.namesort.bam | \
  awk 'BEGIN{OFS="\t"} {printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n", $1,$2,$3,$9,$4,$5,$6,$10}' | \
  gzip -nc > ${run_id}.nodup.tagAlign.gz
else
  bedtools bamtobed -i "${FINAL_NODUP_BAM}" | \
  awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | \
  gzip -nc > ${run_id}.nodup.tagAlign.gz
fi

# TN5 偏移校正
zcat -f ${run_id}.nodup.tagAlign.gz | awk 'BEGIN {OFS = "\t"} {
    if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5}
    if ($2 >= $3) {if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1}}
    print $0
}' | gzip -nc > ${run_id}.tn5.tagAlign.gz
