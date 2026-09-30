#!/bin/bash -e

# ==============================================================================
# Generic Bulk RNA-seq Analysis Pipeline (fastp, STAR, kallisto, RSEM)
#
# 版本: 2.0
# 描述: 此脚本实现了一个标准的RNA-seq分析流程，已将所有物种相关的
#       参考文件路径参数化，使其能够适用于任何细胞系或物种。
#
# 使用方法:
# ./bulkRNA-seq_E4_v2_dUTP.sh -i <fastq_dir> -o <out_dir> -s <sample_id> \
#                     -c <cpus> -m <ram_gb> \
#                     -S <star_index_dir> \
#                     -G <kallisto_gene_idx> \
#                     -T <kallisto_transcript_idx> \
#                     -R <rsem_ref_prefix>
# ==============================================================================

### software vision###
#fastp:0.23.2
#kallisto:0.46.0
#STAR:2.73a
#rsem:1.3.3

# --- 默认参数 改成了1---
fwd_prob=1

# --- 帮助信息函数 ---
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Required Options:"
    echo "  -i <path>   输入FASTQ文件的根目录 (例如: /path/to/project_results)"
    echo "  -o <path>   输出结果的根目录 (例如: /path/to/project_results)"
    echo "  -s <string> 样本ID (例如: SRR123456)"
    echo "  -c <int>    使用的CPU核心数 (线程数)"
    echo "  -m <int>    使用的最大内存 (GB)"
    echo "  -S <path>   STAR索引文件所在的目录"
    echo "  -G <path>   Kallisto 基因 水平的索引文件 (*.idx)"
    echo "  -T <path>   Kallisto 转录本 水平的索引文件 (*.idx)"
    echo "  -R <path>   RSEM参考基因组的前缀 (例如: /path/to/rsem_ref/reference, rsem-prepare-reference的输出)"
    echo ""
    echo "Optional Options:"
    echo "  -h          显示此帮助信息并退出"
    exit 1
}

# --- 解析命令行参数 ---
while getopts "i:o:s:c:m:S:G:T:R:h" opt; do
    case ${opt} in
        i) fastq_dir=$OPTARG ;;
        o) out_dir=$OPTARG ;;
        s) sample_id=$OPTARG ;;
        c) ncpus=$OPTARG ;;
        m) ramGB=$OPTARG ;;
        S) star_index_dir=$OPTARG ;;
        G) kallisto_gene_idx=$OPTARG ;;
        T) kallisto_transcript_idx=$OPTARG ;;
        R) rsem_ref_prefix=$OPTARG ;;
        h) usage ;;
        \?) echo "无效的选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要一个参数." >&2; usage ;;
    esac
done

# --- 检查所有必需的参数是否已提供 ---
if [ -z "$fastq_dir" ] || [ -z "$out_dir" ] || [ -z "$sample_id" ] || \
   [ -z "$ncpus" ] || [ -z "$ramGB" ] || [ -z "$star_index_dir" ] || \
   [ -z "$kallisto_gene_idx" ] || [ -z "$kallisto_transcript_idx" ] || [ -z "$rsem_ref_prefix" ]; then
    echo "错误: 缺少一个或多个必需的参数。"
    usage
fi


# --- 打印运行参数 ---
echo "--- Pipeline Started: $(date) ---"
echo "Sample ID: $sample_id"
echo "FASTQ Dir: $fastq_dir"
echo "Output Dir: $out_dir"
echo "CPUs: $ncpus"
echo "Memory: ${ramGB}GB"
echo "STAR Index: $star_index_dir"
echo "Kallisto Gene Index: $kallisto_gene_idx"
echo "Kallisto Transcript Index: $kallisto_transcript_idx"
echo "RSEM Reference Prefix: $rsem_ref_prefix"
echo "-------------------------------------"


### private environment ###
# 注意: 这一部分可能需要根据你的服务器环境进行调整
eval "$(mamba shell hook --shell bash)"
mamba activate /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4
export PERL5LIB="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/RNAseq_E4/lib/perl5/5.32"
############

### 1. fastp: 质控和过滤 ###
echo "[Step 1/4] Running fastp for quality control..."
mkdir -p "${out_dir}/${sample_id}/fastp"

# 定义文件路径变量，方便后续引用和判断
raw_r1="${fastq_dir}/${sample_id}/reads/${sample_id}_1.fastq.gz"
raw_r2="${fastq_dir}/${sample_id}/reads/${sample_id}_2.fastq.gz"
raw_se="${fastq_dir}/${sample_id}/reads/${sample_id}.fastq.gz"

# 检查是否存在双端文件
if [[ -f "$raw_r1" && -f "$raw_r2" ]]; then
  echo "...Processing Paired-End (PE) data..."
  
  # 如果同时存在没有编号的 SE 文件，则将其删掉
  if [[ -f "$raw_se" ]]; then
    echo "WARNING: Found extra unnumbered file ${sample_id}.fastq.gz, deleting it to proceed with PE analysis."
    rm "$raw_se"
  fi

  paired_end=0
  fastp -g -q 5 -u 50 -n 5 -w "${ncpus}" \
    -i "$raw_r1" \
    -I "$raw_r2" \
    -o "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz" \
    -O "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz" \
    -j "${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json" \
    -h "${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html" \
    -R "${sample_id}_fastp_report"

# 如果不存在双端文件，但存在单端文件
elif [[ -f "$raw_se" ]]; then
  echo "...Processing Single-End (SE) data..."
  paired_end=1
  fastp -g -q 5 -u 50 -n 5 -w "${ncpus}" \
    -i "$raw_se" \
    -o "${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz" \
    -j "${out_dir}/${sample_id}/fastp/${sample_id}_fastp.json" \
    -h "${out_dir}/${sample_id}/fastp/${sample_id}_fastp.html" \
    -R "${sample_id}_fastp_report"

else
  echo "ERROR: Raw sequencing file does not exist or naming is incorrect!"
  echo "Expected PE: $raw_r1 and $raw_r2 OR SE: $raw_se"
  exit 1
fi

### 2. Kallisto & STAR: 定量和比对 ###
echo "[Step 2/4] Running Kallisto for quantification and STAR for alignment..."
mkdir -p "${out_dir}/${sample_id}/kallisto_gene"
mkdir -p "${out_dir}/${sample_id}/kallisto_transcript"
mkdir -p "${out_dir}/${sample_id}/star"

# 检查质控后的文件是否存在
PE_clean=$(find "${out_dir}/${sample_id}/fastp/" -type f \( -name "${sample_id}.clean.R1.fastq.gz" -o -name "${sample_id}.clean.R2.fastq.gz" \))
SE_clean=$(find "${out_dir}/${sample_id}/fastp/" -type f -name "${sample_id}.clean.fastq.gz")

if [[ -n "$PE_clean" && -z "$SE_clean" && "$paired_end" -eq 0 ]]; then
  # --- PE Pipeline ---
  echo "...Running PE analysis for Kallisto and STAR..."
  kallisto quant \
    -t "${ncpus}" --fusion --plaintext --fr-stranded \
    -i "${kallisto_gene_idx}" \
    -o "${out_dir}/${sample_id}/kallisto_gene" \
    "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz" \
    "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz"

  kallisto quant \
    -t "${ncpus}" --fusion --plaintext --fr-stranded \
    -i "${kallisto_transcript_idx}" \
    -o "${out_dir}/${sample_id}/kallisto_transcript" \
    "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz" \
    "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz"

  STAR \
    --genomeDir "${star_index_dir}" \
    --readFilesIn "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R1.fastq.gz" "${out_dir}/${sample_id}/fastp/${sample_id}.clean.R2.fastq.gz" \
    --outFileNamePrefix "${out_dir}/${sample_id}/star/${sample_id}.gsd." \
    --readFilesCommand zcat \
    --runThreadN "${ncpus}" \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM "${ramGB}000000000"

elif [[ -z "$PE_clean" && -n "$SE_clean" && "$paired_end" -eq 1 ]]; then
  # --- SE Pipeline ---
  echo "...Running SE analysis for Kallisto and STAR..."
  # 注意: 单端数据的 -l (平均片段长度) 和 -s (标准差) 是估计值，需要根据实际文库调整
  kallisto quant \
    -t "${ncpus}" --fusion --plaintext --single -l 200 -s 20 --fr-stranded \
    -i "${kallisto_gene_idx}" \
    -o "${out_dir}/${sample_id}/kallisto_gene" \
    "${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz"

  kallisto quant \
    -t "${ncpus}" --fusion --plaintext --single -l 200 -s 20 --fr-stranded \
    -i "${kallisto_transcript_idx}" \
    -o "${out_dir}/${sample_id}/kallisto_transcript" \
    "${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz"

  STAR \
    --genomeDir "${star_index_dir}" \
    --readFilesIn "${out_dir}/${sample_id}/fastp/${sample_id}.clean.fastq.gz" \
    --outFileNamePrefix "${out_dir}/${sample_id}/star/${sample_id}.gsd." \
    --readFilesCommand zcat \
    --runThreadN "${ncpus}" \
    --genomeLoad NoSharedMemory \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --sjdbScore 1 \
    --limitBAMsortRAM "${ramGB}000000000"

else
  echo "ERROR: Clean sequencing files does not exist or has a wrong name!"
  exit 1
fi

### 3. STAR: 生成 bedGraph (可选) ###
echo "[Step 3/4] Generating bedGraph files from BAM..."
bam_file="${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.sortedByCoord.out.bam"
if [ -f "$bam_file" ]; then
  STAR \
    --runMode inputAlignmentsFromBAM \
    --inputBAMfile "${bam_file}" \
    --outWigType bedGraph \
    --outWigStrand Stranded \
    --outFileNamePrefix "${out_dir}/${sample_id}/star/${sample_id}.gsd."
else
  echo "WARNING: Coordinate sorted BAM file not found, skipping bedGraph generation."
fi


### 4. RSEM: 基于比对的定量 ###
echo "[Step 4/4] Running RSEM for alignment-based quantification..."
mkdir -p "${out_dir}/${sample_id}/rsem"
transcriptome_bam_file="${out_dir}/${sample_id}/star/${sample_id}.gsd.Aligned.toTranscriptome.out.bam"

if [ -f "$transcriptome_bam_file" ]; then
  rsem-calculate-expression --version

  if [ "$paired_end" -eq 0 ]; then
    # --- PE RSEM ---
    rsem-calculate-expression \
      --bam --paired-end \
      --estimate-rspd \
      --seed 12345 \
      -p "${ncpus}" \
      --no-bam-output \
      --forward-prob ${fwd_prob} \
      "${transcriptome_bam_file}" \
      "${rsem_ref_prefix}" \
      "${out_dir}/${sample_id}/rsem/${sample_id}_rsem"

  elif [ "$paired_end" -eq 1 ]; then
    # --- SE RSEM ---
    rsem-calculate-expression \
      --bam \
      --estimate-rspd \
      --seed 12345 \
      -p "${ncpus}" \
      --no-bam-output \
      --forward-prob ${fwd_prob} \
      "${transcriptome_bam_file}" \
      "${rsem_ref_prefix}" \
      "${out_dir}/${sample_id}/rsem/${sample_id}_rsem"
  fi
else
  echo "ERROR: Transcriptome BAM file ($transcriptome_bam_file) does not exist! RSEM calculation failed."
  exit 1
fi

echo "--- Pipeline Finished: $(date) ---"
echo "All programs completed successfully for sample ${sample_id}!"