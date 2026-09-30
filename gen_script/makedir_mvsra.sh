#!/bin/bash
## 功能：指定一个细胞系名称和多个 BioProject ID，自动拼接目录名并迁移数据

# 1. 基础源路径
source_root="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/dna_sra"

# 2. 设置当前要处理的细胞系名称 (在此处修改)
cellline="Hela"

# 3. 设置 BioProject ID 列表 (在此处修改)
bioprojects=(
PRJNA358844
PRJNA529767
PRJNA752995
)

# 4. 基础目录结构 
base_structure="00_data/{dna_sra,raw_fastq,tmp,trimmed},01_results/{01_fastqc,02_bam,03_gvcf,04_genomicsdb,05_vcf_raw,06_vcf_filtered,logs,07_vcf_merged},02_jobs/{generated_mapping_jobs,genomicsdb_jobs,sra2fastq},03_logs/{01_mapping_gvcf,02_joint_calling,s2q_log_file}"

# 进入工作根目录
cd /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects || exit

# 5. 循环处理
for bp in "${bioprojects[@]}"; do
  # 拼接出目录名称 (例如: PRJEB10016_H9)
  dir_name="${bp}_${cellline}"

  echo "------------------------------------------------"
  echo "Processing: $dir_name"

  # A. 创建目录结构
  eval "mkdir -p $dir_name/{$base_structure}"
  
  # B. 构建源路径和目标路径
  # 源：/../dna_sra/H9/PRJEB10016
  src_dir="${source_root}/${cellline}/${bp}"
  # 目标：./PRJEB10016_H9/00_data/dna_sra
  target_dir="${dir_name}/00_data/dna_sra"

  # C. 执行剪切操作
  if [ -d "$src_dir" ]; then
    # 检查源目录是否为空
    if [ "$(ls -A "$src_dir")" ]; then
      echo "Moving files from $src_dir to $target_dir"
      mv "$src_dir"/* "$target_dir/"
      
      # D. 删除原目录
      if rmdir "$src_dir" 2>/dev/null; then
        echo "Successfully removed empty source directory: $src_dir"
      else
        echo "Note: Source directory $src_dir was not removed (possibly contains hidden files)."
      fi
    else
      echo "Source directory $src_dir is empty. Skipping move."
    fi
  else
    echo "Warning: Source directory $src_dir does not exist."
  fi
done

echo "------------------------------------------------"
echo "All tasks completed successfully."
