#!/bin/bash
## 用于表观基因组数据的目录初始化与数据转移

# 1. 基础路径配置
source_root="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects/sra_rawdata"
target_root="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects"

# ==========================================================
# 2. 细胞系和项目配置 (每次运行前修改此处！)
# ==========================================================
# 指定当前处理的细胞系名称
cell_line="Hela"
 
# 指定属于该细胞系的 Bioproject ID 列表
bioprojects=(
PRJEB40269
PRJEB59931
PRJEB79721)
# ==========================================================

# 3. 基础目录结构
# 1_result 包含子目录 0_fastq, 1_alignment, 2_tagalign, 3_peak_calling
base_structure="0_data,1_result/{0_fastq,1_alignment,2_tagalign,3_peak_calling},2_jobs,3_logs,tmp"

cd "$target_root" || exit

echo "Current Cell Line: $cell_line"

# 循环处理每一个项目
for bioproject in "${bioprojects[@]}"; do
  # 拼接目标文件夹名称，例如：PRJNA305986_HT1080
  dir_name="${bioproject}_${cell_line}"
  
  echo "------------------------------------------------"
  echo "Processing Project: $bioproject (dir_name: $dir_name)"

  # 4. 创建目标目录结构
  # eval 用于解析大括号扩展
  eval "mkdir -p $dir_name/{$base_structure}"
  
  # 5. 构建源路径和目标路径
  # 源：.../sra_rawdata/HT1080/PRJNA305986
  src_dir="${source_root}/${cell_line}/${bioproject}"
  # 目标：./PRJNA305986_HT1080/0_data/
  target_data_dir="${dir_name}/0_data"

  # 6. 执行数据转移
  if [ -d "$src_dir" ]; then
    # 检查源目录是否为空
    if [ "$(ls -A "$src_dir")" ]; then
      echo "Moving data from ${src_dir} to ${target_root}/${target_data_dir}"
      # 将源目录下所有内容（SRRxxx文件夹、id文件等）移至目标 0_data 下
      mv "$src_dir"/* "$target_data_dir/"
    else
      echo "Source directory $src_dir is empty. Skipping move."
    fi
  else
    echo "Warning: Source directory $src_dir does not exist."
  fi

  # 7. 清理空的源目录
  if [ -d "$src_dir" ]; then
    if rmdir "$src_dir" 2>/dev/null; then
      echo "Successfully removed empty source directory: $src_dir"
    else
      echo "Note: Source directory $src_dir was not removed (it may contain hidden files)."
    fi
  fi
done

# 可选：如果整个细胞系的源目录也空了，则删除
# rmdir "${source_root}/${cell_line}" 2>/dev/null

echo "------------------------------------------------"
echo "All tasks for $cell_line completed."
