#!/bin/bash

# 定义基础路径
RNA_COUNT_RESULT_BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"
EQ_JOBS_BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/EQ_jobs"

# 检查是否提供了正确的参数数量
if [ "$#" -ne 1 ]; then
  echo "用法: $0 <项目目录名>"
  echo "示例: $0 PRJEB41085_HEK293"
  exit 1
fi

prj_name=$1

echo "--- 1. 获取待检查目录列表 ---"
job_scripts_path="$EQ_JOBS_BASE_DIR/$prj_name"

if [ ! -d "$job_scripts_path" ]; then
  echo "错误: 项目脚本目录 '$job_scripts_path' 不存在。"
  exit 1
fi

# 查找所有 .sh 脚本文件，并按名称排序
# 将结果存储为一个字符串，每行一个文件路径
all_sh_files_str=$(find "$job_scripts_path" -maxdepth 1 -type f -name "*.sh" | sort)

if [ -z "$all_sh_files_str" ]; then
  echo "错误: 在 '$job_scripts_path' 中未找到任何 .sh 脚本文件。"
  exit 1
fi

# 从所有 .sh 脚本文件名中提取目录名列表
extracted_dirs=()
while IFS= read -r file_path; do
  dir_name=$(basename "$file_path" .sh)
  extracted_dirs+=("$dir_name")
done <<< "$all_sh_files_str"

# 获取第一个和最后一个脚本对应的目录名，用于显示和确定预期文件数
# 注意：这些变量仅用于信息显示和获取预期文件数，实际循环将遍历整个 extracted_dirs 数组
start_dir_full="${extracted_dirs[0]}"
end_dir_full="${extracted_dirs[$((${#extracted_dirs[@]} - 1))]}"

echo "将检查以下目录 (从 .sh 脚本文件名提取):"
echo "  起始目录名: $start_dir_full"
echo "  结束目录名: $end_dir_full"
echo "  总共要检查的目录数: ${#extracted_dirs[@]}"

echo  "--- 2. 获取预期文件数 ---"
# 检查 tree 命令是否存在
if ! command -v tree >/dev/null 2>&1; then
  echo "错误: 'tree' 命令未找到。请安装它以统计文件数。"
  exit 1
fi

# 检查起始目录是否存在于结果路径中
first_result_dir="$RNA_COUNT_RESULT_BASE_DIR/$start_dir_full"
if [ ! -d "$first_result_dir" ]; then
  echo "错误: 结果目录 '$first_result_dir' 不存在，无法确定预期文件数。"
  exit 1
fi

# 使用 tree 命令获取文件总数
# tree 命令的输出格式通常是 "X directories, Y files" 在最后一行
tree_output=$(tree "$first_result_dir" 2>/dev/null)
files_line=$(echo "$tree_output" | grep -oP '\d+ files$')

if [ -z "$files_line" ]; then
  echo "警告: 无法从 '$first_result_dir' 的 tree 输出中解析文件总数。请手动检查。"
  expected_files=0 # 设置为0，后续检查会失败
else
  expected_files=$(echo "$files_line" | grep -oP '^\d+')
fi

echo "预期文件数: $expected_files"
echo "-------------------------------------"

# 切换到结果目录
cd "$RNA_COUNT_RESULT_BASE_DIR" || { echo "错误: 无法切换到目录 '$RNA_COUNT_RESULT_BASE_DIR'"; exit 1; }

# 预期目录数（保持不变）
expected_dirs=7

echo "--- 检查目录和文件 ---"
echo "预期子目录数 (每个结果目录内): $expected_dirs"
echo "预期文件数 (每个结果目录内): $expected_files"
echo "------------------------"

# 循环遍历从 .sh 脚本文件名中提取的目录列表
for d in "${extracted_dirs[@]}"; do
  if [ ! -d "$d" ]; then
    echo "$d : MISSING DIR"
    continue
  fi
  n_dirs=$(find "$d" -maxdepth 1 -type d | wc -l) # 只统计当前目录下的直接子目录
  n_files=$(find "$d" -type f | wc -l)            # 统计当前目录及其所有子目录下的文件

  # 检查关键文件是否存在
  if [ -f "$d/rsem/${d}_rsem.genes.results" ]; then
    key="rsem.genes.results OK"
  else
    key="MISSING rsem.genes.results"
  fi

  # 打印并标注是否达标
  if [ "$n_dirs" -eq "$expected_dirs" ] && [ "$n_files" -eq "$expected_files" ]; then
    echo "$d : $n_dirs directories, $n_files files | $key | PASS"
  else
    echo "$d : $n_dirs directories, $n_files files | $key | FAIL (expected $expected_dirs dirs, $expected_files files)"
  fi
done