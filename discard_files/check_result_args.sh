#!/bin/bash

# 检查是否提供了正确的参数数量
if [ "$#" -ne 3 ]; then
  echo "用法: $0 <预期文件数> <起始目录名> <结束目录名>"
  echo "示例: $0 25 SRR15559102 SRR15559121"
  echo "示例: $0 27 SRR15111551 SRR15111580"
  exit 1
fi
# cd /gpfs/zhaowm_group/gaoxiaojing/CellLine/projects_results/rnseq
cd /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result

# 从命令行参数获取值
expected_files=$1
start_dir_full=$2
end_dir_full=$3

# 预期目录数（保持不变）
expected_dirs=8

# --- 解析起始和结束目录名，提取公共前缀和数字范围 ---
common_prefix=""
min_len=${#start_dir_full}
if (( ${#end_dir_full} < min_len )); then min_len=${#end_dir_full}; fi

# 找到最长公共前缀
for (( i=0; i<min_len; i++ )); do
  if [[ "${start_dir_full:$i:1}" == "${end_dir_full:$i:1}" ]]; then
    common_prefix+="${start_dir_full:$i:1}"
  else
    break
  fi
done

# 提取数字部分
start_num_str=${start_dir_full#$common_prefix}
end_num_str=${end_dir_full#$common_prefix}

# 验证提取出的部分是否为纯数字
if ! [[ "$start_num_str" =~ ^[0-9]+$ ]] || ! [[ "$end_num_str" =~ ^[0-9]+$ ]]; then
  echo "错误: 无法从 '$start_dir_full' 和 '$end_dir_full' 中提取有效的数字范围。"
  echo "请确保目录名遵循 '前缀+数字' 的模式，且数字部分是连续的。"
  exit 1
fi

start_num=$((10#$start_num_str)) # 强制十进制，避免八进制解释
end_num=$((10#$end_num_str))     # 强制十进制
num_padding_len=${#start_num_str} # 记录原始数字部分的长度，用于补齐前导零

echo "--- 检查目录和文件 ---"
echo "预期目录数: $expected_dirs"
echo "预期文件数: $expected_files"
echo "检查范围: $start_dir_full 到 $end_dir_full"
echo "公共前缀: '$common_prefix'"
echo "数字范围: $start_num 到 $end_num (补齐长度: $num_padding_len)"
echo "------------------------"

# 循环遍历指定范围的目录
for (( i=$start_num; i<=$end_num; i++ )); do
  # 根据原始数字部分的长度补齐前导零，构建完整的目录名
  current_num_padded=$(printf "%0*d" "$num_padding_len" "$i")
  d="$common_prefix$current_num_padded"

  if [ ! -d "$d" ]; then
    echo "$d : MISSING DIR"
    continue
  fi
  n_dirs=$(find "$d" -type d | wc -l)
  n_files=$(find "$d" -type f | wc -l)

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