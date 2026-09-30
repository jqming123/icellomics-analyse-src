#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @File       :mergeRSEM
# @Time       :2024/6/6 16:42
# @Author     :zhoubw
# @Product    :DataSpell
# @Project    :encode
# @Version    :python 3.10.6
# @Description:
# @Usage      :

import os
import pandas as pd
import argparse

def merge_gene_expression(input_dir, sample_map_file, output_file, column_choice):
	# 映射用户输入的列名到实际列名
	column_map = {
		'TPM': 'TPM',
		'transcript_ids': 'transcript_id(s)',
		'length': 'length',
		'effective_length': 'effective_length',
		'expected_count': 'expected_count',
		'FPKM': 'FPKM'
	}

	if column_choice not in column_map:
		raise ValueError(f"Invalid column choice: {column_choice}. Must be one of {list(column_map.keys())}")

	actual_column = column_map[column_choice]

	# 读取技术重复映射文件：无表头，两列，run_id 和 merged_sample
	if not os.path.exists(sample_map_file):
		raise FileNotFoundError(f"The file {sample_map_file} does not exist.")

	sample_map = pd.read_csv(
		sample_map_file,
		sep='\t',
		header=None,
		names=['run_id', 'merged_sample'],
		dtype=str
	)

	sample_map = sample_map.dropna()
	sample_map = sample_map[(sample_map['run_id'].str.strip() != '') & (sample_map['merged_sample'].str.strip() != '')]

	if sample_map.empty:
		raise ValueError("The sample map is empty.")

	# 初始化一个空的DataFrame用于存储结果
	merged_df = pd.DataFrame()
	missing_runs = []      # 缺失或读取失败的 run
	failed_groups = []     # 整个合并样本都失败

	# expected_count 可以做技术重复加和；其他列不建议这样合并
	if actual_column != 'expected_count':
		print(f"Warning: technical replicate merging is primarily intended for expected_count, but current column is {actual_column}.")

	# 按合并后的新样本名分组，sort=False 保留 techrep.tsv 中第一次出现的顺序
	for merged_sample, group_df in sample_map.groupby('merged_sample', sort=False):
		group_merged = None

		for run_id in group_df['run_id']:
			file_path = os.path.join(input_dir, run_id, "rsem", f"{run_id}_rsem.genes.results")

			if not os.path.exists(file_path):
				print(f"Warning: The file {file_path} does not exist. Skipping run {run_id}.")
				missing_runs.append(run_id)
				continue

			try:
				df = pd.read_csv(file_path, sep='\t', usecols=['gene_id', actual_column])
			except Exception as e:
				print(f"Error reading {file_path}: {e}. Skipping run {run_id}.")
				missing_runs.append(run_id)
				continue

			# 组内统一列名，后面按 gene_id 做加和
			df = df.rename(columns={actual_column: run_id})

			if group_merged is None:
				group_merged = df
			else:
				group_merged = pd.merge(group_merged, df, on='gene_id', how='outer')

		# 如果这个新样本对应的所有 run 都失败了
		if group_merged is None:
			failed_groups.append(merged_sample)
			continue

		# 组内技术重复按 gene_id 求和
		value_cols = [col for col in group_merged.columns if col != 'gene_id']
		group_merged[value_cols] = group_merged[value_cols].fillna(0)
		group_merged[merged_sample] = group_merged[value_cols].sum(axis=1)
		group_merged = group_merged[['gene_id', merged_sample]]

		# 再并入总矩阵
		if merged_df.empty:
			merged_df = group_merged
		else:
			merged_df = pd.merge(merged_df, group_merged, on='gene_id', how='outer')

	# 保存结果到一个新的文件
	merged_df.to_csv(output_file, sep='\t', index=False)
	print(f"Merged gene expression matrix saved to {output_file}")

	# 输出不存在的样本ID
	if missing_runs:
		print("The following runs were not found or had errors and were skipped:")
		for run_id in missing_runs:
			print(f"- {run_id}")

	if failed_groups:
		print("The following merged samples had no valid runs and were skipped:")
		for sample in failed_groups:
			print(f"- {sample}")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Merge gene expression matrices.')
	parser.add_argument('-i', '--input_dir', required=True, help='Directory containing sample directories.')
	parser.add_argument('-l', '--sample_list', required=True, help='Two-column TSV without header: run_id and merged sample name.')
	parser.add_argument('-o', '--output_file', required=True, help='Output file to save the merged matrix.')
	parser.add_argument('-c', '--column', required=True, choices=['TPM', 'transcript_ids', 'length', 'effective_length', 'expected_count', 'FPKM'], help='Column to extract from each file.')

	args = parser.parse_args()

	merge_gene_expression(args.input_dir, args.sample_list, args.output_file, args.column)
