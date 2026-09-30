#!/usr/bin/env python3

"""
脚本名称: split_epigen_metadata_by_bioproject.py
描述:
    该脚本用于根据 BioProject ID 拆分表观相关项目（epigen_projects）的测序数据元数据表格。
    它读取指定的 TSV 文件，提取 'BioSample' 和 'Run' 信息。
    注意：脚本会检查对应的项目目录是否存在，仅对已存在的项目目录生成 0_data/sample_run_map.tsv。

输入文件格式 (TSV):
    必须包含列: 'BioProject', 'Run', 'BioSample' (有表头)
    路径模式: /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects/run_spl_prj_tables/{cell_line}_run_spl_prj.tsv

输出结果:
    目录结构: <root_dir>/<BioProject>_<cell_line>/0_data/ （仅当项目目录存在时创建）
    文件名: sample_run_map.tsv (第一列 BioSample，第二列 Run，无表头，Linux 换行符)

使用方法:
    python split_epigen_metadata_by_bioproject.py <cell_line>

参数说明:
    cell_line : 目标细胞系的名称 (例如: HEK293, K562, MRC-5)

日期: 2026-06-03
"""

import pandas as pd
import os
import argparse

def main():
    # --- 1. 设置命令行参数解析 ---
    parser = argparse.ArgumentParser(description="根据 BioProject 拆分表观项目数据表格（仅处理已存在的目录）")
    parser.add_argument("cell_line", help="输入的细胞系名称 (例如: HEK293, K562)")
    args = parser.parse_args()

    cell_line = args.cell_line
    
    # --- 2. 路径配置 ---
    root_dir = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects"
    input_file = os.path.join(root_dir, "run_spl_prj_tables", f"{cell_line}_run_spl_prj.tsv")

    # --- 3. 逻辑执行 ---
    if not os.path.exists(input_file):
        print(f"错误: 找不到输入文件 {input_file}")
        return

    print(f"正在处理表观细胞系项目: {cell_line} ...")
    
    # 读取表格
    df = pd.read_csv(input_file, sep='\t')

    # 检查必要的列是否存在
    required_columns = ['BioProject', 'BioSample', 'Run']
    for col in required_columns:
        if col not in df.columns:
            print(f"错误: 输入文件中缺少必需的列 '{col}'")
            return

    # 按 BioProject 分组处理
    bioprojects = df['BioProject'].unique()

    for prj in bioprojects:
        # 定义核心项目目录路径: <root_dir>/<BioProject>_<cell_line>
        project_dir = os.path.join(root_dir, f"{prj}_{cell_line}")
        
        # --- 关键修改：检查项目主目录是否存在 ---
        if not os.path.exists(project_dir):
            print(f"  [跳过] -> 目录不存在: {project_dir}")
            continue

        # 筛选两列数据，注意顺序：第一列 BioSample，第二列 Run
        subset = df[df['BioProject'] == prj][['BioSample', 'Run']]
        
        # 去除前后的空格并转换为字符串
        subset['BioSample'] = subset['BioSample'].astype(str).str.strip()
        subset['Run'] = subset['Run'].astype(str).str.strip()

        # 目标目录: <project_dir>/0_data
        target_dir = os.path.join(project_dir, "0_data")
        os.makedirs(target_dir, exist_ok=True)

        # 保存文件 (文件名: sample_run_map.tsv，强制指定换行符为 '\n')
        output_path = os.path.join(target_dir, "sample_run_map.tsv")
        subset.to_csv(output_path, sep='\t', index=False, header=False, lineterminator='\n')
        
        print(f"  [成功] -> {output_path} ({len(subset)} 行)")

if __name__ == "__main__":
    main()
