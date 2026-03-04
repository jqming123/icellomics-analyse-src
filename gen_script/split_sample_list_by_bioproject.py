#!/usr/bin/env python3

"""
脚本名称: split_sequencing_data_by_bioproject.py
描述:
    该脚本用于根据 BioProject ID 拆分特定细胞系的测序数据元数据表格。
    它读取指定的 TSV 文件，提取 'Run' 和 'BioSample' 信息，并为每个 BioProject 
    创建标准的目录结构和样本列表文件。

输入文件格式 (TSV):
    必须包含列: 'BioProject', 'Run', 'BioSample' (有表头)
    路径模式: /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/run_spl_prj_tables/{cell_line}_run_spl_prj.tsv

输出结果:
    目录结构: <root_dir>/<BioProject>_<cell_line>/00_data/
    文件名: sample_list.txt (包含 Run 和 BioSample 两列，无表头，Linux 换行符)

使用方法:
    python split_sequencing_data_by_bioproject.py <cell_line>

参数说明:
    cell_line : 目标细胞系的名称 (例如: H9, Vero, HeLa)

依赖项:
    - pandas
    - os
    - argparse

日期: 2026-xx-xx
"""
import pandas as pd
import os
import argparse

def main():
    # --- 1. 设置命令行参数解析 ---
    parser = argparse.ArgumentParser(description="根据 BioProject 拆分测序数据表格")
    parser.add_argument("cell_line", help="输入的细胞系名称 (例如: H9, Vero)")
    args = parser.parse_args()

    cell_line = args.cell_line
    
    # --- 2. 路径配置 ---
    input_file = f"/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/run_spl_prj_tables/{cell_line}_run_spl_prj.tsv"
    root_dir = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects"

    # --- 3. 逻辑执行 ---
    if not os.path.exists(input_file):
        print(f"错误: 找不到输入文件 {input_file}")
        return

    print(f"正在处理细胞系: {cell_line} ...")
    
    # 读取表格
    df = pd.read_csv(input_file, sep='\t')

    # 按 BioProject 分组处理
    bioprojects = df['BioProject'].unique()

    for prj in bioprojects:
        # 筛选两列数据
        subset = df[df['BioProject'] == prj][['Run', 'BioSample']]
        subset['Run'] = subset['Run'].astype(str).str.strip()
        subset['BioSample'] = subset['BioSample'].astype(str).str.strip()

        # 目录: <BioProject>_<细胞系名称>/00_data
        target_dir = os.path.join(root_dir, f"{prj}_{cell_line}", "00_data")
        os.makedirs(target_dir, exist_ok=True)

        # 保存文件(强制指定换行符为 '\n' (Linux 格式))
        output_path = os.path.join(target_dir, "sample_list.txt")
        subset.to_csv(output_path, sep='\t', index=False, header=False, lineterminator='\n')
        
        print(f"  [成功] -> {output_path} ({len(subset)} 行)")

if __name__ == "__main__":
    main()
