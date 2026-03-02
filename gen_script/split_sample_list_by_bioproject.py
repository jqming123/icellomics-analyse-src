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

        # 目录: <BioProject>_<细胞系名称>/00_data
        target_dir = os.path.join(root_dir, f"{prj}_{cell_line}", "00_data")
        os.makedirs(target_dir, exist_ok=True)

        # 保存文件
        output_path = os.path.join(target_dir, "sample_list.txt")
        subset.to_csv(output_path, sep='\t', index=False, header=False)
        
        print(f"  [成功] -> {output_path} ({len(subset)} 行)")

if __name__ == "__main__":
    main()
