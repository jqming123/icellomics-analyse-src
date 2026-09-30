#!/usr/bin/env python3
import sys
import os
import re

def process_group(group_str):
    """处理group字符串：括号、斜杠和空格替换成下划线，合并连续下划线，删除末尾下划线"""
    # 替换括号、斜杠和空格为下划线
    result = re.sub(r'[()/\s]', '_', group_str)
    # 合并连续的两个及以上下划线为一个
    result = re.sub(r'_+', '_', result)
    # 删除末尾的下划线
    result = result.rstrip('_')
    return result

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 split_deg_groups.py <cell_line>")
        print("Example: python3 split_deg_groups.py HT1080")
        sys.exit(1)
    
    cell_line = sys.argv[1]
    
    # 构建文件路径
    input_file = f"/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/groups/{cell_line}/{cell_line}_DEG_groups.tsv"
    output_file = f"/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/groups/{cell_line}/batch_generate_{cell_line}_DEGjobs_input.tsv"
    
    if not os.path.exists(input_file):
        print(f"Error: File not found - {input_file}")
        sys.exit(1)
    
    print(f"Processing DEG groups file: {input_file}")
    print(f"Output directory: {os.path.dirname(input_file)}")
    print("----------------------------------------")
    
    # 读取输入文件
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # 检查是否已经有folder_name和cell_line列
    header = lines[0].rstrip('\n').split('\t')
    has_folder_name = 'folder_name' in header
    has_cell_line = 'cell_line' in header
    
    # 根据表头确定各列的索引
    required_cols = ['Run', 'BioProject', 'group', 'condition']
    col_indices = {}
    for col in required_cols:
        if col in header:
            col_indices[col] = header.index(col)
        else:
            print(f"Error: Required column '{col}' not found in header!")
            sys.exit(1)
    
    # 处理每一行，生成新内容和提取内容
    new_lines = []
    extracted_lines = []
    folder_data = {}  # key: folder_name, value: list of (run, condition)
    
    for i, line in enumerate(lines):
        line = line.rstrip('\n')
        if i == 0:  # 表头
            if not has_folder_name and not has_cell_line:
                new_lines.append(line + '\tfolder_name\tcell_line')
            else:
                new_lines.append(line)
        else:
            parts = line.split('\t')
            if len(parts) > max(col_indices.values()):
                run = parts[col_indices['Run']]
                bioproject = parts[col_indices['BioProject']]
                group = parts[col_indices['group']]
                condition = parts[col_indices['condition']]
                
                # 确定folder_name（始终基于group重新计算，避免历史数据中的脏字符）
                processed_group = process_group(group)
                folder_name = f"{bioproject}_{processed_group}"
                
                # 确定cell_line值
                if has_cell_line:
                    cell_line_val = parts[header.index('cell_line')]
                else:
                    cell_line_val = cell_line
                
                # 添加新列（如果还没有）
                if not has_folder_name and not has_cell_line:
                    new_line = line + f'\t{folder_name}\t{cell_line_val}'
                    new_lines.append(new_line)
                else:
                    new_lines.append(line)
                
                # 提取用于新文件
                line_to_append = f"{folder_name}\t{cell_line_val}"
                if line_to_append not in extracted_lines:
                    extracted_lines.append(line_to_append)
                
                # 收集分组数据
                if folder_name not in folder_data:
                    folder_data[folder_name] = []
                folder_data[folder_name].append((run, condition))
    
    # 写回原文件（添加新列）
    if not has_folder_name or not has_cell_line:
        with open(input_file, 'w') as f:
            f.write('\n'.join(new_lines) + '\n')
    
    # 写入提取的新文件
    with open(output_file, 'w') as f:
        f.write('\n'.join(extracted_lines) + '\n')
    
    # 为每个folder_name创建目录和文件
    cell_line_dir = os.path.dirname(input_file)
    for folder_name, samples in folder_data.items():
        # 创建分组目录
        group_dir = os.path.join(cell_line_dir, folder_name)
        os.makedirs(group_dir, exist_ok=True)
        
        # 创建.tsv文件（实验设计信息）
        tsv_file = os.path.join(group_dir, f"{folder_name}.tsv")
        with open(tsv_file, 'w') as f:
            f.write("sample\tcondition\n")
            for run, condition in samples:
                f.write(f"{run}\t{condition}\n")
        
        # 创建.list文件（样本ID列表）
        list_file = os.path.join(group_dir, f"{folder_name}.list")
        with open(list_file, 'w') as f:
            for run, _ in samples:
                f.write(f"{run}\n")
        
        # 统计样本数
        case_count = sum(1 for _, cond in samples if cond == 'case')
        control_count = sum(1 for _, cond in samples if cond == 'control')
        
        print(f"Created group: {folder_name}")
        print(f"  - Directory: {group_dir}")
        print(f"  - Case samples: {case_count}")
        print(f"  - Control samples: {control_count}")
    
    print("----------------------------------------")
    print(f"Updated: {input_file}")
    print(f"Created: {output_file}")
    print("Processing completed!")

if __name__ == "__main__":
    main()