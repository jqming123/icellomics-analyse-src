#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: GenerateRNAShell_v3.py

功能描述:
    此脚本用于为RNA-seq数据分析批量生成SLURM作业提交脚本。
    它读取一个样本列表文件, 为每个样本动态创建SLURM shell脚本, 
    以执行SRA到FASTQ转换、FASTQ压缩以及调用核心RNA-seq分析流程。
    支持自动检测现有的 fastq.gz 文件以跳过转换步骤。

使用方法:
    python3 GenerateRNAShell_v3.py -s sra_runid_prjid_ref.txt -o <输出目录>

参数:
    -s, --sample_file: 样本列表文件路径 (格式: SRR_ID\tProject_ID\tReference_Name)。
    -o, --output_dir:  生成的SLURM脚本 (.sh) 的输出目录。

依赖:
    - Python 3.x
    - fasterq-dump, pigz, bash
    - SLURM 作业调度系统
    - 核心分析脚本 (bulkRNA-seq_E4_v2.sh)
    - 预构建的参考基因组索引

重要提示:
    请务必在运行前检查同一目录下的config.py文件是否配置
"""

import os
import sys
import argparse
import textwrap

# 从独立的配置文件中导入 CONFIG 字典
try:
    from config import CONFIG
except ImportError:
    print("错误: 找不到配置文件 'config.py'。请确保它与此脚本在同一目录下。")
    sys.exit(1)

# -----------------------------------------------------------------
# 核心逻辑区域: 通常无需修改此部分代码
# -----------------------------------------------------------------

def generate_script_content(srr_name, project_name, ref_name, config):
    """使用f-string模板生成单个shell脚本的内容"""

    # 从配置字典中提取参数
    paths = config['paths']
    slurm = config['slurm_settings']
    params = config['tool_params']
    # 从配置中获取指定物种的参考基因组路径
    try:
        ref = config['reference_genomes'][ref_name]
    except KeyError:
        print(f"严重错误: 在CONFIG中未找到名为 '{ref_name}' 的参考基因组。")
        return None 

    # 构建动态路径
    sra_file_path = os.path.join(paths['sra_data_root'], project_name, srr_name, f"{srr_name}.sra")
    fastq_output_dir = os.path.join(paths['project_results_root'], srr_name, 'reads')
    working_dir = os.path.join(paths['project_results_root'], srr_name)
    temp_dir = os.path.join(paths['project_results_root'], srr_name, 'temp')
    # 确保输出目录存在
    os.makedirs(fastq_output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(slurm['log_dir'], exist_ok=True) # 确保日志目录存在


    # 使用textwrap.dedent来移除多行字符串的公共前导空白，使代码更整洁
    # %x 会被SLURM替换为作业名称(srr_name), %j 会被替换为作业ID
    log_file_path = os.path.join(slurm['log_dir'], project_name, f'{srr_name}_count_%j.log')
    os.makedirs(os.path.dirname(log_file_path), exist_ok=True)

    script_template = textwrap.dedent(f"""\
        #!/bin/bash
        # --- SLURM Directives ---
        #SBATCH --partition={slurm['partition']}
        #SBATCH --job-name={srr_name}_count
        #SBATCH --nodes=1
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task={params['threads']}
        #SBATCH --mem={params['memory_gb']}G
        #SBATCH --time={slurm['time']}
        #SBATCH --output={log_file_path}

        # --- Script Body ---
        echo "开始处理样本: {srr_name}"
        echo "脚本启动时间: $(date)"
        echo "作业ID: $SLURM_JOB_ID"
        echo "在节点: $SLURM_JOB_NODELIST 上运行"

        # 设置执行环境，确保可复现性
        export PATH={paths['custom_bin_path']}:$PATH
        # 设置HOME目录
        export HOME=/home/gaoxiaojing
        # 设置工作目录
        export WORKING_DIR={working_dir}
        mkdir -p "$WORKING_DIR" # 确保目录存在
        
        if ! cd "$WORKING_DIR"; then
            echo "严重错误: 无法切换到工作目录: $WORKING_DIR"
            # 在退出前可以尝试输出更多诊断信息
            echo "当前所在目录是: $(pwd)"
            echo "目标目录权限信息:"
            ls -ld "$WORKING_DIR"
            exit 1
        fi
        echo "已切换到工作目录: $(pwd)"
        
        # --- 调试信息开始 ---
        echo "--- 调试信息 ---"
        echo "当前PATH: $PATH"
        echo "当前HOME目录: $HOME" # 再次打印HOME，确认修改已生效
        echo "当前工作目录 (pwd): $(pwd)"
        
        # 1. Check for existing fastq.gz files
        USE_EXISTING_FASTQ=false
        
        if ls "{fastq_output_dir}"/*.fastq.gz >/dev/null 2>&1; then
            echo "=================================================================="
            echo "[INFO] 检测到在 {fastq_output_dir} 目录中已存在 FASTQ.GZ 文件:"
            ls -la "{fastq_output_dir}"/*.fastq.gz
            echo "将直接跳过 SRA 转换 (fasterq-dump) 和压缩 (pigz) 步骤。"
            echo "=================================================================="
            USE_EXISTING_FASTQ=true
        fi

        if [ "$USE_EXISTING_FASTQ" = "false" ]; then
            echo "未检测到现有的 fastq.gz 文件，尝试使用 SRA 文件进行转换..."
            echo "检查 fasterq-dump 可执行文件路径和版本:"
            which fasterq-dump
            fasterq-dump --version
        
            # SRA to FASTQ conversion
            if [ -f "{sra_file_path}" ]; then
                # 设置临时目录环境变量
                export TMPDIR="{temp_dir}"
                export TMP="$TMPDIR"
                export TEMP="$TMPDIR"
                # 确保临时目录存在且可写
                if [ ! -d "$TMPDIR" ]; then
                    mkdir -p $TMPDIR
                    echo "现在创建临时目录: $TMPDIR"
                fi
                echo "临时目录 ($TMPDIR) 权限检查:"
                ls -ld "$TMPDIR" 
                
                fasterq-dump \\
                    -e {params['threads']} \\
                    --split-3 \\
                    --temp "$TMPDIR" \\
                    -O "{fastq_output_dir}" \\
                    "{sra_file_path}"
                    
                # 检查fasterq-dump是否成功
                if [ $? -eq 0 ]; then
                    echo "fasterq-dump 执行成功"
                    echo "生成的文件列表:"
                    ls -la {fastq_output_dir}/
                else
                    echo "错误: fasterq-dump 执行失败"
                    exit 1
                fi
            else
                echo "=================================================================="
                echo "严重错误: 既未找到 SRA 原始文件: {sra_file_path}"
                echo "          也未在目标目录找到现有的 FASTQ.GZ 文件: {fastq_output_dir}/*.fastq.gz"
                echo "提示: 请将手动下载的 .fastq.gz 文件存放到以下目录后重新提交作业："
                echo "      {fastq_output_dir}/"
                echo "=================================================================="
                exit 1
            fi
            
            # 2. Compress FASTQ files
            echo "步骤2: 使用 pigz 并行压缩FASTQ文件..."
            fastq_files=("{fastq_output_dir}"/*.fastq)

            if [ ${{#fastq_files[@]}} -eq 0 ]; then
                echo "错误: 未找到 FASTQ 文件，无法压缩"
                exit 1
            fi
            
            pigz -p {params['threads']} "${{fastq_files[@]}}"
            
            # 清理临时文件
            echo "清理临时文件..."
            rm -rf {temp_dir}
        fi

        # 3. Run main RNA-seq analysis pipeline (内部成功后会自动清理中间文件)
        echo "步骤3: 运行核心分析流程 bulkRNA-seq_E4_v2.sh..."
        bash {paths['main_script_path']} \\
            -i {paths['project_results_root']} \\
            -o {paths['project_results_root']} \\
            -s {srr_name} \\
            -c {params['threads']} \\
            -m {params['memory_gb']} \\
            -S "{ref['star_index']}" \\
            -G "{ref['kallisto_gene_idx']}" \\
            -T "{ref['kallisto_transcript_idx']}" \\
            -R "{ref['rsem_ref_prefix']}"

        echo "脚本结束时间: $(date)"
        echo "样本 {srr_name} 处理完成。"
    """)
    return script_template


def main():
    """主函数：解析参数并驱动脚本生成"""
    parser = argparse.ArgumentParser(
        description="为RNA-seq分析流程批量生成SLURM作业提交脚本。所有配置均在此脚本顶部修改。",
        # RawTextHelpFormatter 可以在帮助信息中保留换行符
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-s', '--sample_file', required=True,
                        help="输入的样本列表文件路径。\n"
                             "格式: SRR_ID\\tProject_ID\\tReference_Name (tab分隔)\n"
                             "其中 Reference_Name 必须与CONFIG字典中'reference_genomes'的键匹配。\n"
                             "例如:\n"
                             "SRR123456\\tProjectA\\tCriGri-PICRH-1.0\n"
                             "SRR789012\\tProjectB\\thg38_Ensemble")
    parser.add_argument('-o', '--output_dir', required=True, help="生成的.sh脚本的输出目录。")

    args = parser.parse_args()

    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"开始生成脚本，将保存到: {os.path.abspath(args.output_dir)}")
    print(
        f"使用的资源配置：分区(Partition)='{CONFIG['slurm_settings']['partition']}', "
        f"线程数={CONFIG['tool_params']['threads']}, 内存={CONFIG['tool_params']['memory_gb']}G"
    )
    print("将根据样本文件中的第三列选择参考基因组。")
    
    try:
        with open(args.sample_file, 'r') as f:
            for line in f:
                # 忽略空行 and 注释行
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                # 现在需要解析三列
                parts = line.split('\t') 
                if len(parts) != 3:
                    print(f"警告: 跳过格式不正确的行 -> '{line}' (需要3列: SRR_ID, Project_ID, Reference_Name)")
                    continue

                srr_name, project_name, ref_name = parts

                # 调用修改后的函数，传入参考基因组名称
                # (假设 generate_script_content 函数也已相应修改)
                script_content = generate_script_content(srr_name, project_name, ref_name.strip(), CONFIG)
                
                # 检查 generate_script_content 是否成功返回内容 (例如，如果ref_name无效)
                if script_content is None:
                    print(f"\n[致命错误]: 未能在配置中找到 '{ref_name}'。")
                    print(f"发生错误的行为: {line}")
                    print("程序已终止，未生成后续脚本。")
                    sys.exit(1)  # 1 表示非正常退出

                # 写入文件
                script_path = os.path.join(args.output_dir, f"{srr_name}.sh")
                with open(script_path, 'w') as out_f:
                    out_f.write(script_content)
                
                print(f"  -> 已生成: {script_path} (使用参考: {ref_name})")

    except FileNotFoundError:
        print(f"错误: 样本文件 '{args.sample_file}' 未找到。")
        sys.exit(1)

    print("所有脚本生成完毕！")

if __name__ == "__main__":
    main()
