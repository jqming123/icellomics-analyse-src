#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
脚本名称: GenerateRNAShell_v3.py

功能描述:
    此脚本用于为RNA-seq数据分析批量生成SLURM作业提交脚本。
    它读取一个样本列表文件, 为每个样本动态创建SLURM shell脚本, 
    以执行SRA到FASTQ转换、FASTQ压缩以及调用核心RNA-seq分析流程。

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
    请务必在运行前检查并修改脚本顶部的 `CONFIG` 字典以适应您的环境和需求。
"""



import os
import sys
import argparse
import textwrap

# -----------------------------------------------------------------
# 配置区域: 所有可调整的参数都定义在这里
# -----------------------------------------------------------------
# 修改此处的字典值即可调整脚本行为，无需触碰下方的核心逻辑代码。
REF_BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"

CONFIG = {
    # 1. 路径设置
    "paths": {
        "sra_data_root": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_rawdata',
        "project_results_root": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result',
        "main_script_path": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/trs_script/bulkRNA-seq_E4_v2.sh',
        # "custom_bin_path": '/gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/bin:/p300s/zhaowm_group/kongdm/workspace/BIG/sra/.pixi/envs/default/bin:/p300s/zhaowm_group/tangbx/software/sratoolkit.3.1.0-centos_linux64/bin'
        # pixi似乎不应该按下面的用法使用
        # sratoolkit的路径已经在PATH变量里了
        "custom_bin_path": '/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/bin:/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/pixi_0.59.0/trs_env/.pixi/envs/default/bin'
    },

    # 2. SLURM作业调度系统设置
    "slurm_settings": {
        # "partition": 'vmcore128',
        "partition": 'corexd192',
        "time": '7-00:00:00',  # D-HH:MM:SS 格式 
        "log_dir": '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/logs',
    },


    # 3. 工具和资源参数
    "tool_params": {
        "threads": 8,
        "memory_gb": 60
    },

    # 4. 参考基因组设置
    "reference_genomes": {
        "CriGri-PICRH-1.0": { # 中国仓鼠卵巢细胞（已弃用）
            "star_index": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","kallisto.index","GCF_003668045.3_CriGri-PICRH-1.0_genomic.gene.fa.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","kallisto.index","GCF_003668045.3_CriGri-PICRH-1.0_genomic.transcript.fa.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0","rsem.index","reference")
        },
        "CH_Ensemble": { # 中国仓鼠卵巢细胞Ensemble
            "star_index": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensemble","star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensemble","kallisto.index","CriGri-PICRH-1.0.115.gene.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensemble","kallisto.index","CriGri-PICRH-1.0.115.transcript.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR,"CriGri-PICRH-1.0_Ensemble","rsem.index","reference")
        },
        "hg38_Ensemble": { # 人类
            "star_index": os.path.join(REF_BASE_DIR,"hg38_Ensemble","star.index"),
            "kallisto_gene_idx": os.path.join(REF_BASE_DIR,"hg38_Ensemble","kallisto.index","Homo_sapiens.GRCh38.dna_sm.primary_assembly.gene.fa.idx"),
            "kallisto_transcript_idx": os.path.join(REF_BASE_DIR,"hg38_Ensemble","kallisto.index","Homo_sapiens.GRCh38.dna_sm.primary_assembly.transcript.fa.idx"),
            "rsem_ref_prefix": os.path.join(REF_BASE_DIR,"hg38_Ensemble","rsem.index","reference")
        },
        # 在这里添加更多物种...
    }


}


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
        # 可以选择优雅地退出或跳过
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
        
        echo "检查 fasterq-dump 可执行文件路径和版本:"
        which fasterq-dump
        fasterq-dump --version
    
        # 1. SRA to FASTQ conversion
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
            echo "错误: SRA文件不存在: {sra_file_path}"
            exit 1
        fi
        
        # 2. Compress FASTQ files
        echo "步骤2: 使用 pigz 并行压缩FASTQ文件..."
        pigz -p {params['threads']} {fastq_output_dir}/*
        
        # 清理临时文件
        echo "清理临时文件..."
        rm -rf {temp_dir}

        # 3. Run main RNA-seq analysis pipeline
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
    # os.makedirs(args.output_dir, exist_ok=True)

    print(f"开始生成脚本，将保存到: {os.path.abspath(args.output_dir)}")
    print(
        f"使用的资源配置：分区(Partition)='{CONFIG['slurm_settings']['partition']}', "
        f"线程数={CONFIG['tool_params']['threads']}, 内存={CONFIG['tool_params']['memory_gb']}G"
    )
    print("将根据样本文件中的第三列选择参考基因组。")
    
    try:
        with open(args.sample_file, 'r') as f:
            for line in f:
                # 忽略空行和注释行
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
                    print(f"  -> 错误: 未能为 {srr_name} 生成脚本。请检查参考基因组名称 '{ref_name}' 是否在CONFIG中正确定义。")
                    continue

                # 写入文件
                script_path = f"{srr_name}.sh"
                with open(script_path, 'w') as out_f:
                    out_f.write(script_content)

                # 在输出中明确指出使用了哪个参考基因组，便于核对
                print(f"  -> 已生成: {script_path} (使用参考: {ref_name})")

    except FileNotFoundError:
        print(f"错误: 样本文件 '{args.sample_file}' 未找到。")
        sys.exit(1)

    print("所有脚本生成完毕！")

if __name__ == "__main__":
    main()
