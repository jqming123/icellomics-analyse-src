import os
import sys
import argparse
import textwrap

# -----------------------------------------------------------------
# 配置区域: 所有可调整的参数都定义在这里
# -----------------------------------------------------------------
# 修改此处的字典值即可调整脚本行为，无需触碰下方的核心逻辑代码。
CONFIG = {
    # 1. 路径设置
    "paths": {
        "sra_data_root": '/gpfs/zhaowm_group/gaoxiaojing/CellLine/rna',
        "project_results_root": '/gpfs/zhaowm_group/gaoxiaojing/CellLine/projects_results/rnseq',
        "main_script_path": '/gpfs/zhaowm_group/gaoxiaojing/CellLine/script/bulkRNA-seq_E4.sh',
        "custom_bin_path": '/gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/bin:/p300s/zhaowm_group/kongdm/workspace/BIG/sra/.pixi/envs/default/bin:/p300s/zhaowm_group/tangbx/software/sratoolkit.3.1.0-centos_linux64/bin'
    },

    # 2. PBS作业调度系统设置
    "pbs_settings": {
        # "queue": 'q512G',
        "queue": 'core40',
        "walltime": '15:00:00',
        "error_dir": '/gpfs/zhaowm_group/gaoxiaojing/CellLine/sra/error/',
        "output_dir": '/gpfs/zhaowm_group/gaoxiaojing/CellLine/sra/output/',
        "log_dir": '/gpfs/zhaowm_group/gaoxiaojing/CellLine/sra/logs/',
        "scheduler_opts": '#HSCHED -s hschedd=hschedd'
    },

    # 3. 工具和资源参数
    "tool_params": {
        "threads": 8,
        "memory_gb": 60
    }
}


# -----------------------------------------------------------------
# 核心逻辑区域: 通常无需修改此部分代码
# -----------------------------------------------------------------

def generate_script_content(srr_name, project_name, config):
    """使用f-string模板生成单个shell脚本的内容"""

    # 从配置字典中提取参数
    paths = config['paths']
    pbs = config['pbs_settings']
    params = config['tool_params']

    # 构建动态路径
    sra_file_path = os.path.join(paths['sra_data_root'], project_name, srr_name, f"{srr_name}.sra")
    fastq_output_dir = os.path.join(paths['project_results_root'], srr_name, 'reads')
    working_dir = os.path.join(paths['project_results_root'], srr_name)
    temp_dir = os.path.join(paths['project_results_root'], srr_name, 'temp')
    # 确保输出目录存在
    os.makedirs(fastq_output_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)


    # 使用textwrap.dedent来移除多行字符串的公共前导空白，使代码更整洁
    script_template = textwrap.dedent(f"""\
        #!/bin/bash
        # --- PBS Directives ---
        #PBS -q {pbs['queue']}
        #PBS -N {srr_name}
        #PBS -l nodes=1:ppn={params['threads']}
        #PBS -l mem={params['memory_gb']}gb
        #PBS -l walltime={pbs['walltime']}
        #PBS -e {os.path.join(pbs['error_dir'], f'{srr_name}.err')}
        #PBS -o {os.path.join(pbs['output_dir'], f'{srr_name}.out')}
        {pbs['scheduler_opts']}

        # --- Script Body ---
        # 将所有输出重定向到专用日志文件，便于调试
        exec &> {os.path.join(pbs['log_dir'], f'{srr_name}.log')}

        echo "开始处理样本: {srr_name}"
        echo "脚本启动时间: $(date)"

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
        echo "当前工作目录 (pwd): $(pwd)" # **新增：检查当前工作目录**
        
        echo "检查 fasterq-dump 可执行文件路径和版本:"
        which fasterq-dump
        fasterq-dump --version
        
        echo "检查 ~/.ncbi/user-settings.mkfg 文件内容:"
        # 使用 $HOME 而不是硬编码路径，以确保与当前HOME设置一致
        

    
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
        echo "步骤3: 运行核心分析流程 bulkRNA-seq_E4.sh..."
        bash {paths['main_script_path']} \\
            {paths['project_results_root']} \\
            {paths['project_results_root']} \\
            {srr_name} \\
            {params['threads']} \\
            {params['memory_gb']}

        echo "脚本结束时间: $(date)"
        echo "样本 {srr_name} 处理完成。"
    """)
    return script_template


def main():
    """主函数：解析参数并驱动脚本生成"""
    parser = argparse.ArgumentParser(
        description="为RNA-seq分析流程批量生成PBS作业提交脚本。所有配置均在此脚本顶部修改。",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-s', '--sample_file', required=True,
                        help="输入的样本列表文件路径。\n格式: SRR_ID\\tProject_ID (tab分隔)")
    parser.add_argument('-o', '--output_dir', required=True, help="生成的.sh脚本的输出目录。")

    args = parser.parse_args()

    # 确保输出目录存在
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"开始生成脚本，将保存到: {args.output_dir}")
    print(
        f"使用的配置：队列='{CONFIG['pbs_settings']['queue']}', 线程数={CONFIG['tool_params']['threads']}, 内存={CONFIG['tool_params']['memory_gb']}gb")

    try:
        with open(args.sample_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) != 2:
                    print(f"警告: 跳过格式不正确的行 -> '{line}'")
                    continue

                srr_name, project_name = parts

                # 使用全局配置字典生成脚本内容
                script_content = generate_script_content(srr_name, project_name, CONFIG)

                # 写入文件
                script_path = os.path.join(args.output_dir, f"{srr_name}.sh")
                with open(script_path, 'w') as out_f:
                    out_f.write(script_content)

                print(f"  -> 已生成: {script_path}")

    except FileNotFoundError:
        print(f"错误: 样本文件 '{args.sample_file}' 未找到。")
        sys.exit(1)

    print("所有脚本生成完毕！")


if __name__ == "__main__":
    main()
