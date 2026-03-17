#!/bin/bash

# 1. 指定工作目录路径
WORK_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects"

# 2. 定义要处理的项目列表（空格分隔）
# 您可以在括号内添加或删除项目名称
PROJECTS=(
PRJDB10440_HEK293 PRJEB20596_HEK293 PRJEB23952_HEK293 PRJEB55318_HEK293 PRJEB78913_HEK293
)

# 切换到指定的工作目录，如果失败则退出
cd "$WORK_DIR" || { echo "错误: 无法进入目录 $WORK_DIR"; exit 1; }

echo "当前工作目录: $(pwd)"

# 遍历定义的项目列表
for project in "${PROJECTS[@]}"; do
    # 检查该项目目录是否存在
    if [ -d "$project" ]; then
        data_dir="$project/0_data"
        
        # 检查 0_data 目录是否存在
        if [ -d "$data_dir" ]; then
            echo "正在处理项目: $project"
            
            shopt -s nullglob
            # 查找 0_data 目录下所有的 fastq.gz 文件
            for fastq_path in "$data_dir"/*.fastq.gz; do
                filename=$(basename "$fastq_path")
                
                # 使用 sed 删除可能存在的 _1.fastq.gz, _2.fastq.gz 或 .fastq.gz 后缀
                # run_id=$(echo "$filename" | sed -r 's/(_[1][2]?)?\.fastq\.gz$//')
                # 兼容性更好的写法：先去掉 .fastq.gz，再去掉末尾的 _1 或 _2
                run_id=${filename%.fastq.gz}
                run_id=${run_id%_[12]}
                # ------------------------------------

                target_dir="$project/1_result/0_fastq/$run_id/reads"
                
                mkdir -p "$target_dir"
                mv "$fastq_path" "$target_dir/"
                
                echo "  [移动成功] $filename -> $target_dir"
            done
            shopt -u nullglob
        else
            echo "警告: 项目 $project 下未发现 0_data 目录，跳过。"
        fi
    else
        echo "错误: 项目目录 $project 不存在，跳过。"
    fi
done

echo "列表中的所有任务已尝试处理完毕。"
