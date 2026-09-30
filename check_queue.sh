#!/bin/bash

echo "正在获取 Slurm 系统中各个队列（分区）的排队作业和运行作业数量..."
echo ""

PARTITIONS=$(sinfo -h -o "%P" | tr ',' '\n' | sed 's/\*//g' | sort -u)

HEADER_PARTITION="队列 (Partition)  "
HEADER_PENDING="排队作业 (Pending)"
HEADER_RUNNING="运行作业 (Running)"

# 打印分隔符，每个段落使用 20 个横线
echo "+--------------------+--------------------+--------------------+"
printf "| %-18s | %-18s | %-18s |\n" "$HEADER_PARTITION" "$HEADER_PENDING" "$HEADER_RUNNING"
echo "+--------------------+--------------------+--------------------+"

for P in $PARTITIONS; do
    if [ -z "$P" ]; then
        continue
    fi

    PENDING_COUNT=$(squeue -p "$P" -t PENDING -h -o "%i" 2>/dev/null | wc -l)
    RUNNING_COUNT=$(squeue -p "$P" -t RUNNING,COMPLETING -h -o "%i" 2>/dev/null | wc -l)

    printf "| %-18s | %-18s | %-18s |\n" "$P" "$PENDING_COUNT" "$RUNNING_COUNT"
done
echo "+--------------------+--------------------+--------------------+"

echo ""
echo "注意："
echo "  - '排队作业' 指的是状态为 PENDING 的作业。"
echo "  - '运行作业' 指的是状态为 RUNNING 或 COMPLETING 的作业。"
echo "  - 如果某个分区没有作业，其计数会显示为 0。"
# '2>/dev/null' 用于抑制 squeue 在没有作业时可能输出的错误信息。