#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
本脚本用于监控 SLURM 队列中的任务数量，并自动提交新的任务，直到所有任务都被提交。

推荐使用 nohup 命令在后台运行此脚本，并将标准输出和标准错误重定向到日志文件，
以确保脚本在会话断开后仍然继续运行，并方便后续查看运行日志。

示例命令：
nohup python -u autosubmit_v2.py > /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB39258_CHO/03_logs/autosubmit_logs/vep_debug_final_round2.log 2>&1 &

命令解释：
- `nohup`: 运行命令，即使终端关闭，进程也不会停止。
- `python -u`: 以无缓冲模式运行 Python 脚本，这意味着输出会立即写入文件，而不是等到缓冲区满。
- `autosubmit_v2.py`: 要执行的 Python 脚本。
- `>`: 将标准输出重定向到指定文件。
- `/hpcdisk1/.../03_logs/autosubmit_logs/vep_debug_final_round2.log`: 日志文件的路径。
- `2>&1`: 将标准错误（文件描述符2）重定向到与标准输出（文件描述符1）相同的位置。
- `&`: 将命令放到后台运行，释放当前终端。
"""

import os  # 用于处理文件和目录
import time  # 用于添加时间延迟
import subprocess  # 用于执行外部命令
import datetime  # 导入 datetime 模块，用于记录时间

def get_running_tasks_count(username, max_tasks):
    """
    获取指定用户的正在运行的任务数量。

    Args:
        username (str): 要查询的用户名。
        max_tasks(int): 最大任务数

    Returns:
        int: 正在运行的任务数量。 如果执行命令出错，则返回 -1。
    """

    cmd = f"squeue -u {username}"  # 构建查询任务状态的命令

    print(f"Attempting to run command: {cmd}")

    try:
        # 使用 subprocess.run 来同时捕获 stdout 和 stderr
        # capture_output=True 捕获 stdout 和 stderr
        # text=True 将输出解码为文本
        # check=True 会在命令返回非零退出状态时抛出 CalledProcessError(目前没有这个参数)
        result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
        output = result.stdout
        if result.returncode != 0:
            print(f"Warning: squeue returned non-zero exit code {result.returncode}, but attempting to parse output anyway.")
            if result.stderr:
                print(f"squeue stderr:\n{result.stderr}")

        # 改进任务数量的解析逻辑
        lines = output.strip().split('\n')
        job_lines = [line for line in lines if "gaoxiaoj" in line]
        task_count = len(job_lines)
        print(f"running tasks for {username}: {task_count}")
        return task_count

    except FileNotFoundError:
        print("Error: 'squeue' command not found in PATH. Please ensure it's installed and accessible.")
        return -1
    except Exception as e:
        print(f"Unexpected error: {e}")
        return -1

def submit_task(file_name, err_file):
    """
    提交一个任务。

    Args:
        file_name (str): 要提交的任务脚本文件名。
        err_file (str): 错误日志文件名。
    """
    cmd = f"sbatch {file_name}"  # 构建提交任务的命令
    try:
        subprocess.run(cmd, shell=True, check=True)  # 执行命令
        print(f"Successfully submit {file_name}")
    except subprocess.CalledProcessError as e:
        with open(err_file, "a") as f:  # 打开错误日志文件
            f.write(f"Failed to submit download task: {e}")  # 写入错误信息


def main():
    """
    主函数，负责任务的监控和提交。
    """
    # 记录开始时间
    start_time = datetime.datetime.now()
    print(f"Script started at: {start_time}")

    max_tasks=50
    username = "gaoxiaojing"
    script_folder = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/EQ_jobs/PRJNA1008690_HEK293"
    print(f"current script folder: {script_folder}, username:{username}, max tasks: {max_tasks}")

    script_list = [f for f in os.listdir(script_folder) if f.endswith('.sh')]
    script_list.sort() 
    print(f"Found {len(script_list)} .sh scripts to process.")

    while len(script_list) > 0:  # 循环直到所有任务都已提交
        print("Time: ",datetime.datetime.now())
        current_tasks = get_running_tasks_count(username,max_tasks)  # 获取当前正在运行的任务数量
        if 0 <= current_tasks < max_tasks:  # 如果当前任务数量小于 max_tasks，则提交新的任务
            sub_number=min(max_tasks - current_tasks,len(script_list))
            for i in range(0, sub_number):  # 循环提交任务
                current_file_name = script_list[0]  # 获取要提交的脚本文件名
                submit_task(f"{script_folder}/{current_file_name}",f"{script_folder}/{current_file_name}.aterr")  # 提交任务
                script_list.pop(0)  # 从列表中移除已提交的任务
                time.sleep(3)  # 暂停 3 秒
            if len(script_list) > 0:
                print(f"本次提交{sub_number}个任务，还剩{len(script_list)}个任务，三分钟后再尝试提交")
                time.sleep(180)  # 提交完一批任务后，暂停 300 秒
            else:
                print(f"本次提交{sub_number}个任务，还剩{len(script_list)}个任务")

        elif current_tasks >= max_tasks:
            print(f"当前队列中的任务数大于等于{max_tasks}，三分钟后再尝试提交")
            time.sleep(180)  # 如果当前任务数量大于等于 max_tasks，则暂停 300 秒

        elif current_tasks < 0:
            print("上一步出错了！")
            break

    print(f"{script_folder}中的所有任务都已提交！")

    # 记录结束时间
    end_time = datetime.datetime.now()
    print(f"Script finished at: {end_time}")

    # 计算并打印总运行时间
    total_duration = end_time - start_time
    print(f"Total script duration: {total_duration}")

if __name__ == "__main__":
    main()
