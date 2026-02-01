import os  # 用于处理文件和目录
import time  # 用于添加时间延迟
import subprocess  # 用于执行外部命令



def get_running_tasks_count(username):
    """
    获取指定用户的正在运行的任务数量。

    Args:
        username (str): 要查询的用户名。

    Returns:
        int: 正在运行的任务数量。 如果执行命令出错，则返回 20。
    """
    cmd = f"squeue -u {username}"  # 构建查询任务状态的命令
    try:
        output = subprocess.check_output(cmd, shell=True, text=True)  # 执行命令并获取输出
        print(output) # 打印qstat命令的输出，方便调试
        tasks = [item for item in output.strip().split('\n') if "gaoxiao" in item]  # 过滤包含用户名的行，统计任务数量
        task_count = len(tasks)  # 获取任务数量
        print(f"running wget number is {task_count}")  # 打印当前运行的任务数量
        return task_count  # 返回任务数量
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")  # 打印错误信息
        return 20  # 如果执行命令出错，则返回 20，表示超过任务数量限制


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
    except subprocess.CalledProcessError as e:
        with open(err_file, "a") as f:  # 打开错误日志文件
            f.write(f"Failed to submit download task: {e}")  # 写入错误信息


def main():
    """
    主函数，负责任务的监控和提交。
    """
    MAX_TASKS=50
    script_folder = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB9185_CHO/02_jobs/vep_debug_jobs_2/"  # 包含作业脚本的目录
    script_list = os.listdir(script_folder)  # 获取目录下所有脚本的文件名
    while len(script_list) > 0:  # 循环直到所有任务都已提交
        current_tasks = get_running_tasks_count("gaoxiaojing")  # 获取当前正在运行的任务数量
        if current_tasks < MAX_TASKS:  # 如果当前任务数量小于 20，则提交新的任务
            for i in range(0, MAX_TASKS - current_tasks):  # 循环提交任务
                current_file_name = script_list[0]  # 获取要提交的脚本文件名
                submit_task(f"{script_folder}{current_file_name}",f"{script_folder}{current_file_name}.aterr")  # 提交任务
                script_list.pop(0)  # 从列表中移除已提交的任务
                time.sleep(3)  # 暂停 3 秒
            time.sleep(300)  # 提交完一批任务后，暂停 300 秒
        else:
            time.sleep(300)  # 如果当前任务数量大于等于 20，则暂停 300 秒


if __name__ == "__main__":
    main()  # 当脚本直接运行时，执行 main 函数