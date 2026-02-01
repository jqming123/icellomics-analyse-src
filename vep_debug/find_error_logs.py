import os

def find_specific_error_logs(directory_path, log_file_prefix, specific_error_phrase):
    """
    在指定目录中查找包含特定错误短语的日志文件。

    Args:
        directory_path (str): 要搜索的目录路径。
        log_file_prefix (str): 日志文件名的前缀，例如 "vep_debug_chunk_"。
        specific_error_phrase (str): 要查找的特定错误短语（不区分大小写）。

    Returns:
        list: 包含特定错误短语的日志文件名的列表。
    """
    found_error_files = []
    
    # 检查目录是否存在
    if not os.path.isdir(directory_path):
        print(f"错误：目录 '{directory_path}' 不存在。")
        return []

    print(f"正在搜索目录：{directory_path}")
    print(f"日志文件前缀：{log_file_prefix}")
    print(f"正在查找短语：'{specific_error_phrase}'")

    # 准备要查找的短语的小写形式，以便进行不区分大小写的匹配
    target_phrase_lower = specific_error_phrase.lower()

    # 遍历目录中的所有文件和子目录
    for filename in os.listdir(directory_path):
        # 构建完整的文件路径
        file_path = os.path.join(directory_path, filename)

        # 检查是否是符合命名格式的日志文件
        if os.path.isfile(file_path) and \
           filename.startswith(log_file_prefix) and \
           filename.endswith(".log"):
            
            # 标记当前文件是否包含错误
            has_specific_error = False
            try:
                with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                    for line_num, line in enumerate(f, 1):
                        # 将行内容转换为小写，以便进行不区分大小写的匹配
                        lower_line = line.lower()
                        
                        # 检查行中是否包含特定的错误短语
                        if target_phrase_lower in lower_line:
                            print(f"  在文件 '{filename}' 中发现特定错误短语 (行 {line_num}): {line.strip()}")
                            found_error_files.append(filename)
                            has_specific_error = True
                            break # 找到特定短语即可，跳出当前文件的行循环，检查下一个文件
            except Exception as e:
                print(f"  读取文件 '{filename}' 时发生错误: {e}")
    
    return found_error_files

# 定义目标目录和日志文件前缀
target_directory = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB39258_CHO/03_logs/04_vep_annotation_overlap_debug"
log_file_prefix = "vep_debug_chunk_"

# 定义要查找的特定错误短语
specific_error_phrase_to_find = "Died in forked process"

# 调用函数查找错误日志
error_logs = find_specific_error_logs(target_directory, log_file_prefix, specific_error_phrase_to_find)

if error_logs:
    print(f"\n以下日志文件包含短语 '{specific_error_phrase_to_find}'：")
    for log_file in sorted(set(error_logs)): # 使用 set 去重并排序
        print(f"- {log_file}")
else:
    print(f"\n未在指定目录中找到包含短语 '{specific_error_phrase_to_find}' 的日志文件。")