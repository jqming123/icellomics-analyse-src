import os

def get_dir_size(path):
    """递归计算目录大小（字节）"""
    total_size = 0
    try:
        for dirpath, dirnames, filenames in os.walk(path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                if not os.path.islink(fp):
                    total_size += os.path.getsize(fp)
    except Exception as e:
        print(f"无法访问目录 {path}: {e}")
    return total_size

def bytes_to_gb(bytes_size):
    """将字节转换为 GB 格式，保留 3 位小数"""
    return bytes_size / (1024**3)

def main():
    base_dir = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects"
    exclude_dirs = {"CHO_NCBIRef_MergedCache", "dna_sra"}
    suffixes = ["HT1080", "HEK293", "CHO_E"]
    
    # 初始化统计字典
    summary = {suffix: 0 for suffix in suffixes}
    
    if not os.path.exists(base_dir):
        print(f"错误: 路径 {base_dir} 不存在")
        return

    # 打印表头
    print(f"{'Directory Name':<50} | {'Size (GB)':<15}")
    print("-" * 70)

    # 遍历一级子目录
    for entry in sorted(os.listdir(base_dir)):
        # 排除指定目录
        if entry in exclude_dirs:
            continue
        
        full_path = os.path.join(base_dir, entry)
        if not os.path.isdir(full_path):
            continue

        # 匹配后缀
        matched_suffix = None
        for s in suffixes:
            if entry.endswith(s):
                matched_suffix = s
                break
        
        # 即使没有匹配到后缀，如果它是一个项目目录，我们也统计其内部的 dna_sra
        target_subdir = os.path.join(full_path, "00_data", "dna_sra")
        
        if os.path.exists(target_subdir):
            current_size_bytes = get_dir_size(target_subdir)
            current_size_gb = bytes_to_gb(current_size_bytes)
            
            # 打印单个目录统计结果
            print(f"{entry:<50} | {current_size_gb:>10.3f} GB")
            
            if matched_suffix:
                summary[matched_suffix] += current_size_bytes
        else:
            # 如果路径不存在则忽略
            pass

    print("-" * 70)
    print("各后缀汇总统计 (GB):")
    total_all_bytes = 0
    for suffix in suffixes:
        size_bytes = summary[suffix]
        total_all_bytes += size_bytes
        print(f"{suffix:<10}: {bytes_to_gb(size_bytes):>10.3f} GB")
    
    print("-" * 30)
    print(f"{'总计':<10}: {bytes_to_gb(total_all_bytes):>10.3f} GB")

if __name__ == "__main__":
    main()