#!/usr/bin/env python3
"""
基于 agent_DE_config.csv + agent_full_report.md 为 MRC-5 / WI-38 生成 DEG 分组配置。
v3: 解析 MD 文件提取细胞系信息，与 CSV 关键字匹配交叉验证，全面覆盖所有项目。
"""

import csv
import os
import re
import sys

BASE_DIR = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects"
GROUPS_BASE = os.path.join(BASE_DIR, "DEG", "groups")
CSV_PATH = "/hpcdisk1/zhaowm_group/gaoxiaojing/agent_DE_config.csv"
MD_PATH = "/hpcdisk1/zhaowm_group/gaoxiaojing/agent_full_report.md"

MRC5_KEYWORDS = ["MRC5", "MRC-5", "mrc5", "iCAF", "iCAFs"]
WI38_KEYWORDS = ["WI38", "WI-38", "wi38", "WI_38"]

NON_MRC5_CELL_LINES = [
    "A549", "IMR90", "IMR-90", "HEK293", "HEK-293", "293T",
    "U2OS", "HCT116", "Hela", "HeLa", "MCF7", "MCF-7",
    "mouse", "Mouse", "MEF", "mef", "hippocampus", "PBMC",
    "iPSC", "H9", "Detroit", "CCD-18Co", "CS1AN", "GM",
    "W-V", "WI38VA13", "VA13", "Huh7", "HepG2", "BJ",
    "AHLM", "SCLC", "DMS153", "H209",
]

NON_WI38_CELL_LINES = [
    "A549", "IMR90", "IMR-90", "HEK293", "HEK-293", "293T",
    "U2OS", "HCT116", "Hela", "HeLa", "MCF7", "MCF-7",
    "mouse", "Mouse", "MEF", "mef", "MRC5", "MRC-5", "mrc5",
    "iPSC", "H9", "Detroit", "CCD-18Co", "BJ",
    "WI38VA13", "VA13", "W-V",
    "MDAMB436", "SKBR3", "ZR751", "HMEC", "THP1",
    "NCI-H1299", "H1299",
]

SKIP_PATTERNS = [
    r"质量=(weak|moderate|unknown)",
    r"需人工确认",
    r"未生成对比",
    r"^N/A",
    r"无有效分组",
    r"\(无标记实验组\)",
    r"\(无标记对照组\)",
]


def parse_md_cell_lines(md_path):
    """Parse agent_full_report.md to extract bioproject -> cell_line mapping."""
    mapping = {}
    if not os.path.exists(md_path):
        print(f"WARNING: MD file not found: {md_path}")
        return mapping

    with open(md_path, "r") as f:
        content = f.read()

    sections = re.split(r'\n## (PRJ[EN]A?\d+)', content)
    for i in range(1, len(sections), 2):
        bioproject = sections[i].strip()
        section_text = sections[i + 1] if i + 1 < len(sections) else ""

        desc_match = re.search(r'\*\*项目描述\*\*:\s*(.+?)(?:\n|$)', section_text)
        if not desc_match:
            continue
        description = desc_match.group(1)

        is_mrc5 = has_keyword(description, MRC5_KEYWORDS)
        is_wi38 = has_keyword(description, WI38_KEYWORDS)

        if is_mrc5 and not is_wi38:
            mapping[bioproject] = "MRC-5"
        elif is_wi38 and not is_mrc5:
            mapping[bioproject] = "WI-38"
        elif is_mrc5 and is_wi38:
            pass

    return mapping


def is_valid_entry(row):
    exp = row.get("exp_group", "").strip()
    ctrl = row.get("ctrl_group", "").strip()
    if not exp or not ctrl:
        return False
    if not parse_sample_ids(exp) or not parse_sample_ids(ctrl):
        return False
    combined = exp + ctrl + row.get("analysis_content", "")
    for pat in SKIP_PATTERNS:
        if re.search(pat, combined):
            return False
    return True


def parse_sample_ids(text):
    return re.findall(r'[SED]RR\d+', text)


def has_keyword(text, keywords):
    tokens = set(re.split(r'[_\s,;:()\[\]{}]+', text))
    for kw in keywords:
        if kw in tokens:
            return True
    return False


def has_non_target_cell_line(text, non_target_list):
    tokens = set(re.split(r'[_\s,;:()\[\]{}]+', text))
    for cl in non_target_list:
        if cl in tokens:
            return True
        if re.search(r'(?<![a-zA-Z])' + re.escape(cl) + r'(?![a-zA-Z])', text):
            return True
    return False


def classify_entry(row, md_cell_mapping):
    """
    Returns: "MRC-5", "WI-38", or None
    Uses CSV keyword matching as primary, MD file mapping as supplementary.
    """
    bioproject = row.get("bioproject", "").strip()
    folder_name = row.get("folder_name", "").strip()
    analysis_content = row.get("analysis_content", "").strip()
    exp_group = row.get("exp_group", "").strip()
    ctrl_group = row.get("ctrl_group", "").strip()
    exp_vs_ctrl = row.get("exp_vs_ctrl", "").strip()

    full_text = f"{folder_name} {analysis_content} {exp_group} {ctrl_group} {exp_vs_ctrl}"

    is_mrc5_kw = has_keyword(full_text, MRC5_KEYWORDS)
    is_wi38_kw = has_keyword(full_text, WI38_KEYWORDS)

    md_cell = md_cell_mapping.get(bioproject)

    if is_mrc5_kw and is_wi38_kw:
        if md_cell == "MRC-5":
            return "MRC-5"
        elif md_cell == "WI-38":
            return "WI-38"
        return None

    if is_mrc5_kw:
        if has_non_target_cell_line(full_text, NON_MRC5_CELL_LINES):
            mrc5_specific = has_keyword(exp_vs_ctrl, MRC5_KEYWORDS) or \
                            has_keyword(folder_name, MRC5_KEYWORDS)
            if not mrc5_specific:
                return None
        if has_non_target_cell_line(ctrl_group, NON_MRC5_CELL_LINES):
            return None
        if has_non_target_cell_line(exp_group, NON_MRC5_CELL_LINES):
            return None
        return "MRC-5"

    if is_wi38_kw:
        if has_non_target_cell_line(full_text, NON_WI38_CELL_LINES):
            wi38_specific = has_keyword(exp_vs_ctrl, WI38_KEYWORDS) or \
                            has_keyword(folder_name, WI38_KEYWORDS)
            if not wi38_specific:
                return None
        if has_non_target_cell_line(ctrl_group, NON_WI38_CELL_LINES):
            return None
        if has_non_target_cell_line(exp_group, NON_WI38_CELL_LINES):
            return None
        return "WI-38"

    if md_cell == "MRC-5":
        if has_non_target_cell_line(full_text, NON_MRC5_CELL_LINES):
            return None
        if has_non_target_cell_line(ctrl_group, NON_MRC5_CELL_LINES):
            return None
        if has_non_target_cell_line(exp_group, NON_MRC5_CELL_LINES):
            return None
        return "MRC-5"

    if md_cell == "WI-38":
        if has_non_target_cell_line(full_text, NON_WI38_CELL_LINES):
            return None
        if has_non_target_cell_line(ctrl_group, NON_WI38_CELL_LINES):
            return None
        if has_non_target_cell_line(exp_group, NON_WI38_CELL_LINES):
            return None
        return "WI-38"

    return None


def sanitize_folder_name(name):
    return name.replace("/", "_").replace("\\", "_")


def generate_group_files(cell_line, folder_name, exp_group_str, ctrl_group_str):
    safe_folder = sanitize_folder_name(folder_name)
    group_dir = os.path.join(GROUPS_BASE, cell_line, safe_folder)
    os.makedirs(group_dir, exist_ok=True)

    exp_ids = parse_sample_ids(exp_group_str)
    ctrl_ids = parse_sample_ids(ctrl_group_str)

    if not exp_ids or not ctrl_ids:
        return [], []

    tsv_path = os.path.join(group_dir, f"{safe_folder}.tsv")
    list_path = os.path.join(group_dir, f"{safe_folder}.list")

    with open(tsv_path, "w") as f:
        f.write("sample\tcondition\n")
        for sid in exp_ids:
            f.write(f"{sid}\tcase\n")
        for sid in ctrl_ids:
            f.write(f"{sid}\tcontrol\n")

    with open(list_path, "w") as f:
        for sid in exp_ids:
            f.write(f"{sid}\n")
        for sid in ctrl_ids:
            f.write(f"{sid}\n")

    return exp_ids, ctrl_ids


def main():
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: CSV not found: {CSV_PATH}")
        sys.exit(1)

    md_cell_mapping = parse_md_cell_lines(MD_PATH)
    print(f"Parsed {len(md_cell_mapping)} projects from MD file")
    mrc5_from_md = sum(1 for v in md_cell_mapping.values() if v == "MRC-5")
    wi38_from_md = sum(1 for v in md_cell_mapping.values() if v == "WI-38")
    print(f"  MRC-5 from MD: {mrc5_from_md}, WI-38 from MD: {wi38_from_md}")

    all_configs = []
    skipped_count = 0
    csv_mrc5_count = 0
    csv_wi38_count = 0
    md_only_mrc5_count = 0
    md_only_wi38_count = 0

    with open(CSV_PATH, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            bioproject = row.get("bioproject", "").strip()
            folder_name = row.get("folder_name", "").strip()
            analysis_content = row.get("analysis_content", "").strip()
            exp_vs_ctrl = row.get("exp_vs_ctrl", "").strip()
            exp_group = row.get("exp_group", "").strip()
            ctrl_group = row.get("ctrl_group", "").strip()

            if not is_valid_entry(row):
                skipped_count += 1
                continue

            cell_line = classify_entry(row, md_cell_mapping)
            if cell_line is None:
                skipped_count += 1
                continue

            full_text = f"{folder_name} {analysis_content} {exp_group} {ctrl_group} {exp_vs_ctrl}"
            is_csv_mrc5 = has_keyword(full_text, MRC5_KEYWORDS)
            is_csv_wi38 = has_keyword(full_text, WI38_KEYWORDS)

            if cell_line == "MRC-5":
                if is_csv_mrc5:
                    csv_mrc5_count += 1
                else:
                    md_only_mrc5_count += 1
            elif cell_line == "WI-38":
                if is_csv_wi38:
                    csv_wi38_count += 1
                else:
                    md_only_wi38_count += 1

            exp_ids, ctrl_ids = generate_group_files(
                cell_line, folder_name, exp_group, ctrl_group
            )

            if not exp_ids or not ctrl_ids:
                skipped_count += 1
                continue

            all_configs.append({
                "cell_line": cell_line,
                "bioproject": bioproject,
                "folder_name": folder_name,
                "analysis_content": analysis_content,
                "exp_vs_ctrl": exp_vs_ctrl,
                "exp_group": exp_group,
                "ctrl_group": ctrl_group,
                "exp_ids": exp_ids,
                "ctrl_ids": ctrl_ids,
            })

    report_path = os.path.join(BASE_DIR, "DEG", "DE_config_summary.md")
    with open(report_path, "w") as f:
        f.write("# DEG 差异分析配置汇总\n\n")
        f.write(f"生成时间: 2026-05-08 | 基于 agent_DE_config.csv + agent_full_report.md 交叉验证\n\n")
        f.write(f"总计: {len(all_configs)} 个对比 | 跳过无效条目: {skipped_count}\n")
        f.write(f"来源: CSV关键字匹配 MRC-5={csv_mrc5_count}, WI-38={csv_wi38_count} | ")
        f.write(f"MD补充 MRC-5={md_only_mrc5_count}, WI-38={md_only_wi38_count}\n\n")

        for cl in ["MRC-5", "WI-38"]:
            cl_configs = [c for c in all_configs if c["cell_line"] == cl]
            f.write(f"## {cl} ({len(cl_configs)} 个对比)\n\n")
            f.write("| bioproject | 差异分析文件夹名称 | 差异分析的内容 | 实验组 vs 对照组 | 实验组 | 对照组 |\n")
            f.write("|------------|-------------------|---------------|-----------------|--------|--------|\n")
            for cfg in cl_configs:
                exp_label = cfg['exp_vs_ctrl'].split(' vs ')[0] if ' vs ' in cfg['exp_vs_ctrl'] else cfg['exp_vs_ctrl']
                ctrl_label = cfg['exp_vs_ctrl'].split(' vs ')[-1] if ' vs ' in cfg['exp_vs_ctrl'] else ''
                exp_str = f"{exp_label} ({', '.join(cfg['exp_ids'])})"
                ctrl_str = f"{ctrl_label} ({', '.join(cfg['ctrl_ids'])})"
                f.write(f"| {cfg['bioproject']} | {cfg['folder_name']} | {cfg['analysis_content']} | {cfg['exp_vs_ctrl']} | {exp_str} | {ctrl_str} |\n")
            f.write("\n")

    print(f"\nGenerated {len(all_configs)} DEG configs")
    print(f"Skipped {skipped_count} invalid entries")
    print(f"CSV keyword: MRC-5={csv_mrc5_count}, WI-38={csv_wi38_count}")
    print(f"MD supplement: MRC-5={md_only_mrc5_count}, WI-38={md_only_wi38_count}")
    print(f"Report: {report_path}")

    for cl in ["MRC-5", "WI-38"]:
        count = len([c for c in all_configs if c["cell_line"] == cl])
        cl_dir = os.path.join(GROUPS_BASE, cl)
        proj_count = len(os.listdir(cl_dir)) if os.path.exists(cl_dir) else 0
        print(f"  {cl}: {count} comparisons, {proj_count} group folders")


if __name__ == "__main__":
    main()
