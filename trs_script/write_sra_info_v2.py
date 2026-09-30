#!/usr/bin/env python3
"""生成 sra_runid_prjid_ref.txt，并对原始数据进行完整性检查。"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


DEFAULT_RAWDATA_ROOT = (
    "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_rawdata"
)
DEFAULT_FASTQ_BASE_ROOT = (
    "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"
)
DEFAULT_OUTPUT_ROOT = (
    "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/EQ_jobs"
)
DEFAULT_PROJECT_INFO_FILE = "project_info.txt"


def usage_text(prog: str) -> str:
    return f"""用法 1，兼容旧版:
  python {prog} [参数] CELL_LINE REFERENCE_GENOME PROJECT_ID1 PROJECT_ID2 ...

例如:
  python {prog} CHO CH_Ensembl PRJEB30364 PRJNA941080

  python {prog} --skip-sra-check --check-fastq CHO CH_Ensembl PRJEB30364

用法 2，使用项目信息文件:
  python {prog} -i project_info.txt [参数]

参数:
  -i, --input FILE
      指定项目信息文件，文件格式为:
      PROJECT_ID<TAB>CELL_LINE<TAB>REFERENCE_GENOME

  --skip-sra-check
      跳过 .sra 文件数量与 Run ID 数量比较

  --check-fastq
      启用 .fastq.gz 文件存在性检查
      每个 Run 只要在 reads 目录下有至少一个 .fastq.gz 就算通过

  --rawdata-root DIR
      指定原始 SRA 数据根目录
      默认: {DEFAULT_RAWDATA_ROOT}

  --fastq-root DIR
      指定 fastq.gz 结果根目录
      默认: {DEFAULT_FASTQ_BASE_ROOT}

  --output-root DIR
      指定生成表格的输出根目录
      默认: {DEFAULT_OUTPUT_ROOT}

  -h, --help
      显示帮助信息

说明:
  旧版用法中:
    第 1 个参数为细胞系名称（CELL_LINE），例如 CHO
    第 2 个参数 REFERENCE_GENOME，例如 CH_Ensembl
    第 3 个及之后参数为 PROJECT_ID，例如 PRJEB30364

  输出目录名格式: PROJECT_ID_CELL_LINE，例如 PRJEB30364_CHO
  输出文件路径:
    <OUTPUT_ROOT>/PROJECT_ID_CELL_LINE/sra_runid_prjid_ref.txt

  输出表格格式（三列，制表符分隔）:
    SRR_ID<TAB>Project_ID<TAB>Reference_Name
    例如:
    ERR3001877<TAB>PRJEB30364_CHO<TAB>CH_Ensembl

  原始数据目录会优先查找:
    <RAWDATA_ROOT>/PROJECT_ID_CELL_LINE/sra_runid.txt

  例如:
    <RAWDATA_ROOT>/PRJEB30364_CHO/sra_runid.txt
"""


def make_composite_project_id(project_id: str, cell_line: str) -> str:
    """组合项目标识: PROJECT_ID_CELL_LINE，例如 PRJEB30364_CHO。"""
    return f"{project_id}_{cell_line}"


def get_rawdata_base_dir(
    rawdata_root: Path, cell_line: str, project_id: str
) -> Path:
    """自动选择原始数据目录。"""
    composite_id = make_composite_project_id(project_id, cell_line)
    candidates = [
        rawdata_root / composite_id,
        rawdata_root / cell_line / project_id,
        rawdata_root / project_id,
        rawdata_root / cell_line,
    ]

    for candidate in candidates:
        if (candidate / "sra_runid.txt").is_file():
            return candidate

    # 模糊匹配，例如 PRJNA899862_CHO、PRJNA899862_K562 等
    pattern = f"{project_id}_*"
    matched_dirs = sorted(rawdata_root.glob(pattern))
    for matched_dir in matched_dirs:
        if matched_dir.is_dir() and (matched_dir / "sra_runid.txt").is_file():
            return matched_dir

    # 如果都没找到，默认返回实际命名格式，方便后面报错
    return candidates[0]


def read_run_ids(sra_runid_file: Path) -> list[str]:
    run_ids: list[str] = []
    with sra_runid_file.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            first_col = line.split(None, 1)[0]
            if first_col.startswith("#"):
                continue
            run_ids.append(first_col)
    return run_ids


def count_sra_files(rawdata_base_dir: Path) -> int:
    return sum(1 for _ in rawdata_base_dir.rglob("*.sra"))


def count_fastq_in_reads(reads_dir: Path) -> int:
    if not reads_dir.is_dir():
        return 0
    return sum(
        1
        for path in reads_dir.iterdir()
        if path.is_file() and path.name.endswith(".fastq.gz")
    )


def process_project(
    project_id: str,
    cell_line: str,
    reference_genome: str,
    *,
    rawdata_root: Path,
    fastq_base_root: Path,
    output_root: Path,
    check_sra: bool,
    check_fastq: bool,
) -> None:
    composite_project_id = make_composite_project_id(project_id, cell_line)

    print()
    print("================================================")
    print(f"开始处理项目: {composite_project_id}")
    print(f"PROJECT_ID: {project_id}")
    print(f"CELL_LINE: {cell_line}")
    print(f"REFERENCE_GENOME: {reference_genome}")
    print("================================================")

    project_dir = output_root / composite_project_id
    project_dir.mkdir(parents=True, exist_ok=True)

    rawdata_base_dir = get_rawdata_base_dir(rawdata_root, cell_line, project_id)
    original_sra_runid_file = rawdata_base_dir / "sra_runid.txt"

    print("--- 开始数据完整性检查 ---")
    print(f"原始数据目录: {rawdata_base_dir}")
    print(f"Run ID 文件: {original_sra_runid_file}")

    if not rawdata_base_dir.is_dir():
        print(f"警告: 原始数据目录 '{rawdata_base_dir}' 不存在。跳过此项目。")
        return

    if not original_sra_runid_file.is_file():
        print(
            f"警告: 原始 sra_runid 文件 '{original_sra_runid_file}' 不存在。跳过此项目。"
        )
        return

    run_ids = read_run_ids(original_sra_runid_file)
    runid_count = len(run_ids)
    print(f"发现 sra_runid.txt 中的 Run ID 数量: {runid_count}")

    if check_sra:
        print()
        print(f"正在统计 '{rawdata_base_dir}' 及其子目录中的 .sra 文件...")
        sra_file_count = count_sra_files(rawdata_base_dir)
        print(f"发现 .sra 文件数量: {sra_file_count}")
        print(f"发现 Run ID 数量: {runid_count}")
        if sra_file_count == runid_count:
            print("数据完整性检查通过: .sra 文件数量与 Run ID 数量一致。")
        else:
            print("警告: 数据完整性检查失败! .sra 文件数量与 Run ID 数量不匹配。")
            print(
                f"请检查原始数据目录 '{rawdata_base_dir}' "
                f"和文件 '{original_sra_runid_file}'。"
            )
            print("继续处理，但请注意数据可能不完整。")
    else:
        print()
        print("已根据参数 --skip-sra-check 跳过 .sra 文件数量检查。")

    if check_fastq:
        print()
        print("正在检查每个 Run 是否至少存在一个 .fastq.gz 文件...")
        print("fastq.gz 文件位置格式:")
        print(f"{fastq_base_root}/<RUN>/reads/*.fastq.gz")

        runs_with_fastq_count = 0
        total_fastq_file_count = 0
        missing_reads_dir_runs: list[str] = []
        no_fastq_runs: list[str] = []

        for run_id in run_ids:
            reads_dir = fastq_base_root / run_id / "reads"
            if not reads_dir.is_dir():
                missing_reads_dir_runs.append(run_id)
                continue

            current_fastq_count = count_fastq_in_reads(reads_dir)
            total_fastq_file_count += current_fastq_count
            if current_fastq_count > 0:
                runs_with_fastq_count += 1
            else:
                no_fastq_runs.append(run_id)

        print(f"Run ID 总数: {runid_count}")
        print(f"至少有一个 .fastq.gz 的 Run 数量: {runs_with_fastq_count}")
        print(f"发现 .fastq.gz 文件总数: {total_fastq_file_count}")

        if runs_with_fastq_count == runid_count:
            print("数据完整性检查通过: 每个 Run 都至少存在一个 .fastq.gz 文件。")
        else:
            print("警告: 数据完整性检查失败! 存在没有 .fastq.gz 文件的 Run。")
            print("继续处理，但请注意数据可能不完整。")

        if missing_reads_dir_runs:
            print()
            print("以下 Run 缺少 reads 目录:")
            for run_id in missing_reads_dir_runs:
                print(run_id)

        if no_fastq_runs:
            print()
            print("以下 Run 的 reads 目录中没有 .fastq.gz 文件:")
            for run_id in no_fastq_runs:
                print(run_id)
    else:
        print()
        print("未启用 .fastq.gz 文件检查。如需启用，请添加参数 --check-fastq。")

    print("--- 数据完整性检查结束 ---")

    output_file = project_dir / "sra_runid_prjid_ref.txt"
    with output_file.open("w", encoding="utf-8", newline="\n") as handle:
        for run_id in run_ids:
            handle.write(f"{run_id}\t{composite_project_id}\t{reference_genome}\n")

    print()
    print("输出文件路径:")
    print(output_file)
    print()
    print("输出表格格式: SRR_ID\\tProject_ID\\tReference_Name")
    print(output_file.read_text(encoding="utf-8", errors="replace"), end="")

    line_count = len(run_ids)
    print()
    print("最终处理后的 sra_runid_prjid_ref.txt 文件行数:")
    print(line_count)

    print(f"--- 项目 {composite_project_id} 处理完成 ---")


def parse_project_info_file(project_info_file: Path) -> list[tuple[str, str, str]]:
    projects: list[tuple[str, str, str]] = []
    with project_info_file.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line.strip():
                continue

            parts = line.split("\t")
            project_id = parts[0].strip() if parts else ""
            cell_line = parts[1].strip() if len(parts) > 1 else ""
            reference_genome = parts[2].strip() if len(parts) > 2 else ""

            if not project_id:
                continue
            if project_id.startswith("#"):
                continue
            if project_id == "PROJECT_ID":
                continue
            if not cell_line or not reference_genome:
                print("警告: 项目信息不完整，跳过该行:")
                print(f"{project_id} {cell_line} {reference_genome}")
                continue

            projects.append((project_id, cell_line, reference_genome))
    return projects


def build_parser(prog: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description="生成 sra_runid_prjid_ref.txt，并对原始数据进行完整性检查。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=usage_text(prog).split("参数:", 1)[-1],
        add_help=False,
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="project_info_file",
        default=DEFAULT_PROJECT_INFO_FILE,
        help=f"指定项目信息文件（默认: {DEFAULT_PROJECT_INFO_FILE}）",
    )
    parser.add_argument(
        "--skip-sra-check",
        action="store_true",
        help="跳过 .sra 文件数量与 Run ID 数量比较",
    )
    parser.add_argument(
        "--check-fastq",
        action="store_true",
        help="启用 .fastq.gz 文件存在性检查",
    )
    parser.add_argument(
        "--rawdata-root",
        default=DEFAULT_RAWDATA_ROOT,
        help=f"指定原始 SRA 数据根目录（默认: {DEFAULT_RAWDATA_ROOT}）",
    )
    parser.add_argument(
        "--fastq-root",
        default=DEFAULT_FASTQ_BASE_ROOT,
        help=f"指定 fastq.gz 结果根目录（默认: {DEFAULT_FASTQ_BASE_ROOT}）",
    )
    parser.add_argument(
        "--output-root",
        default=DEFAULT_OUTPUT_ROOT,
        help=f"指定生成表格的输出根目录（默认: {DEFAULT_OUTPUT_ROOT}）",
    )
    parser.add_argument("-h", "--help", action="help", help="显示帮助信息")
    parser.add_argument(
        "positional_args",
        nargs="*",
        help="旧版模式: CELL_LINE REFERENCE_GENOME PROJECT_ID ...",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    prog = Path(argv[0]).name if argv else Path(__file__).name
    parser = build_parser(prog)
    args = parser.parse_args(argv[1:] if argv is not None else None)

    rawdata_root = Path(args.rawdata_root)
    fastq_base_root = Path(args.fastq_root)
    output_root = Path(args.output_root)
    check_sra = not args.skip_sra_check
    check_fastq = args.check_fastq

    print("================================================")
    print("开始处理所有项目")
    print("================================================")
    print(f"原始数据根目录: {rawdata_root}")
    print(f"fastq.gz 根目录: {fastq_base_root}")
    print(f"输出根目录: {output_root}")
    print(f".sra 检查: {check_sra}")
    print(f".fastq.gz 检查: {check_fastq}")
    print("================================================")

    positional_args = args.positional_args
    use_project_info_file = "-i" in (argv or sys.argv) or "--input" in (argv or sys.argv)

    if len(positional_args) >= 3:
        cell_line = positional_args[0]
        reference_genome = positional_args[1]
        project_ids = positional_args[2:]

        print("输入模式: 旧版位置参数模式")
        print(f"CELL_LINE: {cell_line}")
        print(f"REFERENCE_GENOME: {reference_genome}")
        print(f"PROJECT_ID 数量: {len(project_ids)}")

        for project_id in project_ids:
            process_project(
                project_id,
                cell_line,
                reference_genome,
                rawdata_root=rawdata_root,
                fastq_base_root=fastq_base_root,
                output_root=output_root,
                check_sra=check_sra,
                check_fastq=check_fastq,
            )
    else:
        project_info_file = Path(args.project_info_file)
        if not use_project_info_file and not project_info_file.is_file():
            print(
                f"错误: 未提供足够的位置参数，也没有找到默认项目信息文件 "
                f"'{project_info_file}'。"
            )
            print()
            print(usage_text(prog))
            return 1

        if not project_info_file.is_file():
            print(f"错误: 项目信息文件 '{project_info_file}' 不存在。")
            return 1

        print("输入模式: project_info 文件模式")
        print(f"项目信息文件: {project_info_file}")

        for project_id, cell_line, reference_genome in parse_project_info_file(
            project_info_file
        ):
            process_project(
                project_id,
                cell_line,
                reference_genome,
                rawdata_root=rawdata_root,
                fastq_base_root=fastq_base_root,
                output_root=output_root,
                check_sra=check_sra,
                check_fastq=check_fastq,
            )

    print()
    print("================================================")
    print("所有项目处理完成!")
    print("================================================")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
