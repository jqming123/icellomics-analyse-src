"""ENCODE-style NRF/PBC library-complexity calculation."""

from __future__ import annotations

import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Iterable


def summarize_multiplicities(multiplicities: Iterable[int]) -> dict[str, float | int | None]:
    values = [value for value in multiplicities if value > 0]
    total = sum(values)
    distinct = len(values)
    once = sum(value == 1 for value in values)
    twice = sum(value == 2 for value in values)
    return {
        "total_fragments": total,
        "distinct_fragments": distinct,
        "one_read_fragments": once,
        "two_read_fragments": twice,
        "nrf": distinct / total if total else None,
        "pbc1": once / distinct if distinct else None,
        "pbc2": once / twice if twice else None,
    }


def compute_library_complexity(
    bam_file: Path,
    layout: str,
    output_file: Path,
    threads: int = 2,
) -> dict[str, object]:
    empty: dict[str, object] = {
        "library_total_fragments": None,
        "library_distinct_fragments": None,
        "library_one_read_fragments": None,
        "library_two_read_fragments": None,
        "nrf": None,
        "pbc1": None,
        "pbc2": None,
        "library_complexity_file": "",
        "warnings": [],
    }
    if not bam_file.exists():
        empty["warnings"] = [f"library_complexity_bam_missing:{bam_file}"]
        return empty
    required = ("samtools", "bedtools", "sort", "uniq")
    missing_tools = [tool for tool in required if shutil.which(tool) is None]
    if missing_tools:
        empty["warnings"] = ["library_complexity_tools_missing:" + ",".join(missing_tools)]
        return empty

    quoted_bam = shlex.quote(str(bam_file))
    if layout == "PE":
        interval_command = (
            f"samtools sort -@ {threads} -n -O BAM {quoted_bam} | "
            "bedtools bamtobed -bedpe -i stdin | "
            "awk 'BEGIN{OFS=\"\\t\"} $1==$4 {print $1,$2,$4,$6,$9,$10}'"
        )
    elif layout == "SE":
        interval_command = (
            f"bedtools bamtobed -i {quoted_bam} | "
            "awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$6}'"
        )
    else:
        empty["warnings"] = [f"library_complexity_layout_unsupported:{layout}"]
        return empty

    command = (
        f"set -o pipefail; {interval_command} | LC_ALL=C sort | uniq -c | "
        "awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} "
        "$1==1{m1++} $1==2{m2++} {m0++;mt+=$1} "
        "END{nrf=(mt?m0/mt:-1);pbc1=(m0?m1/m0:-1);pbc2=(m2?m1/m2:-1); "
        "printf \"%d\\t%d\\t%d\\t%d\\t%.10f\\t%.10f\\t%.10f\\n\",mt,m0,m1,m2,nrf,pbc1,pbc2}'"
    )
    result = subprocess.run(
        command,
        shell=True,
        executable="/bin/bash",
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        empty["warnings"] = [f"library_complexity_failed:{result.stderr.strip() or 'pipeline failed'}"]
        return empty

    fields = result.stdout.strip().split("\t")
    if len(fields) != 7:
        empty["warnings"] = ["library_complexity_invalid_output"]
        return empty
    try:
        values: dict[str, float | int | None] = {
            "total_fragments": int(fields[0]), "distinct_fragments": int(fields[1]),
            "one_read_fragments": int(fields[2]), "two_read_fragments": int(fields[3]),
            "nrf": float(fields[4]) if float(fields[4]) >= 0 else None,
            "pbc1": float(fields[5]) if float(fields[5]) >= 0 else None,
            "pbc2": float(fields[6]) if float(fields[6]) >= 0 else None,
        }
    except ValueError:
        empty["warnings"] = ["library_complexity_invalid_output"]
        return empty
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with output_file.open("w", encoding="utf-8", newline="") as handle:
        handle.write("total\tdistinct\tonce\ttwice\tnrf\tpbc1\tpbc2\n")
        handle.write(
            "{total_fragments}\t{distinct_fragments}\t{one_read_fragments}\t"
            "{two_read_fragments}\t{nrf}\t{pbc1}\t{pbc2}\n".format(**values)
        )
    return {
        "library_total_fragments": values["total_fragments"],
        "library_distinct_fragments": values["distinct_fragments"],
        "library_one_read_fragments": values["one_read_fragments"],
        "library_two_read_fragments": values["two_read_fragments"],
        "nrf": values["nrf"],
        "pbc1": values["pbc1"],
        "pbc2": values["pbc2"],
        "library_complexity_file": str(output_file),
        "warnings": [],
    }
