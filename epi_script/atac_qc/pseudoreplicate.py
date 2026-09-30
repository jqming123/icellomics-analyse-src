"""Deterministic streaming split of tagAlign records into pseudoreplicates."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import Iterator, TextIO


PSEUDOREP_METHOD = "sha256_balanced_stream_split_v2"


def _groups(handle: TextIO, layout: str) -> Iterator[list[str]]:
    if layout == "PE":
        while True:
            first = handle.readline()
            if not first:
                return
            second = handle.readline()
            if not second:
                raise ValueError("PE tagAlign contains an odd number of records")
            yield [first, second]
    else:
        for line in handle:
            yield [line]


def split_pseudoreplicates(
    input_file: Path,
    output_one: Path,
    output_two: Path,
    layout: str,
    seed: int,
) -> tuple[int, int]:
    output_one.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(input_file, "rt", encoding="utf-8") as count_handle:
        line_count = sum(1 for _ in count_handle)
    if layout == "PE" and line_count % 2:
        raise ValueError("PE tagAlign contains an odd number of records")
    group_count = line_count // 2 if layout == "PE" else line_count
    targets = [(group_count + 1) // 2, group_count // 2]
    counts = [0, 0]
    with gzip.open(input_file, "rt", encoding="utf-8") as source, gzip.open(
        output_one, "wt", encoding="utf-8"
    ) as first_out, gzip.open(output_two, "wt", encoding="utf-8") as second_out:
        outputs = (first_out, second_out)
        for index, lines in enumerate(_groups(source, layout)):
            digest = hashlib.sha256(
                f"{seed}\t{index}\t".encode("utf-8") + "".join(lines).encode("utf-8")
            ).digest()
            target = digest[0] & 1
            if counts[target] >= targets[target]:
                target = 1 - target
            outputs[target].writelines(lines)
            counts[target] += 1
    if min(counts) == 0:
        raise ValueError("pseudoreplicate split produced an empty output")
    return counts[0], counts[1]
