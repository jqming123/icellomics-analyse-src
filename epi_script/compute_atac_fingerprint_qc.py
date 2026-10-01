#!/usr/bin/env python3
"""Run deepTools plotFingerprint and record diagnostic QC artifacts."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

from atac_qc.sample_map import group_replicates, read_replicate_records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--biosample", required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--blacklist-bed", type=Path)
    args = parser.parse_args()
    output_dir = args.project_dir / "1_result" / "5_qc" / "fingerprint" / args.biosample
    output_dir.mkdir(parents=True, exist_ok=True)
    status_json = output_dir / f"{args.biosample}.fingerprint.json"
    if shutil.which("plotFingerprint") is None:
        status_json.write_text(json.dumps({
            "status": "NOT_RUN", "warning": "plotFingerprint_not_found",
        }, indent=2) + "\n", encoding="utf-8")
        print("WARNING: plotFingerprint not found; fingerprint QC marked NOT_RUN")
        return 0
    groups = group_replicates(
        read_replicate_records(args.project_dir / "0_data" / "sample_run_map.tsv")
    ).get(args.biosample, {})
    if not groups:
        status_json.write_text(json.dumps({
            "status": "INCOMPLETE", "warning": "biosample_not_found_in_sample_map",
        }, indent=2) + "\n", encoding="utf-8")
        return 0
    bams = [
        args.project_dir / "1_result" / "1_alignment" / "replicates" / args.biosample / rep /
        f"{args.biosample}.{rep}.nodup.bam"
        for rep in groups
    ]
    missing = [str(path) for path in bams if not path.exists()]
    if missing:
        status_json.write_text(json.dumps({
            "status": "INCOMPLETE", "missing_bams": missing,
        }, indent=2) + "\n", encoding="utf-8")
        return 0
    metrics = output_dir / f"{args.biosample}.fingerprint.metrics.tsv"
    plot = output_dir / f"{args.biosample}.fingerprint.png"
    command = [
        "plotFingerprint", "--bamfiles", *map(str, bams), "--labels", *groups.keys(),
        "--plotFile", str(plot), "--outQualityMetrics", str(metrics),
        "--numberOfSamples", "50000", "--skipZeros", "--numberOfProcessors", str(args.threads),
    ]
    if len(bams) >= 2:
        command.extend(["--JSDsample", next(iter(groups))])
    if args.blacklist_bed and args.blacklist_bed.exists():
        command.extend(["--blackListFileName", str(args.blacklist_bed)])
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    payload = {
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "command": command,
        "metrics_file": str(metrics),
        "plot_file": str(plot),
        "stderr": result.stderr.strip(),
        "jsd_reference_replicate": next(iter(groups)) if len(bams) >= 2 else None,
        "interpretation": "Synthetic fingerprint/AUC metrics are sequencing-depth dependent; JS distance needs a declared reference sample.",
    }
    status_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if result.returncode:
        raise SystemExit(result.returncode)
    print(f"Wrote fingerprint QC: {status_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
