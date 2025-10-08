#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch runner + parser for Shapley outputs.

INPUT FILE FORMAT (space/tab-separated), one model per line:
    /path/to/model1.bnet  TGT1  TGT2  TGT3
    /path/to/model2.bnet  P53   MYC   AKT

USAGE:
    python run_shapley_batch.py \
        --input models_and_targets.txt \
        --results results.csv \
        --logs logs \
        [--python-exe python]

What it does:
- For each (model, target):
    Runs:  {python-exe} -m Shapley -e {model} -o {target} -abpki
- Writes raw console output to: {logs}/{model_basename}__{target}.log
- Extracts metrics and appends a row to results.csv

The parser is tolerant to whitespace and fixed phrases as in your example.
If a metric is missing, it’s left blank in the CSV.
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Tuple, List

# --------- Regex patterns (robust to whitespace) ----------
RE_FLOAT = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
RE_INT = r"(\d+)"

PATTERNS = {
    "num_nodes": re.compile(rf"Number of nodes is\s+{RE_INT}", re.IGNORECASE),
    "error_ko": re.compile(rf"Error\s+KO:\s+{RE_FLOAT}", re.IGNORECASE),
    "error_ki": re.compile(rf"Error\s+KI:\s+{RE_FLOAT}", re.IGNORECASE),

    # Two-column metric lines under "Ranking Comparison"
    "kendall": re.compile(rf"Kendall\s+tau-b\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "spearman": re.compile(rf"Spearman\s+rho\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "mean_abs_rank_diff": re.compile(rf"Mean\s+abs\s+rank\s+diff\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "n_common": re.compile(rf"n_common\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),

    # Top-r Overlap
    "top10_jaccard": re.compile(rf"Top-10\s+Jaccard\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "p_at_10": re.compile(rf"Precision@10\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "r_at_10": re.compile(rf"Recall@10\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),

    "top20_jaccard": re.compile(rf"Top-20\s+Jaccard\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "p_at_20": re.compile(rf"Precision@20\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "r_at_20": re.compile(rf"Recall@20\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),

    # Top-group Overlap
    "group_jaccard_1": re.compile(rf"Top-1\s+groups\s+Jaccard\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "group_jaccard_2": re.compile(rf"Top-2\s+groups\s+Jaccard\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),
    "group_jaccard_3": re.compile(rf"Top-3\s+groups\s+Jaccard\s+{RE_FLOAT}\s+{RE_FLOAT}", re.IGNORECASE),

    # Times
    "time_heading": re.compile(rf"TIME\s+process\s+heading:\s+{RE_FLOAT}", re.IGNORECASE),
    "time_sim_orig": re.compile(rf"TIME\s+simulate\s+original\s+network:\s+{RE_FLOAT}", re.IGNORECASE),
    "time_orig_analysis": re.compile(rf"TIME\s+perform\s+original\s+analysis\s+process:\s+{RE_FLOAT}", re.IGNORECASE),
    "time_propagation": re.compile(rf"TIME\s+perform\s+propagation:\s+{RE_FLOAT}", re.IGNORECASE),
}

CSV_COLUMNS = [
    "timestamp",
    "model",
    "model_basename",
    "target",

    "num_nodes",
    "error_ko",
    "error_ki",

    "n_common_ko", "n_common_ki",
    "kendall_ko", "kendall_ki",
    "spearman_ko", "spearman_ki",
    "mean_abs_rank_diff_ko", "mean_abs_rank_diff_ki",

    "top10_jaccard_ko", "top10_jaccard_ki",
    "precision_at_10_ko", "precision_at_10_ki",
    "recall_at_10_ko", "recall_at_10_ki",

    "top20_jaccard_ko", "top20_jaccard_ki",
    "precision_at_20_ko", "precision_at_20_ki",
    "recall_at_20_ko", "recall_at_20_ki",

    "group_jaccard_1_ko", "group_jaccard_1_ki",
    "group_jaccard_2_ko", "group_jaccard_2_ki",
    "group_jaccard_3_ko", "group_jaccard_3_ki",

    "time_heading",
    "time_simulate_original",
    "time_original_analysis",
    "time_propagation",

    "returncode",
    "log_file",
]


def parse_line(line: str) -> Optional[Tuple[str, List[str]]]:
    """
    Parse one line of the input file: model_path target1 target2 target3
    Returns (model, [t1, t2, t3]) or None if invalid.
    """
    if not line.strip() or line.strip().startswith("#"):
        return None
    parts = line.strip().split()
    if len(parts) < 4:
        return None
    model = parts[0]
    targets = parts[1:4]
    return model, targets


def run_command(python_exe: str, model: str, target: str) -> Tuple[int, str]:
    """
    Run the Shapley command and capture stdout+stderr as text.
    """
    cmd = [python_exe, "-m", "Shapley", "-e", model, "-o", target, "-abpki"]
    try:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            check=False,
        )
        return proc.returncode, proc.stdout
    except FileNotFoundError:
        return 127, f"ERROR: Could not execute {' '.join(cmd)}. Is '{python_exe}' on PATH?\n"
    except Exception as exc:
        return 1, f"ERROR: Exception while running {' '.join(cmd)}:\n{exc}\n"


def _extract_pair(match: Optional[re.Match]) -> Tuple[Optional[str], Optional[str]]:
    if not match:
        return None, None
    return match.group(1), match.group(2)


def parse_metrics(text: str) -> Dict[str, Optional[str]]:
    """
    Extract all needed metrics from the command output.
    Returns a dict with string values (or None if missing).
    """
    d: Dict[str, Optional[str]] = {}

    # Singles
    d["num_nodes"] = (m.group(1) if (m := PATTERNS["num_nodes"].search(text)) else None)
    d["error_ko"] = (m.group(1) if (m := PATTERNS["error_ko"].search(text)) else None)
    d["error_ki"] = (m.group(1) if (m := PATTERNS["error_ki"].search(text)) else None)

    # Pairs (KO, KI)
    def put_pair(key_base: str, pat_key: str):
        m = PATTERNS[pat_key].search(text)
        ko, ki = _extract_pair(m)
        d[f"{key_base}_ko"] = ko
        d[f"{key_base}_ki"] = ki

    put_pair("n_common", "n_common")
    put_pair("kendall", "kendall")
    put_pair("spearman", "spearman")
    put_pair("mean_abs_rank_diff", "mean_abs_rank_diff")

    put_pair("top10_jaccard", "top10_jaccard")
    put_pair("precision_at_10", "p_at_10")
    put_pair("recall_at_10", "r_at_10")

    put_pair("top20_jaccard", "top20_jaccard")
    put_pair("precision_at_20", "p_at_20")
    put_pair("recall_at_20", "r_at_20")

    put_pair("group_jaccard_1", "group_jaccard_1")
    put_pair("group_jaccard_2", "group_jaccard_2")
    put_pair("group_jaccard_3", "group_jaccard_3")

    # Times
    d["time_heading"] = (m.group(1) if (m := PATTERNS["time_heading"].search(text)) else None)
    d["time_simulate_original"] = (m.group(1) if (m := PATTERNS["time_sim_orig"].search(text)) else None)
    d["time_original_analysis"] = (m.group(1) if (m := PATTERNS["time_orig_analysis"].search(text)) else None)
    d["time_propagation"] = (m.group(1) if (m := PATTERNS["time_propagation"].search(text)) else None)

    return d


def ensure_parent(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)


def append_csv_row(csv_path: Path, row: Dict[str, Optional[str]]):
    file_exists = csv_path.exists()
    ensure_parent(csv_path)
    with csv_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Path to models+targets text file")
    ap.add_argument("--results", required=True, help="Path to output CSV (appends)")
    ap.add_argument("--logs", default="logs", help="Directory to write per-run logs")
    ap.add_argument("--python-exe", default=sys.executable, help="Python executable to run Shapley (default: current python)")
    args = ap.parse_args()

    in_path = Path(args.input)
    results_path = Path(args.results)
    logs_dir = Path(args.logs)
    logs_dir.mkdir(parents=True, exist_ok=True)

    with in_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()

    jobs = []
    for ln in lines:
        parsed = parse_line(ln)
        if not parsed:
            continue
        model, targets = parsed
        for tgt in targets:
            jobs.append((model, tgt))

    if not jobs:
        print("No valid jobs found in input.", file=sys.stderr)
        sys.exit(2)

    for model, target in jobs:
        model_path = Path(model)
        model_base = model_path.stem
        safe_target = re.sub(r"[^\w\-.]+", "_", target.strip())
        log_file = logs_dir / f"{model_base}__{safe_target}.log"

        rc, out_text = run_command(args.python_exe, str(model_path), target)

        # Write raw log
        ensure_parent(log_file)
        with log_file.open("w", encoding="utf-8") as lf:
            lf.write(out_text)

        # Parse metrics
        metrics = parse_metrics(out_text)

        # Build CSV row
        row = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "model": str(model_path),
            "model_basename": model_base,
            "target": target,
            "returncode": rc,
            "log_file": str(log_file),
        }

        # Fill parsed fields
        for col in CSV_COLUMNS:
            if col in row:
                continue
            row[col] = metrics.get(col)

        append_csv_row(results_path, row)

        # Small console feedback
        ok = "OK" if rc == 0 else f"RC={rc}"
        print(f"[{ok}] {model_base} :: {target}  -> {log_file}")

    print(f"\nDone. Parsed results appended to: {results_path}")
    print(f"Logs saved under: {logs_dir.resolve()}")


if __name__ == "__main__":
    main()
