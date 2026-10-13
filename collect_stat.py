#!/usr/bin/env python3
"""
Collect metrics from Shapley batch log files across several runs.

Usage:
    python collect_stats.py run1_dir run2_dir run3_dir -o report.csv

Outputs:
    report.csv          one row per (run, model, target)
    report_summary.csv  one row per (model, target) with mean/std over runs
"""
import argparse
import csv
import math
import re
import statistics
import sys
from pathlib import Path

NUM = r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|NaN|inf|-inf)"

# column name -> regex (one capture group)
SINGLE = {
    "Relative_RMSE_KO":               rf"Relative_RMSE_KO_is\s+{NUM}",
    "Relative_RMSE_KI":               rf"Relative_RMSE_KI_is\s+{NUM}",
    "NDCG_KO":                        rf"NDCG_KO_is\s+{NUM}",
    "NDCG_KI":                        rf"NDCG_KI_is\s+{NUM}",
    "AVERAGE_CORRECT_ROWS_ALL_NODES": rf"AVERAGE_CORRECT_ROWS_ALL_NODES:\s+{NUM}",
    "NUMBER_OF_INPUT":                rf"NUMBER OF INPUT:\s+{NUM}",
    "NUMBER_OF_NODES":                rf"NUMBER OF NODES:\s+{NUM}",
    "NUMBER_OF_EDGES":                rf"NUMBER OF EDGES:\s+{NUM}",
    "NUMBER_OF_SCCs":                 rf"NUMBER OF SCCs:\s+{NUM}",
    "NUMBER_OF_aFAS":                 rf"NUMBER OF aFAS:\s+{NUM}",
    "TIME_total_Propagation":         rf"TIME total for Propagation:\s+{NUM}",
    # "TIME_total_Simulation":         rf"TIME total for Simulation:\s+{NUM}",
    "TIME_total_FULLY_Simulation":    rf"TIME total for FULLY Simulation:\s+{NUM}",
    "Avg_error_KO":                   rf"Average error for KO values:\s+{NUM}",
    "Avg_error_KI":                   rf"Average error for KI values:\s+{NUM}",
}

# table rows with two values (KO, KI): base column name -> row label regex
PAIRED = {
    "Kendall_tau_b":     r"Kendall tau-b",
    "Spearman_rho":      r"Spearman rho",
    "Top10_Jaccard":     r"Top-10 Jaccard",
    "Top1_groups_Jaccard": r"Top-1\s+groups Jaccard",
}

SINGLE_RE = {k: re.compile(v) for k, v in SINGLE.items()}
PAIRED_RE = {k: re.compile(rf"^\s*{v}\s+{NUM}\s+{NUM}\s*$", re.M) for k, v in PAIRED.items()}

COLUMNS = []
for k in SINGLE:
    COLUMNS.append(k)
    if k == "NDCG_KI":
        pass
for k in PAIRED:
    COLUMNS += [f"{k}_KO", f"{k}_KI"]
# put in a sensible order
ORDER = [
    "Relative_RMSE_KO", "Relative_RMSE_KI", "NDCG_KO", "NDCG_KI",
    "AVERAGE_CORRECT_ROWS_ALL_NODES",
    "Kendall_tau_b_KO", "Kendall_tau_b_KI",
    "Spearman_rho_KO", "Spearman_rho_KI",
    "Top10_Jaccard_KO", "Top10_Jaccard_KI",
    "Top1_groups_Jaccard_KO", "Top1_groups_Jaccard_KI",
    "NUMBER_OF_INPUT", "NUMBER_OF_NODES", "NUMBER_OF_EDGES",
    "NUMBER_OF_SCCs", "NUMBER_OF_aFAS",
    "TIME_total_Propagation", "TIME_total_FULLY_Simulation",
    "Avg_error_KO", "Avg_error_KI",
]
COLUMNS = ORDER

NAME_RE = re.compile(r"^model_(.+?)_(.+)\.log$")


def parse_log(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    vals = {}
    for key, rx in SINGLE_RE.items():
        m = rx.search(text)
        vals[key] = m.group(1) if m else ""
    for key, rx in PAIRED_RE.items():
        m = rx.search(text)
        vals[f"{key}_KO"] = m.group(1) if m else ""
        vals[f"{key}_KI"] = m.group(2) if m else ""
    return vals


def to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folders", nargs="+", help="folders of log files, one per run")
    ap.add_argument("-o", "--output", default="report.csv", help="output CSV path")
    args = ap.parse_args()

    rows = []
    for folder in args.folders:
        fdir = Path(folder)
        if not fdir.is_dir():
            sys.exit(f"Not a directory: {fdir}")
        run = fdir.name
        for log in sorted(fdir.glob("*.log")):
            m = NAME_RE.match(log.name)
            if not m:
                print(f"Skipping unrecognised file name: {log.name}")
                continue
            model, target = m.group(1), m.group(2)
            vals = parse_log(log)

            missing = [k for k in COLUMNS if vals[k] == ""]
            if missing:
                print(f"[{run}/{log.name}] missing: {', '.join(missing)}")
            rows.append({"run": run, "model": model, "target": target, **vals})

    def sort_key(r):
        mk = (0, int(r["model"])) if r["model"].isdigit() else (1, r["model"])
        return (mk, r["target"], r["run"])

    rows.sort(key=sort_key)

    out = Path(args.output)
    header = ["run", "model", "target"] + COLUMNS
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out}")

    # summary: mean / std across runs for each (model, target)
    groups = {}
    for r in rows:
        groups.setdefault((r["model"], r["target"]), []).append(r)

    summary_path = out.with_name(out.stem + "_summary" + out.suffix)
    sum_header = ["model", "target", "n_runs"]
    for c in COLUMNS:
        sum_header += [f"{c}_mean", f"{c}_std"]
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(sum_header)
        for (model, target), grp in groups.items():
            line = [model, target, len(grp)]
            for c in COLUMNS:
                xs = [x for x in (to_float(r[c]) for r in grp) if x is not None]
                clean_xs = [x for x in xs if math.isfinite(x)]
                if clean_xs:
                    line.append(statistics.mean(clean_xs))
                    line.append(statistics.stdev(clean_xs) if len(clean_xs) > 1 else 0.0)
                else:
                    line += ["", ""]
            w.writerow(line)
    print(f"Wrote {len(groups)} (model, target) summaries to {summary_path}")


if __name__ == "__main__":
    main()