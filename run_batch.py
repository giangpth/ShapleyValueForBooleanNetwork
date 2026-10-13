#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path

# ---------------- CONFIG (hard-coded) ----------------
INPUT_LIST = "realbenchmark.txt"            # TAB-separated: path<TAB>out1<TAB>out2<TAB>out3
LOG_DIR    = Path("logs/realbenchmark_Prop_sample_60_3rdRun")  # where the terminal logs are saved 
PROGRAM    = [sys.executable, "-m", "Shapley"]
ECHO_TO_TERMINAL = True                  # also print the log live

# Fixed parameters (edit freely; remove any you don't want)
FIXED_ARGS = [
    "-r", "uniform", # Sampling strategy: all, uniform, antithetic, stratified
    "-p", "60.0", # Sample percentage (0-100)
    "-c", "0.0", # Percentage of nodes to be simulated (0-100)
    "-k", "Shapley", # Weighting scheme: Shapley, Uniform
    "-m", "Propagation", # Method: Propagation, Simulation
    "-t" # Enable Test mode to compare with fully simulation results
]
# -----------------------------------------------------

def run_one(name, expr_path, output, log_path):
    cmd = PROGRAM + ["-n", name, "-e", expr_path, "-o", output] + FIXED_ARGS
    with open(log_path, "w", encoding="utf-8") as log:
        log.write("CMD: " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        for out_line in proc.stdout:
            log.write(out_line)
            if ECHO_TO_TERMINAL:
                print(out_line, end="")
        proc.wait()
        log.write(f"\n[exit code: {proc.returncode}]\n")
    return proc.returncode

def main():
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    with open(INPUT_LIST, encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]

    jobs = []
    for i, line in enumerate(lines, 1):
        parts = [p.strip() for p in line.split("\t") if p.strip()]
        if len(parts) < 2:
            print(f"[line {i}] Skipping malformed line: {line!r}")
            continue
        expr_path, outputs = parts[0], parts[1:]
        stem = Path(expr_path).stem                  # "1" from data/realbenchmark/1.txt
        for out in outputs:                          # one run per output
            jobs.append((f"model_{stem}_{out}", expr_path, out))

    failed = []
    for n, (name, expr_path, out) in enumerate(jobs, 1):
        log_path = LOG_DIR / f"{name}.log"
        print(f"[{n}/{len(jobs)}] {name}")
        rc = run_one(name, expr_path, out, log_path)
        if rc != 0:
            failed.append(name)
            print(f"  -> FAILED (exit {rc}), see {log_path}")

    print(f"\nDone. {len(jobs) - len(failed)} succeeded, {len(failed)} failed.")
    for name in failed:
        print("  failed:", name)

if __name__ == "__main__":
    main()