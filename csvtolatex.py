#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate LaTeX tables from two benchmark CSVs (acyclic & normal).

OUTPUT FILES (exactly four):
  - {outdir}/acyclic_performance.tex     (rotated longtable; may contain 1..N parts)
  - {outdir}/acyclic_time.tex
  - {outdir}/normal_performance.tex      (rotated longtable; may contain 1..N parts)
  - {outdir}/normal_time.tex

Performance tables:
  - Landscape pages (lscape), full-bleed margins (geometry), multipage (longtable)
  - Split into N parts inside the SAME .tex (default: 2), each with its own caption/label
  - First columns: Model, Target; every metric shown as KO/KI pair.
  - NOW includes Top-20 metrics (Top-20 Jaccard, Precision@20, Recall@20).

Time tables:
  - Regular table float (portrait), per-(model,target) rows, includes Total column.

USAGE:
  python make_tables_two_csvs_fullbleed.py \
      --acyclic acyclic.csv \
      --normal  normal.csv \
      --outdir  tables \
      --perf-split 2

LaTeX preamble (required):
  \\usepackage{lscape}
  \\usepackage{longtable}
  \\usepackage{booktabs}
  \\usepackage{geometry}
  \\usepackage{caption}
  \\captionsetup[longtable]{aboveskip=2pt, belowskip=2pt}

CSV expectations:
  - Identifier columns: preferred "model_basename" (fallback "model"), and "target".
  - Performance metrics present as KO/KI pairs (any subset is fine):
      error_ko, error_ki
      kendall_ko, kendall_ki
      spearman_ko, spearman_ki
      mean_abs_rank_diff_ko, mean_abs_rank_diff_ki
      top10_jaccard_ko, top10_jaccard_ki
      precision_at_10_ko, precision_at_10_ki
      recall_at_10_ko, recall_at_10_ki
      top20_jaccard_ko, top20_jaccard_ki
      precision_at_20_ko, precision_at_20_ki
      recall_at_20_ko, recall_at_20_ki
  - Time metrics (any subset is fine):
      time_simulate_original, time_original_analysis, time_propagation, time_heading (optional)
"""

from pathlib import Path
import argparse
import pandas as pd

# ---------------- Config ----------------

# Performance metric groups: (display name, base key) => uses base_ko / base_ki
PERF_COLS = [
    ("Error", "error"),
    ("Kendall $\\tau$-b", "kendall"),
    ("Spearman $\\rho$", "spearman"),
    ("Mean abs rank diff", "mean_abs_rank_diff"),
    ("Top-10 Jaccard", "top10_jaccard"),
    ("Precision@10", "precision_at_10"),
    ("Recall@10", "recall_at_10"),
    # NEW: Top-20 metrics
    ("Top-20 Jaccard", "top20_jaccard"),
    ("Precision@20", "precision_at_20"),
    ("Recall@20", "recall_at_20"),
    # Uncomment to include more:
    # ("Top-1 groups Jaccard", "group_jaccard_1"),
    # ("Top-2 groups Jaccard", "group_jaccard_2"),
    # ("Top-3 groups Jaccard", "group_jaccard_3"),
]

# Time columns: (display name, key)
TIME_COLS = [
    ("Simulate orig", "time_simulate_original"),
    ("Orig analysis", "time_original_analysis"),
    ("Propagation", "time_propagation"),
    # ("Process heading", "time_heading"),
]

ROUND_DEC = 3  # numeric formatting


# ---------------- Helpers ----------------

def _latex_escape(s: str) -> str:
    s = str(s)
    return (s.replace('\\', r'\textbackslash{}')
             .replace('&', r'\&')
             .replace('%', r'\%')
             .replace('$', r'\$')
             .replace('#', r'\#')
             .replace('_', r'\_')
             .replace('{', r'\{')
             .replace('}', r'\}')
             .replace('~', r'\textasciitilde{}')
             .replace('^', r'\textasciicircum{}'))

def _fmt_float(x):
    try:
        return f"{float(x):.{ROUND_DEC}f}"
    except Exception:
        return ""

def _pick_model_col(df: pd.DataFrame) -> str:
    if "model_basename" in df.columns:
        return "model_basename"
    if "model" in df.columns:
        return "model"
    # create one for safety
    df["model_basename"] = "MODEL"
    return "model_basename"

def _pick_target_col(df: pd.DataFrame) -> str | None:
    for c in df.columns:
        if c.lower() == "target":
            return c
    return None

def _ensure_numeric(df: pd.DataFrame, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

def _aggregate_model_target(df: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    """Mean over duplicate (model,target) rows if any."""
    model_col = _pick_model_col(df)
    tgt_col = _pick_target_col(df)
    id_cols = [model_col] + ([tgt_col] if tgt_col else [])
    keep = id_cols + [c for c in value_cols if c in df.columns]
    sub = df[keep].copy()
    _ensure_numeric(sub, keep[len(id_cols):])
    return sub.groupby(id_cols, as_index=False).mean(numeric_only=True)

def _prep_perf_headers(df: pd.DataFrame):
    """Return (headers, present_metric_cols) where headers is list of (disp, base).
       Only keeps metrics where BOTH KO and KI columns exist to keep a rectangular table."""
    headers = []
    present = []
    for disp, base in PERF_COLS:
        ko, ki = f"{base}_ko", f"{base}_ki"
        if ko in df.columns and ki in df.columns:
            headers.append((disp, base))
            present.extend([ko, ki])
    return headers, present

def _sanitize_num_cell(x):
    if x is None:
        return ""
    s = str(x)
    if s.lower() == "nan":
        return r"--"
    try:
        return f"{float(s):.{ROUND_DEC}f}"
    except Exception:
        return s


# ---------------- Performance: rotated full-bleed longtable (split into parts) ----------------

def build_performance_longtable_tex(df: pd.DataFrame, caption_base: str, label_base: str, split_parts: int = 2) -> str:
    """
    Rotated, full-bleed performance tables with Model & Target.
    Uses lscape+geometry correctly so the bottom whitespace disappears.
    Splits into N parts inside one .tex if needed.
    """
    if df.empty:
        return f"% No data for {caption_base}\n"

    model_col = _pick_model_col(df)
    tgt_col = _pick_target_col(df)

    # Which metric pairs are available
    headers, metric_cols = _prep_perf_headers(df)
    if not headers:
        return f"% No performance metrics found for {caption_base}\n"

    # Average duplicate (model,target) runs, sort
    agg = _aggregate_model_target(df, metric_cols)
    if agg.empty:
        return f"% No usable performance rows for {caption_base}\n"
    for c in metric_cols:
        if c in agg.columns:
            agg[c] = agg[c].map(_sanitize_num_cell)
    sort_cols = [model_col] + ([tgt_col] if tgt_col else [])
    agg = agg.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)

    # Col spec: 2 id cols + 2c per metric
    col_spec = "ll" + "cc" * len(headers)

    # Header rows
    top_line = " & ".join(
        [r"\textbf{Model}", r"\textbf{Target}"] +
        [rf"\multicolumn{{2}}{{c}}{{\textbf{{{h}}}}}" for (h, _) in headers]
    ) + r" \\"
    sub_line = " & ".join(["", ""] + [r"\textbf{KO} & \textbf{KI}"] * len(headers)) + r" \\"

    # Split rows into parts
    def chunk_df(frame, parts):
        if parts <= 1 or len(frame) == 0: return [frame]
        n = len(frame); sz = (n + parts - 1) // parts
        return [frame.iloc[i*sz:(i+1)*sz].copy() for i in range(parts)]

    parts = chunk_df(agg, split_parts)
    blocks = []

    for i, part in enumerate(parts, start=1):
        cap = caption_base + (f" (Part {i} of {split_parts})" if split_parts > 1 else "")
        cap = _latex_escape(cap).replace("—", "---")
        lab = f"{label_base}-p{i}"

        lines = []
        lines.append(r"\clearpage")
        # IMPORTANT: set landscape geometry BEFORE the environment
        lines.append(r"\newgeometry{landscape,left=15mm,right=15mm,top=15mm,bottom=15mm,headheight=0pt,headsep=0pt,footskip=0pt,includeheadfoot,heightrounded}")
        lines.append(r"\begin{landscape}")          # from lscape package
        lines.append(r"\thispagestyle{empty}")      # no header/footer artwork
        lines.append(r"\begingroup\flushbottom")    # fill vertical space; no ragged bottom
        lines.append(r"\small")
        # Longtable paddings
        lines.append(r"\setlength{\LTpre}{0pt}")
        lines.append(r"\setlength{\LTpost}{0pt}")
        lines.append(r"\setlength{\LTleft}{0pt}")
        lines.append(r"\setlength{\LTright}{0pt}")
        # tighter cells (adjust if you prefer)
        lines.append(r"\setlength{\tabcolsep}{2.5pt}")
        lines.append(r"\renewcommand{\arraystretch}{0.9}")

        lines.append(rf"\begin{{longtable}}{{{col_spec}}}")
        lines.append(rf"\caption{{{cap}}}\label{{{_latex_escape(lab)}}} \\")
        lines.append(r"\toprule")
        lines.append(top_line)
        lines.append(r"\midrule")
        lines.append(sub_line)
        lines.append(r"\midrule")
        lines.append(r"\endfirsthead")

        lines.append(r"\toprule")
        lines.append(top_line)
        lines.append(r"\midrule")
        lines.append(sub_line)
        lines.append(r"\midrule")
        lines.append(r"\endhead")

        lines.append(r"\midrule")
        lines.append(r"\endfoot")
        lines.append(r"\bottomrule")
        lines.append(r"\endlastfoot")

        # rows
        for _, row in part.iterrows():
            cells = [
                rf"\texttt{{{_latex_escape(row[model_col])}}}",
                rf"\texttt{{{_latex_escape(row[tgt_col])}}}" if tgt_col else r"\texttt{--}",
            ]
            for _disp, base in headers:
                cells += [row.get(f"{base}_ko", ""), row.get(f"{base}_ki", "")]
            lines.append(" & ".join(cells) + r" \\")
        lines.append(r"\end{longtable}")

        lines.append(r"\endgroup")
        lines.append(r"\end{landscape}")
        lines.append(r"\restoregeometry")  # restore portrait layout for next pages
        lines.append(r"\clearpage")

        blocks.append("\n".join(lines))

    return "\n".join(blocks) + "\n"


# ---------------- Time: per-(model,target) table (portrait) ----------------

def build_time_table_tex(df: pd.DataFrame, caption: str, label: str) -> str:
    """Portrait longtable (non-floating) so it always appears where \input is placed."""
    if df.empty:
        return rf"\par\noindent\emph{{No data for { _latex_escape(caption) }.}}\par\n"

    model_col = _pick_model_col(df)
    tgt_col = _pick_target_col(df)
    present = [key for _, key in TIME_COLS if key in df.columns]

    # Small debug note in the .tex to help if nothing matched
    if not present:
        dbg = ", ".join(k for _, k in TIME_COLS)
        return (r"\par\noindent\emph{No time metrics found. "
                rf"Expected one of: { _latex_escape(dbg) }." r"}\par" "\n")

    # Aggregate duplicate (model,target), compute Total
    agg = _aggregate_model_target(df, present)
    if agg.empty:
        return rf"\par\noindent\emph{{No usable time rows for { _latex_escape(caption) }.}}\par\n"

    agg["Total"] = agg[present].sum(axis=1, numeric_only=True)
    for c in present + ["Total"]:
        if c in agg.columns:
            agg[c] = agg[c].map(_fmt_float)

    # Build longtable (portrait)
    col_spec = "ll" + "c" * (len(present) + 1)  # Model, Target, metrics..., Total
    cap = _latex_escape(caption).replace("—", "---")

    # header row
    hdr_cells = [r"\textbf{Model}", r"\textbf{Target}"] + \
                [rf"\textbf{{{_latex_escape(d)}}}" for (d, k) in TIME_COLS if k in present] + \
                [r"\textbf{Total}"]
    top_line = " & ".join(hdr_cells) + r" \\"

    lines = []
    lines.append(r"\clearpage")                    # ensure it starts fresh
    lines.append(r"\begingroup")                   # keep sizing local
    lines.append(r"\small")                        # <- put \small here for time tables
    lines.append(r"\setlength{\tabcolsep}{4pt}")
    lines.append(r"\renewcommand{\arraystretch}{1.0}")
    lines.append(rf"\begin{{longtable}}{{{col_spec}}}")
    lines.append(rf"\caption{{{cap}}}\label{{{_latex_escape(label)}}} \\")
    lines.append(r"\toprule")
    lines.append(top_line)
    lines.append(r"\midrule")
    lines.append(r"\endfirsthead")

    lines.append(r"\toprule")
    lines.append(top_line)
    lines.append(r"\midrule")
    lines.append(r"\endhead")

    lines.append(r"\midrule")
    lines.append(r"\endfoot")
    lines.append(r"\bottomrule")
    lines.append(r"\endlastfoot")

    sort_cols = [model_col] + ([tgt_col] if tgt_col else [])
    for _, row in agg.sort_values(sort_cols, kind="mergesort").iterrows():
        cells = [
            rf"\texttt{{{_latex_escape(row[model_col])}}}",
            rf"\texttt{{{_latex_escape(row[tgt_col])}}}" if tgt_col else r"\texttt{--}",
        ]
        for _, key in TIME_COLS:
            if key in present:
                cells.append(row.get(key, ""))
        cells.append(row.get("Total", ""))
        lines.append(" & ".join(cells) + r" \\")
    lines.append(r"\end{longtable}")
    lines.append(r"\endgroup")
    lines.append(r"\clearpage")
    return "\n".join(lines) + "\n"



# ---------------- Driver ----------------

def build_all_for_csv(csv_path: Path, bench_name: str, outdir: Path, perf_split: int):
    df = pd.read_csv(csv_path)

    # PERFORMANCE (rotated, full-bleed longtable, possibly multiple parts in one file)
    perf_caption = f"{bench_name} benchmark — Performance (per target)"
    perf_label_base = f"tab:{bench_name.lower()}_perf"
    perf_tex = build_performance_longtable_tex(df, perf_caption, perf_label_base, split_parts=perf_split)
    (outdir / f"{bench_name.lower()}_performance.tex").write_text(perf_tex, encoding="utf-8")

    # TIME (portrait table)
    time_caption = f"{bench_name} benchmark — Time (seconds; per target)"
    time_label = f"tab:{bench_name.lower()}_time"
    time_tex = build_time_table_tex(df, time_caption, time_label)
    (outdir / f"{bench_name.lower()}_time.tex").write_text(time_tex, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acyclic", required=True, help="CSV for the acyclic benchmark (perf + time columns)")
    ap.add_argument("--normal", required=True, help="CSV for the normal benchmark (perf + time columns)")
    ap.add_argument("--outdir", default="tables", help="Output directory for .tex tables")
    ap.add_argument("--perf-split", type=int, default=2, help="Split performance table into N parts inside the same .tex")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    build_all_for_csv(Path(args.acyclic), "Acyclic", outdir, args.perf_split)
    build_all_for_csv(Path(args.normal),  "Normal",  outdir, args.perf_split)

    print("Done. Wrote LaTeX tables to:", outdir.resolve())
    print(" - acyclic_performance.tex (rotated longtable, parts inside)")
    print(" - acyclic_time.tex")
    print(" - normal_performance.tex  (rotated longtable, parts inside)")
    print(" - normal_time.tex")


if __name__ == "__main__":
    main()
