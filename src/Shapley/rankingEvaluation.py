from math import isfinite
from statistics import mean
from collections import defaultdict

# ---------- TIE-AWARE BUCKETING ----------

def _close(a, b, atol, rtol):
    return abs(a - b) <= max(atol, rtol * max(abs(a), abs(b)))

def tie_buckets_abs(values, atol=1e-6, rtol=1e-6, rep="max"):
    """
    Build epsilon-tolerant tie buckets by absolute value (descending).
    Returns: list of buckets, each bucket = list[(key, abs_value)]
             Buckets are ordered from highest abs value to lowest.
    """
    items = [(k, abs(v)) for k, v in values.items() if v is not None and isfinite(v)]
    items.sort(key=lambda kv: kv[1], reverse=True)

    buckets, cur = [], []
    for k, val in items:
        if not cur:
            cur = [(k, val)]
            continue
        # compare to representative of current bucket
        if rep == "max":
            ref = max(x[1] for x in cur)
        elif rep == "mean":
            ref = sum(x[1] for x in cur) / len(cur)
        else:
            raise ValueError("rep must be 'max' or 'mean'")
        if _close(val, ref, atol, rtol):
            cur.append((k, val))
        else:
            buckets.append(cur)
            cur = [(k, val)]
    if cur:
        buckets.append(cur)
    return buckets

def ranks_from_buckets(buckets, scheme="dense"):
    """
    Convert buckets -> ranks (1 is best).
    scheme='dense': all in bucket j get rank j (1,2,3,...)
    scheme='average': average competition ranks (1..m ties get (1+...+m)/m, next gets m+1)
    Returns: dict[key] -> float rank
    """
    if scheme not in {"dense", "average"}:
        raise ValueError("scheme must be 'dense' or 'average'")

    ranks = {}
    if scheme == "dense":
        for j, bucket in enumerate(buckets, start=1):
            for k, _ in bucket:
                ranks[k] = float(j)
    else:
        # average competition ranking
        used = 0
        for bucket in buckets:
            start = used + 1
            stop = used + len(bucket)
            avg_rank = (start + stop) / 2.0
            for k, _ in bucket:
                ranks[k] = float(avg_rank)
            used = stop
    return ranks

# ---------- UTILITIES ----------

def intersect_keys(d1, d2):
    """Keep only keys present in both dicts."""
    common = set(d1).intersection(d2)
    return {k: d1[k] for k in common}, {k: d2[k] for k in common}

def pearsonr(x, y):
    """Simple Pearson correlation (no scipy)."""
    mx, my = mean(x), mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denx = (sum((a - mx) ** 2 for a in x) ** 0.5)
    deny = (sum((b - my) ** 2 for b in y) ** 0.5)
    if denx == 0 or deny == 0:
        return 0.0
    return num / (denx * deny)

# ---------- KENDALL TAU-B ----------

def kendall_tau_b(r1, r2):
    """
    r1, r2: dict[key] -> rank (floats allowed, ties allowed)
    Returns tau-b. Uses scipy if available; else O(n^2) fallback.
    """
    keys = list(set(r1).intersection(r2))
    if not keys:
        return 0.0
    x = [r1[k] for k in keys]
    y = [r2[k] for k in keys]

    # Try scipy (fast, robust)
    try:
        from scipy.stats import kendalltau  # type: ignore
        tau, _ = kendalltau(x, y, variant='b')
        return float(tau)
    except Exception:
        pass

    # Fallback: O(n^2) tau-b
    n = len(keys)
    concordant = discordant = 0
    ties_x = ties_y = 0

    for i in range(n - 1):
        for j in range(i + 1, n):
            xi, xj = x[i], x[j]
            yi, yj = y[i], y[j]
            dx = xi - xj
            dy = yi - yj
            if dx == 0 and dy == 0:
                # tied in both -> neither concordant nor discordant, handled via tie terms
                ties_x += 1
                ties_y += 1
            elif dx == 0:
                ties_x += 1
            elif dy == 0:
                ties_y += 1
            else:
                prod = dx * dy
                if prod > 0:
                    concordant += 1
                elif prod < 0:
                    discordant += 1

    num = concordant - discordant
    # tie correction
    denom = ((concordant + discordant + ties_x) * (concordant + discordant + ties_y)) ** 0.5
    return float(num / denom) if denom != 0 else 0.0

# ---------- TOP-OVERLAP METRICS ----------

def expand_to_cover_r(buckets, r):
    """
    Given ordered tie-buckets, take as many whole buckets as needed to reach >= r items.
    Returns the set of keys selected.
    """
    chosen = []
    for bucket in buckets:
        chosen.extend(k for k, _ in bucket)
        if len(chosen) >= r:
            break
    return set(chosen)

def top_r_overlap(buckets1, buckets2, r):
    """
    Tie-aware Top-r overlap by expanding buckets to cover r.
    Returns dict with overlap, jaccard, precision, recall.
    """
    A = expand_to_cover_r(buckets1, r)
    B = expand_to_cover_r(buckets2, r)
    inter = len(A & B)
    union = len(A | B)
    return {
        "overlap": inter / r if r > 0 else 0.0,
        "jaccard": inter / union if union > 0 else 0.0,
        "precision_at_r": inter / len(B) if len(B) > 0 else 0.0,
        "recall_at_r": inter / len(A) if len(A) > 0 else 0.0,
        "A_size": len(A),
        "B_size": len(B),
    }

def top_group_overlap(buckets1, buckets2, k_groups=1):
    """
    Compare the union of the top-k tie-buckets (k_groups >= 1).
    Returns Jaccard and raw sizes.
    """
    A = set(k for bucket in buckets1[:k_groups] for k, _ in bucket)
    B = set(k for bucket in buckets2[:k_groups] for k, _ in bucket)
    inter = len(A & B)
    union = len(A | B)
    return {
        "jaccard": inter / union if union > 0 else 0.0,
        "A_size": len(A),
        "B_size": len(B),
        "intersection": inter,
    }

# ---------- MAIN COMPARISON ----------

def compare_rankings(
    baseline: dict,
    propagated: dict,
    *,
    use_abs=True,
    atol=1e-4,
    rtol=1e-4,
    rep="max",
    rank_scheme="dense",   # 'dense' or 'average'
    top_r_values=(10, 20),
    top_group_ks=(1,2,3),     # e.g., (1,) or (1,2)
):
    """
    Compute tie- & tolerance-aware ranking metrics between two dicts.
    Returns: dict with metrics and intermediate artifacts.
    """
    # 1) Align keys
    b0, p0 = intersect_keys(baseline, propagated)
    if not b0:
        raise ValueError("No overlapping keys between baseline and propagated.")

    # 2) Build buckets (abs vs raw)
    if use_abs:
        b_buckets = tie_buckets_abs(b0, atol=atol, rtol=rtol, rep=rep)
        p_buckets = tie_buckets_abs(p0, atol=atol, rtol=rtol, rep=rep)
    else:
        # If you need signed values, adapt tie_buckets_abs accordingly
        raise NotImplementedError("Signed-value bucketing not implemented; set use_abs=True.")

    # 3) Ranks from buckets
    b_ranks = ranks_from_buckets(b_buckets, scheme=rank_scheme)
    p_ranks = ranks_from_buckets(p_buckets, scheme=rank_scheme)

    # 4) Kendall tau-b
    tau_b = kendall_tau_b(b_ranks, p_ranks)

    # 5) Spearman rho on ranks
    #    (compute over common keys; they already are)
    keys = list(b_ranks.keys() & p_ranks.keys())
    x = [b_ranks[k] for k in keys]
    y = [p_ranks[k] for k in keys]
    spearman_rho = pearsonr(x, y)  # Pearson on rank vectors

    # 6) Mean absolute rank difference (MARD)
    mard = sum(abs(b_ranks[k] - p_ranks[k]) for k in keys) / len(keys)

    # 7) Top-r overlap(s)
    top_r_results = {r: top_r_overlap(b_buckets, p_buckets, r) for r in top_r_values}

    # 8) Top-group overlap(s)
    top_group_results = {k: top_group_overlap(b_buckets, p_buckets, k_groups=k) for k in top_group_ks}

    return {
        "n_common": len(keys),
        "kendall_tau_b": tau_b,
        "spearman_rho_on_ranks": spearman_rho,
        "mean_abs_rank_diff": mard,
        "top_r_overlap": top_r_results,
        "top_group_overlap": top_group_results,
        "meta": {
            "use_abs": use_abs,
            "atol": atol,
            "rtol": rtol,
            "rep": rep,
            "rank_scheme": rank_scheme,
            "top_r_values": top_r_values,
            "top_group_ks": top_group_ks,
            "buckets_baseline_sizes": [len(b) for b in b_buckets],
            "buckets_propagated_sizes": [len(b) for b in p_buckets],
        }
    }


def print_rank_comparison(ko_result, ki_result, name="Ranking Comparison"):
    """
    Print a side-by-side summary of ranking metrics for KO and KI runs.
    Both results should come from compare_rankings(...).
    """
    def fmt(x, prec=3):
        if isinstance(x, (int, float)):
            if x != x:  # NaN check
                return "NaN"
            return f"{x:.{prec}f}"
        return str(x)

    print(f"\n{name}")
    print("=" * len(name))
    print(f"{'Metric':<30} {'KO':>12} {'KI':>12}")
    print("-" * 58)

    # Core global metrics
    metrics = [
        ("n_common", ko_result["n_common"], ki_result["n_common"]),
        ("Kendall tau-b", ko_result["kendall_tau_b"], ki_result["kendall_tau_b"]),
        ("Spearman rho", ko_result["spearman_rho_on_ranks"], ki_result["spearman_rho_on_ranks"]),
        ("Mean abs rank diff", ko_result["mean_abs_rank_diff"], ki_result["mean_abs_rank_diff"]),
    ]
    for name, ko, ki in metrics:
        print(f"{name:<30} {fmt(ko):>12} {fmt(ki):>12}")

    print("-" * 58)

    # Top-r overlap
    print("\nTop-r Overlap:")
    for r in ko_result["top_r_overlap"]:
        ko_t = ko_result["top_r_overlap"][r]
        ki_t = ki_result["top_r_overlap"].get(r, {})
        jacc_ko, jacc_ki = ko_t.get("jaccard", float("nan")), ki_t.get("jaccard", float("nan"))
        prec_ko, prec_ki = ko_t.get("precision_at_r", float("nan")), ki_t.get("precision_at_r", float("nan"))
        rec_ko, rec_ki = ko_t.get("recall_at_r", float("nan")), ki_t.get("recall_at_r", float("nan"))
        print(f"  Top-{r:<2} Jaccard{'':<16}{fmt(jacc_ko):>8} {fmt(jacc_ki):>12}")
        print(f"         Precision@{r:<2}{'':<13}{fmt(prec_ko):>8} {fmt(prec_ki):>12}")
        print(f"         Recall@{r:<2}{'':<16}{fmt(rec_ko):>8} {fmt(rec_ki):>12}")
    print("-" * 58)

    # Top-group overlap
    print("\nTop-group Overlap:")
    for k in ko_result["top_group_overlap"]:
        ko_g = ko_result["top_group_overlap"][k]
        ki_g = ki_result["top_group_overlap"].get(k, {})
        jacc_ko, jacc_ki = ko_g.get("jaccard", float("nan")), ki_g.get("jaccard", float("nan"))
        print(f"  Top-{k:<2} groups Jaccard{'':<9}{fmt(jacc_ko):>8} {fmt(jacc_ki):>12}")
    print("=" * 58)


# Example usage:
# print_rank_comparison(ko_result, ki_result, name="Ranking stability between KO and KI")

# # ---------- EXAMPLE USAGE ----------
# if __name__ == "__main__":
#     baseline = {"A": 1.6667, "B": 0.2, "C": 0.1999999, "D": 1.6668, "E": 0.05}
#     propagated = {"A": 1.66675, "B": 0.21, "C": 0.2, "D": 1.6666, "E": 0.0499}

#     results = compare_rankings(
#         baseline, propagated,
#         use_abs=True,
#         atol=1e-6, rtol=1e-5,  # adjust to your numeric noise level
#         rep="max",
#         rank_scheme="dense",
#         top_r_values=(3, 5),
#         top_group_ks=(1, )
#     )
#     import pprint; pprint.pprint(results)
