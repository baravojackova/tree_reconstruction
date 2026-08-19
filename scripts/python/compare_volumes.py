# -*- coding: utf-8 -*-
# =====================================================================
#  Compare tree-volume results from different methods against a reference.
# ---------------------------------------------------------------------
#  All results live in ONE simple "master" CSV (RESULTS_CSV below), one row
#  per (tree, method). Columns:
#
#     tree, method, total_m3, stem_m3, branch_m3, std_m3
#
#  - tree      : tree ID, e.g. "IND01_054"
#  - method    : a free-text label, e.g. "AdQSM", "TreeQSM mine v2"
#  - total_m3  : total wood volume [m^3]
#  - stem_m3   : stem/trunk volume [m^3]   (may be blank)
#  - branch_m3 : branch volume [m^3]        (may be blank)
#  - std_m3    : standard deviation of the total, if known (may be blank)
#
#  This script:
#     1) reads that CSV,
#     2) for the chosen tree (or ALL trees) prints every method side by side,
#        with its difference from the REFERENCE method (absolute and %),
#     3) across all trees, computes error metrics of every method vs the
#        reference: Bias, MAE, RMSE and CV-RMSE.
#
#  To ADD ANOTHER TREE later: just add more rows to the CSV (same columns)
#  and re-run. Nothing else changes.
#
#  If RESULTS_CSV does not exist yet, this script writes a starter file
#  pre-filled with the IND01_054 numbers we already have, so you can run it
#  immediately and then edit/extend it.
#
#  Dependencies: none (standard library only).
# =====================================================================

import csv
import os
import math

# =====================  PARAMETERS  ==================================
# The master results table. Edit this file to add trees/methods.
RESULTS_CSV = "volume_results.csv"

# Which tree to show in detail: a specific ID like "IND01_054", or "ALL".
SELECT_TREE = "ALL"

# The method used as the "truth" that everything is compared against. Must
# match a value in the 'method' column exactly.
REFERENCE_METHOD = "Reference (destructive)"
# =====================================================================


# Starter content written only if RESULTS_CSV does not exist yet.
STARTER_ROWS = [
    ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3"],
    ["IND01_054", "Reference (destructive)", "1.7169", "1.2557", "0.4612", ""],
    ["IND01_054", "AdQSM (TreesParams)",     "1.6905", "1.1302", "0.5603", ""],
    ["IND01_054", "AdTree calibrated",       "1.9240", "1.3020", "0.6220", ""],
    ["IND01_054", "TreeQSM de Tanago (mean 20)", "3.0838", "1.8383", "1.2455", "0.3113"],
    ["IND01_054", "TreeQSM mine v1",         "3.6806", "1.5317", "2.1488", "0.1915"],
    ["IND01_054", "TreeQSM mine v2",         "3.9744", "1.5737", "2.4007", "0.1923"],
]


def to_float(text):
    """Parse a number, or return None for blank/non-numeric cells."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def load_results(path):
    """Read the master CSV into a list of dict rows with parsed numbers."""
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "tree": r["tree"].strip(),
                "method": r["method"].strip(),
                "total": to_float(r.get("total_m3")),
                "stem": to_float(r.get("stem_m3")),
                "branch": to_float(r.get("branch_m3")),
                "std": to_float(r.get("std_m3")),
            })
    return rows


def fmt(x, width=10, dec=4):
    """Format a number or a blank cell to a fixed width."""
    if x is None:
        return " " * (width - 1) + "-"
    return ("%*.*f" % (width, dec, x))


def print_tree_block(tree, rows):
    """Print the side-by-side method comparison for a single tree."""
    tree_rows = [r for r in rows if r["tree"] == tree]
    ref = next((r for r in tree_rows if r["method"] == REFERENCE_METHOD), None)

    print("=" * 78)
    print("Tree: %s" % tree)
    print("-" * 78)
    print("%-28s %10s %10s %10s %10s" %
          ("method", "total", "stem", "branch", "d(total)"))
    print("-" * 78)

    # reference first (if present), then the rest in file order
    ordered = ([ref] if ref else []) + [r for r in tree_rows if r is not ref]
    for r in ordered:
        if r is None:
            continue
        if ref and ref["total"] is not None and r["total"] is not None:
            d_abs = r["total"] - ref["total"]
            d_pct = d_abs / ref["total"] * 100.0 if ref["total"] else 0.0
            if r is ref:
                dcol = "reference"
            else:
                dcol = "%+.3f (%+.0f%%)" % (d_abs, d_pct)
        else:
            dcol = "-"
        print("%-28s %s %s %s   %s" %
              (r["method"][:28], fmt(r["total"]), fmt(r["stem"]),
               fmt(r["branch"]), dcol))
    if ref is None:
        print("(no reference method '%s' for this tree - showing raw values only)"
              % REFERENCE_METHOD)
    print()


def error_metrics(rows):
    """Across all trees, compute Bias/MAE/RMSE/CV-RMSE of each method vs the
    reference, using only trees where BOTH that method and the reference exist."""
    trees = sorted({r["tree"] for r in rows})
    methods = [m for m in dict.fromkeys(r["method"] for r in rows)
               if m != REFERENCE_METHOD]

    # quick lookup: (tree, method) -> total
    total_of = {(r["tree"], r["method"]): r["total"] for r in rows}

    print("=" * 78)
    print("Error metrics vs '%s'  (total volume, across trees)" % REFERENCE_METHOD)
    print("-" * 78)
    print("%-28s %5s %9s %9s %9s %9s" %
          ("method", "n", "Bias", "MAE", "RMSE", "CV-RMSE%"))
    print("-" * 78)

    for m in methods:
        pairs = []
        for t in trees:
            est = total_of.get((t, m))
            ref = total_of.get((t, REFERENCE_METHOD))
            if est is not None and ref is not None:
                pairs.append((est, ref))
        if not pairs:
            continue
        n = len(pairs)
        errs = [e - r for e, r in pairs]
        refs = [r for _, r in pairs]
        bias = sum(errs) / n
        mae = sum(abs(e) for e in errs) / n
        rmse = math.sqrt(sum(e * e for e in errs) / n)
        cv_rmse = rmse / (sum(refs) / n) * 100.0 if refs else 0.0
        print("%-28s %5d %9.3f %9.3f %9.3f %9.1f" %
              (m[:28], n, bias, mae, rmse, cv_rmse))
    print("(Bias/MAE/RMSE in m^3. n = number of trees compared. With one tree,")
    print(" MAE = RMSE = |error| and CV-RMSE is that error relative to the reference.)")
    print()


# =========================  RUN  =====================================
if not os.path.exists(RESULTS_CSV):
    with open(RESULTS_CSV, "w", encoding="utf-8", newline="") as f:
        csv.writer(f).writerows(STARTER_ROWS)
    print("No results file found - wrote a starter '%s' with the IND01_054 data.\n"
          % RESULTS_CSV)

rows = load_results(RESULTS_CSV)
all_trees = sorted({r["tree"] for r in rows})

if SELECT_TREE.upper() == "ALL":
    for t in all_trees:
        print_tree_block(t, rows)
elif SELECT_TREE in all_trees:
    print_tree_block(SELECT_TREE, rows)
else:
    print("Tree '%s' not found. Available trees: %s" % (SELECT_TREE, ", ".join(all_trees)))
    raise SystemExit

# Error metrics make sense whenever a reference exists (works for 1 or many trees).
error_metrics(rows)
