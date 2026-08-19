# -*- coding: utf-8 -*-
# =====================================================================
#  Compare tree-volume results from different methods against a reference.
# ---------------------------------------------------------------------
#  All results live in ONE simple "master" CSV (RESULTS_CSV below), one row
#  per (tree, method). Columns:
#
#     tree, method, total_m3, stem_m3, branch_m3, std_m3, dbh_m, height_m,
#     taper_cm_per_m
#
#  - tree           : tree ID, e.g. "IND01_054"
#  - method         : a free-text label, e.g. "AdQSM", "TreeQSM mine v2"
#  - total_m3       : total wood volume [m^3]
#  - stem_m3        : stem/trunk volume [m^3]   (may be blank)
#  - branch_m3      : branch volume [m^3]        (may be blank)
#  - std_m3         : standard deviation of the total, if known (may be blank)
#  - dbh_m          : stem diameter at 1.3 m above the tree base [m] (may be blank)
#  - height_m       : total tree height [m]                          (may be blank)
#  - taper_cm_per_m : stem taper between two reference heights [cm/m] (may be blank)
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
    ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3",
     "dbh_m", "height_m", "taper_cm_per_m"],
    ["IND01_054", "Reference (destructive)", "1.7169", "1.2557", "0.4612", "", "", "", ""],
    ["IND01_054", "AdQSM (TreesParams)",     "1.6905", "1.1302", "0.5603", "", "", "", ""],
    ["IND01_054", "AdTree calibrated",       "1.9240", "1.3020", "0.6220", "", "", "", ""],
    ["IND01_054", "TreeQSM de Tanago (mean 20)", "3.0838", "1.8383", "1.2455", "0.3113", "", "", ""],
    ["IND01_054", "TreeQSM mine v1",         "3.6806", "1.5317", "2.1488", "0.1915", "", "", ""],
    ["IND01_054", "TreeQSM mine v2",         "3.9744", "1.5737", "2.4007", "0.1923", "", "", ""],
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
                "dbh": to_float(r.get("dbh_m")),
                "height": to_float(r.get("height_m")),
                "taper": to_float(r.get("taper_cm_per_m")),
            })
    return rows


def fmt(x, width=10, dec=4):
    """Format a number or a blank cell to a fixed width."""
    if x is None:
        return " " * (width - 1) + "-"
    return ("%*.*f" % (width, dec, x))


def pct_diff(value, ref_value):
    """Percent difference of value vs. ref_value, or None if either is
    missing (blank cell in the printed table)."""
    if value is None or ref_value is None or ref_value == 0:
        return None
    return (value - ref_value) / ref_value * 100.0


def fmt_pct(x, width=8):
    """Format a percent difference, or a blank cell if None."""
    if x is None:
        return " " * (width - 1) + "-"
    return "%+*.1f%%" % (width - 1, x)


def print_tree_block(tree, rows):
    """Print the side-by-side method comparison for a single tree."""
    tree_rows = [r for r in rows if r["tree"] == tree]
    ref = next((r for r in tree_rows if r["method"] == REFERENCE_METHOD), None)

    print("=" * 118)
    print("Tree: %s" % tree)
    print("-" * 118)
    print("%-28s %10s %10s %10s %10s   %8s %8s   %8s %8s   %8s %8s" %
          ("method", "total", "stem", "branch", "d(total)",
           "DBH[m]", "d(DBH)", "H[m]", "d(H)", "taper", "d(taper)"))
    print("-" * 118)

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

        dbh_pct = None if r is ref else pct_diff(r["dbh"], ref["dbh"] if ref else None)
        height_pct = None if r is ref else pct_diff(r["height"], ref["height"] if ref else None)
        taper_pct = None if r is ref else pct_diff(r["taper"], ref["taper"] if ref else None)

        print("%-28s %s %s %s   %-16s %s %s   %s %s   %s %s" %
              (r["method"][:28], fmt(r["total"]), fmt(r["stem"]),
               fmt(r["branch"]), dcol,
               fmt(r["dbh"], width=8, dec=3), fmt_pct(dbh_pct),
               fmt(r["height"], width=8, dec=2), fmt_pct(height_pct),
               fmt(r["taper"], width=8, dec=2), fmt_pct(taper_pct)))
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


def field_error_summary(rows, key, label):
    """Print one short line per method: its mean percent error vs. the
    reference for a single field (dbh/height/taper), across all trees where
    both values are known. Kept OUT of the volume RMSE metrics block above
    by design - these fields are not volumes."""
    trees = sorted({r["tree"] for r in rows})
    methods = [m for m in dict.fromkeys(r["method"] for r in rows) if m != REFERENCE_METHOD]
    value_of = {(r["tree"], r["method"]): r[key] for r in rows}

    print("%s error vs '%s':" % (label, REFERENCE_METHOD))
    printed_any = False
    for m in methods:
        diffs = []
        for t in trees:
            d = pct_diff(value_of.get((t, m)), value_of.get((t, REFERENCE_METHOD)))
            if d is not None:
                diffs.append(d)
        if not diffs:
            continue
        printed_any = True
        print("  %-28s %+6.1f%%  (n=%d)" % (m[:28], sum(diffs) / len(diffs), len(diffs)))
    if not printed_any:
        print("  (no method has both a value and a reference value for this field)")
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

# DBH/height/taper are not volumes, so their error is reported separately
# instead of folding them into the RMSE metrics block above.
field_error_summary(rows, "dbh", "DBH")
field_error_summary(rows, "height", "Height")
field_error_summary(rows, "taper", "Taper")
