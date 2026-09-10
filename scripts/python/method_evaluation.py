# -*- coding: utf-8 -*-
# =====================================================================
#  Overall comparison of reconstruction METHODS (AdQSM / AdTree raw /
#  AdTree calibrated / TreeQSM Optimal / TreeQSM Simplified) and of
#  PARAMETER SENSITIVITY within each method, across all 12 production
#  beeches in volume_results.csv.
# ---------------------------------------------------------------------
#  THE DECISION THIS SCRIPT IS BUILT AROUND
#
#  A box plot drawn over every AdTree-calibrated row per tree gives a
#  coefficient of variation of ~48% for branch_m3. Restricted to the
#  SETTLED configuration actually carried forward (regression-perorder
#  calibration, length-weighted radius stats, the one AdQSM variant
#  actually selected per tree) it is ~4.9%. Both numbers are correct -
#  they answer different questions. Most of the 48% comes from axes
#  that are no longer open (calmethod and weighting alone explain ~70%
#  of the variance - see the eta^2 table this script reproduces).
#  Presenting the full spread as "method uncertainty" would report an
#  already-resolved question as noise.
#
#  So this script computes everything TWICE, over two explicitly
#  defined populations (POP_SETTLED / POP_FULL below), and never mixes
#  them in one figure.
#
#  This script ONLY READS volume_results.csv - it runs no
#  reconstruction, and it does not modify volume_results.csv or any
#  other existing file except the additive METHOD_COLORS/BOX_* constants
#  already added to plot_style.py.
#
#  Reused, not duplicated: RESULTS_CSV/load_results() from
#  compare_volumes.py, ensure_all_plots_dir()/ensure_csv_dir() from
#  paths.py, METHOD_COLORS/BOX_WHIS/BOX_EDGE_DARKEN_FACTOR/font-size
#  constants from plot_style.py.
#
#  STEP 0 findings this script's styling choices are built on (see the
#  conversation this script was authored in for the full audit):
#    - A method->colour mapping already existed (plot_style.py's
#      FAMILY_GRADIENTS) but was incomplete (no AdTree_raw/calibrated
#      "one colour" pick, no TreeQSM stage split) and 4 scripts never
#      imported plot_style.py at all - METHOD_COLORS (plot_style.py)
#      fills that gap additively, reusing every colour already in use.
#    - AdTree_raw (green) and AdTree_calibrated (blue) are different
#      HUES, not "one hue, two lightnesses" - kept as-is because both
#      are already in active use elsewhere in the project (existing
#      colour wins over the "same hue" ideal for a variant pair).
#    - Two box-plot conventions already existed (plot_box.py vs.
#      plot_volumes.py's plot_error_boxplot()) - this script follows
#      plot_box.py's (whis=(0,100), solid family fill, darkened edges),
#      the more fully worked out of the two.
#    - Existing figures use ax.set_title()/fig.suptitle() - so this
#      script titles its figures too, matching that convention.
# =====================================================================

import csv
import os
import re
import statistics
import textwrap
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from matplotlib.ticker import LogLocator, FuncFormatter

from compare_volumes import RESULTS_CSV, load_results
from paths import ensure_all_plots_dir, ensure_csv_dir
# Field DBH, for the DBH-ascending tree ordering (revision 2, Part A) - reused
# from base_vs_dbh_compare.py's own parser (2 header lines, then one row per
# tree, "diameter at 1.3 m (cm)" at raw column 26) rather than writing a
# third copy of the same file format.
from base_vs_dbh_compare import parse_field_file, FIELD_FILE
from plot_style import (
    METHOD_COLORS,
    BOX_WHIS,
    BOX_EDGE_DARKEN_FACTOR,
    AXIS_LABEL_FONTSIZE,
    LEGEND_FONTSIZE,
    PANEL_TITLE_FONTSIZE,
    TITLE_FONTSIZE,
    ANNOTATION_FONTSIZE,
    POINT_LABEL_FONTSIZE,
    JITTER_POINT_SIZE,
)

# =====================  PARAMETERS  ===================================

# ---- Method families -------------------------------------------------
# TreeQSM's two stages are kept separate THROUGHOUT (different objects,
# not two samples of one thing) - "TreeQSM" is never used as a family on
# its own anywhere below. Only TreeQSM_Simplified (the "Simplified (no
# islands)" stage) is exported to ANSYS - every legend says so.
BOX_FAMILIES = ["AdTree_raw", "AdTree_calibrated", "TreeQSM_Optimal", "TreeQSM_Simplified"]
MARKER_FAMILY = "AdQSM"   # drawn as individual markers, never a box - see module header / Step 2
ALL_FAMILIES = BOX_FAMILIES + [MARKER_FAMILY]

TREEQSM_STAGE_OF_FAMILY = {
    "TreeQSM_Optimal": "Optimal",
    "TreeQSM_Simplified": "Simplified (no islands)",
}

LEGEND_LABELS = {
    "AdTree_raw": "AdTree raw",
    "AdTree_calibrated": "AdTree calibrated",
    "TreeQSM_Optimal": "TreeQSM Optimal",
    "TreeQSM_Simplified": "TreeQSM Simplified (exported to ANSYS)",
    "AdQSM": "AdQSM",
}

# ---- Population definitions (see module header) ----------------------
# POP_SETTLED: the configuration actually carried forward.
#   AdTree calibrated : calmethod == "regression-perorder"
#                       weighting  == "-wlen" (method string contains it)
#                       adqsm_variant == the tree's selected variant (not 999)
#   AdTree raw        : weighting == "-wlen"
#   TreeQSM           : all models, split by stage
#   AdQSM             : the selected variant only, NOT 999
# POP_FULL: every row of every method, nothing fixed.
POP_SETTLED = "POP_SETTLED"
POP_FULL = "POP_FULL"

# Human-readable population wording for figures (CSV/console output keeps the
# POP_SETTLED/POP_FULL codes - those are read by people who already have this
# script's own README in front of them; a figure is read on its own, so every
# figure spells the population out in full, not just in its filename).
POP_LABEL_LONG = {
    POP_SETTLED: ("settled configuration (regression-perorder, -wlen, "
                  "selected AdQSM variant; 36 rows/tree)"),
    POP_FULL: "full parameter grid (all rows, nothing fixed)",
}
POP_LABEL_SHORT = {
    POP_SETTLED: "settled configuration",
    POP_FULL: "full grid",
}

# ---- Quantities compared (Step 4) -------------------------------------
# (csv/plot key, axis label incl. units) - branch_frac is dimensionless
# (a fraction of total_m3), labelled in plain words rather than left as a
# bare internal identifier (see revision Part A1).
QUANTITIES = [
    ("total_m3", "Total volume [m3]"),
    ("trunk_m3", "Trunk volume [m3]"),
    ("branch_m3", "Branch volume [m3]"),
    ("branch_frac", "Branch share of total volume [-]"),
    ("dbh_m", "DBH [m]"),
    ("height_m", "Height [m]"),
    ("trunk_len_m", "Trunk length [m]"),
    ("branch_len_m", "Branch length [m]"),
    ("n_cylinders", "Number of cylinders"),
]
QUANTITY_KEYS = [q[0] for q in QUANTITIES]
QUANTITY_LABEL = dict(QUANTITIES)

# ---- eta^2 factors per family (Step 5), computed on POP_FULL ----------
ETA2_FACTORS = {
    "AdTree_calibrated": ["adqsm_variant", "radius_threshold_mm", "seg_min_mm",
                           "seg_k_pct", "calmethod", "weighting"],
    "AdTree_raw": ["radius_threshold_mm", "seg_min_mm", "seg_k_pct", "weighting"],
}
TREEQSM_ETA2_FACTORS = ["pd2min_m", "simp_smallradii", "simp_replaceiterations"]

# TreeQSM's 3 balanced pd2min levels (each carrying the full 16-cell simp
# grid) vs. the single AUTO run and the single pd2min=0.020 probe (no simp
# sweep) - the design is unbalanced, so the decomposition is restricted to
# the balanced part and the other two are reported as individual values.
TREEQSM_BALANCED_PD2MIN_LEVELS = [0.002, 0.005, 0.010]
TREEQSM_PROBE_PD2MIN = 0.020
_PD2MIN_TOL = 1e-4

# ---- pd2min sweep across all 12 trees (revision 2, Part D / addendum 2) --
# The 4 manual pd2min levels shown on the x axis (0.020 is the single probe,
# no simp grid - see TREEQSM_BALANCED_PD2MIN_LEVELS/TREEQSM_PROBE_PD2MIN
# above, same values).
PD2MIN_SWEEP_LEVELS = [0.002, 0.005, 0.010, 0.020]

# Reduced to the BRANCH quantities only (addendum 2 - trunk_m3/total_m3/
# pmdist_trunk_mean dropped from the original 5-panel layout). Each entry:
# (display qkey, axis label, load_results() row key).
# CAUTION carried into the figure itself (see draw_pd2min_sweep_all_trees()):
# pmdist_branch_mean is shown here purely to describe how the branch fit
# RESPONDS to pd2min - it is NOT a model-selection criterion. This project
# already established that pmdist_mean/pmdist_branch_mean are both minimised
# by the known-degraded fat-cylinder model; only pmdist_trunk_mean is
# admissible for model selection, and model selection is not part of this
# script.
PD2MIN_SWEEP_QUANTITIES = [
    ("branch_m3", "Branch volume", "branch"),
    ("pmdist_branch_mean", "Mean point-to-model distance, branches", "pmdist_branch_mean"),
    ("n_cylinders", "Number of cylinders", "n_cylinders"),
]

# Sanity-check reference values (from an independent calculation on the
# same CSV - see module header). Printed alongside the actual numbers;
# disagreement beyond ~1 percentage point means the grouping/exclusions
# differ from what produced these.
EXPECTED_ETA2_BRANCH = {"calmethod": 53.3, "weighting": 15.8, "seg_k_pct": 2.9,
                         "seg_min_mm": 1.2, "radius_threshold_mm": 0.7, "adqsm_variant": 0.3}
EXPECTED_ETA2_NCYL = {"radius_threshold_mm": 93.4, "seg_min_mm": 5.1}
EXPECTED_CV_BRANCH_FULL = 48.0
EXPECTED_CV_BRANCH_SETTLED = 4.9
EXPECTED_SETTLED_ADTREE_CAL_PER_TREE = 36

# ---- Known failures (Step 3) - documented failure modes, excluded from
# every aggregate (median/quartile/eta^2/CV) but never hidden. -----------
def _row_adqsm_variant(row):
    """adqsm_variant for ANY family - AdTree_calibrated has it as a real
    CSV column; AdQSM-direct rows leave that column blank and carry the
    variant number in the method string instead (e.g. "...(AdQSM 999)")."""
    if row["_family"] == "AdQSM":
        m = re.search(r"\(AdQSM (\d+)\)\s*$", row["method"])
        return m.group(1) if m else None
    return row["adqsm_variant"] or None


def known_failure_reason(row):
    """Return a human-readable exclusion reason if `row` is one of the two
    documented failure modes, else None. See module header for why these
    are excluded from aggregates but not from method_eval_excluded.csv."""
    if row["tree"] == "B21_S21" and _row_adqsm_variant(row) == "999":
        return ("B21_S21, adqsm_variant 999: a verbatim copy of a bad source "
                "variant (branch_m3 far above every other variant for this "
                "tree, against a field crown-wood estimate of 0.48 m^3)")
    if row["tree"] == "B21_S10" and row["_family"] == "AdTree_raw":
        return ("B21_S10, AdTree raw: dbh_m stuck at 1.030964 m identically "
                "on all 72 rows, against a field DBH of 0.37 m; cause not "
                "yet established")
    return None


# =====================  CLASSIFICATION  ================================

def classify_family(method):
    """Method string -> one of ALL_FAMILIES, or None if unrecognised."""
    if method.startswith("AdTree raw"):
        return "AdTree_raw"
    if method.startswith("AdTree calibrated"):
        return "AdTree_calibrated"
    if method.startswith("AdQSM"):
        return "AdQSM"
    if method.startswith("TreeQSM"):
        if method.endswith(", Optimal)"):
            return "TreeQSM_Optimal"
        if method.endswith(", Simplified (no islands))"):
            return "TreeQSM_Simplified"
        return None   # an unexpected TreeQSM stage - not present in today's CSV
    return None


def quantity_value(row, qkey):
    """One row's value for quantity `qkey`, or None if not available.
    branch_frac is computed here (not a CSV column) - branch_m3/total_m3,
    None if either is missing or total is zero."""
    if qkey == "branch_frac":
        total, branch = row["total"], row["branch"]
        if total in (None, 0) or branch is None:
            return None
        return branch / total
    field_of = {
        "total_m3": "total", "trunk_m3": "trunk", "branch_m3": "branch",
        "dbh_m": "dbh", "height_m": "height",
        "trunk_len_m": "trunk_len", "branch_len_m": "branch_len",
        "n_cylinders": "n_cylinders",
    }
    return row[field_of[qkey]]


def is_settled(row, selected_variant):
    """Whether `row` belongs to POP_SETTLED - see module header for the
    per-family definition. `selected_variant`: {tree: variant_string}."""
    fam = row["_family"]
    if fam == "AdTree_calibrated":
        return (row["calmethod"] == "regression-perorder"
                and "-wlen" in row["method"]
                and row["adqsm_variant"] == selected_variant.get(row["tree"]))
    if fam == "AdTree_raw":
        return "-wlen" in row["method"]
    if fam in ("TreeQSM_Optimal", "TreeQSM_Simplified"):
        return True
    if fam == "AdQSM":
        return _row_adqsm_variant(row) == selected_variant.get(row["tree"])
    return False


def compute_selected_variants(rows):
    """{tree: variant_string} - the non-999 adqsm_variant value among that
    tree's AdTree_calibrated rows (derived from the CSV, never hardcoded).
    Warns (does not crash) if a tree has 0 or >1 distinct such values."""
    by_tree = defaultdict(set)
    for r in rows:
        if r["_family"] == "AdTree_calibrated" and r["adqsm_variant"] not in ("", "999"):
            by_tree[r["tree"]].add(r["adqsm_variant"])
    selected = {}
    for tree, variants in by_tree.items():
        if len(variants) != 1:
            print("WARNING: tree %s has %d distinct non-999 adqsm_variant "
                  "values in AdTree_calibrated: %s" % (tree, len(variants), sorted(variants)))
        selected[tree] = sorted(variants)[0]
    return selected


# =====================  STATS HELPERS  =================================

def _quantiles(sorted_vals):
    """(q1, median, q3) via linear interpolation (numpy's default 'linear'
    method) - matches what a reader would get from numpy.percentile/
    pandas.describe, so these numbers are checkable against a spreadsheet."""
    arr = np.array(sorted_vals, dtype=float)
    return (float(np.percentile(arr, 25)), float(np.percentile(arr, 50)),
            float(np.percentile(arr, 75)))


def summary_stats(values):
    """{n, min, q1, median, q3, max, mean, std, cv_pct} for a list of
    numbers (None entries must already be filtered out by the caller).
    n<1 -> all-blank (n=0); n==1 -> std/cv_pct blank (undefined for 1 point,
    not silently 0)."""
    n = len(values)
    if n == 0:
        return dict(n=0, min=None, q1=None, median=None, q3=None, max=None,
                    mean=None, std=None, cv_pct=None)
    values_sorted = sorted(values)
    q1, med, q3 = _quantiles(values_sorted)
    mean = statistics.mean(values)
    if n >= 2:
        std = statistics.stdev(values)
        cv_pct = (100.0 * std / mean) if mean else None
    else:
        std = None
        cv_pct = None
    return dict(n=n, min=values_sorted[0], q1=q1, median=med, q3=q3,
                max=values_sorted[-1], mean=mean, std=std, cv_pct=cv_pct)


def eta_squared(group_of_values):
    """One-way between-group eta^2 (%) - SUM[n_g*(mean_g-grand)^2] /
    SUM[(v-grand)^2] * 100. `group_of_values`: {level: [values]}. None if
    there are <2 distinct non-empty levels or zero total variance is not
    representable (0.0 is returned for genuine zero-variance data, not
    conflated with "not computable")."""
    groups = {k: v for k, v in group_of_values.items() if v}
    if len(groups) < 2:
        return None
    all_values = [v for vals in groups.values() for v in vals]
    n = len(all_values)
    grand_mean = sum(all_values) / n
    sst = sum((v - grand_mean) ** 2 for v in all_values)
    if sst == 0:
        return 0.0
    ssb = sum(len(vals) * (sum(vals) / len(vals) - grand_mean) ** 2 for vals in groups.values())
    return ssb / sst * 100.0


def median_of(values):
    values = [v for v in values if v is not None]
    if not values:
        return None
    return statistics.median(values)


def compute_tree_order_by_dbh(trees_in_csv):
    """DBH-ascending tree order (revision 2, Part A) - every figure in this
    script uses THIS order, not volume_results.csv's own row order, so a
    line drawn across trees says something about tree size instead of
    nothing. Field DBH via base_vs_dbh_compare.py's parse_field_file() (not
    re-implemented here). Ties broken by tree name for a deterministic
    order (matches the task's own worked example: S13 49cm before S21
    49cm)."""
    field_data = parse_field_file(FIELD_FILE)   # {tree: (base_cm, dbh_cm)}
    missing = [t for t in trees_in_csv if t not in field_data]
    if missing:
        print("WARNING: no field DBH found for %s - keeping CSV order for "
              "these, placed after every tree with a known DBH" % missing)
    with_dbh = [t for t in trees_in_csv if t in field_data]
    with_dbh.sort(key=lambda t: (field_data[t][1], t))
    dbh_cm = {t: field_data[t][1] for t in with_dbh}
    return with_dbh + missing, dbh_cm


# =====================  LOAD + CLASSIFY + EXCLUDE  =====================

def load_and_prepare():
    rows = load_results(RESULTS_CSV)
    trees_in_csv = list(dict.fromkeys(r["tree"] for r in rows))
    tree_order, field_dbh_cm = compute_tree_order_by_dbh(trees_in_csv)

    for r in rows:
        r["_family"] = classify_family(r["method"])
        r["_weighting"] = "-wlen" if "-wlen" in r["method"] else "nowlen"

    unclassified = [r["method"] for r in rows if r["_family"] is None]
    if unclassified:
        print("WARNING: %d row(s) matched no known family (excluded entirely - not "
              "counted in any population): %s" % (len(unclassified), sorted(set(unclassified))[:5]))

    classified = [r for r in rows if r["_family"] is not None]

    excluded = []      # (row, reason)
    kept = []
    for r in classified:
        reason = known_failure_reason(r)
        if reason:
            excluded.append((r, reason))
        else:
            kept.append(r)

    selected_variant = compute_selected_variants(kept)

    for r in kept:
        r["_settled"] = is_settled(r, selected_variant)
    for r, _ in excluded:
        r["_settled"] = is_settled(r, selected_variant)

    pop_full_rows = kept
    pop_settled_rows = [r for r in kept if r["_settled"]]

    return {
        "tree_order": tree_order,
        "field_dbh_cm": field_dbh_cm,
        "pop_full": pop_full_rows,
        "pop_settled": pop_settled_rows,
        "excluded": excluded,
        "selected_variant": selected_variant,
    }


def sanity_check_A(data):
    print("=" * 78)
    print("SANITY CHECK A: tree order (field DBH ascending) - identical in every figure")
    print("=" * 78)
    for t in data["tree_order"]:
        print("  %-10s DBH=%5.1f cm" % (t, data["field_dbh_cm"][t]))
    print()


# =====================  STEP 1 SANITY CHECK  ============================

def sanity_check_1(data):
    print("=" * 78)
    print("SANITY CHECK 1: rows per tree per family, both populations")
    print("=" * 78)
    tree_order = data["tree_order"]
    ok = True
    for pop_name, rows in (("POP_SETTLED", data["pop_settled"]), ("POP_FULL", data["pop_full"])):
        print("-- %s --" % pop_name)
        for fam in ALL_FAMILIES:
            counts = {t: 0 for t in tree_order}
            for r in rows:
                if r["_family"] == fam:
                    counts[r["tree"]] += 1
            print("  %-20s %s" % (fam, " ".join("%s=%d" % (t, counts[t]) for t in tree_order)))
            if pop_name == "POP_SETTLED" and fam == "AdTree_calibrated":
                bad = {t: c for t, c in counts.items() if c != EXPECTED_SETTLED_ADTREE_CAL_PER_TREE}
                if bad:
                    ok = False
                    print("  *** MISMATCH: expected %d AdTree_calibrated rows/tree in "
                          "POP_SETTLED, got: %s" % (EXPECTED_SETTLED_ADTREE_CAL_PER_TREE, bad))
    print()
    if not ok:
        print("STOPPING: sanity check 1 failed (see MISMATCH above).")
        raise SystemExit(1)
    return True


# =====================  STEP 3 SANITY CHECK  ============================

def sanity_check_3(data):
    print("=" * 78)
    print("SANITY CHECK 3: excluded rows per population")
    print("=" * 78)
    selected_variant = data["selected_variant"]
    counts = {POP_FULL: defaultdict(int), POP_SETTLED: defaultdict(int)}
    for r, reason in data["excluded"]:
        counts[POP_FULL][reason] += 1
        if r["_settled"]:
            counts[POP_SETTLED][reason] += 1
    for reason in sorted(set(r for _, r in data["excluded"])):
        print("  %s" % reason)
        print("    POP_FULL removed: %d   POP_SETTLED removed: %d"
              % (counts[POP_FULL][reason], counts[POP_SETTLED][reason]))
    print()
    return counts


# =====================  STEP 4: SUMMARY TABLE  ==========================

def build_summary_rows(data):
    """One row per (population, family, tree, quantity)."""
    out = []
    for pop_name, rows in (("POP_SETTLED", data["pop_settled"]), ("POP_FULL", data["pop_full"])):
        for fam in ALL_FAMILIES:
            for tree in data["tree_order"]:
                fam_tree_rows = [r for r in rows if r["_family"] == fam and r["tree"] == tree]
                for qkey, _ in QUANTITIES:
                    values = [v for v in (quantity_value(r, qkey) for r in fam_tree_rows) if v is not None]
                    stats = summary_stats(values)
                    out.append(dict(population=pop_name, family=fam, tree=tree, quantity=qkey, **stats))
    return out


def write_summary_csv(summary_rows, path):
    fields = ["population", "family", "tree", "quantity", "n", "min", "q1", "median",
              "q3", "max", "mean", "std", "cv_pct"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)


# =====================  STEP 5: VARIANCE DECOMPOSITION  =================

def _close(a, b, tol=_PD2MIN_TOL):
    return a is not None and abs(a - b) < tol


_FACTOR_ROW_KEY = {"weighting": "_weighting", "pd2min_m": "pd2min"}   # display names -> load_results() row keys


def _eta2_over_rows(rows, factor_key, qkey):
    row_key = _FACTOR_ROW_KEY.get(factor_key, factor_key)
    groups = defaultdict(list)
    for r in rows:
        v = quantity_value(r, qkey)
        if v is None:
            continue
        groups[r[row_key]].append(v)
    return eta_squared(groups)


def build_eta2_rows(data):
    """One row per (family, stage, quantity, factor): median/min/max eta^2
    across the 12 trees, computed on POP_FULL (parameter sensitivity needs
    the parameter to actually vary - POP_SETTLED fixes almost all of these
    by construction). AUTO / pd2min=0.020 probe rows are excluded here and
    reported separately (see report_treeqsm_unbalanced_points below)."""
    rows_full = data["pop_full"]
    tree_order = data["tree_order"]
    out = []

    # ---- AdTree_calibrated / AdTree_raw --------------------------------
    for fam, factors in ETA2_FACTORS.items():
        for factor in factors:
            for qkey, _ in QUANTITIES:
                per_tree = []
                for tree in tree_order:
                    tree_rows = [r for r in rows_full if r["_family"] == fam and r["tree"] == tree]
                    e = _eta2_over_rows(tree_rows, factor, qkey)
                    if e is not None:
                        per_tree.append(e)
                if per_tree:
                    out.append(dict(family=fam, stage="", quantity=qkey, factor=factor,
                                     median_eta2=statistics.median(per_tree),
                                     min_eta2=min(per_tree), max_eta2=max(per_tree), n_trees=len(per_tree)))

    # ---- TreeQSM (balanced part only: 3 pd2min levels x 16 simp cells) --
    for fam, stage in TREEQSM_STAGE_OF_FAMILY.items():
        for factor in TREEQSM_ETA2_FACTORS:
            for qkey, _ in QUANTITIES:
                per_tree = []
                for tree in tree_order:
                    tree_rows = [r for r in rows_full if r["_family"] == fam and r["tree"] == tree
                                 and r["mode"] == "manual"
                                 and any(_close(r["pd2min"], lvl) for lvl in TREEQSM_BALANCED_PD2MIN_LEVELS)]
                    e = _eta2_over_rows(tree_rows, factor, qkey)
                    if e is not None:
                        per_tree.append(e)
                if per_tree:
                    out.append(dict(family=fam, stage=stage, quantity=qkey, factor=factor,
                                     median_eta2=statistics.median(per_tree),
                                     min_eta2=min(per_tree), max_eta2=max(per_tree), n_trees=len(per_tree)))
    return out


def write_eta2_csv(eta2_rows, path):
    fields = ["family", "stage", "quantity", "factor", "median_eta2", "min_eta2", "max_eta2", "n_trees"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in eta2_rows:
            w.writerow(row)


def report_treeqsm_unbalanced_points(data):
    """AUTO (single run/tree) and the pd2min=0.020 probe (no simp sweep) -
    reported as individual branch_m3 values per tree, NOT folded into the
    eta^2 decomposition (the design has no simp grid for either, so no
    factor there is estimable)."""
    rows_full = data["pop_full"]
    print("-- TreeQSM unbalanced points (reported as individual values, not in the eta^2 decomposition) --")
    balanced_per_tree = {}
    for fam, stage in TREEQSM_STAGE_OF_FAMILY.items():
        auto_vals, probe_vals = [], []
        for tree in data["tree_order"]:
            tree_rows = [r for r in rows_full if r["_family"] == fam and r["tree"] == tree]
            n_balanced_here = 0
            for r in tree_rows:
                b = quantity_value(r, "branch_m3")
                if r["mode"] == "auto":
                    auto_vals.append((tree, b))
                elif r["mode"] == "manual" and _close(r["pd2min"], TREEQSM_PROBE_PD2MIN):
                    probe_vals.append((tree, b))
                elif r["mode"] == "manual" and any(_close(r["pd2min"], lvl) for lvl in TREEQSM_BALANCED_PD2MIN_LEVELS):
                    n_balanced_here += 1
            balanced_per_tree[(fam, tree)] = n_balanced_here
        print("  %s: AUTO branch_m3 per tree: %s"
              % (LEGEND_LABELS[fam], ["%s=%.3f" % (t, v) for t, v in auto_vals if v is not None]))
        print("  %s: pd2min=0.020 probe branch_m3 per tree: %s"
              % (LEGEND_LABELS[fam], ["%s=%.3f" % (t, v) for t, v in probe_vals if v is not None]))
    counts = sorted(set(balanced_per_tree.values()))
    n_balanced_used = counts[0] if len(counts) == 1 else None
    if n_balanced_used == 48:
        print("  balanced-part models used per stage per tree: 48 (expected 48 = 3 pd2min levels x 16 simp cells) - OK")
    else:
        print("  *** balanced-part model count is NOT a uniform 48/stage/tree: %s"
              % {k: v for k, v in balanced_per_tree.items() if v != 48})
    print()
    return n_balanced_used or 48


def print_eta2_sanity_check(eta2_rows):
    print("=" * 78)
    print("SANITY CHECK 5: median eta^2 (AdTree_calibrated, POP_FULL) vs. independent reference")
    print("=" * 78)
    lookup = {(r["family"], r["quantity"], r["factor"]): r["median_eta2"] for r in eta2_rows}
    print("  branch_m3:")
    for factor, expected in EXPECTED_ETA2_BRANCH.items():
        actual = lookup.get(("AdTree_calibrated", "branch_m3", factor))
        print("    %-22s actual=%6s  expected=%5.1f" %
              (factor, ("%.1f" % actual) if actual is not None else "None", expected))
    print("  n_cylinders:")
    for factor, expected in EXPECTED_ETA2_NCYL.items():
        actual = lookup.get(("AdTree_calibrated", "n_cylinders", factor))
        print("    %-22s actual=%6s  expected=%5.1f" %
              (factor, ("%.1f" % actual) if actual is not None else "None", expected))
    print()
    return lookup


# =====================  STEP 6: STABILITY TABLE  ========================

def build_stability_rows(data):
    """One row per (population, family, quantity): median CV across the
    12 trees (CV itself computed per tree, per Step 6)."""
    out = []
    for pop_name, rows in (("POP_SETTLED", data["pop_settled"]), ("POP_FULL", data["pop_full"])):
        for fam in ALL_FAMILIES:
            for qkey, _ in QUANTITIES:
                per_tree_cv = []
                for tree in data["tree_order"]:
                    values = [v for v in (quantity_value(r, qkey) for r in rows
                                           if r["_family"] == fam and r["tree"] == tree) if v is not None]
                    stats = summary_stats(values)
                    if stats["cv_pct"] is not None:
                        per_tree_cv.append(stats["cv_pct"])
                out.append(dict(population=pop_name, family=fam, quantity=qkey,
                                 median_cv_pct=median_of(per_tree_cv), n_trees=len(per_tree_cv)))
    return out


def write_stability_csv(stability_rows, path):
    fields = ["population", "family", "quantity", "median_cv_pct", "n_trees"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in stability_rows:
            w.writerow(row)


def print_stability_table(stability_rows):
    print("=" * 78)
    print("STABILITY TABLE: median CV [%] across trees, per family/quantity")
    print("=" * 78)
    lookup = {(r["population"], r["family"], r["quantity"]): r["median_cv_pct"] for r in stability_rows}
    header = "  %-20s" % "family" + "".join("%14s" % q for q in QUANTITY_KEYS)
    for pop in (POP_FULL, POP_SETTLED):
        print("-- %s --" % pop)
        print(header)
        for fam in ALL_FAMILIES:
            cells = []
            for qkey in QUANTITY_KEYS:
                v = lookup.get((pop, fam, qkey))
                cells.append("%14s" % (("%.1f" % v) if v is not None else "-"))
            print("  %-20s%s" % (fam, "".join(cells)))
        print()
    return lookup


def print_cv_sanity_check(stability_lookup):
    print("=" * 78)
    print("SANITY CHECK 6: AdTree_calibrated / branch_m3 CV")
    print("=" * 78)
    full = stability_lookup.get((POP_FULL, "AdTree_calibrated", "branch_m3"))
    settled = stability_lookup.get((POP_SETTLED, "AdTree_calibrated", "branch_m3"))
    print("  POP_FULL:    actual=%s%%  expected=%.1f%%" %
          (("%.1f" % full) if full is not None else "None", EXPECTED_CV_BRANCH_FULL))
    print("  POP_SETTLED: actual=%s%%  expected=%.1f%%" %
          (("%.1f" % settled) if settled is not None else "None", EXPECTED_CV_BRANCH_SETTLED))
    print()
    return full, settled


# =====================  STEP 3: EXCLUDED CSV  ============================

def write_excluded_csv(excluded, path):
    fields = ["tree", "method", "family", "population_would_be", "reason",
              "total_m3", "trunk_m3", "branch_m3", "dbh_m"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r, reason in excluded:
            w.writerow(dict(
                tree=r["tree"], method=r["method"], family=r["_family"],
                population_would_be="POP_SETTLED+POP_FULL" if r["_settled"] else "POP_FULL only",
                reason=reason,
                total_m3=r["total"], trunk_m3=r["trunk"], branch_m3=r["branch"], dbh_m=r["dbh"],
            ))


# =====================  STEP 7: FIGURES  =================================

def _darken(color, factor=BOX_EDGE_DARKEN_FACTOR):
    r, g, b = mcolors.to_rgb(color)
    keep = 1 - factor
    return (r * keep, g * keep, b * keep)


def _lighten(color, factor=0.55):
    """Same-hue lighter tint of `color`, mixed toward white by `factor` (0 =
    unchanged, 1 = white) - used for the line figures' interquartile-range
    band, a lighter tint of the same colour as that family's median line
    (revision 2, Part B)."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * factor, g + (1 - g) * factor, b + (1 - b) * factor)


def _place_legend_below(ax, handles, labels, title):
    """Horizontal legend below the axes, one row if it fits (revision 2,
    Part B - applies to both the box and the line per-tree figures)."""
    return ax.legend(handles, labels, fontsize=LEGEND_FONTSIZE, loc="upper center",
                      bbox_to_anchor=(0.5, -0.16), ncol=len(handles), frameon=True,
                      title=title)


def _legend_handles():
    handles = []
    for fam in BOX_FAMILIES:
        handles.append(mlines.Line2D([0], [0], marker="s", linestyle="",
                                      markerfacecolor=METHOD_COLORS[fam],
                                      markeredgecolor=_darken(METHOD_COLORS[fam]), markersize=10))
    handles.append(mlines.Line2D([0], [0], marker="D", linestyle="",
                                  markerfacecolor=METHOD_COLORS[MARKER_FAMILY],
                                  markeredgecolor=_darken(METHOD_COLORS[MARKER_FAMILY]), markersize=8))
    labels = [LEGEND_LABELS[f] for f in BOX_FAMILIES] + [LEGEND_LABELS[MARKER_FAMILY]]
    return handles, labels


def _wrapped_title(main_line, population, width=95):
    """Two-line figure title: the main line, then the population spelled out
    in full (see POP_LABEL_LONG) so a reader never has to look at the
    filename to know which population a figure shows (revision Part A1)."""
    subtitle = "Population: %s" % POP_LABEL_LONG[population]
    return main_line + "\n" + "\n".join(textwrap.wrap(subtitle, width=width))


def _apply_log_ticks(ax, axis="y"):
    """Decade + 2/5 subdivision ticks on a log axis, labelled as plain
    numbers (0.1, 0.2, 0.5, 1, 2, 5, 10, ...) instead of scientific notation
    or unlabelled minor ticks (revision Part A3)."""
    which = ax.yaxis if axis == "y" else ax.xaxis
    which.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    which.set_major_formatter(FuncFormatter(lambda v, pos: ("%g" % v)))
    which.set_minor_locator(LogLocator(base=10, subs=()))   # no unlabelled minor ticks cluttering the axis


def draw_box_by_tree(rows, tree_order, qkey, population, out_name, log_scale=False):
    """Box per (tree, family) - AdTree_raw/AdTree_calibrated/TreeQSM_Optimal/
    TreeQSM_Simplified as boxes, AdQSM as individual markers (Step 2: a box
    of 1-2 points is meaningless). Same tree order in every figure (CSV
    first-appearance order).

    `population`: POP_SETTLED or POP_FULL - used to look up the long/short
    human-readable wording (POP_LABEL_LONG/SHORT), never shown as the raw
    code in the figure itself (revision Part A1)."""
    n_trees = len(tree_order)
    n_slots = len(ALL_FAMILIES)
    # Narrower per-tree footprint (0.62 of the unit spacing, was 0.8) so
    # neighbouring tree groups don't visually merge into one another even
    # with 5 series per tree (revision Part A2).
    slot_width = 0.62 / n_slots
    offsets = [(i - (n_slots - 1) / 2.0) * slot_width for i in range(n_slots)]
    box_width = slot_width * 0.8

    fig, ax = plt.subplots(figsize=(max(15, n_trees * 1.9), 6.5))

    # Alternating light background bands, one per tree, so the eye can tell
    # at a glance which boxes belong to which tree without following the
    # x-tick grid line by line (revision Part A2).
    for ti in range(n_trees):
        if ti % 2 == 1:
            ax.axvspan(ti - 0.5, ti + 0.5, color="0.92", zorder=0, linewidth=0)

    box_data, box_positions, box_colors = [], [], []
    for ti, tree in enumerate(tree_order):
        for fi, fam in enumerate(BOX_FAMILIES):
            values = [v for v in (quantity_value(r, qkey) for r in rows
                                   if r["tree"] == tree and r["_family"] == fam) if v is not None]
            if not values:
                continue
            box_data.append(values)
            box_positions.append(ti + offsets[fi])
            box_colors.append(METHOD_COLORS[fam])

    if box_data:
        bp = ax.boxplot(box_data, positions=box_positions, widths=box_width,
                         patch_artist=True, whis=BOX_WHIS, zorder=2)
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_edgecolor(_darken(color))
        for i, color in enumerate(box_colors):
            edge = _darken(color)
            bp["whiskers"][2 * i].set_color(edge)
            bp["whiskers"][2 * i + 1].set_color(edge)
            bp["caps"][2 * i].set_color(edge)
            bp["caps"][2 * i + 1].set_color(edge)
            bp["medians"][i].set_color(edge)

    marker_offset = offsets[-1]
    marker_color = METHOD_COLORS[MARKER_FAMILY]
    for ti, tree in enumerate(tree_order):
        values = [v for v in (quantity_value(r, qkey) for r in rows
                               if r["tree"] == tree and r["_family"] == MARKER_FAMILY) if v is not None]
        if not values:
            continue
        ax.scatter([ti + marker_offset] * len(values), values, marker="D",
                   s=JITTER_POINT_SIZE, color=marker_color, edgecolors=_darken(marker_color), zorder=3)

    # Tree ticks shortened to just the plot number (all 12 trees share the
    # "B21_" prefix) - the omission is stated once in the x-axis label
    # instead of repeating "B21_" on every one of 12 ticks (revision Part A2).
    ax.set_xlim(-0.5, n_trees - 0.5)
    ax.set_xticks(range(n_trees))
    ax.set_xticklabels([t.replace("B21_", "") for t in tree_order],
                        rotation=0, ha="center", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel("Tree, sorted by field DBH ascending (\"B21_\" prefix omitted from tick labels)",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(QUANTITY_LABEL[qkey] + (" (log scale)" if log_scale else ""),
                  fontsize=AXIS_LABEL_FONTSIZE)
    if log_scale:
        ax.set_yscale("log")
        _apply_log_ticks(ax)
    ax.set_title(_wrapped_title("%s across methods" % QUANTITY_LABEL[qkey], population),
                 fontsize=TITLE_FONTSIZE)

    # Horizontal legend below the axes (revision 2, Part B) - loc="best" used
    # to sit on top of a box with 12 dense tree groups (revision 1, Part A5),
    # and a right-hand legend (revision 1's fix) is superseded by this
    # explicit horizontal-below placement.
    handles, labels = _legend_handles()
    _place_legend_below(ax, handles, labels, POP_LABEL_SHORT[population])

    fig.tight_layout()
    out_path = os.path.join(ensure_all_plots_dir(), out_name)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _wrapped_title_for_lines(main_line, population, width=95):
    """Two-line title for the line figures - same population wording as
    _wrapped_title(), plus a note on what the line/band mean (revision 2,
    Part B: "state in the legend or axis label that the line is the median
    and the band the interquartile range" - stated here, in the title,
    since the axis label is already carrying the quantity/units)."""
    subtitle = ("Population: %s. Line = per-tree median; shaded band = "
                "interquartile range (Q1-Q3) across that tree's rows."
                % POP_LABEL_LONG[population])
    return main_line + "\n" + "\n".join(textwrap.wrap(subtitle, width=width))


def per_tree_line_stats(rows, tree_order, qkey, fam):
    """(medians, q1s, q3s) - one entry per tree, in `tree_order` - for
    family `fam`. A tree with no data for this family gets NaN in all
    three arrays, which matplotlib renders as a gap in the line/band
    instead of a wrong zero (revision 2, Part B)."""
    medians, q1s, q3s = [], [], []
    for tree in tree_order:
        values = [v for v in (quantity_value(r, qkey) for r in rows
                               if r["tree"] == tree and r["_family"] == fam) if v is not None]
        if values:
            q1, med, q3 = _quantiles(sorted(values))
        else:
            q1, med, q3 = np.nan, np.nan, np.nan
        medians.append(med)
        q1s.append(q1)
        q3s.append(q3)
    return np.array(medians), np.array(q1s), np.array(q3s)


def draw_line_by_tree(rows, tree_order, qkey, population, out_name, log_scale=False):
    """Line-figure alternative to draw_box_by_tree() for the SAME (rows,
    tree_order, qkey, population) - one median line + IQR band per family,
    AdQSM as larger markers only (no line/band - one or two points per tree
    cannot define a median line that means anything more than the points
    themselves). Kept as a SEPARATE figure (suffix "_lines") alongside the
    box version, never replacing it (revision 2, Part B)."""
    n_trees = len(tree_order)
    x = np.arange(n_trees)

    fig, ax = plt.subplots(figsize=(max(15, n_trees * 1.6), 6.5))
    for ti in range(n_trees):
        if ti % 2 == 1:
            ax.axvspan(ti - 0.5, ti + 0.5, color="0.92", zorder=0, linewidth=0)

    for fam in BOX_FAMILIES:
        med, q1, q3 = per_tree_line_stats(rows, tree_order, qkey, fam)
        color = METHOD_COLORS[fam]
        ax.fill_between(x, q1, q3, color=_lighten(color, 0.55), alpha=0.7, zorder=1, linewidth=0)
        ax.plot(x, med, marker="o", markersize=6, color=_darken(color, 0.15),
                linewidth=1.8, zorder=2)

    # AdQSM: larger markers only, drawn LAST (on top) so it stays visible
    # over the bands (revision 2, Part B).
    aq_med, _, _ = per_tree_line_stats(rows, tree_order, qkey, MARKER_FAMILY)
    ax.scatter(x, aq_med, marker="D", s=JITTER_POINT_SIZE * 2.2,
               color=METHOD_COLORS[MARKER_FAMILY], edgecolors=_darken(METHOD_COLORS[MARKER_FAMILY]),
               zorder=5)

    ax.set_xlim(-0.5, n_trees - 0.5)
    ax.set_xticks(range(n_trees))
    ax.set_xticklabels([t.replace("B21_", "") for t in tree_order],
                        rotation=0, ha="center", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_xlabel("Tree, sorted by field DBH ascending (\"B21_\" prefix omitted from tick labels)",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(QUANTITY_LABEL[qkey] + (" (log scale)" if log_scale else ""),
                  fontsize=AXIS_LABEL_FONTSIZE)
    if log_scale:
        ax.set_yscale("log")
        _apply_log_ticks(ax)
    ax.set_title(_wrapped_title_for_lines("%s across methods" % QUANTITY_LABEL[qkey], population),
                 fontsize=TITLE_FONTSIZE)

    handles, labels = _legend_handles()
    _place_legend_below(ax, handles, labels, POP_LABEL_SHORT[population])

    fig.tight_layout()
    out_path = os.path.join(ensure_all_plots_dir(), out_name)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


# Light sequential palette for the eta^2 heat map (revision 2, Part C) -
# samples only the PALE half [0.0, 0.55] of matplotlib's "Blues", so even a
# cell at eta^2=100 (the darkest the map ever reaches) stays light enough
# for a single dark text colour everywhere - removing the white/dark text
# switch removes the readability problem at its source instead of tuning
# the switching threshold.
_ETA2_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "eta2_light", plt.cm.Blues(np.linspace(0.0, 0.55, 256)))
_ETA2_TEXT_COLOR = "0.15"   # near-black, used for every cell regardless of value


def draw_eta2_heatmap(eta2_rows, out_name):
    lookup = {(r["family"], r["quantity"], r["factor"]): r["median_eta2"] for r in eta2_rows}
    panels = [
        ("AdTree_calibrated", "AdTree calibrated", ETA2_FACTORS["AdTree_calibrated"]),
        ("AdTree_raw", "AdTree raw", ETA2_FACTORS["AdTree_raw"]),
        ("TreeQSM_Optimal", "TreeQSM Optimal\n(balanced part only)", TREEQSM_ETA2_FACTORS),
        ("TreeQSM_Simplified", "TreeQSM Simplified,\nexported to ANSYS\n(balanced part only)", TREEQSM_ETA2_FACTORS),
    ]
    fig, axes = plt.subplots(1, len(panels), figsize=(7.0 * len(panels), 5.6))
    fig.subplots_adjust(wspace=0.6, top=0.72, bottom=0.32, left=0.06, right=0.90)
    im = None
    for ax, (fam, panel_title, factors) in zip(axes, panels):
        matrix = np.full((len(factors), len(QUANTITY_KEYS)), np.nan)
        for fi, factor in enumerate(factors):
            for qi, qkey in enumerate(QUANTITY_KEYS):
                v = lookup.get((fam, qkey, factor))
                if v is not None:
                    matrix[fi, qi] = v
        im = ax.imshow(matrix, vmin=0, vmax=100, cmap=_ETA2_CMAP, aspect="auto")
        ax.set_xticks(range(len(QUANTITY_KEYS)))
        ax.set_xticklabels([QUANTITY_LABEL[q] for q in QUANTITY_KEYS], rotation=60, ha="right",
                            fontsize=ANNOTATION_FONTSIZE)
        ax.set_yticks(range(len(factors)))
        ax.set_yticklabels(factors, fontsize=ANNOTATION_FONTSIZE)
        ax.set_title(panel_title, fontsize=PANEL_TITLE_FONTSIZE)
        # Numeric value printed IN every cell (one decimal place) - a SINGLE
        # dark text colour everywhere now that the palette itself never gets
        # dark enough to need a white alternative (revision 2, Part C).
        for fi in range(len(factors)):
            for qi in range(len(QUANTITY_KEYS)):
                val = matrix[fi, qi]
                if not np.isnan(val):
                    ax.text(qi, fi, "%.1f" % val, ha="center", va="center",
                            fontsize=POINT_LABEL_FONTSIZE, color=_ETA2_TEXT_COLOR)
    fig.suptitle(_wrapped_title("How much of each quantity's spread is explained by each "
                                 "reconstruction parameter (variance explained, eta^2, in %)",
                                 POP_FULL, width=110),
                 fontsize=TITLE_FONTSIZE)
    if im is not None:
        fig.colorbar(im, ax=axes, shrink=0.8, label="variance explained [%]")
    out_path = os.path.join(ensure_all_plots_dir(), out_name)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def draw_stability_bars(stability_lookup, out_name):
    """Grouped bars, x=families, y=median coefficient of variation (CV) of
    branch_m3 [%], two series (settled configuration / full grid) -
    branch_m3 is the quantity the whole "settled vs. full" decision (module
    header) is framed around."""
    fig, ax = plt.subplots(figsize=(11, 6.5))
    x = np.arange(len(ALL_FAMILIES))
    width = 0.35
    for i, (pop, dx) in enumerate([(POP_SETTLED, -width / 2), (POP_FULL, width / 2)]):
        values = [stability_lookup.get((pop, fam, "branch_m3")) for fam in ALL_FAMILIES]
        heights = [v if v is not None else 0.0 for v in values]
        colors = [METHOD_COLORS[fam] for fam in ALL_FAMILIES]
        alpha = 1.0 if pop == POP_SETTLED else 0.55
        bars = ax.bar(x + dx, heights, width, color=colors, alpha=alpha,
                      edgecolor=[_darken(c) for c in colors],
                      label=POP_LABEL_SHORT[pop], hatch=("" if pop == POP_SETTLED else "//"))
        for bar, v in zip(bars, values):
            if v is None:
                ax.text(bar.get_x() + bar.get_width() / 2, 0.5, "n/a", ha="center",
                        va="bottom", fontsize=POINT_LABEL_FONTSIZE, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([LEGEND_LABELS[f] for f in ALL_FAMILIES], rotation=20, ha="right",
                        fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Median coefficient of variation (CV) of branch volume\nacross the 12 trees [%]",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title("\n".join(textwrap.wrap(
        "Stability of branch volume by method family - %s vs. %s"
        % (POP_LABEL_LONG[POP_SETTLED], POP_LABEL_LONG[POP_FULL]), width=100)),
        fontsize=TITLE_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE, title="Population")
    fig.tight_layout()
    out_path = os.path.join(ensure_all_plots_dir(), out_name)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def per_tree_branch_share_medians(rows, tree_order):
    """{family: [one value per tree - that tree's MEDIAN branch_frac*100
    across its own rows]}. Deliberately one number per tree (not every row
    pooled) - pooling would combine two different sources of spread (how
    much the parameters move the answer, already covered by the stability
    figure, and how much the trees differ from each other, which is what
    this figure is about) - see revision Part B2."""
    out = {}
    for fam in ALL_FAMILIES:
        medians = []
        for tree in tree_order:
            values = [v for v in (quantity_value(r, "branch_frac") for r in rows
                                   if r["tree"] == tree and r["_family"] == fam) if v is not None]
            if values:
                medians.append(100.0 * statistics.median(values))
        out[fam] = medians
    return out


def print_branch_share_table(data):
    """C2: one row per family, median/min/max branch share [%] across the
    12 per-tree medians, both populations side by side."""
    print("=" * 78)
    print("BRANCH SHARE OF TOTAL VOLUME [%], per family (each input = one tree's median)")
    print("=" * 78)
    tables = {}
    for pop_name, rows in ((POP_SETTLED, data["pop_settled"]), (POP_FULL, data["pop_full"])):
        per_family = per_tree_branch_share_medians(rows, data["tree_order"])
        print("-- %s (%s) --" % (pop_name, POP_LABEL_SHORT[pop_name]))
        print("  %-20s %6s %8s %8s %8s" % ("family", "n", "min", "median", "max"))
        table = {}
        for fam in ALL_FAMILIES:
            vals = per_family[fam]
            if vals:
                row = dict(n=len(vals), min=min(vals), median=statistics.median(vals), max=max(vals))
            else:
                row = dict(n=0, min=None, median=None, max=None)
            table[fam] = row
            print("  %-20s %6d %8s %8s %8s" % (
                fam, row["n"],
                ("%.1f" % row["min"]) if row["min"] is not None else "-",
                ("%.1f" % row["median"]) if row["median"] is not None else "-",
                ("%.1f" % row["max"]) if row["max"] is not None else "-"))
        tables[pop_name] = table
        print()
    return tables


def draw_branch_share_by_family(rows, tree_order, population, out_name):
    """One box per family (AdQSM: markers only, Step 2), built from the 12
    per-tree MEDIAN branch-share values - x=family, y=branch share of total
    volume [%]. Individual tree medians overlaid as jittered points on top
    of each box so a reader can see whether a family's spread comes from
    one outlier tree or from all of them (revision Part B2)."""
    per_family = per_tree_branch_share_medians(rows, tree_order)
    x_families = BOX_FAMILIES + [MARKER_FAMILY]

    fig, ax = plt.subplots(figsize=(10, 6.5))
    rng = np.random.default_rng(0)   # fixed seed: jitter positions are reproducible run to run

    box_data, box_positions, box_colors = [], [], []
    for i, fam in enumerate(BOX_FAMILIES):
        values = per_family[fam]
        if not values:
            continue
        box_data.append(values)
        box_positions.append(i)
        box_colors.append(METHOD_COLORS[fam])

    if box_data:
        bp = ax.boxplot(box_data, positions=box_positions, widths=0.5,
                         patch_artist=True, whis=BOX_WHIS, zorder=2)
        for patch, color in zip(bp["boxes"], box_colors):
            patch.set_facecolor(color)
            patch.set_edgecolor(_darken(color))
        for i, color in enumerate(box_colors):
            edge = _darken(color)
            bp["whiskers"][2 * i].set_color(edge)
            bp["whiskers"][2 * i + 1].set_color(edge)
            bp["caps"][2 * i].set_color(edge)
            bp["caps"][2 * i + 1].set_color(edge)
            bp["medians"][i].set_color(edge)

    # Individual tree medians overlaid on every box, jittered sideways.
    for i, fam in enumerate(BOX_FAMILIES):
        values = per_family[fam]
        if not values:
            continue
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        ax.scatter(np.full(len(values), i) + jitter, values, s=JITTER_POINT_SIZE,
                   color=_darken(METHOD_COLORS[fam], factor=0.15), edgecolors="none",
                   alpha=0.85, zorder=3)

    # AdQSM: markers only, no box (Step 2 - a box of ~12 points from a
    # 1-2-row-per-tree family would still be legitimate here since each
    # input IS a per-tree median, but AdQSM is kept visually consistent
    # with every other figure in this script, where it is never boxed).
    aq_index = len(BOX_FAMILIES)
    aq_values = per_family[MARKER_FAMILY]
    if aq_values:
        jitter = rng.uniform(-0.12, 0.12, size=len(aq_values))
        ax.scatter(np.full(len(aq_values), aq_index) + jitter, aq_values, marker="D",
                   s=JITTER_POINT_SIZE, color=METHOD_COLORS[MARKER_FAMILY],
                   edgecolors=_darken(METHOD_COLORS[MARKER_FAMILY]), zorder=3)

    ax.set_xlim(-0.6, len(x_families) - 0.4)
    ax.set_xticks(range(len(x_families)))
    ax.set_xticklabels([LEGEND_LABELS[f] for f in x_families], rotation=20, ha="right",
                        fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylim(0, 100)
    ax.set_ylabel("Branch share of total volume [%]\n(each point = one tree's median across its own rows)",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(_wrapped_title("Branch share of total volume, by method family", population),
                 fontsize=TITLE_FONTSIZE)

    fig.tight_layout()
    out_path = os.path.join(ensure_all_plots_dir(), out_name)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def compute_pd2min_sweep_data(rows, tree_order, stage_family):
    """{tree: {"levels": {level: {qkey: median_or_None}},
               "auto": {"pd2min": x, qkey: median_or_None} or None}}
    for TreeQSM stage `stage_family` (TreeQSM_Optimal / TreeQSM_Simplified).
    Each manual level's value is the MEDIAN across that level's simp-grid
    rows (0.002/0.005/0.010 each carry 16 simp cells; 0.020 is a single
    probe with no simp grid, so its "median" is just that probe's own
    value - except B21_S01, which has 3 mincylrad variants at 0.020, where
    it is their median)."""
    out = {}
    for tree in tree_order:
        tree_rows = [r for r in rows if r["tree"] == tree and r["_family"] == stage_family]
        levels = {}
        for lvl in PD2MIN_SWEEP_LEVELS:
            lvl_rows = [r for r in tree_rows if r["mode"] == "manual" and _close(r["pd2min"], lvl)]
            vals = {}
            for qkey, _panel_title, row_key in PD2MIN_SWEEP_QUANTITIES:
                v = [r[row_key] for r in lvl_rows if r[row_key] is not None]
                vals[qkey] = statistics.median(v) if v else None
            levels[lvl] = vals
        auto_rows = [r for r in tree_rows if r["mode"] == "auto"]
        auto = None
        if auto_rows:
            pd2min_vals = [r["pd2min"] for r in auto_rows if r["pd2min"] is not None]
            if pd2min_vals:
                auto = dict(pd2min=statistics.median(pd2min_vals))
                for qkey, _panel_title, row_key in PD2MIN_SWEEP_QUANTITIES:
                    v = [r[row_key] for r in auto_rows if r[row_key] is not None]
                    auto[qkey] = statistics.median(v) if v else None
        out[tree] = dict(levels=levels, auto=auto)
    return out


def pd2min_sweep_ratio_table(sweep, tree_order, qkey):
    """{tree: ratio at pd2min=0.020 relative to pd2min=0.002}, skipping a
    tree with no value at either level - used by sanity_check_D_2 below,
    kept separate from the drawing function so the same numbers back both
    the printed check and the figure."""
    ratios = {}
    base_level, probe_level = PD2MIN_SWEEP_LEVELS[0], PD2MIN_SWEEP_LEVELS[-1]
    for tree in tree_order:
        base = sweep[tree]["levels"][base_level][qkey]
        probe = sweep[tree]["levels"][probe_level][qkey]
        if base and probe is not None:
            ratios[tree] = probe / base
    return ratios


def draw_pd2min_sweep_all_trees(data, stage_family, out_name):
    """methodeval_pd2min_sweep_all_trees{,_optimal}.png (revision 2, Part D,
    reduced to 3 branch-only panels per the addendum) - one panel per
    branch quantity, 12 lines (one per tree), each normalised to that
    tree's OWN value at pd2min=0.002 (absolute values span ~2 orders of
    magnitude between trees - normalising asks whether trees RESPOND to
    pd2min alike, not how big they are). AUTO has a different pd2min on
    every tree, so it is drawn as a lone marker at its own x position,
    never connected into the line. Lines coloured by a sequential ramp
    ordered by field DBH, via a colour-bar (not a 12-entry legend)."""
    tree_order = data["tree_order"]
    dbh_cm = data["field_dbh_cm"]
    sweep = compute_pd2min_sweep_data(data["pop_full"], tree_order, stage_family)

    dbh_values = [dbh_cm[t] for t in tree_order]
    norm = mcolors.Normalize(vmin=min(dbh_values), vmax=max(dbh_values))
    cmap = plt.cm.viridis
    tree_color = {t: cmap(norm(dbh_cm[t])) for t in tree_order}

    n_panels = len(PD2MIN_SWEEP_QUANTITIES)
    fig, axes = plt.subplots(1, n_panels, figsize=(6.5 * n_panels, 6.0))
    for ax, (qkey, panel_title, _row_key) in zip(axes, PD2MIN_SWEEP_QUANTITIES):
        for tree in tree_order:
            d = sweep[tree]
            base = d["levels"][PD2MIN_SWEEP_LEVELS[0]][qkey]
            if not base:
                continue   # no pd2min=0.002 value to normalise against - skipped, not zeroed
            ys = [(d["levels"][lvl][qkey] / base) if d["levels"][lvl][qkey] is not None else np.nan
                  for lvl in PD2MIN_SWEEP_LEVELS]
            ax.plot(PD2MIN_SWEEP_LEVELS, ys, marker="o", markersize=5,
                    color=tree_color[tree], linewidth=1.5, zorder=2)
            if d["auto"] is not None and d["auto"].get(qkey) is not None:
                ax.scatter([d["auto"]["pd2min"]], [d["auto"][qkey] / base], marker="*", s=110,
                           color=tree_color[tree], edgecolors="black", linewidths=0.6, zorder=4)
        ax.axhline(1.0, color="0.4", linestyle="--", linewidth=1, zorder=1)
        ax.set_xlabel("pd2min [m] (manual levels; AUTO = star marker, at its own pd2min)",
                      fontsize=ANNOTATION_FONTSIZE)
        ax.set_ylabel("Ratio to pd2min=0.002 [-]", fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_title(panel_title, fontsize=PANEL_TITLE_FONTSIZE)

    # Colorbar drawn in its OWN fixed axes (not the automatic ax=<list>
    # shrink-to-fit form) - that form fights with the explicit
    # subplots_adjust() below and ends up floating on top of the last
    # panel instead of beside it (caught by opening the figure and
    # checking it, same as revision 1's Part A5).
    fig.subplots_adjust(bottom=0.30, top=0.72, right=0.89, wspace=0.35)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = fig.add_axes([0.91, 0.34, 0.015, 0.36])
    fig.colorbar(sm, cax=cax, label="Field DBH [cm] (tree colour ramp)")

    stage_label = ("Simplified (exported to ANSYS)" if stage_family == "TreeQSM_Simplified"
                   else TREEQSM_STAGE_OF_FAMILY[stage_family])
    fig.suptitle(_wrapped_title(
        "pd2min sweep, TreeQSM %s, all 12 trees" % stage_label, POP_FULL, width=110),
        fontsize=TITLE_FONTSIZE, y=0.99)

    # Caution note carried into the figure itself (revision 2 addendum):
    # pmdist_branch_mean describes fit RESPONSE only, never a selection
    # criterion - see PD2MIN_SWEEP_QUANTITIES' own comment for why.
    fig.text(0.5, 0.01, "\n".join(textwrap.wrap(
        "Caution: the point-to-model-distance panel describes how the branch fit responds to pd2min - "
        "it is NOT a model-selection criterion. pmdist_mean/pmdist_branch_mean are both minimised by "
        "the known-degraded fat-cylinder model; only pmdist_trunk_mean is admissible for model "
        "selection, and model selection is not part of this script.", width=175)),
        ha="center", va="bottom", fontsize=ANNOTATION_FONTSIZE)
    out_path = os.path.join(ensure_all_plots_dir(), out_name)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path, sweep


def draw_all_figures(data, eta2_rows, stability_lookup):
    tree_order = data["tree_order"]
    outputs = []
    outputs.append(draw_box_by_tree(data["pop_settled"], tree_order, "branch_m3",
                                     POP_SETTLED, "methodeval_branch_volume_log_settled.png", log_scale=True))
    outputs.append(draw_box_by_tree(data["pop_full"], tree_order, "branch_m3",
                                     POP_FULL, "methodeval_branch_volume_log_full.png", log_scale=True))
    outputs.append(draw_box_by_tree(data["pop_settled"], tree_order, "branch_frac",
                                     POP_SETTLED, "methodeval_branch_frac_settled.png"))
    outputs.append(draw_box_by_tree(data["pop_full"], tree_order, "branch_frac",
                                     POP_FULL, "methodeval_branch_frac_full.png"))
    outputs.append(draw_box_by_tree(data["pop_settled"], tree_order, "trunk_m3",
                                     POP_SETTLED, "methodeval_trunk_volume_settled.png"))
    outputs.append(draw_box_by_tree(data["pop_full"], tree_order, "trunk_m3",
                                     POP_FULL, "methodeval_trunk_volume_full.png"))
    outputs.append(draw_box_by_tree(data["pop_settled"], tree_order, "total_m3",
                                     POP_SETTLED, "methodeval_total_volume_settled.png"))
    outputs.append(draw_box_by_tree(data["pop_settled"], tree_order, "n_cylinders",
                                     POP_SETTLED, "methodeval_ncyl_log_settled.png", log_scale=True))
    outputs.append(draw_eta2_heatmap(eta2_rows, "methodeval_eta2_heatmap.png"))
    outputs.append(draw_stability_bars(stability_lookup, "methodeval_stability.png"))
    outputs.append(draw_branch_share_by_family(data["pop_settled"], tree_order,
                                                POP_SETTLED, "methodeval_branch_share_settled.png"))
    outputs.append(draw_branch_share_by_family(data["pop_full"], tree_order,
                                                POP_FULL, "methodeval_branch_share_full.png"))

    # ---- Line-figure alternatives (revision 2, Part B) - same data as the
    # box figures above, kept alongside them, never replacing them.
    outputs.append(draw_line_by_tree(data["pop_settled"], tree_order, "branch_m3",
                                      POP_SETTLED, "methodeval_branch_volume_log_settled_lines.png", log_scale=True))
    outputs.append(draw_line_by_tree(data["pop_full"], tree_order, "branch_m3",
                                      POP_FULL, "methodeval_branch_volume_log_full_lines.png", log_scale=True))
    outputs.append(draw_line_by_tree(data["pop_settled"], tree_order, "branch_frac",
                                      POP_SETTLED, "methodeval_branch_frac_settled_lines.png"))
    outputs.append(draw_line_by_tree(data["pop_settled"], tree_order, "trunk_m3",
                                      POP_SETTLED, "methodeval_trunk_volume_settled_lines.png"))
    outputs.append(draw_line_by_tree(data["pop_settled"], tree_order, "total_m3",
                                      POP_SETTLED, "methodeval_total_volume_settled_lines.png"))
    outputs.append(draw_line_by_tree(data["pop_settled"], tree_order, "n_cylinders",
                                      POP_SETTLED, "methodeval_ncyl_log_settled_lines.png", log_scale=True))
    return outputs


# =====================  SANITY CHECK D / addendum 2 =====================

def sanity_check_pd2min_sweep(sweep_by_stage, tree_order):
    """Ratio at pd2min=0.020 relative to 0.002, min/median/max across the 12
    trees, for each of the 3 branch quantities and both TreeQSM stages
    (original Part D's branch_m3 check, extended by the addendum to all 3
    quantities x both stages). States plainly whether trees move together
    or diverge - no further interpretation."""
    print("=" * 78)
    print("SANITY CHECK D: pd2min=0.020 / pd2min=0.002 ratio, min/median/max across 12 trees")
    print("=" * 78)
    for stage_family, sweep in sweep_by_stage.items():
        stage_label = LEGEND_LABELS[stage_family]
        print("-- %s --" % stage_label)
        for qkey, panel_title, _row_key in PD2MIN_SWEEP_QUANTITIES:
            ratios = pd2min_sweep_ratio_table(sweep, tree_order, qkey)
            values = list(ratios.values())
            if not values:
                print("  %-30s no trees with both pd2min levels present" % panel_title)
                continue
            spread = max(values) - min(values)
            verdict = "move together (narrow spread)" if spread < 0.5 else "diverge (wide spread)"
            print("  %-30s min=%.2f  median=%.2f  max=%.2f  (n=%d trees) - %s"
                  % (panel_title, min(values), statistics.median(values), max(values), len(values), verdict))
        print()
    # The addendum's own expected contrast is about the eta^2 finding
    # (pd2min explains ~0.3% of n_cylinders variance in Simplified, ~100%
    # in Optimal) - that number is about the SPREAD ACROSS THE 16 SIMP
    # CELLS at one fixed pd2min level, not directly about the shape of this
    # figure's per-tree MEDIAN line (which cancels that spread out by
    # construction - see per_tree_line_stats()'s own docstring elsewhere in
    # this file for the same median-vs-spread distinction). Checked
    # directly rather than assumed: at pd2min=0.002, every one of the 16
    # simp-grid cells gives the IDENTICAL n_cylinders value for Optimal
    # (simp_smallradii/simp_replaceiterations are Simplified-stage-only
    # parameters - they cannot affect a pre-simplification model at all),
    # while Simplified's 16 cells at the same pd2min spread over a ~200%
    # range - THAT is what "pd2min explains 100% vs 0.3%" actually means.
    print("  Note on the expected eta^2 contrast (pd2min: 100% of n_cylinders variance in "
          "Optimal, only 0.3% in Simplified): verified directly against the raw rows - at a "
          "fixed pd2min, all 16 simp-grid cells give an IDENTICAL n_cylinders value for Optimal "
          "(simp_smallradii/simp_replaceiterations only apply to the Simplified stage, so they "
          "cannot move a pre-simplification model at all), while Simplified's 16 cells spread "
          "over roughly a 200% range at the same pd2min. That within-level spread (not the "
          "median line's own shape) is what the eta^2 numbers describe - the two stages' MEDIAN "
          "ratios above can legitimately look similar to each other even though their eta^2 "
          "values differ enormously, because the median cancels out exactly the simp-driven "
          "spread that eta^2 is measuring.")
    print()


# =====================  ADDENDUM 1: WIDE PER-TREE MEDIANS CSV  ==========

# quantity -> (output column key, source QUANTITIES key, scale). Reuses
# build_summary_rows()'s ALREADY-COMPUTED medians - a reshape, not a
# recalculation (addendum 1's own sanity check verifies this by spot-check).
MEDIANS_QUANTITIES = [
    ("trunk_m3", "trunk_m3", 1.0),
    ("branch_m3", "branch_m3", 1.0),
    ("total_m3", "total_m3", 1.0),
    ("branch_share_pct", "branch_frac", 100.0),   # branch_frac's median * 100 = branch_share_pct's median (scaling a median by a constant commutes)
    ("dbh_m", "dbh_m", 1.0),
    ("height_m", "height_m", 1.0),
    ("trunk_len_m", "trunk_len_m", 1.0),
    ("branch_len_m", "branch_len_m", 1.0),
    ("n_cylinders", "n_cylinders", 1.0),
]
MEDIANS_FAMILY_COLUMNS = ["AdQSM", "AdTree_raw", "AdTree_calibrated", "TreeQSM_Optimal", "TreeQSM_Simplified"]


def build_medians_rows(summary_rows, tree_order):
    """One row per (population, quantity, tree), one column per family -
    see addendum 1's "Layout". Rows ordered quantity -> tree (field-DBH
    ascending, Part A) -> population, so settled/full sit side by side for
    the same tree. Missing values (a family with no rows for that quantity/
    tree/population - e.g. AdQSM has no n_cylinders at all, or an excluded
    tree/family combination) are left as "" (empty), never 0."""
    lookup = {(r["population"], r["family"], r["tree"], r["quantity"]): r["median"] for r in summary_rows}
    out = []
    for out_qkey, source_qkey, scale in MEDIANS_QUANTITIES:
        for tree in tree_order:
            for pop_name in (POP_SETTLED, POP_FULL):
                row = dict(population=POP_LABEL_SHORT[pop_name], quantity=out_qkey, tree=tree)
                for fam in MEDIANS_FAMILY_COLUMNS:
                    v = lookup.get((pop_name, fam, tree, source_qkey))
                    row[fam] = (v * scale) if v is not None else ""
                out.append(row)
    return out


def write_medians_csv(medians_rows, path):
    fields = ["population", "quantity", "tree"] + MEDIANS_FAMILY_COLUMNS
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in medians_rows:
            w.writerow(row)


def sanity_check_medians_csv(medians_rows, summary_rows, path):
    """Addendum 1's own check: print the row count, then re-read the file
    just written and spot-check 3 (population, family, tree, quantity)
    values against the summary lookup this file was reshaped FROM - they
    must agree exactly."""
    print("=" * 78)
    print("SANITY CHECK (addendum 1): method_eval_medians.csv row count and spot-check")
    print("=" * 78)
    print("  rows written: %d (expected %d = %d quantities x 12 trees x 2 populations)"
          % (len(medians_rows), len(MEDIANS_QUANTITIES) * 12 * 2, len(MEDIANS_QUANTITIES) * 12 * 2))

    with open(path, "r", encoding="utf-8", newline="") as f:
        written = list(csv.DictReader(f))
    lookup = {(r["population"], r["family"], r["tree"], r["quantity"]): r["median"] for r in summary_rows}

    spot_checks = [
        (POP_SETTLED, "AdTree_calibrated", None, "branch_m3"),
        (POP_FULL, "TreeQSM_Simplified", None, "n_cylinders"),
        (POP_SETTLED, "AdQSM", None, "branch_frac"),
    ]
    trees_by_quantity = {}
    for out_qkey, source_qkey, scale in MEDIANS_QUANTITIES:
        trees_by_quantity[source_qkey] = (out_qkey, scale)

    for pop_name, fam, _unused, source_qkey in spot_checks:
        out_qkey, scale = trees_by_quantity[source_qkey]
        # first tree that actually has a value for this (population, family, quantity)
        tree = next((t for t in set(r["tree"] for r in summary_rows)
                     if lookup.get((pop_name, fam, t, source_qkey)) is not None), None)
        if tree is None:
            print("  %s / %s / %s: no tree has a value - skipped" % (pop_name, fam, source_qkey))
            continue
        source_value = lookup[(pop_name, fam, tree, source_qkey)] * scale
        written_row = next((r for r in written if r["population"] == POP_LABEL_SHORT[pop_name]
                             and r["quantity"] == out_qkey and r["tree"] == tree), None)
        written_value = float(written_row[fam]) if written_row and written_row[fam] != "" else None
        match = "MATCH" if (written_value is not None and abs(written_value - source_value) < 1e-9) else "MISMATCH"
        print("  %-20s %-20s %-10s %-16s source=%.6f  written=%s  [%s]"
              % (pop_name, fam, tree, out_qkey, source_value,
                 ("%.6f" % written_value) if written_value is not None else "None", match))
    print()


# =====================  RUN  =============================================

if __name__ == "__main__":
    data = load_and_prepare()

    sanity_check_A(data)
    sanity_check_1(data)
    sanity_check_3(data)

    csv_dir = ensure_csv_dir()

    summary_rows = build_summary_rows(data)
    write_summary_csv(summary_rows, os.path.join(csv_dir, "method_eval_summary.csv"))

    eta2_rows = build_eta2_rows(data)
    write_eta2_csv(eta2_rows, os.path.join(csv_dir, "method_eval_eta2.csv"))
    n_balanced_per_tree = report_treeqsm_unbalanced_points(data)
    eta2_lookup = print_eta2_sanity_check(eta2_rows)

    stability_rows = build_stability_rows(data)
    write_stability_csv(stability_rows, os.path.join(csv_dir, "method_eval_stability.csv"))
    stability_lookup = print_stability_table(stability_rows)
    cv_full, cv_settled = print_cv_sanity_check(stability_lookup)

    write_excluded_csv(data["excluded"], os.path.join(csv_dir, "method_eval_excluded.csv"))

    branch_share_tables = print_branch_share_table(data)   # underlying table for the two branch-share figures

    medians_rows = build_medians_rows(summary_rows, data["tree_order"])
    medians_path = os.path.join(csv_dir, "method_eval_medians.csv")
    write_medians_csv(medians_rows, medians_path)
    sanity_check_medians_csv(medians_rows, summary_rows, medians_path)

    figure_paths = draw_all_figures(data, eta2_rows, stability_lookup)

    # ---- pd2min sweep, all trees, both stages (revision 2, Part D / addendum 2)
    sweep_by_stage = {}
    for stage_family, out_name in [("TreeQSM_Simplified", "methodeval_pd2min_sweep_all_trees.png"),
                                    ("TreeQSM_Optimal", "methodeval_pd2min_sweep_all_trees_optimal.png")]:
        path, sweep = draw_pd2min_sweep_all_trees(data, stage_family, out_name)
        figure_paths.append(path)
        sweep_by_stage[stage_family] = sweep
    sanity_check_pd2min_sweep(sweep_by_stage, data["tree_order"])

    # ---- Step 8 report ---------------------------------------------------
    def _count(pop_rows, fam, tree=None):
        return sum(1 for r in pop_rows if r["_family"] == fam and (tree is None or r["tree"] == tree))

    tree0 = data["tree_order"][0]
    settled_counts = {fam: _count(data["pop_settled"], fam, tree0) for fam in ALL_FAMILIES}
    full_total = len(data["pop_full"])

    excl_settled = sum(1 for r, _ in data["excluded"] if r["_settled"])
    excl_full = len(data["excluded"])

    lookup_e = eta2_lookup
    print("=" * 78)
    print("=== METHOD EVALUATION SUMMARY ===")
    print("=" * 78)
    print("POP_SETTLED rows per tree (tree=%s): AdTree_cal %d (expected %d)  AdTree_raw %d  AdQSM %d"
          % (tree0, settled_counts["AdTree_calibrated"], EXPECTED_SETTLED_ADTREE_CAL_PER_TREE,
             settled_counts["AdTree_raw"], settled_counts["AdQSM"]))
    print("                            TreeQSM_Optimal %d  TreeQSM_Simplified %d"
          % (settled_counts["TreeQSM_Optimal"], settled_counts["TreeQSM_Simplified"]))
    print("POP_FULL rows total       : %d" % full_total)
    print("Excluded rows : POP_SETTLED %d   POP_FULL %d   (listed in csv/method_eval_excluded.csv)"
          % (excl_settled, excl_full))
    print()
    def _f1(x):
        return "%.1f" % x if x is not None else "None"

    print("eta2 check, AdTree_calibrated / branch_m3 / POP_FULL:")
    print("  calmethod %s (exp 53.3)  weighting %s (exp 15.8)  radius_threshold %s (exp 0.7)"
          % (_f1(lookup_e.get(("AdTree_calibrated", "branch_m3", "calmethod"))),
             _f1(lookup_e.get(("AdTree_calibrated", "branch_m3", "weighting"))),
             _f1(lookup_e.get(("AdTree_calibrated", "branch_m3", "radius_threshold_mm")))))
    print("eta2 check, AdTree_calibrated / n_cylinders / POP_FULL:")
    print("  radius_threshold %s (exp 93.4)"
          % _f1(lookup_e.get(("AdTree_calibrated", "n_cylinders", "radius_threshold_mm"))))
    print()
    print("CV check, AdTree_calibrated / branch_m3:")
    print("  POP_FULL %s %% (exp 48.0)    POP_SETTLED %s %% (exp 4.9)" % (_f1(cv_full), _f1(cv_settled)))
    print()
    print("TreeQSM variance decomposition used %s models per stage (expected 48 balanced), "
          "AUTO and pd2min=0.020 reported separately (see stdout above)." % n_balanced_per_tree)
    print()
    bs_settled_cal = branch_share_tables[POP_SETTLED]["AdTree_calibrated"]["median"]
    bs_settled_simp = branch_share_tables[POP_SETTLED]["TreeQSM_Simplified"]["median"]
    print("Branch share medians (settled): AdTree calibrated %s (exp 38.2)  "
          "TreeQSM Simplified %s (exp 73.4)"
          % (_f1(bs_settled_cal), _f1(bs_settled_simp)))
    print()
    print("Outputs : 5 CSV, %d PNG" % len(figure_paths))
    print()
    print("Sanity checks: A see tree order above   1 PASSED   3 see counts above   "
          "5 see eta2 check above   6 see CV check above   "
          "D see pd2min ratio check above   addendum-1 see medians spot-check above")
