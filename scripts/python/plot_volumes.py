# -*- coding: utf-8 -*-
# =====================================================================
#  Draw and save charts that summarize volume_results.csv (the shared
#  master results table produced by ply_to_geom.py, qsm_volume_mean.py,
#  reference_volume.py, ... and printed in detail by compare_volumes.py).
# ---------------------------------------------------------------------
#  This script makes THREE PNG charts every time you run it:
#
#   a) total_volume_by_tree.png
#        Bar chart of total_m3 per method, with one GROUP of bars per tree.
#        The reference method's bar is drawn in a different colour and
#        labelled "(reference)" in the legend, so it's easy to spot.
#
#   b) error_boxplot.png
#        Box plot of the percentage error of each method's total_m3 vs. the
#        reference, one box per method, built from the SAME trees used for
#        that method's box (percent error computed the same way as
#        pct_diff() in compare_volumes.py - imported from there, not
#        re-implemented). With fewer than 3 trees a box plot isn't very
#        informative (each "box" is really just 1-2 points) - the script
#        still draws it, but prints a warning to the console explaining this.
#
#   c) error_metrics_bar.png
#        Bar chart of Bias / MAE / RMSE per method, using the exact same
#        calculation as compare_volumes.py's error_metrics() table
#        (via compute_error_metrics(), imported from compare_volumes.py -
#        so the numbers here can never drift out of sync with the printed
#        table).
#
#  All three PNGs are written into a "plots" subfolder next to this script
#  (created automatically if it doesn't exist yet). The path of each saved
#  file is printed to the console.
#
#  Dependencies: matplotlib   (install: pip install matplotlib)
#  This script also IMPORTS a few things from compare_volumes.py, which
#  must live in the same folder:
#     RESULTS_CSV, REFERENCE_METHOD, load_results, pct_diff, compute_error_metrics
# =====================================================================

import os
import matplotlib.pyplot as plt

from compare_volumes import (
    RESULTS_CSV,
    REFERENCE_METHOD,
    load_results,
    pct_diff,
    compute_error_metrics,
)

# =====================  PARAMETERS  ==================================
# Folder (relative to this script's working directory) where the PNG
# charts are saved. Created automatically if it doesn't exist.
PLOTS_DIR = "plots"
# =====================================================================


def ensure_plots_dir():
    """Create PLOTS_DIR if it doesn't exist yet, and return its path."""
    if not os.path.isdir(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)   # makedirs (not mkdir) also creates parent folders if needed
    return PLOTS_DIR


def save_and_report(fig, filename):
    """Save one matplotlib figure into PLOTS_DIR, close it (frees memory),
    and print the full path so you know where to look for it."""
    out_path = os.path.join(ensure_plots_dir(), filename)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


# ----------------------------------------------------------------------
# a) Bar chart: total_m3 per method, grouped by tree.
# ----------------------------------------------------------------------
def plot_total_volume_by_tree(rows):
    trees = sorted({r["tree"] for r in rows})
    # dict.fromkeys() keeps the methods in the order they first appear in
    # the CSV (a plain set() would print them in a random order every run).
    methods = list(dict.fromkeys(r["method"] for r in rows))

    # Quick lookup: (tree, method) -> total_m3 (may be missing -> None).
    total_of = {(r["tree"], r["method"]): r["total"] for r in rows}

    n_methods = len(methods)
    # Each tree gets a group of bars 0.8 units wide, split evenly between
    # the methods, with a small gap (the leftover 0.2) to the next tree's group.
    bar_width = 0.8 / max(n_methods, 1)

    fig, ax = plt.subplots(figsize=(max(6, 1.6 * len(trees) * n_methods), 6))

    for i, method in enumerate(methods):
        heights = [total_of.get((t, method)) or 0.0 for t in trees]  # missing -> 0 (drawn as no bar)
        # x position of this method's bar within each tree's group
        x = [tree_idx + i * bar_width for tree_idx in range(len(trees))]
        is_ref = (method == REFERENCE_METHOD)
        ax.bar(x, heights, width=bar_width,
               color="tab:orange" if is_ref else None,
               label=method + (" (reference)" if is_ref else ""))

    # Put the x-axis tick for each tree in the middle of its group of bars.
    group_centers = [tree_idx + bar_width * (n_methods - 1) / 2.0 for tree_idx in range(len(trees))]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(trees)
    ax.set_xlabel("Tree")
    ax.set_ylabel("Total volume [m^3]")
    ax.set_title("Total volume by method, per tree")
    # Legend outside the plot area (to the right) so it doesn't cover bars.
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    save_and_report(fig, "total_volume_by_tree.png")


# ----------------------------------------------------------------------
# b) Box plot: percentage error vs. reference, per method, across trees.
# ----------------------------------------------------------------------
def plot_error_boxplot(rows):
    trees = sorted({r["tree"] for r in rows})
    methods = [m for m in dict.fromkeys(r["method"] for r in rows) if m != REFERENCE_METHOD]
    total_of = {(r["tree"], r["method"]): r["total"] for r in rows}

    if len(trees) < 3:
        print("NOTE: only %d tree(s) currently in %s - a box plot only becomes "
              "meaningful with 3+ trees (with fewer, each 'box' is really just "
              "1-2 points). Drawing it anyway so you can see the layout."
              % (len(trees), RESULTS_CSV))

    data = []     # one list of % errors per method (only methods with >=1 value)
    labels = []
    for m in methods:
        errors_pct = []
        for t in trees:
            d = pct_diff(total_of.get((t, m)), total_of.get((t, REFERENCE_METHOD)))
            if d is not None:
                errors_pct.append(d)
        if errors_pct:
            data.append(errors_pct)
            labels.append(m)

    if not data:
        print("No method has both a total_m3 value and a reference value - skipping box plot.")
        return

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(labels)), 6))
    ax.boxplot(data, labels=labels)
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)  # 0% error = perfect match
    ax.set_ylabel("Error vs. reference [%]")
    ax.set_title("Total-volume error distribution by method (across trees)")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    save_and_report(fig, "error_boxplot.png")


# ----------------------------------------------------------------------
# c) Bar chart: Bias / MAE / RMSE per method (same numbers as
#    compare_volumes.py's printed "Error metrics" table).
# ----------------------------------------------------------------------
def plot_error_metrics_bar(rows):
    metrics = compute_error_metrics(rows)   # reuses compare_volumes.py's own calculation
    if not metrics:
        print("No method could be compared to the reference - skipping error-metrics chart.")
        return

    methods = [m["method"] for m in metrics]
    bias = [m["bias"] for m in metrics]
    mae = [m["mae"] for m in metrics]
    rmse = [m["rmse"] for m in metrics]

    x = list(range(len(methods)))
    width = 0.25   # 3 bars per method (bias, mae, rmse), each this wide

    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(methods)), 6))
    ax.bar([xi - width for xi in x], bias, width=width, label="Bias")
    ax.bar(x,                        mae,  width=width, label="MAE")
    ax.bar([xi + width for xi in x], rmse, width=width, label="RMSE")
    ax.axhline(0.0, color="gray", linestyle="-", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=30, ha="right")
    ax.set_ylabel("Volume error [m^3]")
    ax.set_title("Error metrics vs. reference '%s'" % REFERENCE_METHOD)
    ax.legend()
    fig.tight_layout()
    save_and_report(fig, "error_metrics_bar.png")


# =========================  RUN  =====================================
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit(
            "'%s' not found - run compare_volumes.py (or one of the volume "
            "scripts, e.g. ply_to_geom.py) first so it gets created." % RESULTS_CSV
        )

    rows = load_results(RESULTS_CSV)

    plot_total_volume_by_tree(rows)
    plot_error_boxplot(rows)
    plot_error_metrics_bar(rows)
