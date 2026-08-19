# -*- coding: utf-8 -*-
# =====================================================================
#  Draw and save charts that summarize volume_results.csv (the shared
#  master results table produced by ply_to_geom.py, qsm_volume_mean.py,
#  reference_volume.py, ... and printed in detail by compare_volumes.py).
# ---------------------------------------------------------------------
#  This script makes several PNG charts every time you run it:
#
#   a) total_volume_by_tree.png
#        Bar chart of total_m3 per method, with one GROUP of bars per tree.
#        The reference method's bar is drawn in a different colour and
#        labelled "(reference)" in the legend, so it's easy to spot.
#        Always drawn, no matter how many trees are in the CSV.
#
#   b) error_boxplot.png
#        Box plot of the percentage error of each method's total_m3 vs. the
#        reference, one box per method, built from the SAME trees used for
#        that method's box (percent error computed the same way as
#        pct_diff() in compare_volumes.py - imported from there, not
#        re-implemented). Needs at least 2 trees to say anything meaningful
#        (comparing across trees is the whole point) - with fewer, the RUN
#        section below skips it and prints why, instead of drawing something
#        misleading.
#
#   c) error_metrics_bar.png
#        Bar chart of Bias / MAE / RMSE per method, using the exact same
#        calculation as compare_volumes.py's error_metrics() table
#        (via compute_error_metrics(), imported from compare_volumes.py -
#        so the numbers here can never drift out of sync with the printed
#        table). Same "needs >=2 trees" rule as (b), for the same reason.
#
#   d) tree_overview_<tree>.png (one file PER tree in the CSV)
#        A single figure with a 2x2 grid of bar charts for ONE tree: total
#        volume, DBH, height and taper - one bar per method in each subplot.
#        Unlike (b)/(c), this makes sense even with just ONE tree in the
#        CSV (it doesn't compare across trees, only across methods for the
#        SAME tree), so it is always drawn for every tree found in the CSV.
#
#  All PNGs are written into a "plots" subfolder next to this script
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
import matplotlib.patches as mpatches   # used to build the shared legend in plot_tree_overview()

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
# d) One tree's overview: 2x2 grid of bar charts (total volume, DBH,
#    height, taper), one bar per method, for a SINGLE tree. Makes sense
#    even with only 1 tree in the CSV (unlike the boxplot/RMSE charts,
#    which compare ACROSS trees).
# ----------------------------------------------------------------------
def plot_tree_overview(rows, tree):
    tree_rows = [r for r in rows if r["tree"] == tree]
    if not tree_rows:
        print("No rows for tree '%s' - skipping tree overview." % tree)
        return

    # Methods in the order they first appear for THIS tree (dict.fromkeys()
    # trick again, see plot_total_volume_by_tree). Keeping a stable order
    # means each method gets the same colour every time you re-run this.
    methods = list(dict.fromkeys(r["method"] for r in tree_rows))
    # One row dict per method, for quick lookups below (assumes at most one
    # row per (tree, method) pair, which is how upsert_result() keeps the CSV).
    row_of = {r["method"]: r for r in tree_rows}
    ref_row = row_of.get(REFERENCE_METHOD)   # None if this tree has no reference row

    # Pick ONE fixed colour per method and reuse it in all 4 subplots AND
    # the shared legend below, so e.g. "AdQSM" is always the same colour.
    # The reference method always gets the same highlight colour used in
    # plot_total_volume_by_tree, for visual consistency between the charts.
    # matplotlib's default cycle includes an orange that is the same colour
    # as "tab:orange" (used for the reference highlight below) - compare via
    # to_rgba() rather than the raw value, since depending on matplotlib
    # version/theme the cycle may hold hex strings OR (r, g, b) tuples.
    import matplotlib.colors as mcolors
    default_colors = [c for c in plt.rcParams["axes.prop_cycle"].by_key()["color"]
                       if mcolors.to_rgba(c) != mcolors.to_rgba("tab:orange")]
    color_of = {}
    next_color_i = 0
    for m in methods:
        if m == REFERENCE_METHOD:
            color_of[m] = "tab:orange"
        else:
            color_of[m] = default_colors[next_color_i % len(default_colors)]
            next_color_i += 1

    # Which field (as returned by load_results()) goes in which subplot,
    # and what to title that subplot.
    fields = [
        ("total",  "Total volume [m^3]"),
        ("dbh",    "DBH [m]"),
        ("height", "Height [m]"),
        ("taper",  "Taper [cm/m]"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    for ax, (field_key, subplot_title) in zip(axes.flat, fields):
        # Skip methods with no value (None) for THIS field entirely, instead
        # of trying to draw a bar for them (which would crash matplotlib).
        present_methods = [m for m in methods if row_of[m][field_key] is not None]
        values = [row_of[m][field_key] for m in present_methods]
        colors = [color_of[m] for m in present_methods]

        x_positions = list(range(len(present_methods)))
        ax.bar(x_positions, values, color=colors)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(present_methods, rotation=30, ha="right", fontsize=8)
        ax.set_title(subplot_title)

        if not present_methods:
            # Nothing to plot for this field at all (e.g. no method has a
            # taper value yet) - say so instead of leaving a blank mystery panel.
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            continue

        # Percent-difference label above every non-reference bar, using the
        # same pct_diff() calculation compare_volumes.py uses for its table.
        ref_value = ref_row[field_key] if ref_row is not None else None
        for xi, m in zip(x_positions, present_methods):
            if m == REFERENCE_METHOD:
                continue
            d_pct = pct_diff(row_of[m][field_key], ref_value)
            if d_pct is None:   # no reference value to compare against - leave blank
                continue
            ax.annotate("%+.0f%%" % d_pct,
                        xy=(xi, row_of[m][field_key]),
                        xytext=(0, 3), textcoords="offset points",   # 3 points above the bar top
                        ha="center", va="bottom", fontsize=7)

    fig.suptitle("Tree overview: %s" % tree, fontsize=14)

    # ONE legend for the whole figure (not one per subplot, which would just
    # repeat the same method names 4 times) - built from coloured squares
    # ("patches") rather than real bar handles, since not every method has a
    # bar in every subplot.
    legend_handles = [
        mpatches.Patch(color=color_of[m],
                       label=m + (" (reference)" if m == REFERENCE_METHOD else ""))
        for m in methods
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=min(3, len(methods)), fontsize=8, bbox_to_anchor=(0.5, -0.02))

    # rect leaves room at the top for suptitle and at the bottom for the legend.
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    save_and_report(fig, "tree_overview_%s.png" % tree)


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
    ax.boxplot(data, tick_labels=labels)
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
    all_trees = sorted({r["tree"] for r in rows})

    # Always makes sense, regardless of how many trees are in the CSV.
    plot_total_volume_by_tree(rows)

    # One overview PNG per tree currently in the CSV - loops over whatever
    # is actually there, so new trees get their own chart automatically the
    # next time you run this, with nothing to edit here.
    for tree in all_trees:
        plot_tree_overview(rows, tree)

    # The boxplot and RMSE/Bias/MAE charts compare methods ACROSS trees, so
    # they need at least 2 trees to say anything real; with only 1, skip
    # them (instead of drawing a single degenerate point/bar) and say why.
    if len(all_trees) >= 2:
        plot_error_boxplot(rows)
        plot_error_metrics_bar(rows)
    else:
        print("Only %d tree in %s - skipping error_boxplot.png and "
              "error_metrics_bar.png (both compare methods ACROSS trees, so "
              "they need at least 2 trees to be meaningful). Add more trees "
              "to the CSV and re-run to get them." % (len(all_trees), RESULTS_CSV))
