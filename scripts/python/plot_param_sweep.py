# -*- coding: utf-8 -*-
# =====================================================================
#  Plot how TreeQSM's reconstruction metrics respond to ONE swept
#  parameter (PatchDiam2Min, "pd2min") across the small number of
#  DISTINCT reconstructions actually run for one tree.
#
#  WHY this is a separate script from plot_box.py, not a mode of it:
#  volume_results.csv holds many repeated exports of the SAME underlying
#  .mat model across simp_smallradii / simp_replaceiterations / the
#  Optimal-vs-Simplified pair - after deduplicating on the actual
#  reconstruction parameters (RECON_KEY below), each distinct
#  reconstruction collapses to exactly ONE row. plot_box.py's box plots
#  describe a DISTRIBUTION of several values per group - with a single
#  value per group after deduplication, a box would degenerate to a
#  flat line with no distribution to show. This script instead plots
#  ONE POINT per reconstruction, connected in sweep order, which is the
#  shape this data actually has.
#
#  Deliberately excluded (per explicit instruction, not an oversight):
#   - no Spearman/Pearson correlation coefficients anywhere - n=4 on the
#     clean SERIES sweep is far too small for a coefficient to mean
#     anything.
#   - no box plots (see above).
#   - no unit conversion / derived columns beyond what's specified here
#     - e.g. dbh is plotted in METRES exactly as load_results() returns
#     it, never converted to cm.
# =====================================================================

import os
import statistics

import matplotlib.pyplot as plt

from compare_volumes import RESULTS_CSV, load_results, to_float
from plot_volumes import ensure_plots_dir, PLOTS_DIR, treeqsm_pd_token, _pd_field_token
from plot_style import (
    FAMILY_GRADIENTS,
    family_shades,
    DEFAULT_LINEWIDTH,
    DEFAULT_MARKERSIZE,
    TITLE_FONTSIZE,
    PANEL_TITLE_FONTSIZE,
    AXIS_LABEL_FONTSIZE,
    LEGEND_FONTSIZE,
    ANNOTATION_FONTSIZE,
)

# to_float/_pd_field_token are imported for parity with the conventions
# named for this script (compare_volumes.py's/plot_volumes.py's own
# import lists) - not called directly here: load_results() already
# returns parsed floats for every field this script reads, and
# treeqsm_pd_token() already calls _pd_field_token() internally for each
# of the 3 PatchDiam values, so there is no separate direct call needed.
_ = to_float, _pd_field_token

# =====================  PARAMETERS  ===================================
SELECT_TREE = "B21_S01"
BRANCH_FILTER = "none"
METHOD_PREFIX = "TreeQSM mine"
MODEL_VARIANT_SUFFIX = "Optimal)"

# RECON_KEY / SWEEP_PARAM / SERIES_BASELINE use load_results()'s OWN
# dict keys, confirmed for real against its actual output (see this
# project's Phase A/C0 diagnostics) - NOT re-derived or guessed here.
#
# IMPORTANT ASYMMETRY - do not "fix" this later by mistake:
# load_results()'s CSV-header -> dict-key mapping is a hardcoded
# per-column dict, not one consistent rule. The 5 reconstruction-key
# columns below happen to have their "_m" unit suffix STRIPPED
# (pd1_m -> pd1, mincylrad_m -> mincylrad, ...), while the 3 pmdist_*
# columns in PANEL_FIELDS further down keep their FULL CSV header name
# unchanged (pmdist_mean -> pmdist_mean - nothing to strip, there's no
# unit suffix on those to begin with). Both are correct; they are simply
# different columns that were named differently by hand.
RECON_KEY = ["pd1", "pd2min", "pd2max", "mincylrad", "mode"]
# raw CSV headers: pd1_m, pd2min_m, pd2max_m, mincylrad_m, mode

SWEEP_PARAM = "pd2min"          # raw CSV header: pd2min_m

SERIES_BASELINE = {"pd1": 0.07, "pd2max": 0.10,
                   "mincylrad": 0.0025, "mode": "manual"}

# PANEL_FIELDS: (load_results() dict key, axis label) for the 8 metric
# columns - dict keys confirmed for real against load_results()'s actual
# output before this file was written (see the project's Phase A/C0
# diagnostics: all 8 present, none missing).
PANEL_FIELDS = [
    ("pmdist_mean", "pmdist_mean [m]"),
    ("pmdist_trunk_mean", "pmdist_trunk_mean [m]"),
    ("pmdist_branch_mean", "pmdist_branch_mean [m]"),
    ("dbh", "DBH [m]"),
    ("trunk", "Trunk volume [m^3]"),
    ("branch", "Branch volume [m^3]"),
    ("total", "Total volume [m^3]"),
    ("n_cylinders", "Number of cylinders"),
]

X_SCALE = "linear"
NCOLS = 4
SHOW_SPREAD_ANNOTATION = True

# LABEL_SERIES_POINTS: within SERIES, every reconstruction-key value
# except the swept one (SWEEP_PARAM) is fixed by SERIES_BASELINE - so a
# SERIES point's run tag is fully determined by its own x position
# already, and the label adds no information while causing overlap at
# closely spaced pd2min values. Default False.
LABEL_SERIES_POINTS = False
# LABEL_OFF_SERIES_POINTS: OFF_SERIES points differ from each other (and
# from SERIES) in mincylrad/mode, neither of which appears on any axis -
# so their run tag carries real information the plot has no other way
# to show. Default True.
LABEL_OFF_SERIES_POINTS = True

# SPREAD_IN_TITLE: True (default) appends the spread to each panel's own
# title (e.g. "Branch volume [m^3]  (spread 84.3 %)") instead of a
# separate floating annotation box. Reason: a fixed-corner annotation
# box collides with data points whose position varies per panel; the
# title is always free regardless of where the data sits. The separate-
# annotation code (used when this is False) is kept intact below, not
# deleted. SHOW_SPREAD_ANNOTATION (above) remains the master on/off
# switch for BOTH forms - when it's False, neither the title suffix nor
# the separate annotation box is drawn.
SPREAD_IN_TITLE = True

# OFF_SERIES_LABEL_ALTERNATE: two OFF_SERIES points can share the same
# SWEEP_PARAM x value (e.g. pd2min=0.020 for both mincylrad=0.00125 and
# mincylrad=0.05595) - in whichever panels their y values also land
# close together, their labels collide. When True (default), points are
# grouped by their exact x value and, within any group of 2+, the label
# offset alternates above/below (+-OFF_SERIES_LABEL_DY_PT) instead of
# every label using the same fixed offset - see _off_series_offsets()
# below. A group of 1 (nothing to disambiguate) keeps today's placement
# unchanged. When False, every OFF_SERIES label keeps today's single
# fixed offset - that code path is kept intact, not deleted.
OFF_SERIES_LABEL_ALTERNATE = True
OFF_SERIES_LABEL_DY_PT = 6      # vertical offset in points, for the alternating case

# OFF_SERIES_LABEL_MODE: OFF_SERIES points sit at pd2min 0.020/0.025,
# close to the right edge of the x range - no INLINE offset direction
# fixes label overflow there (right spills out of the axes, left runs
# over the data). "legend" (default) avoids the problem entirely by
# identifying each OFF_SERIES point via a distinct marker+colour in the
# shared figure legend instead of text inside the panels - see
# _off_series_legend_style() below. "inline" keeps the EXACT previous
# behaviour (single shared marker/colour, text via LABEL_OFF_SERIES_
# POINTS/OFF_SERIES_LABEL_ALTERNATE/OFF_SERIES_LABEL_DY_PT) - that code
# path is kept intact, not deleted. "none" draws neither: OFF_SERIES
# points still appear on the plot (never remove data points), but with
# one shared marker/colour and one shared, generic legend entry - no
# per-point text and no per-point legend entry either.
OFF_SERIES_LABEL_MODE = "legend"     # "legend" | "inline" | "none"

# Marker shape per OFF_SERIES point in "legend" mode, cycled in the SAME
# deterministic order as _off_series_offsets() already uses (sort by
# mincylrad, then mode) - see _off_series_legend_style(). If there are
# more OFF_SERIES points than entries here, the list is cycled with one
# console warning rather than crashing.
OFF_SERIES_MARKERS = ["X", "^", "s"]
# =====================================================================


def _run_tag(row):
    """Short, unambiguous label for one reconstruction, built from
    RECON_KEY's own values - NEVER from the CSV "method" string, which
    (after deduplication) still carries the simp_smallradii/
    simp_replaceiterations tokens of whichever export happened to be
    FIRST in the file for that reconstruction - an arbitrary and
    misleading choice once collapsed to one row per reconstruction.

    Reuses treeqsm_pd_token() (plot_volumes.py, which itself calls
    _pd_field_token()) for the PatchDiam part, and the exact same
    mcr_tag/mode-shortening convention shorten_method_label()
    (plot_volumes.py) already uses for MinCylRad/mode - not new
    formatting logic. NOTE: treeqsm_pd_token()'s own "clean" PatchDiam
    format is NOT zero-padded (e.g. pd1=0.07 -> "7", not "07") - that is
    its existing, established behaviour elsewhere in this project, kept
    as-is here rather than hand-building a differently-padded format."""
    pd_token = treeqsm_pd_token(row["pd1"], row["pd2min"], row["pd2max"])
    mincylrad = row["mincylrad"]
    if mincylrad is None or abs(mincylrad - 0.0025) < 1e-9:
        mcr_tag = ""
    else:
        mcr_tag = "_mcr%03d" % round(mincylrad * 10000)
    mode_short = {"manual": "man", "auto": "aut"}.get(row["mode"], row["mode"])
    return "%s%s_%s" % (pd_token, mcr_tag, mode_short)


def _off_series_offsets(off_series_rows, sweep_param, alternate, dy_pt):
    """Compute a per-row (dx, dy) annotate() offset, IN POINTS (for
    textcoords="offset points" - display units, so the offset behaves
    the same regardless of each panel's own y data range/scale), for
    every OFF_SERIES row - so two rows sharing the same SWEEP_PARAM x
    value get their labels pushed apart vertically instead of landing on
    top of each other in whichever panels they also land close in y.

    alternate=False: every row keeps today's single fixed offset
    (4, -10) - unchanged behaviour, kept intact for that case.

    alternate=True: rows are grouped by their exact SWEEP_PARAM value.
    A group of 1 (nothing to disambiguate) keeps the SAME (4, -10)
    placement. A group of 2+ is sorted DETERMINISTICALLY by
    (mincylrad, mode) - so the figure is reproducible run to run - and
    offsets alternate (4, +dy_pt), (4, -dy_pt), (4, +dy_pt), ... in that
    sorted order.

    Returns {id(row): (dx, dy)}, keyed by object identity - the row
    dicts from load_results() aren't meant to be hashed by value here,
    and every row passed in is already a distinct dict object."""
    offsets = {}
    if not alternate:
        for r in off_series_rows:
            offsets[id(r)] = (4, -10)
        return offsets

    groups = {}
    for r in off_series_rows:
        groups.setdefault(r[sweep_param], []).append(r)

    for group_rows in groups.values():
        if len(group_rows) == 1:
            offsets[id(group_rows[0])] = (4, -10)
            continue
        ordered = sorted(group_rows, key=lambda r: (
            r["mincylrad"] if r["mincylrad"] is not None else float("-inf"),
            r["mode"] or "",
        ))
        for i, r in enumerate(ordered):
            dy = dy_pt if i % 2 == 0 else -dy_pt
            offsets[id(r)] = (4, dy)
    return offsets


def _off_series_legend_style(off_series_rows, markers):
    """For OFF_SERIES_LABEL_MODE == "legend": assign each OFF_SERIES row
    its own marker (cycled from `markers`) and a colour-index (0..n-1,
    for indexing into a same-length family_shades(...) list built by the
    caller), in the SAME deterministic (mincylrad, mode) order
    _off_series_offsets() already uses within a group - so the
    marker/colour assignment is reproducible run to run, independent of
    load_results()'s own row order.

    If there are more rows than markers, `markers` is cycled (modulo)
    and a single console warning is printed - never raises.

    Returns {id(row): (marker, color_idx)}."""
    ordered = sorted(off_series_rows, key=lambda r: (
        r["mincylrad"] if r["mincylrad"] is not None else float("-inf"),
        r["mode"] or "",
    ))
    if markers and len(ordered) > len(markers):
        print("WARNING: %d OFF_SERIES point(s) but only %d marker(s) in "
              "OFF_SERIES_MARKERS - cycling the marker list." % (len(ordered), len(markers)))
    style = {}
    for i, r in enumerate(ordered):
        marker = markers[i % len(markers)] if markers else "X"
        style[id(r)] = (marker, i)
    return style


def _spread_pct(values):
    """(max-min)/median*100 over `values` - the ONLY summary statistic
    this script ever computes; no correlation coefficient anywhere (n=4
    on the SERIES sweep is far too small for one to mean anything)."""
    vmin, vmax = min(values), max(values)
    median = statistics.median(values)
    return (vmax - vmin) / median * 100 if median else float("nan")


def run():
    rows = load_results(RESULTS_CSV)

    # ---- pipeline step 2: filter (tree / branch_filter / method prefix+suffix) ----
    step_tree = [r for r in rows if r["tree"] == SELECT_TREE]
    step_bf = [r for r in step_tree if r["branch_filter"] == BRANCH_FILTER]
    step_prefix = [r for r in step_bf if r["method"].startswith(METHOD_PREFIX)]
    filtered = [r for r in step_prefix if r["method"].endswith(MODEL_VARIANT_SUFFIX)]
    print("Filter pipeline: tree=%d -> branch_filter=%d -> prefix=%d -> suffix=%d"
          % (len(step_tree), len(step_bf), len(step_prefix), len(filtered)))
    if len(filtered) != 40:
        print("WARNING: expected 40 rows after filtering, got %d - continuing anyway."
              % len(filtered))

    # ---- pipeline step 3: deduplicate on RECON_KEY, keep first ----------------
    seen = set()
    dedup = []
    for r in filtered:
        key = tuple(r[k] for k in RECON_KEY)
        if key not in seen:
            seen.add(key)
            dedup.append(r)
    print("After dedup on RECON_KEY: %d rows" % len(dedup))
    if len(dedup) != 7:
        print("WARNING: expected 7 rows after dedup, got %d - continuing anyway." % len(dedup))

    # ---- pipeline step 4: split SERIES / OFF_SERIES ----------------------------
    def matches_baseline(r):
        return all(r[k] == v for k, v in SERIES_BASELINE.items())

    SERIES = sorted([r for r in dedup if matches_baseline(r)], key=lambda r: r[SWEEP_PARAM])
    OFF_SERIES = [r for r in dedup if not matches_baseline(r)]
    print("SERIES: %d rows, OFF_SERIES: %d rows" % (len(SERIES), len(OFF_SERIES)))
    if len(SERIES) != 4 or len(OFF_SERIES) != 3:
        print("WARNING: expected SERIES=4 / OFF_SERIES=3, got SERIES=%d / OFF_SERIES=%d - "
              "continuing anyway." % (len(SERIES), len(OFF_SERIES)))

    if not SERIES:
        raise SystemExit("plot_param_sweep.py: no SERIES rows matched SERIES_BASELINE - nothing to plot.")

    # ---- build the small-multiples figure ---------------------------------------
    n_panels = len(PANEL_FIELDS)
    nrows = -(-n_panels // NCOLS)  # ceil division, no numpy needed for this
    fig, axes = plt.subplots(nrows, NCOLS, figsize=(4.5 * NCOLS, 4 * nrows), sharex=True)
    axes_flat = list(axes.flatten()) if n_panels > 1 else [axes]

    series_x = [r[SWEEP_PARAM] for r in SERIES]
    off_x = [r[SWEEP_PARAM] for r in OFF_SERIES]
    # Colors from the shared palette (plot_style.py's FAMILY_GRADIENTS) -
    # SERIES uses TreeQSM's own family's deepest stop (this is TreeQSM
    # data); OFF_SERIES uses the Reference family's flat highlight stop,
    # so it reads as visually distinct/"off the main sweep" rather than
    # just another TreeQSM shade.
    series_color = FAMILY_GRADIENTS["TreeQSM"][-1]
    off_color = FAMILY_GRADIENTS["Reference"][1]

    # Computed ONCE (not per panel) - the offset only depends on which
    # OFF_SERIES rows share an x value, not on any panel's own y data.
    off_series_offsets = _off_series_offsets(
        OFF_SERIES, SWEEP_PARAM, OFF_SERIES_LABEL_ALTERNATE, OFF_SERIES_LABEL_DY_PT)

    # OFF_SERIES_LABEL_MODE == "legend": one marker+colour per point,
    # computed ONCE (not per panel) - colours sampled from the shared
    # palette (plot_style.py's family_shades()), not hardcoded, one per
    # OFF_SERIES row.
    off_series_style = {}
    off_series_colors = []
    if OFF_SERIES_LABEL_MODE == "legend" and OFF_SERIES:
        off_series_style = _off_series_legend_style(OFF_SERIES, OFF_SERIES_MARKERS)
        off_series_colors = family_shades("Reference", len(OFF_SERIES))

    for panel_idx, (field_key, axis_label) in enumerate(PANEL_FIELDS):
        ax = axes_flat[panel_idx]
        series_y = [r[field_key] for r in SERIES]

        ax.plot(series_x, series_y, marker="o", markersize=DEFAULT_MARKERSIZE,
                linewidth=DEFAULT_LINEWIDTH, color=series_color,
                label="SERIES (%s sweep)" % SWEEP_PARAM, zorder=2)
        if LABEL_SERIES_POINTS:
            # Off by default - see LABEL_SERIES_POINTS' own comment in the
            # PARAMETERS block above for why (redundant with x position,
            # causes overlap at closely spaced pd2min values). Code kept
            # intact, not deleted, for whoever flips the flag back on.
            for r, x, y in zip(SERIES, series_x, series_y):
                ax.annotate(_run_tag(r), (x, y), fontsize=ANNOTATION_FONTSIZE,
                            textcoords="offset points", xytext=(4, 4))

        if OFF_SERIES:
            off_y = [r[field_key] for r in OFF_SERIES]
            # OFF_SERIES markers are NEVER joined by a line (scatter,
            # not plot, in every mode below) - they vary a DIFFERENT
            # parameter (mincylrad, or auto mode) at this SAME pd2min
            # value, so drawing a line through them would misleadingly
            # suggest they belong to the pd2min sweep.
            if OFF_SERIES_LABEL_MODE == "legend":
                # One scatter call PER POINT, each with its own marker/
                # colour/label - no inline text at all. Avoids the
                # inline-label overflow problem entirely (OFF_SERIES
                # points sit close to the right edge of the pd2min
                # range, so no inline offset direction stays inside the
                # axes) by moving identification into the shared figure
                # legend instead.
                for r in OFF_SERIES:
                    marker, color_idx = off_series_style[id(r)]
                    ax.scatter([r[SWEEP_PARAM]], [r[field_key]], marker=marker,
                               s=DEFAULT_MARKERSIZE ** 2 * 3, color=off_series_colors[color_idx],
                               label=_run_tag(r), zorder=3)
            elif OFF_SERIES_LABEL_MODE == "inline":
                # EXACT previous behaviour, kept intact - single shared
                # marker/colour, one shared legend entry, per-point text
                # via LABEL_OFF_SERIES_POINTS/OFF_SERIES_LABEL_ALTERNATE/
                # OFF_SERIES_LABEL_DY_PT.
                ax.scatter(off_x, off_y, marker="X", s=DEFAULT_MARKERSIZE ** 2 * 3,
                           color=off_color, label="OFF_SERIES (other parameter varies)",
                           zorder=3)
                if LABEL_OFF_SERIES_POINTS:
                    # On by default - these points differ in mincylrad/
                    # mode, neither shown on any axis, so the label is
                    # the only place that information appears (see
                    # LABEL_OFF_SERIES_POINTS' own comment above).
                    for r, x, y in zip(OFF_SERIES, off_x, off_y):
                        dx, dy = off_series_offsets[id(r)]
                        ax.annotate(_run_tag(r), (x, y), fontsize=ANNOTATION_FONTSIZE,
                                    textcoords="offset points", xytext=(dx, dy))
            else:   # "none" - points still drawn (never remove data), no
                    # per-point text, no per-point legend entry either -
                    # one shared marker/colour/legend entry, same style
                    # as the "inline" scatter above but with nothing else.
                ax.scatter(off_x, off_y, marker="X", s=DEFAULT_MARKERSIZE ** 2 * 3,
                           color=off_color, label="OFF_SERIES (other parameter varies)",
                           zorder=3)

        ax.set_xscale(X_SCALE)
        ax.set_xlabel("%s [m]" % SWEEP_PARAM, fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel(axis_label, fontsize=AXIS_LABEL_FONTSIZE)

        # Spread text: same underlying number (_spread_pct(series_y)),
        # just given a different HOME depending on SPREAD_IN_TITLE - see
        # that constant's own comment in the PARAMETERS block above for
        # why the title is the default. Both branches are gated behind
        # SHOW_SPREAD_ANNOTATION, the master on/off switch for either form.
        if SHOW_SPREAD_ANNOTATION and SPREAD_IN_TITLE:
            spread_pct = _spread_pct(series_y)
            ax.set_title("%s  (spread %.1f %%)" % (axis_label, spread_pct), fontsize=PANEL_TITLE_FONTSIZE)
        else:
            ax.set_title(axis_label, fontsize=PANEL_TITLE_FONTSIZE)

        if SHOW_SPREAD_ANNOTATION and not SPREAD_IN_TITLE:
            # Separate-annotation form, kept intact (not deleted) for
            # SPREAD_IN_TITLE = False. matplotlib autoscales each panel's
            # y-axis independently, so a genuinely FLAT quantity (e.g.
            # dbh, ~1.6% spread) would otherwise fill the whole panel
            # height and look like a dramatic trend, indistinguishable
            # from a panel with a real large spread (e.g. branch, ~84%)
            # unless the actual magnitude is stated numerically - that's
            # what this annotation is for, computed over SERIES rows only.
            spread_pct = _spread_pct(series_y)
            ax.text(0.02, 0.98, "spread = %.1f %%" % spread_pct, transform=ax.transAxes,
                    ha="left", va="top", fontsize=ANNOTATION_FONTSIZE,
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="lightgray"))

    # Any leftover panel slots (NCOLS doesn't evenly divide len(PANEL_FIELDS))
    # are hidden rather than left showing empty axes.
    for extra_idx in range(n_panels, len(axes_flat)):
        axes_flat[extra_idx].axis("off")

    # ONE legend for the whole figure (every panel repeats the same set
    # of handles - 1 SERIES + 1 shared OFF_SERIES entry in "inline"/
    # "none" mode, or 1 SERIES + one entry PER OFF_SERIES point in
    # "legend" mode) - built from panel 0's own handles/labels, same
    # pattern plot_volumes.py's plot_tree_overview() uses for its shared
    # legend. ncol=2 (unchanged) already wraps to a second row on its own
    # once there are more than 2 entries (4, in "legend" mode with today's
    # data) - no font shrinking involved.
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=LEGEND_FONTSIZE,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("%s: %s sweep (TreeQSM, branch_filter=%r)" % (SELECT_TREE, SWEEP_PARAM, BRANCH_FILTER),
                 fontsize=TITLE_FONTSIZE)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    ensure_plots_dir()   # creates PLOTS_DIR if it doesn't exist yet
    out_path = os.path.join(PLOTS_DIR, "param_sweep_%s_%s.png" % (SELECT_TREE, SWEEP_PARAM))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)

    # ---- spread table, printed for review (not interpreted here) --------------
    print()
    print("Spread (%%) per panel field, over the %d SERIES rows only:" % len(SERIES))
    for field_key, _axis_label in PANEL_FIELDS:
        series_y = [r[field_key] for r in SERIES]
        print("  %-20s spread=%.2f %%" % (field_key, _spread_pct(series_y)))


if __name__ == "__main__":
    run()
