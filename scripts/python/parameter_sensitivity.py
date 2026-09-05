# -*- coding: utf-8 -*-
# =====================================================================
#  Parameter sensitivity analysis: how do AdTree's (adqsm_variant/
#  radius_threshold_mm/seg_*) or TreeQSM's (mode/pd1/pd2min/pd2max/simp_*)
#  input parameters affect the reconstructed outputs, and which combination
#  lands closest to the reference? Three purpose-built charts, replacing
#  ad-hoc box-plot squinting:
#    1) facet grid    - one panel per FACET_PARAM value, METRIC vs.
#                        X_PARAM, lines coloured by LINE_PARAM.
#    2) ranked deviation dot plot - every individual row, sorted by %
#                        deviation from the resolved reference - doubles
#                        as a "which combo is closest" leaderboard.
#    3) correlation heatmap - Spearman correlation between each structured
#                        parameter and each output metric. Fast numeric
#                        screening only - this UNDERSTATES non-monotonic
#                        categorical effects (e.g. AdQSM variant identity
#                        isn't a smooth numeric axis, it's a label) - the
#                        facet grid and ranked plot are the primary tools,
#                        this is a supplement, not a replacement.
#  A 4th option (parallel coordinates) was tried and rejected - too few
#  discrete levels per axis (7x4x3) produced an unreadable "spaghetti" plot
#  with no extra insight over the other three - not implemented here.
# ---------------------------------------------------------------------
#  FAMILY ("adtree"/"treeqsm") and BRANCH_FILTER ("10cm"/"none") are NEVER
#  mixed within one chart - branch_filter alone drives which reference row
#  get_reference() resolves to (the destructive reference for "10cm", the
#  dynamically-resolved smallest-numbered AdQSM variant for "none"),
#  independent of FAMILY: an AdTree analysis and a TreeQSM analysis on the
#  same branch_filter compare against the exact same reference row.
#
#  Dependencies: matplotlib, numpy (install: pip install matplotlib numpy)
# =====================================================================

import math
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines

from compare_volumes import RESULTS_CSV, REFERENCE_METHOD, load_results, resolve_reference_method_none
from plot_volumes import FAMILY_GRADIENTS, classify_family, ensure_plots_dir, shorten_method_label, TREE_MARKERS
from plot_box import parse_treeqsm_method, TREEQSM_REF_LINE_COLOR

# =====================  PARAMETERS  ===================================
SELECT_TREE = "IND01_054"
# NOT read anywhere in the current RUN section - REFERENCE_TREES drives
# every chart now, including single-tree runs (set REFERENCE_TREES to a
# one-element list for that). Kept only in case a future direct/REPL call
# to one of this file's per-tree functions (plot_facet_grid(),
# plot_ranked_deviation(), plot_correlation_heatmap()) wants a default
# `tree=` value to reference - otherwise vestigial.
REFERENCE_TREES = ["IND01_054", "IND03_088", "IND07_083"]
# The 3 trees with a real destructive reference - used by the
# multi-tree ranked-deviation and correlation-heatmap functions below,
# and looped over for per-tree facet grids in the RUN section.

BRANCH_FILTERS_TO_RUN = ["10cm", "none"]
# Both branch_filter modes run automatically in one script execution.
# Set to a single-element list (e.g. ["10cm"]) to only generate one.
# "10cm" (vs. destructive reference) tells us about accuracy; "none"
# (fully unfiltered, vs. AdQSM) is closer to the actual geometry that
# will go into ANSYS for the production beech trees (which won't be
# filtered to >=10cm) - both useful for different reasons. This drives
# which reference row is used per branch_filter (see get_reference()).
FAMILY = "adtree"          # "adtree" or "treeqsm" - which parameter family to analyze;
                             # NEVER mixed in one chart, same principle as BRANCH_FILTER
METRIC = "n_cylinders"       # one output field at a time for the facet grid / ranked plot
                             # (correlation heatmap still checks all of METRICS_FOR_CORR at once)

FACET_Y_MODE = "absolute"   # "absolute" = raw METRIC value on the facet grid's y-axis (default,
                              # unchanged behaviour) - reference line drawn at the reference's own
                              # absolute value. "relative" = % deviation from the resolved reference
                              # (same formula plot_ranked_deviation() uses) - reference line then
                              # sits at 0% instead. Only affects plot_facet_grid() - plot_ranked_deviation()
                              # is already relative-only, plot_correlation_heatmap() is unaffected either way.

# --- adtree family settings ---
ADTREE_FACET_PARAM = "adqsm_variant"
ADTREE_X_PARAM = "radius_threshold_mm"
ADTREE_LINE_PARAM = "seg_min_mm"
ADTREE_PARAMS_FOR_CORR = ["adqsm_variant", "radius_threshold_mm", "seg_min_mm", "seg_max_mm", "seg_k_pct"]

# Dedicated diagnostic view (Bára's request): how does seg_min_mm
# actually shape branch_m3 deviation, given ANOVA found a real but
# non-monotonic effect Spearman couldn't see? Faceted by
# radius_threshold_mm (only 4 distinct values - readable panel
# count), x-axis = seg_min_mm, lines coloured by adqsm_variant.
SEGMIN_VIEW_METRIC = "branch_m3"
SEGMIN_VIEW_FACET_PARAM = "radius_threshold_mm"
SEGMIN_VIEW_X_PARAM = "seg_min_mm"
SEGMIN_VIEW_LINE_PARAM = "adqsm_variant"

# --- treeqsm family settings ---
TREEQSM_FACET_PARAM = "mode"
TREEQSM_X_PARAM = "simp_smallradii"
TREEQSM_LINE_PARAM = "simp_replaceiterations"
TREEQSM_PARAMS_FOR_CORR = ["pd1_m", "pd2min_m", "pd2max_m", "simp_maxorder", "simp_smallradii", "simp_replaceiterations"]
# NOTE: "mode" is a string ("manual"/"auto"), not numeric - deliberately
# NOT in TREEQSM_PARAMS_FOR_CORR above (Spearman needs an ordinal/numeric
# axis), but it's fine as TREEQSM_FACET_PARAM (the facet grid just needs
# distinct groups, not an ordering).

METRICS_FOR_CORR = ["total_m3", "trunk_m3", "branch_m3", "dbh_m"]

CORR_HEATMAP_MODES = ["raw", "pct_dev"]
# "raw" = original behavior: pools each tree's raw metric values
# directly (scale NOT normalized across trees - kept for direct
# comparison against the newer mode, not because it's recommended).
# "pct_dev" = pools each tree's % deviation from ITS OWN reference
# (via _pct_dev_pooled(), same normalization tornado/ranked-deviation
# already use) - removes the cross-tree scale confound the "raw"
# mode has. Both run automatically; set to a single-element list to
# only generate one.

TORNADO_METRICS = ["total_m3", "trunk_m3", "branch_m3", "n_cylinders"]
# Metrics shown as tornado diagrams (see plot_tornado_multi_tree()).
# The first three compare against each tree's destructive reference
# (same % formula as the ranked-deviation/facet-grid charts); n_cylinders
# has no reference row, so it instead compares against each tree's
# OWN MEAN n_cylinders across its selected rows - "how far from this
# tree's typical cylinder count", not "how far from a truth value".

# --- TreeQSM simp_smallradii/simp_replaceiterations sensitivity -------
# Confirmed directly from volume_results.csv (not assumed): within the
# "Simplified (no islands)" stage ONLY (the one stage that actually goes
# into ANSYS - "Optimal" and "Filtered <10cm" are other TreeQSM stages,
# not analyzed here), simp_smallradii/simp_replaceiterations have a
# genuine 6-combo sweep - (0.005,0),(0.005,1),(0.005,2),(0.010,0),
# (0.015,0),(0.020,0) - identical for mode="manual" and mode="auto",
# identical across all 3 reference trees. pd1_m/pd2min_m/pd2max_m are
# CONSTANT within each (tree, mode) pair (no sweep) - deliberately NOT
# included here; they're handled separately in a later manual-vs-auto
# reliability comparison task.
TREEQSM_STAGE_FOR_SENSITIVITY = "Simplified (no islands)"
TREEQSM_SIMP_PARAMS = ["simp_smallradii", "simp_replaceiterations"]
TREEQSM_SIMP_METRICS = ["total_m3", "trunk_m3", "branch_m3", "n_cylinders", "dbh_m"]
RUN_TREEQSM_SIMP_ANALYSIS = True   # independent of FAMILY - this
# block always targets family="treeqsm" regardless of what FAMILY is
# currently set to elsewhere in this file, since it's a standalone
# supplementary analysis, not part of the generic per-FAMILY loop.
# =====================================================================


# PARAMETERS-block field name (written in the more-familiar CSV-header
# spelling, e.g. "total_m3"/"pd1_m") -> load_results()'s own row dict key
# (which strips those "_m3"/"_m" suffixes, e.g. "total"/"pd1") - confirmed
# by inspecting a live load_results() row before writing this script.
# adqsm_variant/radius_threshold_mm/seg_*/mode/simp_* already match their
# row-dict key 1:1, so only the volume/DBH/PD fields need translating;
# _rowkey()'s .get() default handles those pass-through cases.
_FIELD_TO_ROWKEY = {
    "total_m3": "total", "trunk_m3": "trunk", "branch_m3": "branch", "std_m3": "std",
    "dbh_m": "dbh", "height_m": "height", "taper_cm_per_m": "taper",
    "trunk_len_m": "trunk_len", "branch_len_m": "branch_len",
    "pd1_m": "pd1", "pd2min_m": "pd2min", "pd2max_m": "pd2max", "mincylrad_m": "mincylrad",
}


def _rowkey(field_name):
    """Translate one PARAMETERS-block field name to its load_results() row
    dict key - see _FIELD_TO_ROWKEY above."""
    return _FIELD_TO_ROWKEY.get(field_name, field_name)


def get_reference(rows, branch_filter):
    """Returns (ref_row, ref_method, ref_color) - the SAME reference every
    other chart in this project compares against for this branch_filter:
    the destructive field reference (pink, FAMILY_GRADIENTS["Reference"][1])
    for "10cm", the dynamically-resolved smallest-numbered AdQSM variant
    (teal, plot_box.py's own TREEQSM_REF_LINE_COLOR - imported, not
    duplicated as a hex literal) for "none".

    `rows` should already be filtered to ONE tree before calling this -
    resolve_reference_method_none() picks the smallest variant PRESENT IN
    `rows`, so passing multi-tree rows here would resolve across trees,
    which no other chart in this project does.
    """
    if branch_filter == "10cm":
        method = REFERENCE_METHOD
        color = FAMILY_GRADIENTS["Reference"][1]
    else:
        method = resolve_reference_method_none(rows)
        color = TREEQSM_REF_LINE_COLOR
    ref_row = next((r for r in rows if r["method"] == method), None) if method else None
    return ref_row, method, color


def select_rows(rows, family, tree, branch_filter):
    """Rows for ONE family/tree/branch_filter - the population every chart
    in this script analyzes.

    "adtree": "AdTree calibrated ..." rows only (mirrors plot_box.py's
    assign_adtree_groups() default scope) - "AdTree raw" rows are excluded
    since they don't vary by the AdQSM-dependent parameters this script
    analyzes (adqsm_variant has no effect on raw AdTree geometry at all -
    see adtree_reconstruct_compare.py's STEP 5 fix, which correctly leaves
    adqsm_variant blank on those rows). Can be revisited if a future
    version should include raw AdTree too.

    "treeqsm": "TreeQSM mine (...)" rows that ALSO have structured params
    (mode != "") - the SAME check plot_box.py's assign_treeqsm_groups()
    already uses to exclude older pre-params_*.csv rows.

    The AdQSM (TreesParams) reference row itself is deliberately NEVER
    included here for either family - it's used ONLY as get_reference()'s
    comparison target, never as an analyzed parameter row (it has no
    structured AdTree OR TreeQSM columns to plot in the first place).
    """
    base = [r for r in rows if r["tree"] == tree and r["branch_filter"] == branch_filter]
    if family == "adtree":
        return [r for r in base if r["method"].startswith("AdTree calibrated")]
    if family == "treeqsm":
        return [r for r in base if parse_treeqsm_method(r["method"]) is not None and r["mode"]]
    raise ValueError("family must be 'adtree' or 'treeqsm', got %r" % family)


def filter_treeqsm_stage(rows, stage_substring):
    """Rows whose method contains `stage_substring` (e.g. 'Simplified
    (no islands)') - a plain substring check against the full method
    string, since TreeQSM stage names appear verbatim inside it (see
    parse_treeqsm_method-style parsing already used elsewhere in this
    file for the exact shape). Does not otherwise filter by tree/
    branch_filter/mode - callers still go through select_rows() for that."""
    return [r for r in rows if stage_substring in r["method"]]


def _family_name_for(rows):
    """classify_family() needs one representative method string - any row
    in an already family-filtered list works, they all share one family by
    construction (select_rows() only ever returns rows from ONE family)."""
    if not rows:
        return None
    return classify_family(rows[0]["method"], REFERENCE_METHOD)


def _gradient_color_map(values, family_name):
    """{value: color}, spreading `values` evenly across
    FAMILY_GRADIENTS[family_name] - the SAME t = i/(n-1)-across-a-gradient
    formula plot_box.py's build_group_color_map() uses, applied directly
    here (not by calling that function, which expects a groups/
    family_of_group shape this script doesn't build) rather than inventing
    a second colour formula. Falls back to flat grey for an unrecognized
    family rather than crashing."""
    stops = FAMILY_GRADIENTS.get(family_name, ["#9a9a9a", "#9a9a9a"])
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "%s_gradient" % (family_name or "unknown").lower().replace(" ", "_"), stops)
    n = len(values)
    return {v: cmap(i / (n - 1) if n > 1 else 0.5) for i, v in enumerate(values)}


def _fmt_num(v):
    """Display a numeric value without a spurious trailing ".0" for whole
    numbers (e.g. radius_threshold_mm=5.0 -> "5") - purely cosmetic, for
    the "closest combination" print line below."""
    if v is None:
        return "None"
    if isinstance(v, float) and v == int(v):
        return str(int(v))
    return str(v)


def _output_dir(tree):
    """plots/<tree>/ - same per-tree subfolder convention
    adtree_reconstruct_compare.py's own FIGURES_DIR already established,
    built on top of ensure_plots_dir() (plots/) rather than duplicating it."""
    out_dir = os.path.join(ensure_plots_dir(), tree)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def _save(fig, family, branch_filter, metric, chart, tree=None):
    """Save into plots/<tree>/sensitivity_<family>_<branch_filter>_<metric>_<chart>.png
    when `tree` is given (single-tree charts, unchanged from before this
    tree parameter existed), or directly into plots/ as
    sensitivity_multitree_<family>_<branch_filter>_<metric>_<chart>.png
    when `tree` is None (the pooled-across-trees functions) - the
    "multitree" filename tag keeps pooled and per-tree outputs from ever
    colliding or overwriting each other."""
    if tree is not None:
        out_dir = _output_dir(tree)
        out_path = os.path.join(out_dir, "sensitivity_%s_%s_%s_%s.png" % (family, branch_filter, metric, chart))
    else:
        out_path = os.path.join(
            ensure_plots_dir(), "sensitivity_multitree_%s_%s_%s_%s.png" % (family, branch_filter, metric, chart))
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


# ----------------------------------------------------------------------
# 1) Facet grid: one panel per facet_param value, metric vs. x_param,
#    lines coloured by line_param.
# ----------------------------------------------------------------------
def plot_facet_grid(rows, family, branch_filter, metric, facet_param, x_param, line_param, tree, filename_tag=None):
    sel = select_rows(rows, family, tree, branch_filter)
    if not sel:
        print("No rows for family=%r, tree=%r, branch_filter=%r - skipping facet grid."
              % (family, tree, branch_filter))
        return

    tree_rows = [r for r in rows if r["tree"] == tree and r["branch_filter"] == branch_filter]
    ref_row, ref_method, ref_color = get_reference(tree_rows, branch_filter)
    metric_key, facet_key, x_key, line_key = (
        _rowkey(metric), _rowkey(facet_param), _rowkey(x_param), _rowkey(line_param))
    ref_val = ref_row[metric_key] if ref_row is not None else None

    # FACET_Y_MODE == "relative": same "% deviation from reference" formula
    # plot_ranked_deviation() uses, applied per-point here instead - guarded
    # the same way that function guards it (no reference value, or a
    # reference value of 0 -> skip the whole chart rather than raising
    # ZeroDivisionError or plotting against a missing reference).
    relative = FACET_Y_MODE == "relative"
    if relative and (ref_val is None or ref_val == 0):
        print("FACET_Y_MODE='relative' but no usable reference value for '%s' (branch_filter='%s') "
              "- skipping facet grid." % (metric, branch_filter))
        return

    def y_of(row):
        v = row[metric_key]
        if v is None:
            return None
        return 100.0 * (v - ref_val) / ref_val if relative else v

    facet_values = sorted({r[facet_key] for r in sel if r[facet_key] not in (None, "")})
    line_values = sorted({r[line_key] for r in sel if r[line_key] not in (None, "")})
    if not facet_values:
        print("No non-blank '%s' values in this selection - skipping facet grid." % facet_param)
        return

    family_name = _family_name_for(sel)
    line_colors = _gradient_color_map(line_values, family_name)

    n = len(facet_values)
    ncols = min(4, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.2 * ncols, 3.6 * nrows), squeeze=False)
    axes_flat = list(axes.flat)
    for ax in axes_flat[n:]:
        ax.axis("off")

    legend_handles = {}
    for ax, fval in zip(axes_flat, facet_values):
        panel_rows = [r for r in sel if r[facet_key] == fval]
        for lval in line_values:
            sub = sorted(
                (r for r in panel_rows
                 if r[line_key] == lval and r[metric_key] is not None and r[x_key] is not None),
                key=lambda r: r[x_key])
            if not sub:
                continue
            xs = [r[x_key] for r in sub]
            ys = [y_of(r) for r in sub]
            line, = ax.plot(xs, ys, marker="o", markersize=4, color=line_colors[lval], label=str(lval))
            legend_handles.setdefault(lval, line)
        # Reference line: at 0% in relative mode (the reference IS the
        # zero-point by definition), at its own absolute value otherwise -
        # unaffected either way when ref_val is None (nothing to draw).
        if ref_val is not None:
            ax.axhline(0.0 if relative else ref_val, linestyle="--", color=ref_color, linewidth=1.2)
        ax.set_title("%s = %s" % (facet_param, fval), fontsize=10)
        ax.set_xlabel(x_param, fontsize=8)
        ax.set_ylabel(("%s (%% dev. from reference)" % metric) if relative else metric, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(False)

    if legend_handles:
        ordered_vals = [v for v in line_values if v in legend_handles]
        fig.legend([legend_handles[v] for v in ordered_vals], [str(v) for v in ordered_vals],
                   title=line_param, loc="lower center",
                   ncol=min(len(ordered_vals), 6), fontsize=8, bbox_to_anchor=(0.5, -0.02))

    fig.suptitle(
        "%s parameter sensitivity: %s  (%s, branch_filter='%s')\nReference: %s"
        % (family, metric, tree, branch_filter, ref_method), fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 0.90))
    # Mode suffix on THIS chart's filename only (ranked/heatmap have no
    # such toggle) - so absolute/relative runs of the same
    # family/branch_filter/metric don't silently overwrite each other.
    # filename_tag (optional): appended when a metric/x_param/facet_param
    # combination reuses an existing METRIC value (e.g. the SEGMIN_VIEW_*
    # combination also uses "branch_m3") - without it, that combination's
    # PNG would silently collide with/overwrite a differently-configured
    # facet grid using the same metric. None (default) = today's exact
    # filename, unchanged.
    chart = ("facetgrid_%s_%s" % (FACET_Y_MODE, filename_tag) if filename_tag is not None
             else "facetgrid_%s" % FACET_Y_MODE)
    _save(fig, family, branch_filter, metric, chart, tree=tree)


# ----------------------------------------------------------------------
# 2) Ranked deviation dot plot: every individual row, sorted by |% dev|,
#    dashed line at 0% - also a practical "which combo is closest" leaderboard.
# ----------------------------------------------------------------------
def plot_ranked_deviation(rows, family, branch_filter, metric, tree):
    sel = select_rows(rows, family, tree, branch_filter)
    if not sel:
        print("No rows for family=%r, tree=%r, branch_filter=%r - skipping ranked deviation plot."
              % (family, tree, branch_filter))
        return

    tree_rows = [r for r in rows if r["tree"] == tree and r["branch_filter"] == branch_filter]
    ref_row, ref_method, ref_color = get_reference(tree_rows, branch_filter)
    metric_key = _rowkey(metric)
    if ref_row is None or ref_row[metric_key] is None:
        print("No reference value for '%s' (branch_filter='%s') - skipping ranked deviation plot."
              % (metric, branch_filter))
        return
    ref_val = ref_row[metric_key]
    if ref_val == 0:
        print("Reference value for '%s' is 0 - percent deviation is undefined, skipping ranked deviation plot."
              % metric)
        return

    scored = [(100.0 * (r[metric_key] - ref_val) / ref_val, r) for r in sel if r[metric_key] is not None]
    if not scored:
        print("No row has a value for '%s' - skipping ranked deviation plot." % metric)
        return
    scored.sort(key=lambda s: abs(s[0]))

    family_name = _family_name_for(sel)
    stops = FAMILY_GRADIENTS.get(family_name, ["#5a92d6"])
    base_color = stops[len(stops) // 2]

    # Shortened labels - same shorten_method_label() plot_volumes.py's own
    # charts (and plot_box.py's point labels) already use everywhere else,
    # instead of the long raw method string.
    labels = [shorten_method_label(s[1]["method"]) for s in scored]
    values = [s[0] for s in scored]
    y = list(range(len(scored)))

    fig, ax = plt.subplots(figsize=(9, max(4, 0.22 * len(scored))))
    ax.axvline(0.0, linestyle="--", color=ref_color, linewidth=1.2, label="Reference: %s" % ref_method)
    ax.scatter(values, y, color=base_color, zorder=3, s=25)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.invert_yaxis()   # closest match (top of the sorted list) at the TOP of the chart
    ax.set_xlabel("%% deviation from reference (%s)" % metric)
    ax.set_title("%s ranked by deviation from reference: %s\n(%s, branch_filter='%s')"
                 % (family, metric, tree, branch_filter))
    ax.legend(fontsize=8, loc="best")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    _save(fig, family, branch_filter, metric, "ranked", tree=tree)

    best_row = scored[0][1]
    if family == "adtree":
        detail = "adqsm_variant=%s, radius_threshold_mm=%s, seg_k_pct=%s" % (
            best_row.get("adqsm_variant"), _fmt_num(best_row.get("radius_threshold_mm")),
            _fmt_num(best_row.get("seg_k_pct")))
    elif family == "treeqsm":
        detail = "mode=%s, simp_smallradii=%s, simp_replaceiterations=%s" % (
            best_row.get("mode"), _fmt_num(best_row.get("simp_smallradii")),
            _fmt_num(best_row.get("simp_replaceiterations")))
    else:
        detail = ""
    print("Closest to reference (%s): %s  [%s] -> %+.2f%%"
          % (metric, shorten_method_label(best_row["method"]), detail, values[0]))


# ----------------------------------------------------------------------
# 3) Correlation heatmap: Spearman correlation between each structured
#    parameter and each output metric - fast numeric screening ONLY, see
#    the module docstring for why this understates categorical effects.
# ----------------------------------------------------------------------
def _spearman(x, y):
    """Spearman rank correlation between two equal-length numeric lists -
    computed as Pearson correlation on each list's own ranks (no SciPy
    dependency needed; mathematically identical to Spearman's rho)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or np.all(x == x[0]) or np.all(y == y[0]):
        return None
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def _eta_squared(values, group_labels):
    """One-way ANOVA effect size: what fraction of the total variance in
    `values` is explained by grouping them according to `group_labels`
    (one label per value, e.g. that row's parameter value - equality-
    based grouping, no numeric casting needed, unlike tornado's min/max
    ordering which specifically needs numeric comparison). Returns None
    if fewer than 2 distinct groups, or if every value is identical
    (SS_total == 0, effect size undefined). eta_squared = SS_between /
    SS_total, in [0, 1] - higher means that parameter's grouping explains
    more of the spread in `values`. Main-effects only (this parameter
    alone) - no interaction terms, consistent with this file's existing
    Spearman heatmap being a numeric-screening supplement, not a full
    factorial ANOVA."""
    values = np.asarray(values, dtype=float)
    groups = {}
    for v, g in zip(values, group_labels):
        groups.setdefault(g, []).append(v)
    if len(groups) < 2:
        return None
    grand_mean = values.mean()
    ss_total = float(np.sum((values - grand_mean) ** 2))
    if ss_total == 0:
        return None
    ss_between = sum(len(vs) * (np.mean(vs) - grand_mean) ** 2 for vs in groups.values())
    return float(ss_between / ss_total)


def plot_correlation_heatmap(rows, family, branch_filter, params, metrics, tree):
    sel = select_rows(rows, family, tree, branch_filter)
    if not sel:
        print("No rows for family=%r, tree=%r, branch_filter=%r - skipping correlation heatmap."
              % (family, tree, branch_filter))
        return

    matrix = np.full((len(params), len(metrics)), np.nan)
    for i, p in enumerate(params):
        p_key = _rowkey(p)
        for j, m in enumerate(metrics):
            m_key = _rowkey(m)
            pairs = [(r[p_key], r[m_key]) for r in sel
                     if r[p_key] not in (None, "") and r[m_key] is not None]
            if len(pairs) < 2:
                continue
            xs, ys = zip(*pairs)
            rho = _spearman(xs, ys)
            if rho is not None:
                matrix[i, j] = rho

    fig, ax = plt.subplots(figsize=(1.3 * len(metrics) + 2, 0.6 * len(params) + 2))
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=30, ha="right")
    ax.set_yticks(range(len(params)))
    ax.set_yticklabels(params)
    for i in range(len(params)):
        for j in range(len(metrics)):
            v = matrix[i, j]
            if not np.isnan(v):
                ax.text(j, i, "%.2f" % v, ha="center", va="center",
                        color="white" if abs(v) > 0.6 else "black", fontsize=8)
    fig.colorbar(im, ax=ax, label="Spearman rho")
    ax.set_title(
        "%s parameter/output correlation (Spearman)  (%s, branch_filter='%s')\n"
        "NOTE: understates non-monotonic/categorical effects (e.g. AdQSM variant identity) -\n"
        "the facet grid and ranked plot are the primary tools, this is a numeric-screening supplement."
        % (family, tree, branch_filter), fontsize=9)
    fig.tight_layout()
    _save(fig, family, branch_filter, "allmetrics", "corrheatmap", tree=tree)


def _pct_dev_pooled(rows, family, branch_filter, metric, trees, allow_own_mean_fallback=True,
                     reference_rows=None):
    """For each tree in `trees`, resolve that tree's own baseline for
    `metric`: get_reference()'s row value if one exists and is
    non-zero (same as every other %-deviation chart in this file),
    otherwise (only when `allow_own_mean_fallback` is True) the MEAN of
    that tree's own select_rows() population for this metric (only
    reached for metrics with no reference row, e.g. n_cylinders -
    printed once as a NOTE so the fallback is visible, not silent).
    `allow_own_mean_fallback=False` restores the older "no reference ->
    skip this tree" behavior instead - plot_ranked_deviation_multi_tree()
    below passes False specifically so gaining this fallback doesn't
    change ITS existing behavior for reference-less metrics (it was never
    designed to plot "deviation from own mean", only "deviation from
    reference" - plot_tornado_multi_tree() is the one chart meant to
    handle both).

    `reference_rows`: where to look up the reference row FROM, separate
    from `rows` (the analyzed population). None (default) = look for the
    reference inside `rows` itself, exactly as before this parameter
    existed. Needed because `rows` can be pre-filtered to a SUBSET that
    excludes the reference row entirely - e.g. the TreeQSM
    simp_smallradii/simp_replaceiterations sensitivity analysis passes a
    `rows` already restricted to the "Simplified (no islands)" stage,
    which never includes the destructive-reference/AdQSM row (neither is
    itself a "Simplified (no islands)" TreeQSM row) - without a separate
    `reference_rows` pointing at the FULL dataset, get_reference() would
    silently find nothing there and this function would incorrectly fall
    back to (or skip past) the "own mean" path even for metrics that
    genuinely have a usable reference. Returns (pooled, baseline_kind):
      - pooled: a flat list of (tree, row, pct_deviation) triples pooled
        across all `trees` - the same shape
        plot_ranked_deviation_multi_tree() already built before this
        helper existed (reused here, not reinvented).
      - baseline_kind: "reference" or "own_mean" - whichever baseline
        path was actually used. Every metric in this file uses ONE
        consistent path across every tree (a metric either always has a
        reference row, or never does), so a single flag describes the
        whole pooled result.
    Skips a tree entirely (with a printed reason) if no baseline could be
    resolved at all (e.g. empty selection, reference missing/zero with
    the fallback disabled, or both reference and own mean unusable)."""
    metric_key = _rowkey(metric)
    pooled = []
    baseline_kind = None
    ref_source = reference_rows if reference_rows is not None else rows
    for t in trees:
        sel = select_rows(rows, family, t, branch_filter)
        if not sel:
            print("No rows for family=%r, tree=%r, branch_filter=%r - skipping this tree."
                  % (family, t, branch_filter))
            continue
        tree_rows = [r for r in ref_source if r["tree"] == t and r["branch_filter"] == branch_filter]
        ref_row, ref_method, ref_color = get_reference(tree_rows, branch_filter)

        baseline_val = None
        this_kind = None
        if ref_row is not None and ref_row[metric_key] is not None and ref_row[metric_key] != 0:
            baseline_val = ref_row[metric_key]
            this_kind = "reference"
        elif allow_own_mean_fallback:
            own_vals = [r[metric_key] for r in sel if r[metric_key] is not None]
            if own_vals:
                own_mean = sum(own_vals) / len(own_vals)
                if own_mean != 0:
                    baseline_val = own_mean
                    this_kind = "own_mean"
                    print("NOTE: no usable reference value for '%s' (tree=%r) - using this tree's own "
                          "mean (%.4f) as the baseline instead." % (metric, t, own_mean))

        if baseline_val is None:
            print("No baseline (reference or own mean) could be resolved for '%s' (tree=%r, "
                  "branch_filter='%s') - skipping this tree." % (metric, t, branch_filter))
            continue

        baseline_kind = baseline_kind or this_kind
        for r in sel:
            if r[metric_key] is not None:
                pooled.append((t, r, 100.0 * (r[metric_key] - baseline_val) / baseline_val))

    return pooled, baseline_kind


# ----------------------------------------------------------------------
# 2b) Ranked deviation, POOLED across several trees - each tree's rows are
#     scored against THAT tree's own reference (same formula as
#     plot_ranked_deviation()), then every tree's scored rows are combined
#     into ONE ranking. Added alongside plot_ranked_deviation() (not a
#     replacement) - see this file's own module docstring/RUN section for
#     which one is actually used where.
# ----------------------------------------------------------------------
def plot_ranked_deviation_multi_tree(rows, family, branch_filter, metric, trees, filename_tag=None,
                                      reference_rows=None):
    # _pct_dev_pooled() is the same pooling logic this function used to
    # build inline, extracted so plot_tornado_multi_tree() below can reuse
    # it too. allow_own_mean_fallback=False keeps THIS function's own
    # behavior byte-for-byte unchanged by that extraction: it was always
    # "deviation from reference, skip a tree/metric with none" (never
    # "deviation from own mean") - plot_tornado_multi_tree() is the one
    # chart that actually wants that fallback (e.g. for n_cylinders).
    # reference_rows: see _pct_dev_pooled()'s own docstring - None
    # (default) preserves today's exact behavior for every existing call
    # site (reference looked up inside `rows` itself).
    pooled, _baseline_kind = _pct_dev_pooled(rows, family, branch_filter, metric, trees,
                                              allow_own_mean_fallback=False,
                                              reference_rows=reference_rows)
    if not pooled:
        print("No row (across %d trees) has a value for '%s' - skipping pooled ranked deviation plot."
              % (len(trees), metric))
        return
    pooled.sort(key=lambda p: abs(p[2]))

    # tree_marker_map built from `trees` (sorted, for a stable assignment)
    # - same {tree: TREE_MARKERS[i % len(TREE_MARKERS)]} pattern
    # plot_volumes.py's plot_error_boxplot() uses for its own per-tree
    # marker legend.
    tree_marker_map = {t: TREE_MARKERS[i % len(TREE_MARKERS)] for i, t in enumerate(sorted(trees))}
    # Fixed by branch_filter alone (see get_reference()) - computed
    # directly rather than depending on a specific tree's own
    # get_reference() call, since every tree resolves to the same colour.
    ref_color = FAMILY_GRADIENTS["Reference"][1] if branch_filter == "10cm" else TREEQSM_REF_LINE_COLOR

    family_name = _family_name_for([r for (_, r, _) in pooled])
    stops = FAMILY_GRADIENTS.get(family_name, ["#5a92d6"])
    base_color = stops[len(stops) // 2]

    # Shortened labels, same as plot_ranked_deviation() - PLUS a " [tree]"
    # suffix, since the same method/parameter combo can now appear once per
    # tree, and without the suffix two rows for the same combo (different
    # trees) would show identical y-axis text.
    labels = ["%s [%s]" % (shorten_method_label(r["method"]), t) for (t, r, _) in pooled]
    values = [p for (_, _, p) in pooled]
    y = list(range(len(pooled)))

    fig, ax = plt.subplots(figsize=(9, max(4, 0.22 * len(pooled))))
    ax.axvline(0.0, linestyle="--", color=ref_color, linewidth=1.2, label="Reference (0%)")
    for i, (t, r, p) in enumerate(pooled):
        ax.plot(p, i, marker=tree_marker_map[t], color=base_color, markersize=6,
                markeredgewidth=0, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.invert_yaxis()   # closest match (top of the sorted list) at the TOP of the chart
    ax.set_xlabel("%% deviation from reference (%s)" % metric)
    ax.set_title("%s ranked by deviation from reference: %s\n(%d trees pooled, branch_filter='%s')"
                 % (family, metric, len(trees), branch_filter))
    ax.grid(axis="x", alpha=0.3)

    # TWO separate legends: the reference axvline (kept inside the plot
    # area, same as plot_ranked_deviation()'s own legend), plus a second
    # one mapping marker SHAPE -> tree name, placed OUTSIDE the axes (to
    # the right) so it never overlaps the dots - same approach as
    # plot_volumes.py's plot_error_boxplot() tree legend. ax.add_artist()
    # is needed because a second ax.legend() call would otherwise replace
    # (not add to) the first.
    ref_legend = ax.legend(fontsize=8, loc="best")
    ax.add_artist(ref_legend)
    tree_legend_handles = [
        mlines.Line2D([0], [0], marker=tree_marker_map[t], color="none",
                      markerfacecolor="#555555", markeredgewidth=0,
                      markersize=8, label=t)
        for t in sorted(trees)
    ]
    ax.legend(handles=tree_legend_handles, title="Tree", loc="upper left",
              bbox_to_anchor=(1.01, 1.0), borderaxespad=0)

    fig.tight_layout()
    # filename_tag (optional, same convention as plot_facet_grid()'s own
    # filename_tag): None (default) preserves today's exact filename;
    # given, appends "_<filename_tag>" so e.g. the TreeQSM simp_* analysis
    # (simp_manual/simp_auto/simp_pooled) can never collide with/overwrite
    # a plain call using the same metric/branch_filter/family.
    chart = "ranked_%s" % filename_tag if filename_tag is not None else "ranked"
    _save(fig, family, branch_filter, metric, chart, tree=None)

    best_tree, best_row, best_pct = pooled[0]
    if family == "adtree":
        detail = "adqsm_variant=%s, radius_threshold_mm=%s, seg_k_pct=%s" % (
            best_row.get("adqsm_variant"), _fmt_num(best_row.get("radius_threshold_mm")),
            _fmt_num(best_row.get("seg_k_pct")))
    elif family == "treeqsm":
        detail = "mode=%s, simp_smallradii=%s, simp_replaceiterations=%s" % (
            best_row.get("mode"), _fmt_num(best_row.get("simp_smallradii")),
            _fmt_num(best_row.get("simp_replaceiterations")))
    else:
        detail = ""
    print("Closest to reference (%s): %s  [tree=%s]  [%s] -> %+.2f%%"
          % (metric, shorten_method_label(best_row["method"]), best_tree, detail, best_pct))


# ----------------------------------------------------------------------
# 3b) Correlation heatmap, POOLED across several trees - TWO modes (see
#     CORR_HEATMAP_MODES in the PARAMETERS block), kept SIDE BY SIDE (not
#     one replacing the other), same project convention as e.g. the
#     calmethod min5mm/regression-perorder comparison:
#       "raw"     - the RAW (param_value, metric_value) pairs from every
#                   tree are concatenated FIRST, then ONE Spearman rho is
#                   computed on the combined set (NOT a per-tree rho
#                   averaged afterwards). This is the ORIGINAL behavior -
#                   UNCHANGED code path, still pools raw values across
#                   trees of different absolute scale.
#       "pct_dev" - the same idea, but on (param_value, pct_deviation)
#                   pairs instead - pct_deviation comes from
#                   _pct_dev_pooled() (same normalization tornado/
#                   ranked-deviation already use), which removes the
#                   cross-tree SCALE confound the "raw" mode has (3 trees
#                   of very different absolute volume/cylinder-count are
#                   no longer pooled as if they were one population).
# ----------------------------------------------------------------------
def plot_correlation_heatmap_multi_tree(rows, family, branch_filter, params, metrics, trees, mode="raw",
                                         filename_tag=None, reference_rows=None):
    if mode not in ("raw", "pct_dev"):
        raise ValueError("mode must be 'raw' or 'pct_dev', got %r" % mode)

    if mode == "raw":
        # UNCHANGED from before CORR_HEATMAP_MODES existed - do not alter
        # this code path (see this function's own docstring above).
        pooled_sel = []
        for t in trees:
            sel = select_rows(rows, family, t, branch_filter)
            if not sel:
                print("No rows for family=%r, tree=%r, branch_filter=%r - contributing nothing from this "
                      "tree to the pooled correlation heatmap." % (family, t, branch_filter))
            pooled_sel.extend(sel)

        if not pooled_sel:
            print("No rows (across %d trees) for family=%r, branch_filter=%r - skipping pooled correlation "
                  "heatmap." % (len(trees), family, branch_filter))
            return

        matrix = np.full((len(params), len(metrics)), np.nan)
        for i, p in enumerate(params):
            p_key = _rowkey(p)
            for j, m in enumerate(metrics):
                m_key = _rowkey(m)
                pairs = [(r[p_key], r[m_key]) for r in pooled_sel
                         if r[p_key] not in (None, "") and r[m_key] is not None]
                if len(pairs) < 2:
                    continue
                xs, ys = zip(*pairs)
                rho = _spearman(xs, ys)
                if rho is not None:
                    matrix[i, j] = rho
    else:
        # pct_dev: _pct_dev_pooled() called ONCE PER METRIC (not once per
        # param - mirrors the raw path's "once per tree" call pattern,
        # just keyed by metric instead since the pooled pct_deviation
        # values are metric-specific), then EVERY param's column reuses
        # that same per-metric pooled list.
        # reference_rows: see _pct_dev_pooled()'s own docstring - None
        # (default) preserves today's exact behavior for every existing
        # call site (reference looked up inside `rows` itself).
        pooled_by_metric = {}
        for m in metrics:
            pooled_m, _baseline_kind = _pct_dev_pooled(rows, family, branch_filter, m, trees,
                                                        reference_rows=reference_rows)
            pooled_by_metric[m] = pooled_m

        if not any(pooled_by_metric.values()):
            print("No rows (across %d trees) for family=%r, branch_filter=%r - skipping pooled pct_dev "
                  "correlation heatmap." % (len(trees), family, branch_filter))
            return

        matrix = np.full((len(params), len(metrics)), np.nan)
        for i, p in enumerate(params):
            p_key = _rowkey(p)
            for j, m in enumerate(metrics):
                pairs = [(row[p_key], pct) for (_, row, pct) in pooled_by_metric[m]
                         if row[p_key] not in (None, "")]
                if len(pairs) < 2:
                    continue
                xs, ys = zip(*pairs)
                rho = _spearman(xs, ys)
                if rho is not None:
                    matrix[i, j] = rho

    fig, ax = plt.subplots(figsize=(1.3 * len(metrics) + 2, 0.6 * len(params) + 2))
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=30, ha="right")
    ax.set_yticks(range(len(params)))
    ax.set_yticklabels(params)
    for i in range(len(params)):
        for j in range(len(metrics)):
            v = matrix[i, j]
            if not np.isnan(v):
                ax.text(j, i, "%.2f" % v, ha="center", va="center",
                        color="white" if abs(v) > 0.6 else "black", fontsize=8)
    fig.colorbar(im, ax=ax, label="Spearman rho")
    # mode_line makes the raw-vs-pct_dev difference visible ON THE CHART
    # itself, not just in the filename/code.
    mode_line = ("(scale-normalized: % deviation from own reference)" if mode == "pct_dev"
                 else "(raw pooled values - NOT scale-normalized across trees)")
    ax.set_title(
        "%s parameter/output correlation (Spearman)  (%d trees pooled, branch_filter='%s')\n"
        "%s\n"
        "NOTE: understates non-monotonic/categorical effects (e.g. AdQSM variant identity) -\n"
        "the facet grid and ranked plot are the primary tools, this is a numeric-screening supplement."
        % (family, len(trees), branch_filter, mode_line), fontsize=9)
    fig.tight_layout()
    # Filename tags the mode ("_raw"/"_pctdev") - this DOES change the
    # "raw" mode's filename vs. before CORR_HEATMAP_MODES existed
    # (previously plain "corrheatmap", now "corrheatmap_raw") - the one
    # deliberate exception to "don't change existing filenames" in this
    # task, unavoidable now that two variants exist side by side.
    chart = "corrheatmap_raw" if mode == "raw" else "corrheatmap_pctdev"
    # filename_tag (optional, distinct from `mode`) - same convention as
    # plot_facet_grid()'s own filename_tag: None (default) preserves the
    # filename exactly as above; given, appends "_<filename_tag>".
    if filename_tag is not None:
        chart = "%s_%s" % (chart, filename_tag)
    _save(fig, family, branch_filter, "allmetrics", chart, tree=None)


def _fmt_param_value(param, value):
    """Display one structured parameter's raw VALUE (not a %% result) for
    the tornado diagram's end-of-bar annotations. Casts to float first
    (same "numeric-string" tolerance _spearman()'s callers rely on for
    values like adqsm_variant's "04") so _fmt_num() can drop a spurious
    trailing ".0" - then appends the column's own unit suffix ("mm" for
    every "*_mm" column, "%%" for "*_pct") so e.g. radius_threshold_mm=5.0
    reads as "5mm" on the chart, not a bare "5"."""
    try:
        num = float(value)
    except (TypeError, ValueError):
        num = None
    text = _fmt_num(num) if num is not None else str(value)
    if param.endswith("_mm"):
        return text + "mm"
    if param.endswith("_pct"):
        return text + "%"
    return text


# ----------------------------------------------------------------------
# 4) Tornado diagram, POOLED across several trees - for each structured
#    parameter, how far apart is the pooled mean %% deviation between its
#    MINIMUM and MAXIMUM present value? Bars sorted by that spread,
#    largest first (classic tornado ordering) - a quick "which parameter
#    actually moves the outcome" screening tool, complementary to the
#    correlation heatmap (monotonic-only) and facet grid (full detail,
#    one tree at a time).
# ----------------------------------------------------------------------
def plot_tornado_multi_tree(rows, family, branch_filter, metric, params, trees, filename_tag=None,
                             reference_rows=None):
    # reference_rows: see _pct_dev_pooled()'s own docstring - None
    # (default) preserves today's exact behavior for every existing call
    # site (reference looked up inside `rows` itself).
    pooled, baseline_kind = _pct_dev_pooled(rows, family, branch_filter, metric, trees,
                                             reference_rows=reference_rows)
    if not pooled:
        print("No row (across %d trees) has a value for '%s' - skipping tornado diagram."
              % (len(trees), metric))
        return

    bars = []   # (param, low_val, high_val, low_mean, high_mean)
    for p in params:
        p_key = _rowkey(p)
        raw_vals = {row[p_key] for (_, row, _) in pooled if row[p_key] not in (None, "")}
        if len(raw_vals) < 2:
            print("Fewer than 2 distinct values of '%s' in the pooled selection - skipping this "
                  "parameter in the tornado diagram." % p)
            continue
        # Cast to float for ORDERING only (same tolerance _spearman()'s
        # callers already rely on for numeric-string values like
        # adqsm_variant's "04") - grouping below still matches on the
        # ORIGINAL (uncast) value, not the float.
        ordered_vals = sorted(raw_vals, key=lambda v: float(v))
        min_val, max_val = ordered_vals[0], ordered_vals[-1]
        low_group = [pct for (_, row, pct) in pooled if row[p_key] == min_val]
        high_group = [pct for (_, row, pct) in pooled if row[p_key] == max_val]
        low_mean = sum(low_group) / len(low_group)
        high_mean = sum(high_group) / len(high_group)
        bars.append((p, min_val, max_val, low_mean, high_mean))

    if not bars:
        print("No parameter had 2+ distinct values in the pooled selection - skipping tornado "
              "diagram for '%s'." % metric)
        return

    # Classic tornado ordering: largest |high - low| spread first (most
    # impactful parameter), reversed onto the chart below via
    # invert_yaxis() so it lands at the TOP, not the bottom.
    bars.sort(key=lambda b: abs(b[4] - b[3]), reverse=True)

    family_name = _family_name_for([r for (_, r, _) in pooled])
    stops = FAMILY_GRADIENTS.get(family_name, ["#5a92d6"])
    base_color = stops[len(stops) // 2]

    fig, ax = plt.subplots(figsize=(9, max(3, 0.6 * len(bars) + 1)))
    y = list(range(len(bars)))

    # x-axis padding computed from the bars' actual data range (NOT
    # ax.margins(), which only pads relative to plotted artists and could
    # still leave an end-value annotation sitting right against the left
    # spine - close enough to visually collide with the y-axis tick text,
    # e.g. seg_min_mm's bar in the n_cylinders tornado chart) - set BEFORE
    # the annotations are placed below so there's always room on both sides.
    all_bar_vals = [v for (_, _, _, low_mean, high_mean) in bars for v in (low_mean, high_mean)]
    data_min, data_max = min(all_bar_vals), max(all_bar_vals)
    padding = 0.12 * (data_max - data_min) if data_max != data_min else 1.0
    ax.set_xlim(data_min - padding, data_max + padding)

    for i, (p, min_val, max_val, low_mean, high_mean) in enumerate(bars):
        left = min(low_mean, high_mean)
        width = abs(high_mean - low_mean)
        ax.barh(i, width, left=left, height=0.6, color=base_color, alpha=0.75)
        # Annotate each bar's two ends with the ACTUAL parameter values
        # (not the %% result) - whichever end (low_mean/high_mean) is
        # smaller gets the smaller-valued annotation, regardless of the
        # low/high VALUE's own sign or magnitude. Offset via
        # textcoords="offset points" (not a data-scale fraction) so the
        # text never overlaps the bar fill regardless of the chart's x range.
        left_val, right_val = (min_val, max_val) if low_mean <= high_mean else (max_val, min_val)
        ax.annotate(_fmt_param_value(p, left_val), xy=(left, i), xytext=(-6, 0),
                    textcoords="offset points", ha="right", va="center", fontsize=8)
        ax.annotate(_fmt_param_value(p, right_val), xy=(left + width, i), xytext=(6, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=8)
    ax.axvline(0.0, linestyle="--", color="#999999", linewidth=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([p for (p, *_rest) in bars])
    ax.invert_yaxis()   # most impactful parameter (sorted first) at the TOP
    # (x-axis limits/padding already set above, before the annotations were drawn)

    # xlabel reflects whichever baseline path _pct_dev_pooled() ACTUALLY
    # used (threaded back via baseline_kind), rather than re-deriving it
    # from the metric name a second time here - the two can never drift
    # out of sync this way.
    xlabel = ("% deviation from tree's own mean" if baseline_kind == "own_mean"
              else "% deviation from reference")
    ax.set_xlabel(xlabel)
    title = ("%s parameter impact (tornado): %s\n(%d trees pooled, branch_filter='%s')"
             % (family, metric, len(trees), branch_filter))
    if metric == "n_cylinders":
        title += "\nBaseline: each tree's own mean n_cylinders (no destructive reference for this field)"
    ax.set_title(title)
    fig.tight_layout()
    # filename_tag (optional, same convention as plot_facet_grid()'s own
    # filename_tag): None (default) preserves today's exact filename.
    chart = "tornado_%s" % filename_tag if filename_tag is not None else "tornado"
    _save(fig, family, branch_filter, metric, chart, tree=None)


# ----------------------------------------------------------------------
# 5) One-way ANOVA / variance decomposition (eta-squared) - for each
#    structured parameter, what FRACTION of the total variance in the
#    pooled %-deviation values is explained by that parameter's grouping
#    alone? Complementary to the tornado diagram (which only compares the
#    parameter's min vs. max value) and the Spearman heatmap (monotonic
#    relationships only) - this catches a parameter that matters a lot
#    but non-monotonically (e.g. a categorical like adqsm_variant, or a
#    U-shaped numeric effect neither of the other two charts would rank
#    highly). pct_dev-ONLY (no "raw" variant) - pooling raw values across
#    trees of different scale would bias eta-squared the same way it
#    biased the raw correlation heatmap (see that chart's own rationale
#    comment) - so this chart always normalizes via _pct_dev_pooled()
#    first, same as the tornado diagram.
# ----------------------------------------------------------------------
def plot_anova_multi_tree(rows, family, branch_filter, metric, params, trees, filename_tag=None,
                           reference_rows=None):
    # reference_rows: see _pct_dev_pooled()'s own docstring - None
    # (default) preserves today's exact behavior for every existing call
    # site (reference looked up inside `rows` itself).
    pooled, baseline_kind = _pct_dev_pooled(rows, family, branch_filter, metric, trees,
                                             reference_rows=reference_rows)
    if not pooled:
        print("No row (across %d trees) has a value for '%s' - skipping ANOVA chart."
              % (len(trees), metric))
        return

    bars = []   # (param, eta_squared)
    for p in params:
        p_key = _rowkey(p)
        vals, labels = [], []
        for (_, row, pct) in pooled:
            v = row.get(p_key)
            if v in (None, ""):
                continue
            vals.append(pct)
            labels.append(v)
        eta_sq = _eta_squared(vals, labels)
        if eta_sq is None:
            print("Could not compute eta-squared for '%s' (fewer than 2 groups, or zero variance) - "
                  "skipping this parameter in the ANOVA chart." % p)
            continue
        bars.append((p, eta_sq))

    if not bars:
        print("No parameter had a usable eta-squared in the pooled selection - skipping ANOVA "
              "chart for '%s'." % metric)
        return

    # Most-explaining parameter first, same visual convention as the
    # tornado diagram (largest effect at the TOP via invert_yaxis() below).
    bars.sort(key=lambda b: b[1], reverse=True)

    fig, ax = plt.subplots(figsize=(9, max(3, 0.6 * len(bars) + 1)))
    y = list(range(len(bars)))
    # Flat neutral color - this isn't a per-method chart (no color_map
    # needed) - "#5a92d6" is the same "AdTree calibrated" family colour
    # already used elsewhere in this file (e.g. plot_ranked_deviation()'s
    # own fallback base_color), not a new hex value.
    bar_color = "#5a92d6"
    ax.barh(y, [eta_sq for (_, eta_sq) in bars], height=0.6, color=bar_color, alpha=0.85)
    for i, (p, eta_sq) in enumerate(bars):
        ax.annotate("%.2f" % eta_sq, xy=(eta_sq, i), xytext=(6, 0),
                    textcoords="offset points", ha="left", va="center", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([p for (p, _eta_sq) in bars])
    ax.invert_yaxis()   # most-explaining parameter (sorted first) at the TOP
    ax.set_xlim(0, 1)   # eta-squared is bounded in [0, 1]
    ax.margins(x=0.1)   # headroom so the outward-offset annotations never get clipped
    ax.set_xlabel("Eta² (share of variance explained by this parameter)")

    title = ("%s parameter variance explained (ANOVA eta²): %s\n(%d trees pooled, branch_filter='%s')"
             % (family, metric, len(trees), branch_filter))
    if metric == "n_cylinders":
        title += "\nBaseline: each tree's own mean n_cylinders (no destructive reference for this field)"
    ax.set_title(title)
    fig.tight_layout()
    # filename_tag (optional, same convention as plot_facet_grid()'s own
    # filename_tag): None (default) preserves today's exact filename.
    chart = "anova_%s" % filename_tag if filename_tag is not None else "anova"
    _save(fig, family, branch_filter, metric, chart, tree=None)


# =========================  RUN  =====================================
# Calls all plotting functions once per BRANCH_FILTERS_TO_RUN entry, for
# the currently-set FAMILY/METRIC - does NOT loop over FAMILY or METRIC
# automatically (re-run by hand with different PARAMETER values for
# those; a batch-all-combinations mode can be a follow-up once this
# version's been reviewed).
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit("'%s' not found - run compare_volumes.py first so it gets created." % RESULTS_CSV)

    all_rows = load_results(RESULTS_CSV)

    if FAMILY == "adtree":
        facet_param, x_param, line_param = ADTREE_FACET_PARAM, ADTREE_X_PARAM, ADTREE_LINE_PARAM
        params_for_corr = ADTREE_PARAMS_FOR_CORR
    elif FAMILY == "treeqsm":
        facet_param, x_param, line_param = TREEQSM_FACET_PARAM, TREEQSM_X_PARAM, TREEQSM_LINE_PARAM
        params_for_corr = TREEQSM_PARAMS_FOR_CORR
    else:
        raise SystemExit("FAMILY must be 'adtree' or 'treeqsm', got %r" % FAMILY)

    # Outer loop over BRANCH_FILTERS_TO_RUN (see PARAMETERS block) - every
    # plotting function already takes branch_filter as an explicit
    # parameter, so this is a RUN-section-only change; no function
    # signatures needed touching for this part.
    for branch_filter in BRANCH_FILTERS_TO_RUN:
        # Facet grid stays per-tree - looped over REFERENCE_TREES so all 3
        # reference trees get their own facet grid PNG (plots/<tree>/...).
        for tree in REFERENCE_TREES:
            plot_facet_grid(all_rows, FAMILY, branch_filter, METRIC, facet_param, x_param, line_param, tree=tree)

        # Dedicated seg_min_mm-vs-branch_m3 diagnostic view (see
        # SEGMIN_VIEW_* in the PARAMETERS block) - only meaningful for
        # "adtree" (TreeQSM has no seg_min_mm), guarded so a future
        # TreeQSM run doesn't crash or silently do something meaningless.
        if FAMILY == "adtree":
            for tree in REFERENCE_TREES:
                plot_facet_grid(all_rows, FAMILY, branch_filter, SEGMIN_VIEW_METRIC,
                                 SEGMIN_VIEW_FACET_PARAM, SEGMIN_VIEW_X_PARAM,
                                 SEGMIN_VIEW_LINE_PARAM, tree=tree,
                                 filename_tag="segminview")

        # Ranked deviation / correlation heatmap now use the POOLED-across-
        # trees versions here - the single-tree plot_ranked_deviation()/
        # plot_correlation_heatmap() functions above are left untouched and
        # still callable directly (e.g. from a REPL) with one tree.
        plot_ranked_deviation_multi_tree(all_rows, FAMILY, branch_filter, METRIC, REFERENCE_TREES)
        # dbh_m specifically - answers "which adqsm_variant gets closest to
        # the REAL measured DBH" (the "10cm" mode's destructive-reference
        # dbh_m IS that real measurement; "none" mode's DBH reference is
        # only AdQSM's own reported DBH, not independent - generated here
        # for completeness, but "10cm" is the one that actually answers
        # this question).
        plot_ranked_deviation_multi_tree(all_rows, FAMILY, branch_filter, "dbh_m", REFERENCE_TREES)
        for corr_mode in CORR_HEATMAP_MODES:
            plot_correlation_heatmap_multi_tree(all_rows, FAMILY, branch_filter, params_for_corr,
                                                 METRICS_FOR_CORR, REFERENCE_TREES, mode=corr_mode)

        for metric in TORNADO_METRICS:
            plot_tornado_multi_tree(all_rows, FAMILY, branch_filter, metric, params_for_corr, REFERENCE_TREES)

        # Reuses TORNADO_METRICS (not a separate ANOVA_METRICS list) - the
        # same 4 metrics make sense for both charts, and a second near-
        # duplicate list would just be another thing to keep in sync.
        for metric in TORNADO_METRICS:
            plot_anova_multi_tree(all_rows, FAMILY, branch_filter, metric, params_for_corr, REFERENCE_TREES)

    # =====================================================================
    # TreeQSM simp_smallradii/simp_replaceiterations sensitivity - see
    # TREEQSM_STAGE_FOR_SENSITIVITY/TREEQSM_SIMP_PARAMS/TREEQSM_SIMP_METRICS
    # in the PARAMETERS block. Always targets family="treeqsm", independent
    # of the FAMILY-driven loops above - a standalone supplementary
    # analysis, not part of the generic per-FAMILY loop. Run separately for
    # mode="manual", mode="auto", and both pooled together, since mode is a
    # potential confound (same reasoning as the raw-vs-pct_dev correlation
    # heatmap's own rationale for why an unacknowledged confound can
    # distort a pooled reading).
    # =====================================================================
    if RUN_TREEQSM_SIMP_ANALYSIS:
        stage_rows = filter_treeqsm_stage(all_rows, TREEQSM_STAGE_FOR_SENSITIVITY)
        mode_variants = [
            ("simp_pooled", stage_rows),
            ("simp_manual", [r for r in stage_rows if r["mode"] == "manual"]),
            ("simp_auto",   [r for r in stage_rows if r["mode"] == "auto"]),
        ]
        for branch_filter in BRANCH_FILTERS_TO_RUN:
            for tag, variant_rows in mode_variants:
                for metric in TREEQSM_SIMP_METRICS:
                    # reference_rows=all_rows (the FULL, unfiltered dataset)
                    # - NOT variant_rows - so get_reference() can actually
                    # find "Reference (destructive)" ("10cm") or the
                    # resolved AdQSM variant ("none"); neither is itself a
                    # "Simplified (no islands)" TreeQSM row, so it would be
                    # invisible if we only looked inside variant_rows (see
                    # _pct_dev_pooled()'s own docstring for why this
                    # parameter exists). The ANALYZED population stays
                    # variant_rows, unchanged.
                    plot_ranked_deviation_multi_tree(variant_rows, "treeqsm", branch_filter,
                                                      metric, REFERENCE_TREES, filename_tag=tag,
                                                      reference_rows=all_rows)
                    plot_tornado_multi_tree(variant_rows, "treeqsm", branch_filter, metric,
                                             TREEQSM_SIMP_PARAMS, REFERENCE_TREES, filename_tag=tag,
                                             reference_rows=all_rows)
                    # "mode" only added to ANOVA's own param list for the
                    # pooled variant - it would be a single, non-varying
                    # group within simp_manual/simp_auto alone, which
                    # _eta_squared() already handles by returning
                    # None/skipping, but there's no reason to even attempt
                    # it there.
                    anova_params = TREEQSM_SIMP_PARAMS + (["mode"] if tag == "simp_pooled" else [])
                    plot_anova_multi_tree(variant_rows, "treeqsm", branch_filter, metric,
                                           anova_params, REFERENCE_TREES, filename_tag=tag,
                                           reference_rows=all_rows)
                plot_correlation_heatmap_multi_tree(variant_rows, "treeqsm", branch_filter,
                                                     TREEQSM_SIMP_PARAMS, TREEQSM_SIMP_METRICS,
                                                     REFERENCE_TREES, mode="pct_dev", filename_tag=tag,
                                                     reference_rows=all_rows)
