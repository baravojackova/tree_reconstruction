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

from compare_volumes import RESULTS_CSV, REFERENCE_METHOD, load_results, resolve_reference_method_none
from plot_volumes import FAMILY_GRADIENTS, classify_family, ensure_plots_dir, shorten_method_label
from plot_box import parse_treeqsm_method, TREEQSM_REF_LINE_COLOR

# =====================  PARAMETERS  ===================================
SELECT_TREE = "IND01_054"
BRANCH_FILTER = "10cm"     # "10cm" or "none" - NEVER both in one run/chart;
                             # this drives which reference row is used (see get_reference())
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


def _output_dir():
    """plots/<SELECT_TREE>/ - same per-tree subfolder convention
    adtree_reconstruct_compare.py's own FIGURES_DIR already established,
    built on top of ensure_plots_dir() (plots/) rather than duplicating it."""
    out_dir = os.path.join(ensure_plots_dir(), SELECT_TREE)
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def _save(fig, family, branch_filter, metric, chart):
    out_path = os.path.join(
        _output_dir(), "sensitivity_%s_%s_%s_%s.png" % (family, branch_filter, metric, chart))
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


# ----------------------------------------------------------------------
# 1) Facet grid: one panel per facet_param value, metric vs. x_param,
#    lines coloured by line_param.
# ----------------------------------------------------------------------
def plot_facet_grid(rows, family, branch_filter, metric, facet_param, x_param, line_param):
    sel = select_rows(rows, family, SELECT_TREE, branch_filter)
    if not sel:
        print("No rows for family=%r, tree=%r, branch_filter=%r - skipping facet grid."
              % (family, SELECT_TREE, branch_filter))
        return

    tree_rows = [r for r in rows if r["tree"] == SELECT_TREE and r["branch_filter"] == branch_filter]
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
        % (family, metric, SELECT_TREE, branch_filter, ref_method), fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 0.90))
    # Mode suffix on THIS chart's filename only (ranked/heatmap have no
    # such toggle) - so absolute/relative runs of the same
    # family/branch_filter/metric don't silently overwrite each other.
    _save(fig, family, branch_filter, metric, "facetgrid_%s" % FACET_Y_MODE)


# ----------------------------------------------------------------------
# 2) Ranked deviation dot plot: every individual row, sorted by |% dev|,
#    dashed line at 0% - also a practical "which combo is closest" leaderboard.
# ----------------------------------------------------------------------
def plot_ranked_deviation(rows, family, branch_filter, metric):
    sel = select_rows(rows, family, SELECT_TREE, branch_filter)
    if not sel:
        print("No rows for family=%r, tree=%r, branch_filter=%r - skipping ranked deviation plot."
              % (family, SELECT_TREE, branch_filter))
        return

    tree_rows = [r for r in rows if r["tree"] == SELECT_TREE and r["branch_filter"] == branch_filter]
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
                 % (family, metric, SELECT_TREE, branch_filter))
    ax.legend(fontsize=8, loc="best")
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    _save(fig, family, branch_filter, metric, "ranked")

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


def plot_correlation_heatmap(rows, family, branch_filter, params, metrics):
    sel = select_rows(rows, family, SELECT_TREE, branch_filter)
    if not sel:
        print("No rows for family=%r, tree=%r, branch_filter=%r - skipping correlation heatmap."
              % (family, SELECT_TREE, branch_filter))
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
        % (family, SELECT_TREE, branch_filter), fontsize=9)
    fig.tight_layout()
    _save(fig, family, branch_filter, "allmetrics", "corrheatmap")


# =========================  RUN  =====================================
# Calls all three plotting functions once, for the currently-set FAMILY/
# BRANCH_FILTER/METRIC - does NOT loop over every combination automatically
# in this first version (re-run by hand with different PARAMETER values;
# a batch-all-combinations mode can be a follow-up once this version's
# been reviewed).
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

    plot_facet_grid(all_rows, FAMILY, BRANCH_FILTER, METRIC, facet_param, x_param, line_param)
    plot_ranked_deviation(all_rows, FAMILY, BRANCH_FILTER, METRIC)
    plot_correlation_heatmap(all_rows, FAMILY, BRANCH_FILTER, params_for_corr, METRICS_FOR_CORR)
