# -*- coding: utf-8 -*-
# =====================================================================
#  Compare reconstruction APPROACHES (AdTree raw/calibrated, TreeQSM
#  manual/auto, TreeQSM de Tanago, other AdQSM variants) against the
#  appropriate reference, in 3 separate blocks - each block pairs ONE
#  TreeQSM stage with the reference that's actually fair to compare it
#  against (see BLOCKS below). This is deliberately NOT a calibration-
#  method decision (see calmethod_decision_summary.py for that) - AdTree
#  here uses ONLY the settled calmethod="regression-perorder"
#  (AdTree_Calmethod_Decision.docx), never "min5mm".
# ---------------------------------------------------------------------
#  WHY 3 separate blocks, not one pooled comparison: TreeQSM's "Optimal",
#  "Filtered <10cm", and "Simplified (no islands)" are different
#  PROCESSING STAGES of the same reconstruction, each meaningful against
#  a DIFFERENT reference:
#    1) Filtered_10cm_vs_Destructive - the one stage that's diameter-
#       restricted to match the destructive field reference's own
#       >=10cm measurement floor - the fair, ground-truth accuracy check.
#    2) Optimal_vs_AdQSM - full/unfiltered TreeQSM vs. AdQSM (a
#       reconstruction, not ground truth) - a methods-vs-each-other cross
#       check, not an accuracy claim (same caveat as compare_volumes.py's
#       own "none" mode).
#    3) Simplified_vs_AdQSM - the stage that actually goes into ANSYS,
#       also vs. AdQSM, restricted to ONE simp_smallradii/
#       simp_replaceiterations combo (SIMP_BASELINE below) so this block
#       compares reconstruction APPROACHES, not simplification settings
#       (that sensitivity is parameter_sensitivity.py's job).
#
#  AdTree is pooled across EVERY radius_threshold_mm/adqsm_variant combo
#  present, same convention calmethod_decision_summary.py already
#  established - the reconstruction APPROACH is under test, not any one
#  particular parameter choice.
#
#  CRITICAL: every group's rows are filtered from `rows_bf` (branch_filter-
#  filtered ONLY), but each tree's REFERENCE row is looked up from that
#  SAME rows_bf, never from an already-group/stage-filtered subset - a
#  stage filter (e.g. "Simplified (no islands)") would silently exclude
#  the reference row too (it's neither AdTree nor that TreeQSM stage),
#  which is exactly the bug this project hit in parameter_sensitivity.py
#  earlier this session. See verify_references_or_die() below, which
#  checks this explicitly before any chart/CSV output is produced.
#
#  Reuses (imports, does not duplicate): compare_volumes.py's
#  load_results/pct_diff/filter_by_branch_filter/REFERENCE_METHOD/
#  RESULTS_CSV/resolve_reference_method_none; calmethod_decision_summary.py's
#  error_stats_from_triples()/tree_spread_std() (the exact same Bias/MAE/
#  RMSE/CV-RMSE math, not reimplemented); parameter_sensitivity.py's
#  REFERENCE_TREES/_rowkey().
#
#  Dependencies: matplotlib (install: pip install matplotlib)
# =====================================================================

import csv
import os
import re

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

from compare_volumes import (
    RESULTS_CSV,
    REFERENCE_METHOD,
    load_results,
    pct_diff,
    filter_by_branch_filter,
    resolve_reference_method_none,
)
from calmethod_decision_summary import error_stats_from_triples, tree_spread_std
from parameter_sensitivity import REFERENCE_TREES, _rowkey
from plot_volumes import (
    FAMILY_GRADIENTS,
    classify_family,
    build_method_color_map,
    shorten_method_label,
    ANNOTATION_FONTSIZE,
)

# =====================  PARAMETERS  ==================================
OUTPUT_CSV = "reconstruction_method_decision_summary.csv"
BEST_PICKS_CSV = "reconstruction_method_best_picks.csv"
PLOTS_DIR = "plots"
METRICS = ["total_m3", "trunk_m3", "branch_m3", "dbh_m", "n_cylinders"]
CHART_METRICS = ["total_m3", "trunk_m3", "branch_m3", "n_cylinders"]   # only these get bar charts, all 5 go in the CSV
SIMP_BASELINE = {"simp_smallradii": 0.005, "simp_replaceiterations": 0.0}   # block 3 only, compared as floats
_SIMP_EPS = 1e-6   # float-equality tolerance for matching SIMP_BASELINE

# "Best single configuration" selection - only meaningful for the two
# groups that actually POOL MULTIPLE candidate rows per tree (AdTree
# calibrated across every radius_threshold_mm/adqsm_variant combo, AdQSM
# across every "other" variant) - the other groups already have exactly
# one row per tree, so there's nothing to select between.
BEST_OF_GROUPS = ["AdTree calibrated (regression-perorder)", "AdQSM (other variants)"]
SINGLE_VALUE_GROUPS = ["AdTree raw", "TreeQSM manual", "TreeQSM auto", "TreeQSM de Tanago"]
COMBINED_SCORE_METRICS = ["total_m3", "trunk_m3", "branch_m3"]

# "Realistic" AdTree pick (Part A) - a NON-cherry-picked configuration,
# achievable on a tree with NO known volume reference (unlike the
# "(best-of-tree)" oracle picks above, which search for the closest
# volume match using knowledge of the answer). radius_threshold_mm and
# seg_* are FIXED, defensible defaults - NOT searched/cherry-picked.
# Matches a commonly-occurring combo already present in the data
# ("seg100-500-k50"), not an arbitrary new choice. adqsm_variant is
# instead chosen via the already-validated DBH-matching method (see
# select_variant_by_dbh()) - achievable on a new tree since it only
# needs a field-measured DBH, not a full destructive volume reference.
REALISTIC_RADIUS_THRESHOLD_MM = 5
REALISTIC_SEG_MIN_MM = 100
REALISTIC_SEG_MAX_MM = 500
REALISTIC_SEG_K_PCT = 50
REALISTIC_GROUPS = ["AdTree calibrated (realistic)"]   # present in ALL 3 blocks (unlike "AdQSM other variants")

BLOCKS = [
    {"name": "Filtered_10cm_vs_Destructive", "branch_filter": "10cm",
     "treeqsm_stage": "Filtered <10cm",
     "de_tanago_method": "TreeQSM de Tanago (mean, Filtered<10cm)",
     "include_adqsm_other": False, "apply_simp_baseline": False},
    {"name": "Optimal_vs_AdQSM", "branch_filter": "none",
     "treeqsm_stage": "Optimal",
     "de_tanago_method": "TreeQSM de Tanago (mean)",
     "include_adqsm_other": True, "apply_simp_baseline": False},
    {"name": "Simplified_vs_AdQSM", "branch_filter": "none",
     "treeqsm_stage": "Simplified (no islands)",
     "de_tanago_method": None,
     "include_adqsm_other": True, "apply_simp_baseline": True},
]
# =====================================================================


def ensure_plots_dir():
    if not os.path.isdir(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)
    return PLOTS_DIR


# ----------------------------------------------------------------------
# Group selectors - each takes rows ALREADY filtered to one branch_filter
# (rows_bf) and returns the matching subset. None of these touch the
# reference row's own lookup - that's handled separately (see
# resolve_ref_methods_for_block()/build_ref_of_tree() below), specifically
# so a group filter can never silently strip the reference out too.
# ----------------------------------------------------------------------
def adtree_raw(rows):
    return [r for r in rows if r["method"].startswith("AdTree raw")]


def adtree_calibrated(rows):
    """ONLY calmethod='regression-perorder' - the settled decision
    (AdTree_Calmethod_Decision.docx). Unlike calmethod_decision_summary.py
    (which deliberately compares min5mm vs. regression-perorder), every
    other script in this project should treat this as the one AdTree
    calibration to use."""
    return [r for r in rows if r["method"].startswith("AdTree calibrated") and r["calmethod"] == "regression-perorder"]


def _matches_simp_baseline(row):
    sr, ri = row.get("simp_smallradii"), row.get("simp_replaceiterations")
    if sr is None or ri is None:
        return False
    return (abs(sr - SIMP_BASELINE["simp_smallradii"]) < _SIMP_EPS
            and abs(ri - SIMP_BASELINE["simp_replaceiterations"]) < _SIMP_EPS)


def treeqsm_manual(rows, stage_substr, apply_simp_baseline=False):
    out = [r for r in rows if "TreeQSM mine" in r["method"] and r["mode"] == "manual" and stage_substr in r["method"]]
    if apply_simp_baseline:
        out = [r for r in out if _matches_simp_baseline(r)]
    return out


def treeqsm_auto(rows, stage_substr, apply_simp_baseline=False):
    out = [r for r in rows if "TreeQSM mine" in r["method"] and r["mode"] == "auto" and stage_substr in r["method"]]
    if apply_simp_baseline:
        out = [r for r in out if _matches_simp_baseline(r)]
    return out


def de_tanago(rows, exact_method_string):
    if exact_method_string is None:
        return []
    return [r for r in rows if r["method"] == exact_method_string]


def adqsm_other_variants(rows, reference_method):
    return [r for r in rows if r["method"].startswith("AdQSM (TreesParams)") and r["method"] != reference_method]


# ----------------------------------------------------------------------
# Reference resolution - ALWAYS from rows_bf (branch_filter-filtered
# only), never from a group/stage-filtered subset. See this module's own
# header comment for why that distinction matters.
# ----------------------------------------------------------------------
def resolve_ref_methods_for_block(rows_bf, branch_filter):
    """{tree: reference_method} for this branch_filter - REFERENCE_METHOD
    (same for every tree) for '10cm'; resolved PER TREE via
    resolve_reference_method_none() for 'none' (variant numbering differs
    per tree - never assume one shared value, same rule as everywhere
    else in this project)."""
    result = {}
    for t in REFERENCE_TREES:
        if branch_filter == "10cm":
            result[t] = REFERENCE_METHOD
        else:
            tree_rows = [r for r in rows_bf if r["tree"] == t]
            result[t] = resolve_reference_method_none(tree_rows)
    return result


def build_ref_of_tree(rows_bf, ref_method_by_tree, metric_key):
    """{tree: reference value for metric_key} - looked up in rows_bf
    (never a group-filtered subset)."""
    ref_of_tree = {}
    for t, ref_method in ref_method_by_tree.items():
        if not ref_method:
            continue
        row = next((r for r in rows_bf if r["tree"] == t and r["method"] == ref_method), None)
        if row is not None and row[metric_key] is not None:
            ref_of_tree[t] = row[metric_key]
    return ref_of_tree


def build_ref_row_of_tree(rows_bf, ref_method_by_tree):
    """{tree: full reference row dict} - looked up in rows_bf (never a
    group-filtered subset), same reference-resolution rule as
    build_ref_of_tree() above (which only extracts ONE metric's value;
    this returns the WHOLE row, needed for combined_pct_score()'s
    simultaneous multi-metric comparison)."""
    ref_row_of_tree = {}
    for t, ref_method in ref_method_by_tree.items():
        if not ref_method:
            continue
        row = next((r for r in rows_bf if r["tree"] == t and r["method"] == ref_method), None)
        if row is not None:
            ref_row_of_tree[t] = row
    return ref_row_of_tree


def combined_pct_score(row, ref_row):
    """Average of |% deviation from ref_row| across total_m3, trunk_m3,
    branch_m3 for one candidate row. Returns None if ANY of the three
    values is missing on either row (ref_row or row) - a row can't be
    scored on a partial set of the three metrics, skip it entirely (the
    row excluded from selection, not crashing) rather than scoring on
    1-2 of the 3."""
    if row is None or ref_row is None:
        return None
    abs_pcts = []
    for metric in COMBINED_SCORE_METRICS:
        key = _rowkey(metric)
        est = row.get(key)
        ref = ref_row.get(key)
        if est is None or ref is None:
            return None
        pct = pct_diff(est, ref)
        if pct is None:
            return None
        abs_pcts.append(abs(pct))
    return sum(abs_pcts) / len(abs_pcts)


def select_best_combined(candidate_rows, ref_row):
    """Among candidate_rows (already restricted to one tree, one group,
    one block), return the single row with the LOWEST
    combined_pct_score(row, ref_row), or None if no row could be scored
    at all (printed why - e.g. every candidate missing branch_m3)."""
    scored = [(combined_pct_score(row, ref_row), row) for row in candidate_rows]
    scored = [(score, row) for score, row in scored if score is not None]
    if not scored:
        print("    No candidate could be scored (missing total_m3/trunk_m3/branch_m3 on "
              "the candidate and/or the reference row) - no winner selected.")
        return None
    scored.sort(key=lambda sr: sr[0])
    return scored[0][1]


def get_real_dbh(rows, tree):
    """Return the destructive reference's dbh_m for `tree` (from the
    branch_filter='10cm' 'Reference (destructive)' row) - the one real,
    physically-measured DBH value, reused across all 3 blocks regardless
    of that block's own branch_filter (a physical property of the tree,
    not something that depends on which TreeQSM stage/reconstruction
    mode is being analyzed)."""
    rows_10cm = filter_by_branch_filter(rows, "10cm")
    row = next((r for r in rows_10cm if r["tree"] == tree and r["method"] == REFERENCE_METHOD), None)
    if row is None or row.get("dbh") is None:
        print("  WARNING: no destructive-reference dbh_m found for tree=%r - cannot resolve a "
              "real DBH for the realistic AdTree pick." % tree)
        return None
    return row["dbh"]


def select_variant_by_dbh(rows, tree, real_dbh_m):
    """Among 'AdQSM (TreesParams) (AdQSM XX)' rows for `tree`
    (branch_filter='none'), return the variant label (e.g. '08') whose
    own dbh_m is closest to `real_dbh_m`. Prints which variant was
    chosen and the resulting |difference| in cm. Returns None (with a
    clear printed reason) if no AdQSM row has a dbh_m value for this
    tree."""
    if real_dbh_m is None:
        print("  Cannot select a variant by DBH for tree=%r - no real DBH available." % tree)
        return None
    rows_none = filter_by_branch_filter(rows, "none")
    candidates = [
        r for r in rows_none
        if r["tree"] == tree and r["method"].startswith("AdQSM (TreesParams)") and r.get("dbh") is not None
    ]
    if not candidates:
        print("  No AdQSM (TreesParams) row with a dbh_m value found for tree=%r - cannot select "
              "a variant by DBH." % tree)
        return None
    best = min(candidates, key=lambda r: abs(r["dbh"] - real_dbh_m))
    variant = _adqsm_variant_of(best)
    diff_cm = abs(best["dbh"] - real_dbh_m) * 100.0
    print("  DBH-match for tree=%s: chosen variant=%s (variant dbh=%.4fm vs real dbh=%.4fm, |diff|=%.2fcm)"
          % (tree, variant, best["dbh"], real_dbh_m, diff_cm))
    return variant


def find_realistic_adtree_row(rows_bf, tree, chosen_variant):
    """The single 'AdTree calibrated' row for `tree` matching ALL of
    REALISTIC_RADIUS_THRESHOLD_MM/chosen_variant/REALISTIC_SEG_*/
    calmethod='regression-perorder' - float comparisons use the same
    tolerance convention (_SIMP_EPS) already used elsewhere in this
    project's structured-param matching (_matches_simp_baseline()).
    Returns None if no such row exists for this tree/block (the caller
    prints exactly which tree/block/combo is missing)."""
    candidates = [
        r for r in rows_bf
        if r["tree"] == tree
        and r["method"].startswith("AdTree calibrated")
        and r["calmethod"] == "regression-perorder"
        and r.get("adqsm_variant") == chosen_variant
        and r.get("radius_threshold_mm") is not None
        and abs(r["radius_threshold_mm"] - REALISTIC_RADIUS_THRESHOLD_MM) < _SIMP_EPS
        and r.get("seg_min_mm") is not None
        and abs(r["seg_min_mm"] - REALISTIC_SEG_MIN_MM) < _SIMP_EPS
        and r.get("seg_max_mm") is not None
        and abs(r["seg_max_mm"] - REALISTIC_SEG_MAX_MM) < _SIMP_EPS
        and r.get("seg_k_pct") is not None
        and abs(r["seg_k_pct"] - REALISTIC_SEG_K_PCT) < _SIMP_EPS
    ]
    if not candidates:
        return None
    if len(candidates) > 1:
        print("  WARNING: %d rows matched the realistic AdTree combo for tree=%r - using the first one."
              % (len(candidates), tree))
    return candidates[0]


def build_groups_for_block(rows_bf, block, ref_method_by_tree):
    """{group_name: group_rows} for this block - insertion order is the
    chart/table display order."""
    groups = {
        "AdTree raw": adtree_raw(rows_bf),
        "AdTree calibrated (regression-perorder)": adtree_calibrated(rows_bf),
        "TreeQSM manual": treeqsm_manual(rows_bf, block["treeqsm_stage"], block["apply_simp_baseline"]),
        "TreeQSM auto": treeqsm_auto(rows_bf, block["treeqsm_stage"], block["apply_simp_baseline"]),
    }
    if block["de_tanago_method"] is not None:
        groups["TreeQSM de Tanago"] = de_tanago(rows_bf, block["de_tanago_method"])
    if block["include_adqsm_other"]:
        # reference_method differs PER TREE in "none" mode - exclude each
        # row's OWN tree's reference method, not one shared value.
        groups["AdQSM (other variants)"] = [
            r for r in rows_bf
            if r["method"].startswith("AdQSM (TreesParams)") and r["method"] != ref_method_by_tree.get(r["tree"])
        ]
    return groups


def gather_group_errors(group_rows, ref_of_tree, metric_key):
    """Same pooling pattern as calmethod_decision_summary.py's
    gather_calmethod_errors(), generalized to any pre-filtered
    `group_rows` (the AdTree-calmethod-specific filtering there is now
    done by the group-selector functions above, before this is called) -
    pools every row in `group_rows` against ref_of_tree[tree] for the
    SAME tree. Returns (triples, per_tree_mean_pct), identical shape to
    that script's own function."""
    triples = []
    pct_by_tree = {}
    for r in group_rows:
        est = r[metric_key]
        ref = ref_of_tree.get(r["tree"])
        if est is None or ref is None:
            continue
        triples.append((r["tree"], est, ref))
        pct = pct_diff(est, ref)
        if pct is not None:
            pct_by_tree.setdefault(r["tree"], []).append(pct)
    per_tree_mean_pct = {t: sum(vals) / len(vals) for t, vals in pct_by_tree.items()}
    return triples, per_tree_mean_pct


def verify_references_or_die(all_rows):
    """MANDATORY bug-prevention check (see this module's own header
    comment for the incident that motivated it) - resolves each block's
    per-tree reference method from rows_bf (branch_filter-filtered only)
    and confirms it's a REAL reference ("Reference (destructive)" for
    "10cm", an "AdQSM (TreesParams) ..." variant for "none"), never None
    and never a silently-wrong fallback. Raises SystemExit immediately if
    any tree/block fails, BEFORE any chart/CSV output is produced."""
    print("=" * 78)
    print("STEP 6 - mandatory reference-resolution check (must pass before continuing)")
    print("-" * 78)
    ok = True
    for block in BLOCKS:
        rows_bf = filter_by_branch_filter(all_rows, block["branch_filter"])
        ref_method_by_tree = resolve_ref_methods_for_block(rows_bf, block["branch_filter"])
        for t in REFERENCE_TREES:
            method = ref_method_by_tree.get(t)
            print("block=%-30s tree=%-12s -> resolved reference: %s" % (block["name"], t, method))
            if not method:
                ok = False
                print("  FATAL: no reference method resolved for this tree/block.")
            elif block["branch_filter"] == "10cm" and method != "Reference (destructive)":
                ok = False
                print("  FATAL: expected 'Reference (destructive)', got %r." % method)
            elif block["branch_filter"] == "none" and not method.startswith("AdQSM (TreesParams)"):
                ok = False
                print("  FATAL: expected an 'AdQSM (TreesParams) ...' variant, got %r." % method)
    print()
    if not ok:
        raise SystemExit("Reference resolution check FAILED (see FATAL lines above) - "
                          "aborting before computing any chart/CSV output.")


def write_summary_csv(summary_rows):
    fieldnames = (
        ["block", "metric", "group", "n_obs", "bias", "mae", "rmse", "cv_rmse_pct",
         "mean_pct_error", "tree_spread_std_pct"]
        + ["%s_pct" % t for t in REFERENCE_TREES]
    )
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)
    print("Wrote:", OUTPUT_CSV)


def write_best_picks_csv(best_picks):
    """best_picks: {(block_name, tree, group_name): (winning_row, score)}
    -> a SEPARATE CSV from OUTPUT_CSV, one row per winning pick, with the
    winning row's OWN values for all 5 METRICS (not just the 3 used to
    select it) so the exact picks are preserved for later reference."""
    fieldnames = ["block", "tree", "group", "winning_method", "combined_pct_score"] + METRICS
    with open(BEST_PICKS_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for (block_name, tree, group_name), (row, score) in best_picks.items():
            out = dict(block=block_name, tree=tree, group=group_name,
                       winning_method=row["method"], combined_pct_score=score)
            for m in METRICS:
                out[m] = row.get(_rowkey(m))
            w.writerow(out)
    print("Wrote:", BEST_PICKS_CSV)


# Y-axis label per metric - same dict-lookup-with-fallback convention as
# plot_volumes.py's plot_tree_overview() FIELD_UNITS handling (m3 for
# volumes, a metric-specific override for the one non-volume field that
# needs its own wording, e.g. "" / a plain count for n_cylinders - here
# "Cylinder count" instead of a nonsensical "Error [m^3]").
_METRIC_YLABELS = {
    "n_cylinders": "Cylinder count",
}


def plot_block_metric_chart(summary_rows, block, metric, group_order):
    by_key = {(r["block"], r["metric"], r["group"]): r for r in summary_rows}
    rows_for_metric = [by_key.get((block["name"], metric, g)) for g in group_order]

    abs_bias = [abs(r["bias"]) if r and r["bias"] is not None else 0.0 for r in rows_for_metric]
    rmse = [r["rmse"] if r and r["rmse"] is not None else 0.0 for r in rows_for_metric]

    x = list(range(len(group_order)))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, 1.8 * len(group_order)), 6))
    ax.bar([xi - width / 2 for xi in x], abs_bias, width=width, label="|Bias|", color="#8fe0cf")
    ax.bar([xi + width / 2 for xi in x], rmse, width=width, label="RMSE", color="#2f8f8a")
    ax.set_xticks(x)
    ax.set_xticklabels(group_order, rotation=20, ha="right")
    ax.set_ylabel(_METRIC_YLABELS.get(metric, "Error [m^3]" if metric.endswith("_m3") else "Error"))
    ax.set_title("%s: %s\n(branch_filter='%s')" % (block["name"], metric, block["branch_filter"]))
    ax.legend()
    fig.tight_layout()

    out_path = os.path.join(ensure_plots_dir(), "reconstruction_decision_%s_%s.png" % (block["name"], metric))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


def _group_family_representative(group_name):
    """Strip a bestof-chart group's DISPLAY name down to a representative
    METHOD-PREFIX string classify_family() can match via startswith() -
    e.g. 'AdTree calibrated (regression-perorder) (best-of-tree)' ->
    'AdTree calibrated', 'TreeQSM manual' -> 'TreeQSM', 'AdQSM (other
    variants) (best-of-tree)' -> 'AdQSM'. Used ONLY to pick a per-GROUP
    tick-label colour (via build_method_color_map()) - never to look up
    an actual row."""
    name = group_name
    if name.endswith(" (best-of-tree)"):
        name = name[: -len(" (best-of-tree)")]
    for prefix in ("AdTree raw", "AdTree calibrated", "TreeQSM", "AdQSM"):
        if name.startswith(prefix):
            return prefix
    return name


def _build_group_tick_color_map(group_order):
    """{group_display_name: colour} - same "build a color map, then
    tick_label.set_color(...)" pattern plot_error_metrics_bar() (in
    plot_volumes.py) already uses, applied to our GROUP names instead of
    individual method strings. reference_method=None (no method in
    `representatives` is ever an actual reference row, so nothing should
    get the isolated "Reference" pink highlight here)."""
    representatives = [_group_family_representative(g) for g in group_order]
    color_of_repr = build_method_color_map(representatives, None)
    return {g: color_of_repr.get(r, "#333333") for g, r in zip(group_order, representatives)}


def _fmt_num_short(v):
    """Same "no spurious trailing .0" cosmetic rule this project's other
    scripts already use (see parameter_sensitivity.py's _fmt_num())."""
    if v is None:
        return "?"
    return str(int(v)) if float(v) == int(v) else str(v)


def _adqsm_variant_of(row):
    """The AdQSM variant number for one row. For AdTree rows,
    row["adqsm_variant"] is already parsed (load_results()) and used
    directly. For "AdQSM (TreesParams) (AdQSM XX)" rows themselves,
    that field is BLANK in the underlying CSV (confirmed directly -
    adqsm_variant is only ever populated on AdTree calibration rows,
    never on AdQSM's own reference rows) - the only place the variant
    number actually exists for those rows is inside the method string
    itself, so this is the one place in this module that parses it from
    there, as a fallback only when the field is empty."""
    variant = row.get("adqsm_variant")
    if variant:
        return variant
    m = re.search(r"AdQSM (\d+)\)$", row.get("method", ""))
    return m.group(1) if m else "?"


def _format_winner_short(base_group, row):
    """Compact winning-configuration text for one (group, tree) pick -
    e.g. 'r10mm/AdQSM08/seg50-500-k50' for either AdTree group (oracle
    best-of-tree OR Part A's realistic pick - same format for both,
    matched by prefix since both start with 'AdTree calibrated'),
    'AdQSM05' for the AdQSM group. Pulled from the row's own already-
    parsed fields (radius_threshold_mm/seg_*, and adqsm_variant for
    AdTree rows) - see _adqsm_variant_of() for the one AdQSM-row
    exception."""
    if base_group.startswith("AdTree calibrated"):
        return "r%smm/AdQSM%s/seg%s-%s-k%s" % (
            _fmt_num_short(row.get("radius_threshold_mm")), _adqsm_variant_of(row),
            _fmt_num_short(row.get("seg_min_mm")), _fmt_num_short(row.get("seg_max_mm")),
            _fmt_num_short(row.get("seg_k_pct")))
    if base_group == "AdQSM (other variants)":
        return "AdQSM%s" % _adqsm_variant_of(row)
    return row.get("method", "")


def _best_picks_lookup_group(group_display_name):
    """Map a chart/table DISPLAY group name back to the group key stored
    in `best_picks` - "X (best-of-tree)" -> "X" (the oracle groups,
    whose display name adds that suffix), "AdTree calibrated
    (realistic)" -> itself unchanged (Part A's realistic group's display
    name already IS its best_picks key - no suffix to strip, it was
    never an oracle pick). Returns None for a SINGLE_VALUE_GROUPS member
    (no per-tree selection, not in best_picks at all)."""
    if group_display_name.endswith(" (best-of-tree)"):
        return group_display_name[: -len(" (best-of-tree)")]
    if group_display_name in REALISTIC_GROUPS:
        return group_display_name
    return None


def plot_bestof_pooled_chart(summary_rows, block, metric, bestof_group_order, best_picks):
    """Same grouped-bar style as plot_block_metric_chart(), but reading
    the "(best-of-tree)" stats rows (Step 5) for the two groups that HAVE
    one, alongside the as-is single-value groups for comparison - a
    SEPARATE chart, does not touch/replace plot_block_metric_chart()'s
    own pooled-average chart. Each group's x-axis tick label is tinted
    with its family colour, and a text box (only when a BEST_OF_GROUPS
    group is actually present) lists all 3 trees' winning configurations,
    since one pooled bar can't show 3 different per-tree winners."""
    by_key = {(r["block"], r["metric"], r["group"]): r for r in summary_rows}
    rows_for_metric = [by_key.get((block["name"], metric, g)) for g in bestof_group_order]

    abs_bias = [abs(r["bias"]) if r and r["bias"] is not None else 0.0 for r in rows_for_metric]
    rmse = [r["rmse"] if r and r["rmse"] is not None else 0.0 for r in rows_for_metric]

    x = list(range(len(bestof_group_order)))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, 1.8 * len(bestof_group_order)), 6))
    ax.bar([xi - width / 2 for xi in x], abs_bias, width=width, label="|Bias|", color="#8fe0cf")
    ax.bar([xi + width / 2 for xi in x], rmse, width=width, label="RMSE", color="#2f8f8a")
    ax.set_xticks(x)
    ax.set_xticklabels(bestof_group_order, rotation=20, ha="right")
    color_of_group = _build_group_tick_color_map(bestof_group_order)
    for tick_label, g in zip(ax.get_xticklabels(), bestof_group_order):
        tick_label.set_color(color_of_group.get(g, "#333333"))
    ax.set_ylabel(_METRIC_YLABELS.get(metric, "Error [m^3]" if metric.endswith("_m3") else "Error"))
    ax.set_title("%s: %s (best-of-tree selection)\n(branch_filter='%s')" % (block["name"], metric, block["branch_filter"]))
    ax.legend(loc="upper left")

    # Winner-list text box - upper RIGHT (legend is pinned upper LEFT
    # above), only built when at least one group with a real per-tree
    # selection (oracle best-of-tree OR Part A's realistic pick) is
    # actually present in this block (skipped entirely otherwise).
    lines = []
    for g in bestof_group_order:
        lookup_group = _best_picks_lookup_group(g)
        if lookup_group is None:
            continue
        lines.append("%s:" % g)
        for t in REFERENCE_TREES:
            entry = best_picks.get((block["name"], t, lookup_group))
            if entry is None:
                lines.append("  %s: (no pick)" % t)
                continue
            lines.append("  %s: %s" % (t, _format_winner_short(lookup_group, entry[0])))
    if lines:
        ax.text(0.98, 0.98, "\n".join(lines), transform=ax.transAxes, ha="right", va="top",
                 fontsize=ANNOTATION_FONTSIZE, family="monospace",
                 bbox=dict(boxstyle="round", facecolor="white", edgecolor="#999999", alpha=0.9))

    fig.tight_layout()

    out_path = os.path.join(ensure_plots_dir(), "reconstruction_decision_%s__bestof_pooled_%s.png" % (block["name"], metric))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


# Tree shade order (lightest -> darkest), CONSISTENT across every group
# in plot_bestof_pertree_chart() - IND01_054 = lightest, IND03_088 =
# middle, IND07_083 = darkest. Same t = i/(n-1)-across-a-gradient formula
# build_method_color_map() already uses elsewhere in this project.
_TREE_SHADE_T = {t: (i / (len(REFERENCE_TREES) - 1) if len(REFERENCE_TREES) > 1 else 0.5)
                 for i, t in enumerate(REFERENCE_TREES)}

# Generic grey swatches for the shade-legend (Step B4) - NOT tied to any
# one family's own gradient (that would wrongly imply the legend only
# applies to one group), just illustrating "lighter = which tree".
_SHADE_LEGEND_GREYS = ["#cccccc", "#888888", "#333333"]


def _group_tree_color(group_name, tree):
    """One bar's colour in plot_bestof_pertree_chart(): the GROUP's
    family gradient (FAMILY_GRADIENTS, via the same representative-
    prefix classification _build_group_tick_color_map() already uses),
    at the SHADE position for this specific tree (_TREE_SHADE_T) - e.g.
    "AdTree calibrated ..." groups all render in shades of blue, with
    the shade itself telling trees apart."""
    family = _group_family_representative(group_name)
    stops = FAMILY_GRADIENTS.get(family, ["#9a9a9a", "#9a9a9a"])
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "%s_gradient" % (family or "unknown").lower().replace(" ", "_"), stops)
    return cmap(_TREE_SHADE_T.get(tree, 0.5))


def _build_shade_legend_handles():
    """Proxy legend handles (Step B4) explaining the lightest/middle/
    darkest -> tree convention - a SEPARATE legend from the (removed)
    per-tree colour legend, since colour now encodes family+shade, not
    tree alone."""
    tree_shade_words = []
    n = len(REFERENCE_TREES)
    for i in range(n):
        if i == 0:
            tree_shade_words.append("Lightest")
        elif i == n - 1:
            tree_shade_words.append("Darkest")
        else:
            tree_shade_words.append("Middle")
    return [
        mpatches.Patch(facecolor=_SHADE_LEGEND_GREYS[min(i, len(_SHADE_LEGEND_GREYS) - 1)],
                       edgecolor="#666666", label="%s: %s" % (word, t))
        for i, (t, word) in enumerate(zip(REFERENCE_TREES, tree_shade_words))
    ]


def plot_bestof_pertree_chart(groups, best_picks, ref_of_tree, block, metric, bestof_group_order):
    """THREE bars per group (one per reference tree) instead of one
    aggregated bar - each tree's own % deviation from ITS OWN reference,
    for that tree's single best-pick row (oracle best-of-tree, Part A's
    realistic pick) or single as-is row (SINGLE_VALUE_GROUPS) - so per-
    tree variation is visible rather than collapsed into one pooled
    number. x-axis tick labels are plain black (unlike the pooled
    chart, which keeps its own family-coloured tick labels); each bar
    is coloured by its GROUP's family, shaded by WHICH TREE it is
    (_group_tree_color()); a winner-list text box (same style as
    plot_bestof_pooled_chart()'s own) replaces the old per-bar rotated
    annotations."""
    metric_key = _rowkey(metric)

    def row_for(group_name, tree):
        lookup_group = _best_picks_lookup_group(group_name)
        if lookup_group is not None:
            entry = best_picks.get((block["name"], tree, lookup_group))
            return (entry[0] if entry else None), lookup_group
        tree_rows = [r for r in groups.get(group_name, []) if r["tree"] == tree]
        return (tree_rows[0] if tree_rows else None), None

    def value_for(group_name, tree):
        row, _lookup_group = row_for(group_name, tree)
        if row is None or row.get(metric_key) is None or ref_of_tree.get(tree) is None:
            return None
        return pct_diff(row[metric_key], ref_of_tree[tree])

    x = list(range(len(bestof_group_order)))
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(6, 2.6 * len(bestof_group_order)), 6.5))
    for i, tree in enumerate(REFERENCE_TREES):
        offset = (i - (len(REFERENCE_TREES) - 1) / 2.0) * width
        values = [value_for(g, tree) or 0.0 for g in bestof_group_order]
        xs = [xi + offset for xi in x]
        colors = [_group_tree_color(g, tree) for g in bestof_group_order]
        ax.bar(xs, values, width=width, color=colors)
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)
    ax.set_xticks(x)
    ax.set_xticklabels(bestof_group_order, rotation=20, ha="right")   # plain black - Step B1

    # Winner-list text box - same style/positioning as
    # plot_bestof_pooled_chart()'s own box (upper right), listing every
    # group with a real per-tree selection.
    lines = []
    for g in bestof_group_order:
        lookup_group = _best_picks_lookup_group(g)
        if lookup_group is None:
            continue
        lines.append("%s:" % g)
        for t in REFERENCE_TREES:
            entry = best_picks.get((block["name"], t, lookup_group))
            if entry is None:
                lines.append("  %s: (no pick)" % t)
                continue
            lines.append("  %s: %s" % (t, _format_winner_short(lookup_group, entry[0])))
    if lines:
        ax.text(0.98, 0.98, "\n".join(lines), transform=ax.transAxes, ha="right", va="top",
                 fontsize=ANNOTATION_FONTSIZE, family="monospace",
                 bbox=dict(boxstyle="round", facecolor="white", edgecolor="#999999", alpha=0.9))

    ax.set_ylabel("%% deviation from reference (%s)" % metric)
    ax.set_title("%s: %s per-tree best-pick deviation\n(branch_filter='%s')" % (block["name"], metric, block["branch_filter"]))
    # Shade-convention legend (Step B4) - upper LEFT, so it never
    # overlaps the winner-list box (upper right, see above).
    ax.legend(handles=_build_shade_legend_handles(), title="Bar shade = Tree", loc="upper left")
    fig.tight_layout()

    out_path = os.path.join(ensure_plots_dir(), "reconstruction_decision_%s__bestof_pertree_%s.png" % (block["name"], metric))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


def append_best_pick_stats(summary_rows, block, rows_bf, ref_method_by_tree, group_display_name, winning_rows):
    """Bias/MAE/RMSE/tree-spread stats for a set of per-tree winning
    rows (exactly one per tree, already selected elsewhere) - appended
    to `summary_rows` under `group_display_name`, reusing
    gather_group_errors()/error_stats_from_triples()/tree_spread_std()
    exactly as the pooled per-group loop does (Step 5's own math, not
    reimplemented). Shared by both the oracle best-of-tree groups and
    Part A's realistic-pick group so this stats-computation shape is
    written once, not duplicated three times."""
    if not winning_rows:
        return
    for metric in METRICS:
        metric_key = _rowkey(metric)
        ref_of_tree = build_ref_of_tree(rows_bf, ref_method_by_tree, metric_key)
        triples, per_tree_mean_pct = gather_group_errors(winning_rows, ref_of_tree, metric_key)
        stats = error_stats_from_triples(triples)
        spread = tree_spread_std(per_tree_mean_pct, REFERENCE_TREES)
        row = dict(
            block=block["name"], metric=metric, group=group_display_name,
            n_obs=stats["n_obs"], bias=stats["bias"], mae=stats["mae"],
            rmse=stats["rmse"], cv_rmse_pct=stats["cv_rmse"],
            mean_pct_error=(sum(per_tree_mean_pct[t] for t in REFERENCE_TREES if t in per_tree_mean_pct)
                            / len([t for t in REFERENCE_TREES if t in per_tree_mean_pct])
                            if any(t in per_tree_mean_pct for t in REFERENCE_TREES) else None),
            tree_spread_std_pct=spread,
        )
        for t in REFERENCE_TREES:
            row["%s_pct" % t] = per_tree_mean_pct.get(t)
        summary_rows.append(row)


def print_block_summaries(summary_rows):
    print("=" * 78)
    print("PER-BLOCK SUMMARY (total_m3 RMSE comparison)")
    print("-" * 78)
    for block in BLOCKS:
        rows_this = [r for r in summary_rows
                     if r["block"] == block["name"] and r["metric"] == "total_m3" and r["rmse"] is not None]
        if not rows_this:
            print("%s: no usable total_m3 data." % block["name"])
            continue
        best = min(rows_this, key=lambda r: r["rmse"])
        print("%s: lowest total_m3 RMSE = '%s' (RMSE=%.4f, bias=%.4f, n=%d)"
              % (block["name"], best["group"], best["rmse"], best["bias"], best["n_obs"]))
        interesting = {r["group"]: r for r in rows_this
                       if r["group"] in ("AdTree calibrated (regression-perorder)", "TreeQSM manual", "TreeQSM auto")}
        if interesting:
            ranked = sorted(interesting.items(), key=lambda kv: kv[1]["rmse"])
            print("   AdTree calibrated / TreeQSM manual / TreeQSM auto, ranked by RMSE: "
                  + ", ".join("%s (RMSE=%.4f)" % (k, v["rmse"]) for k, v in ranked))
    print()


# =========================  RUN  =====================================
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit("'%s' not found - run compare_volumes.py first so it gets created." % RESULTS_CSV)

    all_rows = load_results(RESULTS_CSV)

    verify_references_or_die(all_rows)

    summary_rows = []
    group_order_by_block = {}
    # Captured per-block (purely additive - the pooled loop's own logic
    # below is completely unchanged) so the best-of-tree block further
    # down can reuse rows_bf/ref_method_by_tree/groups for ALL 3 blocks
    # without recomputing them.
    rows_bf_by_block = {}
    ref_method_by_tree_by_block = {}
    groups_by_block = {}
    for block in BLOCKS:
        rows_bf = filter_by_branch_filter(all_rows, block["branch_filter"])
        ref_method_by_tree = resolve_ref_methods_for_block(rows_bf, block["branch_filter"])
        groups = build_groups_for_block(rows_bf, block, ref_method_by_tree)
        group_order_by_block[block["name"]] = list(groups.keys())
        rows_bf_by_block[block["name"]] = rows_bf
        ref_method_by_tree_by_block[block["name"]] = ref_method_by_tree
        groups_by_block[block["name"]] = groups

        for metric in METRICS:
            metric_key = _rowkey(metric)
            ref_of_tree = build_ref_of_tree(rows_bf, ref_method_by_tree, metric_key)
            for group_name, group_rows in groups.items():
                triples, per_tree_mean_pct = gather_group_errors(group_rows, ref_of_tree, metric_key)
                if not triples:
                    print("NOTE: no data for block=%r metric=%r group=%r." % (block["name"], metric, group_name))
                stats = error_stats_from_triples(triples)
                spread = tree_spread_std(per_tree_mean_pct, REFERENCE_TREES)
                row = dict(
                    block=block["name"], metric=metric, group=group_name,
                    n_obs=stats["n_obs"], bias=stats["bias"], mae=stats["mae"],
                    rmse=stats["rmse"], cv_rmse_pct=stats["cv_rmse"],
                    mean_pct_error=(sum(per_tree_mean_pct[t] for t in REFERENCE_TREES if t in per_tree_mean_pct)
                                    / len([t for t in REFERENCE_TREES if t in per_tree_mean_pct])
                                    if any(t in per_tree_mean_pct for t in REFERENCE_TREES) else None),
                    tree_spread_std_pct=spread,
                )
                for t in REFERENCE_TREES:
                    row["%s_pct" % t] = per_tree_mean_pct.get(t)
                summary_rows.append(row)

    # =====================================================================
    # Best single configuration ("best-of-tree") selection - additive,
    # separate from the pooled-group computation above. Only meaningful
    # for BEST_OF_GROUPS (AdTree calibrated, AdQSM other variants) - the
    # other groups already have exactly one row per tree.
    # =====================================================================
    print("=" * 78)
    print("BEST-OF-TREE SELECTION (single combined total/trunk/branch score per tree)")
    print("-" * 78)
    best_picks = {}   # (block_name, tree, group_name) -> (winning_row, score)
    for block in BLOCKS:
        rows_bf = rows_bf_by_block[block["name"]]
        ref_method_by_tree = ref_method_by_tree_by_block[block["name"]]
        groups = groups_by_block[block["name"]]
        ref_row_of_tree = build_ref_row_of_tree(rows_bf, ref_method_by_tree)

        for group_name in BEST_OF_GROUPS:
            if group_name not in groups:
                continue
            group_rows = groups[group_name]
            for t in REFERENCE_TREES:
                candidates = [r for r in group_rows if r["tree"] == t]
                ref_row = ref_row_of_tree.get(t)
                if ref_row is None:
                    print("  [%s] %s / %s -> no reference row available, skipping selection."
                          % (block["name"], t, group_name))
                    continue
                best_row = select_best_combined(candidates, ref_row)
                if best_row is None:
                    print("  [%s] %s / %s -> NO WINNER (%d candidate(s), none scoreable)"
                          % (block["name"], t, group_name, len(candidates)))
                    continue
                score = combined_pct_score(best_row, ref_row)
                print("[%s] %s / %s -> BEST: %s" % (block["name"], t, group_name, best_row["method"]))
                best_picks[(block["name"], t, group_name)] = (best_row, score)
    print()

    # =====================================================================
    # Part A - "realistic" AdTree pick: radius_threshold_mm/seg_* FIXED
    # (never searched), adqsm_variant chosen via DBH-matching against the
    # one real, physically-measured DBH - achievable on a NEW tree with
    # no known volume reference, unlike the oracle best-of-tree pick
    # above. chosen_variant is resolved ONCE PER TREE (the real DBH is a
    # physical property of the tree, independent of branch_filter mode)
    # and then applied identically across all 3 BLOCKS.
    # =====================================================================
    print("=" * 78)
    print("REALISTIC ADTREE PICK (fixed radius_threshold_mm/seg_*, variant chosen by DBH match)")
    print("-" * 78)
    chosen_variant_by_tree = {}
    for t in REFERENCE_TREES:
        real_dbh = get_real_dbh(all_rows, t)
        chosen_variant_by_tree[t] = select_variant_by_dbh(all_rows, t, real_dbh)
    print()

    for block in BLOCKS:
        rows_bf = rows_bf_by_block[block["name"]]
        ref_method_by_tree = ref_method_by_tree_by_block[block["name"]]
        ref_row_of_tree = build_ref_row_of_tree(rows_bf, ref_method_by_tree)
        for t in REFERENCE_TREES:
            chosen_variant = chosen_variant_by_tree[t]
            if chosen_variant is None:
                print("  [%s] %s / AdTree calibrated (realistic) -> skipped, no DBH-matched variant."
                      % (block["name"], t))
                continue
            row = find_realistic_adtree_row(rows_bf, t, chosen_variant)
            if row is None:
                print("  [%s] %s / AdTree calibrated (realistic) -> MISSING: no row matches "
                      "radius_threshold_mm=%s, adqsm_variant=%s, seg=%s-%s-k%s, calmethod=regression-perorder."
                      % (block["name"], t, REALISTIC_RADIUS_THRESHOLD_MM, chosen_variant,
                         REALISTIC_SEG_MIN_MM, REALISTIC_SEG_MAX_MM, REALISTIC_SEG_K_PCT))
                continue
            ref_row = ref_row_of_tree.get(t)
            score = combined_pct_score(row, ref_row)
            print("[%s] %s / AdTree calibrated (realistic) -> %s" % (block["name"], t, row["method"]))
            best_picks[(block["name"], t, "AdTree calibrated (realistic)")] = (row, score)
    print()

    write_best_picks_csv(best_picks)
    print()

    # Step 5: Bias/MAE/RMSE/tree-spread stats for the per-tree best-picks
    # (both the oracle best-of-tree groups AND Part A's realistic-pick
    # group) - added as ADDITIONAL rows to summary_rows before the ONE
    # write_summary_csv() call below, so both pooled and best-of/realistic
    # results land in the same OUTPUT_CSV, clearly distinguished.
    for block in BLOCKS:
        rows_bf = rows_bf_by_block[block["name"]]
        ref_method_by_tree = ref_method_by_tree_by_block[block["name"]]
        groups = groups_by_block[block["name"]]
        for group_name in BEST_OF_GROUPS:
            if group_name not in groups:
                continue
            winning_rows = [best_picks[(block["name"], t, group_name)][0]
                             for t in REFERENCE_TREES if (block["name"], t, group_name) in best_picks]
            append_best_pick_stats(summary_rows, block, rows_bf, ref_method_by_tree,
                                    group_name + " (best-of-tree)", winning_rows)
        for group_name in REALISTIC_GROUPS:
            winning_rows = [best_picks[(block["name"], t, group_name)][0]
                             for t in REFERENCE_TREES if (block["name"], t, group_name) in best_picks]
            append_best_pick_stats(summary_rows, block, rows_bf, ref_method_by_tree,
                                    group_name, winning_rows)

    write_summary_csv(summary_rows)
    print()

    for block in BLOCKS:
        for metric in CHART_METRICS:
            plot_block_metric_chart(summary_rows, block, metric, group_order_by_block[block["name"]])
    print()

    # Best-of-tree / realistic-pick charts (Step 6 / Part A Step A6) -
    # additive, alongside (not replacing) the pooled-average charts just
    # saved above.
    for block in BLOCKS:
        groups = groups_by_block[block["name"]]
        bestof_group_order = (
            [g for g in SINGLE_VALUE_GROUPS if g in groups]
            + [g + " (best-of-tree)" for g in BEST_OF_GROUPS if g in groups]
            + list(REALISTIC_GROUPS)
        )
        for metric in CHART_METRICS:
            plot_bestof_pooled_chart(summary_rows, block, metric, bestof_group_order, best_picks)
            ref_of_tree = build_ref_of_tree(rows_bf_by_block[block["name"]],
                                             ref_method_by_tree_by_block[block["name"]], _rowkey(metric))
            plot_bestof_pertree_chart(groups, best_picks, ref_of_tree, block, metric, bestof_group_order)
    print()

    print_block_summaries(summary_rows)
