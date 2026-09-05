# -*- coding: utf-8 -*-
# =====================================================================
#  Decide between AdTree's two radius-calibration methods -
#  calmethod="min5mm" vs. calmethod="regression-perorder" - as the
#  PRIMARY calibration this project uses going forward, using an
#  OBJECTIVE, DECOMPOSED comparison (total/trunk/branch volume computed
#  and reported SEPARATELY, never just "total" alone).
# ---------------------------------------------------------------------
#  WHY decomposed, not just total: an earlier investigation (see
#  CHANGELOG_adtree.md) found that closeness in TOTAL volume alone can
#  hide ERROR CANCELLATION between trunk over-correction and branch
#  under-reconstruction - two calibration methods can land on a similar
#  total_m3 number for completely different (and not equally good)
#  reasons. Reporting trunk/branch separately (and checking their bias
#  SIGNS against each other - see the cancellation check below) is the
#  only way to catch that.
#
#  WHICH branch_filter MODE IS "PRIMARY":
#    - "10cm" (vs. the REAL destructive field reference) is the PRIMARY
#      evidence here - it's the only mode with an actual ground-truth
#      measurement behind it (see compare_volumes.py's own header for
#      why branch_filter exists at all).
#    - "none" (vs. AdQSM) is SECONDARY / confirmatory ONLY. Both AdTree
#      calibration methods are themselves DERIVED FROM AdQSM data (that's
#      what "calibration" means here - rescaling AdTree radii to match
#      AdQSM's own per-branch-order values), so comparing them against
#      AdQSM in "none" mode is partially circular: it can confirm a
#      total-consistency story, but it can never be treated as an
#      independent accuracy check the way "10cm" can.
#
#  This script pools EVERY radius_threshold_mm/adqsm_variant combination
#  present for a given calmethod (rather than picking just one combo) -
#  the calibration METHOD is the thing under test here, not any one
#  particular radius threshold or AdQSM variant choice.
#
#  Produces: OUTPUT_CSV (one row per field x branch_filter x calmethod),
#  three PNG charts (10cm mode only, one per field), and a plain-language
#  printed recommendation at the end.
#
#  Dependencies: matplotlib (install: pip install matplotlib). Reuses
#  load_results/pct_diff/filter_by_branch_filter/REFERENCE_METHOD/
#  RESULTS_CSV/resolve_reference_method_none from compare_volumes.py -
#  see that file for what those actually do, not re-explained here.
# =====================================================================

import csv
import math
import os

import matplotlib.pyplot as plt

from compare_volumes import (
    RESULTS_CSV,
    REFERENCE_METHOD,
    load_results,
    pct_diff,
    filter_by_branch_filter,
    resolve_reference_method_none,
)

# =====================  PARAMETERS  ==================================
OUTPUT_CSV = "calmethod_decision_summary.csv"
PLOTS_DIR = "plots"   # same convention as plot_volumes.py

FIELDS = ["total", "trunk", "branch"]
CALMETHODS = ["min5mm", "regression-perorder"]

# The 3 reference trees with a real destructive reference row - see
# compare_volumes.py's REFERENCE_METHOD. Hard-coded as individual CSV
# columns further down (one "<tree>_pct" column each) rather than a
# generic loop, since a fixed 3-tree comparison is what this decision is
# actually based on right now. If a 4th reference tree is EVER added
# later, this list AND the per-tree CSV columns in build_summary_rows()
# below would both need updating by hand - not attempted generically here.
REFERENCE_TREES = ["IND01_054", "IND03_088", "IND07_083"]

BRANCH_FILTERS = ["10cm", "none"]

# Cancellation check (Step 5) threshold - a bias below this (in m^3) is
# treated as "negligible", regardless of its sign, so two tiny opposite-
# sign biases don't get flagged as "cancellation" when neither is large
# enough to matter.
CANCELLATION_MIN_ABS_BIAS_M3 = 0.05
# =====================================================================


def gather_calmethod_errors(rows, reference_method, field, calmethod):
    """For rows already filtered to one branch_filter mode: pool EVERY
    'AdTree calibrated ...' row whose calmethod matches `calmethod`
    (across all radius_threshold_mm/adqsm_variant combinations present),
    matched against `reference_method`'s own row for the SAME tree.

    Returns:
      - triples: a flat list of (tree, est, ref) - the raw pooled
        observations, one per matched combination x tree.
      - per_tree_mean_pct: {tree: mean % error across that tree's
        combinations} - for the cross-tree consistency check. A tree
        with zero matching combinations for this calmethod is simply
        absent from this dict (not inserted as a fake 0).

    Only rows where method.startswith("AdTree calibrated") AND
    row["calmethod"] == calmethod are included - raw/uncalibrated AdTree
    rows and every other method (TreeQSM, AdQSM, the reference itself)
    are excluded by construction.
    """
    ref_of_tree = {r["tree"]: r[field] for r in rows if r["method"] == reference_method}

    triples = []
    pct_by_tree = {}   # tree -> list of % errors, one per matching combination
    for r in rows:
        if not r["method"].startswith("AdTree calibrated"):
            continue
        if r["calmethod"] != calmethod:
            continue
        est = r[field]
        ref = ref_of_tree.get(r["tree"])
        if est is None or ref is None:
            continue
        triples.append((r["tree"], est, ref))
        pct = pct_diff(est, ref)
        if pct is not None:
            pct_by_tree.setdefault(r["tree"], []).append(pct)

    per_tree_mean_pct = {t: sum(vals) / len(vals) for t, vals in pct_by_tree.items()}
    return triples, per_tree_mean_pct


def error_stats_from_triples(triples):
    """Bias/MAE/RMSE/CV-RMSE% from a flat (tree, est, ref) triple list -
    the EXACT SAME formulas as compare_volumes.py's compute_error_metrics()
    (bias = mean(est-ref), mae = mean(|est-ref|), rmse = sqrt(mean((est-ref)^2)),
    cv_rmse = rmse / mean(ref) * 100) - not reinvented, just applied to a
    pooled-across-combinations triple list instead of compute_error_metrics()'s
    own one-row-per-(tree,method) shape (that function can't be called
    directly here: it keys lookups by (tree, method) in a dict, which would
    silently DROP all but one combination per tree - exactly the pooling
    this script needs to keep)."""
    n = len(triples)
    if n == 0:
        return dict(n_obs=0, bias=None, mae=None, rmse=None, cv_rmse=None)
    errs = [est - ref for _, est, ref in triples]
    refs = [ref for _, _, ref in triples]
    bias = sum(errs) / n
    mae = sum(abs(e) for e in errs) / n
    rmse = math.sqrt(sum(e * e for e in errs) / n)
    cv_rmse = rmse / (sum(refs) / n) * 100.0 if refs else 0.0
    return dict(n_obs=n, bias=bias, mae=mae, rmse=rmse, cv_rmse=cv_rmse)


def tree_spread_std(per_tree_mean_pct, trees):
    """Population standard deviation (not sample - dividing by n, not
    n-1) of the per-tree mean %% errors across `trees` - a small value
    means this calmethod behaves CONSISTENTLY across the 3 reference
    trees, a large value means it's sensitive to which tree you look at.
    Trees missing from per_tree_mean_pct (zero matching combinations)
    are simply excluded, same as gather_calmethod_errors() leaves them out."""
    vals = [per_tree_mean_pct[t] for t in trees if t in per_tree_mean_pct]
    if not vals:
        return None
    mean = sum(vals) / len(vals)
    return math.sqrt(sum((v - mean) ** 2 for v in vals) / len(vals))


def ensure_plots_dir():
    if not os.path.isdir(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)
    return PLOTS_DIR


def build_summary_rows(rows, reference_method_10cm, reference_method_none):
    """One dict per (field, branch_filter, calmethod) - see the module
    header / STEP 4 in the task this script was written for. Also returns
    a parallel {(field, branch_filter, calmethod): triples} dict so the
    Step 6 charts and Step 7 recommendation can reuse the SAME pooled
    observations without recomputing them."""
    rows_by_filter = {bf: filter_by_branch_filter(rows, bf) for bf in BRANCH_FILTERS}
    reference_method_of = {"10cm": reference_method_10cm, "none": reference_method_none}

    summary_rows = []
    triples_by_key = {}
    for field in FIELDS:
        for branch_filter in BRANCH_FILTERS:
            reference_method = reference_method_of[branch_filter]
            filtered = rows_by_filter[branch_filter]
            for calmethod in CALMETHODS:
                triples, per_tree_mean_pct = gather_calmethod_errors(
                    filtered, reference_method, field, calmethod)
                triples_by_key[(field, branch_filter, calmethod)] = triples
                stats = error_stats_from_triples(triples)
                spread = tree_spread_std(per_tree_mean_pct, REFERENCE_TREES)

                row = dict(
                    field=field,
                    branch_filter=branch_filter,
                    calmethod=calmethod,
                    n_obs=stats["n_obs"],
                    bias_m3=stats["bias"],
                    mae_m3=stats["mae"],
                    rmse_m3=stats["rmse"],
                    cv_rmse_pct=stats["cv_rmse"],
                    mean_pct_error=(sum(per_tree_mean_pct[t] for t in REFERENCE_TREES if t in per_tree_mean_pct)
                                    / len([t for t in REFERENCE_TREES if t in per_tree_mean_pct])
                                    if any(t in per_tree_mean_pct for t in REFERENCE_TREES) else None),
                    tree_spread_std_pct=spread,
                )
                for t in REFERENCE_TREES:
                    row["%s_pct" % t] = per_tree_mean_pct.get(t)
                summary_rows.append(row)
    return summary_rows, triples_by_key


def write_summary_csv(summary_rows):
    fieldnames = (
        ["field", "branch_filter", "calmethod", "n_obs", "bias_m3", "mae_m3", "rmse_m3",
         "cv_rmse_pct", "mean_pct_error", "tree_spread_std_pct"]
        + ["%s_pct" % t for t in REFERENCE_TREES]
    )
    # tree_spread_std_pct is the POPULATION standard deviation (divide by n,
    # not n-1) of the 3 reference trees' own mean %% errors - noted here in
    # the header comment (not a separate column) since that's what the
    # number in that column actually is.
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in summary_rows:
            w.writerow(row)
    print("Wrote:", OUTPUT_CSV)


def print_cancellation_check(summary_rows):
    """STEP 5: trunk vs. branch bias side by side for branch_filter='10cm',
    per calmethod - re-prints the Step-4 rows in this specific layout,
    does NOT recompute anything."""
    by_key = {(r["field"], r["branch_filter"], r["calmethod"]): r for r in summary_rows}

    print("=" * 78)
    print("CANCELLATION CHECK (branch_filter='10cm' - trunk vs. branch bias sign)")
    print("-" * 78)
    for calmethod in CALMETHODS:
        trunk_row = by_key.get(("trunk", "10cm", calmethod))
        branch_row = by_key.get(("branch", "10cm", calmethod))
        trunk_bias = trunk_row["bias_m3"] if trunk_row else None
        branch_bias = branch_row["bias_m3"] if branch_row else None
        print("calmethod=%-22s trunk_bias=%8s m3   branch_bias=%8s m3" % (
            calmethod,
            "%.4f" % trunk_bias if trunk_bias is not None else "-",
            "%.4f" % branch_bias if branch_bias is not None else "-"))
        if trunk_bias is None or branch_bias is None:
            print("  (missing data - cannot evaluate)")
            continue
        opposite_signs = (trunk_bias > 0) != (branch_bias > 0)
        both_meaningful = abs(trunk_bias) > CANCELLATION_MIN_ABS_BIAS_M3 and abs(branch_bias) > CANCELLATION_MIN_ABS_BIAS_M3
        if opposite_signs and both_meaningful:
            print("  WARNING: possible error cancellation")
        else:
            print("  OK: no sign of trunk/branch error cancellation for this calmethod")
    print()


def plot_calmethod_decision_chart(summary_rows, field):
    """STEP 6: one grouped bar chart for `field`, branch_filter='10cm' only
    (the primary evidence - see module header). 2 x-axis groups (one per
    calmethod), each showing |Bias| and RMSE side by side - same 2-bars-
    per-group visual style as plot_volumes.py's plot_error_metrics_bar(),
    just with CALMETHODS instead of full method strings as the x-axis
    categories."""
    by_key = {(r["field"], r["branch_filter"], r["calmethod"]): r for r in summary_rows}
    rows_for_field = [by_key.get((field, "10cm", cm)) for cm in CALMETHODS]

    abs_bias = [abs(r["bias_m3"]) if r and r["bias_m3"] is not None else 0.0 for r in rows_for_field]
    rmse = [r["rmse_m3"] if r and r["rmse_m3"] is not None else 0.0 for r in rows_for_field]

    x = list(range(len(CALMETHODS)))
    width = 0.35

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.bar([xi - width / 2 for xi in x], abs_bias, width=width, label="|Bias|", color="#8fe0cf")
    ax.bar([xi + width / 2 for xi in x], rmse, width=width, label="RMSE", color="#2f8f8a")
    ax.set_xticks(x)
    ax.set_xticklabels(CALMETHODS)
    ax.set_ylabel("Volume error [m^3]")
    ax.set_title("%s: min5mm vs. regression-perorder\n(branch_filter='10cm')"
                  % {"total": "Total volume", "trunk": "Trunk volume", "branch": "Branch volume"}[field])
    ax.legend()
    fig.tight_layout()

    out_path = os.path.join(ensure_plots_dir(), "calmethod_decision_%s.png" % field)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


def print_recommendation(summary_rows):
    """STEP 7: plain-language final recommendation, citing the actual
    numbers from `summary_rows` (never asserted without a number next to
    it)."""
    by_key = {(r["field"], r["branch_filter"], r["calmethod"]): r for r in summary_rows}

    branch_min5mm = by_key.get(("branch", "10cm", "min5mm"))
    branch_reg = by_key.get(("branch", "10cm", "regression-perorder"))
    trunk_min5mm = by_key.get(("trunk", "10cm", "min5mm"))
    trunk_reg = by_key.get(("trunk", "10cm", "regression-perorder"))

    print("=" * 78)
    print("RECOMMENDATION")
    print("-" * 78)

    # a) branch RMSE comparison
    if branch_min5mm and branch_reg and branch_min5mm["rmse_m3"] is not None and branch_reg["rmse_m3"] is not None:
        rmse_min5mm, rmse_reg = branch_min5mm["rmse_m3"], branch_reg["rmse_m3"]
        if rmse_reg < rmse_min5mm:
            lower_name, lower_rmse, higher_rmse = "regression-perorder", rmse_reg, rmse_min5mm
        else:
            lower_name, lower_rmse, higher_rmse = "min5mm", rmse_min5mm, rmse_reg
        diff = higher_rmse - lower_rmse
        print("a) Branch volume (10cm mode): '%s' has lower RMSE (%.4f m3 vs. %.4f m3,"
              " a difference of %.4f m3)." % (lower_name, lower_rmse, higher_rmse, diff))
    else:
        print("a) Branch volume (10cm mode): insufficient data to compare RMSE.")

    # b) trunk comparison
    if trunk_min5mm and trunk_reg and trunk_min5mm["bias_m3"] is not None and trunk_reg["bias_m3"] is not None:
        bias_min5mm, bias_reg = trunk_min5mm["bias_m3"], trunk_reg["bias_m3"]
        rmse_min5mm_t, rmse_reg_t = trunk_min5mm["rmse_m3"], trunk_reg["rmse_m3"]
        bias_diff = abs(bias_min5mm - bias_reg)
        rmse_diff = abs(rmse_min5mm_t - rmse_reg_t) if rmse_min5mm_t is not None and rmse_reg_t is not None else None
        print("b) Trunk volume (10cm mode): bias min5mm=%.4f m3 vs. regression-perorder=%.4f m3"
              " (difference %.4f m3); RMSE %.4f m3 vs. %.4f m3%s." % (
                  bias_min5mm, bias_reg, bias_diff, rmse_min5mm_t, rmse_reg_t,
                  " (difference %.4f m3)" % rmse_diff if rmse_diff is not None else ""))
        if bias_diff < 0.05:
            print("   -> As expected, trunk numbers do NOT differ meaningfully between the two"
                  " calmethods (trunk calibration is shared).")
        else:
            print("   -> UNEXPECTED: trunk numbers differ by more than a negligible amount between"
                  " the two calmethods - this contradicts the assumption that trunk calibration"
                  " is shared and should be investigated before trusting this comparison.")
    else:
        print("b) Trunk volume (10cm mode): insufficient data to compare.")

    # c) tree consistency (spread)
    if branch_min5mm and branch_reg and branch_min5mm["tree_spread_std_pct"] is not None \
            and branch_reg["tree_spread_std_pct"] is not None:
        spread_min5mm = branch_min5mm["tree_spread_std_pct"]
        spread_reg = branch_reg["tree_spread_std_pct"]
        more_consistent = "regression-perorder" if spread_reg < spread_min5mm else "min5mm"
        print("c) Branch volume (10cm mode) cross-tree consistency: tree_spread_std_pct"
              " min5mm=%.2f%% vs. regression-perorder=%.2f%% - '%s' is more consistent across"
              " the %d reference trees." % (spread_min5mm, spread_reg, more_consistent, len(REFERENCE_TREES)))
    else:
        print("c) Branch volume (10cm mode) cross-tree consistency: insufficient data to compare.")

    # d) final call
    print("-" * 78)
    if branch_min5mm and branch_reg and branch_min5mm["rmse_m3"] is not None and branch_reg["rmse_m3"] is not None \
            and branch_min5mm["tree_spread_std_pct"] is not None and branch_reg["tree_spread_std_pct"] is not None:
        rmse_winner = "regression-perorder" if branch_reg["rmse_m3"] < branch_min5mm["rmse_m3"] else "min5mm"
        consistency_winner = "regression-perorder" if branch_reg["tree_spread_std_pct"] < branch_min5mm["tree_spread_std_pct"] else "min5mm"
        if rmse_winner == consistency_winner:
            print("d) RECOMMENDATION: use calmethod='%s' - it has both the lower branch RMSE"
                  " and the more consistent cross-tree behavior, with no meaningful trunk-side"
                  " difference to offset that advantage." % rmse_winner)
        else:
            print("d) RECOMMENDATION: the evidence is MIXED - '%s' has the lower branch RMSE but"
                  " '%s' is more consistent across the 3 reference trees. No single calmethod wins"
                  " on both criteria here; a decision would need to weigh accuracy vs. consistency"
                  " explicitly rather than being read off this table alone." % (rmse_winner, consistency_winner))
    else:
        print("d) RECOMMENDATION: insufficient data to make a call.")
    print()


# =========================  RUN  =====================================
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit("'%s' not found - run compare_volumes.py first so it gets created." % RESULTS_CSV)

    rows = load_results(RESULTS_CSV)

    # REFERENCE_METHOD_NONE resolved dynamically, same pattern as
    # compare_volumes.py's own RUN section (and plot_volumes.py's, since
    # a previous session) - printed once, before the field/branch_filter
    # loops below use it.
    rows_none = filter_by_branch_filter(rows, "none")
    reference_method_none = resolve_reference_method_none(rows_none)
    if reference_method_none is None:
        raise SystemExit("No 'AdQSM (TreesParams) (AdQSM XX)' row found in branch_filter='none' data - "
                          "cannot resolve a 'none'-mode reference.")
    print("Resolved REFERENCE_METHOD_NONE for this run: %s" % reference_method_none)
    print()

    summary_rows, _triples_by_key = build_summary_rows(rows, REFERENCE_METHOD, reference_method_none)
    write_summary_csv(summary_rows)
    print()

    print_cancellation_check(summary_rows)

    for field in FIELDS:
        plot_calmethod_decision_chart(summary_rows, field)
    print()

    print_recommendation(summary_rows)
