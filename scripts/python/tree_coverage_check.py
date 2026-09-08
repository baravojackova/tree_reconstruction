# -*- coding: utf-8 -*-
# =====================================================================
#  Pure READ / print-only diagnostic: for tree X, which parameter
#  combinations already have rows in volume_results.csv, for which
#  branch_filter modes, and is a reference row present at all?
#
#  WHY: Bara is about to reconstruct 2 more reference trees (B21_S01,
#  Ave_B) and wants to quickly check each new tree's coverage against
#  IND01_054's, to keep the parameter grid comparable across trees before
#  running the planned multi-tree sensitivity comparison. This script
#  never writes anything - no backup/staged-confirmation needed beyond the
#  normal "review before running" pass.
#
#  Reused, not reimplemented: load_results()/RESULTS_CSV (compare_volumes.py),
#  select_rows()/get_reference() (parameter_sensitivity.py),
#  parse_treeqsm_method() (plot_box.py).
# ---------------------------------------------------------------------
#  Dependencies: none beyond the project's other scripts (no matplotlib/numpy
#  needed - this is plain-text console output only).
# =====================================================================

import csv
import itertools
import os
import re

from compare_volumes import RESULTS_CSV, load_results
from parameter_sensitivity import select_rows, get_reference
from plot_box import parse_treeqsm_method
from paths import ensure_csv_dir

# =====================  PARAMETERS  ===================================
SELECT_TREE = "ALL"   # a specific tree name, or "ALL" to report on every tree present in the CSV
# =====================================================================


def _fmt_val(v):
    """Format one value for the CSV export's "values" cell - integer-
    looking floats without a trailing ".0" (e.g. 8.0 -> "8"), everything
    else via plain str()."""
    if isinstance(v, float) and v == int(v):
        return str(int(v))
    return str(v)


def _csv_row(tree, family, branch_filter, mode, parameter, values):
    """One long-format CSV row: multiple values joined with '; ' in a
    single cell (not one row per value) - keeps the CSV short and each
    parameter easy to scan on one line in Excel."""
    return {"tree": tree, "family": family, "branch_filter": branch_filter or "",
            "mode": mode or "", "parameter": parameter,
            "values": "; ".join(_fmt_val(v) for v in values)}


def _cartesian_gap_report(label, axis_names, axis_values, actual_combos):
    """Print total-possible vs. actual combination counts for the cartesian
    product of `axis_values` (one sorted-unique-value list per axis), and
    list the SPECIFIC missing combinations (not just a count) when fewer
    than the full product actually appears in `actual_combos` (a set of
    tuples, same shape/order as axis_values).

    Skips the check entirely (no data to cross, not a "gap") when any axis
    has zero values - that means this bucket is empty, already reported by
    the caller's own row-count line above."""
    if any(len(vals) == 0 for vals in axis_values):
        print("    (skipping coverage-gap check - at least one axis has no values)")
        return

    total_possible = 1
    for vals in axis_values:
        total_possible *= len(vals)

    print("    %s: %d of %d possible combinations present"
          % (label, len(actual_combos), total_possible))
    if len(actual_combos) < total_possible:
        all_combos = set(itertools.product(*axis_values))
        missing = sorted(all_combos - actual_combos, key=lambda c: [str(x) for x in c])
        print("    MISSING (%d):" % len(missing))
        for combo in missing:
            print("      %s" % ", ".join("%s=%s" % (name, val) for name, val in zip(axis_names, combo)))
    else:
        print("    (no gaps)")


def report_reference_availability(rows, tree):
    print("\n[1] Reference availability")
    csv_rows = []
    for bf in ("10cm", "none"):
        tree_rows = [r for r in rows if r["tree"] == tree and r["branch_filter"] == bf]
        ref_row, ref_method, ref_color = get_reference(tree_rows, bf)
        if ref_row is not None:
            print("  branch_filter=%-5s -> %s" % (bf, ref_method))
            value = ref_method
        else:
            detail = ("resolved method %r has no matching row" % ref_method) if ref_method \
                else "no candidate method resolved at all"
            print("  branch_filter=%-5s -> NONE (%s)" % (bf, detail))
            value = "NONE"
        csv_rows.append(_csv_row(tree, "reference", bf, "", "resolved_method", [value]))
    return csv_rows


def report_adtree_calibrated(rows, tree, branch_filter):
    sel = select_rows(rows, "adtree", tree, branch_filter)
    print("  branch_filter=%s: %d row(s)" % (branch_filter, len(sel)))
    if not sel:
        return []

    variants = sorted({r["adqsm_variant"] for r in sel if r["adqsm_variant"] not in (None, "")})
    thresholds = sorted({r["radius_threshold_mm"] for r in sel if r["radius_threshold_mm"] is not None})
    segs = sorted({(r["seg_min_mm"], r["seg_max_mm"], r["seg_k_pct"]) for r in sel
                   if r["seg_min_mm"] is not None and r["seg_max_mm"] is not None and r["seg_k_pct"] is not None})
    print("    adqsm_variant values: %s" % variants)
    print("    radius_threshold_mm values: %s" % thresholds)
    print("    (seg_min_mm, seg_max_mm, seg_k_pct) combos: %s" % segs)

    actual = {(r["adqsm_variant"], r["radius_threshold_mm"],
               (r["seg_min_mm"], r["seg_max_mm"], r["seg_k_pct"])) for r in sel}
    _cartesian_gap_report("AdTree calibrated coverage",
                           ["adqsm_variant", "radius_threshold_mm", "seg(min,max,k)"],
                           [variants, thresholds, segs], actual)

    # CSV export: same underlying value-lists as above (variants/thresholds
    # already are; seg_mins/seg_maxs/seg_ks are the marginal per-field
    # projections of the SAME `segs` combos already printed, not a fresh
    # scan of `sel`) - re-expressed in adtree_reconstruct_compare.py's own
    # PARAMETERS-block units (meters/ratio, not the CSV's mm/percent
    # storage) so a value can be pasted straight into RADIUS_THRESHOLDS/
    # SEG_LEN_MIN/SEG_LEN_MAX/SEG_LEN_K without manual conversion.
    seg_mins = sorted({s[0] for s in segs})
    seg_maxs = sorted({s[1] for s in segs})
    seg_ks = sorted({s[2] for s in segs})
    return [
        _csv_row(tree, "adtree_calibrated", branch_filter, "", "ADQSM_VARIANTS", variants),
        _csv_row(tree, "adtree_calibrated", branch_filter, "", "RADIUS_THRESHOLDS_m",
                  [t / 1000.0 for t in thresholds]),
        _csv_row(tree, "adtree_calibrated", branch_filter, "", "SEG_LEN_MIN_m", [s / 1000.0 for s in seg_mins]),
        _csv_row(tree, "adtree_calibrated", branch_filter, "", "SEG_LEN_MAX_m", [s / 1000.0 for s in seg_maxs]),
        _csv_row(tree, "adtree_calibrated", branch_filter, "", "SEG_LEN_K", [k / 100.0 for k in seg_ks]),
    ]


def report_adtree_raw(rows, tree, branch_filter):
    # select_rows() only supports "adtree" (AdTree calibrated) and "treeqsm" -
    # it deliberately EXCLUDES "AdTree raw" rows from the "adtree" family
    # (see its own docstring: raw rows don't vary by the AdQSM-dependent
    # parameters that family analyzes). Filtered locally here instead -
    # same base tree/branch_filter pattern select_rows() itself uses.
    sel = [r for r in rows if r["tree"] == tree and r["branch_filter"] == branch_filter
           and r["method"].startswith("AdTree raw")]
    print("  branch_filter=%s: %d row(s)" % (branch_filter, len(sel)))
    if not sel:
        return []

    thresholds = sorted({r["radius_threshold_mm"] for r in sel if r["radius_threshold_mm"] is not None})
    segs = sorted({(r["seg_min_mm"], r["seg_max_mm"], r["seg_k_pct"]) for r in sel
                   if r["seg_min_mm"] is not None and r["seg_max_mm"] is not None and r["seg_k_pct"] is not None})
    print("    radius_threshold_mm values: %s" % thresholds)
    print("    (seg_min_mm, seg_max_mm, seg_k_pct) combos: %s" % segs)

    actual = {(r["radius_threshold_mm"], (r["seg_min_mm"], r["seg_max_mm"], r["seg_k_pct"])) for r in sel}
    _cartesian_gap_report("AdTree raw coverage",
                           ["radius_threshold_mm", "seg(min,max,k)"],
                           [thresholds, segs], actual)

    # Same conversion-to-PARAMETERS-block-units convention as
    # report_adtree_calibrated() above - no ADQSM_VARIANTS row here (raw
    # rows correctly have none, see STEP 5's fix).
    seg_mins = sorted({s[0] for s in segs})
    seg_maxs = sorted({s[1] for s in segs})
    seg_ks = sorted({s[2] for s in segs})
    return [
        _csv_row(tree, "adtree_raw", branch_filter, "", "RADIUS_THRESHOLDS_m", [t / 1000.0 for t in thresholds]),
        _csv_row(tree, "adtree_raw", branch_filter, "", "SEG_LEN_MIN_m", [s / 1000.0 for s in seg_mins]),
        _csv_row(tree, "adtree_raw", branch_filter, "", "SEG_LEN_MAX_m", [s / 1000.0 for s in seg_maxs]),
        _csv_row(tree, "adtree_raw", branch_filter, "", "SEG_LEN_K", [k / 100.0 for k in seg_ks]),
    ]


def report_treeqsm(rows, tree, branch_filter):
    sel = select_rows(rows, "treeqsm", tree, branch_filter)
    print("  branch_filter=%s: %d row(s)" % (branch_filter, len(sel)))
    if not sel:
        return []

    modes = sorted({r["mode"] for r in sel if r["mode"]})
    pd_combos = sorted({(r["pd1"], r["pd2min"], r["pd2max"]) for r in sel
                         if r["pd1"] is not None and r["pd2min"] is not None and r["pd2max"] is not None})
    simp_combos = sorted({(r["simp_maxorder"], r["simp_smallradii"], r["simp_replaceiterations"]) for r in sel
                           if r["simp_maxorder"] is not None and r["simp_smallradii"] is not None
                           and r["simp_replaceiterations"] is not None})
    stages = sorted({parse_treeqsm_method(r["method"])[1] for r in sel
                      if parse_treeqsm_method(r["method"]) is not None})
    print("    mode values: %s" % modes)
    print("    (pd1, pd2min, pd2max) combos: %s" % pd_combos)
    print("    (simp_maxorder, simp_smallradii, simp_replaceiterations) combos: %s" % simp_combos)
    print("    stage values (reported separately, NOT part of the gap check - "
          "different stages of the same run are expected, not a gap): %s" % stages)

    # mode and (pd1, pd2min, pd2max) are COUPLED, not independent axes -
    # choosing a mode determines where the PD values come from (auto always
    # produces its own searched PD, manual always uses the fixed
    # man_PD1/2Min/2Max) - so treating them as two separate cartesian axes
    # produces false-positive "gaps" for pairings that can never
    # legitimately exist (e.g. auto mode with manual's fixed PD values -
    # confirmed as a real false positive on IND01_054's own data before
    # this fix). Combined into ONE recon_key axis instead, built from the
    # PAIRINGS actually observed in the data (not modes x pd_combos).
    # simp_maxorder/simp_smallradii/simp_replaceiterations show no such
    # coupling with mode/PD in the current data (vary independently in
    # every observed row), so they stay a separate axis for now - revisit
    # if future tree data reveals a similar coupling there.
    recon_keys = sorted({(r["mode"], r["pd1"], r["pd2min"], r["pd2max"]) for r in sel
                          if r["mode"] and r["pd1"] is not None and r["pd2min"] is not None
                          and r["pd2max"] is not None})
    print("    recon (mode, pd1, pd2min, pd2max) combos actually paired together: %s" % recon_keys)

    actual = {((r["mode"], r["pd1"], r["pd2min"], r["pd2max"]),
               (r["simp_maxorder"], r["simp_smallradii"], r["simp_replaceiterations"])) for r in sel}
    _cartesian_gap_report("TreeQSM coverage (recon x simp)",
                           ["recon(mode,pd1,pd2min,pd2max)", "simp(maxorder,smallradii,replaceiter)"],
                           [recon_keys, simp_combos], actual)

    # CSV export: per-mode marginal unique values for each of the 6
    # structured TreeQSM params - one CSV row per (mode, parameter),
    # derived by filtering the SAME `sel` rows by mode (not a fresh
    # CSV re-read), matching the coupling insight above (mode/PD move
    # together, so values are naturally grouped by mode here too).
    csv_rows = []
    for m in modes:
        mode_rows = [r for r in sel if r["mode"] == m]
        field_labels = [
            ("pd1", "pd1"), ("pd2min", "pd2min"), ("pd2max", "pd2max"),
            ("simp_maxorder", "simp_maxorder"), ("simp_smallradii", "simp_smallradii"),
            ("simp_replaceiterations", "simp_replaceiterations"),
        ]
        for field_key, param_name in field_labels:
            vals = sorted({r[field_key] for r in mode_rows if r[field_key] is not None})
            csv_rows.append(_csv_row(tree, "treeqsm", branch_filter, m, param_name, vals))
    return csv_rows


_ADQSM_DIRECT_RE = re.compile(r"^AdQSM \(TreesParams\) \(AdQSM (\d+)\)$")


def report_adqsm_direct(rows, tree):
    # AdQSM-direct rows are deliberately excluded from select_rows()'s
    # adtree/treeqsm families (see parameter_sensitivity.py's own
    # select_rows() docstring: they have no structured AdTree OR TreeQSM
    # columns to plot/analyze) - matched here by method-string prefix
    # instead, same regex resolve_reference_method_none() uses internally.
    sel = [r for r in rows if r["tree"] == tree and r["method"].startswith("AdQSM (TreesParams)")]
    variants = sorted({m.group(1) for m in (_ADQSM_DIRECT_RE.match(r["method"]) for r in sel) if m})
    print("\n[5] AdQSM-direct coverage")
    print("  adqsm_variant values with a row: %s" % variants)
    return [_csv_row(tree, "adqsm_direct", "", "", "adqsm_variant", variants)]


def print_tree_report(rows, tree):
    print("=" * 100)
    print("TREE: %s" % tree)
    print("=" * 100)

    csv_rows = []
    csv_rows.extend(report_reference_availability(rows, tree))

    print("\n[2] AdTree calibrated coverage")
    for bf in ("10cm", "none"):
        csv_rows.extend(report_adtree_calibrated(rows, tree, bf))

    print("\n[3] AdTree raw coverage")
    for bf in ("10cm", "none"):
        csv_rows.extend(report_adtree_raw(rows, tree, bf))

    print("\n[4] TreeQSM coverage")
    for bf in ("10cm", "none"):
        csv_rows.extend(report_treeqsm(rows, tree, bf))

    csv_rows.extend(report_adqsm_direct(rows, tree))
    print()
    return csv_rows


def write_coverage_csv(csv_rows, out_path):
    """Write the same information the console printout shows, in long
    format (one row per parameter, multiple values joined with '; ' in one
    cell) - a lightweight, disposable diagnostic export, not part of the
    shared volume_results.csv conventions, but still routed under csv/ like
    every other CSV this project writes (see ensure_csv_dir() at the call
    site below) - out_path is already a full path by the time it gets here."""
    fieldnames = ["tree", "family", "branch_filter", "mode", "parameter", "values"]
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)
    print("Saved:", out_path)


# =========================  RUN  =====================================
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit("'%s' not found - run compare_volumes.py first so it gets created." % RESULTS_CSV)

    all_rows = load_results(RESULTS_CSV)
    trees = sorted({r["tree"] for r in all_rows}) if SELECT_TREE.upper() == "ALL" else [SELECT_TREE]

    all_csv_rows = []
    for i, tree in enumerate(trees):
        if i > 0:
            print("\n" + "#" * 100 + "\n")
        all_csv_rows.extend(print_tree_report(all_rows, tree))

    out_name = "coverage_all_trees.csv" if SELECT_TREE.upper() == "ALL" else "coverage_%s.csv" % SELECT_TREE
    out_path = os.path.join(ensure_csv_dir(), out_name)
    write_coverage_csv(all_csv_rows, out_path)
