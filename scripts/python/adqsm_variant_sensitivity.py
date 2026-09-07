# -*- coding: utf-8 -*-
# =====================================================================
#  Reproduces (and saves) the ad hoc AdQSM-variant-sensitivity analysis:
#  for a dense grid of AdQSM segmentation-setting variants on one tree,
#  overlay AdQSM's own total volume against AdTree's CALIBRATED total
#  volume (per variant) on one log-scale chart, so the "close cluster
#  near the reference DBH vs. wildly inflated outliers" split (already
#  found once, ad hoc, for B21_S01: 44 close / 17 outliers out of 61
#  variants) can be reproduced for ANY tree with the same data available,
#  and reviewed later from the saved CSV/PNG instead of re-deriving it by
#  hand each time.
#
#  WHY AdQSM and AdTree-calibrated total volume track each other so
#  closely: AdTree's calibrated trunk radius is taken ENTIRELY from
#  AdQSM's own taper.txt curve at each height (see
#  adtree_reconstruct_compare.py's trunk_radius_func comment) - so a bad/
#  inflated AdQSM taper curve for a given variant propagates directly
#  into AdTree's calibrated result for that SAME variant, not just
#  AdQSM's own row. This script's whole point is to make that shared
#  failure mode visible per variant, for both series at once.
#
#  Data sources (both already produced by other scripts in this repo -
#  this script does not reconstruct anything itself):
#    - RESULTS_CSV (volume_results.csv, written by adtree_reconstruct_
#      compare.py) - "AdQSM (TreesParams) (AdQSM NNN)" rows and "AdTree
#      calibrated ... (AdQSM NNN)" rows (both calmethod=regression-
#      perorder and, optionally, calref=min5mm).
#    - TAPER_SUMMARY_CSV (taper_curve_compare_summary.csv, written by
#      taper_curve_compare.py) - dbh_cm/pct_diff_vs_reference per variant,
#      computed from AdQSM's raw taper.txt curve against the measured
#      reference DBH.
#
#  Both files must already exist for TREE_NAME (run adtree_reconstruct_
#  compare.py and taper_curve_compare.py for it first if not).
#
#  Dependencies: numpy, matplotlib (install: pip install numpy matplotlib)
# =====================================================================

import csv
import os

import numpy as np
import matplotlib.pyplot as plt

# Reuse (do not re-implement): load_results() already parses
# volume_results.csv into row dicts with every column this script needs
# (total/trunk/method/calmethod/adqsm_variant), including the structured
# AdTree columns (calmethod, adqsm_variant) - see compare_volumes.py's own
# header comment for the full column list.
from compare_volumes import load_results

# Reuse (do not re-implement): the one place in this codebase that
# already knows how to pull an AdQSM variant number out of EITHER an
# AdQSM (TreesParams) row (parsed from the method string, since
# adqsm_variant is blank on those rows) OR an AdTree calibrated row
# (already-structured adqsm_variant column) - see its own docstring in
# reconstruction_method_decision_summary.py for why both cases are
# needed. Imported despite the leading underscore - it is a plain
# module-level function, not a class attribute, so it is perfectly
# importable; it is just named as an internal helper within that file.
from reconstruction_method_decision_summary import _adqsm_variant_of

# Reuse (do not re-implement): the plots/ folder convention every other
# per-tree chart in this codebase already builds on (plots/<tree>/...).
from plot_volumes import ensure_plots_dir

# Shared visual style (colors/sizes) - see plot_style.py's own header.
# This chart previously had no colors/sizes in common with the rest of
# the project at all (default matplotlib color cycle, markersize=3,
# linewidth=1, default fontsizes) - PANEL_TITLE_FONTSIZE/AXIS_LABEL_FONTSIZE
# below are used UNCHANGED (this chart's title/axis labels match every
# other chart's sizing); CHART_LINEWIDTH/CHART_MARKERSIZE/
# CHART_LEGEND_FONTSIZE (STYLE block below) are deliberately BIGGER than
# the project defaults, per Bara's request for this chart specifically.
from plot_style import FAMILY_GRADIENTS, PANEL_TITLE_FONTSIZE, AXIS_LABEL_FONTSIZE

# =====================  PARAMETERS  ===================================
TREE_NAME = "B21_S01"

# True (default): also read/plot/export the secondary "[calref=min5mm]"
# AdTree-calibrated series alongside the primary "[calmethod=regression-
# perorder]" one. False: omit it entirely (chart, table, and console
# summary all skip it) - e.g. when COMPUTE_CALREF_MIN5MM was False for
# this tree's adtree_reconstruct_compare.py run and no calref=min5mm rows
# exist in RESULTS_CSV for it at all.
INCLUDE_CALREF_MIN5MM = True

# Same path convention already used by every other script in this repo
# for these two files (adtree_reconstruct_compare.py/compare_volumes.py/
# etc. for RESULTS_CSV; taper_curve_compare.py for TAPER_SUMMARY_CSV) -
# both are plain relative paths, resolved from wherever this script is
# actually run from (normally scripts/python/).
RESULTS_CSV = "volume_results.csv"
TAPER_SUMMARY_CSV = "taper_curve_compare_summary.csv"

# The "close to reference" vs. "outlier" split threshold used by the
# console summary (Step 6) - a variant's |pct_diff_vs_reference| (already
# computed by taper_curve_compare.py against the measured reference DBH)
# below this is "close", at/above it is an "outlier". 10 (%) matches the
# cutoff named in this task; expose it as a parameter rather than a
# hard-coded literal so it can be tuned later without hunting through the
# function bodies below.
PCT_DIFF_OUTLIER_THRESHOLD_PCT = 10.0
# =====================================================================

# =====================  STYLE  ========================================
# This chart's own emphasis tier - heavier than plot_style.py's project-
# wide DEFAULT_LINEWIDTH(1)/DEFAULT_MARKERSIZE(6, though this chart never
# used that one either)/LEGEND_FONTSIZE(8), per Bara's explicit request
# for THIS chart specifically (many overlapping thin lines were hard to
# read). PANEL_TITLE_FONTSIZE/AXIS_LABEL_FONTSIZE are used directly from
# plot_style.py, unmodified, in plot_variant_sensitivity() below - no
# local override needed for those two.
CHART_LINEWIDTH = 5        # thicker than the project default (1) - Bara's request
CHART_MARKERSIZE = 10       # larger than the project default (this chart's own old value, 3)
CHART_LEGEND_FONTSIZE = 12  # larger than the project default (8)
AXIS_LABEL_FONTSIZE = 12
ANNOTATION_FONTSIZE = 12
TITLE_FONTSIZE = 14
# =====================================================================


def read_taper_summary(path, tree_name):
    """Read TAPER_SUMMARY_CSV (taper_curve_compare.py's own output - no
    reusable reader function exists there to import, it only WRITES this
    file in its own RUN section, so a local implementation is needed
    here) into {variant: {"dbh_cm": float, "pct_diff_vs_reference": float}}
    for `tree_name` only. Raises SystemExit if the file is missing
    entirely (nothing to merge against) - a missing PER-VARIANT row
    inside an existing file is handled separately, gracefully, by the
    caller (see missing_taper below), not here."""
    if not os.path.exists(path):
        raise SystemExit(
            "'%s' not found - run taper_curve_compare.py for %r first "
            "(with it included in TREES_TO_RUN)." % (path, tree_name))
    with open(path, "r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    by_variant = {}
    for r in rows:
        if r["tree"] != tree_name:
            continue
        by_variant[r["variant"]] = {
            "dbh_cm": float(r["dbh_cm"]) if r.get("dbh_cm") not in (None, "") else None,
            "pct_diff_vs_reference": (float(r["pct_diff_vs_reference"])
                                       if r.get("pct_diff_vs_reference") not in (None, "") else None),
        }
    return by_variant


def _select_family_rows(rows, tree_name, family):
    """Return {variant: row} for one method family ("adqsm" /
    "regperorder" / "calref_min5mm") for `tree_name`, keyed by
    _adqsm_variant_of(row) (shared with reconstruction_method_decision_
    summary.py, so both scripts agree on what a row's variant number is).

    Filtering is done on volume_results.csv's own STRUCTURED columns
    (calmethod) rather than regex-parsing the method string, wherever a
    structured column actually exists for that family - the same
    "structured over text-parsed" preference this codebase already
    applies elsewhere (see plot_box.py's TREEQSM_PARAM_SHORT_NAMES
    comment for the general rationale). AdQSM's own "(TreesParams)" rows
    are the one exception: they have no structured equivalent (calmethod/
    adqsm_variant are blank on them by design - see _adqsm_variant_of()'s
    own docstring), so the method-string prefix check already used
    elsewhere in this codebase (e.g. reconstruction_method_decision_
    summary.py's adqsm_other_variants()) is the only option there.

    SAFETY CHECK: raises SystemExit (does not silently pick one or
    average) if more than one row matches the SAME (variant, family) for
    this tree - e.g. multiple radius_threshold_mm/seg_* variants present
    for that (variant, calmethod) pair - since averaging/picking one
    would silently produce a chart that doesn't correspond to any single
    real reconstruction.
    """
    if family == "adqsm":
        picked = [r for r in rows if r["tree"] == tree_name
                  and r["method"].startswith("AdQSM (TreesParams)")]
    elif family == "regperorder":
        picked = [r for r in rows if r["tree"] == tree_name and r["calmethod"] == "regression-perorder"]
    elif family == "calref_min5mm":
        picked = [r for r in rows if r["tree"] == tree_name and r["calmethod"] == "min5mm"]
    else:
        raise ValueError("unknown family %r" % family)

    grouped = {}
    for r in picked:
        grouped.setdefault(_adqsm_variant_of(r), []).append(r)

    ambiguous = {variant: group for variant, group in grouped.items() if len(group) > 1}
    if ambiguous:
        lines = ["  variant %s: %d matching rows -> %s"
                 % (variant, len(group), [g["method"] for g in group])
                 for variant, group in sorted(ambiguous.items(), key=lambda kv: int(kv[0]))]
        raise SystemExit(
            "adqsm_variant_sensitivity.py: ambiguous rows for family=%r, tree=%r - "
            "more than one row matches (variant, method family) for %d variant(s). Add "
            "a filter (e.g. a RADIUS_THRESHOLD_MM/SEG_VARIANT parameter) to disambiguate "
            "before this script can build a trustworthy chart:\n%s"
            % (family, tree_name, len(ambiguous), "\n".join(lines)))

    return {variant: group[0] for variant, group in grouped.items()}


def _blank_if_none(value):
    """CSV convention already used throughout this codebase: a missing
    numeric value is written as "" (empty cell), never the literal string
    "None" - keeps the exported CSV machine-readable."""
    return "" if value is None else value


def build_merged_rows(tree_name):
    """Build one row dict per AdQSM variant present for `tree_name` in
    RESULTS_CSV's "AdQSM (TreesParams)" rows (the driving variant list -
    see the task's own "one row per AdQSM variant" framing), merged with
    the matching AdTree-calibrated rows and the taper_curve_compare_
    summary.csv row for that same variant. A variant missing from one of
    the OTHER sources (taper summary / regperorder / calref) does not
    raise - it is reported once as a batched warning and left blank in
    that row, same "degrade gracefully, never crash on an optional value"
    convention this codebase uses everywhere else (e.g. compare_volumes.
    py's load_results(), plot_box.py's None handling)."""
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit("'%s' not found." % RESULTS_CSV)
    rows = load_results(RESULTS_CSV)

    adqsm_by_variant = _select_family_rows(rows, tree_name, "adqsm")
    if not adqsm_by_variant:
        raise SystemExit("No 'AdQSM (TreesParams)' rows found for tree=%r in %r." % (tree_name, RESULTS_CSV))
    regperorder_by_variant = _select_family_rows(rows, tree_name, "regperorder")
    calref_by_variant = _select_family_rows(rows, tree_name, "calref_min5mm") if INCLUDE_CALREF_MIN5MM else {}

    taper_by_variant = read_taper_summary(TAPER_SUMMARY_CSV, tree_name)

    merged_rows = []
    missing_taper, missing_regperorder, missing_calref = [], [], []
    for variant in sorted(adqsm_by_variant, key=int):
        adqsm_row = adqsm_by_variant[variant]
        taper = taper_by_variant.get(variant)
        regp_row = regperorder_by_variant.get(variant)
        calref_row = calref_by_variant.get(variant) if INCLUDE_CALREF_MIN5MM else None

        if taper is None:
            missing_taper.append(variant)
        if regp_row is None:
            missing_regperorder.append(variant)
        if INCLUDE_CALREF_MIN5MM and calref_row is None:
            missing_calref.append(variant)

        row = {
            "variant": variant,
            "dbh_cm": _blank_if_none(taper["dbh_cm"] if taper is not None else None),
            "pct_diff_vs_reference": _blank_if_none(taper["pct_diff_vs_reference"] if taper is not None else None),
            "adqsm_total_m3": _blank_if_none(adqsm_row["total"]),
            "adqsm_trunk_m3": _blank_if_none(adqsm_row["trunk"]),
            "adtree_regperorder_total_m3": _blank_if_none(regp_row["total"] if regp_row is not None else None),
            "adtree_regperorder_trunk_m3": _blank_if_none(regp_row["trunk"] if regp_row is not None else None),
        }
        if INCLUDE_CALREF_MIN5MM:
            row["adtree_calref_total_m3"] = _blank_if_none(calref_row["total"] if calref_row is not None else None)
            row["adtree_calref_trunk_m3"] = _blank_if_none(calref_row["trunk"] if calref_row is not None else None)
        merged_rows.append(row)

    if missing_taper:
        print("WARNING: %d/%d variant(s) have no taper_curve_compare_summary.csv row for %r - "
              "dbh_cm/pct_diff_vs_reference left blank: %s"
              % (len(missing_taper), len(merged_rows), tree_name, missing_taper))
    if missing_regperorder:
        print("WARNING: %d/%d variant(s) have no '[calmethod=regression-perorder]' AdTree row for %r - "
              "adtree_regperorder_* left blank: %s"
              % (len(missing_regperorder), len(merged_rows), tree_name, missing_regperorder))
    if INCLUDE_CALREF_MIN5MM and missing_calref:
        print("WARNING: %d/%d variant(s) have no '[calref=min5mm]' AdTree row for %r - "
              "adtree_calref_* left blank: %s"
              % (len(missing_calref), len(merged_rows), tree_name, missing_calref))

    return merged_rows


def plot_variant_sensitivity(tree_name, merged_rows):
    """Log-scale line chart of total volume vs. AdQSM variant number,
    one line per series (AdQSM, AdTree regression-perorder, and -
    optionally - AdTree calref=min5mm). Saved as plots/<tree>/
    adqsm_variant_sensitivity_<tree>.png, the same plots/<tree>/ per-tree
    subfolder convention taper_curve_compare.py/parameter_sensitivity.py
    already use (ensure_plots_dir(), reused from plot_volumes.py, is the
    SAME helper those build on)."""
    x = [int(r["variant"]) for r in merged_rows]

    def series(key):
        return [r[key] if r[key] != "" else np.nan for r in merged_rows]

    # Colors: same families used everywhere else in the project
    # (FAMILY_GRADIENTS, plot_style.py) instead of matplotlib's generic
    # default color cycle - "AdQSM" gets that family's deepest/last shade
    # (a single solid line, so the boldest stop reads best); the two
    # AdTree-calibrated series share the "AdTree calibrated" family but
    # use its two deepest stops so they stay visually grouped together
    # while remaining distinguishable from each other - the deeper one
    # goes to regression-perorder since that's the PRIMARY calibration
    # method project-wide (calref=min5mm is the secondary/backup one).
    adqsm_color = FAMILY_GRADIENTS["AdQSM"][-1]
    regperorder_color = FAMILY_GRADIENTS["AdTree calibrated"][-1]
    calref_color = FAMILY_GRADIENTS["AdTree calibrated"][-2]

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(x, series("adqsm_total_m3"), marker="o", markersize=CHART_MARKERSIZE,
            linewidth=CHART_LINEWIDTH, color=adqsm_color, label="AdQSM")
    ax.plot(x, series("adtree_regperorder_total_m3"), marker="o", markersize=CHART_MARKERSIZE,
            linewidth=CHART_LINEWIDTH, color=regperorder_color,
            label="AdTree calibrated (regression-perorder)")
    if INCLUDE_CALREF_MIN5MM:
        ax.plot(x, series("adtree_calref_total_m3"), marker="o", markersize=CHART_MARKERSIZE,
                linewidth=CHART_LINEWIDTH, color=calref_color,
                label="AdTree calibrated (calref=min5mm)")

    ax.set_yscale("log")
    ax.set_xlabel("AdQSM variant", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Total volume [m^3] (log scale)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title("%s: total volume vs. AdQSM variant (%d variants)" % (tree_name, len(merged_rows)),
                 fontsize=PANEL_TITLE_FONTSIZE)
    ax.legend(fontsize=CHART_LEGEND_FONTSIZE, loc="best")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()

    out_dir = os.path.join(ensure_plots_dir(), tree_name)
    os.makedirs(out_dir, exist_ok=True)   # only ever creates this tree's own subfolder - never touches any other file/folder under plots/
    png_path = os.path.join(out_dir, "adqsm_variant_sensitivity_%s.png" % tree_name)
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", png_path)
    return out_dir, png_path


def write_merged_csv(out_dir, tree_name, merged_rows):
    """Save the full merged table next to the PNG (out_dir, i.e.
    plots/<tree>/) as adqsm_variant_sensitivity_<tree>.csv, sorted by
    variant (merged_rows already is, from build_merged_rows())."""
    fieldnames = ["variant", "dbh_cm", "pct_diff_vs_reference",
                  "adqsm_total_m3", "adqsm_trunk_m3",
                  "adtree_regperorder_total_m3", "adtree_regperorder_trunk_m3"]
    if INCLUDE_CALREF_MIN5MM:
        fieldnames += ["adtree_calref_total_m3", "adtree_calref_trunk_m3"]

    csv_path = os.path.join(out_dir, "adqsm_variant_sensitivity_%s.csv" % tree_name)
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(merged_rows)
    print("Saved:", csv_path)
    return csv_path


def print_cluster_summary(tree_name, merged_rows):
    """Console summary (Step 6): how many variants are "close" to the
    measured reference DBH (|pct_diff_vs_reference| < threshold) vs.
    "outliers" (>= threshold), and the min/max/mean AdQSM total volume
    within each group - mirrors the ad hoc cluster split already found
    for B21_S01 (44 close / 17 outliers on a dbh<65cm-ish cutoff; this
    uses pct_diff_vs_reference instead, so exact counts may differ
    slightly)."""
    have_pct = [r for r in merged_rows if r["pct_diff_vs_reference"] != ""]
    skipped = [r for r in merged_rows if r["pct_diff_vs_reference"] == ""]
    close_group = [r for r in have_pct if abs(r["pct_diff_vs_reference"]) < PCT_DIFF_OUTLIER_THRESHOLD_PCT]
    outlier_group = [r for r in have_pct if abs(r["pct_diff_vs_reference"]) >= PCT_DIFF_OUTLIER_THRESHOLD_PCT]

    def vol_stats(group):
        vols = [r["adqsm_total_m3"] for r in group if r["adqsm_total_m3"] != ""]
        if not vols:
            return None
        return min(vols), max(vols), sum(vols) / len(vols)

    print()
    print("=" * 78)
    print("%s: %d/%d variants within +/-%.0f%% of the reference DBH ('close'), "
          "%d 'outlier' (threshold on pct_diff_vs_reference)"
          % (tree_name, len(close_group), len(have_pct), PCT_DIFF_OUTLIER_THRESHOLD_PCT, len(outlier_group)))
    if skipped:
        print("  (%d variant(s) excluded from this split - no taper_curve_compare_summary.csv match: %s)"
              % (len(skipped), [r["variant"] for r in skipped]))
    close_stats = vol_stats(close_group)
    outlier_stats = vol_stats(outlier_group)
    if close_stats:
        print("  close cluster   AdQSM total_m3: min=%.4f  max=%.4f  mean=%.4f" % close_stats)
    else:
        print("  close cluster   AdQSM total_m3: (no variants with a valid value)")
    if outlier_stats:
        print("  outlier cluster AdQSM total_m3: min=%.4f  max=%.4f  mean=%.4f" % outlier_stats)
    else:
        print("  outlier cluster AdQSM total_m3: (no variants with a valid value)")
    print("=" * 78)


def run(tree_name):
    print("=" * 78)
    print("Tree: %s  (INCLUDE_CALREF_MIN5MM=%s)" % (tree_name, INCLUDE_CALREF_MIN5MM))
    print("-" * 78)
    merged_rows = build_merged_rows(tree_name)
    out_dir, _png_path = plot_variant_sensitivity(tree_name, merged_rows)
    write_merged_csv(out_dir, tree_name, merged_rows)
    print_cluster_summary(tree_name, merged_rows)


# =========================  RUN  =====================================
if __name__ == "__main__":
    run(TREE_NAME)
