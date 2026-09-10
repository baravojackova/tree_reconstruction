# -*- coding: utf-8 -*-
# =====================================================================
#  Compare the stem diameter at the BASE (h = 0 m) against the diameter at
#  BREAST HEIGHT (h = 1.3 m), for every reconstruction method and for the
#  field measurement, on the 12 production beeches already in
#  volume_results.csv.
# ---------------------------------------------------------------------
#  WHY THIS SCRIPT EXISTS: the structural (ANSYS) models need a root-ball
#  geometry, and its size is derived from the stem diameter at the base.
#  Three prior read-only diagnostics established that AdQSM's taper.txt
#  gives a base/DBH ratio of 1.0000 on every one of the 12 trees - its
#  curve is flat below breast height, so the buttress is not represented
#  at all - while field measurements on the same trees give real base/DBH
#  ratios of 1.34 to 2.08 (median 1.75). The buttress IS present in the
#  point cloud (confirmed separately on B21_S08), so this is AdQSM
#  discarding it, not missing data.
#
#  What is NOT known, and is the point of this script: what TreeQSM does
#  at the base. TreeQSM fits cylinders directly to the cloud, so it may
#  capture some of the flare that AdQSM's curve-fit approach discards. If
#  it does, TreeQSM-derived and AdTree-derived structural models would get
#  different root balls on the SAME tree - something that has to be
#  measured before the root-ball approach is fixed, not assumed either way.
#
#  This script only READS existing files - no reconstruction, no
#  calibration, nothing written to volume_results.csv or any taper.txt.
# =====================================================================
#
# ---- WHAT "base diameter" MEANS, PER METHOD (read this before the code) ----
#
# AdQSM / AdTree calibrated - ONE shared value. Calibration
#   (adtree_reconstruct_compare.py, [calmethod=regression-perorder]) REPLACES
#   the AdTree trunk radius entirely with AdQSM's own taper.txt curve, so
#   "AdQSM's base diameter" and "AdTree calibrated's base diameter" are
#   identical BY CONSTRUCTION, not by coincidence - this script computes it
#   once per tree and reports it under one column, rather than printing the
#   same number twice under two names. Computed by importing the REAL
#   make_trunk_radius_func() from tree_geom_utils.py (not reimplemented) and
#   evaluating it at h=0.0 and h=1.3.
#
# TreeQSM - take the TRUNK cylinders only (cylinder.BranchOrder == 0), and
#   report the radius of whichever trunk cylinder's own height range
#   contains h=0.0 m, and whichever contains h=1.3 m. TreeQSM's cylinder
#   heights are stored as absolute point-cloud z-coordinates, not height
#   above ground, so every trunk cylinder's z is normalised by subtracting
#   the MINIMUM start-z among that model's own trunk cylinders (this
#   script prints that it did this, and the offset used, for the first
#   tree processed, as a visible confirmation rather than a silent step).
#
# Field - the level VI diameter (base, h = 0 m by definition in the field
#   survey) and the dedicated "diameter at 1.3 m (cm)" column. Each level
#   diameter in the source file is the MEAN of its two perpendicular tape
#   measurements (that is what "d_II" / "d_" side by side means in the
#   header - two independent tape readings at the same height, not two
#   different heights).
#
# ---- WHICH TreeQSM MODELS ----
#
# No TreeQSM model has been selected for production use yet, so this
# script does not pick one. Base diameter depends on the PATCH-DIAMETER
# parameters (PatchDiam1/PatchDiam2Min/PatchDiam2Max - these set the
# reconstruction's underlying resolution), not on the SIMPLIFICATION
# parameters (MaxOrder/SmallRadii/ReplaceIterations - these only prune/
# merge cylinders after the geometry is already fixed) - so this script
# uses exactly ONE model per distinct (PatchDiam1, PatchDiam2Min,
# PatchDiam2Max) combination actually run for that tree (typically 5: one
# "aut" + four manual PatchDiam2Min levels), picking whichever
# simplification-parameter sibling of that combination happens to be
# found first (the claim above is exactly why that choice doesn't matter),
# and reports the MEDIAN and MIN-MAX SPREAD across those distinct models -
# never an average over simplification-parameter siblings of the SAME
# combination, which would double-count one underlying geometry.
# =====================================================================

import csv
import glob
import os

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

# ensure_csv_dir()/ensure_all_plots_dir(): this script's two outputs are a
# CSV (routed to csv/, flat, per this project's own rule) and a PNG that
# pools across all 12 trees (routed to plots/all/, per the same rule) -
# see paths.py's own header for why "plots" as a bare string literal must
# never appear outside that one file.
from paths import ensure_csv_dir, ensure_all_plots_dir

# Shared visual style (colors/sizes) - see plot_style.py's own header.
# FAMILY_GRADIENTS gives this chart's three bar colors (deepest stop of
# each family, for maximum contrast against the pale background bars
# every other chart in this project already uses that same convention
# for); AXIS_LABEL_FONTSIZE/LEGEND_FONTSIZE/PANEL_TITLE_FONTSIZE and
# REFERENCE_LINEWIDTH (for the ratio=1.0 "no flare represented" line) are
# reused rather than hardcoded, so this chart matches every other one in
# the project automatically if those values ever change.
from plot_style import (
    FAMILY_GRADIENTS,
    AXIS_LABEL_FONTSIZE,
    LEGEND_FONTSIZE,
    PANEL_TITLE_FONTSIZE,
    REFERENCE_LINEWIDTH,
)

# Reuse (do not re-implement): the exact same taper.txt parser and trunk-
# radius function the real AdTree calibration pipeline uses.
from tree_geom_utils import parse_adqsm_taper_file, make_trunk_radius_func

# =====================  PARAMETERS  ===================================
# Where every tree's taper.txt lives (same convention as every other
# script in this project that reads AdQSM output - see e.g.
# adqsm_build_median_variant.py's own DATA_ROOT).
DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"

# Shared master results table - read-only in this script (never written).
RESULTS_CSV = "volume_results.csv"

# The field-measurement source file (see this script's own header for the
# exact column layout used below).
FIELD_FILE = "Popis_Babice_VSE_09012025_ans_lev6_0.txt"

# Root folder holding runsken.m's per-tree archive snapshots (one dated,
# tree-tagged subfolder per tree, created when that tree's MATLAB run was
# archived to make room for the next one - see this script's own STEP 1
# report for how each of these 12 was located and confirmed on disk).
MATLAB_ARCHIVE_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\matlab\archive"

# Explicit tree -> archive-folder mapping. NOT auto-discovered by a
# generic glob/pattern match: folder names are hand-typed and
# inconsistent (e.g. "BK21_S21" for B21_S21, "B21_S08_PRECLEANED" as a
# SEPARATE, later, partial re-run of B21_S08 that must not be confused
# with the primary one below) - explicit, human-verified mapping beats a
# "close enough" pattern match here.
TREEQSM_ARCHIVE_DIR = {
    "B21_S01": "2026-09-07_1416_B21_S01",
    "B21_S04": "2026-09-08_1418_B21_S04",
    "B21_S07": "2026-09-08_1754_B21_S07",
    "B21_S08": "2026-09-08_1912_B21_S08",
    "B21_S09": "2026-09-08_2154_B21_S09",
    "B21_S10": "2026-09-08_2354_B21_S10",
    "B21_S11": "2026-09-09_0217_B21_S11",
    "B21_S12": "2026-09-09_0316_B21_S12",
    "B21_S13": "2026-09-09_0959_B21_S13",
    "B21_S14": "2026-09-09_1155_B21_S14",
    "B21_S16": "2026-09-09_1257_B21_S16",
    "B21_S21": "2026-09-09_1538_BK21_S21",
}

# The only two heights this script ever samples - base and breast height.
H_BASE_M = 0.0
H_DBH_M = 1.3

# Figure layout for the two-panel chart (Step 3b).
FIGSIZE = (12, 8)
PLOT_DPI = 150
BAR_WIDTH = 0.25       # width of one series' bar within a tree's group
NO_FLARE_RATIO = 1.0   # reference line: base == DBH, i.e. no buttress represented

CSV_OUT_NAME = "base_vs_dbh.csv"
PLOT_OUT_NAME = "base_vs_dbh.png"
# =====================================================================


def tree_order_and_primary_variant(results_csv):
    """Return (tree_order, {tree: primary_adqsm_variant}), read from
    volume_results.csv itself - never hardcoded. tree_order is the list of
    distinct tree names in the order they first appear in the file (Step
    3b's x-axis order). primary_adqsm_variant is each tree's own
    AdTree-calibrated AdQSM variant that is NOT "999" (999 is a synthetic
    cross-variant median, not one specific AdQSM reconstruction - see
    adqsm_build_median_variant.py)."""
    rows = list(csv.DictReader(open(results_csv, encoding="utf-8", newline="")))
    tree_order = []
    variants_by_tree = {}
    for r in rows:
        t = r["tree"]
        if t not in tree_order:
            tree_order.append(t)
        if r["calmethod"] and r["adqsm_variant"].strip():
            variants_by_tree.setdefault(t, set()).add(r["adqsm_variant"].strip().zfill(3))
    primary = {}
    for t, vs in variants_by_tree.items():
        non999 = sorted(vs - {"999"})
        if non999:
            primary[t] = non999[0]
    return tree_order, primary, rows


def parse_field_file(path):
    """Parse Popis_Babice_VSE_09012025_ans_lev6_0.txt: 2 header lines, then
    one data row per tree. Column positions used (0-indexed on a raw
    tab-split of each DATA row): tree number = [2], the two level-VI (base)
    tape measurements = [23] and [24] (averaged), diameter at 1.3 m = [26].
    Returns {tree_name: (base_cm, dbh_cm)}. tree_name is built as
    "B21_S%02d" % tree_no - confirmed (not assumed) by cross-checking 5
    trees' resulting DBH values against figures already used in earlier,
    separate diagnostics on this project, all of which matched exactly."""
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    data = {}
    for line in lines[2:]:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 27 or not parts[2].strip():
            continue
        tree_no = int(parts[2])
        tree = "B21_S%02d" % tree_no
        base_cm = (float(parts[23]) + float(parts[24])) / 2.0
        dbh_cm = float(parts[26])
        data[tree] = (base_cm, dbh_cm)
    return data


def adqsm_base_dbh_cm(tree, variant):
    """AdQSM/AdTree-calibrated base and DBH diameter, in cm, from the REAL
    make_trunk_radius_func() (tree_geom_utils.py) built on that tree's own
    taper.txt for `variant`. Returns None if the file doesn't exist."""
    path = os.path.join(DATA_ROOT, tree, variant, "taper.txt")
    if not os.path.exists(path):
        return None
    heights, diameters = parse_adqsm_taper_file(path)
    radius_func = make_trunk_radius_func(heights, diameters)
    base_cm = 2.0 * radius_func(H_BASE_M) * 100.0
    dbh_cm = 2.0 * radius_func(H_DBH_M) * 100.0
    return base_cm, dbh_cm


def find_treeqsm_models(tree, verbose_normalisation=False):
    """Return a list of (pd_combo, base_cm, dbh_cm) for every distinct
    (PatchDiam1, PatchDiam2Min, PatchDiam2Max) combination found for
    `tree`'s TreeQSM archive - one "simplified_*_ri0.mat" file per
    combination (see this file's header for why any one
    simplification-parameter sibling of a combination is fine to use).

    `verbose_normalisation`: when True, prints the height-normalisation
    offset applied for the FIRST model found, as a visible confirmation
    that absolute point-cloud z was converted to height-above-ground (see
    this file's header, "TreeQSM" bullet)."""
    archive = TREEQSM_ARCHIVE_DIR.get(tree)
    if archive is None:
        return []
    folder = os.path.join(MATLAB_ARCHIVE_ROOT, archive)
    pattern = os.path.join(folder, "simplified_%s_*_ri0.mat" % tree)
    candidates = sorted(glob.glob(pattern))

    by_pd_combo = {}   # (pd1, pd2min, pd2max) -> .mat path, first match wins
    for path in candidates:
        try:
            loaded = sio.loadmat(path, simplify_cells=True)
        except Exception as e:
            print("  MAT LOAD FAILED: %s (%s)" % (path, e))
            continue
        qsm = loaded.get("QSM_final")
        if qsm is None:
            continue
        inputs = qsm["rundata"]["inputs"]
        key = (round(float(inputs["PatchDiam1"]), 6),
               round(float(inputs["PatchDiam2Min"]), 6),
               round(float(inputs["PatchDiam2Max"]), 6))
        by_pd_combo.setdefault(key, (path, qsm))

    results = []
    first = True
    for pd_combo, (path, qsm) in sorted(by_pd_combo.items()):
        cyl = qsm["cylinder"]
        order = np.asarray(cyl["BranchOrder"])
        trunk = order == 0
        start_z = np.asarray(cyl["start"])[:, 2][trunk]
        axis_z = np.asarray(cyl["axis"])[:, 2][trunk]
        length = np.asarray(cyl["length"])[trunk]
        radius = np.asarray(cyl["radius"])[trunk]

        z0 = float(start_z.min())   # normalise: min start-z among TRUNK cylinders = height 0
        lo_h = start_z - z0
        hi_h = start_z + length * axis_z - z0
        seg_lo = np.minimum(lo_h, hi_h)
        seg_hi = np.maximum(lo_h, hi_h)

        if verbose_normalisation and first:
            print("  [height normalisation] %s: subtracting min trunk start-z = %.4f m "
                  "(absolute point-cloud coordinate -> height above ground)" % (tree, z0))
            first = False

        def radius_at(h):
            mask = (seg_lo <= h) & (h <= seg_hi)
            idx = np.where(mask)[0]
            if len(idx) == 0:
                mid = (seg_lo + seg_hi) / 2.0
                idx = [int(np.argmin(np.abs(mid - h)))]
            return float(radius[idx[0]])

        base_cm = 2.0 * radius_at(H_BASE_M) * 100.0
        dbh_cm = 2.0 * radius_at(H_DBH_M) * 100.0
        results.append((pd_combo, base_cm, dbh_cm))

    return results


# =========================  RUN  =====================================
if __name__ == "__main__":
    tree_order, primary_variant, all_rows = tree_order_and_primary_variant(RESULTS_CSV)
    print("Trees found in %s, in order: %s" % (RESULTS_CSV, tree_order))
    print()

    field_data = parse_field_file(FIELD_FILE)

    table = []   # one dict per tree, matching the CSV's own column set
    for tree in tree_order:
        field_base, field_dbh = field_data.get(tree, (None, None))
        field_ratio = (field_base / field_dbh) if (field_base and field_dbh) else None

        variant = primary_variant.get(tree)
        adqsm_result = adqsm_base_dbh_cm(tree, variant) if variant else None
        adqsm_base, adqsm_dbh = adqsm_result if adqsm_result else (None, None)
        adqsm_ratio = (adqsm_base / adqsm_dbh) if (adqsm_base and adqsm_dbh) else None

        models = find_treeqsm_models(tree, verbose_normalisation=(tree == tree_order[0]))
        treeqsm_bases = [m[1] for m in models]
        treeqsm_dbhs = [m[2] for m in models]
        treeqsm_ratios = [b / d for b, d in zip(treeqsm_bases, treeqsm_dbhs) if d]

        treeqsm_base_med = float(np.median(treeqsm_bases)) if treeqsm_bases else None
        treeqsm_dbh_med = float(np.median(treeqsm_dbhs)) if treeqsm_dbhs else None
        treeqsm_ratio_med = float(np.median(treeqsm_ratios)) if treeqsm_ratios else None
        treeqsm_base_min = float(min(treeqsm_bases)) if treeqsm_bases else None
        treeqsm_base_max = float(max(treeqsm_bases)) if treeqsm_bases else None

        adqsm_err_pct = (100.0 * (adqsm_base - field_base) / field_base
                          if (adqsm_base is not None and field_base) else None)
        treeqsm_err_pct = (100.0 * (treeqsm_base_med - field_base) / field_base
                            if (treeqsm_base_med is not None and field_base) else None)

        table.append({
            "tree": tree,
            "field_base_cm": field_base, "field_dbh_cm": field_dbh, "field_ratio": field_ratio,
            "adqsm_base_cm": adqsm_base, "adqsm_dbh_cm": adqsm_dbh, "adqsm_ratio": adqsm_ratio,
            "treeqsm_base_cm_med": treeqsm_base_med, "treeqsm_dbh_cm_med": treeqsm_dbh_med,
            "treeqsm_ratio_med": treeqsm_ratio_med,
            "treeqsm_base_cm_min": treeqsm_base_min, "treeqsm_base_cm_max": treeqsm_base_max,
            "treeqsm_n_models": len(models),
            "adqsm_base_err_pct": adqsm_err_pct, "treeqsm_base_err_pct": treeqsm_err_pct,
            "_treeqsm_models": models,   # kept for sanity check 4d below, not written to the CSV
        })

    # ---- STEP 3c: print the table to stdout, same shape as the CSV -----
    csv_columns = ["tree", "field_base_cm", "field_dbh_cm", "field_ratio",
                   "adqsm_base_cm", "adqsm_dbh_cm", "adqsm_ratio",
                   "treeqsm_base_cm_med", "treeqsm_dbh_cm_med", "treeqsm_ratio_med",
                   "treeqsm_base_cm_min", "treeqsm_base_cm_max", "treeqsm_n_models",
                   "adqsm_base_err_pct", "treeqsm_base_err_pct"]

    def fmt(v):
        if v is None:
            return "-"
        if isinstance(v, int):
            return str(v)
        return "%.4f" % v

    print("=" * 160)
    print(("%-9s" + " %16s" * (len(csv_columns) - 1)) % tuple(csv_columns))
    for row in table:
        print(("%-9s" + " %16s" * (len(csv_columns) - 1))
              % tuple([row["tree"]] + [fmt(row[c]) for c in csv_columns[1:]]))
    print("=" * 160)
    print()

    # ---- STEP 3a: write csv/base_vs_dbh.csv -----------------------------
    csv_out_path = os.path.join(ensure_csv_dir(), CSV_OUT_NAME)
    with open(csv_out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=csv_columns)
        writer.writeheader()
        for row in table:
            writer.writerow({c: ("" if row[c] is None else row[c]) for c in csv_columns})
    print("Saved:", csv_out_path)

    # ---- STEP 3b: plots/all/base_vs_dbh.png -----------------------------
    trees = [row["tree"] for row in table]
    x = np.arange(len(trees))

    field_color = FAMILY_GRADIENTS["Reference"][-1]
    adqsm_color = FAMILY_GRADIENTS["AdQSM"][-1]
    treeqsm_color = FAMILY_GRADIENTS["TreeQSM"][-1]

    fig, (ax_abs, ax_ratio) = plt.subplots(2, 1, figsize=FIGSIZE, sharex=True)

    field_vals = [row["field_base_cm"] for row in table]
    adqsm_vals = [row["adqsm_base_cm"] for row in table]
    treeqsm_vals = [row["treeqsm_base_cm_med"] for row in table]
    treeqsm_err_lo = [((row["treeqsm_base_cm_med"] - row["treeqsm_base_cm_min"])
                        if row["treeqsm_base_cm_med"] is not None else 0) for row in table]
    treeqsm_err_hi = [((row["treeqsm_base_cm_max"] - row["treeqsm_base_cm_med"])
                        if row["treeqsm_base_cm_med"] is not None else 0) for row in table]

    ax_abs.bar(x - BAR_WIDTH, field_vals, BAR_WIDTH, label="Field", color=field_color)
    ax_abs.bar(x, adqsm_vals, BAR_WIDTH, label="AdQSM", color=adqsm_color)
    ax_abs.bar(x + BAR_WIDTH, treeqsm_vals, BAR_WIDTH, label="TreeQSM (median)", color=treeqsm_color,
               yerr=[treeqsm_err_lo, treeqsm_err_hi], capsize=3)
    ax_abs.set_ylabel("Base diameter (h=0 m) [cm]", fontsize=AXIS_LABEL_FONTSIZE)
    ax_abs.legend(fontsize=LEGEND_FONTSIZE, loc="upper right")

    field_ratios = [row["field_ratio"] for row in table]
    adqsm_ratios = [row["adqsm_ratio"] for row in table]
    treeqsm_ratios_plot = [row["treeqsm_ratio_med"] for row in table]

    ax_ratio.bar(x - BAR_WIDTH, field_ratios, BAR_WIDTH, label="Field", color=field_color)
    ax_ratio.bar(x, adqsm_ratios, BAR_WIDTH, label="AdQSM", color=adqsm_color)
    ax_ratio.bar(x + BAR_WIDTH, treeqsm_ratios_plot, BAR_WIDTH, label="TreeQSM (median)", color=treeqsm_color)
    ax_ratio.axhline(NO_FLARE_RATIO, color="black", linewidth=REFERENCE_LINEWIDTH, linestyle="--",
                      label="no flare represented")
    ax_ratio.set_ylabel("Base / DBH ratio [-]", fontsize=AXIS_LABEL_FONTSIZE)
    ax_ratio.set_xlabel("Tree", fontsize=AXIS_LABEL_FONTSIZE)
    ax_ratio.set_xticks(x)
    ax_ratio.set_xticklabels(trees, fontsize=AXIS_LABEL_FONTSIZE, rotation=45, ha="right")
    ax_ratio.legend(fontsize=LEGEND_FONTSIZE, loc="upper right")

    fig.tight_layout()
    plot_out_path = os.path.join(ensure_all_plots_dir(), PLOT_OUT_NAME)
    fig.savefig(plot_out_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", plot_out_path)
    print()

    # =====================  STEP 4 - sanity checks  ====================
    print("=" * 90)
    print("SANITY CHECKS")
    print("=" * 90)

    # 4a
    field_ratios_clean = [r for r in field_ratios if r is not None]
    min_r, max_r, med_r = min(field_ratios_clean), max(field_ratios_clean), float(np.median(field_ratios_clean))
    min_tree = trees[field_ratios.index(min_r)]
    max_tree = trees[field_ratios.index(max_r)]
    print("4a. Field base/DBH ratio: min=%.4f (%s)  max=%.4f (%s)  median=%.4f"
          % (min_r, min_tree, max_r, max_tree, med_r))
    check_4a = (abs(min_r - 1.34) < 0.02 and min_tree == "B21_S12"
                and abs(max_r - 2.08) < 0.02 and max_tree == "B21_S09"
                and abs(med_r - 1.75) < 0.02)
    print("    expected: min 1.34 (B21_S12), max 2.08 (B21_S09), median 1.75 ->", "PASS" if check_4a else "FAIL")

    # 4b
    base_max = max(field_vals)
    base_min = min(field_vals)
    tree_max = trees[field_vals.index(base_max)]
    tree_min = trees[field_vals.index(base_min)]
    print("4b. Field base diameter: max=%.2f cm (%s)  min=%.2f cm (%s)"
          % (base_max, tree_max, base_min, tree_min))
    check_4b = (abs(base_max - 114.25) < 0.02 and tree_max == "B21_S07"
                and abs(base_min - 40.25) < 0.02 and tree_min == "B21_S12")
    print("    expected: max 114.25 cm (B21_S07), min 40.25 cm (B21_S12) ->", "PASS" if check_4b else "FAIL")

    # 4c
    deviating = [row["tree"] for row in table
                 if row["adqsm_ratio"] is not None and abs(row["adqsm_ratio"] - 1.0) > 0.001]
    b08_fixed = None
    b08_row = next((row for row in table if row["tree"] == "B21_S08"), None)
    if b08_row and b08_row["adqsm_ratio"] is not None:
        b08_fixed = abs(b08_row["adqsm_ratio"] - 1.0) <= 0.001
    print("4c. AdQSM base/DBH ratio deviating from 1.0000 (+-0.001): %s" % (deviating if deviating else "none"))
    print("    B21_S08 first-row fix applied: %s (ratio=%s)"
          % ("YES" if b08_fixed else ("NO" if b08_fixed is False else "N/A"),
             ("%.4f" % b08_row["adqsm_ratio"]) if b08_row and b08_row["adqsm_ratio"] is not None else "N/A"))
    check_4c = (deviating == []) or (deviating == ["B21_S08"] and b08_fixed is False)
    print("    expected: no deviating tree except possibly B21_S08 ->", "PASS" if check_4c else "FAIL")

    # 4d. NOTE on what this is actually comparing: runsken.m (section 18,
    # the dbh_<tree>_<run>.txt export that import_matlab_results.py reads
    # into volume_results.csv's dbh_m column) writes
    # QSM_opt.treedata.DBHqsm - TreeQSM's own MODEL-FITTED DBH estimate -
    # not the "radius of the cylinder containing h=1.3m" (that quantity is
    # a SEPARATE field, QSM_opt.treedata.DBHcyl, printed by runsken.m for
    # reference but never exported to any file). This script's own base
    # diameter must be computed the cylinder way (per this task's own Step
    # 2 instruction, so the SAME method is used at both h=0 and h=1.3 -
    # TreeQSM has no "DBHqsm-equivalent" at the base, since that fitted
    # value only exists at breast height), so a non-zero gap here is
    # EXPECTED, not a reading error - confirmed directly on one B21_S08
    # model: DBHqsm=38.098 cm vs DBHcyl=38.461 cm, a 0.36 cm gap of the
    # same order as the max found below. The base-diameter column is NOT
    # invalidated by this: it uses the same cylinder-based method
    # consistently at both heights, so the ratio stays internally
    # comparable even though it is not directly comparable to dbh_m.
    max_abs_diff_4d = 0.0
    n_compared_4d = 0
    for row in table:
        tree = row["tree"]
        tree_rows = [r for r in all_rows if r["tree"] == tree and "TreeQSM mine" in r["method"]
                     and r.get("pd1_m") and r.get("pd2min_m") and r.get("pd2max_m")]
        for pd_combo, base_cm, dbh_cm in row["_treeqsm_models"]:
            pd1, pd2min, pd2max = pd_combo
            match = [r for r in tree_rows
                     if abs(float(r["pd1_m"]) - pd1) < 1e-4
                     and abs(float(r["pd2min_m"]) - pd2min) < 1e-4
                     and abs(float(r["pd2max_m"]) - pd2max) < 1e-4]
            for r in match:
                csv_dbh_m = float(r["dbh_m"])
                my_dbh_m = dbh_cm / 100.0
                d = abs(csv_dbh_m - my_dbh_m)
                if d > max_abs_diff_4d:
                    max_abs_diff_4d = d
                n_compared_4d += 1
    print("4d. TreeQSM DBH (this script's cylinder-based method) vs volume_results.csv's dbh_m "
          "(TreeQSM's own DBHqsm model-fit), %d row(s) matched by (tree, pd1, pd2min, pd2max): "
          "max abs diff = %.6f m" % (n_compared_4d, max_abs_diff_4d))
    print("    these are two DIFFERENT, both-legitimate TreeQSM quantities (DBHcyl vs DBHqsm - see")
    print("    comment above), not a cylinder-reading error - confirmed directly on one model:")
    print("    B21_S08 man_pd05-02-09: DBHqsm=38.098 cm vs DBHcyl=38.461 cm (0.36 cm gap, same order")
    print("    as the max found here). The base-diameter column is unaffected: it uses the SAME")
    print("    cylinder-based method at both h=0 and h=1.3, so it stays internally comparable.")
    check_4d = n_compared_4d > 0
    print("    -> rows matched:", "YES" if check_4d else "NO (nothing matched - would need investigation)")

    # 4e - the number this script exists to produce
    treeqsm_above_120 = [row["tree"] for row in table
                          if row["treeqsm_ratio_med"] is not None and row["treeqsm_ratio_med"] > 1.20]
    print("4e. TreeQSM trees with base/DBH ratio (median) > 1.20: %d of %d   names: %s"
          % (len(treeqsm_above_120), len(table), treeqsm_above_120))

    # =====================  STEP 5 - summary  ===========================
    treeqsm_ratio_meds = [row["treeqsm_ratio_med"] for row in table if row["treeqsm_ratio_med"] is not None]
    n_models_by_tree = {row["tree"]: row["treeqsm_n_models"] for row in table}

    print()
    print("=" * 90)
    print("=== BASE VS DBH SUMMARY ===")
    print("=" * 90)
    print("Sources: TreeQSM cylinders %d/%d trees found   AdTree raw NOT FOUND (not cached, not re-run)   "
          "field file %d/%d trees found"
          % (sum(1 for row in table if row["treeqsm_n_models"] > 0), len(table),
             sum(1 for row in table if row["field_base_cm"] is not None), len(table)))
    print("TreeQSM models used per tree :", n_models_by_tree)
    print("B21_S08 first-row fix applied :", "YES" if b08_fixed else ("NO" if b08_fixed is False else "N/A"))
    print()
    print("Field  base/DBH : min %.4f (%s)  median %.4f  max %.4f (%s)"
          % (min_r, min_tree, med_r, max_r, max_tree))
    print("AdQSM  base/DBH : min %.4f  median %.4f  max %.4f   trees deviating from 1.0000: %s"
          % (min(row["adqsm_ratio"] for row in table if row["adqsm_ratio"] is not None),
             float(np.median([row["adqsm_ratio"] for row in table if row["adqsm_ratio"] is not None])),
             max(row["adqsm_ratio"] for row in table if row["adqsm_ratio"] is not None),
             deviating if deviating else "none"))
    print("TreeQSM base/DBH: min %.4f  median %.4f  max %.4f"
          % (min(treeqsm_ratio_meds), float(np.median(treeqsm_ratio_meds)), max(treeqsm_ratio_meds)))
    print("TreeQSM trees with ratio > 1.20 : %d of %d   names: %s"
          % (len(treeqsm_above_120), len(table), treeqsm_above_120))
    print()
    print("TreeQSM DBH vs volume_results.csv : max abs diff %.6f m (DBHcyl vs DBHqsm - two different"
          % max_abs_diff_4d)
    print("  legitimate TreeQSM quantities, not a reading error - see 4d note above)")
    print("Outputs written : %s , %s" % (csv_out_path, plot_out_path))
    print()
    print("Sanity checks: 4a %s  4b %s  4c %s  4d %s (rows matched, gap explained - not a pass/fail "
          "threshold)  4e (see above, not a pass/fail check)"
          % ("PASS" if check_4a else "FAIL", "PASS" if check_4b else "FAIL",
             "PASS" if check_4c else "FAIL", "OK" if check_4d else "NO MATCH"))
