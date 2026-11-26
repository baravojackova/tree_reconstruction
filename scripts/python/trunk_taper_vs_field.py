# -*- coding: utf-8 -*-
# =====================================================================
#  Compare the reconstructed trunk TAPER CURVE (not just base/DBH) against
#  the field-measured taper, at all six field measurement heights, on the
#  12 production beeches already in volume_results.csv.
# ---------------------------------------------------------------------
#  WHY THIS SCRIPT EXISTS: base_vs_dbh_compare.py compared exactly two
#  heights (0 m and 1.3 m) and found AdQSM's base/DBH ratio is 1.0000 on
#  every tree while the field gives 1.34-2.08. Two heights are not enough
#  to say anything about the trunk as a whole. The field dataset has SIX
#  measured heights per tree (0 m up to between 4.6 m and 9.8 m depending
#  on the tree) - an independent reference for the whole lower trunk that
#  has not been used yet. This script only READS existing files - no
#  reconstruction, no calibration, nothing written to volume_results.csv
#  or any taper.txt.
#
#  IMPORTANT NAMING PRECISION, preserved throughout this file's output
#  (column names, figure legends): the AdTree-CALIBRATED trunk is not
#  AdTree's own measurement. Calibration (adtree_reconstruct_compare.py,
#  [calmethod=regression-perorder]) REPLACES the AdTree trunk radius
#  ENTIRELY with AdQSM's own taper.txt curve (see
#  AdTree_Calmethod_Decision.docx S5) - so "AdQSM's taper" and "AdTree
#  calibrated's taper" are identical by construction, not by coincidence.
#  This script therefore has ONE column/series for both, named
#  "adqsm_taper" (never "adtree") - printing a separate "AdTree" column
#  or line would wrongly imply two independent measurements exist.
#
#  REUSE NOTE (see also this script's own Step 1 sanity-check print):
#  base_vs_dbh_compare.py's tree_order_and_primary_variant() IS imported
#  and reused directly below - no duplication there. Its
#  parse_field_file() and find_treeqsm_models(), however, are NOT reused:
#  both hard-code exactly two heights (0.0 m and 1.3 m via that file's own
#  H_BASE_M/H_DBH_M module constants) baked into their return values, with
#  no way to ask either one for a DIFFERENT, per-tree list of six
#  arbitrary heights without changing what they return - and that file is
#  finished/reviewed and must not be edited. The .mat-loading and
#  pd-combo-discovery LOOP inside find_treeqsm_models(), and the tab-
#  column field-file parsing logic inside parse_field_file(), are
#  therefore duplicated below in a more general form (arbitrary height
#  list, not two fixed heights). FLAGGED HERE for a later, separate
#  cleanup pass: if a third script ever needs this again, both files'
#  height-sampling logic should be factored into one shared, height-
#  parameterised helper instead of a third copy.
# =====================================================================

import csv
import glob
import hashlib
import os

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ensure_csv_dir()/ensure_all_plots_dir(): two CSVs (routed to csv/, flat)
# and two PNGs that pool across all 12 trees (routed to plots/all/) - see
# paths.py's own header for why "plots" as a bare string literal must
# never appear outside that one file.
from paths import ensure_csv_dir, ensure_all_plots_dir
import paths

# Shared visual style - see plot_style.py's own header. Reused so this
# chart's colors/sizes match every other chart in the project (and
# base_vs_dbh_compare.py specifically, which uses the SAME three family
# colors for the same three sources).
from plot_style import (
    FAMILY_GRADIENTS,
    AXIS_LABEL_FONTSIZE,
    LEGEND_FONTSIZE,
    PANEL_TITLE_FONTSIZE,
    DEFAULT_LINEWIDTH,
)

# Reuse (do not re-implement): the exact same taper.txt parser and trunk-
# radius function the real AdTree calibration pipeline uses - imported
# directly from tree_geom_utils.py, the same way base_vs_dbh_compare.py
# itself does (not a duplication - both scripts call the one real
# implementation).
from tree_geom_utils import parse_adqsm_taper_file, make_trunk_radius_func

# Reuse (do not re-implement): tree order + each tree's own primary
# (non-999) AdQSM variant, read from volume_results.csv itself - fully
# reusable as-is, no height assumption baked into it, so imported
# directly rather than duplicated.
from base_vs_dbh_compare import (
    tree_order_and_primary_variant,
    parse_field_file,
    TREEQSM_ARCHIVE_DIR,
    MATLAB_ARCHIVE_ROOT,
    DATA_ROOT,
    RESULTS_CSV,
)

# =====================  PARAMETERS  ===================================
# NOTE: this is a DATA file - Git ignores it, so on a new computer (e.g. the
# Mac) it must be copied by hand into scripts/python/ before running.
FIELD_FILE = str(paths.PYTHON_DIR / "Popis_Babice_VSE_09012025_ans_lev6_0.txt")

# Field survey's own level labels, in the file's own column order (level I
# = highest measured point, level VI = base, h = 0 m by definition) -
# height columns 7-12, diameter column PAIRS 13-24 (two perpendicular
# tape measurements per level - each level's diameter is their mean).
FIELD_LEVELS = ["I", "II", "III", "IV", "V", "VI"]
FIELD_HEIGHT_COLS = [7, 8, 9, 10, 11, 12]
FIELD_DIAM_COL_PAIRS = [(13, 14), (15, 16), (17, 18), (19, 20), (21, 22), (23, 24)]

H_DBH_M = 1.3   # breast height - the boundary of the "above base" error figure

# Fine height grid used ONLY for the continuous curves drawn in the
# small-multiples figure (Step 3c) - NOT used for any reported number;
# every number in the two CSVs comes from evaluating each method at the
# tree's own exact field heights (see compute_*_diam_at_heights() below).
CURVE_GRID_STEP_M = 0.1

FIGSIZE_GRID = (16, 11)     # 3x4 small-multiples figure
FIGSIZE_BAR = (10, 6)
PLOT_DPI = 150
BAR_WIDTH = 0.32

TAPER_CSV_NAME = "trunk_taper_vs_field.csv"
ERROR_CSV_NAME = "trunk_taper_error_summary.csv"
GRID_PLOT_NAME = "trunk_taper_vs_field.png"
ERROR_PLOT_NAME = "trunk_taper_error.png"
# =====================================================================


# =====================  EXTENSION PARAMETERS  ==========================
#  Per-pd2min breakdown (follow-up task). Nothing above this block was
#  changed - see the follow-up task's own Hard rule 5 ("the four current
#  outputs must come out byte-identical") and Step 1's checksum
#  requirement. Everything from here to the end of the file is additive:
#  new constants, new functions, and new code appended after the
#  original script's own last line.
# =====================================================================

# The four manual pd2min levels (Step 4a's x-axis / Step 3's CSV rows).
# The AUTO model is NEVER one of these - its own pd2min is tree-specific
# (see Background: "the AUTO model has a tree-specific pd2min and must
# not be pooled with the manual levels").
PD2MIN_LEVELS_MANUAL = [0.002, 0.005, 0.010, 0.020]

# Tolerance when matching a model's own (PatchDiam1, PatchDiam2Min,
# PatchDiam2Max) - read from rundata.inputs, never from a filename tag -
# against volume_results.csv's own pd1_m/pd2min_m/pd2max_m columns
# (Step 3's n_cylinders/volume lookup).
PD2MIN_MATCH_TOL_M = 1e-5

# SANITY CHECK 2 tolerance, in percentage points. The median of 5
# PER-MODEL err_above_base_pct values cannot equal the original script's
# per-tree value exactly: that value is mean-abs-error of the
# PER-HEIGHT MEDIAN diameter across models, while this is the MEDIAN,
# across models, of each model's OWN mean-abs-error. Median and
# mean-abs-error do not commute, so a small gap between the two is
# expected - this bounds how large a gap is still "the same underlying
# number, decomposed" rather than a sign of an actual bug.
SANITY2_TOL_PP = 1.0

# Step 1's own baseline: sha256 of each of the four ORIGINAL outputs,
# captured by running this script BEFORE this extension existed (see the
# follow-up task's Step 1). Hardcoded so every future run can self-check
# hard rule 5 ("byte-identical") without a separate baseline file.
BASELINE_SHA256 = {
    TAPER_CSV_NAME: "bae061692dbba82a2668cb2e474fba3984d576c19ed4a505090f90b3944fea1e",
    ERROR_CSV_NAME: "023cb098018e6be6dd40df7473bd11a39487687a7483717b4294506e521c9687",
    GRID_PLOT_NAME: "360c5899b892a56d65a638db990d64b7a6c5cf94f85b047d0a92ca6d7f42be6b",
    ERROR_PLOT_NAME: "91faf381efd42e0c17a20787fed499397f5e4330302240db2396bbad999dadc0",
}

PER_MODEL_CSV_NAME = "trunk_taper_per_pd2min.csv"
ERROR_PER_PD2MIN_PLOT_NAME = "trunk_taper_error_per_pd2min.png"
COST_PLOT_NAME = "pd2min_volume_vs_cylinders.png"
FIGSIZE_PD2MIN_ERR = (11, 7)
FIGSIZE_COST = (10, 9)
# =====================================================================


def parse_field_file_all_levels(path):
    """Parse Popis_Babice_VSE_09012025_ans_lev6_0.txt into ALL SIX levels
    per tree (base_vs_dbh_compare.py's parse_field_file() only keeps level
    VI and the dedicated 1.3 m column - see this file's header for why
    that function isn't reused here). Returns
    {tree: [(level_label, height_m, diam_cm), ...]}, sorted ascending by
    height (level VI / h=0 first). tree_name is built as "B21_S%02d" %
    tree_no, same convention as base_vs_dbh_compare.py, for the same
    reason (confirmed by cross-checking against known field values in an
    earlier, separate diagnostic)."""
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    data = {}
    for line in lines[2:]:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 25 or not parts[2].strip():
            continue
        tree_no = int(parts[2])
        tree = "B21_S%02d" % tree_no
        rows = []
        for level, hcol, (d1col, d2col) in zip(FIELD_LEVELS, FIELD_HEIGHT_COLS, FIELD_DIAM_COL_PAIRS):
            height_m = float(parts[hcol])
            diam_cm = (float(parts[d1col]) + float(parts[d2col])) / 2.0
            rows.append((level, height_m, diam_cm))
        rows.sort(key=lambda r: r[1])
        data[tree] = rows
    return data


def adqsm_diam_at_heights(tree, variant, heights_m):
    """AdQSM/AdTree-calibrated diameter, in cm, at every height in
    `heights_m`, from the REAL make_trunk_radius_func() built on that
    tree's own taper.txt for `variant`. Returns None if the file doesn't
    exist, else a list of diameters (cm), same order/length as
    `heights_m`."""
    path = os.path.join(DATA_ROOT, tree, variant, "taper.txt")
    if not os.path.exists(path):
        return None
    heights, diameters = parse_adqsm_taper_file(path)
    radius_func = make_trunk_radius_func(heights, diameters)
    return [2.0 * radius_func(h) * 100.0 for h in heights_m]


def _load_treeqsm_trunk_models(tree):
    """Discover one simplified_*_ri0.mat per distinct (PatchDiam1,
    PatchDiam2Min, PatchDiam2Max) combination for `tree` (SAME selection
    rule as base_vs_dbh_compare.py's find_treeqsm_models() - see this
    file's header for why the loading loop itself is duplicated rather
    than imported). Returns a list of dicts, one per model:
    {"pd_combo": (pd1,pd2min,pd2max), "seg_lo": array, "seg_hi": array,
     "radius": array} - seg_lo/seg_hi are each TRUNK cylinder's own
    height-above-ground range (min start-z among trunk cylinders
    subtracted, exactly like base_vs_dbh_compare.py does), radius is that
    cylinder's own radius [m]."""
    archive = TREEQSM_ARCHIVE_DIR.get(tree)
    if archive is None:
        return []
    folder = os.path.join(MATLAB_ARCHIVE_ROOT, archive)
    pattern = os.path.join(folder, "simplified_%s_*_ri0.mat" % tree)
    candidates = sorted(glob.glob(pattern))

    by_pd_combo = {}
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
        by_pd_combo.setdefault(key, qsm)

    models = []
    for pd_combo, qsm in sorted(by_pd_combo.items()):
        cyl = qsm["cylinder"]
        order = np.asarray(cyl["BranchOrder"])
        trunk = order == 0
        start_z = np.asarray(cyl["start"])[:, 2][trunk]
        axis_z = np.asarray(cyl["axis"])[:, 2][trunk]
        length = np.asarray(cyl["length"])[trunk]
        radius = np.asarray(cyl["radius"])[trunk]

        z0 = float(start_z.min())
        lo_h = start_z - z0
        hi_h = start_z + length * axis_z - z0
        seg_lo = np.minimum(lo_h, hi_h)
        seg_hi = np.maximum(lo_h, hi_h)

        models.append({"pd_combo": pd_combo, "seg_lo": seg_lo, "seg_hi": seg_hi, "radius": radius})
    return models


def _radius_at_height(model, h):
    """The radius [m] of whichever trunk cylinder's own height range
    contains `h`, for one model dict from _load_treeqsm_trunk_models()."""
    mask = (model["seg_lo"] <= h) & (h <= model["seg_hi"])
    idx = np.where(mask)[0]
    if len(idx) == 0:
        mid = (model["seg_lo"] + model["seg_hi"]) / 2.0
        idx = [int(np.argmin(np.abs(mid - h)))]
    return float(model["radius"][idx[0]])


def treeqsm_diam_at_heights(models, heights_m):
    """For a list of models (from _load_treeqsm_trunk_models()), return,
    for EACH height in heights_m: (median_cm, min_cm, max_cm) across all
    models - the exact "median across the five models of the radius of
    the containing cylinder" rule this task specifies, generalised from
    base_vs_dbh_compare.py's fixed two-height version to any height list."""
    result = []
    for h in heights_m:
        diams_cm = [2.0 * _radius_at_height(m, h) * 100.0 for m in models]
        result.append((float(np.median(diams_cm)), float(min(diams_cm)), float(max(diams_cm))))
    return result


def treeqsm_curve_cm(models, height_grid_m):
    """Smooth TreeQSM taper curve for the FIGURE ONLY (Step 3c) - not used
    for any reported number. For each model, builds a height/diameter
    curve from that model's own trunk cylinders (midpoint height, own
    diameter), linearly interpolated (clamped outside the model's own
    range, same convention as make_trunk_radius_func()), then returns the
    MEDIAN curve across models at `height_grid_m` - matching the
    "TreeQSM (median)" convention base_vs_dbh_compare.py's own figure
    already uses."""
    curves = []
    for m in models:
        mid_h = (m["seg_lo"] + m["seg_hi"]) / 2.0
        order = np.argsort(mid_h)
        mid_h_sorted = mid_h[order]
        diam_cm_sorted = 2.0 * m["radius"][order] * 100.0
        curves.append(np.interp(height_grid_m, mid_h_sorted, diam_cm_sorted))
    return np.median(np.array(curves), axis=0)


def mean_abs_pct_error(diffs_pct):
    return float(np.mean(np.abs(diffs_pct))) if diffs_pct else None


# =====================  EXTENSION FUNCTIONS  ============================

def _load_treeqsm_trunk_models_with_mode(tree):
    """SAME selection/dedup rule and SAME trunk-cylinder geometry as
    _load_treeqsm_trunk_models() above - a PARALLEL read of the same
    files, not a modification of that function (Hard rule 5: the four
    original outputs must stay byte-identical, so that function is never
    touched - this duplicates its short loading loop rather than risk
    it). Additionally keeps, per model:
      "mode"    - "auto"/"manual", from the run tag's own "aut"/"man"
                  prefix (the archive's own naming convention - see
                  import_matlab_results.py's identical "aut"/"man" split)
                  - cross-checked below against that run's own
                  params_<tree>_<run>.csv "mode" column, so this is
                  confirmed, not assumed.
    pd1/pd2min/pd2max themselves are NEVER re-derived from a filename -
    "pd_combo" (as in the original function) always comes straight from
    rundata.inputs, exactly like the original."""
    archive = TREEQSM_ARCHIVE_DIR.get(tree)
    if archive is None:
        return []
    folder = os.path.join(MATLAB_ARCHIVE_ROOT, archive)
    pattern = os.path.join(folder, "simplified_%s_*_ri0.mat" % tree)
    candidates = sorted(glob.glob(pattern))

    by_pd_combo = {}
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

    models = []
    for pd_combo, (path, qsm) in sorted(by_pd_combo.items()):
        cyl = qsm["cylinder"]
        order = np.asarray(cyl["BranchOrder"])
        trunk = order == 0
        start_z = np.asarray(cyl["start"])[:, 2][trunk]
        axis_z = np.asarray(cyl["axis"])[:, 2][trunk]
        length = np.asarray(cyl["length"])[trunk]
        radius = np.asarray(cyl["radius"])[trunk]

        z0 = float(start_z.min())
        lo_h = start_z - z0
        hi_h = start_z + length * axis_z - z0
        seg_lo = np.minimum(lo_h, hi_h)
        seg_hi = np.maximum(lo_h, hi_h)

        run_tag = os.path.basename(path)[len("simplified_%s_" % tree):-len(".mat")]
        mode = "auto" if run_tag.startswith("aut") else "manual"
        params_path = os.path.join(folder, "params_%s_%s.csv" % (tree, run_tag))
        if os.path.exists(params_path):
            with open(params_path, "r", encoding="utf-8-sig", newline="") as pf:
                prow = next(csv.DictReader(pf), None)
            csv_mode = (prow.get("mode") or "").strip() if prow else ""
            if csv_mode and csv_mode != mode:
                print("  WARNING: mode mismatch for %s: filename says %s, params csv says %s"
                      % (run_tag, mode, csv_mode))

        models.append({"pd_combo": pd_combo, "seg_lo": seg_lo, "seg_hi": seg_hi,
                        "radius": radius, "mode": mode})
    return models


def _diam_cm_at_heights_for_model(model, heights_m):
    """diameter [cm] at every height in heights_m, for ONE model dict -
    reuses _radius_at_height() exactly as-is (same "containing cylinder,
    else nearest midpoint" rule as the median-across-models version
    above), just not yet reduced across models."""
    return [2.0 * _radius_at_height(model, h) * 100.0 for h in heights_m]


def _pd2min_level_of_row(row):
    """Bucket one trunk_taper_per_pd2min.csv row into one of
    PD2MIN_LEVELS_MANUAL (nearest value - the float read from
    rundata.inputs is already exact to the intended level, this only
    guards float noise), or "AUTO" for an auto-mode row (its own pd2min
    is tree-specific - must never be pooled with the manual levels, per
    the follow-up task's own Background)."""
    if row["mode"] == "auto":
        return "AUTO"
    return min(PD2MIN_LEVELS_MANUAL, key=lambda L: abs(L - row["pd2min_m"]))


def _to_float_or_none(x):
    """Parse a raw (string) volume_results.csv cell to float, None for
    blank/missing - same convention as compare_volumes.py's to_float(),
    duplicated here as a two-line helper rather than adding a new
    cross-file import for one function."""
    if x is None or (isinstance(x, str) and not x.strip()):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _match_volume_results_row(all_rows, tree, mode, pd_combo):
    """Find volume_results.csv's own "TreeQSM mine (..., Simplified (no
    islands))" row at sr=0.010/ri=1 whose own (pd1_m, pd2min_m, pd2max_m)
    match `pd_combo` NUMERICALLY (never by reconstructing/parsing a
    filename pd-tag - pd tags are not unique across trees, exactly the
    trap a previous diagnostic on this project already found and had to
    work around). Returns the raw CSV row dict, or None if no sr010/ri1
    export of this pd-combo exists in the archive (see Step 3's own
    report for how often that happens - never silently substituted with
    a different sr/ri)."""
    pd1, pd2min, pd2max = pd_combo
    for r in all_rows:
        if r["tree"] != tree or (r.get("mode") or "").strip() != mode:
            continue
        if "sr010_ri1" not in r["method"] or "Simplified (no islands)" not in r["method"]:
            continue
        rp1 = _to_float_or_none(r.get("pd1_m"))
        rp2 = _to_float_or_none(r.get("pd2min_m"))
        rp3 = _to_float_or_none(r.get("pd2max_m"))
        if rp1 is None or rp2 is None or rp3 is None:
            continue
        if (abs(rp1 - pd1) <= PD2MIN_MATCH_TOL_M and abs(rp2 - pd2min) <= PD2MIN_MATCH_TOL_M
                and abs(rp3 - pd2max) <= PD2MIN_MATCH_TOL_M):
            return r
    return None


# =========================  RUN  =====================================
if __name__ == "__main__":
    tree_order, primary_variant, all_rows = tree_order_and_primary_variant(RESULTS_CSV)
    print("Trees, in order:", tree_order)

    field_all_levels = parse_field_file_all_levels(FIELD_FILE)

    # ---- STEP 1 sanity check: loaders imported vs re-written -----------
    print()
    print("=" * 90)
    print("SANITY CHECK 1 - loader provenance")
    print("=" * 90)
    print("tree_order_and_primary_variant() : IMPORTED from base_vs_dbh_compare.py (reused as-is)")
    print("AdQSM taper (parse_adqsm_taper_file / make_trunk_radius_func) : IMPORTED from tree_geom_utils.py")
    print("parse_field_file_all_levels() : RE-WRITTEN (base_vs_dbh_compare.py's parse_field_file() only")
    print("  returns level VI + the 1.3 m column, hard-coded - see this file's header note)")
    print("TreeQSM .mat loader (_load_treeqsm_trunk_models) : RE-WRITTEN (find_treeqsm_models() hard-codes")
    print("  exactly two heights (0.0/1.3 m) into its return value - see this file's header note)")
    print("  -> flagged for a later, separate cleanup pass, as instructed.")
    print()

    taper_rows = []      # for csv/trunk_taper_vs_field.csv
    error_rows = []       # for csv/trunk_taper_error_summary.csv
    treeqsm_models_by_tree = {}
    max_field_height_by_tree = {}

    n_models_report = {}
    for tree in tree_order:
        variant = primary_variant.get(tree)
        levels = field_all_levels.get(tree)
        if levels is None:
            print("MISSING: field data for %s" % tree)
            continue

        heights_m = [lvl[1] for lvl in levels]
        field_diams_cm = [lvl[2] for lvl in levels]
        max_field_height_by_tree[tree] = max(heights_m)

        adqsm_diams = adqsm_diam_at_heights(tree, variant, heights_m) if variant else None
        if adqsm_diams is None:
            print("MISSING: AdQSM taper.txt for %s / %s" % (tree, variant))
            adqsm_diams = [None] * len(heights_m)

        models = _load_treeqsm_trunk_models(tree)
        treeqsm_models_by_tree[tree] = models
        n_models_report[tree] = len(models)
        treeqsm_stats = treeqsm_diam_at_heights(models, heights_m) if models else \
            [(None, None, None)] * len(heights_m)

        adqsm_diffs_all, adqsm_diffs_above = [], []
        treeqsm_diffs_all, treeqsm_diffs_above = [], []

        for (level, h, field_d), adqsm_d, (tq_med, tq_min, tq_max) in zip(levels, adqsm_diams, treeqsm_stats):
            adqsm_pct = (100.0 * (adqsm_d - field_d) / field_d) if (adqsm_d is not None and field_d) else None
            treeqsm_pct = (100.0 * (tq_med - field_d) / field_d) if (tq_med is not None and field_d) else None

            taper_rows.append({
                "tree": tree, "level": level, "height_m": h,
                "field_d_cm": field_d, "adqsm_d_cm": adqsm_d,
                "treeqsm_d_cm_med": tq_med, "treeqsm_d_cm_min": tq_min, "treeqsm_d_cm_max": tq_max,
                "adqsm_diff_pct": adqsm_pct, "treeqsm_diff_pct": treeqsm_pct,
            })

            if adqsm_pct is not None:
                adqsm_diffs_all.append(adqsm_pct)
                if h > 0.0:
                    adqsm_diffs_above.append(adqsm_pct)
            if treeqsm_pct is not None:
                treeqsm_diffs_all.append(treeqsm_pct)
                if h > 0.0:
                    treeqsm_diffs_above.append(treeqsm_pct)

        error_rows.append({
            "tree": tree, "n_heights": len(heights_m), "max_field_height_m": max(heights_m),
            "adqsm_err_all_pct": mean_abs_pct_error(adqsm_diffs_all),
            "adqsm_err_above_base_pct": mean_abs_pct_error(adqsm_diffs_above),
            "treeqsm_err_all_pct": mean_abs_pct_error(treeqsm_diffs_all),
            "treeqsm_err_above_base_pct": mean_abs_pct_error(treeqsm_diffs_above),
        })

    print("TreeQSM models found per tree:", n_models_report, "(expected 5 for all 12)")
    print()

    # ---- STEP 3e: print both tables to stdout ---------------------------
    taper_cols = ["tree", "level", "height_m", "field_d_cm", "adqsm_d_cm",
                  "treeqsm_d_cm_med", "treeqsm_d_cm_min", "treeqsm_d_cm_max",
                  "adqsm_diff_pct", "treeqsm_diff_pct"]
    error_cols = ["tree", "n_heights", "max_field_height_m",
                  "adqsm_err_all_pct", "adqsm_err_above_base_pct",
                  "treeqsm_err_all_pct", "treeqsm_err_above_base_pct"]

    def fmt(v):
        if v is None:
            return "-"
        if isinstance(v, int):
            return str(v)
        return "%.4f" % v

    print("=" * 150)
    print("csv/trunk_taper_vs_field.csv (long format, one row per tree x level)")
    print("=" * 150)
    print(("%-9s %-6s" + " %10s" * (len(taper_cols) - 2)) % tuple(taper_cols))
    for row in taper_rows:
        print(("%-9s %-6s" + " %10s" * (len(taper_cols) - 2))
              % tuple([row["tree"], row["level"]] + [fmt(row[c]) for c in taper_cols[2:]]))
    print()

    print("=" * 150)
    print("csv/trunk_taper_error_summary.csv (one row per tree)")
    print("=" * 150)
    print(("%-9s" + " %20s" * (len(error_cols) - 1)) % tuple(error_cols))
    for row in error_rows:
        print(("%-9s" + " %20s" * (len(error_cols) - 1))
              % tuple([row["tree"]] + [fmt(row[c]) for c in error_cols[1:]]))
    print()

    # ---- STEP 3a/3b: write the two CSVs ----------------------------------
    taper_csv_path = os.path.join(ensure_csv_dir(), TAPER_CSV_NAME)
    with open(taper_csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=taper_cols)
        writer.writeheader()
        for row in taper_rows:
            writer.writerow({c: ("" if row[c] is None else row[c]) for c in taper_cols})
    print("Saved:", taper_csv_path)

    error_csv_path = os.path.join(ensure_csv_dir(), ERROR_CSV_NAME)
    with open(error_csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=error_cols)
        writer.writeheader()
        for row in error_rows:
            writer.writerow({c: ("" if row[c] is None else row[c]) for c in error_cols})
    print("Saved:", error_csv_path)
    print()

    # ---- STEP 3c: 3x4 small multiples -----------------------------------
    adqsm_color = FAMILY_GRADIENTS["AdQSM"][-1]
    treeqsm_color = FAMILY_GRADIENTS["TreeQSM"][-1]
    field_color = FAMILY_GRADIENTS["Reference"][-1]

    fig, axes = plt.subplots(3, 4, figsize=FIGSIZE_GRID, sharex=False, sharey=False)
    axes_flat = axes.flatten()

    for i, tree in enumerate(tree_order):
        ax = axes_flat[i]
        levels = field_all_levels.get(tree)
        if levels is None:
            continue
        heights_m = [lvl[1] for lvl in levels]
        field_diams_cm = [lvl[2] for lvl in levels]
        max_h = max_field_height_by_tree[tree]

        # shade the buttress zone (below breast height) so it's visually
        # distinct from the round part of the stem
        ax.axvspan(0.0, H_DBH_M, color="#f0f0f0", zorder=0)

        grid = np.arange(0.0, max_h + 1e-9, CURVE_GRID_STEP_M)
        variant = primary_variant.get(tree)
        adqsm_curve = adqsm_diam_at_heights(tree, variant, grid) if variant else None
        if adqsm_curve is not None:
            ax.plot(grid, adqsm_curve, color=adqsm_color, linewidth=DEFAULT_LINEWIDTH * 2, label="adqsm_taper")

        models = treeqsm_models_by_tree.get(tree)
        if models:
            treeqsm_curve = treeqsm_curve_cm(models, grid)
            ax.plot(grid, treeqsm_curve, color=treeqsm_color, linewidth=DEFAULT_LINEWIDTH * 2, label="treeqsm")

        ax.scatter(heights_m, field_diams_cm, color=field_color, zorder=5, s=30, label="field")

        ax.set_title(tree, fontsize=PANEL_TITLE_FONTSIZE)
        ax.set_xlabel("Height [m]", fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel("Diameter [cm]", fontsize=AXIS_LABEL_FONTSIZE)

    for j in range(len(tree_order), len(axes_flat)):
        axes_flat[j].axis("off")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=LEGEND_FONTSIZE, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    grid_plot_path = os.path.join(ensure_all_plots_dir(), GRID_PLOT_NAME)
    fig.savefig(grid_plot_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", grid_plot_path)

    # ---- STEP 3d: error bar chart (err_above_base only) ------------------
    trees_x = [row["tree"] for row in error_rows]
    x = np.arange(len(trees_x))
    adqsm_err_vals = [row["adqsm_err_above_base_pct"] for row in error_rows]
    treeqsm_err_vals = [row["treeqsm_err_above_base_pct"] for row in error_rows]

    fig2, ax2 = plt.subplots(figsize=FIGSIZE_BAR)
    ax2.bar(x - BAR_WIDTH / 2, adqsm_err_vals, BAR_WIDTH, label="adqsm_taper", color=adqsm_color)
    ax2.bar(x + BAR_WIDTH / 2, treeqsm_err_vals, BAR_WIDTH, label="treeqsm", color=treeqsm_color)
    ax2.set_xticks(x)
    ax2.set_xticklabels(trees_x, fontsize=AXIS_LABEL_FONTSIZE, rotation=45, ha="right")
    ax2.set_ylabel("Mean absolute %% error\n(5 heights above base, h=0 m excluded)", fontsize=AXIS_LABEL_FONTSIZE)
    ax2.set_xlabel("Tree", fontsize=AXIS_LABEL_FONTSIZE)
    ax2.legend(fontsize=LEGEND_FONTSIZE)
    fig2.tight_layout()
    error_plot_path = os.path.join(ensure_all_plots_dir(), ERROR_PLOT_NAME)
    fig2.savefig(error_plot_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig2)
    print("Saved:", error_plot_path)
    print()

    # =====================  STEP 4 - sanity checks  ====================
    print("=" * 90)
    print("SANITY CHECKS")
    print("=" * 90)

    # 4a. Uses parse_field_file() (IMPORTED from base_vs_dbh_compare.py,
    # not re-derived) for the ratio, because "diameter at 1.3 m (cm)" is
    # its OWN dedicated column [26] in the source file - it is NOT
    # generally one of the six level heights (e.g. B21_S01's nearest
    # level to 1.3 m is actually at h=1.16 m), so approximating it from
    # the nearest level would silently use the wrong denominator. The
    # base value itself (level VI) is cross-checked against this script's
    # own parse_field_file_all_levels() below, as an extra confirmation
    # the two parses agree.
    field_base_dbh = parse_field_file(FIELD_FILE)
    b07_base = field_base_dbh["B21_S07"][0]
    b12_base = field_base_dbh["B21_S12"][0]
    ratios = {t: base / dbh for t, (base, dbh) in field_base_dbh.items()}
    min_t = min(ratios, key=ratios.get)
    max_t = max(ratios, key=ratios.get)
    # cross-check: level VI from the all-levels parser must match the
    # dedicated base/dbh parser's own base value, for every tree
    level_vi_mismatch = [t for t, rows in field_all_levels.items()
                          if abs(next(lvl[2] for lvl in rows if lvl[1] == 0.0) - field_base_dbh[t][0]) > 1e-6]
    if level_vi_mismatch:
        print("  WARNING: level-VI cross-check mismatch for:", level_vi_mismatch)
    print("4a. Field base diam: B21_S07=%.2f cm  B21_S12=%.2f cm | ratio min=%.4f (%s)  max=%.4f (%s)"
          % (b07_base, b12_base, ratios[min_t], min_t, ratios[max_t], max_t))
    check_4a = (abs(b07_base - 114.25) < 0.02 and abs(b12_base - 40.25) < 0.02
                and abs(ratios["B21_S12"] - 1.3417) < 0.001 and abs(ratios["B21_S09"] - 2.0764) < 0.001)
    print("    expected: 114.25/40.25 cm, ratio min 1.3417 (B21_S12) max 2.0764 (B21_S09) ->",
          "PASS" if check_4a else "FAIL")

    # 4b
    def diffs_for(tree):
        rows_for_tree = [r for r in taper_rows if r["tree"] == tree]
        rows_for_tree.sort(key=lambda r: r["height_m"])
        return [r["adqsm_diff_pct"] for r in rows_for_tree]

    b04_expected = [-49.5, -13.1, -9.6, -8.2, -10.8, -12.0]
    b16_expected = [-33.6, -12.5, -6.4, -2.9, 0.2, 3.6]
    b04_actual = diffs_for("B21_S04")
    b16_actual = diffs_for("B21_S16")
    print("4b. B21_S04 adqsm_diff_pct: actual  =", ["%.1f" % v for v in b04_actual])
    print("             expected =", b04_expected)
    print("    B21_S16 adqsm_diff_pct: actual  =", ["%.1f" % v for v in b16_actual])
    print("             expected =", b16_expected)
    check_4b = (all(abs(a - e) < 1.0 for a, e in zip(b04_actual, b04_expected))
                and all(abs(a - e) < 1.0 for a, e in zip(b16_actual, b16_expected)))
    print("    ->", "PASS" if check_4b else "FAIL")

    # 4c
    max_diff_4c = 0.0
    for tree in tree_order:
        variant = primary_variant.get(tree)
        if not variant:
            continue
        row_match = [r for r in all_rows if r["tree"] == tree and r["calmethod"]
                     and r["adqsm_variant"].strip().zfill(3) == variant]
        if not row_match:
            continue
        csv_dbh_m = float(row_match[0]["dbh_m"])
        my_row = next((r for r in taper_rows if r["tree"] == tree and abs(r["height_m"] - H_DBH_M) < 1e-6), None)
        if my_row is None or my_row["adqsm_d_cm"] is None:
            continue
        d = abs(csv_dbh_m - my_row["adqsm_d_cm"] / 100.0)
        max_diff_4c = max(max_diff_4c, d)
    print("4c. adqsm_d_cm @ 1.3m vs volume_results.csv dbh_m: max abs diff = %.6f m" % max_diff_4c)
    check_4c = max_diff_4c < 0.0001
    print("    expected: < 0.1 mm ->", "PASS" if check_4c else "FAIL")

    # 4d
    max_diff_4d = 0.0
    for tree in tree_order:
        variant = primary_variant.get(tree)
        if not variant:
            continue
        row_match = [r for r in all_rows if r["tree"] == tree and r["calmethod"]
                     and r["adqsm_variant"].strip().zfill(3) == variant]
        if not row_match:
            continue
        csv_dbh_m = float(row_match[0]["dbh_m"])
        my_row = next((r for r in taper_rows if r["tree"] == tree and abs(r["height_m"] - H_DBH_M) < 1e-6), None)
        if my_row is None or my_row["treeqsm_d_cm_med"] is None:
            continue
        d = abs(csv_dbh_m - my_row["treeqsm_d_cm_med"] / 100.0)
        max_diff_4d = max(max_diff_4d, d)
    print("4d. treeqsm_d_cm_med @ 1.3m vs volume_results.csv dbh_m: max abs diff = %.6f m "
          "(expected ~DBHcyl vs DBHqsm gap, already explained in base_vs_dbh_compare.py)" % max_diff_4d)
    check_4d = max_diff_4d < 0.02
    print("    ->", "PASS (within known DBHcyl/DBHqsm gap)" if check_4d else "FAIL (larger than expected)")

    # 4e - headline result
    adqsm_above_meds = [row["adqsm_err_above_base_pct"] for row in error_rows
                         if row["adqsm_err_above_base_pct"] is not None]
    treeqsm_above_meds = [row["treeqsm_err_above_base_pct"] for row in error_rows
                           if row["treeqsm_err_above_base_pct"] is not None]
    adqsm_med_above = float(np.median(adqsm_above_meds))
    treeqsm_med_above = float(np.median(treeqsm_above_meds))
    n_treeqsm_lower_above = sum(1 for row in error_rows
                                 if row["treeqsm_err_above_base_pct"] is not None
                                 and row["adqsm_err_above_base_pct"] is not None
                                 and row["treeqsm_err_above_base_pct"] < row["adqsm_err_above_base_pct"])
    print()
    print("4e. mean absolute %% error, five heights above the base, median across 12 trees:")
    print("      adqsm_taper : %.4f %%" % adqsm_med_above)
    print("      treeqsm     : %.4f %%" % treeqsm_med_above)
    print("    trees where treeqsm has the lower error : %d of %d" % (n_treeqsm_lower_above, len(error_rows)))
    print()

    # =====================  STEP 5 - summary  ===========================
    n_heights_per_tree = {t: len(field_all_levels[t]) for t in tree_order if t in field_all_levels}
    min_max_h = [max_field_height_by_tree[t] for t in tree_order if t in max_field_height_by_tree]

    adqsm_all_meds_list = [row["adqsm_err_all_pct"] for row in error_rows if row["adqsm_err_all_pct"] is not None]
    treeqsm_all_meds_list = [row["treeqsm_err_all_pct"] for row in error_rows if row["treeqsm_err_all_pct"] is not None]
    n_treeqsm_lower_all = sum(1 for row in error_rows
                               if row["treeqsm_err_all_pct"] is not None and row["adqsm_err_all_pct"] is not None
                               and row["treeqsm_err_all_pct"] < row["adqsm_err_all_pct"])

    print("=" * 90)
    print("=== TRUNK TAPER VS FIELD SUMMARY ===")
    print("=" * 90)
    print("Loaders : tree_order_and_primary_variant() imported from base_vs_dbh_compare.py;")
    print("          field-file and TreeQSM .mat loaders re-written (see header note)")
    print("TreeQSM models per tree :", n_models_report, "(expected 5)")
    print("Field heights per tree  :", n_heights_per_tree, "(expected 6)")
    print("Comparison valid up to  : %.2f m (min across trees) .. %.2f m (max)"
          % (min(min_max_h), max(min_max_h)))
    print()
    print("err_above_base, median across 12 trees : adqsm %.4f %%   treeqsm %.4f %%"
          % (adqsm_med_above, treeqsm_med_above))
    print("err_all       , median across 12 trees : adqsm %.4f %%   treeqsm %.4f %%"
          % (float(np.median(adqsm_all_meds_list)), float(np.median(treeqsm_all_meds_list))))
    print("trees where treeqsm lower (above base) : %d of %d" % (n_treeqsm_lower_above, len(error_rows)))
    print("trees where treeqsm lower (all heights): %d of %d" % (n_treeqsm_lower_all, len(error_rows)))
    print()
    print("Outputs : %s , %s ," % (taper_csv_path, error_csv_path))
    print("          %s , %s" % (grid_plot_path, error_plot_path))
    print()
    print("Sanity checks: 4a %s  4b %s  4c %s  4d %s"
          % ("PASS" if check_4a else "FAIL", "PASS" if check_4b else "FAIL",
             "PASS" if check_4c else "FAIL", "PASS" if check_4d else "FAIL"))

    # =====================================================================
    #  EXTENSION - per-pd2min breakdown (follow-up task, Steps 1-5).
    #  Nothing above this line was changed by this extension - this block
    #  only READS variables the original script already computed above
    #  (tree_order, field_all_levels, error_rows, all_rows, field_base_dbh)
    #  and adds new outputs alongside the original four.
    # =====================================================================

    # ---- SANITY CHECK 1: the four original outputs, byte-identical -------
    print()
    print("=" * 90)
    print("SANITY CHECK 1 - original outputs unchanged (Hard rule 5)")
    print("=" * 90)
    check_1 = True
    for name, path in ((TAPER_CSV_NAME, taper_csv_path), (ERROR_CSV_NAME, error_csv_path),
                        (GRID_PLOT_NAME, grid_plot_path), (ERROR_PLOT_NAME, error_plot_path)):
        actual_hash = hashlib.sha256(open(path, "rb").read()).hexdigest()
        expected_hash = BASELINE_SHA256[name]
        ok = (actual_hash == expected_hash)
        check_1 = check_1 and ok
        print("  %-32s : %s" % (name, "PASS" if ok else "FAIL (hash changed since baseline)"))
    print()

    print("=" * 90)
    print("EXTENSION: per-pd2min breakdown")
    print("=" * 90)

    per_model_rows = []
    sanity2_rows = []           # (tree, median_of_5, existing_value, diff)
    missing_volume_match = []   # (tree, mode, pd2min) - no sr010/ri1 CSV row found

    for tree in tree_order:
        levels = field_all_levels.get(tree)
        if levels is None:
            continue
        heights_m = [lvl[1] for lvl in levels]
        dbh_cm = field_base_dbh.get(tree, (None, None))[1]

        models_meta = _load_treeqsm_trunk_models_with_mode(tree)
        per_model_err_above = []

        for m in models_meta:
            diam_cm = _diam_cm_at_heights_for_model(m, heights_m)
            diffs_pct, diffs_above = [], []
            for (level, h, field_d), dcm in zip(levels, diam_cm):
                pct = (100.0 * (dcm - field_d) / field_d) if field_d else None
                diffs_pct.append(pct)
                if h > 0.0 and pct is not None:
                    diffs_above.append(pct)
            err_above = mean_abs_pct_error(diffs_above)
            err_all = mean_abs_pct_error([d for d in diffs_pct if d is not None])
            per_model_err_above.append(err_above)

            pd1, pd2min, pd2max = m["pd_combo"]
            vol_row = _match_volume_results_row(all_rows, tree, m["mode"], m["pd_combo"])
            if vol_row is None:
                missing_volume_match.append((tree, m["mode"], round(pd2min, 4)))
                n_cyl = trunk_m3 = branch_m3 = total_m3 = None
            else:
                n_cyl = _to_float_or_none(vol_row.get("n_cylinders"))
                trunk_m3 = _to_float_or_none(vol_row.get("trunk_m3"))
                branch_m3 = _to_float_or_none(vol_row.get("branch_m3"))
                total_m3 = _to_float_or_none(vol_row.get("total_m3"))

            row = {"tree": tree, "field_dbh_cm": dbh_cm, "mode": m["mode"],
                   "pd1_m": pd1, "pd2min_m": pd2min, "pd2max_m": pd2max,
                   "n_cylinders": n_cyl, "trunk_m3": trunk_m3, "branch_m3": branch_m3, "total_m3": total_m3,
                   "err_above_base_pct": err_above, "err_all_pct": err_all}
            for i, d in enumerate(diffs_pct, start=1):
                row["diff_h%d_pct" % i] = d
            per_model_rows.append(row)

        if per_model_err_above and all(v is not None for v in per_model_err_above):
            med = float(np.median(per_model_err_above))
            existing = next((r["treeqsm_err_above_base_pct"] for r in error_rows if r["tree"] == tree), None)
            if existing is not None:
                sanity2_rows.append((tree, med, existing, abs(med - existing)))

    models_per_tree_ext = {t: sum(1 for r in per_model_rows if r["tree"] == t) for t in tree_order}
    print("Models per tree (extension loader) :", models_per_tree_ext, "(expected 5)")
    print()

    # ---- SANITY CHECK 2 ----------------------------------------------------
    print("SANITY CHECK 2 - median(err_above_base_pct across the 5 models) vs. the ORIGINAL")
    print("  per-tree value in trunk_taper_error_summary.csv (median-then-error vs.")
    print("  error-then-median do not commute - a small gap is expected, see module header)")
    for tree, med, existing, diff in sanity2_rows:
        print("  %-9s median_of_5_models=%.4f  existing=%.4f  diff=%.4f pp" % (tree, med, existing, diff))
    max_diff_2 = max((d for _, _, _, d in sanity2_rows), default=None)
    check_2 = (max_diff_2 is not None and max_diff_2 <= SANITY2_TOL_PP)
    print("  max diff = %.4f pp (tolerance %.1f pp) -> %s"
          % (max_diff_2 if max_diff_2 is not None else -1.0, SANITY2_TOL_PP, "PASS" if check_2 else "FAIL"))
    print()

    if missing_volume_match:
        levels_missing = sorted({(mode, lvl) for _, mode, lvl in missing_volume_match})
        print("NOTE: no volume_results.csv row at sr=0.010/ri=1 'Simplified (no islands)' for:")
        for mode, lvl in levels_missing:
            n = sum(1 for t, mo, l in missing_volume_match if mo == mode and l == lvl)
            print("  mode=%-6s pd2min=%.4f m : missing for %d/%d trees" % (mode, lvl, n, len(tree_order)))
        print("  (n_cylinders/trunk_m3/branch_m3/total_m3 left BLANK in the new CSV for these rows -")
        print("   never substituted from a different sr/ri, per this task's own hard rules)")
        print()

    # ---- Step 3: write the new CSV ------------------------------------------
    per_model_cols = (["tree", "field_dbh_cm", "mode", "pd1_m", "pd2min_m", "pd2max_m",
                        "n_cylinders", "trunk_m3", "branch_m3", "total_m3",
                        "err_above_base_pct", "err_all_pct"]
                       + ["diff_h%d_pct" % i for i in range(1, 7)])

    def _dbh_sort_key(row):
        dbh = row["field_dbh_cm"]
        return (dbh if dbh is not None else 1e9, row["tree"], row["pd2min_m"] or 0.0)

    per_model_rows_sorted = sorted(per_model_rows, key=_dbh_sort_key)

    per_model_csv_path = os.path.join(ensure_csv_dir(), PER_MODEL_CSV_NAME)
    with open(per_model_csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=per_model_cols)
        writer.writeheader()
        for row in per_model_rows_sorted:
            writer.writerow({c: ("" if row.get(c) is None else row[c]) for c in per_model_cols})
    print("Saved:", per_model_csv_path)
    print()

    # DBH-ascending tree order for the two NEW figures (project convention -
    # method_evaluation.py's compute_tree_order_by_dbh() - kept SEPARATE
    # from `tree_order` above, which is CSV first-appearance order and is
    # used by every ORIGINAL output; never touched, per Hard rule 5).
    dbh_tree_order = sorted([t for t in tree_order if t in field_base_dbh],
                             key=lambda t: (field_base_dbh[t][1], t))
    dbh_tree_order += [t for t in tree_order if t not in field_base_dbh]

    dbh_cmap = plt.cm.viridis(np.linspace(0, 1, len(dbh_tree_order)))
    tree_color = {t: dbh_cmap[i] for i, t in enumerate(dbh_tree_order)}

    # ---- Step 4a: accuracy figure -------------------------------------------
    fig3, ax3 = plt.subplots(figsize=FIGSIZE_PD2MIN_ERR)
    per_level_values = {L: [] for L in PD2MIN_LEVELS_MANUAL}

    for tree in dbh_tree_order:
        rows_t = [r for r in per_model_rows if r["tree"] == tree]
        manual_pts = sorted(
            [(r["pd2min_m"], r["err_above_base_pct"]) for r in rows_t
             if r["mode"] == "manual" and r["err_above_base_pct"] is not None])
        if manual_pts:
            xs, ys = zip(*manual_pts)
            ax3.plot(xs, ys, marker="o", color=tree_color[tree], linewidth=DEFAULT_LINEWIDTH,
                      markersize=5, label=tree.replace("B21_", ""))
            for r in rows_t:
                if r["mode"] == "manual" and r["err_above_base_pct"] is not None:
                    per_level_values[_pd2min_level_of_row(r)].append(r["err_above_base_pct"])

        auto_row = next((r for r in rows_t if r["mode"] == "auto"), None)
        if auto_row is not None and auto_row["err_above_base_pct"] is not None:
            ax3.scatter([auto_row["pd2min_m"]], [auto_row["err_above_base_pct"]],
                        marker="*", s=110, color=tree_color[tree], zorder=5,
                        edgecolors="black", linewidths=0.5)

    median_xs = PD2MIN_LEVELS_MANUAL
    median_ys = [float(np.median(per_level_values[L])) if per_level_values[L] else None for L in median_xs]
    ax3.plot(median_xs, median_ys, color="black", linewidth=DEFAULT_LINEWIDTH * 3,
              marker="o", markersize=7, label="median (12 trees)", zorder=10)

    ax3.set_xscale("log")
    ax3.set_xticks(PD2MIN_LEVELS_MANUAL)
    ax3.set_xticklabels(["%.3f" % L for L in PD2MIN_LEVELS_MANUAL])
    ax3.set_xlabel("PatchDiam2Min [m]  (log scale; ★ = AUTO model, plotted at its own\n"
                    "tree-specific pd2min, unconnected - never pooled with the 4 manual levels)",
                    fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_ylabel("Mean absolute %% error\n(5 heights above base, h=0 m excluded)",
                    fontsize=AXIS_LABEL_FONTSIZE)
    ax3.set_title("Trunk taper error vs. field, by pd2min level", fontsize=PANEL_TITLE_FONTSIZE)
    ax3.grid(True, which="both", alpha=0.2)
    ax3.legend(fontsize=LEGEND_FONTSIZE - 2, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig3.tight_layout()
    err_pd2min_plot_path = os.path.join(ensure_all_plots_dir(), ERROR_PER_PD2MIN_PLOT_NAME)
    fig3.savefig(err_pd2min_plot_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig3)
    print("Saved:", err_pd2min_plot_path)

    # ---- Step 4b: cost figure -----------------------------------------------
    fig4, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=FIGSIZE_COST, sharex=True)

    for tree in dbh_tree_order:
        rows_by_level = {_pd2min_level_of_row(r): r for r in per_model_rows
                          if r["tree"] == tree and r["mode"] == "manual"}
        base = rows_by_level.get(0.002)
        if base is None or base["trunk_m3"] is None:
            continue   # no sr010/ri1 export at pd2min=0.002 - nothing to normalise against
                       # (not expected to happen - see the missing-match note above)

        xs, trunk_ratio, branch_ratio, cyl_ratio = [], [], [], []
        for L in PD2MIN_LEVELS_MANUAL:
            r = rows_by_level.get(L)
            if r is None or r["trunk_m3"] is None:
                continue   # this tree/level has no sr010/ri1 export - gap left in the line,
                           # not filled from a different sr/ri (see missing-match note above)
            xs.append(L)
            trunk_ratio.append(r["trunk_m3"] / base["trunk_m3"])
            branch_ratio.append(r["branch_m3"] / base["branch_m3"] if r["branch_m3"] and base["branch_m3"] else None)
            cyl_ratio.append(r["n_cylinders"] / base["n_cylinders"] if r["n_cylinders"] and base["n_cylinders"] else None)

        ax_top.plot(xs, trunk_ratio, marker="o", linestyle="-", color=tree_color[tree], linewidth=DEFAULT_LINEWIDTH)
        ax_top.plot(xs, branch_ratio, marker="s", linestyle="--", color=tree_color[tree], linewidth=DEFAULT_LINEWIDTH)
        ax_bot.plot(xs, cyl_ratio, marker="o", linestyle="-", color=tree_color[tree], linewidth=DEFAULT_LINEWIDTH)

    ax_top.axhline(1.0, color="grey", linewidth=1, linestyle=":", zorder=0)
    ax_bot.axhline(1.0, color="grey", linewidth=1, linestyle=":", zorder=0)
    ax_top.set_ylabel("Volume, relative to\nthat tree's pd2min=0.002 value", fontsize=AXIS_LABEL_FONTSIZE)
    ax_bot.set_ylabel("n_cylinders, relative to\nthat tree's pd2min=0.002 value", fontsize=AXIS_LABEL_FONTSIZE)
    ax_bot.set_xlabel("PatchDiam2Min [m]", fontsize=AXIS_LABEL_FONTSIZE)
    ax_bot.set_xticks(PD2MIN_LEVELS_MANUAL)
    ax_bot.set_xticklabels(["%.3f" % L for L in PD2MIN_LEVELS_MANUAL])
    ax_top.set_title("Cost of pd2min, normalised per tree (Simplified, sr=0.010, ri=1)",
                      fontsize=PANEL_TITLE_FONTSIZE)
    ax_top.grid(True, alpha=0.2)
    ax_bot.grid(True, alpha=0.2)

    # ONE shared legend for the whole figure (project convention - see the
    # original script's own small-multiples figure above): line STYLE
    # (trunk vs branch) and per-tree COLOR combined into one handle list.
    style_handles = [mlines.Line2D([0], [0], color="black", linestyle="-", marker="o", label="trunk_m3"),
                      mlines.Line2D([0], [0], color="black", linestyle="--", marker="s", label="branch_m3")]
    tree_handles = [mlines.Line2D([0], [0], color=tree_color[t], linestyle="-", marker="o",
                                   label=t.replace("B21_", "")) for t in dbh_tree_order]
    fig4.legend(handles=style_handles + tree_handles, fontsize=LEGEND_FONTSIZE - 2,
                ncol=1, loc="center left", bbox_to_anchor=(1.0, 0.5))
    fig4.tight_layout(rect=(0, 0, 0.86, 1))
    cost_plot_path = os.path.join(ensure_all_plots_dir(), COST_PLOT_NAME)
    fig4.savefig(cost_plot_path, dpi=PLOT_DPI, bbox_inches="tight")
    plt.close(fig4)
    print("Saved:", cost_plot_path)
    print()

    # =====================  STEP 5 (extension) - the answer  ================
    best_count = {L: 0 for L in PD2MIN_LEVELS_MANUAL}
    best_count["AUTO"] = 0
    per_level_err_all_trees = {L: [] for L in PD2MIN_LEVELS_MANUAL}
    per_level_err_all_trees["AUTO"] = []

    for tree in tree_order:
        rows_t = [r for r in per_model_rows if r["tree"] == tree]
        per_level = {}
        for r in rows_t:
            if r["err_above_base_pct"] is None:
                continue
            lvl = _pd2min_level_of_row(r)
            per_level[lvl] = r["err_above_base_pct"]
            per_level_err_all_trees[lvl].append(r["err_above_base_pct"])
        if not per_level:
            continue
        best = min(per_level, key=per_level.get)
        best_count[best] += 1

    level_medians = {L: (float(np.median(per_level_err_all_trees[L])) if per_level_err_all_trees[L] else None)
                      for L in PD2MIN_LEVELS_MANUAL}
    auto_median = float(np.median(per_level_err_all_trees["AUTO"])) if per_level_err_all_trees["AUTO"] else None

    manual_medians_known = {L: v for L, v in level_medians.items() if v is not None}
    best_median_level = min(manual_medians_known, key=manual_medians_known.get) if manual_medians_known else None
    worst_median_level = max(manual_medians_known, key=manual_medians_known.get) if manual_medians_known else None
    spread_pp = (manual_medians_known[worst_median_level] - manual_medians_known[best_median_level]) \
        if best_median_level is not None else None

    # Cylinder-count / trunk-volume ratio AT the best-MEDIAN level, relative
    # to pd2min=0.002 - one ratio per tree that has an sr010/ri1 export at
    # BOTH 0.002 and the best-median level (never substituted, same rule as
    # everywhere else in this extension - see the missing-match note above
    # for which trees/levels that excludes).
    cyl_ratios_at_best, trunk_ratios_at_best = [], []
    if best_median_level is not None:
        for tree in tree_order:
            rows_by_level = {_pd2min_level_of_row(r): r for r in per_model_rows
                              if r["tree"] == tree and r["mode"] == "manual"}
            base = rows_by_level.get(0.002)
            at_best = rows_by_level.get(best_median_level)
            if base and at_best and base["n_cylinders"] and at_best["n_cylinders"]:
                cyl_ratios_at_best.append(at_best["n_cylinders"] / base["n_cylinders"])
            if base and at_best and base["trunk_m3"] and at_best["trunk_m3"]:
                trunk_ratios_at_best.append(at_best["trunk_m3"] / base["trunk_m3"])

    print("=" * 90)
    print("=== TAPER ERROR PER pd2min ===")
    print("=" * 90)
    print("Models per tree :", models_per_tree_ext, "(expected 5)")
    print()
    print("err_above_base_pct, median across 12 trees:")
    for L in PD2MIN_LEVELS_MANUAL:
        v = level_medians[L]
        print("  pd2min %.3f : %s" % (L, ("%.4f %%" % v) if v is not None else "no data"))
    print("  AUTO         : %s" % (("%.4f %%" % auto_median) if auto_median is not None else "no data"))
    print()
    print("Best level per tree (count):")
    print("  0.002 : %d trees    0.005 : %d trees" % (best_count[0.002], best_count[0.005]))
    print("  0.010 : %d trees    0.020 : %d trees    AUTO : %d trees"
          % (best_count[0.010], best_count[0.020], best_count["AUTO"]))
    n_distinct_best = sum(1 for v in best_count.values() if v > 0)
    if n_distinct_best > 1:
        print("  -> TREES DISAGREE on the best level (%d of 5 candidate levels win for at least" % n_distinct_best)
        print("     one tree each) - the medians above do NOT reflect a consensus choice; see the")
        print("     per-tree breakdown in trunk_taper_per_pd2min.csv / the figure for the full picture.")
    print()
    if spread_pp is not None:
        print("Difference between the best (%.3f) and worst (%.3f) MEDIAN level : %.4f percentage points"
              % (best_median_level, worst_median_level, spread_pp))
    else:
        print("Difference between the best and worst median level : no data")
    print()
    if cyl_ratios_at_best:
        print("Cylinder count at the best level (pd2min=%.3f), relative to pd2min=0.002, median : %.4f (n=%d trees)"
              % (best_median_level, float(np.median(cyl_ratios_at_best)), len(cyl_ratios_at_best)))
    else:
        print("Cylinder count at the best level, relative to pd2min=0.002, median : no data")
    if trunk_ratios_at_best:
        print("Trunk volume  at the best level (pd2min=%.3f), relative to pd2min=0.002, median : %.4f (n=%d trees)"
              % (best_median_level, float(np.median(trunk_ratios_at_best)), len(trunk_ratios_at_best)))
    else:
        print("Trunk volume  at the best level, relative to pd2min=0.002, median : no data")
    print()
    print("Outputs : %s" % per_model_csv_path)
    print("          %s , %s" % (err_pd2min_plot_path, cost_plot_path))
    print()
    print("Sanity checks: 1 %s  2 %s" % ("PASS" if check_1 else "FAIL", "PASS" if check_2 else "FAIL"))
