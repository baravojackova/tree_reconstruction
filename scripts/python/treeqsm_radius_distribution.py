# -*- coding: utf-8 -*-
# =====================================================================
#  DIAGNOSTIC (read-only) - does TreeQSM have enough branch length below
#  5/10 mm radius for an AdTree-style radius threshold to make the two
#  families comparable, or would it change almost nothing?
# ---------------------------------------------------------------------
#  WHY THIS SCRIPT EXISTS: volume_results.csv already showed that
#  n_cylinders means something different in each family - AdTree's radius
#  threshold DISCARDS thin cylinders (branch length 2021 m -> 506 m at
#  5 mm -> 20 mm on B21_S07), while TreeQSM's sr/ri simplification MERGES
#  cylinders along a branch without changing branch length at all (4373 m,
#  identical across all 16 sr/ri cells). Applying AdTree's threshold to
#  TreeQSM at export time has been proposed as a fix - but TreeQSM's
#  MinCylRad already acts as a floor on cylinder radius, so the fix only
#  works if a meaningful share of TreeQSM's branch length actually sits
#  below 5-10 mm. This script measures that, directly from the .mat
#  cylinder data, and does nothing else - no threshold is applied to any
#  export, no file used elsewhere in the pipeline is modified.
#
#  HARD RULES THIS SCRIPT OBEYS (see the task this was written for):
#   1. Never deletes/moves/renames/reorganises anything under
#      scripts/matlab/archive/ - read-only, os.path.exists() and
#      scipy.io.loadmat() are the only operations performed on it.
#   2. Runs no reconstruction (no MATLAB, no runsken.m, no AdTree/AdQSM).
#   3. Never writes to volume_results.csv, any taper.txt, or any existing
#      .py file - this module reads volume_results.csv via
#      compare_volumes.load_results() (import, not reimplementation) and
#      writes only its own three declared outputs.
#   4. Outputs: this script, one CSV, one PNG - nothing else.
#   5. No threshold is ever applied to write geometry anywhere - every
#      "cut at N mm" number below is a measurement (a sum over the cylinders
#      that WOULD survive), not an action.
#
#  STEP 1 - ARCHIVE FOLDER RESOLUTION, worth reading before the code:
#  archive folder names follow scripts/matlab/archive/<timestamp>_<TREE>/,
#  auto-discovered below by regex per tree. Eleven of twelve trees resolve
#  to exactly one such folder. B21_S21's does NOT - it is named
#  "...BK21_S21" (typo: "BK21" instead of "B21"), so the regex finds zero
#  matches for "B21_S21". This exact inconsistency is not a new finding -
#  it is already hand-verified and documented in this repo, in
#  base_vs_dbh_compare.py's TREEQSM_ARCHIVE_DIR + its own comment ("folder
#  names are hand-typed and inconsistent (e.g. "BK21_S21" for B21_S21)").
#  Every filename INSIDE that folder is correctly tagged "B21_S21"
#  (confirmed: branchlen_B21_S21_..., dbh_B21_S21_..., etc.) - only the
#  folder's own name has the typo. Given that prior, independent
#  verification already exists in this codebase, FOLDER_NAME_EXCEPTIONS
#  below applies it explicitly (asked of and approved by the user for
#  this run) rather than silently either trusting or discarding the regex
#  result.
# =====================================================================

import os
import re
import csv

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

from paths import ensure_csv_dir, ensure_all_plots_dir
from compare_volumes import RESULTS_CSV, load_results
from base_vs_dbh_compare import parse_field_file, FIELD_FILE
from plot_style import TREE_MARKERS, DEFAULT_LINEWIDTH, LEGEND_FONTSIZE, \
    AXIS_LABEL_FONTSIZE, TITLE_FONTSIZE, PANEL_TITLE_FONTSIZE

# =====================  PARAMETERS  ====================================

# ---- Archive location and model selection -----------------------------
# scripts/matlab/archive, relative to scripts/python (this project's own
# launch-directory convention - see paths.py's own header comment).
ARCHIVE_ROOT = os.path.join("..", "matlab", "archive")

# One TreeQSM archive folder per tree, matching
# <YYYY-MM-DD_HHMM>_<TREE> exactly (Step 1). Auto-discovered by regex,
# never hand-typed, so a future re-run catches a folder rename/move.
FOLDER_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}_\d{4}_(?P<tree>.+)$")

# The one documented exception to strict regex matching - see module
# header. Applied ON TOP OF the regex result, and reported as such in
# Step 1's output, never silently.
FOLDER_NAME_EXCEPTIONS = {
    "B21_S21": "2026-09-09_1538_BK21_S21",
}

# tree -> (pd-tag at pd2min=0.002, pd-tag at pd2min=0.010), both at
# sr=0.010, ri=1, mo9, Simplified - given by the task, not derived.
MODEL_TABLE = {
    "B21_S01": ("man_pd07-002-10_mo9_sr010_ri1",   "man_pd07-01-10_mo9_sr010_ri1"),
    "B21_S04": ("man_pd05-002-10_mo9_sr010_ri1",   "man_pd05-01-10_mo9_sr010_ri1"),
    "B21_S07": ("man_pd08-002-14_mo9_sr010_ri1",   "man_pd08-01-14_mo9_sr010_ri1"),
    "B21_S08": ("man_pd05-002-09_mo9_sr010_ri1",   "man_pd05-01-09_mo9_sr010_ri1"),
    "B21_S09": ("man_pd05-002-08_mo9_sr010_ri1",   "man_pd05-01-08_mo9_sr010_ri1"),
    "B21_S10": ("man_pd05-002-09_mo9_sr010_ri1",   "man_pd05-01-09_mo9_sr010_ri1"),
    "B21_S11": ("man_pd06-002-11_mo9_sr010_ri1",   "man_pd06-01-11_mo9_sr010_ri1"),
    "B21_S12": ("man_pd04-002-05_mo9_sr010_ri1",   "man_pd04-01-05_mo9_sr010_ri1"),
    "B21_S13": ("man_pd06-002-10_mo9_sr010_ri1",   "man_pd06-01-10_mo9_sr010_ri1"),
    "B21_S14": ("man_pd05-002-09_mo9_sr010_ri1",   "man_pd05-01-09_mo9_sr010_ri1"),
    "B21_S16": ("man_pd045-002-05_mo9_sr010_ri1",  "man_pd045-01-05_mo9_sr010_ri1"),
    "B21_S21": ("man_pd066-002-099_mo9_sr010_ri1", "man_pd066-01-099_mo9_sr010_ri1"),
}
PD2MIN_OF_TAG_SLOT = (0.002, 0.010)   # matches MODEL_TABLE's tuple order

# ---- Radius thresholds (Step 3 / Step 4) -------------------------------
RADIUS_THRESHOLDS_MM = [2.5, 5, 7.5, 10, 15, 20]      # Step 3 sweep
ADTREE_COMPARE_THRESHOLDS_MM = [5, 10, 20]            # Step 4 side-by-side

# ---- Hard-floor detection ----------------------------------------------
# A radius value counts as "a visible floor" if this many percent or more
# of a model's branch cylinders sit within FLOOR_CLUSTER_TOL_MM of one
# single rounded radius value (rounded to FLOOR_ROUND_MM precision).
FLOOR_ROUND_MM = 0.01
FLOOR_SHARE_PCT_THRESHOLD = 1.0
# Separately: how close (in mm) a cylinder radius must be to
# rundata.inputs.MinCylRad to count as "sitting on the MinCylRad floor
# itself" (distinct from the "most common rounded value" check above -
# see Step 3's report for why both are checked).
MINCYLRAD_TOL_MM = 0.05

# ---- AdTree settled-configuration filter (Step 4) ----------------------
# Matches method_evaluation.py's POP_SETTLED definition for AdTree_calibrated
# (calmethod == "regression-perorder", "-wlen" weighting, the tree's one
# non-999 AdQSM variant) plus the actual exported segmentation config
# "seg50-500-k80-wlen" (see export_geom_ansys.py's NPZ_FILE constant) -
# not re-derived here, just pinned to the same settled values.
ADTREE_CALMETHOD = "regression-perorder"
ADTREE_SEG_TAG = "seg50-500-k80"

# ---- Sanity-check anchors (Step 2), pd2min = 0.002 ----------------------
SANITY_ANCHORS = {
    "B21_S07": {"n_cylinders": 19220, "branch_m3": 10.725, "branch_len_m": 4373.008},
    "B21_S16": {"n_cylinders": 2611,  "branch_m3": 0.8189},
    "B21_S01": {"n_cylinders": 21864, "branch_m3": 10.4292},
}
SANITY_TOL_REL = 0.02   # 2% relative tolerance (see Step 2 report for why
                         # branch_m3 needs a tolerance at all, not exact
                         # match, even when reading the right file)

# ---- Outputs (Hard rule 4 - exactly these three, nothing else) ---------
CSV_OUT_NAME = "treeqsm_radius_distribution.csv"
PLOT_OUT_NAME = "treeqsm_radius_distribution.png"

FIGSIZE = (13, 6)
PLOT_DPI = 150


# =====================  STEP 1 - ARCHIVE RESOLUTION  ====================

def resolve_archive_folders():
    """Auto-discover one archive folder per MODEL_TABLE tree by regex
    (Step 1), print every tree's match count (Sanity check 1), apply the
    one documented exception on top, and return {tree: folder_name}.
    Trees with zero matches AND no exception stay out of `resolved`
    (reported as MISSING, not silently skipped)."""
    all_dirs = [d for d in os.listdir(ARCHIVE_ROOT)
                if os.path.isdir(os.path.join(ARCHIVE_ROOT, d))]

    print("--- STEP 1: archive folder resolution (SANITY CHECK 1) ---")
    resolved = {}
    for tree in MODEL_TABLE:
        pattern = re.compile(r"^\d{4}-\d{2}-\d{2}_\d{4}_%s$" % re.escape(tree))
        matches = sorted(d for d in all_dirs if pattern.match(d))
        if len(matches) == 1:
            resolved[tree] = matches[0]
            print("  %-10s -> 1 match : %s" % (tree, matches[0]))
        elif len(matches) > 1:
            print("  %-10s -> %d matches (STOPPING - ambiguous):" % (tree, len(matches)))
            for m in matches:
                print("       %s" % m)
            raise SystemExit("Multiple archive folders matched %s - resolve by hand "
                              "before re-running (see Step 1's own instructions)." % tree)
        else:
            if tree in FOLDER_NAME_EXCEPTIONS:
                exc = FOLDER_NAME_EXCEPTIONS[tree]
                exists = os.path.isdir(os.path.join(ARCHIVE_ROOT, exc))
                print("  %-10s -> 0 regex matches; MISSING under strict pattern "
                      "matching. Using documented exception (user-approved, see "
                      "module header) -> %s (exists on disk: %s)"
                      % (tree, exc, exists))
                if exists:
                    resolved[tree] = exc
            else:
                print("  %-10s -> MISSING (no folder matches the pattern)" % tree)

    print("Archive folders matched : %d / %d trees" % (len(resolved), len(MODEL_TABLE)))
    return resolved


# =====================  MODEL LOADING  ==================================

def load_model(tree, folder, tag):
    """Load one simplified_<tree>_<tag>.mat (QSM_final) plus its sibling
    branchlen_<tree>_<tag>.txt. Returns a dict of everything Steps 2-4
    need, or None if the .mat file does not exist (reported, not raised -
    Step 1's own MISSING/STOP protocol already covers absent archives;
    this only guards an unexpectedly absent single file within a present
    one)."""
    mat_path = os.path.join(ARCHIVE_ROOT, folder, "simplified_%s_%s.mat" % (tree, tag))
    if not os.path.exists(mat_path):
        print("    MISSING FILE: %s" % mat_path)
        return None

    d = sio.loadmat(mat_path, struct_as_record=False, squeeze_me=True)
    qsm = d["QSM_final"]
    cyl = qsm.cylinder

    radius = np.asarray(cyl.radius, dtype=float)
    length = np.asarray(cyl.length, dtype=float)
    branch_order = np.asarray(cyl.BranchOrder, dtype=float)
    branch_mask = branch_order > 0

    branchlen_path = os.path.join(ARCHIVE_ROOT, folder, "branchlen_%s_%s.txt" % (tree, tag))
    branch_len_txt = None
    if os.path.exists(branchlen_path):
        with open(branchlen_path, "r") as f:
            branch_len_txt = float(f.read().strip())

    return {
        "tree": tree, "folder": folder, "tag": tag,
        "radius": radius, "length": length, "branch_mask": branch_mask,
        "n_cylinders": len(radius),
        "min_cyl_rad_input_m": float(qsm.rundata.inputs.MinCylRad),
        # treedata.BranchVolume/1000: what volume_results.csv's branch_m3
        # actually is (see Step 2's report) - NOT the naive
        # sum(pi*r^2*l) over branch cylinders, which runs ~1-1.3% lower
        # here (that sum is still computed below, for Step 3's per-
        # cylinder measurement, where it is the right quantity - see
        # module header / Step 2 report for why these differ).
        "branch_m3_treedata": float(qsm.treedata.BranchVolume) / 1000.0,
        "branch_len_txt_m": branch_len_txt,
    }


# =====================  STEP 3 - THE MEASUREMENT  =======================

def cylinder_measurements(model):
    """Per-cylinder branch (BranchOrder>0) radius/length/volume arrays,
    plus the cumulative-below-threshold fractions (count/length/volume)
    at every RADIUS_THRESHOLDS_MM value, plus the percentile/floor
    diagnostics Step 3 asks for. Pure measurement - nothing written
    anywhere, no threshold applied to any export."""
    r = model["radius"][model["branch_mask"]]
    l = model["length"][model["branch_mask"]]
    vol = np.pi * r ** 2 * l

    total_len = l.sum()
    total_vol = vol.sum()
    n = len(r)

    fractions = {}
    for t_mm in RADIUS_THRESHOLDS_MM:
        below = r < (t_mm / 1000.0)
        fractions[t_mm] = {
            "count": below.sum() / n if n else None,
            "length": l[below].sum() / total_len if total_len else None,
            "volume": vol[below].sum() / total_vol if total_vol else None,
        }

    # Hard-floor detection: (a) share of cylinders within MINCYLRAD_TOL_MM
    # of rundata.inputs.MinCylRad itself, (b) share sitting on the single
    # most common rounded radius value (a floor need not equal
    # MinCylRad - simplification/replacement can introduce its own,
    # different, round-number floor - see Step 3's report for what was
    # actually found).
    mincylrad_mm = model["min_cyl_rad_input_m"] * 1000.0
    near_mincylrad = np.isclose(r * 1000.0, mincylrad_mm, atol=MINCYLRAD_TOL_MM)
    r_rounded_mm = np.round(r * 1000.0, 2)
    vals, counts = np.unique(r_rounded_mm, return_counts=True)
    top_i = np.argmax(counts)
    mode_value_mm, mode_count = float(vals[top_i]), int(counts[top_i])

    return {
        "branch_cyl_n": n,
        "branch_len_cylsum_m": total_len,
        "branch_vol_cylsum_m3": total_vol,
        "fractions": fractions,
        "radius_min_mm": float(r.min()) * 1000.0 if n else None,
        "radius_p1_mm": float(np.percentile(r, 1)) * 1000.0 if n else None,
        "radius_p5_mm": float(np.percentile(r, 5)) * 1000.0 if n else None,
        "radius_median_mm": float(np.percentile(r, 50)) * 1000.0 if n else None,
        "near_mincylrad_share_pct": 100.0 * near_mincylrad.sum() / n if n else None,
        "mode_radius_mm": mode_value_mm,
        "mode_radius_share_pct": 100.0 * mode_count / n if n else None,
    }


# =====================  STEP 4 - ADTREE SETTLED FIGURES  ================

def adtree_settled_branch_len(results_rows, tree):
    """{radius_threshold_mm: branch_len_m} for tree's AdTree-calibrated
    settled row (Step 4) - calmethod==regression-perorder, "-wlen"
    weighting, ADTREE_SEG_TAG segmentation, and the tree's one non-999
    AdQSM variant (matches method_evaluation.py's is_settled() for
    AdTree_calibrated, reused by definition rather than re-imported to
    avoid pulling in that module's whole POP_SETTLED machinery for one
    lookup)."""
    variants = sorted({
        r["adqsm_variant"] for r in results_rows
        if r["tree"] == tree and r["method"].startswith("AdTree calibrated")
        and r["adqsm_variant"] not in ("", "999")
    })
    if len(variants) != 1:
        print("    WARNING: %s has %d non-999 AdQSM variants for AdTree "
              "calibrated (expected 1): %s" % (tree, len(variants), variants))
    variant = variants[0] if variants else None

    out = {}
    for r in results_rows:
        if (r["tree"] != tree or not r["method"].startswith("AdTree calibrated")
                or r["calmethod"] != ADTREE_CALMETHOD
                or "-wlen" not in r["method"]
                or ADTREE_SEG_TAG not in r["method"]
                or r["adqsm_variant"] != variant):
            continue
        t_mm = r["radius_threshold_mm"]
        if t_mm in ADTREE_COMPARE_THRESHOLDS_MM:
            out[t_mm] = r["branch_len"]
    return out


# =====================  MAIN  ============================================

def main():
    resolved_folders = resolve_archive_folders()

    print()
    print("--- STEP 2/3: loading + verifying + measuring each model ---")
    results_rows = load_results(RESULTS_CSV)
    field_dbh_cm = {t: dbh for t, (_base, dbh) in parse_field_file(FIELD_FILE).items()}

    csv_rows = []
    sanity_lines = []
    models_loaded = 0
    adtree_cache = {}

    for tree in resolved_folders:
        for slot, tag in enumerate(MODEL_TABLE[tree]):
            pd2min = PD2MIN_OF_TAG_SLOT[slot]
            folder = resolved_folders[tree]
            full_path = os.path.join(ARCHIVE_ROOT, folder, "simplified_%s_%s.mat" % (tree, tag))
            print("  %s pd2min=%.3f -> %s" % (tree, pd2min, full_path))

            model = load_model(tree, folder, tag)
            if model is None:
                continue
            models_loaded += 1
            meas = cylinder_measurements(model)

            # ---- Step 2: cross-check against volume_results.csv --------
            csv_method = "TreeQSM mine (%s, Simplified (no islands))" % tag
            csv_row = next((r for r in results_rows
                             if r["tree"] == tree and r["method"] == csv_method), None)
            csv_n = csv_row["n_cylinders"] if csv_row else None
            csv_branch_m3 = csv_row["branch"] if csv_row else None
            csv_branch_len = csv_row["branch_len"] if csv_row else None
            csv_mincylrad = csv_row["mincylrad"] if csv_row else None

            n_ok = (csv_n is not None and model["n_cylinders"] == int(csv_n))
            bm3_ok = (csv_branch_m3 is not None and abs(model["branch_m3_treedata"] - csv_branch_m3)
                      <= SANITY_TOL_REL * csv_branch_m3)
            blen_ok = (csv_branch_len is not None and model["branch_len_txt_m"] is not None
                       and abs(model["branch_len_txt_m"] - csv_branch_len) <= SANITY_TOL_REL * csv_branch_len)

            if tree in SANITY_ANCHORS:
                a = SANITY_ANCHORS[tree]
                if pd2min == 0.002:
                    line = ("  ANCHOR %-8s pd2min 0.002 : n_cylinders actual=%d expected=%d | "
                            "branch_m3 actual=%.4f expected=%.4f" %
                            (tree, model["n_cylinders"], a["n_cylinders"],
                             model["branch_m3_treedata"], a["branch_m3"]))
                    if "branch_len_m" in a:
                        line += (" | branch_len_m(txt) actual=%.3f expected=%.3f" %
                                 (model["branch_len_txt_m"], a["branch_len_m"]))
                    sanity_lines.append(line)

            # ---- assemble one CSV output row ---------------------------
            row = {
                "tree": tree, "pd2min_m": pd2min, "archive_folder": folder, "model_tag": tag,
                "field_dbh_cm": field_dbh_cm.get(tree),
                "n_cylinders": model["n_cylinders"], "n_cylinders_csv": csv_n, "n_cylinders_match": n_ok,
                "branch_m3_treedata": round(model["branch_m3_treedata"], 6),
                "branch_m3_csv": csv_branch_m3, "branch_m3_match_within_2pct": bm3_ok,
                "branch_len_txt_m": model["branch_len_txt_m"], "branch_len_csv_m": csv_branch_len,
                "branch_len_txt_match_csv": blen_ok,
                "branch_len_cylsum_m": round(meas["branch_len_cylsum_m"], 3),
                "branch_vol_cylsum_m3": round(meas["branch_vol_cylsum_m3"], 6),
                "min_cyl_rad_input_mm": round(model["min_cyl_rad_input_m"] * 1000.0, 4),
                "min_cyl_rad_csv_mm": (round(csv_mincylrad * 1000.0, 4) if csv_mincylrad is not None else None),
                "branch_radius_min_mm": round(meas["radius_min_mm"], 4),
                "branch_radius_p1_mm": round(meas["radius_p1_mm"], 4),
                "branch_radius_p5_mm": round(meas["radius_p5_mm"], 4),
                "branch_radius_median_mm": round(meas["radius_median_mm"], 4),
                "near_mincylrad_share_pct": round(meas["near_mincylrad_share_pct"], 4),
                "mode_radius_mm": meas["mode_radius_mm"],
                "mode_radius_share_pct": round(meas["mode_radius_share_pct"], 4),
                "hard_floor_detected": meas["mode_radius_share_pct"] >= FLOOR_SHARE_PCT_THRESHOLD,
            }
            for t_mm in RADIUS_THRESHOLDS_MM:
                f = meas["fractions"][t_mm]
                key = ("%.1f" % t_mm).replace(".", "_")
                row["frac_count_below_%smm" % key] = round(f["count"], 6) if f["count"] is not None else None
                row["frac_length_below_%smm" % key] = round(f["length"], 6) if f["length"] is not None else None
                row["frac_volume_below_%smm" % key] = round(f["volume"], 6) if f["volume"] is not None else None

            # ---- Step 4: AdTree settled figures + the cut comparison ---
            if tree not in adtree_cache:
                adtree_cache[tree] = adtree_settled_branch_len(results_rows, tree)
            adtree_len = adtree_cache[tree]
            for t_mm in ADTREE_COMPARE_THRESHOLDS_MM:
                row["adtree_r%dmm_branch_len_m" % t_mm] = adtree_len.get(t_mm)
                frac_below = meas["fractions"][t_mm]["length"]
                cut_len = (meas["branch_len_cylsum_m"] * (1.0 - frac_below)
                           if frac_below is not None else None)
                row["treeqsm_cut%dmm_branch_len_m" % t_mm] = round(cut_len, 3) if cut_len is not None else None
                if cut_len is not None and adtree_len.get(t_mm):
                    row["ratio_cut%dmm_vs_adtree" % t_mm] = round(cut_len / adtree_len[t_mm], 3)
                else:
                    row["ratio_cut%dmm_vs_adtree" % t_mm] = None

            csv_rows.append(row)

    print()
    print("--- SANITY CHECK 2 anchors ---")
    for line in sanity_lines:
        print(line)

    # =====================  WRITE CSV  ===================================
    csv_path = os.path.join(ensure_csv_dir(), CSV_OUT_NAME)  # ensure_csv_dir() returns the csv folder
    fieldnames = list(csv_rows[0].keys()) if csv_rows else []
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        # Trees ordered by field DBH ascending (project convention - see
        # method_evaluation.py's compute_tree_order_by_dbh()), ties by name.
        def _dbh_key(row):
            dbh = row["field_dbh_cm"]
            return (dbh if dbh is not None else 1e9, row["tree"], row["pd2min_m"])
        for row in sorted(csv_rows, key=_dbh_key):
            writer.writerow(row)
    print("\nWrote %s (%d rows)" % (csv_path, len(csv_rows)))

    # =====================  FIGURE  =======================================
    draw_figure(csv_rows, field_dbh_cm)

    # =====================  STEP 5 - REPORT  ==============================
    print_final_report(resolved_folders, models_loaded, csv_rows, sanity_lines)


def draw_figure(csv_rows, field_dbh_cm):
    """Cumulative branch LENGTH vs. cylinder radius, one line per tree, one
    panel per pd2min level (Step 4's figure spec) - built directly from
    the per-cylinder radius/length arrays (re-loaded here per tree/tag,
    since csv_rows only holds the already-bucketed fractions at 6 fixed
    thresholds - the figure needs the full curve, not just those 6
    points)."""
    tree_order = sorted(MODEL_TABLE.keys(),
                         key=lambda t: (field_dbh_cm.get(t, 1e9), t))
    resolved = {row["tree"]: row["archive_folder"] for row in csv_rows}
    colors = plt.cm.tab20(np.linspace(0, 1, len(tree_order)))
    color_of = {t: colors[i] for i, t in enumerate(tree_order)}

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, sharey=True)
    pd2min_levels = [0.002, 0.010]

    for ax, pd2min in zip(axes, pd2min_levels):
        slot = PD2MIN_OF_TAG_SLOT.index(pd2min)
        for ti, tree in enumerate(tree_order):
            if tree not in resolved:
                continue
            tag = MODEL_TABLE[tree][slot]
            model = load_model(tree, resolved[tree], tag)
            if model is None:
                continue
            r_mm = model["radius"][model["branch_mask"]] * 1000.0
            l = model["length"][model["branch_mask"]]
            order = np.argsort(r_mm)
            r_sorted = r_mm[order]
            cum_len_frac = np.cumsum(l[order]) / l.sum()
            ax.plot(r_sorted, cum_len_frac, color=color_of[tree],
                    marker=TREE_MARKERS[ti % len(TREE_MARKERS)], markevery=0.1,
                    markersize=4, linewidth=DEFAULT_LINEWIDTH,
                    label=tree.replace("B21_", ""))

        ax.set_xscale("log")
        ax.set_xlim(1, 100)
        ax.set_ylim(0, 1)
        for t_mm in ADTREE_COMPARE_THRESHOLDS_MM:
            ax.axvline(t_mm, color="grey", linestyle="--", linewidth=1, zorder=0)
        ax.set_xlabel("cylinder radius [mm] (log scale)", fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_title("pd2min = %.3f m" % pd2min, fontsize=PANEL_TITLE_FONTSIZE)
        ax.grid(True, which="both", alpha=0.25)

    axes[0].set_ylabel("cumulative fraction of branch length", fontsize=AXIS_LABEL_FONTSIZE)
    axes[1].legend(fontsize=LEGEND_FONTSIZE - 2, loc="lower right", ncol=2,
                    title="tree (DBH asc.)")
    fig.suptitle("TreeQSM: cumulative branch length below cylinder radius\n"
                 "(dashed lines: 5 / 10 / 20 mm, AdTree's own threshold values)",
                 fontsize=TITLE_FONTSIZE)
    fig.tight_layout(rect=(0, 0, 1, 0.92))

    out_dir = ensure_all_plots_dir()
    out_path = os.path.join(out_dir, PLOT_OUT_NAME)
    fig.savefig(out_path, dpi=PLOT_DPI)
    plt.close(fig)
    print("Wrote %s" % out_path)


def print_final_report(resolved_folders, models_loaded, csv_rows, sanity_lines):
    n_trees_matched = len(resolved_folders)
    multiples = []   # SANITY CHECK 1 already raises SystemExit on multiples,
                      # so this is always empty by the time we get here.

    mincylrad_values = sorted({row["min_cyl_rad_input_mm"] for row in csv_rows})

    # Hard-floor summary across all loaded models.
    floor_rows = [r for r in csv_rows if r["hard_floor_detected"]]
    if floor_rows:
        floor_val = np.median([r["mode_radius_mm"] for r in floor_rows])
        floor_share = np.median([r["mode_radius_share_pct"] for r in floor_rows])
        floor_txt = "YES  at ~%.2f mm, holding ~%.2f %% of branch cylinders (median over %d/%d models flagged)" \
                    % (floor_val, floor_share, len(floor_rows), len(csv_rows))
    else:
        floor_txt = "NO (no model has >= %.1f%% of branch cylinders on one rounded radius value)" \
                    % FLOOR_SHARE_PCT_THRESHOLD

    def median_frac(key):
        vals = [r[key] for r in csv_rows if r.get(key) is not None]
        return 100.0 * np.median(vals) if vals else None

    ratios = [r["ratio_cut10mm_vs_adtree"] for r in csv_rows if r.get("ratio_cut10mm_vs_adtree") is not None]

    print()
    print("=== TREEQSM RADIUS DISTRIBUTION ===")
    print("Archive folders matched : %d / 12 trees   (multiples: %s)"
          % (n_trees_matched, multiples if multiples else "none"))
    print("Models loaded           : %d / %d" % (models_loaded, len(csv_rows)))
    print("MinCylRad values found  : %s mm" % mincylrad_values)
    print()
    print("Sanity check 2 (S07 / S16 / S01):")
    for line in sanity_lines:
        print(" ", line.strip())
    print()
    print("Hard floor detected : %s" % floor_txt)
    print()
    print("Fraction of branch LENGTH below threshold, median across %d models:" % len(csv_rows))
    print("   2.5 mm %.1f %%    5 mm %.1f %%    10 mm %.1f %%    20 mm %.1f %%" % (
        median_frac("frac_length_below_2_5mm") or 0, median_frac("frac_length_below_5_0mm") or 0,
        median_frac("frac_length_below_10_0mm") or 0, median_frac("frac_length_below_20_0mm") or 0))
    print("Fraction of branch VOLUME below threshold, median across %d models:" % len(csv_rows))
    print("   2.5 mm %.1f %%    5 mm %.1f %%    10 mm %.1f %%    20 mm %.1f %%" % (
        median_frac("frac_volume_below_2_5mm") or 0, median_frac("frac_volume_below_5_0mm") or 0,
        median_frac("frac_volume_below_10_0mm") or 0, median_frac("frac_volume_below_20_0mm") or 0))
    print()
    if ratios:
        print("TreeQSM branch length after a 10 mm cut, as a multiple of AdTree r10mm:")
        print("   min %.2f   median %.2f   max %.2f" % (min(ratios), np.median(ratios), max(ratios)))
    else:
        print("TreeQSM branch length after a 10 mm cut, as a multiple of AdTree r10mm: no data")

    print()
    print("NOTE on branch_len_m: volume_results.csv's branch_len_m column (and this "
          "script's branch_len_txt_m/branch_len_csv_m columns) come from "
          "branchlen_<tree>_<run>.txt, which runsken.m writes from QSM_opt.treedata."
          "BranchLength - the OPTIMAL stage, BEFORE simplification - not from the "
          "Simplified model these .mat files actually contain (see runsken.m around "
          "its trunklen_*/branchlen_*.txt export). This script's own branch-length "
          "measurements (branch_len_cylsum_m, the Step 3/4 threshold fractions, and "
          "the figure) are all computed from the Simplified-stage cylinders directly, "
          "which is why branch_len_cylsum_m runs below branch_len_txt_m/branch_len_csv_m "
          "here - both numbers are reported in the CSV so this is checkable, not hidden.")
    print()
    print("NOTE on B21_S21: included via the documented BK21_S21 folder-name exception "
          "(user-approved this run) - see module header and base_vs_dbh_compare.py's "
          "TREEQSM_ARCHIVE_DIR for the prior, independent verification this reuses.")


if __name__ == "__main__":
    main()
