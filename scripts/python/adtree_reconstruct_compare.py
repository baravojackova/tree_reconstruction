# -*- coding: utf-8 -*-
# =====================================================================
#  Calibrate an AdTree skeleton (.ply) against AdQSM, compare it to the
#  other reconstruction methods, and save the calibrated geometry.
# ---------------------------------------------------------------------
#  This is step 1 of a 2-step pipeline (see export_geom_ansys.py for step 2):
#
#   1) Reads a binary .ply (vertices with x,y,z,radius and edges).
#   2) Merges coincident vertices -> turns the disconnected skeleton into
#      a single connected tree with real branch points.
#   3) Roots the tree at the base (lowest z) and finds a unique parent
#      for every node (this removes any loops).
#   4) Smooths the centerline (Laplacian smoothing along branches) to
#      remove zig-zag noise from the AdTree skeleton.
#   5) Prunes thin twigs below a chosen radius threshold.
#   6) Resamples the dense skeleton into segments whose length adapts to
#      the local branch radius (keeps the curvature but drastically
#      reduces the element count).
#   7) Shifts the x,y coordinates to the origin (z is left untouched,
#      since it is the height above ground).
#   8) Optionally CALIBRATES the cylinder radii against AdQSM data: a taper
#      curve for the trunk, and a per-order-group log-log power-law
#      regression fit for the rest ([calmethod=regression-perorder], the
#      adopted primary method - see CHANGELOG_adtree.md), replacing the
#      AdTree radii before writing. A fixed-threshold median-ratio
#      calibration ("[calref=min5mm]") is also computed and written as a
#      secondary/backup reference point, alongside (not instead of) the
#      primary method.
#   9) Prints/upserts volume, DBH, height, taper, trunk/branch length
#      comparisons into the shared volume_results.csv (see compare_volumes.py).
#   10) SAVES the final calibrated geometry (xyz, cyl, root, RECENTER_XY, ...)
#       to a "calib_<tree>_r<mm>mm<variant>.npz" file per threshold/variant -
#       it does NOT write geom_*.txt itself any more (see the CHANGE note
#       right before the np.savez(...) call in the RUN section below for why,
#       and export_geom_ansys.py for the script that actually writes it).
#   11) Optionally shows/saves a 3D preview of the reduced model so you
#       can inspect it before importing into ANSYS.
#
#  All the actual geometry/calibration MATH lives in tree_geom_utils.py -
#  this script only holds the PARAMETERS and the RUN section that calls
#  those functions in the right order and prints/saves the results.
#
#  Dependencies: numpy, scipy, matplotlib   (install: pip install numpy scipy matplotlib)
# =====================================================================

import os
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from tree_geom_utils import (
    read_ply, merge_vertices, smooth_centerline,
    parse_adqsm_taper_file, make_trunk_radius_func, parse_adqsm_branch_file,
    parse_adqsm_params_file, print_adqsm_branch_file_sample, report_adqsm_thin_branch,
    compute_order_calibration_factors, apply_order_calibration_factors,
    parse_adqsm_branch_file_raw, build_quantile_matched_pairs, fit_radius_regression,
    group_orders_for_fitting, apply_radius_regression_per_order, plot_radius_regression_per_order,
    convert, cylinder_metrics, volume_stats,
    report_thin_branch_volume, upsert_result, stem_diameter_at_height,
    _fmt_dbh, _fmt_taper, print_volume_stats, raw_skeleton_stats,
    report_volume, plot_model,
)
# PLOTS_DIR: FIGURES_DIR below is built from this instead of a local bare
# "plots" literal, so paths.py stays the single place that name is
# defined project-wide (see paths.py's own header for why).
# ensure_tree_geom_dir(): the plot_model() call below (tree-SHAPE render
# PNG, derived from a geom_*.txt-style stub) is routed to plots/<tree>/
# geom/ - a DIFFERENT artifact from the adtree_adqsm_radius_regression_
# perorder_*.png diagnostic saved via FIGURES_DIR above, which stays
# directly under plots/<tree>/, not geom/ (see paths.py's own
# ensure_tree_geom_dir() docstring for this same "two different things
# both named 'geom'" distinction).
from paths import PLOTS_DIR, ensure_tree_geom_dir

# =====================  PARAMETERS  ===================================
# PRINT_TIMING: print wall-clock elapsed seconds (time.perf_counter()) for
# the initial read_ply+merge_vertices+smoothing stage, each convert() call
# (labelled by call site), each AdQSM variant iteration, each
# RADIUS_THRESHOLDS iteration, and the whole run - added to cost a
# SEG_LEN sweep before building it. Purely additive instrumentation - set
# to False to silence it without touching anything else.
PRINT_TIMING = True

# Tree ID for THIS run. This is the ONLY thing you need to change to switch
# trees - it names this tree's row in the shared results table (RESULTS_CSV,
# see upsert_result() calls below) AND builds AdQSM_DIR/AdTree_DIR/INPUT_PLY
# right below it automatically, so those don't need editing separately.
TREE_NAME = "B21_S08"

# Base folder holding every tree's data, one subfolder per tree named after
# TREE_NAME (e.g. ".../data/IND07_083/..."). Change this only if you move the
# whole "data" folder somewhere else - it does NOT depend on which tree you're
# processing.
DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"

# Directory holding this tree's source input files (the AdTree skeleton .ply,
# and the AdQSM taper/branch/params exports) - built from DATA_ROOT + TREE_NAME
# above. INPUT_PLY and the ADQSM_*_FILE paths below are all resolved relative
# to these two directories.
#
# --- AdQSM variant(s) ---------------------------------------------------
# AdQSM can be reconstructed several times with different settings, each
# saved in its own subfolder (e.g. ".../05", ".../08" - the folder name is
# just whatever you called that reconstruction run). You can either:
#   (1) point at ONE such folder with AdQSM_DIR (simple, old behaviour), or
#   (2) list SEVERAL subfolder names in ADQSM_VARIANTS to process all of
#       them in a single run of this script (similar to how RADIUS_THRESHOLDS
#       lets you try several radius thresholds in one run).
#
# Case (1) - single variant (default, still works exactly as before):
AdQSM_DIR = os.path.join(DATA_ROOT, TREE_NAME, "05")

# Case (2) - several variants. Leave ADQSM_VARIANTS empty/None to use only
# AdQSM_DIR above (case 1). To use several variants instead, set BOTH:
#   ADQSM_BASE_DIR = os.path.join(DATA_ROOT, TREE_NAME, "AdQSM")
#   ADQSM_VARIANTS = ["05", "08"]
# Each name in ADQSM_VARIANTS must be a subfolder of ADQSM_BASE_DIR that
# contains its own taper.txt, BranchStructure.txt and TreesParams.txt.
ADQSM_BASE_DIR = os.path.join(DATA_ROOT, TREE_NAME)
ADQSM_VARIANTS = ["060","999"]

AdTree_DIR = os.path.join(DATA_ROOT, TREE_NAME)

# input skeleton from AdTree - filename follows the "<TREE_NAME> - Cloud_skeleton.ply"
# convention used for every tree, so it's derived from TREE_NAME too.
INPUT_PLY = os.path.join(AdTree_DIR, "%s_noplate_clean_skeleton.ply" % TREE_NAME)

# Radius threshold (in METERS). You can give several values -> several variants.
# Remove all branches whose radius is below this threshold. The trunk (branch order 0) is never removed, even if its radius is below the threshold.
# Example of a single variant:   RADIUS_THRESHOLDS = [0.010]
# Example of several variants:   RADIUS_THRESHOLDS = [0.010, 0.020, 0.030]
RADIUS_THRESHOLDS = [0.005,0.01,0.015,0.02]        # 0.030 m = 30 mm radius (60 mm diameter)

# Fixed reference threshold(s) (in METERS) used to build the
# "[calref=minXmm]" calibration factors (see the "FIXED calibration
# reference set" block in the RUN section below) - one factors dict per
# entry, each producing its own labelled family of rows. Deliberately a
# SEPARATE list from RADIUS_THRESHOLDS (which controls final-model
# pruning, not calibration), NOT derived from it - if it were derived,
# changing RADIUS_THRESHOLDS would silently change what "[calref=minXmm]"
# is actually calibrated against, breaking comparability with rows
# already stored in volume_results.csv. Each entry's label is
# "[calref=min%dmm]" % round(ref_thr * 1000) - e.g. 0.005 ->
# "[calref=min5mm]".
#
# Decision (see CHANGELOG_adtree.md): per-order regression
# ([calmethod=regression-perorder]) is the PRIMARY calibration method
# going forward; calref=min5mm is kept only as a SECONDARY/backup
# reference point, so this list was collapsed back to its single min5mm
# entry (it previously also tested min2mm/min3mm/min4mm during the
# investigation - see the changelog). The list-based mechanism below is
# left as-is (not reverted to old single-value code) since it already
# handles one entry with no extra complexity.
CALIBRATION_REF_THRESHOLDS_MM = [0.005]

# Whether to compute/write the calref=min5mm secondary calibration variant
# at all. OFF by default since per-order regression
# ([calmethod=regression-perorder]) is the PRIMARY method - turn True only
# when you want to re-verify calref=min5mm against fresh data. When False,
# no calref=min5mm rows are computed or upserted this run, but any EXISTING
# calref=min5mm rows already in volume_results.csv from past runs are left
# untouched (last-verified snapshot) until this is re-enabled and re-run.
COMPUTE_CALREF_MIN5MM = True  # secondary/backup calibration check
# (fixed-reference median-ratio method) - OFF by default since
# per-order regression is the primary method; turn True only when
# you want to re-verify calref=min5mm against fresh data.

# Minimum quantile-matched pairs (see build_quantile_matched_pairs()) a
# branch order must have before it gets its OWN regression fit for the
# "[calmethod=regression-perorder]" calibration variant (see
# group_orders_for_fitting() in tree_geom_utils.py). Sparser orders are
# merged with the next order(s), walked in ascending order, until this
# minimum is reached. Tune this after seeing the printed raw per-order
# pair counts if 15 turns out too strict/loose for a given tree.
#TODO: revise fittign function
MIN_PAIRS_PER_ORDER = 15

# Adaptive segment length used for resampling (in METERS). The target length at
# a given point is SEG_LEN_K * local_radius, clamped to [SEG_LEN_MIN, SEG_LEN_MAX]:
# thick branches (trunk) get long segments, thin twigs get short/fine ones.
# If SEG_LEN_MIN == SEG_LEN_MAX, this reduces to the old constant-length
# resampling (that fixed value, regardless of radius).
# Set SEG_LEN_MIN to 0 or None to disable resampling entirely (keep every point).
# TODO: Check if seglen mim influence the regresion - see changes doc 
# CLAMPED DIAMETER [mm] = 2 × SEG_LEN_MIN [mm] / SEG_LEN_K
# LENGTH ≈ SEG_LEN_MIN
#
# SEG_LEN_MIN_LIST / SEG_LEN_K_LIST: swept as a grid - the RUN section below
# loops over every (SEG_LEN_MIN, SEG_LEN_K) combination in turn, instead of
# requiring a manual edit + rerun per combination. SEG_LEN_MIN/SEG_LEN_K (no
# _LIST suffix) are no longer set here as scalars - the sweep loop assigns
# the CURRENT combination's values into those exact names each iteration, so
# the three convert() call sites elsewhere in this file need no changes to
# their argument lists.
#
# Lower bound (0.035 m): sits just above the AdTree skeleton's own vertex
# spacing (measured median post-smoothing edge length = 0.0276 m, see
# scratch_skeleton_diag.py) - below that, resampling discards nothing, so a
# smaller SEG_LEN_MIN would just be a duplicate run of the same geometry,
# not a distinct one.
# Upper bound (0.065 m): keeps the resampling floor below the trunk radius of
# even the thinnest tree in the set (25 cm measured diameter, ~0.139 m raw
# AdTree radius) at the smallest K (0.5) in the grid, so the trunk always
# stays in the linear (non-floor) region of local_seg_len() for every
# combination in this grid.
SEG_LEN_MIN_LIST = [0.035, 0.05, 0.065]   # metres  -- STEP 1 clean invariance run: case B, same session as case A
SEG_LEN_K_LIST   = [0.5, 0.8, 1.0]    # dimensionless
SEG_LEN_MAX      = 0.5                     # longest allowed segment (m), for the trunk - unchanged, still scalar

# --- Short suffix identifying THIS resampling configuration -------------
# WHY THIS EXISTS: SEG_LEN_MIN/MAX/K (just above) control how densely the
# skeleton gets resampled, so changing any of them produces a DIFFERENT
# final cylinder model for the same tree/threshold/AdQSM variant. Without
# tagging results with which resampling setting produced them, re-running
# this script after tweaking SEG_LEN_* would silently OVERWRITE the
# previous run's row in RESULTS_CSV (upsert_result() replaces any existing
# row with the same tree+method - see tree_geom_utils.py) and its .npz file
# on disk - even though the two runs describe genuinely different
# geometries, not a correction of the same one. SEG_VARIANT_SUFFIX is
# appended to every AdTree method name AND output filename further down, so
# different resampling settings coexist side by side instead of clobbering
# each other.
#
# Format: "_seg{min_mm}-{max_mm}-k{k*100}". Millimetres (not metres) and
# k*100 (not the raw 0..1 fraction) are used purely to keep the suffix a
# short string of whole numbers with no decimal points, which would
# otherwise need extra escaping to stay safe inside both filenames and CSV
# cells. Example: SEG_LEN_MIN=0.01, SEG_LEN_MAX=0.3, SEG_LEN_K=0.5 ->
# "_seg10-300-k50".
#
# Now computed INSIDE the sweep loop in the RUN section below (not here at
# module level), once per (SEG_LEN_MIN, SEG_LEN_K) combination, using the
# same format string plus a "-wlen" tag when RADIUS_STAT_WEIGHTING ==
# "length" (see that parameter, further below) - "" (no tag) when "none",
# keeping this exact format string byte-identical in that default case.

# Shortest permissible cylinder (in METERS). Shorter ones (e.g. at branch
# points) are not created, to avoid degenerate zero-length beams in ANSYS.
# At the end of the resampling step, any cylinder shorter than this is merged into its parent cylinder. Set to 0 or None to disable this check entirely.
MIN_CYL_LEN = 0.005               # 0.1 mm

# Rounding used when merging coincident points (number of decimal places).
# 5 is a safe default.
MERGE_DECIMALS = 5

# Shift x,y to zero (removes a large georeferencing offset)? Recommended True.
# NOTE: the actual shift is applied later, by write_geom() in
# export_geom_ansys.py - this script only carries the flag through into the
# .npz so that step 2 shifts (or doesn't) exactly the way this run intended.
RECENTER_XY = True

# --- Centerline smoothing (applied BEFORE pruning/resampling) --------
# Laplacian smoothing removes zig-zag noise from the AdTree skeleton. Each
# smoothing pass moves every free point along a branch toward the average of
# its two neighbours by SMOOTH_ALPHA. The root, junctions (degree >= 3), and
# branch tips (leaves) are never moved, so topology and branch endpoints stay
# fixed. Radii are not affected.
# TODO: Sensitivity analysis of SMOOTH_ITERS/SMOOTH_ALPHA - how much smoothing is too much?
SMOOTH_ITERS = 5                  # number of smoothing passes (0 = off)
SMOOTH_ALPHA = 0.5                # 0..1 strength per pass

# --- Radius calibration against AdQSM (optional) ----------------------
# Replaces the AdTree skeleton radii with AdQSM-calibrated radii on the final
# (pruned + resampled) cylinders. The trunk (branch order 0) is taken from the
# AdQSM taper curve; every other branch order is scaled by a single factor
# (AdQSM median radius / AdTree median radius) computed for that order. Set
# CALIBRATE_RADII = False to skip this step entirely and use the raw AdTree
# radii, exactly like before this feature existed.
CALIBRATE_RADII = True

# RADIUS_STAT_WEIGHTING_LIST: each entry is "none" | "length", passed
# straight through to build_quantile_matched_pairs() and
# compute_order_calibration_factors() (tree_geom_utils.py - see both for the
# full explanation). "none" is BYTE-IDENTICAL to this project's behaviour
# before this parameter existed: each AdTree cylinder counts as one sample
# in the per-order quantile/median regardless of its length. "length"
# instead weights each cylinder by its own length, so a densely-resampled
# region (many short cylinders) doesn't over-represent itself relative to a
# sparsely-resampled one covering the same physical branch length - see
# CHANGELOG_adtree.md for the investigation that motivated this (order-1
# median AdTree radius shifted 39% across two SEG_LEN settings when
# unweighted, 0.15% when length-weighted, while total length/volume per
# order agreed to within 1% either way). Only the AdTree side is ever
# weighted; see the two functions' own docstrings for why the AdQSM side is
# not.
#
# A LIST, not a single value: the unweighted variant is being kept as a
# comparison baseline while further beech trees are processed, so every run
# needs to compute BOTH variants side by side (two rows per method/threshold
# in RESULTS_CSV, distinguished by SEG_VARIANT_SUFFIX's "-wlen" tag - see
# below), rather than requiring two separate runs (one per weighting) per
# tree. The RUN section loops over this list at the point where the fit is
# computed and applied - NOT around convert() - since geometry does not
# depend on weighting; see the "RADIUS_STAT_WEIGHTING_LIST loop" comment
# further down for exactly where and why.
RADIUS_STAT_WEIGHTING_LIST = ["none", "length"]

# Measured trunk diameter at breast height (1.3 m), in METERS. If given, the
# taper curve is rescaled so its value at 1.3 m matches this measurement.
# Set to None to use the taper curve exactly as read from ADQSM_TAPER_FILE.
FIELD_DBH = None

# --- Thin-branch (< 10 cm) volume diagnostic ---------------------------
# The de Tanago field reference only measured branches down to a 10 cm
# TAPER diameter (see AdQSM.pdf Appendix A) - the same cut-off already
# applied to TreeQSM in runsken.m (section 17b, "Cylinders, cut-off 10 cm").
# This block prints, for BOTH AdTree calibrated cylinders and the raw
# AdQSM BranchStructure.txt, how much volume (or - for AdQSM, if no branch
# length is available - just count/share) sits below this cut-off. The
# point is to check whether AdQSM/AdTree calibrated agreeing closely with
# the reference is a REAL result, or just because they structurally put
# little volume into branches that thin anyway (unlike TreeQSM, which
# reconstructs all of them, thin or not).
THIN_BRANCH_CUT_CM = 10.0

# Print the first few raw rows of AdQSM's BranchStructure.txt (with column
# names/indices) before using it, so you can SEE its actual columns instead
# of trusting the parsing code blindly. Purely diagnostic; safe to turn off
# once you've checked your file's layout.
PRINT_ADQSM_BRANCH_SAMPLE = True

# Also write a second CSV row per (radius threshold, AdQSM variant) with the
# AdTree-calibrated volume restricted to cylinders >= THIN_BRANCH_CUT_CM
# (similar in spirit to TreeQSM's "...Filtered..." rows), so it can be
# compared directly in compare_volumes.py / plot_volumes.py. Set False to skip.
WRITE_THIN_BRANCH_FILTERED_ROW = False


# Shared master results table (see compare_volumes.py). When CALIBRATE_RADII
# is True, each generated threshold variant upserts its own row into this CSV
# so results from all methods live in one place.
RESULTS_CSV = "volume_results.csv"

# Output folders, so the working directory doesn't fill up with dozens of
# .npz/.png files mixed in with the scripts. NPZ_DIR is this script's own
# new folder; FIGURES_DIR reuses the project's shared PLOTS_DIR convention
# (paths.py - already used by plot_box.py/plot_volumes.py for their own
# charts), grouped under a per-tree subfolder. export_geom_ansys.py has its
# OWN matching NPZ_DIR parameter (see that file) - keep both in sync by
# hand if this one ever changes.
NPZ_DIR = "npz"
FIGURES_DIR = os.path.join(PLOTS_DIR, TREE_NAME)

# Reference heights [m] used for DBH (lower) and the taper metric (lower/
# upper). DBH is the stem diameter at TAPER_H_LOWER (1.3 m = breast height).
TAPER_H_LOWER = 1.3    # lower reference height [m]
TAPER_H_UPPER = 10.0   # upper reference height [m]
# =====================================================================

# --- Build the list of AdQSM variants to actually process --------------
# You should NOT need to edit this block - edit AdQSM_DIR (single variant)
# or ADQSM_BASE_DIR + ADQSM_VARIANTS (several variants) above instead.
#
# The RUN section further below loops over ADQSM_VARIANT_LIST. Each entry
# is a tuple (variant_label, taper_file, branch_file, params_file):
#   - variant_label is None when you used the simple single-AdQSM_DIR case
#     (case 1 above) - in that case output filenames/method names get NO
#     extra suffix, exactly like before this feature existed.
#   - variant_label is the variant's subfolder name (e.g. "05") when you
#     used ADQSM_VARIANTS (case 2) - in that case it's appended to output
#     filenames and to method names in RESULTS_CSV, so the variants don't
#     overwrite each other and stay distinguishable in the results table.
if ADQSM_VARIANTS:
    ADQSM_VARIANT_LIST = []
    for variant_name in ADQSM_VARIANTS:
        variant_dir = os.path.join(ADQSM_BASE_DIR, variant_name)
        ADQSM_VARIANT_LIST.append((
            variant_name,
            os.path.join(variant_dir, "taper.txt"),            # taper curve
            os.path.join(variant_dir, "BranchStructure.txt"),  # per-branch table
            os.path.join(variant_dir, "TreesParams.txt"),      # whole-tree params
        ))
else:
    ADQSM_VARIANT_LIST = [(
        None,
        os.path.join(AdQSM_DIR, "taper.txt"),
        os.path.join(AdQSM_DIR, "BranchStructure.txt"),
        os.path.join(AdQSM_DIR, "TreesParams.txt"),
    )]
# =====================================================================

# --- Output file name(s) --------------------------------------------
# The FINAL geom_*.txt name that export_geom_ansys.py will eventually write
# (this script itself only writes the intermediate .npz - see the RUN
# section) is no longer a fixed/hand-edited name here. Instead it's built
# automatically from the same ingredients as npz_name below (TREE_NAME,
# threshold, variant_suffix, SEG_VARIANT_SUFFIX), just with "geom_" instead
# of "calib_" and ".txt" instead of ".npz" - see the `out = ...` line in the
# RUN section. This guarantees the geom_*.txt name always matches the
# calib_*.npz it came from (e.g. calib_IND01_054_r5mm_seg100-500-k50.npz ->
# geom_IND01_054_r5mm_seg100-500-k50.txt), so you can tell at a glance which
# .npz produced which geom_*.txt, and different thresholds/variants/segment
# settings never silently overwrite each other's exported file.
# =====================================================================

# =====================  VISUALIZATION  ================================
# Show an interactive 3D preview of the reduced beam model (one window per
# radius threshold) so you can check it before importing into ANSYS.
SHOW_PLOT = False

# Also save the preview as a PNG next to each output file (True/False).
# The PNG name is derived from the output file name (.txt -> .png).
SAVE_PLOT_PNG = True
# =====================================================================


# =========================  RUN  =====================================
if PRINT_TIMING:
    _run_start = time.perf_counter()

os.makedirs(NPZ_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

if PRINT_TIMING:
    _prep_start = time.perf_counter()

print("Reading:", INPUT_PLY)
xyz, rad, edges = read_ply(INPUT_PLY)
print("  vertices: %d, edges: %d" % (len(xyz), len(edges)))

xyz, rad, edges = merge_vertices(xyz, rad, edges, MERGE_DECIMALS)
ncomp, _ = connected_components(
    csr_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
               shape=(len(xyz), len(xyz))), directed=False)
print("  after merging: %d points, %d edges, connected components: %d" % (len(xyz), len(edges), ncomp))
xyz_raw = xyz.copy()   # AdTree geometry/radii before smoothing, kept as the calibration baseline

smooth_root = int(np.argmin(xyz[:, 2]))
xyz = smooth_centerline(xyz, edges, smooth_root, SMOOTH_ITERS, SMOOTH_ALPHA)
print("  centerline smoothing: %d passes, alpha=%.2f" % (SMOOTH_ITERS, SMOOTH_ALPHA))

if PRINT_TIMING:
    print("[TIMING] read_ply + merge_vertices + smoothing: %.2f s" % (time.perf_counter() - _prep_start))

print("\n%-12s %-12s %-12s %-12s" % ("threshold", "cylinders", "length [m]", "file"))
print("Volume verification below is computed from the exact cylinders that will be "
      "written to each geom file (pi * radius^2 * length per cylinder).\n")
z_base = float(xyz[:, 2].min())   # tree base; DBH/height/taper are measured from here
multiple_variants = len(ADQSM_VARIANT_LIST) > 1   # True only if you used ADQSM_VARIANTS (case 2 above)

# --- SEG_LEN sweep bookkeeping ------------------------------------------
# Unconditional (NOT gated by PRINT_TIMING): counts how many
# (AdQSM variant, SEG_LEN_MIN, SEG_LEN_K) combinations this run actually
# executes, and how long the whole sweep takes, for the end-of-run summary
# printed after the WHOLE RUN line below.
_n_seg_combinations = 0
_sweep_run_start = time.perf_counter()

# Outer loop: one pass per AdQSM variant (just one pass, using the plain
# AdQSM_DIR, unless you filled in ADQSM_VARIANTS). Everything inside this
# loop (loading AdQSM data, calibrating, saving results) is repeated once
# per variant, so several reconstructions can be compared side by side in
# the same RESULTS_CSV without overwriting each other.
for variant_label, taper_file, branch_file, params_file in ADQSM_VARIANT_LIST:
    if PRINT_TIMING:
        _variant_start = time.perf_counter()

    # variant_suffix/variant_method_suffix are "" when there is only one
    # variant (so filenames/method names look exactly like before this
    # feature existed); otherwise they tag the variant name onto them.
    variant_suffix = ("_adqsm%s" % variant_label) if variant_label else ""
    variant_method_suffix = (" (AdQSM %s)" % variant_label) if variant_label else ""

    if CALIBRATE_RADII:
        print("\nLoading AdQSM calibration data%s..."
              % (" (variant: %s)" % variant_label if variant_label else ""))
        taper_heights, taper_diameters = parse_adqsm_taper_file(taper_file)
        # trunk_radius_func: every calibrated trunk (order 0) cylinder's
        # radius comes ENTIRELY from THIS curve, at the AdTree cylinder's
        # own height - never from AdTree's own measured radius. See
        # apply_order_calibration_factors()/apply_radius_regression_per_order()
        # in tree_geom_utils.py, both of which look up trunk radius here
        # regardless of calmethod. That means a bad taper.txt directly and
        # fully determines a bad calibrated trunk volume, with no other
        # signal to catch it - _reject_taper_spikes() (called from inside
        # parse_adqsm_taper_file() above) only drops an isolated single-row
        # spike, not a whole implausible SECTION of the curve.
        #
        # Confirmed case: data/IND07_083/04/taper.txt reports an almost
        # exactly constant diameter (~0.708 m) from 0-22.6 m height - real
        # trunk taper should decrease measurably over 22+ metres, so this
        # is not physically plausible - then turns erratic from 23.6-27.6 m
        # before the already-handled 28.6 m spike. This inflates that
        # tree's calibrated trunk volume well past AdQSM's own official
        # TrunkVolume (TreesParams.txt) for the same tree - identically for
        # both calmethod=min5mm and calmethod=regression-perorder, since
        # neither reads AdTree's own trunk radius at all. See
        # TODO_investigations.md (item 8) for the follow-up: checking other
        # trees' taper.txt files for the same kind of implausibly flat or
        # erratic section before trusting their calibrated trunk volumes.
        trunk_radius_func = make_trunk_radius_func(taper_heights, taper_diameters, FIELD_DBH)
        adqsm_median_by_order = parse_adqsm_branch_file(branch_file)
        # raw_diam_by_order (Task: regression calibration method): the RAW
        # per-order AdQSM radius lists behind adqsm_median_by_order's medians
        # - needed for build_quantile_matched_pairs() below, which matches
        # whole DISTRIBUTIONS rather than single median points.
        raw_diam_by_order = parse_adqsm_branch_file_raw(branch_file)
        print("  %s: %d height/diameter rows, %.1f-%.1f m"
              % (taper_file, len(taper_heights), taper_heights.min(), taper_heights.max()))
        if FIELD_DBH is not None:
            print("  taper curve rescaled so radius at 1.3 m = FIELD_DBH/2 = %.4f m" % (FIELD_DBH / 2.0))
        print("  %s: AdQSM median radius by order: %s"
              % (branch_file,
                 ", ".join("%d=%.4f m" % (o, r) for o, r in sorted(adqsm_median_by_order.items()))))

        raw_stats = raw_skeleton_stats(xyz_raw, rad, edges)
        print("  raw PLY skeleton baseline (all tree edges, AdTree radii):")
        print_volume_stats("raw skeleton (AdTree)", raw_stats)

        # ---- (task 1) upsert the AdQSM reference itself into RESULTS_CSV ----
        # This is AdQSM's OWN reported volume (straight from TreesParams.txt),
        # not anything derived from the AdTree skeleton - it does not depend on
        # RADIUS_THRESHOLDS, so it's only written once per variant (here),
        # outside the threshold loop below.
        adqsm_ref = parse_adqsm_params_file(params_file)
        if adqsm_ref is not None:
            print_volume_stats("(d) AdQSM reference (TreesParams)", adqsm_ref)
            if adqsm_ref.get("height") is not None:
                print("      AdQSM TreeHeight: %.2f m" % adqsm_ref["height"])
            # branch_filter = "none": AdQSM's own TreesParams.txt totals are its
            # full reconstruction, not restricted to any diameter cut-off.
            upsert_result(RESULTS_CSV, TREE_NAME,
                          "AdQSM (TreesParams)%s" % variant_method_suffix,
                          adqsm_ref["total_vol"], adqsm_ref["trunk_vol"], adqsm_ref["branch_vol"], None,
                          adqsm_ref.get("dbh"), adqsm_ref.get("height"), None,
                          # trunk_len/branch_len: already in this dict, straight from
                          # TrunkLength/BranchLength in TreesParams.txt (see parse_adqsm_params_file).
                          adqsm_ref.get("trunk_len"), adqsm_ref.get("branch_len"),
                          branch_filter="none")
        else:
            print("  (no TreesParams.txt reference found at %s - skipping that row)" % params_file)

        # ---- (Part B) thin-branch diagnostic straight from BranchStructure.txt ----
        if PRINT_ADQSM_BRANCH_SAMPLE:
            print("\n  BranchStructure.txt column check (first 10 rows):")
            print_adqsm_branch_file_sample(branch_file, n=10)
        # NOTE: report_adqsm_thin_branch() returns its >=cut_cm-kept cylinder-
        # approximation totals (see tree_geom_utils.py), but they are
        # DELIBERATELY NOT upserted into RESULTS_CSV any more - REVERTED
        # after the approximation was shown to badly overestimate: for all
        # three AdQSM variants tested, the >=10cm-filtered approximation
        # exceeded even AdQSM's own OFFICIAL, UNFILTERED whole-tree
        # BranchVolume from TreesParams.txt - a logically impossible result
        # (a filtered subset can never exceed its own unfiltered total), so
        # it cannot be trusted as a comparison row. See CHANGELOG_adtree.md.
        # The return value is still printed to console (unchanged) for
        # manual reference; only the upsert_result() call that used to write
        # an "AdQSM (BranchStructure, cyl. approx.)" row was removed.
        if WRITE_THIN_BRANCH_FILTERED_ROW:
            report_adqsm_thin_branch(branch_file, cut_cm=THIN_BRANCH_CUT_CM, params_file=params_file)

    # ---- SEG_LEN sweep: run every (SEG_LEN_MIN, SEG_LEN_K) combination ----
    # Nesting choice: this sweep sits INSIDE the AdQSM-variant loop but
    # AFTER the AdQSM-data-loading block above (taper/branch/params file
    # parsing, the AdQSM-reference upsert, the thin-branch sample print) -
    # none of that depends on SEG_LEN_MIN/SEG_LEN_K (see PHASE A survey,
    # finding A4), so it stays outside the sweep and runs ONCE per AdQSM
    # variant, not once per sweep combination - the same per-variant cost
    # as before this feature existed. Everything from here down to the end
    # of the RADIUS_THRESHOLDS loop DOES depend on SEG_LEN_MIN/SEG_LEN_K
    # (ref_root_0/ref_cyl_0/ref_order_0, factors_by_ref, and order_to_ab all
    # do - see A4), so all of it is wrapped in the sweep, once per
    # combination.
    for SEG_LEN_MIN in SEG_LEN_MIN_LIST:
        for SEG_LEN_K in SEG_LEN_K_LIST:
            # SEG_VARIANT_SUFFIX_BASE: same format string as before this
            # feature existed (see the PARAMETERS block comment) - computed
            # here, once per (SEG_LEN_MIN, SEG_LEN_K) combination. Carries NO
            # weighting information - the per-weighting "-wlen" tag is
            # appended further below, once per RADIUS_STAT_WEIGHTING_LIST
            # entry, giving each weighting its own SEG_VARIANT_SUFFIX.
            SEG_VARIANT_SUFFIX_BASE = "_seg%d-%d-k%d" % (
                round(SEG_LEN_MIN * 1000),
                round(SEG_LEN_MAX * 1000),
                round(SEG_LEN_K * 100),
            )

            # _WEIGHTING_ITER: what the two "RADIUS_STAT_WEIGHTING_LIST loop"
            # blocks below (fit COMPUTATION here, fit APPLICATION inside the
            # RADIUS_THRESHOLDS loop) iterate over. When CALIBRATE_RADII is
            # False there is no fit at all - raw AdTree cylinders are
            # exported as-is, exactly as before RADIUS_STAT_WEIGHTING_LIST
            # existed - so this is the single-element [None] sentinel
            # instead of the real list, keeping that path's SEG_VARIANT_SUFFIX
            # (= SEG_VARIANT_SUFFIX_BASE, no tag), single export, and absence
            # of upserts byte-identical to before.
            _WEIGHTING_ITER = RADIUS_STAT_WEIGHTING_LIST if CALIBRATE_RADII else [None]

            _n_seg_combinations += 1
            print("[seg sweep] variant=%s  SEG_LEN_MIN=%.3f  SEG_LEN_K=%.2f  base_suffix=%s"
                  % (variant_label if variant_label else "(single)", SEG_LEN_MIN, SEG_LEN_K, SEG_VARIANT_SUFFIX_BASE))

            if CALIBRATE_RADII:

                # ---- FIXED calibration reference set: calref=min5mm (secondary) ---
                # Originally diagnosed: calibrate_cylinder_radii()'s self-referencing
                # behaviour computed each order's AdTree median from the SAME,
                # already-pruned cylinder set it was calibrating, so a higher
                # RADIUS_THRESHOLDS value mechanically shrank the factor and
                # over-rescaled even the thick, never-pruned cylinders of that order
                # - see CHANGELOG_adtree.md (Steps 1-3) for the full investigation,
                # including calref=unpruned and calref=min2/3/4mm, since removed.
                #
                # DECISION (CHANGELOG_adtree.md, Step 7): per-order regression
                # ([calmethod=regression-perorder], below) is the PRIMARY calibration
                # method going forward. This fixed 5mm-reference factor set
                # (calref=min5mm) is kept only as a SECONDARY/backup reference point
                # - computed ONCE per AdQSM variant (not inside the RADIUS_THRESHOLDS
                # loop below), reused for every threshold.
                #
                # ref_cyl_0/ref_order_0 (the fully unpruned reference set, thr=0.0)
                # is kept too - NOT for calref (calref=unpruned was removed), but
                # because build_quantile_matched_pairs() below (for the per-order
                # regression) still needs it as its fixed AdTree reference population.
                #
                # Both convert() calls below (unpruned reference, calref
                # reference) run EXACTLY ONCE per (AdQSM variant, SEG_LEN
                # combination) regardless of how many entries
                # RADIUS_STAT_WEIGHTING_LIST has - geometry does not depend on
                # weighting (see that parameter's own PARAMETERS-block
                # comment), only the STATISTIC computed over it does. The
                # RADIUS_STAT_WEIGHTING_LIST loop sits AFTER both convert()
                # calls, around compute_order_calibration_factors()/
                # build_quantile_matched_pairs()/the regression fit instead -
                # see that loop's own comment below.
                if PRINT_TIMING:
                    _t0 = time.perf_counter()
                ref_root_0, ref_cyl_0, ref_order_0 = convert(
                    xyz, rad, edges, 0.0, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
                if PRINT_TIMING:
                    print("  [TIMING] convert() [unpruned reference, thr=0.0]: %.2f s" % (time.perf_counter() - _t0))

                # ref_geom_by_thr: {ref_thr: (ref_root, ref_cyl, ref_order)} -
                # the calref reference GEOMETRY, one convert() call per
                # CALIBRATION_REF_THRESHOLDS_MM entry (just min5mm), computed
                # ONCE here, weighting-independent. compute_order_calibration_
                # factors() (which DOES depend on weighting) is called once
                # per weighting further below, reusing this cached geometry
                # instead of re-running convert().
                ref_geom_by_thr = {}
                if COMPUTE_CALREF_MIN5MM:
                    for ref_thr in CALIBRATION_REF_THRESHOLDS_MM:
                        if PRINT_TIMING:
                            _t0 = time.perf_counter()
                        ref_geom_by_thr[ref_thr] = convert(
                            xyz, rad, edges, ref_thr, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
                        if PRINT_TIMING:
                            print("  [TIMING] convert() [calref thr=%.3f]: %.2f s" % (ref_thr, time.perf_counter() - _t0))

                # ---- RADIUS_STAT_WEIGHTING_LIST loop: fit COMPUTATION -------------
                # Placed HERE - around compute_order_calibration_factors()/
                # build_quantile_matched_pairs()/the per-order regression fit -
                # NOT around either convert() call above: the reference
                # geometry (ref_cyl_0/ref_order_0, ref_geom_by_thr) is fixed
                # once per SEG_LEN combination regardless of weighting; only
                # the per-order STATISTIC each weighting computes over that
                # same fixed geometry differs. Results are cached per
                # weighting (factors_by_ref_by_weighting, order_to_ab_by_
                # weighting, seg_variant_suffix_by_weighting) for the
                # RADIUS_THRESHOLDS loop below, which APPLIES them - see that
                # loop's matching RADIUS_STAT_WEIGHTING_LIST loop.
                factors_by_ref_by_weighting = {}
                order_to_ab_by_weighting = {}
                seg_variant_suffix_by_weighting = {}
                for RADIUS_STAT_WEIGHTING in _WEIGHTING_ITER:
                    # SEG_VARIANT_SUFFIX: SEG_VARIANT_SUFFIX_BASE plus a
                    # "-wlen" tag when RADIUS_STAT_WEIGHTING == "length" - ""
                    # (no tag) when "none", so that case's suffix/filenames/
                    # method-names stay BYTE-IDENTICAL to before this
                    # parameter existed. A length-weighted run cannot collide
                    # with (or silently overwrite) an unweighted run's row in
                    # RESULTS_CSV or its exported calib_*.npz/geom_*.txt -
                    # both are named from this exact suffix (see `out`/
                    # `npz_name` in the RADIUS_THRESHOLDS loop below), and the
                    # FINAL exported cylinders genuinely differ between the
                    # two weightings, so this is not just a label collision,
                    # it would be two different geometries fighting over one
                    # filename. Chosen over tagging only the calmethod label
                    # because the calmethod tag alone never reaches the
                    # geom_*.txt/calib_*.npz filenames (those are keyed on
                    # SEG_VARIANT_SUFFIX, not calmethod) - extending
                    # SEG_VARIANT_SUFFIX fixes both at once.
                    SEG_VARIANT_SUFFIX = SEG_VARIANT_SUFFIX_BASE + (
                        "-wlen" if RADIUS_STAT_WEIGHTING == "length" else "")
                    seg_variant_suffix_by_weighting[RADIUS_STAT_WEIGHTING] = SEG_VARIANT_SUFFIX
                    print("  [weighting=%s] suffix=%s" % (RADIUS_STAT_WEIGHTING, SEG_VARIANT_SUFFIX))

                    # factors_by_ref: {ref_thr: factors_dict} - one fixed factors dict per
                    # CALIBRATION_REF_THRESHOLDS_MM entry (just min5mm now), computed
                    # ONCE here per weighting (not inside the RADIUS_THRESHOLDS loop below).
                    factors_by_ref = {}
                    if COMPUTE_CALREF_MIN5MM:
                        for ref_thr in CALIBRATION_REF_THRESHOLDS_MM:
                            ref_root, ref_cyl, ref_order = ref_geom_by_thr[ref_thr]
                            factors_by_ref[ref_thr] = compute_order_calibration_factors(
                                ref_cyl, ref_order, adqsm_median_by_order,
                                xyz=xyz, weighting=RADIUS_STAT_WEIGHTING)
                            print("  Fixed calibration factors - reference set 'min%dmm' (thr=%.3f), "
                                  "reused for every RADIUS_THRESHOLDS run:" % (round(ref_thr * 1000), ref_thr))
                            for o in sorted(factors_by_ref[ref_thr]):
                                print("    order %d : factor = %.3f" % (o, factors_by_ref[ref_thr][o]))
                    factors_by_ref_by_weighting[RADIUS_STAT_WEIGHTING] = factors_by_ref

                    # ---- PRIMARY calibration method: per-order (grouped) regression ---
                    # Adopted as the primary calibration method (CHANGELOG_adtree.md,
                    # Step 7), after the investigation found order-dependent bias a
                    # single global fit could not capture (order 1's own ratio ~1.6 vs.
                    # ~2.1-2.35 for every other order). group_orders_for_fitting() merges
                    # sparse orders together (walking ascending, greedy upward merge) so
                    # every group still has >= MIN_PAIRS_PER_ORDER pairs for a stable
                    # two-parameter fit; each group then gets its own (a, b) via
                    # fit_radius_regression() on that group's own pooled quantile-matched
                    # pairs (build_quantile_matched_pairs(), reusing the fixed
                    # ref_cyl_0/ref_order_0 reference population above). Computed ONCE
                    # per AdQSM variant/weighting here, reused for every RADIUS_THRESHOLDS
                    # value below.
                    adtree_matched, adqsm_matched, order_labels_matched = build_quantile_matched_pairs(
                        ref_cyl_0, ref_order_0, raw_diam_by_order,
                        xyz=xyz, weighting=RADIUS_STAT_WEIGHTING)

                    order_to_group = group_orders_for_fitting(order_labels_matched, MIN_PAIRS_PER_ORDER)

                    group_fits = []    # [(group_orders_tuple, a, b), ...] - for the diagnostic plot
                    order_to_ab = {}   # {order: (a, b)} - for apply_radius_regression_per_order()
                    for group_orders in dict.fromkeys(order_to_group.values()):   # de-duplicated, first-seen order
                        group_mask = np.isin(order_labels_matched, list(group_orders))
                        print("  Per-order regression: fitting group orders=%s (n_pairs=%d)..."
                              % (str(group_orders), int(group_mask.sum())))
                        g_a, g_b = fit_radius_regression(adtree_matched[group_mask], adqsm_matched[group_mask])
                        group_fits.append((group_orders, g_a, g_b))
                        for o in group_orders:
                            order_to_ab[o] = (g_a, g_b)
                    order_to_ab_by_weighting[RADIUS_STAT_WEIGHTING] = order_to_ab

                    # order1_merge_note: short text for plot_radius_regression_per_order()'s
                    # on-plot annotation (bottom-left corner) - group_orders_for_fitting()
                    # already printed the loud console warning above when this applies;
                    # this just makes the same fact visible on the PNG itself.
                    order1_merge_note = None
                    group_of_1 = order_to_group.get(1)
                    if group_of_1 is not None and len(group_of_1) > 1:
                        order1_merge_note = ("order 1 MERGED with order(s) %s\n(see CHANGELOG_adtree.md)"
                                              % [o for o in group_of_1 if o != 1])

                    regression_perorder_plot_path = plot_radius_regression_per_order(
                        adtree_matched, adqsm_matched, order_labels_matched, group_fits,
                        TREE_NAME, variant_label, order1_merge_note=order1_merge_note,
                        filename_suffix=SEG_VARIANT_SUFFIX, plots_dir=FIGURES_DIR)
                    print("  Saved per-order regression diagnostic plot: %s" % regression_perorder_plot_path)
            else:
                seg_variant_suffix_by_weighting = {None: SEG_VARIANT_SUFFIX_BASE}

            # Inner loop: one pass per radius threshold (same as before this feature
            # existed), now repeated for each AdQSM variant AND each SEG_LEN sweep
            # combination above.
            for thr in RADIUS_THRESHOLDS:
                if PRINT_TIMING:
                    _thr_start = time.perf_counter()
                    _t0 = time.perf_counter()
                root, cyl, cyl_order = convert(xyz, rad, edges, thr, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
                if PRINT_TIMING:
                    print("  [TIMING] convert() [RADIUS_THRESHOLDS thr=%.3f]: %.2f s" % (thr, time.perf_counter() - _t0))
                # This single convert() call above is the ONLY per-threshold
                # geometry computation - it runs exactly once here regardless
                # of how many RADIUS_STAT_WEIGHTING_LIST entries there are;
                # `cyl`/`cyl_order` are read (never mutated) by every
                # weighting pass in the RADIUS_STAT_WEIGHTING_LIST loop below.

                # Height of the pruned model: z-range of the nodes actually used by these
                # cylinders. Unaffected by radius calibration (geometry doesn't change).
                node_ids = sorted({idx for a, b, r, pid in cyl for idx in (a, b)})
                height_m = float(xyz[node_ids, 2].max() - xyz[node_ids, 2].min()) if node_ids else None

                if CALIBRATE_RADII:
                    # orig_stats/raw_dbh/raw_taper/orig_thin describe the RAW
                    # (uncalibrated) cylinders straight from convert() above -
                    # weighting-independent (weighting only affects the FIT,
                    # not the raw AdTree geometry), so computed ONCE per
                    # threshold, outside the RADIUS_STAT_WEIGHTING_LIST loop
                    # below, instead of once per weighting.
                    orig_lengths, orig_radii = cylinder_metrics(xyz, cyl)
                    orig_stats = volume_stats(orig_lengths, orig_radii, np.asarray(cyl_order))
                    # DBH/taper of the UNCALIBRATED (raw AdTree) trunk, before radii are replaced.
                    raw_dbh = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_LOWER)
                    raw_d_upper = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_UPPER)
                    raw_taper = ((raw_dbh - raw_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                                 if raw_dbh is not None and raw_d_upper is not None else None)

                    # ---- thin-branch diagnostic on the RAW (uncalibrated) cylinders ----
                    # source_label="AdTree raw" makes this printout visually distinct
                    # from the calibrated one further below (same function, same cut_cm,
                    # different cylinder set) - see report_thin_branch_volume()'s
                    # docstring in tree_geom_utils.py for why the label exists.
                    if WRITE_THIN_BRANCH_FILTERED_ROW:
                        orig_thin = report_thin_branch_volume(orig_lengths, orig_radii, cyl_order,
                                                               cut_cm=THIN_BRANCH_CUT_CM, source_label="AdTree raw")
                    else:
                        orig_thin = None

                # ---- RADIUS_STAT_WEIGHTING_LIST loop: fit APPLICATION -------------
                # Same list (via _WEIGHTING_ITER) as the fit-COMPUTATION loop
                # above - this pass APPLIES each weighting's cached
                # factors_by_ref/order_to_ab to THIS threshold's raw
                # `cyl`/`cyl_order` (computed once, just above, OUTSIDE this
                # loop) and exports/upserts the result. `cyl`/`cyl_order` are
                # only READ here, never reassigned, so the SAME raw geometry
                # is reused unmutated by every weighting pass; each pass's
                # calibrated cylinders go into a fresh local `final_cyl`
                # instead of overwriting `cyl` (the old single-weighting code
                # used to reassign `cyl = regperorder_cyl` here - now that
                # would corrupt the next weighting pass's "raw" input).
                for RADIUS_STAT_WEIGHTING in _WEIGHTING_ITER:
                    SEG_VARIANT_SUFFIX = seg_variant_suffix_by_weighting[RADIUS_STAT_WEIGHTING]
                    if CALIBRATE_RADII:
                        print("  [weighting=%s]" % RADIUS_STAT_WEIGHTING)

                    # `out` is the geom_*.txt name step 2 (export_geom_ansys.py) will
                    # eventually write - computed here (once per weighting, alongside
                    # the threshold/variant it belongs to) and carried inside the .npz
                    # below.
                    #
                    # Built with the EXACT SAME ingredients (and in the same order) as
                    # npz_name further down - just "geom_"/".txt" instead of
                    # "calib_"/".npz" - so every geom_*.txt name matches the calib_*.npz
                    # it was exported from at a glance, e.g.:
                    #   calib_IND01_054_r5mm_seg100-500-k50.npz
                    #   geom_IND01_054_r5mm_seg100-500-k50.txt
                    # This replaces the old fixed OUTPUT_NAME/OUTPUT_PATTERN constants -
                    # every combination of tree/threshold/variant/segment-settings/
                    # weighting now gets its own name automatically, so nothing can
                    # silently overwrite a previous run's exported file.
                    out = "geom_%s_r%dmm%s%s.txt" % (TREE_NAME, round(thr * 1000), variant_suffix, SEG_VARIANT_SUFFIX)

                    # final_cyl: the raw AdTree cylinders by default (used as-is when
                    # CALIBRATE_RADII is False - identical to `cyl`, matching the old
                    # behaviour exactly), reassigned to the calibrated regperorder_cyl
                    # below when CALIBRATE_RADII is True.
                    final_cyl = cyl

                    if CALIBRATE_RADII:
                        factors_by_ref = factors_by_ref_by_weighting[RADIUS_STAT_WEIGHTING]
                        order_to_ab = order_to_ab_by_weighting[RADIUS_STAT_WEIGHTING]

                        # ---- SECONDARY calibration variant: calref=min5mm -----------------
                        # Apply THIS weighting's FIXED factors_by_ref[...] dict (computed
                        # once per AdQSM variant/weighting, above the RADIUS_THRESHOLDS
                        # loop, from the fixed min5mm reference cylinder set) to THIS
                        # threshold's still-RAW `cyl`/`cyl_order`. Kept as a
                        # secondary/backup reference point alongside the primary
                        # per-order regression rows below (see CHANGELOG_adtree.md,
                        # Step 7).
                        fixedref_data = {}
                        if COMPUTE_CALREF_MIN5MM:
                            fixedref_variants = [
                                ("min%dmm" % round(ref_thr * 1000), factors_by_ref[ref_thr])
                                for ref_thr in CALIBRATION_REF_THRESHOLDS_MM
                            ]
                            for ref_name, ref_factors in fixedref_variants:
                                fr_new_r = apply_order_calibration_factors(
                                    xyz, cyl, cyl_order, trunk_radius_func, ref_factors)
                                fr_cyl = [(a, b, float(fr_new_r[i]), pid) for i, (a, b, r, pid) in enumerate(cyl)]
                                fr_lengths, fr_radii = cylinder_metrics(xyz, fr_cyl)
                                fr_stats = volume_stats(fr_lengths, fr_radii, np.asarray(cyl_order))
                                if WRITE_THIN_BRANCH_FILTERED_ROW:
                                    fr_thin = report_thin_branch_volume(
                                        fr_lengths, fr_radii, cyl_order, cut_cm=THIN_BRANCH_CUT_CM,
                                        source_label="AdTree calibrated [calref=%s]" % ref_name)
                                else:
                                    fr_thin = None
                                fr_dbh = stem_diameter_at_height(xyz, fr_cyl, cyl_order, z_base, TAPER_H_LOWER)
                                fr_d_upper = stem_diameter_at_height(xyz, fr_cyl, cyl_order, z_base, TAPER_H_UPPER)
                                fr_taper = ((fr_dbh - fr_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                                            if fr_dbh is not None and fr_d_upper is not None else None)
                                fixedref_data[ref_name] = dict(stats=fr_stats, thin=fr_thin, dbh=fr_dbh,
                                                                taper=fr_taper, n_cylinders=len(fr_cyl))

                        # ---- PRIMARY calibration variant: per-order regression -------------
                        # Applies THIS weighting's order_to_ab (computed ONCE per AdQSM
                        # variant/weighting above, via group_orders_for_fitting() + one
                        # fit_radius_regression() call per group) to THIS threshold's
                        # still-RAW cyl/cyl_order. This is the ADOPTED PRIMARY calibration
                        # method (CHANGELOG_adtree.md, Step 7) - its cylinders
                        # (regperorder_cyl) become `final_cyl`, used for the exported
                        # .npz/geom_*.txt below, replacing the OLD, buggy self-referencing
                        # calibrate_cylinder_radii() (removed).
                        #
                        # Sanity guard (order_to_ab coverage): order_to_ab was built from
                        # order_labels_matched, i.e. only orders that (a) survive in the
                        # UNPRUNED reference set AND (b) have their own entry in AdQSM's
                        # BranchStructure.txt (build_quantile_matched_pairs() skips - and
                        # prints a warning for - any order present in only one of the
                        # two). A given threshold's cyl_order can only ever be a SUBSET
                        # of the unpruned reference's orders (pruning removes cylinders,
                        # it never invents a new order), so this gap can only matter if
                        # AdQSM's own table is missing an order AdTree has - checked
                        # explicitly here (not just left to
                        # apply_radius_regression_per_order()'s internal per-cylinder
                        # warning) so a coverage gap is visible immediately, per
                        # threshold, instead of only inside a buried per-cylinder
                        # warning.
                        cyl_orders_present = set(np.asarray(cyl_order).tolist()) - {0}
                        missing_from_order_to_ab = sorted(cyl_orders_present - set(order_to_ab))
                        if missing_from_order_to_ab:
                            print("  WARNING: order_to_ab has NO fit for order(s) %s present in this "
                                  "threshold's cyl_order (missing from AdQSM's own BranchStructure.txt, "
                                  "or otherwise skipped by build_quantile_matched_pairs) - those cylinders "
                                  "will be left UNSCALED by apply_radius_regression_per_order() below."
                                  % missing_from_order_to_ab)

                        regperorder_new_r = apply_radius_regression_per_order(
                            xyz, cyl, cyl_order, trunk_radius_func, order_to_ab)
                        regperorder_cyl = [(a, b, float(regperorder_new_r[i]), pid)
                                            for i, (a, b, r, pid) in enumerate(cyl)]
                        regperorder_lengths, regperorder_radii = cylinder_metrics(xyz, regperorder_cyl)
                        regperorder_stats = volume_stats(regperorder_lengths, regperorder_radii, np.asarray(cyl_order))
                        if WRITE_THIN_BRANCH_FILTERED_ROW:
                            regperorder_thin = report_thin_branch_volume(
                                regperorder_lengths, regperorder_radii, cyl_order, cut_cm=THIN_BRANCH_CUT_CM,
                                source_label="AdTree calibrated [calmethod=regression-perorder]")
                        else:
                            regperorder_thin = None
                        regperorder_dbh = stem_diameter_at_height(xyz, regperorder_cyl, cyl_order, z_base, TAPER_H_LOWER)
                        regperorder_d_upper = stem_diameter_at_height(xyz, regperorder_cyl, cyl_order, z_base, TAPER_H_UPPER)
                        regperorder_taper = ((regperorder_dbh - regperorder_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                                              if regperorder_dbh is not None and regperorder_d_upper is not None else None)
                        regperorder_n_cylinders = len(regperorder_cyl)

                        # `final_cyl` becomes the PRIMARY-calibrated (per-order regression)
                        # cylinders for THIS weighting - everything below this point (the
                        # .npz save, the "processed, CALIBRATED" report, plot_model()) uses
                        # this. `cyl` itself is left untouched so the NEXT weighting pass
                        # (if any) still sees the same raw geometry.
                        final_cyl = regperorder_cyl

                    # ---- CHANGE vs. the old single-file ply_to_geom.py: save the final
                    # (possibly calibrated) geometry to an .npz file INSTEAD OF calling
                    # write_geom()/writing geom_*.txt directly. Reason: this script does
                    # calibration + comparison + printing, which you may want to re-run
                    # or tweak (e.g. different AdQSM variant) without re-exporting to
                    # ANSYS every time, and conversely you may want to re-export to
                    # ANSYS without redoing the whole calibration. Splitting the
                    # pipeline here lets export_geom_ansys.py
                    # do ONLY the second half, fast, from already-calibrated data.
                    #
                    # What goes into the .npz (so it can be reloaded with NO information
                    # loss - i.e. write_geom() on the reloaded data produces a BIT-IDENTICAL
                    # geom_*.txt to what the old single-file script would have written):
                    #   xyz          : (N,3) float64 - ALL node coordinates (root, and every
                    #                  node any cylinder in `final_cyl` references by index).
                    #                  This is the SAME xyz array used above throughout
                    #                  calibration - not cropped/renumbered, so final_cyl's
                    #                  (a, b) indices stay valid after reloading.
                    #   cyl           : (n_cyl,4) float64 - one row per cylinder, columns
                    #                  [a, b, radius, parent_cyl_id] (a/b/parent are node/
                    #                  cylinder INDICES, stored as float64 for a uniform
                    #                  array; export_geom_ansys.py casts them back to int).
                    #   cyl_order     : (n_cyl,) int - branch order per cylinder (0=trunk,
                    #                  >=1=branch) - write_geom() now writes this as the
                    #                  11th geom_*.txt column (see tree_geom_utils.py),
                    #                  and it's also kept here so nothing is lost if you
                    #                  want to recompute volume_stats()/report_volume()
                    #                  etc. from the .npz later.
                    #   root          : the root node index (scalar) - write_geom() needs
                    #                  it to compute the x,y recentring offset.
                    #   recenter_xy   : the RECENTER_XY flag used for THIS run (scalar bool).
                    #   geom_filename : the `out` filename computed above - so
                    #                  export_geom_ansys.py writes the SAME geom_*.txt name
                    #                  this script would have used, without recomputing
                    #                  the tree/threshold/variant-suffix naming logic.
                    #   tree_name, variant_label, threshold_m : just metadata, so you can
                    #                  tell which run produced a given .npz file later.
                    # SEG_VARIANT_SUFFIX added here too (after variant_suffix, same "end
                    # of the name" placement as the method names above) - a run with
                    # different SEG_LEN_MIN/MAX/K/weighting settings now writes a
                    # DIFFERENT .npz file on disk instead of silently overwriting the
                    # previous run's one.
                    npz_name = os.path.join(
                        NPZ_DIR, "calib_%s_r%dmm%s%s.npz" % (TREE_NAME, round(thr * 1000), variant_suffix, SEG_VARIANT_SUFFIX))
                    cyl_array = np.array([(a, b, r, pid) for a, b, r, pid in final_cyl], dtype=np.float64)
                    np.savez(npz_name,
                             xyz=xyz,
                             cyl=cyl_array,
                             cyl_order=np.asarray(cyl_order, dtype=np.int64),
                             root=np.array(root),
                             recenter_xy=np.array(RECENTER_XY),
                             geom_filename=np.array(out),
                             tree_name=np.array(TREE_NAME),
                             variant_label=np.array(variant_label if variant_label else ""),
                             threshold_m=np.array(thr))

                    total_len = sum(float(np.linalg.norm(xyz[b] - xyz[a])) for a, b, _, _ in final_cyl)
                    print("%-12s %-12d %-12.1f %-12s" % ("%d mm" % round(thr * 1000), len(final_cyl), total_len, npz_name))
                    report_volume(xyz, final_cyl, thr)   # uses the (possibly calibrated) radii above

                    if CALIBRATE_RADII:
                        # cal_stats/cal_dbh/cal_taper/cal_thin (used in the report block
                        # below) are exactly regperorder_stats/regperorder_dbh/
                        # regperorder_taper/regperorder_thin computed above for THIS
                        # weighting - aliased under their original names purely so the
                        # report block below (predating the multi-method investigation)
                        # doesn't need renaming.
                        cal_stats, cal_dbh, cal_taper, cal_thin = (
                            regperorder_stats, regperorder_dbh, regperorder_taper, regperorder_thin)

                        # ---- upsert both the uncalibrated and calibrated rows for this threshold ----
                        # "AdTree raw" does NOT depend on AdQSM at all, so it gets no variant
                        # suffix - it's simply re-written (with identical VALUES) for every
                        # variant AND for every weighting in this loop (raw AdTree radii
                        # don't depend on weighting either), which is harmless since
                        # upsert_result overwrites by (tree, method), not duplicates - the
                        # per-weighting SEG_VARIANT_SUFFIX just gives each weighting its own
                        # otherwise-identical row instead of colliding.
                        # branch_filter = "none": raw AdTree radii, no diameter cut-off applied.
                        upsert_result(RESULTS_CSV, TREE_NAME,
                                      "AdTree raw r%dmm%s" % (round(thr * 1000), SEG_VARIANT_SUFFIX),
                                      orig_stats["total_vol"], orig_stats["trunk_vol"], orig_stats["branch_vol"], None,
                                      raw_dbh, height_m, raw_taper,
                                      # trunk_len/branch_len: already in this dict (volume_stats()
                                      # computes them the same way as trunk_vol/branch_vol).
                                      orig_stats["trunk_len"], orig_stats["branch_len"],
                                      branch_filter="none",
                                      # n_cylinders (Task B): total cylinder count for this
                                      # threshold's reconstruction. Raw and calibrated share
                                      # the exact same count - calibration only replaces
                                      # radii, it never adds, removes, or splits cylinders,
                                      # so len(cyl) here is identical to len(final_cyl) at
                                      # every calibrated row below.
                                      n_cylinders=len(cyl),
                                      # adqsm_variant=None (not variant_label): raw AdTree
                                      # geometry never touches AdQSM at all, so it doesn't
                                      # actually depend on which variant happened to be
                                      # active during this loop iteration - the method
                                      # name itself already confirms this (no variant
                                      # suffix). Passing variant_label here would make
                                      # assign_adtree_groups() (plot_box.py) accidentally
                                      # lump raw rows into a calibrated row's group
                                      # whenever their radius_threshold_mm/seg_* happen to
                                      # match (see STEP 5's fix). radius_threshold_mm/
                                      # seg_min_mm/seg_max_mm/seg_k_pct are kept - raw
                                      # AdTree DOES genuinely depend on the pruning
                                      # threshold and resampling settings.
                                      adqsm_variant=None, radius_threshold_mm=round(thr * 1000),
                                      seg_min_mm=round(SEG_LEN_MIN * 1000), seg_max_mm=round(SEG_LEN_MAX * 1000),
                                      seg_k_pct=round(SEG_LEN_K * 100))

                        if WRITE_THIN_BRANCH_FILTERED_ROW:
                            # Same idea as the calibrated (>=10cm only) rows further below,
                            # but for the RAW (uncalibrated) cylinders instead - uses
                            # orig_thin (computed once per threshold, above, BEFORE this
                            # weighting loop). This row does NOT depend on which AdQSM
                            # variant or weighting is active (raw AdTree radii never touch
                            # AdQSM at all - same reasoning as the plain "AdTree raw" row
                            # above), so it gets no variant suffix either.
                            # DBH/taper reuse raw_dbh/raw_taper (the UNCALIBRATED trunk's
                            # own values), not cal_dbh/cal_taper, to stay consistent with
                            # "this row describes the raw model, not the calibrated one."
                            # SEG_VARIANT_SUFFIX at the very end again, same rule as above.
                            upsert_result(RESULTS_CSV, TREE_NAME,
                                          "AdTree raw r%dmm (>=%.0fcm only)%s"
                                          % (round(thr * 1000), THIN_BRANCH_CUT_CM, SEG_VARIANT_SUFFIX),
                                          orig_thin["total_vol_kept"], orig_thin["trunk_vol_kept"],
                                          orig_thin["branch_vol_kept"], None,
                                          raw_dbh, height_m, raw_taper,
                                          # Same fix as the calibrated row above, using the "raw"
                                          # (uncalibrated) cylinder set's kept lengths instead.
                                          orig_thin["trunk_len_kept"], orig_thin["branch_len_kept"],
                                          branch_filter="10cm",
                                          # n_cylinders (Task B): same idea as the calibrated
                                          # row above, using orig_thin's "n_cyl_kept" (the raw/
                                          # uncalibrated cylinder set's filtered count) instead.
                                          n_cylinders=orig_thin["n_cyl_kept"],
                                          # adqsm_variant=None - same reasoning as the
                                          # plain "AdTree raw" row above (STEP 5 fix).
                                          adqsm_variant=None, radius_threshold_mm=round(thr * 1000),
                                          seg_min_mm=round(SEG_LEN_MIN * 1000), seg_max_mm=round(SEG_LEN_MAX * 1000),
                                          seg_k_pct=round(SEG_LEN_K * 100))

                        # ---- SECONDARY calibration variant: calref=min5mm - upsert 2 rows --
                        # (one "none"/full, one "(>=10cm only)") for the retained secondary
                        # reference (fixedref_data: just "min5mm" now - see
                        # CALIBRATION_REF_THRESHOLDS_MM's definition above and
                        # CHANGELOG_adtree.md, Step 7). Method-name tag "[calref=minXmm]"
                        # sits right after variant_method_suffix, same position the
                        # "(AdQSM 08)" variant tag already occupies. Iterates
                        # fixedref_data's own keys rather than a hard-coded tuple, so this
                        # stays correct even if CALIBRATION_REF_THRESHOLDS_MM ever grows
                        # again, with no second place to keep in sync.
                        for ref_name in fixedref_data:
                            fr = fixedref_data[ref_name]
                            upsert_result(RESULTS_CSV, TREE_NAME,
                                          "AdTree calibrated r%dmm%s [calref=%s]%s"
                                          % (round(thr * 1000), variant_method_suffix, ref_name, SEG_VARIANT_SUFFIX),
                                          fr["stats"]["total_vol"], fr["stats"]["trunk_vol"], fr["stats"]["branch_vol"], None,
                                          fr["dbh"], height_m, fr["taper"],
                                          fr["stats"]["trunk_len"], fr["stats"]["branch_len"],
                                          branch_filter="none",
                                          n_cylinders=fr["n_cylinders"],
                                          adqsm_variant=variant_label, radius_threshold_mm=round(thr * 1000),
                                          seg_min_mm=round(SEG_LEN_MIN * 1000), seg_max_mm=round(SEG_LEN_MAX * 1000),
                                          seg_k_pct=round(SEG_LEN_K * 100), calmethod=ref_name)
                            if WRITE_THIN_BRANCH_FILTERED_ROW:
                                fr_thin = fr["thin"]
                                upsert_result(RESULTS_CSV, TREE_NAME,
                                              "AdTree calibrated r%dmm%s [calref=%s] (>=%.0fcm only)%s"
                                              % (round(thr * 1000), variant_method_suffix, ref_name,
                                                 THIN_BRANCH_CUT_CM, SEG_VARIANT_SUFFIX),
                                              fr_thin["total_vol_kept"], fr_thin["trunk_vol_kept"],
                                              fr_thin["branch_vol_kept"], None,
                                              fr["dbh"], height_m, fr["taper"],
                                              fr_thin["trunk_len_kept"], fr_thin["branch_len_kept"],
                                              branch_filter="10cm",
                                              n_cylinders=fr_thin["n_cyl_kept"],
                                              adqsm_variant=variant_label, radius_threshold_mm=round(thr * 1000),
                                              seg_min_mm=round(SEG_LEN_MIN * 1000), seg_max_mm=round(SEG_LEN_MAX * 1000),
                                              seg_k_pct=round(SEG_LEN_K * 100), calmethod=ref_name)

                        # ---- PRIMARY calibration variant: per-order regression - upsert 2 rows --
                        # Mirrors the "none"/"(>=10cm only)" pattern above exactly, using
                        # the regperorder_* values computed earlier (from
                        # apply_radius_regression_per_order() with this weighting's
                        # order_to_ab) - the ADOPTED PRIMARY calibration method
                        # (CHANGELOG_adtree.md, Step 7).
                        upsert_result(RESULTS_CSV, TREE_NAME,
                                      "AdTree calibrated r%dmm%s [calmethod=regression-perorder]%s"
                                      % (round(thr * 1000), variant_method_suffix, SEG_VARIANT_SUFFIX),
                                      regperorder_stats["total_vol"], regperorder_stats["trunk_vol"],
                                      regperorder_stats["branch_vol"], None,
                                      regperorder_dbh, height_m, regperorder_taper,
                                      regperorder_stats["trunk_len"], regperorder_stats["branch_len"],
                                      branch_filter="none",
                                      n_cylinders=regperorder_n_cylinders,
                                      adqsm_variant=variant_label, radius_threshold_mm=round(thr * 1000),
                                      seg_min_mm=round(SEG_LEN_MIN * 1000), seg_max_mm=round(SEG_LEN_MAX * 1000),
                                      seg_k_pct=round(SEG_LEN_K * 100), calmethod="regression-perorder")
                        if WRITE_THIN_BRANCH_FILTERED_ROW:
                            upsert_result(RESULTS_CSV, TREE_NAME,
                                          "AdTree calibrated r%dmm%s [calmethod=regression-perorder] (>=%.0fcm only)%s"
                                          % (round(thr * 1000), variant_method_suffix,
                                             THIN_BRANCH_CUT_CM, SEG_VARIANT_SUFFIX),
                                          regperorder_thin["total_vol_kept"], regperorder_thin["trunk_vol_kept"],
                                          regperorder_thin["branch_vol_kept"], None,
                                          regperorder_dbh, height_m, regperorder_taper,
                                          regperorder_thin["trunk_len_kept"], regperorder_thin["branch_len_kept"],
                                          branch_filter="10cm",
                                          n_cylinders=regperorder_thin["n_cyl_kept"],
                                          adqsm_variant=variant_label, radius_threshold_mm=round(thr * 1000),
                                          seg_min_mm=round(SEG_LEN_MIN * 1000), seg_max_mm=round(SEG_LEN_MAX * 1000),
                                          seg_k_pct=round(SEG_LEN_K * 100), calmethod="regression-perorder")

                        print("  DBH (at %.1f m)   : raw AdTree = %s   |   calibrated = %s"
                              % (TAPER_H_LOWER, _fmt_dbh(raw_dbh), _fmt_dbh(cal_dbh)))
                        print("  Taper (%.1f-%.1f m): raw AdTree = %s   |   calibrated = %s"
                              % (TAPER_H_LOWER, TAPER_H_UPPER, _fmt_taper(raw_taper), _fmt_taper(cal_taper)))
                        print("  Height (pruned model): %s" % (("%.2f m" % height_m) if height_m is not None else "n/a"))

                        print("  Volume comparison (a) raw skeleton vs. (b) processed/AdTree vs. (c) processed/calibrated:")
                        print_volume_stats("(a) raw skeleton (AdTree)", raw_stats)
                        print_volume_stats("(b) processed, AdTree radii", orig_stats)
                        print_volume_stats("(c) processed, CALIBRATED [calmethod=regression-perorder]", cal_stats)
                    print()

                    if SHOW_PLOT or SAVE_PLOT_PNG:
                        # `out` itself is NOT touched here (see the NPZ_DIR/FIGURES_DIR
                        # comment above) - it's also stored verbatim as `geom_filename`
                        # inside the .npz below, for export_geom_ansys.py to read back
                        # later as the bare (no-folder) name it should write.
                        # ensure_tree_geom_dir(TREE_NAME) (plots/<tree>/geom/) is
                        # prefixed ONLY at this call site, purely to steer where
                        # plot_model() derives its PNG path from (out.txt -> out.png) -
                        # this is the tree-SHAPE render PNG, routed to plots/<tree>/
                        # geom/ specifically (NOT plots/<tree>/ directly, which is
                        # where the adtree_adqsm_radius_regression_perorder_*.png
                        # diagnostic above goes via FIGURES_DIR - a different artifact).
                        plot_model(xyz, final_cyl, root, RECENTER_XY, thr,
                                   os.path.join(ensure_tree_geom_dir(TREE_NAME), out), SHOW_PLOT, SAVE_PLOT_PNG)

                if PRINT_TIMING:
                    print("  [TIMING] RADIUS_THRESHOLDS iteration thr=%.3f (total): %.2f s"
                          % (thr, time.perf_counter() - _thr_start))

    if PRINT_TIMING:
        print("[TIMING] AdQSM variant '%s' (total): %.2f s"
              % (variant_label, time.perf_counter() - _variant_start))

if PRINT_TIMING:
    print("[TIMING] WHOLE RUN (total): %.2f s" % (time.perf_counter() - _run_start))

print("[seg sweep] %d combination(s) run - total wall-clock time: %.2f s"
      % (_n_seg_combinations, time.perf_counter() - _sweep_run_start))
