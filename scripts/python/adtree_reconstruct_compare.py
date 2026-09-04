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

# =====================  PARAMETERS  ===================================
# Tree ID for THIS run. This is the ONLY thing you need to change to switch
# trees - it names this tree's row in the shared results table (RESULTS_CSV,
# see upsert_result() calls below) AND builds AdQSM_DIR/AdTree_DIR/INPUT_PLY
# right below it automatically, so those don't need editing separately.
TREE_NAME = "IND01_054"

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
ADQSM_VARIANTS = ["04", "05", "06","07", "08", "09","10"]   # e.g. ["05", "08"]

AdTree_DIR = os.path.join(DATA_ROOT, TREE_NAME)

# input skeleton from AdTree - filename follows the "<TREE_NAME> - Cloud_skeleton.ply"
# convention used for every tree, so it's derived from TREE_NAME too.
INPUT_PLY = os.path.join(AdTree_DIR, "%s - Cloud_skeleton.ply" % TREE_NAME)

# Radius threshold (in METERS). You can give several values -> several variants.
# Remove all branches whose radius is below this threshold. The trunk (branch order 0) is never removed, even if its radius is below the threshold.
# Example of a single variant:   RADIUS_THRESHOLDS = [0.010]
# Example of several variants:   RADIUS_THRESHOLDS = [0.010, 0.020, 0.030]
RADIUS_THRESHOLDS = [0.005,
                     0.010,
                     0.015,
                     0.020]        # 0.030 m = 30 mm radius (60 mm diameter)

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
COMPUTE_CALREF_MIN5MM = False   # secondary/backup calibration check
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
SEG_LEN_MIN = 0.01                # shortest allowed segment (m), for thin twigs
SEG_LEN_MAX = 0.5                # longest allowed segment (m), for the trunk
SEG_LEN_K = 0.5                   # target length = SEG_LEN_K * local_radius

# --- Build a short suffix identifying THIS resampling configuration -----
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
# cells. Example: SEG_LEN_MIN=0.01, SEG_LEN_MAX=0.3, SEG_LEN_K=0.5 (the
# current values above) -> "_seg10-300-k50".
SEG_VARIANT_SUFFIX = "_seg%d-%d-k%d" % (
    round(SEG_LEN_MIN * 1000),   # metres -> whole millimetres
    round(SEG_LEN_MAX * 1000),   # metres -> whole millimetres
    round(SEG_LEN_K * 100),      # e.g. 0.5 -> 50
)

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
WRITE_THIN_BRANCH_FILTERED_ROW = True


# Shared master results table (see compare_volumes.py). When CALIBRATE_RADII
# is True, each generated threshold variant upserts its own row into this CSV
# so results from all methods live in one place.
RESULTS_CSV = "volume_results.csv"

# Output folders, so the working directory doesn't fill up with dozens of
# .npz/.png files mixed in with the scripts. NPZ_DIR is this script's own
# new folder; FIGURES_DIR reuses the project's existing "plots/" convention
# (already used by plot_box.py/plot_volumes.py for their own charts),
# grouped under a per-tree subfolder. export_geom_ansys.py has its OWN
# matching NPZ_DIR parameter (see that file) - keep both in sync by hand if
# this one ever changes.
NPZ_DIR = "npz"
FIGURES_DIR = os.path.join("plots", TREE_NAME)

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
os.makedirs(NPZ_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

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

print("\n%-12s %-12s %-12s %-12s" % ("threshold", "cylinders", "length [m]", "file"))
print("Volume verification below is computed from the exact cylinders that will be "
      "written to each geom file (pi * radius^2 * length per cylinder).\n")
z_base = float(xyz[:, 2].min())   # tree base; DBH/height/taper are measured from here
multiple_variants = len(ADQSM_VARIANT_LIST) > 1   # True only if you used ADQSM_VARIANTS (case 2 above)

# Outer loop: one pass per AdQSM variant (just one pass, using the plain
# AdQSM_DIR, unless you filled in ADQSM_VARIANTS). Everything inside this
# loop (loading AdQSM data, calibrating, saving results) is repeated once
# per variant, so several reconstructions can be compared side by side in
# the same RESULTS_CSV without overwriting each other.
for variant_label, taper_file, branch_file, params_file in ADQSM_VARIANT_LIST:
    # variant_suffix/variant_method_suffix are "" when there is only one
    # variant (so filenames/method names look exactly like before this
    # feature existed); otherwise they tag the variant name onto them.
    variant_suffix = ("_adqsm%s" % variant_label) if variant_label else ""
    variant_method_suffix = (" (AdQSM %s)" % variant_label) if variant_label else ""

    if CALIBRATE_RADII:
        print("\nLoading AdQSM calibration data%s..."
              % (" (variant: %s)" % variant_label if variant_label else ""))
        taper_heights, taper_diameters = parse_adqsm_taper_file(taper_file)
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
        report_adqsm_thin_branch(branch_file, cut_cm=THIN_BRANCH_CUT_CM, params_file=params_file)

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
        ref_root_0, ref_cyl_0, ref_order_0 = convert(
            xyz, rad, edges, 0.0, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)

        # factors_by_ref: {ref_thr: factors_dict} - one fixed factors dict per
        # CALIBRATION_REF_THRESHOLDS_MM entry (just min5mm now), computed
        # ONCE here (not inside the RADIUS_THRESHOLDS loop below).
        factors_by_ref = {}
        if COMPUTE_CALREF_MIN5MM:
            for ref_thr in CALIBRATION_REF_THRESHOLDS_MM:
                ref_root, ref_cyl, ref_order = convert(
                    xyz, rad, edges, ref_thr, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
                factors_by_ref[ref_thr] = compute_order_calibration_factors(
                    ref_cyl, ref_order, adqsm_median_by_order)
                print("  Fixed calibration factors - reference set 'min%dmm' (thr=%.3f), "
                      "reused for every RADIUS_THRESHOLDS run:" % (round(ref_thr * 1000), ref_thr))
                for o in sorted(factors_by_ref[ref_thr]):
                    print("    order %d : factor = %.3f" % (o, factors_by_ref[ref_thr][o]))

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
        # per AdQSM variant here, reused for every RADIUS_THRESHOLDS value
        # below.
        adtree_matched, adqsm_matched, order_labels_matched = build_quantile_matched_pairs(
            ref_cyl_0, ref_order_0, raw_diam_by_order)

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
            TREE_NAME, variant_label, order1_merge_note=order1_merge_note)
        print("  Saved per-order regression diagnostic plot: %s" % regression_perorder_plot_path)

    # Inner loop: one pass per radius threshold (same as before this feature
    # existed), now repeated for each AdQSM variant above.
    for thr in RADIUS_THRESHOLDS:
        root, cyl, cyl_order = convert(xyz, rad, edges, thr, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
        # `out` is the geom_*.txt name step 2 (export_geom_ansys.py) will
        # eventually write - computed here (once, alongside the threshold/
        # variant it belongs to) and carried inside the .npz below.
        #
        # Built with the EXACT SAME ingredients (and in the same order) as
        # npz_name further down - just "geom_"/".txt" instead of
        # "calib_"/".npz" - so every geom_*.txt name matches the calib_*.npz
        # it was exported from at a glance, e.g.:
        #   calib_IND01_054_r5mm_seg100-500-k50.npz
        #   geom_IND01_054_r5mm_seg100-500-k50.txt
        # This replaces the old fixed OUTPUT_NAME/OUTPUT_PATTERN constants -
        # every combination of tree/threshold/variant/segment-settings now
        # gets its own name automatically, so nothing can silently overwrite
        # a previous run's exported file.
        out = "geom_%s_r%dmm%s%s.txt" % (TREE_NAME, round(thr * 1000), variant_suffix, SEG_VARIANT_SUFFIX)

        # Height of the pruned model: z-range of the nodes actually used by these
        # cylinders. Unaffected by radius calibration (geometry doesn't change).
        node_ids = sorted({idx for a, b, r, pid in cyl for idx in (a, b)})
        height_m = float(xyz[node_ids, 2].max() - xyz[node_ids, 2].min()) if node_ids else None

        if CALIBRATE_RADII:
            orig_lengths, orig_radii = cylinder_metrics(xyz, cyl)
            orig_stats = volume_stats(orig_lengths, orig_radii, np.asarray(cyl_order))
            # DBH/taper of the UNCALIBRATED (raw AdTree) trunk, before radii are replaced.
            raw_dbh = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_LOWER)
            raw_d_upper = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_UPPER)
            raw_taper = ((raw_dbh - raw_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                         if raw_dbh is not None and raw_d_upper is not None else None)

            # ---- thin-branch diagnostic on the RAW (uncalibrated) cylinders ----
            # Computed here (BEFORE `cyl`'s radii get overwritten by the
            # per-order regression calibration below) since orig_lengths/
            # orig_radii are only valid for the CURRENT (raw AdTree) radii
            # at this point.
            # source_label="AdTree raw" makes this printout visually distinct
            # from the calibrated one further below (same function, same cut_cm,
            # different cylinder set) - see report_thin_branch_volume()'s
            # docstring in tree_geom_utils.py for why the label exists.
            orig_thin = report_thin_branch_volume(orig_lengths, orig_radii, cyl_order,
                                                   cut_cm=THIN_BRANCH_CUT_CM, source_label="AdTree raw")

            # ---- SECONDARY calibration variant: calref=min5mm -----------------
            # Apply the FIXED factors_by_ref[...] dict (computed once per
            # AdQSM variant, above the RADIUS_THRESHOLDS loop, from the fixed
            # min5mm reference cylinder set) to THIS threshold's still-RAW
            # `cyl`/`cyl_order` - captured here BEFORE the primary
            # (per-order regression) calibration below overwrites `cyl`'s
            # radii. Kept as a secondary/backup reference point alongside the
            # primary per-order regression rows below (see
            # CHANGELOG_adtree.md, Step 7).
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
                    fr_thin = report_thin_branch_volume(
                        fr_lengths, fr_radii, cyl_order, cut_cm=THIN_BRANCH_CUT_CM,
                        source_label="AdTree calibrated [calref=%s]" % ref_name)
                    fr_dbh = stem_diameter_at_height(xyz, fr_cyl, cyl_order, z_base, TAPER_H_LOWER)
                    fr_d_upper = stem_diameter_at_height(xyz, fr_cyl, cyl_order, z_base, TAPER_H_UPPER)
                    fr_taper = ((fr_dbh - fr_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                                if fr_dbh is not None and fr_d_upper is not None else None)
                    fixedref_data[ref_name] = dict(stats=fr_stats, thin=fr_thin, dbh=fr_dbh,
                                                    taper=fr_taper, n_cylinders=len(fr_cyl))

            # ---- PRIMARY calibration variant: per-order regression -------------
            # Applies order_to_ab (computed ONCE per AdQSM variant above, via
            # group_orders_for_fitting() + one fit_radius_regression() call
            # per group) to THIS threshold's still-RAW cyl/cyl_order. This is
            # the ADOPTED PRIMARY calibration method (CHANGELOG_adtree.md,
            # Step 7) - its cylinders (regperorder_cyl) become the final `cyl`
            # used for the exported .npz/geom_*.txt below, replacing the OLD,
            # buggy self-referencing calibrate_cylinder_radii() (removed).
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
            regperorder_thin = report_thin_branch_volume(
                regperorder_lengths, regperorder_radii, cyl_order, cut_cm=THIN_BRANCH_CUT_CM,
                source_label="AdTree calibrated [calmethod=regression-perorder]")
            regperorder_dbh = stem_diameter_at_height(xyz, regperorder_cyl, cyl_order, z_base, TAPER_H_LOWER)
            regperorder_d_upper = stem_diameter_at_height(xyz, regperorder_cyl, cyl_order, z_base, TAPER_H_UPPER)
            regperorder_taper = ((regperorder_dbh - regperorder_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                                  if regperorder_dbh is not None and regperorder_d_upper is not None else None)
            regperorder_n_cylinders = len(regperorder_cyl)

            # `cyl` now becomes the PRIMARY-calibrated (per-order regression)
            # cylinders - everything below this point (the .npz save, the
            # "processed, CALIBRATED" report, plot_model()) uses this, same
            # as the OLD self-referencing calibrate_cylinder_radii() call
            # used to reassign `cyl` here (see CHANGELOG_adtree.md, Step 7,
            # for why that call was removed).
            cyl = regperorder_cyl

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
        #                  node any cylinder in `cyl` references by index).
        #                  This is the SAME xyz array used above throughout
        #                  calibration - not cropped/renumbered, so cyl's
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
        # different SEG_LEN_MIN/MAX/K settings now writes a DIFFERENT .npz
        # file on disk instead of silently overwriting the previous run's one.
        npz_name = os.path.join(
            NPZ_DIR, "calib_%s_r%dmm%s%s.npz" % (TREE_NAME, round(thr * 1000), variant_suffix, SEG_VARIANT_SUFFIX))
        cyl_array = np.array([(a, b, r, pid) for a, b, r, pid in cyl], dtype=np.float64)
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

        total_len = sum(float(np.linalg.norm(xyz[b] - xyz[a])) for a, b, _, _ in cyl)
        print("%-12s %-12d %-12.1f %-12s" % ("%d mm" % round(thr * 1000), len(cyl), total_len, npz_name))
        report_volume(xyz, cyl, thr)   # uses the (possibly calibrated) radii above

        if CALIBRATE_RADII:
            # cal_stats/cal_dbh/cal_taper/cal_thin (used in the report block
            # below) are exactly regperorder_stats/regperorder_dbh/
            # regperorder_taper/regperorder_thin computed above - `cyl` was
            # already reassigned to regperorder_cyl right after that block
            # (the PRIMARY calibration method, see CHANGELOG_adtree.md, Step
            # 7), so recomputing them here would just repeat the same
            # numbers. Aliased under their original names purely so the
            # report block below (predating the multi-method investigation)
            # doesn't need renaming.
            cal_stats, cal_dbh, cal_taper, cal_thin = (
                regperorder_stats, regperorder_dbh, regperorder_taper, regperorder_thin)

            # ---- upsert both the uncalibrated and calibrated rows for this threshold ----
            # "AdTree raw" does NOT depend on AdQSM at all, so it gets no variant
            # suffix - it's simply re-written (with identical values) for every
            # variant, which is harmless since upsert_result overwrites by
            # (tree, method), not duplicates.
            # branch_filter = "none": raw AdTree radii, no diameter cut-off applied.
            # SEG_VARIANT_SUFFIX appended at the very end of the method name
            # (see where it's built, next to SEG_LEN_MIN/MAX/K above) - keeps
            # results from a different resampling setting as a SEPARATE row
            # instead of overwriting this one.
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
                          # so len(cyl) here is identical to len(cyl) at every
                          # calibrated row below.
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
                # orig_thin (computed earlier from orig_lengths/orig_radii,
                # BEFORE any calibration replaced `cyl`'s radii). This row does
                # NOT depend on which AdQSM variant is active (raw AdTree radii
                # never touch AdQSM at all - same reasoning as the plain
                # "AdTree raw" row above), so it gets no variant suffix either.
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
                              seg_k_pct=round(SEG_LEN_K * 100))
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
                                  seg_k_pct=round(SEG_LEN_K * 100))

            # ---- PRIMARY calibration variant: per-order regression - upsert 2 rows --
            # Mirrors the "none"/"(>=10cm only)" pattern above exactly, using
            # the regperorder_* values computed earlier (from
            # apply_radius_regression_per_order() with this variant's
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
                          seg_k_pct=round(SEG_LEN_K * 100))
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
                              seg_k_pct=round(SEG_LEN_K * 100))

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
            # later as the bare (no-folder) name it should write. FIGURES_DIR
            # is prefixed ONLY at this call site, purely to steer where
            # plot_model() derives its PNG path from (out.txt -> out.png).
            plot_model(xyz, cyl, root, RECENTER_XY, thr, os.path.join(FIGURES_DIR, out), SHOW_PLOT, SAVE_PLOT_PNG)
