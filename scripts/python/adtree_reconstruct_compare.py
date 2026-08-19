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
#   8) Optionally CALIBRATES the cylinder radii against AdQSM data (a taper
#      curve for the trunk, and median branch diameters per branch order for
#      the rest), replacing the AdTree radii before writing.
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
    calibrate_cylinder_radii, convert, cylinder_metrics, volume_stats,
    report_thin_branch_volume, upsert_result, stem_diameter_at_height,
    _fmt_dbh, _fmt_taper, print_volume_stats, raw_skeleton_stats,
    report_volume, plot_model,
)

# =====================  PARAMETERS  ===================================
# Directory holding this tree's source input files (the AdTree skeleton .ply,
# and the AdQSM taper/branch/params exports). Edit this to point at a
# different tree's folder; INPUT_PLY and the ADQSM_*_FILE paths below are all
# resolved relative to it.
#
# --- AdQSM variant(s) ---------------------------------------------------
# AdQSM can be reconstructed several times with different settings, each
# saved in its own subfolder (e.g. ".../AdQSM/05", ".../AdQSM/08" - the
# folder name is just whatever you called that reconstruction run). You can
# either:
#   (1) point at ONE such folder with AdQSM_DIR (simple, old behaviour), or
#   (2) list SEVERAL subfolder names in ADQSM_VARIANTS to process all of
#       them in a single run of this script (similar to how RADIUS_THRESHOLDS
#       lets you try several radius thresholds in one run).
#
# Case (1) - single variant (default, still works exactly as before):
AdQSM_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data\IND01_54\AdQSM\05"

# Case (2) - several variants. Leave ADQSM_VARIANTS empty/None to use only
# AdQSM_DIR above (case 1). To use several variants instead, set BOTH:
#   ADQSM_BASE_DIR = r"C:\...\data\IND01_54\AdQSM"
#   ADQSM_VARIANTS = ["05", "08"]
# Each name in ADQSM_VARIANTS must be a subfolder of ADQSM_BASE_DIR that
# contains its own taper.txt, BranchStructure.txt and TreesParams.txt.
ADQSM_BASE_DIR = None
ADQSM_VARIANTS = None   # e.g. ["05", "08"]

AdTree_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data\IND01_54"

INPUT_PLY = os.path.join(AdTree_DIR, "IND01_054 - Cloud_skeleton.ply")   # input skeleton from AdTree

# Radius threshold (in METERS). You can give several values -> several variants.
# Example of a single variant:   RADIUS_THRESHOLDS = [0.010]
# Example of several variants:   RADIUS_THRESHOLDS = [0.010, 0.020, 0.030]
RADIUS_THRESHOLDS = [0.005]        # 0.030 m = 30 mm radius (60 mm diameter)

# Adaptive segment length used for resampling (in METERS). The target length at
# a given point is SEG_LEN_K * local_radius, clamped to [SEG_LEN_MIN, SEG_LEN_MAX]:
# thick branches (trunk) get long segments, thin twigs get short/fine ones.
# If SEG_LEN_MIN == SEG_LEN_MAX, this reduces to the old constant-length
# resampling (that fixed value, regardless of radius).
# Set SEG_LEN_MIN to 0 or None to disable resampling entirely (keep every point).
SEG_LEN_MIN = 0.01                # shortest allowed segment (m), for thin twigs
SEG_LEN_MAX = 0.3                # longest allowed segment (m), for the trunk
SEG_LEN_K = 0.5                   # target length = SEG_LEN_K * local_radius

# Shortest permissible cylinder (in METERS). Shorter ones (e.g. at branch
# points) are not created, to avoid degenerate zero-length beams in ANSYS.
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

# Tree ID used for this tree's rows in the shared results table (RESULTS_CSV).
TREE_NAME = "IND01_054"

# Shared master results table (see compare_volumes.py). When CALIBRATE_RADII
# is True, each generated threshold variant upserts its own row into this CSV
# so results from all methods live in one place.
RESULTS_CSV = "volume_results.csv"

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
# These name the FINAL geom_*.txt that export_geom_ansys.py will eventually
# write (this script itself only writes the intermediate .npz - see the
# RUN section). Keeping the naming here (instead of in export_geom_ansys.py)
# means the geom_*.txt name is decided once, alongside the threshold/variant
# it belongs to, and carried inside the .npz - so step 2 never has to guess it.
#
# OUTPUT_NAME is used as-is when RADIUS_THRESHOLDS has exactly ONE value.
OUTPUT_NAME = "geom_r30_optim.txt"
# OUTPUT_PATTERN is used instead when RADIUS_THRESHOLDS has SEVERAL values,
# so the generated files don't overwrite each other. {r} is replaced by the
# radius threshold in mm (e.g. geom_r10.txt).
OUTPUT_PATTERN = "geom_r{r}.txt"
# =====================================================================

# =====================  VISUALIZATION  ================================
# Show an interactive 3D preview of the reduced beam model (one window per
# radius threshold) so you can check it before importing into ANSYS.
SHOW_PLOT = True

# Also save the preview as a PNG next to each output file (True/False).
# The PNG name is derived from the output file name (.txt -> .png).
SAVE_PLOT_PNG = False
# =====================================================================


# =========================  RUN  =====================================
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
multiple_thresholds = len(RADIUS_THRESHOLDS) > 1
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
        report_adqsm_thin_branch(branch_file, cut_cm=THIN_BRANCH_CUT_CM, params_file=params_file)

    # Inner loop: one pass per radius threshold (same as before this feature
    # existed), now repeated for each AdQSM variant above.
    for thr in RADIUS_THRESHOLDS:
        root, cyl, cyl_order = convert(xyz, rad, edges, thr, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
        # `out` is the geom_*.txt name step 2 (export_geom_ansys.py) will
        # eventually write - computed here (once, alongside the threshold/
        # variant it belongs to) and carried inside the .npz below.
        out = OUTPUT_PATTERN.format(r=int(round(thr * 1000))) if multiple_thresholds else OUTPUT_NAME
        if variant_suffix:
            # insert the variant tag just before the file extension, e.g.
            # "geom_r30_optim.txt" -> "geom_r30_optim_adqsm05.txt"
            name_part, ext_part = os.path.splitext(out)
            out = name_part + variant_suffix + ext_part

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

            new_radii, factors = calibrate_cylinder_radii(
                xyz, cyl, cyl_order, trunk_radius_func, adqsm_median_by_order)
            cyl = [(a, b, float(new_radii[i]), pid) for i, (a, b, r, pid) in enumerate(cyl)]

        # ---- CHANGE vs. the old single-file ply_to_geom.py: save the final
        # (possibly calibrated) geometry to an .npz file INSTEAD OF calling
        # write_geom()/writing geom_*.txt directly. Reason: this script does
        # calibration + comparison + printing, which you may want to re-run
        # or tweak (e.g. different AdQSM variant) without re-exporting to
        # ANSYS every time, and conversely you may want to re-export to
        # ANSYS (different OUTPUT_NAME, etc.) without redoing the whole
        # calibration. Splitting the pipeline here lets export_geom_ansys.py
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
        #                  >=1=branch) - not needed by write_geom() itself, but
        #                  kept so nothing is lost if you want to recompute
        #                  volume_stats()/report_volume() etc. from the .npz later.
        #   root          : the root node index (scalar) - write_geom() needs
        #                  it to compute the x,y recentring offset.
        #   recenter_xy   : the RECENTER_XY flag used for THIS run (scalar bool).
        #   geom_filename : the `out` filename computed above - so
        #                  export_geom_ansys.py writes the SAME geom_*.txt name
        #                  this script would have used, without recomputing
        #                  OUTPUT_NAME/OUTPUT_PATTERN/variant-suffix logic.
        #   tree_name, variant_label, threshold_m : just metadata, so you can
        #                  tell which run produced a given .npz file later.
        npz_name = "calib_%s_r%dmm%s.npz" % (TREE_NAME, round(thr * 1000), variant_suffix)
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
            print("  AdQSM calibration factors (order: AdQSM_median_radius / AdTree_median_radius):")
            for o in sorted(factors):
                print("    order %d : factor = %.3f" % (o, factors[o]))
            cal_lengths, cal_radii = cylinder_metrics(xyz, cyl)
            cal_stats = volume_stats(cal_lengths, cal_radii, np.asarray(cyl_order))

            # ---- (Part A) thin-branch diagnostic on the CALIBRATED cylinders ----
            cal_thin = report_thin_branch_volume(cal_lengths, cal_radii, cyl_order,
                                                  cut_cm=THIN_BRANCH_CUT_CM)

            # DBH/taper of the CALIBRATED trunk, using the now-replaced radii.
            cal_dbh = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_LOWER)
            cal_d_upper = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_UPPER)
            cal_taper = ((cal_dbh - cal_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                         if cal_dbh is not None and cal_d_upper is not None else None)

            # ---- upsert both the uncalibrated and calibrated rows for this threshold ----
            # "AdTree raw" does NOT depend on AdQSM at all, so it gets no variant
            # suffix - it's simply re-written (with identical values) for every
            # variant, which is harmless since upsert_result overwrites by
            # (tree, method), not duplicates.
            # branch_filter = "none": raw AdTree radii, no diameter cut-off applied.
            upsert_result(RESULTS_CSV, TREE_NAME, "AdTree raw r%dmm" % round(thr * 1000),
                          orig_stats["total_vol"], orig_stats["trunk_vol"], orig_stats["branch_vol"], None,
                          raw_dbh, height_m, raw_taper,
                          # trunk_len/branch_len: already in this dict (volume_stats()
                          # computes them the same way as trunk_vol/branch_vol).
                          orig_stats["trunk_len"], orig_stats["branch_len"],
                          branch_filter="none")
            # "AdTree calibrated" DOES depend on which AdQSM variant it was
            # calibrated against, so it gets the variant suffix (when set).
            # NOTE: radius calibration only rescales/replaces RADII, never the
            # cylinder geometry itself, so cal_stats's lengths equal orig_stats's -
            # calibration cannot change how much length was reconstructed.
            # branch_filter = "none": calibrated radii, but still the full
            # reconstruction - no diameter cut-off applied (see the THIRD row
            # right below for the "10cm"-filtered counterpart of this one).
            upsert_result(RESULTS_CSV, TREE_NAME,
                          "AdTree calibrated r%dmm%s" % (round(thr * 1000), variant_method_suffix),
                          cal_stats["total_vol"], cal_stats["trunk_vol"], cal_stats["branch_vol"], None,
                          cal_dbh, height_m, cal_taper,
                          cal_stats["trunk_len"], cal_stats["branch_len"],
                          branch_filter="none")
            # Optional THIRD row: same tree/DBH/height/taper, but total/stem/branch
            # volume restricted to cylinders >= THIN_BRANCH_CUT_CM (like TreeQSM's
            # "...Filtered..." rows) - lets compare_volumes.py/plot_volumes.py show
            # the reference vs. an apples-to-apples "same cut-off" comparison.
            if WRITE_THIN_BRANCH_FILTERED_ROW:
                # branch_filter = "10cm": this row IS the diameter-cut-off variant
                # (its name already says "(>=10cm only)") - trunk_len/branch_len
                # are left as None here since report_thin_branch_volume() doesn't
                # currently track filtered LENGTH, only filtered volume (see the
                # summary from when this row was added).
                upsert_result(RESULTS_CSV, TREE_NAME,
                              "AdTree calibrated r%dmm%s (>=%.0fcm only)"
                              % (round(thr * 1000), variant_method_suffix, THIN_BRANCH_CUT_CM),
                              cal_thin["total_vol_kept"], cal_thin["trunk_vol_kept"],
                              cal_thin["branch_vol_kept"], None,
                              cal_dbh, height_m, cal_taper,
                              branch_filter="10cm")

            print("  DBH (at %.1f m)   : raw AdTree = %s   |   calibrated = %s"
                  % (TAPER_H_LOWER, _fmt_dbh(raw_dbh), _fmt_dbh(cal_dbh)))
            print("  Taper (%.1f-%.1f m): raw AdTree = %s   |   calibrated = %s"
                  % (TAPER_H_LOWER, TAPER_H_UPPER, _fmt_taper(raw_taper), _fmt_taper(cal_taper)))
            print("  Height (pruned model): %s" % (("%.2f m" % height_m) if height_m is not None else "n/a"))

            print("  Volume comparison (a) raw skeleton vs. (b) processed/AdTree vs. (c) processed/calibrated:")
            print_volume_stats("(a) raw skeleton (AdTree)", raw_stats)
            print_volume_stats("(b) processed, AdTree radii", orig_stats)
            print_volume_stats("(c) processed, CALIBRATED", cal_stats)
        print()

        if SHOW_PLOT or SAVE_PLOT_PNG:
            # `out` is only used here to derive the PNG name (out.txt -> out.png) -
            # plot_model() never writes `out` itself, so this is unaffected by
            # geom_*.txt no longer being written directly by this script.
            plot_model(xyz, cyl, root, RECENTER_XY, thr, out, SHOW_PLOT, SAVE_PLOT_PNG)
