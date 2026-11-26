# -*- coding: utf-8 -*-
# =====================================================================
#  READ-ONLY diagnostic: does AdQSM's branch order o correspond to the
#  SAME physical category as AdTree's branch order o, or does AdQSM's
#  own fork-detection lag behind AdTree's by some roughly constant shift
#  k as order increases?
#
#  HYPOTHESIS: AdQSM misses some forks going up the crown, so a branch
#  AdTree labels order o+k is labelled order o by AdQSM, with k growing
#  with order. If true, AdQSM's per-order median diameter should align
#  with AdTree's median diameter at a SHIFTED order o+k, not at the same
#  order o - tested here by comparing AdQSM order o against AdTree order
#  o+k for k = 0..6, and checking whether some k > 0 makes the two
#  reconstructions' per-order diameters agree more CONSISTENTLY across
#  orders (lower std of log-ratio), not just closer on average (lower
#  mean of log-ratio) - the latter could just be absorbing an overall
#  scale difference between the two reconstructions, not a labelling fix.
#
#  Does NOT import adtree_reconstruct_compare.py: that module has NO
#  `if __name__ == "__main__":` guard - its RUN section executes
#  unconditionally on import. INPUT_PLY/MERGE_DECIMALS/SMOOTH_ITERS/
#  SMOOTH_ALPHA are hardcoded below to match that file's own values
#  (same convention as scratch_skeleton_diag.py/scratch_calib_sensitivity.py),
#  not imported.
#
#  Writes NOTHING to disk: no CSV upserts, no plots, no .npz files, no
#  calibration of any kind. Console output only. Calls convert() exactly
#  ONCE, with thr=0.0 (no pruning) - no loop, no RADIUS_THRESHOLDS.
#
#  Dependencies: numpy
# =====================================================================

import csv
import os

import numpy as np

from tree_geom_utils import read_ply, merge_vertices, smooth_centerline, convert

# paths.py (same folder as this script) knows where the project lives on
# THIS computer, so no folder path has to be typed in by hand here.
import paths

# =====================  PARAMETERS  ===================================
# TREE_NAME / ADQSM_VARIANT / INPUT_PLY are edited TOGETHER when switching
# trees, run once per tree (no multi-tree loop here - see STEP 3's own
# per-tree runs). INPUT_PLY is a literal, explicit path - NOT built from
# TREE_NAME by string formatting - because the AdTree skeleton filename
# does not follow one convention across trees (see the four candidates
# below: "_noplate_clean_skeleton", " - Cloud_skeleton_adtree", " - Cloud_
# skeleton", " - Cloud_nobutt_skeleton"). A constructed name would silently
# pick the wrong file, or fail outright, on the next tree. The differing
# base/butt treatment these names imply is known, documented and
# intentional - it affects order 0 only, which this script never compares
# (the per-tree order list derived at runtime in run() excludes order 0).
#
# ADQSM_VARIANT choice: order/parent/length/height/angle columns in
# BranchStructure.txt are byte-identical across every HS variant within a
# tree (confirmed for B21_S01, re-confirmed for all three IND trees in the
# STEP 1 inventory), so the variant affects diameters only, never this
# alignment test. "06" (mid-range of 04-10) is used for the three IND
# trees, which have no "999" variant (that median-of-variants build exists
# only for B21_S01). B21_S01 itself keeps "999" so this script reproduces
# the earlier single-tree alignment run exactly.
#
# --- B21_S01 (control) ---
# TREE_NAME = "B21_S01"
# ADQSM_VARIANT = "999"
# INPUT_PLY = (r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction"
#              r"\data\B21_S01\B21_S01_noplate_clean_skeleton.ply")

# --- IND01_054 ---
# TREE_NAME = "IND01_054"
# ADQSM_VARIANT = "06"
# INPUT_PLY = (r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction"
#              r"\data\IND01_054\IND01_054 - Cloud_skeleton_adtree.ply")

# --- IND03_088 ---
# TREE_NAME = "IND03_088"
# ADQSM_VARIANT = "06"
# INPUT_PLY = (r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction"
#              r"\data\IND03_088\IND03_088 - Cloud_skeleton.ply")

# --- IND07_083 - ACTIVE ---
TREE_NAME = "IND07_083"
ADQSM_VARIANT = "06"
INPUT_PLY = str(paths.tree_data_dir("IND07_083") / "IND07_083 - Cloud_nobutt_skeleton.ply")

BRANCH_FILE = os.path.join(str(paths.tree_data_dir(TREE_NAME)), ADQSM_VARIANT, "BranchStructure.txt")

# Copied verbatim from adtree_reconstruct_compare.py's own PARAMETERS block -
# NOT imported (see header comment above for why).
MERGE_DECIMALS = 5
SMOOTH_ITERS = 5
SMOOTH_ALPHA = 0.5

# The single convert() call this script makes - thr=0.0 (unpruned), same
# SEG_LEN_MIN/MAX/K/MIN_CYL_LEN as the "0.065/0.5" case used throughout the
# length-weighting investigation. Identical for every tree, so the
# cross-tree comparison in STEP 3 is apples-to-apples on resampling
# density.
SEG_LEN_MIN = 0.065
SEG_LEN_MAX = 0.5
SEG_LEN_K = 0.5
MIN_CYL_LEN = 0.005

K_SHIFTS = list(range(0, 7))   # k = 0..6
# ADQSM_ORDERS is NOT a fixed constant here - a hardcoded range(1, 8) would
# silently drop IND01_054/IND03_088's order 8 (their own maximum - see STEP
# 1c), computing those trees' test from 7 pairs instead of 8. Derived at
# runtime instead, from whatever orders are actually present in that tree's
# own BranchStructure.txt (excluding order 0, the trunk) - see run(), which
# also prints the resulting list so it's visible which orders went in.

# RESULTS_SUMMARY_PATH: one row appended per k, per run, so all four trees'
# results can be compared without re-running anything. A scratch file in
# this session's own scratchpad directory - NOT volume_results.csv, NOT
# under plots/ - purely this diagnostic's own accumulating record.
# NOTE: this path is deliberately left hard-coded. It points into a temporary
# scratchpad folder of the Claude session OUTSIDE the project (not into
# tree_reconstruction), so paths.py has no equivalent for it.
RESULTS_SUMMARY_PATH = (
    r"C:\Users\Spravce\AppData\Local\Temp\claude\c--Users-Spravce-Documents-BARA-01-Skeny-Babice-tree-reconstruction"
    r"\72accd12-932f-436d-b4d0-4dc437d3c140\scratchpad\order_alignment_summary.csv")
# =====================================================================


def append_results_summary(path, tree_name, adqsm_variant, results):
    """Append one row per (k, tree, adqsm_variant) to RESULTS_SUMMARY_PATH -
    tree, adqsm_variant, k, n_pairs, mean_abs_log_ratio, std_log_ratio.
    Writes the header only if the file doesn't exist yet, so repeated runs
    (one per tree) accumulate into a single comparable table instead of
    overwriting each other."""
    write_header = not os.path.exists(path)
    with open(path, "a", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["tree", "adqsm_variant", "k", "n_pairs",
                              "mean_abs_log_ratio", "std_log_ratio"])
        for k, mean_abs_log, std_log, n_pairs in results:
            writer.writerow([
                tree_name, adqsm_variant, k, n_pairs,
                "" if mean_abs_log is None else "%.6f" % mean_abs_log,
                "" if std_log is None else "%.6f" % std_log,
            ])


def weighted_quantile(values, weights, quantiles):
    """Plain cumulative-weight interpolation (no new dependency) - same
    implementation as scratch_calib_sensitivity.py's own weighted_quantile()
    (itself later ported into tree_geom_utils.py as part of the
    length-weighting fix) - confirmed correct there, reused verbatim here
    rather than reimplemented."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    sorter = np.argsort(values)
    values = values[sorter]
    weights = weights[sorter]
    cum_weights = np.cumsum(weights) - 0.5 * weights
    cum_weights = cum_weights / weights.sum()
    return np.interp(quantiles, cum_weights, values)


def parse_adqsm_branch_structure(path):
    """Per-order (diameter[m], length[m]) lists from BranchStructure.txt -
    same column convention as tree_geom_utils.py's parse_adqsm_branch_file_raw()
    (order=col0, diameter=col2 in metres) plus length=col5 (see the file's own
    header row: order/parent/diameter(m)/volume(L)/area(m^2)/length(m)/...)."""
    diam_by_order = {}
    length_by_order = {}
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 7:
                continue
            try:
                order = int(parts[0])
                diameter = float(parts[2])
                length = float(parts[5])
            except ValueError:
                continue
            diam_by_order.setdefault(order, []).append(diameter)
            length_by_order.setdefault(order, []).append(length)
    return diam_by_order, length_by_order


def run():
    print("=" * 90)
    print("Reading:", INPUT_PLY)
    xyz, rad, edges = read_ply(INPUT_PLY)
    print("  vertices: %d, edges: %d" % (len(xyz), len(edges)))

    xyz, rad, edges = merge_vertices(xyz, rad, edges, MERGE_DECIMALS)
    print("  after merging: %d points, %d edges" % (len(xyz), len(edges)))

    smooth_root = int(np.argmin(xyz[:, 2]))
    xyz = smooth_centerline(xyz, edges, smooth_root, SMOOTH_ITERS, SMOOTH_ALPHA)
    print("  centerline smoothing: %d passes, alpha=%.2f" % (SMOOTH_ITERS, SMOOTH_ALPHA))

    print("\nconvert() - single call, thr=0.0, SEG_LEN_MIN=%.3f, SEG_LEN_MAX=%.3f, SEG_LEN_K=%.2f, MIN_CYL_LEN=%.3f"
          % (SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN))
    root, cyl, cyl_order = convert(xyz, rad, edges, 0.0, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K, MIN_CYL_LEN)
    print("  total cylinders: %d" % len(cyl))

    # ---- STEP 2: AdTree side, per-order length-weighted diameter stats ----
    radii = np.array([c[2] for c in cyl], dtype=float)
    lengths = np.array([float(np.linalg.norm(xyz[b] - xyz[a])) for a, b, r, pid in cyl])
    orders_arr = np.asarray(cyl_order)

    print("\n" + "=" * 90)
    print("STEP 2: AdTree side (length-weighted per-order diameter stats)")
    print("=" * 90)
    print("%-6s %10s %14s %16s %16s %12s %12s" % (
        "order", "n_cyl", "total_len_m", "wmed_diam_m", "wmed_diam_cm", "wp25_cm", "wp75_cm"))

    adtree_stats = {}
    max_order = int(orders_arr.max())
    tot_cyl = 0
    tot_len = 0.0
    for o in range(0, max_order + 1):
        mask = orders_arr == o
        n = int(mask.sum())
        if n == 0:
            continue
        r_o = radii[mask]
        l_o = lengths[mask]
        diam_o = 2.0 * r_o
        wq = weighted_quantile(diam_o, l_o, [0.25, 0.5, 0.75])
        total_len_o = float(l_o.sum())
        adtree_stats[o] = dict(n=n, total_len=total_len_o, wmedian_diam_m=float(wq[1]),
                                wp25_cm=float(wq[0]) * 100, wp75_cm=float(wq[2]) * 100)
        tot_cyl += n
        tot_len += total_len_o
        print("%-6d %10d %14.2f %16.5f %16.3f %12.3f %12.3f" % (
            o, n, total_len_o, wq[1], wq[1] * 100, wq[0] * 100, wq[2] * 100))

    print("\nTOTAL cylinders: %d   TOTAL length: %.2f m" % (tot_cyl, tot_len))

    # These specific expected values are B21_S01's own known-good numbers
    # from the original single-tree alignment run - only meaningful as a
    # cross-check when TREE_NAME is actually B21_S01. Gated so running this
    # script for an IND tree doesn't print misleading "NO MATCH" lines
    # against numbers that were never expected to apply to it.
    if TREE_NAME == "B21_S01":
        print("\nCross-check against Phase A's known values (B21_S01 only):")
        checks = [
            (0, 220, 27.84),
            (5, 43717, 2687.92),
            (7, 54247, 3240.94),
        ]
        for o, exp_n, exp_len in checks:
            got_n = adtree_stats[o]["n"]
            got_len = adtree_stats[o]["total_len"]
            match_n = "match" if got_n == exp_n else "NO MATCH"
            match_len = "match" if abs(got_len - exp_len) < 0.01 else "NO MATCH"
            print("  order %d: n expected=%d got=%d (%s)   total_len expected=%.2f got=%.2f (%s)"
                  % (o, exp_n, got_n, match_n, exp_len, got_len, match_len))
        totals_match_n = "match" if tot_cyl == 319584 else "NO MATCH"
        totals_match_len = "match" if abs(tot_len - 19371.60) < 0.01 else "NO MATCH"
        print("  totals: n expected=319584 got=%d (%s)   total_len expected=19371.60 got=%.2f (%s)"
              % (tot_cyl, totals_match_n, tot_len, totals_match_len))
    else:
        print("\n(Cross-check against known values skipped - those numbers are B21_S01-specific.)")

    # ---- STEP 3: alignment test ----
    adqsm_diam_by_order, adqsm_length_by_order = parse_adqsm_branch_structure(BRANCH_FILE)
    adqsm_median_diam = {o: float(np.median(v)) for o, v in adqsm_diam_by_order.items()}
    adqsm_orders = sorted(o for o in adqsm_diam_by_order if o != 0)

    print("\n" + "=" * 90)
    print("STEP 3: alignment test - AdQSM order o vs. AdTree order o+k")
    print("=" * 90)
    print("AdQSM branch orders present (excl. trunk, order 0): %s" % adqsm_orders)

    results = []
    for k in K_SHIFTS:
        rows = []
        for o in adqsm_orders:
            if o not in adqsm_median_diam:
                continue
            shifted = o + k
            if shifted not in adtree_stats:
                continue
            adqsm_d = adqsm_median_diam[o]
            adtree_d = adtree_stats[shifted]["wmedian_diam_m"]
            ratio = adtree_d / adqsm_d
            rows.append((o, shifted, adqsm_d, adtree_d, ratio))

        if not rows:
            print("\nk=%d: no valid (o, o+k) pairs - skipping" % k)
            results.append((k, None, None, 0))
            continue

        log_ratios = np.array([np.log(r[4]) for r in rows])
        mean_abs_log = float(np.mean(np.abs(log_ratios)))
        std_log = float(np.std(log_ratios))
        results.append((k, mean_abs_log, std_log, len(rows)))

        print("\nk=%d (%d pairs):" % (k, len(rows)))
        print("  %-4s %-8s %14s %14s %10s" % ("o", "o+k", "AdQSM_diam_m", "AdTree_diam_m", "ratio"))
        for o, shifted, adqsm_d, adtree_d, ratio in rows:
            print("  %-4d %-8d %14.5f %14.5f %10.4f" % (o, shifted, adqsm_d, adtree_d, ratio))
        print("  mean|log(ratio)| = %.5f   std(log(ratio)) = %.5f" % (mean_abs_log, std_log))

    print("\n" + "=" * 90)
    print("SUMMARY TABLE - all k, sorted by k")
    print("=" * 90)
    print("%-6s %10s %12s %14s %8s" % ("k", "n_pairs", "mean|logR|", "std(logR)", ""))
    for k, mean_abs_log, std_log, n_pairs in results:
        if mean_abs_log is None:
            print("%-6d %10d %12s %14s" % (k, 0, "n/a", "n/a"))
        else:
            print("%-6d %10d %12.5f %14.5f" % (k, n_pairs, mean_abs_log, std_log))

    valid_results = [r for r in results if r[1] is not None]
    best_mean = min(valid_results, key=lambda r: r[1])
    best_std = min(valid_results, key=lambda r: r[2])
    k0 = next(r for r in results if r[0] == 0)

    print("\nk minimising mean|log(ratio)|: k=%d (mean=%.5f)" % (best_mean[0], best_mean[1]))
    print("k minimising std(log(ratio)):  k=%d (std=%.5f)" % (best_std[0], best_std[2]))
    if k0[1] is not None:
        print("k=0 for reference: mean=%.5f, std=%.5f" % (k0[1], k0[2]))
    print("=" * 90)

    append_results_summary(RESULTS_SUMMARY_PATH, TREE_NAME, ADQSM_VARIANT, results)
    print("\nAppended %d row(s) to %s" % (len(results), RESULTS_SUMMARY_PATH))


if __name__ == "__main__":
    run()
