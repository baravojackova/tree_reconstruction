# -*- coding: utf-8 -*-
# =====================================================================
#  READ-ONLY diagnostic: reproduce adtree_reconstruct_compare.py's
#  input-preparation stage EXACTLY (read_ply -> merge_vertices ->
#  smooth_centerline), then STOP before convert() - no resampling is
#  ever run here - to measure two numbers needed before choosing a
#  SEG_LEN_MIN / SEG_LEN_K sweep grid for the beech production trees:
#
#   (A) skeleton edge-length spacing, to check whether a candidate
#       SEG_LEN_MIN is still above the skeleton's own resolution or
#       already below it (in which case resampling keeps every vertex
#       and that grid value is a duplicate run, not a distinct one).
#   (B) the RAW, UNCALIBRATED AdTree radius on the trunk - the
#       resampling clamp (local_seg_len(), tree_geom_utils.py) uses the
#       AdTree radius, not the measured DBH and not AdQSM's calibrated
#       radius, so the break diameter must be compared against THIS
#       radius, not either of those other two.
#
#  Does NOT import adtree_reconstruct_compare.py: that module has NO
#  `if __name__ == "__main__":` guard (confirmed by inspection) - its
#  entire RUN section (real file reads, convert(), file writes) runs
#  unconditionally at import time, which this read-only diagnostic must
#  never trigger. MERGE_DECIMALS/SMOOTH_ITERS/SMOOTH_ALPHA below are
#  therefore hardcoded to match that file's own values, with the exact
#  source lines cited in each comment, rather than imported.
#
#  Writes NOTHING to disk - console output only. Does not call
#  convert() (tree_geom_utils.py) at all.
#
#  Dependencies: numpy (install: pip install numpy)
# =====================================================================

import os

import numpy as np

from tree_geom_utils import read_ply, merge_vertices, smooth_centerline

# =====================  PARAMETERS  ===================================
TREE_NAME = "B21_S01"

DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"
# Same convention as adtree_reconstruct_compare.py:104 (AdTree_DIR =
# os.path.join(DATA_ROOT, TREE_NAME)) and :108 (INPUT_PLY).
AdTree_DIR = os.path.join(DATA_ROOT, TREE_NAME)
INPUT_PLY = os.path.join(AdTree_DIR, "%s_noplate_clean_skeleton.ply" % TREE_NAME)

# Copied verbatim from adtree_reconstruct_compare.py:207/222/223 (NOT
# imported - see this file's own header comment for why).
MERGE_DECIMALS = 5
SMOOTH_ITERS = 5
SMOOTH_ALPHA = 0.5

DBH_HEIGHT_M = 1.3
DBH_BAND_M = 0.10        # +/- band around DBH_HEIGHT_M

SEG_LEN_MIN_GRID = [0.01, 0.05, 0.10, 0.20]
SEG_LEN_K_GRID = [0.2, 0.5, 1.0]
# =====================================================================


def _percentiles(values, ps=(0, 1, 5, 25, 50, 75, 95, 99, 100)):
    """{p: value} for each percentile in `ps` (0/100 are min/max)."""
    return {p: float(np.percentile(values, p)) for p in ps}


def _print_distribution(label, values, ps=(0, 1, 5, 25, 50, 75, 95, 99, 100)):
    values = np.asarray(values, dtype=np.float64)
    pct = _percentiles(values, ps)
    names = {0: "min", 1: "p1", 5: "p5", 25: "p25", 50: "median", 75: "p75",
             95: "p95", 99: "p99", 100: "max"}
    parts = ["count=%d" % len(values)]
    for p in ps:
        parts.append("%s=%.5f" % (names.get(p, "p%d" % p), pct[p]))
    parts.append("mean=%.5f" % float(np.mean(values)))
    print("  %s: %s" % (label, "  ".join(parts)))


def run():
    print("=" * 90)
    print("Reading:", INPUT_PLY)
    xyz, rad, edges = read_ply(INPUT_PLY)
    n_v_before, n_e_before = len(xyz), len(edges)

    # ---- 2a: vertex/edge counts before/after merge_vertices() -----------------
    xyz, rad, edges = merge_vertices(xyz, rad, edges, MERGE_DECIMALS)
    n_v_after, n_e_after = len(xyz), len(edges)
    print("2a. Vertex/edge counts:")
    print("  before merge_vertices(): %d vertices, %d edges" % (n_v_before, n_e_before))
    print("  after  merge_vertices(): %d vertices, %d edges" % (n_v_after, n_e_after))
    print()

    def edge_lengths(xyz_arr):
        d = xyz_arr[edges[:, 0]] - xyz_arr[edges[:, 1]]
        return np.sqrt((d ** 2).sum(axis=1))

    # ---- 2b: edge length distribution, before and after smoothing -------------
    lengths_before_smooth = edge_lengths(xyz)
    print("2b. Edge length distribution [m]:")
    _print_distribution("BEFORE smoothing", lengths_before_smooth)

    smooth_root = int(np.argmin(xyz[:, 2]))
    xyz = smooth_centerline(xyz, edges, smooth_root, SMOOTH_ITERS, SMOOTH_ALPHA)
    print("  (applied smooth_centerline(): %d passes, alpha=%.2f - moves xyz only, "
          "rad is untouched, per tree_geom_utils.py's own docstring)" % (SMOOTH_ITERS, SMOOTH_ALPHA))
    lengths_after_smooth = edge_lengths(xyz)
    _print_distribution("AFTER  smoothing", lengths_after_smooth)
    print()

    # ---- 2c: % of edges shorter than each SEG_LEN_MIN_GRID value ---------------
    print("2c. %% of skeleton edges SHORTER than each SEG_LEN_MIN candidate "
          "(post-smoothing lengths) - answers (A):")
    n_edges = len(lengths_after_smooth)
    for seg_min in SEG_LEN_MIN_GRID:
        pct = 100.0 * np.count_nonzero(lengths_after_smooth < seg_min) / n_edges
        print("  SEG_LEN_MIN=%.2f m: %6.2f %% of edges are already shorter than this value"
              % (seg_min, pct))
    print()

    # ---- 2d: radius distribution over all vertices -----------------------------
    print("2d. Radius distribution over all skeleton vertices [m] (RAW AdTree radii):")
    n_bad = int(np.count_nonzero((rad <= 0) | np.isnan(rad)))
    _print_distribution("radius", rad, ps=(0, 5, 25, 50, 75, 95, 100))
    print("  vertices with radius <= 0 or NaN: %d" % n_bad)
    print()

    # ---- 2e: trunk radius at breast height -------------------------------------
    z_base = float(xyz[:, 2].min())   # same convention as adtree_reconstruct_compare.py:365 (post-smoothing xyz)
    height = xyz[:, 2] - z_base
    band_mask = (height >= DBH_HEIGHT_M - DBH_BAND_M) & (height <= DBH_HEIGHT_M + DBH_BAND_M)
    n_band = int(np.count_nonzero(band_mask))
    print("2e. Trunk radius at breast height (%.2f +/- %.2f m above the lowest skeleton vertex):"
          % (DBH_HEIGHT_M, DBH_BAND_M))
    print("  NOTE: this is the RAW AdTree skeleton radius - NOT calibrated against AdQSM, "
          "and NOT the measured/destructive-reference DBH.")
    print("  vertices in band: %d" % n_band)
    if n_band == 0:
        print("  (no vertices in this height band - cannot report radius stats or a trunk candidate)")
    else:
        band_rad = rad[band_mask]
        r_min, r_med, r_max = float(band_rad.min()), float(np.median(band_rad)), float(band_rad.max())
        print("  radius min/median/max [m] : %.5f / %.5f / %.5f" % (r_min, r_med, r_max))
        print("  diameter min/median/max [cm]: %.3f / %.3f / %.3f" % (2 * r_min * 100, 2 * r_med * 100, 2 * r_max * 100))
        trunk_idx = int(np.argmax(band_rad))
        trunk_radius = float(band_rad[trunk_idx])
        print("  TRUNK CANDIDATE (largest radius in band): radius=%.5f m (diameter=%.3f cm) "
              "- RAW AdTree, uncalibrated" % (trunk_radius, 2 * trunk_radius * 100))
    print()

    # ---- 2f: break-point table ---------------------------------------------------
    print("2f. Break-point table (sorted by break diameter ascending):")
    print("  Sanity check: SEG_LEN_MIN=0.10, SEG_LEN_K=0.5 must give break radius=0.20 m, "
          "break diameter=0.40 m.")
    combos = []
    for seg_min in SEG_LEN_MIN_GRID:
        for seg_k in SEG_LEN_K_GRID:
            break_radius = seg_min / seg_k
            break_diam_m = 2 * break_radius
            pct_below = 100.0 * np.count_nonzero(rad < break_radius) / len(rad)
            combos.append((seg_min, seg_k, break_radius, break_diam_m, pct_below))
    combos.sort(key=lambda c: c[3])

    sanity_ok = None
    for seg_min, seg_k, break_radius, break_diam_m, pct_below in combos:
        flag = ""
        if abs(seg_min - 0.10) < 1e-12 and abs(seg_k - 0.5) < 1e-12:
            sanity_ok = abs(break_radius - 0.20) < 1e-9 and abs(break_diam_m - 0.40) < 1e-9
            flag = "  <-- SANITY CHECK ROW"
        print("  SEG_LEN_MIN=%.2f  SEG_LEN_K=%.2f  ->  break radius=%.4f m  "
              "break diameter=%.4f m (%.2f cm)  %% vertices below break radius=%6.2f %%%s"
              % (seg_min, seg_k, break_radius, break_diam_m, break_diam_m * 100, pct_below, flag))
    if sanity_ok is None:
        print("  SANITY CHECK: SEG_LEN_MIN=0.10/SEG_LEN_K=0.5 combination not found in the grids above - cannot check.")
    elif sanity_ok:
        print("  SANITY CHECK PASSED: SEG_LEN_MIN=0.10/SEG_LEN_K=0.5 -> break radius=0.20 m, break diameter=0.40 m, as expected.")
    else:
        print("  SANITY CHECK FAILED: SEG_LEN_MIN=0.10/SEG_LEN_K=0.5 did NOT give break radius=0.20 m / "
              "break diameter=0.40 m - the formula is wrong, see values printed above.")
    print()

    # ---- 2g: trunk on floor or linear region, per combination -------------------
    print("2g. Trunk candidate vs. break radius, per (SEG_LEN_MIN, SEG_LEN_K) combination:")
    if n_band == 0:
        print("  (no trunk candidate available from 2e - skipping)")
    else:
        for seg_min, seg_k, break_radius, break_diam_m, _pct_below in combos:
            on_floor = trunk_radius < break_radius
            print("  SEG_LEN_MIN=%.2f  SEG_LEN_K=%.2f  break_radius=%.4f m  trunk_radius=%.5f m  -> %s"
                  % (seg_min, seg_k, break_radius, trunk_radius,
                     "ON THE FLOOR (trunk_radius < break_radius)" if on_floor else "linear region"))
    print("=" * 90)


if __name__ == "__main__":
    run()
