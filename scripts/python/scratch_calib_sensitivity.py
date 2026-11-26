# -*- coding: utf-8 -*-
# =====================================================================
#  READ-ONLY diagnostic: reproduce adtree_reconstruct_compare.py's
#  pipeline from read_ply() through the per-order regression fit AND the
#  calref=min5mm scalar-factor computation (everything between its first
#  convert() call and the start of its RADIUS_THRESHOLDS loop), for TWO
#  SEG_LEN settings, to distinguish two candidate explanations for why
#  calibrated volume is highly sensitive to SEG_LEN_MIN/K while raw
#  (uncalibrated) volume is not:
#
#   H1 - the FIT itself changes (build_quantile_matched_pairs() pairs
#        AdTree cylinders with AdQSM branch rows; a different AdTree
#        cylinder count/distribution could shift the matched quantiles,
#        hence the fitted (a, b)).
#   H2 - the fit is stable and something in the APPLY step (which is NOT
#        a pure multiplicative rescale - see fit_radius_regression()'s
#        power-law form) is responsible.
#
#  Does NOT import adtree_reconstruct_compare.py: that module has NO
#  `if __name__ == "__main__":` guard (confirmed in the prior
#  scratch_skeleton_diag.py task) - its RUN section executes
#  unconditionally on import. INPUT_PLY/MERGE_DECIMALS/SMOOTH_ITERS/
#  SMOOTH_ALPHA are therefore hardcoded below to match that file's own
#  values (cited inline), not imported.
#
#  Writes NOTHING to disk: no CSV upserts, no plots, no .npz files.
#  Console output only.
#
#  Dependencies: numpy, scipy (via tree_geom_utils's fit_radius_regression)
# =====================================================================

import os

import numpy as np

from tree_geom_utils import (
    read_ply, merge_vertices, smooth_centerline,
    parse_adqsm_branch_file, parse_adqsm_branch_file_raw,
    convert, build_quantile_matched_pairs, fit_radius_regression,
    group_orders_for_fitting, compute_order_calibration_factors,
)

# paths.py (same folder as this script) knows where the project lives on
# THIS computer, so no folder path has to be typed in by hand here.
import paths

# =====================  PARAMETERS  ===================================
TREE_NAME = "B21_S01"
ADQSM_VARIANT = "999"

# (SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K) - the two extremes of the observed
# 9-combination sweep's volume spread (regression-perorder total: 13.246 m3
# at MIN=0.035/K=1.0 vs. 10.228 m3 at MIN=0.065/K=0.5 - see the task's
# OBSERVATION TO EXPLAIN).
SEG_CASES = [
    (0.035, 0.5, 1.0),
    (0.065, 0.5, 0.5),
]

MIN_CYL_LEN = 0.005          # adtree_reconstruct_compare.py PARAMETERS block
MIN_PAIRS_PER_ORDER = 15     # adtree_reconstruct_compare.py PARAMETERS block

# The single calref reference threshold used for the "[calref=min5mm]"
# scalar-factor path in adtree_reconstruct_compare.py
# (CALIBRATION_REF_THRESHOLDS_MM = [0.005] there - copied here, not
# imported, same reasoning as the constants below).
CALREF_THRESHOLD_M = 0.005

DATA_ROOT = str(paths.DATA_DIR)
# Same convention as adtree_reconstruct_compare.py:113/117 (AdTree_DIR,
# INPUT_PLY) and :110-111/314-315 (ADQSM_BASE_DIR, BranchStructure.txt path).
AdTree_DIR = os.path.join(DATA_ROOT, TREE_NAME)
INPUT_PLY = os.path.join(AdTree_DIR, "%s_noplate_clean_skeleton.ply" % TREE_NAME)
ADQSM_BASE_DIR = os.path.join(DATA_ROOT, TREE_NAME)
BRANCH_FILE = os.path.join(ADQSM_BASE_DIR, ADQSM_VARIANT, "BranchStructure.txt")

# Copied verbatim from adtree_reconstruct_compare.py:216 (MERGE_DECIMALS)
# and :231-232 (SMOOTH_ITERS/SMOOTH_ALPHA) - NOT imported, since that file
# has no __main__ guard (see header comment above).
MERGE_DECIMALS = 5
SMOOTH_ITERS = 5
SMOOTH_ALPHA = 0.5

PROBE_RADII_M = [0.005, 0.01, 0.02, 0.05, 0.10, 0.20]
PROBE_ORDERS = [1, 2, 3]

LENGTH_WEIGHT_QUANTILES = [0.10, 0.50, 0.90]   # p10/p50/p90, as fractions
# =====================================================================


def weighted_quantile(values, weights, quantiles):
    """Plain cumulative-weight interpolation (no new dependency): sort by
    value, place each sample at the MIDPOINT of its own cumulative-weight
    interval (Hazen-type), then linearly interpolate the requested
    `quantiles` (fractions in [0, 1]) against that cumulative-weight axis
    via np.interp(). With all weights equal this reduces to the same
    linear-interpolation quantile convention np.percentile() uses, so it is
    directly comparable to the unweighted np.percentile() calls elsewhere in
    this script."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    sorter = np.argsort(values)
    values = values[sorter]
    weights = weights[sorter]
    cum_weights = np.cumsum(weights) - 0.5 * weights
    cum_weights = cum_weights / weights.sum()
    return np.interp(quantiles, cum_weights, values)


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

    print("\nLoading AdQSM data:", BRANCH_FILE)
    adqsm_median_by_order = parse_adqsm_branch_file(BRANCH_FILE)
    raw_diam_by_order = parse_adqsm_branch_file_raw(BRANCH_FILE)
    print("  AdQSM median radius by order: %s"
          % ", ".join("%d=%.4f m" % (o, r) for o, r in sorted(adqsm_median_by_order.items())))

    # Collected per case, for the side-by-side tables printed after both
    # cases have run.
    case_results = []

    for seg_min, seg_max, seg_k in SEG_CASES:
        label = "MIN=%.3f, MAX=%.3f, K=%.2f" % (seg_min, seg_max, seg_k)
        print("\n" + "=" * 90)
        print("CASE:", label)
        print("=" * 90)

        # ---- convert() call #1: unpruned reference (thr=0.0) --------------
        # Same call as adtree_reconstruct_compare.py's ref_root_0/ref_cyl_0/
        # ref_order_0 (its first convert() call).
        ref_root_0, ref_cyl_0, ref_order_0 = convert(
            xyz, rad, edges, 0.0, seg_min, seg_max, seg_k, MIN_CYL_LEN)

        # 2a. Total cylinder count + per-order counts in the reference set.
        ref_orders_arr = np.asarray(ref_order_0)
        print("\n2a. Reference set (thr=0.0) cylinder counts:")
        print("  total: %d" % len(ref_cyl_0))
        order_counts_2a = {}
        for o in sorted(set(ref_orders_arr.tolist())):
            n = int((ref_orders_arr == o).sum())
            order_counts_2a[o] = n
            print("    order %d : %d" % (o, n))

        # ---- HYPOTHESIS TEST (length-weighting): per-order length/volume/
        # weighted-vs-unweighted radius quantiles, over the SAME reference
        # cylinder set (ref_cyl_0/ref_order_0) used above - not the AdQSM-
        # matched subset, so this reproduces build_quantile_matched_pairs()'s
        # own per-order np.percentile() call (unweighted side) exactly.
        ref_radii_0 = np.array([c[2] for c in ref_cyl_0], dtype=float)
        ref_lengths_0 = np.array(
            [float(np.linalg.norm(xyz[b] - xyz[a])) for a, b, r, pid in ref_cyl_0])
        len_stats_by_order = {}
        for o in sorted(set(ref_orders_arr.tolist())):
            mask = ref_orders_arr == o
            radii_o = ref_radii_0[mask]
            lengths_o = ref_lengths_0[mask]
            unweighted_q = np.percentile(radii_o, [q * 100 for q in LENGTH_WEIGHT_QUANTILES])
            weighted_q = weighted_quantile(radii_o, lengths_o, LENGTH_WEIGHT_QUANTILES)
            len_stats_by_order[o] = dict(
                n=int(mask.sum()),
                total_length=float(lengths_o.sum()),
                mean_length=float(lengths_o.mean()) if mask.sum() else float("nan"),
                volume=float(np.sum(np.pi * radii_o ** 2 * lengths_o)),
                unweighted_q=unweighted_q,
                weighted_q=weighted_q,
            )

        # ---- convert() call #2: calref reference (thr=CALREF_THRESHOLD_M) -
        # Same call as adtree_reconstruct_compare.py's ref_root/ref_cyl/
        # ref_order inside "if COMPUTE_CALREF_MIN5MM: for ref_thr in
        # CALIBRATION_REF_THRESHOLDS_MM:".
        calref_root, calref_cyl, calref_order = convert(
            xyz, rad, edges, CALREF_THRESHOLD_M, seg_min, seg_max, seg_k, MIN_CYL_LEN)
        calref_factors = compute_order_calibration_factors(
            calref_cyl, calref_order, adqsm_median_by_order)

        # 2b. AdQSM population size, total and per order - must be IDENTICAL
        # across cases (see STEP 1d): it is read straight from
        # BranchStructure.txt, independent of SEG_LEN_*.
        print("\n2b. AdQSM population size (total and per order) - "
              "must be identical across cases:")
        adqsm_counts = {o: len(v) for o, v in raw_diam_by_order.items()}
        print("  total: %d" % sum(adqsm_counts.values()))
        for o in sorted(adqsm_counts):
            print("    order %d : %d" % (o, adqsm_counts[o]))

        # ---- build_quantile_matched_pairs() + per-order regression fit ----
        adtree_matched, adqsm_matched, order_labels_matched = build_quantile_matched_pairs(
            ref_cyl_0, ref_order_0, raw_diam_by_order)

        # 2c. Number of matched pairs per order.
        print("\n2c. Matched pairs per order:")
        pair_counts = {}
        for o in sorted(set(order_labels_matched.tolist())):
            n = int((order_labels_matched == o).sum())
            pair_counts[o] = n
            print("    order %d : %d" % (o, n))

        order_to_group = group_orders_for_fitting(order_labels_matched, MIN_PAIRS_PER_ORDER)

        order_to_ab = {}
        for group_orders in dict.fromkeys(order_to_group.values()):
            group_mask = np.isin(order_labels_matched, list(group_orders))
            g_a, g_b = fit_radius_regression(adtree_matched[group_mask], adqsm_matched[group_mask])
            for o in group_orders:
                order_to_ab[o] = (g_a, g_b)

        # 2d. Per-order AdTree/AdQSM radius quantiles (p10/p50/p90) among the
        # matched pairs.
        quantiles_by_order = {}
        for o in sorted(set(order_labels_matched.tolist())):
            mask = order_labels_matched == o
            at = adtree_matched[mask]
            aq = adqsm_matched[mask]
            quantiles_by_order[o] = {
                "adtree_p10": float(np.percentile(at, 10)),
                "adtree_p50": float(np.percentile(at, 50)),
                "adtree_p90": float(np.percentile(at, 90)),
                "adqsm_p10": float(np.percentile(aq, 10)),
                "adqsm_p50": float(np.percentile(aq, 50)),
                "adqsm_p90": float(np.percentile(aq, 90)),
            }

        case_results.append(dict(
            label=label, seg_min=seg_min, seg_max=seg_max, seg_k=seg_k,
            order_counts=order_counts_2a, adqsm_counts=adqsm_counts,
            pair_counts=pair_counts, order_to_group=order_to_group,
            order_to_ab=order_to_ab, quantiles_by_order=quantiles_by_order,
            calref_factors=calref_factors, len_stats_by_order=len_stats_by_order,
        ))

    # =====================================================================
    # Side-by-side comparisons across the two cases
    # =====================================================================
    print("\n" + "=" * 90)
    print("SIDE-BY-SIDE COMPARISON")
    print("=" * 90)

    c0, c1 = case_results[0], case_results[1]

    print("\n2d. Matched-pair radius quantiles per order, both cases side by side:")
    all_orders = sorted(set(c0["quantiles_by_order"]) | set(c1["quantiles_by_order"]))
    for o in all_orders:
        q0 = c0["quantiles_by_order"].get(o)
        q1 = c1["quantiles_by_order"].get(o)
        print("  order %d:" % o)
        print("    case A (%s): n_pairs=%s" % (c0["label"], c0["pair_counts"].get(o)))
        if q0:
            print("      AdTree p10/p50/p90 = %.5f / %.5f / %.5f m"
                  % (q0["adtree_p10"], q0["adtree_p50"], q0["adtree_p90"]))
            print("      AdQSM  p10/p50/p90 = %.5f / %.5f / %.5f m"
                  % (q0["adqsm_p10"], q0["adqsm_p50"], q0["adqsm_p90"]))
        print("    case B (%s): n_pairs=%s" % (c1["label"], c1["pair_counts"].get(o)))
        if q1:
            print("      AdTree p10/p50/p90 = %.5f / %.5f / %.5f m"
                  % (q1["adtree_p10"], q1["adtree_p50"], q1["adtree_p90"]))
            print("      AdQSM  p10/p50/p90 = %.5f / %.5f / %.5f m"
                  % (q1["adqsm_p10"], q1["adqsm_p50"], q1["adqsm_p90"]))
        if q0:
            shift_pct = 100.0 * (q1["adtree_p50"] - q0["adtree_p50"]) / q0["adtree_p50"] if q1 and q0["adtree_p50"] else None
            if shift_pct is not None:
                print("      AdTree p50 shift (B vs A): %.2f %%" % shift_pct)

    print("\n2e. Fitted (a, b) per order, both cases side by side, + ratio of a:")
    all_orders_ab = sorted(set(c0["order_to_ab"]) | set(c1["order_to_ab"]))
    for o in all_orders_ab:
        ab0 = c0["order_to_ab"].get(o)
        ab1 = c1["order_to_ab"].get(o)
        if ab0 and ab1:
            a0, b0 = ab0
            a1, b1 = ab1
            print("  order %d: case A a=%.5f b=%.4f  |  case B a=%.5f b=%.4f  |  a ratio (B/A) = %.4f"
                  % (o, a0, b0, a1, b1, a1 / a0 if a0 else float("nan")))
        else:
            print("  order %d: case A=%s  case B=%s" % (o, ab0, ab1))

    print("\n2f. Calibrated radius from each case's fit, at probe radii, for orders %s:" % PROBE_ORDERS)
    for o in PROBE_ORDERS:
        ab0 = c0["order_to_ab"].get(o)
        ab1 = c1["order_to_ab"].get(o)
        print("  order %d:" % o)
        if not ab0 or not ab1:
            print("    (missing fit for this order in at least one case: A=%s B=%s)" % (ab0, ab1))
            continue
        a0, b0 = ab0
        a1, b1 = ab1
        for r in PROBE_RADII_M:
            rc0 = a0 * r ** b0
            rc1 = a1 * r ** b1
            ratio = rc1 / rc0 if rc0 else float("nan")
            print("    r=%.3f m -> case A: %.5f m   case B: %.5f m   ratio (B/A) = %.4f"
                  % (r, rc0, rc1, ratio))

    print("\n2g. calref=min5mm scalar factor per order, both cases side by side, + ratio:")
    all_orders_cf = sorted(set(c0["calref_factors"]) | set(c1["calref_factors"]))
    for o in all_orders_cf:
        f0 = c0["calref_factors"].get(o)
        f1 = c1["calref_factors"].get(o)
        if f0 is not None and f1 is not None:
            print("  order %d: case A factor=%.5f  case B factor=%.5f  ratio (B/A) = %.4f"
                  % (o, f0, f1, f1 / f0 if f0 else float("nan")))
        else:
            print("  order %d: case A=%s  case B=%s" % (o, f0, f1))

    # =====================================================================
    # HYPOTHESIS TEST: does length-weighting the per-order radius quantiles
    # (over the full reference cylinder set, not the AdQSM-matched subset)
    # remove the cross-case shift that the UNWEIGHTED quantiles show?
    # =====================================================================
    print("\n" + "=" * 90)
    print("HYPOTHESIS TEST: length-weighted vs. unweighted radius quantiles")
    print("=" * 90)

    all_orders_len = sorted(set(c0["len_stats_by_order"]) | set(c1["len_stats_by_order"]))

    print("\nA. Per-order AdTree radius quantiles (p10/p50/p90), unweighted vs. "
          "length-weighted, both cases + B/A ratio per quantile:")
    for o in all_orders_len:
        s0 = c0["len_stats_by_order"].get(o)
        s1 = c1["len_stats_by_order"].get(o)
        print("  order %d:" % o)
        if not s0 or not s1:
            print("    (missing in at least one case: A=%s B=%s)" % (bool(s0), bool(s1)))
            continue
        uq0, uq1 = s0["unweighted_q"], s1["unweighted_q"]
        wq0, wq1 = s0["weighted_q"], s1["weighted_q"]
        print("    UNWEIGHTED   case A p10/p50/p90 = %.5f / %.5f / %.5f m"
              % (uq0[0], uq0[1], uq0[2]))
        print("    UNWEIGHTED   case B p10/p50/p90 = %.5f / %.5f / %.5f m"
              % (uq1[0], uq1[1], uq1[2]))
        print("    UNWEIGHTED   B/A ratio p10/p50/p90 = %.4f / %.4f / %.4f"
              % tuple(uq1[i] / uq0[i] if uq0[i] else float("nan") for i in range(3)))
        print("    LEN-WEIGHTED case A p10/p50/p90 = %.5f / %.5f / %.5f m"
              % (wq0[0], wq0[1], wq0[2]))
        print("    LEN-WEIGHTED case B p10/p50/p90 = %.5f / %.5f / %.5f m"
              % (wq1[0], wq1[1], wq1[2]))
        print("    LEN-WEIGHTED B/A ratio p10/p50/p90 = %.4f / %.4f / %.4f"
              % tuple(wq1[i] / wq0[i] if wq0[i] else float("nan") for i in range(3)))

    print("\nB. Per-order cylinder count, total length, mean length, both cases:")
    print("   (sanity check: total_length B/A ratio should be close to 1 - same "
          "geometry, different discretisation)")
    for o in all_orders_len:
        s0 = c0["len_stats_by_order"].get(o)
        s1 = c1["len_stats_by_order"].get(o)
        if not s0 or not s1:
            print("  order %d: missing in at least one case" % o)
            continue
        len_ratio = s1["total_length"] / s0["total_length"] if s0["total_length"] else float("nan")
        print("  order %d:" % o)
        print("    case A: n=%d  total_length=%.2f m  mean_length=%.4f m"
              % (s0["n"], s0["total_length"], s0["mean_length"]))
        print("    case B: n=%d  total_length=%.2f m  mean_length=%.4f m"
              % (s1["n"], s1["total_length"], s1["mean_length"]))
        print("    total_length B/A ratio = %.4f" % len_ratio)

    print("\nC. Per-order volume contribution sum(pi*r^2*L), both cases + B/A ratio:")
    print("   (control: should be close to 1 for every order - confirms the "
          "geometry itself is invariant)")
    for o in all_orders_len:
        s0 = c0["len_stats_by_order"].get(o)
        s1 = c1["len_stats_by_order"].get(o)
        if not s0 or not s1:
            print("  order %d: missing in at least one case" % o)
            continue
        vol_ratio = s1["volume"] / s0["volume"] if s0["volume"] else float("nan")
        print("  order %d: case A volume=%.6f m3   case B volume=%.6f m3   B/A ratio = %.4f"
              % (o, s0["volume"], s1["volume"], vol_ratio))

    print("\n" + "=" * 90)


if __name__ == "__main__":
    run()
