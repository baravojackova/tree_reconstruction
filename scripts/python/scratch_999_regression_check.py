# -*- coding: utf-8 -*-
# =====================================================================
#  THROWAWAY regression-check script (explicitly requested to be left in
#  place, not deleted, by the task it was written for).
#
#  Does NOT modify adqsm_build_median_variant.py or tree_geom_utils.py -
#  imports both as MODULES (importing adqsm_build_median_variant.py this
#  way does NOT execute its run() RUN section, which only fires under
#  `if __name__ == "__main__":` at the bottom of that file) and calls
#  their existing functions directly, unchanged.
#
#  PURPOSE: for a set of trees whose "999" median variant already exists
#  (built under the OLD, per-variant spike-filtering code, before the
#  cross-variant union-based rewrite), recompute taper.txt with the NEW
#  logic into scratch_999_regression/<TREE>/taper.txt - a folder created
#  next to this script, OUTSIDE data/ entirely - and compare it,
#  height-by-height, against the EXISTING data/<TREE>/999/taper.txt.
#
#  Writes ONLY under scratch_999_regression/. Never writes into data/,
#  never modifies any existing taper.txt/BranchStructure.txt/
#  TreesParams.txt (those are only ever opened here in "r" mode). Safe to
#  re-run any time - each run just overwrites its own scratch output.
# =====================================================================

import os

import numpy as np

import adqsm_build_median_variant as m
from tree_geom_utils import parse_adqsm_taper_file, _reject_taper_spikes

SCRATCH_DIR = "scratch_999_regression"

# Trees to check + their ACTUAL source-variant list on disk, confirmed by
# a prior read-only inventory (NOT assumed to be identical for every
# tree - see that inventory's own findings).
#
# B21_S01 also has an EXTRA "099" variant folder alongside the usual 7 -
# confirmed by the user to be a leftover, NOT part of the original 999
# build, so it is deliberately left out of this list (and never read by
# this script at all - _variant_dir()/build_median_taper() below only
# ever touch the 7 folders actually listed for each tree).
TREES = {
    "B21_S01": ["040", "050", "060", "070", "080", "090", "100"],
    "B21_S04": ["040", "050", "060", "070", "080", "090", "100"],
    "B21_S07": ["040", "050", "060", "070", "080", "090", "100"],
    "B21_S08": ["040", "050", "060", "070", "080", "090", "100"],
    "B21_S09": ["040", "050", "060", "070", "080", "090", "100"],
}


def _ratio_map(tree_name, variant):
    """{height: ratio} for every INTERIOR row of one source variant's raw
    (unfiltered) taper.txt, where ratio = diameter[i] / max(neighbour
    diameters) - the SAME quantity _reject_taper_spikes() compares
    against SPIKE_FACTOR. Computed here ONLY for this report's printed
    numbers - never used to decide anything: the actual pass/fail
    decision always comes from calling _reject_taper_spikes() itself
    (see check_tree() below), not from this helper, so there is still
    only ONE definition of "what counts as a spike" in the codebase."""
    path = os.path.join(m.DATA_ROOT, tree_name, variant, "taper.txt")
    h, d = parse_adqsm_taper_file(path, raw=True)
    ratios = {}
    for i in range(1, len(d) - 1):
        neighbor_max = max(d[i - 1], d[i + 1])
        ratios[h[i]] = (d[i] / neighbor_max) if neighbor_max > 0 else float("nan")
    return ratios


def check_tree(tree_name, source_variants):
    print("=" * 90)
    print("Tree:", tree_name, " source variants:", source_variants)
    print("-" * 90)

    out_dir = os.path.join(SCRATCH_DIR, tree_name)
    os.makedirs(out_dir, exist_ok=True)

    # ---- recompute with the NEW logic, into scratch (never data/) -----
    ref_h, median_diam, _per_variant, new_path = m.build_median_taper(tree_name, source_variants, out_dir)
    print("  (scratch taper.txt written to: %s)" % new_path)

    # ---- per-variant failing heights + ratio, for this report ---------
    failing_per_variant = {}
    ratio_per_variant = {}
    for v in source_variants:
        path = os.path.join(m.DATA_ROOT, tree_name, v, "taper.txt")
        h, d = parse_adqsm_taper_file(path, raw=True)
        failing = _reject_taper_spikes(h, d, path, factor=m.SPIKE_FACTOR, detect_only=True)
        failing_per_variant[v] = set(failing.tolist())
        ratio_per_variant[v] = _ratio_map(tree_name, v)

    union_heights = set()
    for v in source_variants:
        union_heights |= failing_per_variant[v]

    print()
    print("  Union of dropped heights under NEW logic (SPIKE_FACTOR=%.3g): %s"
          % (m.SPIKE_FACTOR, sorted(union_heights) if union_heights else "(none)"))
    for h in sorted(union_heights):
        failed_in = [v for v in source_variants if h in failing_per_variant[v]]
        kept_in = [v for v in source_variants if v not in failed_in]
        ratios_str = ", ".join("%s=%.4g" % (v, ratio_per_variant[v].get(h, float("nan"))) for v in source_variants)
        print("    height %.4g m: FAILED in %s; kept in %s" % (h, failed_in, kept_in if kept_in else "(none)"))
        print("      ratios (all variants): %s" % ratios_str)
        if m.TAPER_H_LOWER <= h <= m.TAPER_H_UPPER:
            print("      *** inside TAPER_H_LOWER..TAPER_H_UPPER (%.1f-%.1f m) ***"
                  % (m.TAPER_H_LOWER, m.TAPER_H_UPPER))

    # ---- load the EXISTING (old-code) result for comparison -----------
    old_path = os.path.join(m.DATA_ROOT, tree_name, "999", "taper.txt")
    h_old, d_old = parse_adqsm_taper_file(old_path, raw=True)
    h_new, d_new = ref_h, median_diam

    old_map = dict(zip(h_old.tolist(), d_old.tolist()))
    new_map = dict(zip(h_new.tolist(), d_new.tolist()))
    all_heights = sorted(set(old_map) | set(new_map))

    changed_rows = []   # (height, old_value, new_value) - present in both, value differs
    dropped_now = []    # present in old, missing in new
    added_now = []      # present in new, missing in old (should be impossible; checked anyway)
    for h in all_heights:
        in_old, in_new = h in old_map, h in new_map
        if in_old and in_new:
            ov, nv = old_map[h], new_map[h]
            if not np.isclose(ov, nv):
                changed_rows.append((h, ov, nv))
        elif in_old and not in_new:
            dropped_now.append(h)
        elif in_new and not in_old:
            added_now.append(h)

    identical = not changed_rows and not dropped_now and not added_now
    print()
    print("  RESULT: %s" % ("IDENTICAL" if identical else "DIFFERENT"))

    in_window_flags = []
    if not identical:
        if changed_rows:
            print("  Changed heights (old -> new, abs diff, pct diff):")
            for h, ov, nv in changed_rows:
                absd = nv - ov
                pct = (absd / ov * 100.0) if ov != 0 else float("inf")
                in_win = m.TAPER_H_LOWER <= h <= m.TAPER_H_UPPER
                in_window_flags.append(in_win)
                flag = "  *** INSIDE 1.3-10.0 m ***" if in_win else ""
                print("    h=%.4g m: old=%.6f  new=%.6f  absdiff=%.6f  pctdiff=%.3f%%%s"
                      % (h, ov, nv, absd, pct, flag))
        if dropped_now:
            print("  Heights dropped under NEW logic that were NOT dropped before:")
            for h in dropped_now:
                in_win = m.TAPER_H_LOWER <= h <= m.TAPER_H_UPPER
                in_window_flags.append(in_win)
                flag = "  *** INSIDE 1.3-10.0 m ***" if in_win else ""
                print("    h=%.4g m%s" % (h, flag))
        if added_now:
            print("  Heights present under NEW logic that were ABSENT before (unexpected - should not "
                  "be possible when SPIKE_FACTOR only decreased):")
            for h in added_now:
                in_window_flags.append(m.TAPER_H_LOWER <= h <= m.TAPER_H_UPPER)
                print("    h=%.4g m" % h)

    # ---- diameter @ 1.3 m, old and new, regardless of identical --------
    i_old = int(np.argmin(np.abs(h_old - m.TAPER_H_LOWER)))
    i_new = int(np.argmin(np.abs(h_new - m.TAPER_H_LOWER)))
    print()
    print("  diameter @ 1.3 m: OLD=%.6f m   NEW=%.6f m" % (d_old[i_old], d_new[i_new]))
    print("=" * 90)
    print()

    any_in_window = any(in_window_flags)
    return identical, len(changed_rows), len(dropped_now) + len(added_now), len(union_heights), any_in_window


if __name__ == "__main__":
    summary = []
    for tree_name, source_variants in TREES.items():
        identical, n_changed, n_row_set_diff, n_dropped_union, in_window = check_tree(tree_name, source_variants)
        summary.append((tree_name, identical, n_changed, n_row_set_diff, n_dropped_union, in_window))

    print("=" * 90)
    print("SUMMARY")
    print("%-10s %-11s %-16s %-16s %-20s" % ("tree", "identical", "heights_changed", "heights_dropped", "change_in_1.3-10.0m"))
    for tree_name, identical, n_changed, n_row_set_diff, n_dropped_union, in_window in summary:
        print("%-10s %-11s %-16d %-16d %-20s" % (tree_name, str(identical), n_changed, n_dropped_union, str(in_window)))
