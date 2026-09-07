# -*- coding: utf-8 -*-
# =====================================================================
#  Diagnostic: overlay EVERY available AdQSM variant's taper.txt curve
#  for a given tree on one chart, so DBH/trunk-volume differences between
#  variants can be inspected visually and numerically - reusable later on
#  the 12 production beech trees, not just the 3 reference trees, so
#  variant discovery is AUTOMATIC (scan the tree's data folder for
#  subfolders containing a taper.txt), never a hand-typed list.
# ---------------------------------------------------------------------
#  WHY this exists: trunk (order 0) cylinder radius is calibrated ENTIRELY
#  from AdQSM's own taper.txt curve (see adtree_reconstruct_compare.py's
#  own comment near trunk_radius_func, and TODO_investigations.md item 8
#  for a confirmed case of a bad taper.txt inflating calibrated trunk
#  volume). Different AdQSM variants can report meaningfully different
#  taper curves for the SAME tree - this script makes that difference
#  visible at a glance, and prints the numbers (DBH, TrunkVolume) needed
#  to judge which variant's export looks most trustworthy.
#
#  This script does NOT rescale anything (field_dbh=None is passed to
#  make_trunk_radius_func() for every variant) - the whole point is to
#  compare variants' RAW curves against each other; rescaling any one of
#  them to match a measured DBH would hide exactly the difference this
#  script exists to show. A measured DBH, if given (MEASURED_DBH_M), is
#  drawn only as a horizontal reference line - it never feeds back into
#  any curve shown here.
#
#  Reuses (imports, does not duplicate) parse_adqsm_taper_file(),
#  make_trunk_radius_func(), parse_adqsm_params_file() from
#  tree_geom_utils.py. CBH is NOT part of that function's returned dict
#  (only TrunkVolume/BranchVolume/DBH/etc. are), so this script parses
#  the "CBH:" token directly out of TreesParams.txt itself, the same
#  simple "Key: value" regex idiom parse_adqsm_params_file() already uses
#  internally - not a duplication of that function's own parsing (which
#  covers many other fields this script doesn't need), just the one
#  extra token it doesn't expose.
#
#  Dependencies: numpy, matplotlib (install: pip install numpy matplotlib)
# =====================================================================

import csv
import os
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines   # proxy handle for the one combined "Crown Base Height" legend entry

from tree_geom_utils import parse_adqsm_taper_file, make_trunk_radius_func, parse_adqsm_params_file

# Shared visual style (colors/sizes) - see plot_style.py's own header.
# family_shades() extends FAMILY_GRADIENTS["AdQSM"]'s 4 hand-picked stops
# into as many distinct shades as this chart needs (one per AdQSM variant,
# up to ~21 on B21_S01) - see its own docstring in plot_style.py.
from plot_style import (
    family_shades,
    LEGEND_FONTSIZE,
    AXIS_LABEL_FONTSIZE,
    PANEL_TITLE_FONTSIZE,
    DEFAULT_LINEWIDTH,
)

# =====================  PARAMETERS  ===================================
DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"

# Trees to process in this run - just add a name to extend this to a
# production beech tree once its data folder exists; nothing else needs
# to change.
TREES_TO_RUN = ["B21_S01"]

# Optional: real field-measured DBH per tree, in METERS, for a
# horizontal reference line on the chart - e.g. {"IND07_083": 0.75}.
# Leave a tree out of this dict (or leave the dict empty) to skip the
# reference line for that tree - this is purely a visual aid, it does
# NOT feed into any calibration (unlike adtree_reconstruct_compare.py's
# own FIELD_DBH, which DOES rescale the taper curve - this script never
# rescales anything, it shows each variant's RAW curve so they stay
# comparable to each other).
MEASURED_DBH_M = {"B21_S01": 0.56}

SUMMARY_CSV_PATH = "taper_curve_compare_summary.csv"
# One row per (tree, variant) across every tree in THIS run's
# TREES_TO_RUN - overwritten fresh each run (not an upsert/append
# across sessions like volume_results.csv - this is a lightweight
# diagnostic export, not the master results table).

DEFAULT_LINEWIDTH = 6
# =====================================================================


def discover_variants(tree_dir):
    """Return sorted variant labels (subfolder names) under `tree_dir`
    that contain a taper.txt file - e.g. ["04","05","06","07","08","09",
    "10"]. Sorted numerically if every label is a plain integer string,
    else sorted as plain strings (never crash on a non-numeric variant
    folder name - just fall back to string sort and print a note)."""
    if not os.path.isdir(tree_dir):
        return []
    labels = [
        name for name in sorted(os.listdir(tree_dir))
        if os.path.isfile(os.path.join(tree_dir, name, "taper.txt"))
    ]
    try:
        return sorted(labels, key=int)
    except ValueError:
        print("  NOTE: not every variant folder name under %s is a plain integer - "
              "falling back to plain string sort." % tree_dir)
        return sorted(labels)


def _extract_cbh(params_path):
    """Pull the "CBH: <value>" token directly out of a TreesParams.txt
    file, in METERS - same simple "Key: value" regex idiom
    parse_adqsm_params_file() uses internally, just for the one token
    that function's own returned dict doesn't expose. Returns None if
    the file is missing or has no CBH token (never raises)."""
    if not os.path.exists(params_path):
        return None
    with open(params_path, "r", encoding="latin-1") as f:
        for line in f:
            if "CBH:" not in line:
                continue
            m = re.search(r"CBH:\s*([-\d.eE+]+)", line)
            if m:
                return float(m.group(1))
    return None


def plot_taper_variants(tree_name):
    tree_dir = os.path.join(DATA_ROOT, tree_name)
    variants = discover_variants(tree_dir)
    if not variants:
        print("WARNING: no variant folder with a taper.txt found under %s - skipping %s."
              % (tree_dir, tree_name))
        return

    print("=" * 90)
    print("Tree: %s  (%d variants found: %s)" % (tree_name, len(variants), ", ".join(variants)))
    print("-" * 90)

    fig, ax = plt.subplots(figsize=(10, 7))
    table_rows = []   # (variant, dbh_cm, trunk_vol_or_None, cbh_or_None)
    cbh_drawn = False   # True once at least one variant's CBH line has actually been drawn -
                          # drives whether the single combined legend entry gets added at all below

    # One distinct AdQSM-family shade per variant, instead of leaving color
    # unset (matplotlib's default cycle only has ~10 colors, so it used to
    # repeat once past variant #10) - see family_shades()'s own docstring
    # in plot_style.py for why a continuous colormap is used instead of
    # just cycling FAMILY_GRADIENTS["AdQSM"]'s 4 raw stops.
    variant_colors = family_shades("AdQSM", len(variants))

    for variant_idx, variant in enumerate(variants):
        variant_dir = os.path.join(tree_dir, variant)
        taper_path = os.path.join(variant_dir, "taper.txt")
        params_path = os.path.join(variant_dir, "TreesParams.txt")

        heights, diameters = parse_adqsm_taper_file(taper_path)
        # field_dbh=None: NO rescaling - see this module's own header
        # comment for why (comparing variants' RAW curves is the point).
        trunk_radius_func = make_trunk_radius_func(heights, diameters, field_dbh=None)
        dbh_variant = 2.0 * trunk_radius_func(1.3)

        params = parse_adqsm_params_file(params_path)
        if params is None:
            print("  (no TreesParams.txt reference found at %s - TrunkVolume/CBH shown as 'n/a' "
                  "for this variant)" % params_path)
            trunk_vol = None
        else:
            trunk_vol = params["trunk_vol"]
        cbh = _extract_cbh(params_path)

        table_rows.append((variant, dbh_variant * 100.0, trunk_vol, cbh))

        legend_label = "Variant %s (DBH=%.1fcm, TrunkVolume=%s)" % (
            variant, dbh_variant * 100.0, ("%.2fm3" % trunk_vol) if trunk_vol is not None else "n/a")
        line, = ax.plot(heights, diameters, marker="o", markersize=3, linewidth=DEFAULT_LINEWIDTH,
                         color=variant_colors[variant_idx], label=legend_label)

        # CBH: thin dashed vertical line in THIS variant's own colour, low
        # alpha - the diagnostic that revealed the crown-base connection
        # for IND07_083 (see TODO_investigations.md) - kept visible but
        # unobtrusive, one per variant that actually has a CBH value.
        # Deliberately NOT given a `label=` here (would add one near-
        # duplicate legend row per variant, up to 21 on B21_S01) - a
        # single combined legend entry is added once, below the loop,
        # via a neutral proxy handle instead (see cbh_drawn/mlines.Line2D
        # below) - kept color-neutral rather than reusing this variant's
        # own colour, since a first-variant-coloured swatch would wrongly
        # suggest CBH is tied to that one variant's colour specifically,
        # when in fact every variant draws its own differently-coloured
        # CBH line.
        if cbh is not None:
            ax.axvline(cbh, color=line.get_color(), linestyle="--", linewidth=DEFAULT_LINEWIDTH, alpha=0.4)
            cbh_drawn = True

    ax.axvline(1.3, color="gray", linestyle="--", linewidth=DEFAULT_LINEWIDTH, label="DBH height")

    if tree_name in MEASURED_DBH_M:
        ax.axhline(MEASURED_DBH_M[tree_name], color="red", linestyle=":", linewidth=DEFAULT_LINEWIDTH,
                   label="Measured DBH (field)")

    ax.set_xlabel("Height [m]", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Diameter [m]", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title("%s: AdQSM taper.txt curves across %d variants" % (tree_name, len(variants)),
                 fontsize=PANEL_TITLE_FONTSIZE)

    # One combined legend entry for every per-variant CBH line drawn above
    # (cbh_drawn), instead of one per variant - a neutral gray dashed proxy
    # handle (matching "DBH height"'s own neutral style just above), NOT
    # any one variant's own colour, so the swatch doesn't wrongly look
    # tied to a single variant. Only added when at least one CBH line was
    # actually drawn (cbh_drawn), so a tree with no CBH data at all gets
    # no unused legend row. Every other legend entry (each variant's own
    # curve, "DBH height", "Measured DBH (field)") comes from
    # get_legend_handles_labels() unchanged, in its normal order.
    handles, labels = ax.get_legend_handles_labels()
    if cbh_drawn:
        cbh_proxy = mlines.Line2D([], [], color="gray", linestyle="--", linewidth=DEFAULT_LINEWIDTH, alpha=0.4,
                                   label="Crown Base Height (CBH)")
        handles.append(cbh_proxy)
        labels.append(cbh_proxy.get_label())
    ax.legend(handles=handles, labels=labels, fontsize=LEGEND_FONTSIZE, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()

    out_dir = os.path.join("plots", tree_name)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "taper_curve_compare_%s.png" % tree_name)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)

    # Console table - % difference vs. MEASURED_DBH_M for this tree if
    # given, else vs. the SMALLEST-numbered variant present (mirroring
    # resolve_reference_method_none()'s "smallest variant present"
    # convention used elsewhere in this project, for consistency).
    if tree_name in MEASURED_DBH_M:
        ref_dbh_cm = MEASURED_DBH_M[tree_name] * 100.0
        ref_label = "measured DBH (%.1fcm)" % ref_dbh_cm
    else:
        ref_dbh_cm = table_rows[0][1]   # variants list is already sorted -> smallest-numbered first
        ref_label = "smallest variant present (%s, %.1fcm)" % (table_rows[0][0], ref_dbh_cm)

    print()
    print("%-10s %-10s %-15s %-12s" % ("variant", "DBH_cm", "TrunkVolume_m3", "%%diff_vs_%s" % ref_label))
    # result_rows: one dict per variant, for the combined multi-tree summary
    # CSV written once by the RUN section below - built alongside the
    # console printout using the exact SAME pct_diff/ref_label already
    # computed for it (not recomputed a second time), so the two can never
    # drift apart. Purely additive - every print/plot/save above/below is
    # unchanged.
    result_rows = []
    for variant, dbh_cm, trunk_vol, cbh in table_rows:
        pct_diff = 100.0 * (dbh_cm - ref_dbh_cm) / ref_dbh_cm if ref_dbh_cm else float("nan")
        print("%-10s %-10.2f %-15s %+.2f%%" % (
            variant, dbh_cm, ("%.4f" % trunk_vol) if trunk_vol is not None else "n/a", pct_diff))
        result_rows.append({
            "tree": tree_name,
            "variant": variant,
            "dbh_cm": dbh_cm,
            # "" (not the console table's "n/a") for a missing value - CSVs
            # should stay machine-readable, "n/a" is a display-only choice.
            "trunk_vol_m3": trunk_vol if trunk_vol is not None else "",
            "cbh_m": cbh if cbh is not None else "",
            "pct_diff_vs_reference": pct_diff,
            "reference_label": ref_label,
        })
    print()

    return result_rows


# =========================  RUN  =====================================
if __name__ == "__main__":
    all_rows = []
    n_trees_with_rows = 0   # trees that actually contributed rows, NOT len(TREES_TO_RUN) -
                              # a skipped tree (no variant folders found, plot_taper_variants()
                              # returns None) must not be counted as one of the "trees" below
    for tree in TREES_TO_RUN:
        rows = plot_taper_variants(tree)
        if rows:
            all_rows.extend(rows)
            n_trees_with_rows += 1

    if all_rows:
        with open(SUMMARY_CSV_PATH, "w", encoding="utf-8", newline="") as f:
            fieldnames = ["tree", "variant", "dbh_cm", "trunk_vol_m3", "cbh_m",
                          "pct_diff_vs_reference", "reference_label"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        print("Saved summary CSV: %s (%d rows across %d trees)" %
              (SUMMARY_CSV_PATH, len(all_rows), n_trees_with_rows))
    else:
        print("No rows to write - every tree in TREES_TO_RUN was skipped (see WARNINGs above).")
