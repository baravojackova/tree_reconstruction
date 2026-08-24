# -*- coding: utf-8 -*-
# =====================================================================
#  STANDALONE diagnostic: how much of AdQSM's branch LENGTH (and, via the
#  reliable cylinder-volume approximation, branch VOLUME) sits in very thin
#  branches, binned by diameter.
# ---------------------------------------------------------------------
#  WHY THIS SCRIPT EXISTS: AdQSM's branch_len_m in volume_results.csv is
#  much larger than AdTree's/TreeQSM's. One suspected reason is that AdQSM
#  reconstructs a huge number of very thin twigs that AdTree/TreeQSM simply
#  cannot resolve from the point cloud (lidar resolution limit). This script
#  checks that hypothesis directly from AdQSM's own BranchStructure.txt: it
#  bins every branch by diameter and reports what share of the total LENGTH
#  and total VOLUME each diameter class holds. If a small share of volume
#  sits in a huge share of length, that supports the "many thin twigs"
#  explanation.
#
#  This script does NOT touch volume_results.csv or any other shared
#  pipeline file - it only reads BranchStructure.txt/TreesParams.txt and
#  prints/plots a report. Safe to run any time, on its own.
#
#  Volume is computed with pi*r^2*length per branch (the same reliable
#  cylinder approximation used elsewhere in this project, e.g.
#  report_thin_branch_volume() in tree_geom_utils.py) - NOT the file's own
#  "volume(L)" column, which adtree_reconstruct_compare.py's
#  report_adqsm_thin_branch() already found to be off by orders of
#  magnitude from AdQSM's own official totals and therefore untrustworthy.
# =====================================================================

import os

import numpy as np
import matplotlib
matplotlib.use("Agg")   # draw without opening a window (this script only saves a PNG)
import matplotlib.pyplot as plt

# tree_geom_utils.py already has the exact BranchStructure.txt header/column
# helpers used by adtree_reconstruct_compare.py - reuse them here so this
# diagnostic can never disagree with the main pipeline about column meaning.
from tree_geom_utils import (
    _read_adqsm_branch_header, _find_adqsm_column,
    print_adqsm_branch_file_sample, parse_adqsm_params_file,
)

# =====================  PARAMETERS  ===================================
# Same tree/AdQSM-variant this project is currently focused on (see
# TREE_NAME/AdQSM_DIR in adtree_reconstruct_compare.py). Change these two
# lines to point at a different tree/variant - this is a one-off diagnostic,
# not part of the multi-variant pipeline, so the path is just hard-coded
# here instead of rebuilding ADQSM_VARIANTS machinery.
DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"
TREE_NAME = "IND01_054"
ADQSM_DIR = os.path.join(DATA_ROOT, TREE_NAME, "05")

# The two AdQSM export files this script reads.
BRANCH_FILE = os.path.join(ADQSM_DIR, "BranchStructure.txt")   # per-branch table (diameter, length, ...)
PARAMS_FILE = os.path.join(ADQSM_DIR, "TreesParams.txt")       # whole-tree totals, used for the DBH unit check

# Manual column override: leave as None to auto-detect (by name, case-
# insensitive substring match on "diam"/"length"); set to a column INDEX
# (e.g. 2) if auto-detection fails or picks the wrong column - the header
# printed in STEP 1 below shows every column's index so you can pick one.
DIAMETER_COL = None
LENGTH_COL = None

# Diameter bin edges in CENTIMETRES. Every branch is put into the bin
# [edges[i], edges[i+1]) it falls into; the last edge (1000) is just a large
# number so the final bin ("20+") catches everything above 20 cm.
DIAMETER_BINS_CM = [0, 1, 2, 5, 10, 20, 1000]

# Where the summary bar chart is saved - same "plots" folder plot_volumes.py
# already uses, so all diagnostic PNGs for this project live in one place.
PLOTS_DIR = "plots"
PLOT_FILENAME = "adqsm_thin_branch_length_share.png"
# =====================================================================


# =========================  STEP 1  ===================================
# Print the raw header + first 10 data rows (with column index) BEFORE
# using the file, so the diameter/length columns can be visually confirmed
# instead of trusted blindly. Reuses the exact same diagnostic function
# adtree_reconstruct_compare.py already uses for this (PRINT_ADQSM_BRANCH_SAMPLE).
print("=" * 70)
print("STEP 1: raw BranchStructure.txt sample (%s)" % BRANCH_FILE)
print("=" * 70)
print_adqsm_branch_file_sample(BRANCH_FILE, n=10)


# =========================  STEP 2  ===================================
# Parse the whole table into plain Python lists, one list per column, using
# the same TAB-split / "first field is an integer order" rule the existing
# AdQSM parsers in tree_geom_utils.py already use, so this reads the file
# exactly the way the rest of the project does.
print("\n" + "=" * 70)
print("STEP 2: parsing BranchStructure.txt")
print("=" * 70)

header_cols = _read_adqsm_branch_header(BRANCH_FILE)   # column names, or None if no header row
if not header_cols:
    raise SystemExit("No header row found in %s - cannot auto-detect columns. "
                      "Set DIAMETER_COL/LENGTH_COL manually and adapt this script." % BRANCH_FILE)
print("Header columns found: %s" % ", ".join("[%d] %s" % (i, n) for i, n in enumerate(header_cols)))

# rows: one list of raw string fields per valid data row (same "len(parts)
# >= 7 and first field is an int" validity rule as parse_adqsm_branch_file()
# in tree_geom_utils.py, so this keeps exactly the rows the main pipeline
# would also use).
rows = []
with open(BRANCH_FILE, "r", encoding="latin-1") as f:
    for line in f:
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 7:
            continue
        try:
            int(parts[0])   # data rows start with an integer branch order
        except ValueError:
            continue
        rows.append(parts)
print("Parsed %d branch rows." % len(rows))

# Auto-detect the diameter/length columns by name, unless the user overrode
# them above. _find_adqsm_column() does a case-insensitive substring match
# and returns None if there's no match at all.
diam_idx = DIAMETER_COL if DIAMETER_COL is not None else _find_adqsm_column(header_cols, ["diam"])
length_idx = LENGTH_COL if LENGTH_COL is not None else _find_adqsm_column(header_cols, ["length"])

# Count how many columns actually contain "diam"/"length" in their name, so
# an ambiguous match (0 or >1 candidates) is caught instead of silently
# picking a possibly-wrong column.
diam_candidates = [i for i, n in enumerate(header_cols) if "diam" in n.lower()]
length_candidates = [i for i, n in enumerate(header_cols) if "length" in n.lower()]

if DIAMETER_COL is None and len(diam_candidates) != 1:
    raise SystemExit("Could not uniquely auto-detect the diameter column (found %d candidates: %s). "
                      "Available columns: %s. Set DIAMETER_COL manually at the top of this script."
                      % (len(diam_candidates), diam_candidates, header_cols))
if LENGTH_COL is None and len(length_candidates) != 1:
    raise SystemExit("Could not uniquely auto-detect the length column (found %d candidates: %s). "
                      "Available columns: %s. Set LENGTH_COL manually at the top of this script."
                      % (len(length_candidates), length_candidates, header_cols))

print("Using diameter column: [%d] %s" % (diam_idx, header_cols[diam_idx]))
print("Using length column  : [%d] %s" % (length_idx, header_cols[length_idx]))

# Pull the two chosen columns out of `rows` into numpy float arrays.
diam_raw = np.array([float(r[diam_idx]) for r in rows])     # raw diameter values, UNIT NOT YET CONFIRMED
length_m = np.array([float(r[length_idx]) for r in rows])   # branch length, already in metres (header says "length(m)")


# =========================  STEP 3  ===================================
# Confirm whether the diameter column is in cm, m, or mm by comparing its
# max value against the tree's known DBH (from TreesParams.txt).
print("\n" + "=" * 70)
print("STEP 3: diameter column unit check")
print("=" * 70)
print("Raw diameter column: min=%.4f  max=%.4f  mean=%.4f" % (diam_raw.min(), diam_raw.max(), diam_raw.mean()))

params = parse_adqsm_params_file(PARAMS_FILE)
if params is None or params.get("dbh") is None:
    print("WARNING: could not read DBH from %s - cannot confirm units. "
          "Continuing and treating the diameter column as CENTIMETRES (the visually "
          "most likely case) - CHECK THE PRINTED MIN/MAX/MEAN ABOVE YOURSELF." % PARAMS_FILE)
    diam_to_m = lambda d: d / 100.0   # assume cm -> m
else:
    dbh_m = params["dbh"]   # DBH in METRES (parse_adqsm_params_file already converts TreeDBH cm -> m)
    print("Reference DBH (TreesParams.txt): %.4f m" % dbh_m)
    ratio = diam_raw.max() / dbh_m   # how many "raw units" fit in one metre of DBH
    # ratio ~100  -> raw column is in CENTIMETRES (1 m = 100 cm)
    # ratio ~1    -> raw column is already in METRES
    # ratio ~1000 -> raw column is in MILLIMETRES (1 m = 1000 mm)
    if 30 <= ratio <= 300:
        print("-> max diameter is ~%.0fx the DBH: column looks like CENTIMETRES. Using cm -> m (/100)." % ratio)
        diam_to_m = lambda d: d / 100.0
    elif 0.3 <= ratio <= 3:
        print("!" * 70)
        print("WARNING: max diameter is only ~%.1fx the DBH (in metres) - the diameter" % ratio)
        print("column looks like it is ALREADY IN METRES, not centimetres as usually")
        print("assumed for this kind of AdQSM export. Using the column AS-IS (no /100),")
        print("but please double-check the printed min/max/mean above yourself.")
        print("!" * 70)
        diam_to_m = lambda d: d
    elif 300 <= ratio <= 3000:
        print("!" * 70)
        print("WARNING: max diameter is ~%.0fx the DBH - the diameter column looks like" % ratio)
        print("it is in MILLIMETRES, not centimetres. Using mm -> m (/1000), but please")
        print("double-check the printed min/max/mean above yourself.")
        print("!" * 70)
        diam_to_m = lambda d: d / 1000.0
    else:
        print("!" * 70)
        print("WARNING: diameter/DBH ratio (%.2f) does not clearly match cm, m, or mm. "
              "Continuing and treating the diameter column as CENTIMETRES (the visually "
              "most likely case) - CHECK THE PRINTED MIN/MAX/MEAN ABOVE YOURSELF." % ratio)
        print("!" * 70)
        diam_to_m = lambda d: d / 100.0

diam_m = diam_to_m(diam_raw)   # diameter column converted to metres, ready for volume math
diam_cm = diam_m * 100.0       # ... and to centimetres, for the binning in STEP 4


# =========================  STEP 4  ===================================
# Bin every branch by its (now confirmed-unit) diameter and sum length and
# cylinder-approximation volume per bin.
print("\n" + "=" * 70)
print("STEP 4: length/volume share by diameter bin")
print("=" * 70)

volume_m3 = np.pi * (diam_m / 2.0) ** 2 * length_m   # per-branch cylinder volume approx. (pi*r^2*length)
total_length = float(length_m.sum())
total_volume = float(volume_m3.sum())

bin_labels = []      # human-readable "lo-hi" range per bin, for the printed/plotted table
n_branches_list = []
length_sum_list = []
length_pct_list = []
volume_sum_list = []
volume_pct_list = []

for i in range(len(DIAMETER_BINS_CM) - 1):
    lo, hi = DIAMETER_BINS_CM[i], DIAMETER_BINS_CM[i + 1]
    # last bin (hi == 1000, our "infinity" placeholder) is labelled "lo+" instead of "lo-1000"
    label = ("%g+" % lo) if hi >= 1000 else ("%g-%g" % (lo, hi))
    mask = (diam_cm >= lo) & (diam_cm < hi)

    n_branches = int(mask.sum())
    length_sum = float(length_m[mask].sum())
    volume_sum = float(volume_m3[mask].sum())
    length_pct = (length_sum / total_length * 100.0) if total_length else 0.0
    volume_pct = (volume_sum / total_volume * 100.0) if total_volume else 0.0

    bin_labels.append(label)
    n_branches_list.append(n_branches)
    length_sum_list.append(length_sum)
    length_pct_list.append(length_pct)
    volume_sum_list.append(volume_sum)
    volume_pct_list.append(volume_pct)

# Print one row per bin, columns as requested.
print("%-12s %-12s %-14s %-12s %-16s %-12s"
      % ("diam [cm]", "n_branches", "length_sum[m]", "length_pct[%]", "volume_sum[m3]", "volume_pct[%]"))
for i in range(len(bin_labels)):
    print("%-12s %-12d %-14.2f %-12.1f %-16.5f %-12.1f"
          % (bin_labels[i], n_branches_list[i], length_sum_list[i], length_pct_list[i],
             volume_sum_list[i], volume_pct_list[i]))
print("%-12s %-12d %-14.2f %-12.1f %-16.5f %-12.1f"
      % ("TOTAL", sum(n_branches_list), total_length, 100.0, total_volume, 100.0))


# =========================  STEP 5  ===================================
# Bar chart: two bars per diameter bin (length_pct vs volume_pct), so the
# "thin twigs = lots of length, little volume" mismatch is visible at a glance.
print("\n" + "=" * 70)
print("STEP 5: plotting")
print("=" * 70)

if not os.path.isdir(PLOTS_DIR):
    os.makedirs(PLOTS_DIR)   # create the shared plots folder if this is the first script to run

x = np.arange(len(bin_labels))   # one x position per diameter bin
width = 0.35                     # width of each bar, so two bars fit side by side per bin

fig, ax = plt.subplots(figsize=(9, 5))
ax.bar(x - width / 2, length_pct_list, width, label="length share [%]")
ax.bar(x + width / 2, volume_pct_list, width, label="volume share [%]")
ax.set_xticks(x)
ax.set_xticklabels(bin_labels)
ax.set_xlabel("diameter bin [cm]")
ax.set_ylabel("share of total [%]")
ax.set_title("AdQSM %s: branch length vs. volume share by diameter" % TREE_NAME)
ax.legend()
fig.tight_layout()

plot_path = os.path.join(PLOTS_DIR, PLOT_FILENAME)
fig.savefig(plot_path, dpi=150)
plt.close(fig)   # free the figure's memory now that it's saved
print("Saved: %s" % plot_path)
