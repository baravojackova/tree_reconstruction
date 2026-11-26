"""
plot_stiffness_comparison.py
------------------------------
Compares experimental pulling-test stiffness (grey boxplots, per tree per
year) against FEM model stiffness results (colored point overlay, per tree
per scenario).

Produces eight figures, two variants of each of two quantities, in a full
and a simplified form (built by the shared build_stiffness_figure()):
  - stiffness_comparison_rotation_all.png / _trimmed.png                  (full: cam/nocam/combined)
  - stiffness_comparison_rotation_combined_all.png / _trimmed.png         (VARIANT 1: combined only)
  - stiffness_comparison_displacement_all.png / _trimmed.png              (full: Pt3/Pt4/combined)
  - stiffness_comparison_displacement_combined_all.png / _trimmed.png     (VARIANT 1: combined only)

The full figures have 3 side-by-side "stiffness" subplots (two "single
sensor" views and one "combined" view - median of the two sensors, pull by
pull); the "_combined" ones keep only the "combined" panel, for a less
dense write-up figure while the per-sensor comparison is still needed
elsewhere. Both share the same tree x-axis, the same FEM point overlay,
and the same PLUS grouped-bar deviation subplot underneath (relative
deviation %, FEM vs experimental median, one bar per model - see
draw_deviation_bars()) - built from ONE deviation_table per quantity, so a
full and its simplified sibling can never show different numbers.  The
stiffness subplots and stiffness_deviation.csv are IDENTICAL between the
"_all" and "_trimmed" files - the two variants differ ONLY in which
models' bars are drawn in the deviation subplot, per
DEVIATION_EXCLUDE_MODELS (its own comment explains why: a known-bad
parameter setting on the treeqsm_pd002 runs produces deviations extreme
enough - down to -2240% - to make the "_all" deviation subplot unreadable
for every other model, hence "_trimmed"; "_all" still exists because how
the fault manifests is itself worth seeing).

MODEL_INCLUDE (None by default) restricts every figure above to a chosen
subset of models - a FILTER ON WHAT IS DRAWN ONLY (stiffness_deviation.csv
always keeps every row for every model). When set, every filename gets
MODEL_INCLUDE_SUFFIX appended (e.g. "..._all_noacel.png") so a filtered run
can never overwrite the full one - see resolve_model_filter() and
MODEL_INCLUDE's own comment.

How to use it:
  - Edit the values in the "CONFIGURATION" section below.
  - Run the script from a terminal in VS Code:  python plot_stiffness_comparison.py
  - Check the eight saved PNG files and the console summary.
"""

import re
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pathlib import Path
from datetime import datetime
import sys

# This script lives in scripts/python/ansys_eval, but paths.py lives one
# folder up (scripts/python). Python only looks for imports next to the
# script, so we first add the folder above to its search list.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folder this script lives in. Spelled out in full on purpose: a bare
# filename is resolved against the folder the TERMINAL is in, not the
# folder this script lives in.
WORK_FOLDER = str(paths.ANSYS_EVAL_DIR)

# Input data and output files live in their own subfolders of WORK_FOLDER.
DATA_FOLDER = Path(WORK_FOLDER) / "data"
OUT_FOLDER = Path(WORK_FOLDER) / "out"

EXPERIMENTAL_CSV_PATH = DATA_FOLDER / "dynatree_data_merge_0926.csv"
FEM_CSV_PATH = OUT_FOLDER / "slopes_combined.csv"          # Stage 2 output

DAYS_TO_PROCESS = ["2021-03-22", "2022-04-05"]      # list of measurement days to include; each day's year becomes one boxplot group
MEASUREMENT = "M01"
TYPE_FILTER = "normal"
REDUCTION_LEVEL_FILTER = "none"    # None = no filter, a string = exact match, a list = match any value in the list (same filter convention as filter_and_plot.py)
TREE_FILTER = None                 # None = all trees, otherwise a tree ID or list of tree IDs

# The experimental dataset labels trees "BK07"; the reconstruction project
# labels the same tree "B21_S07". The number is identical in both, so the
# conversion is purely cosmetic - but without it, no experimental row ever
# matches a FEM row and every plot comes out empty.
TREE_ID_PREFIX = "B21_S"

# True  = only draw trees that have at least one FEM result, so the figure
#         is not mostly empty while the run plan is incomplete.
# False = draw every tree present in the experimental data.
ONLY_TREES_WITH_FEM = True

# Set once by main() (after TREES is resolved against the actual data),
# read by with_filter_suffix() below. Stays "" - meaning "no suffix" -
# whenever TREES is "all", which is what keeps every filename unchanged by
# this particular filter (GRAVITY_VARIANT/MODEL_INCLUDE add their own
# suffixes independently - see with_filter_suffix()).
TREES_SUFFIX = ""

# Set once by main() before any figure is built, read by
# add_provenance_footer() via build_stiffness_figure(). See
# check_input_freshness() for how PROVENANCE_NEWEST_INPUT_MTIME is computed.
PROVENANCE_NEWEST_INPUT_MTIME = None
PROVENANCE_FILTERS_TEXT = ""


def gravity_variant_filename_suffix(variant: str) -> str:
    """
    Filename fragment for GRAVITY_VARIANT. "" for "with_gravity" ONLY - the
    81 unsuffixed PNGs/CSVs already on disk were all written before
    GRAVITY_VARIANT existed, under what was effectively a gravity-on
    default, so "with_gravity" is the one setting that keeps filenames
    unchanged. "no_acel" (the new default) and "both" both get an explicit
    suffix, so a gravity-off figure can never silently overwrite one of
    those gravity-on files - see GRAVITY_VARIANT's own comment.
    """
    return {"no_acel": "gravity-noacel", "both": "gravity-both"}.get(variant, "")


def with_filter_suffix(path) -> Path:
    """
    Append every ACTIVE **drawing** filter's suffix to path, before its
    extension, in a fixed order (gravity variant, then model filter, then
    tree filter) - e.g. "stiffness_comparison_rotation_all.png" ->
    "..._gravity-noacel_noacel_trees-07-12.png". Used for FIGURE
    (PNG) paths ONLY - see csv_filter_suffix() for stiffness_deviation.csv,
    whose content GRAVITY_VARIANT/MODEL_INCLUDE never narrow (draw-only
    filters), so neither may appear in ITS name. Returns path unchanged
    only when GRAVITY_VARIANT="with_gravity" AND nothing else is filtered -
    see gravity_variant_filename_suffix()'s own comment for why
    "with_gravity" (not "all=nothing filtered") is the unsuffixed case now.
    """
    parts = []
    gravity_suffix = gravity_variant_filename_suffix(GRAVITY_VARIANT)
    if gravity_suffix:
        parts.append(gravity_suffix)
    if MODEL_INCLUDE is not None:
        parts.append(MODEL_INCLUDE_SUFFIX)
    if TREES_SUFFIX:
        parts.append(TREES_SUFFIX)
    if not parts:
        return Path(path)
    root, ext = str(path).rsplit(".", 1)
    return Path(f"{root}_{'_'.join(parts)}.{ext}")


def csv_filter_suffix(path) -> Path:
    """
    Append ONLY the TREES suffix to path - used for stiffness_deviation.csv.
    Unlike the PNGs (with_filter_suffix()), stiffness_deviation.csv's ROW
    CONTENT is filtered by TREES but NEVER by MODEL_INCLUDE or
    GRAVITY_VARIANT (both are draw-only
    filters - see fem_df vs. fem_df_drawing in main()), so neither may
    appear in this file's name. This is what lets
    method_deviation_overview.py always read ONE canonical
    stiffness_deviation.csv regardless of what THIS script's own
    MODEL_INCLUDE/GRAVITY_VARIANT happen to be set to - an earlier version
    of Stage 7 briefly broke that guarantee by reusing with_filter_suffix()
    here too (see GRAVITY_VARIANT's own comment for the postmortem).
    """
    if not TREES_SUFFIX:
        return Path(path)
    root, ext = str(path).rsplit(".", 1)
    return Path(f"{root}_{TREES_SUFFIX}.{ext}")


# Each quantity ("rotation", "displacement") is written as TWO figures -
# see stiffness_plot_path() and DEVIATION_EXCLUDE_MODELS below for why.
def stiffness_plot_path(quantity: str, variant: str) -> Path:
    """
    quantity: "rotation" or "displacement". variant: "all" or "trimmed".
    -> .../out/stiffness_comparison_<quantity>_<variant>[_<suffix>].png
    The suffix is built by with_filter_suffix() from whichever of
    GRAVITY_VARIANT/MODEL_INCLUDE/TREES are active - see that function's comment.
    """
    return with_filter_suffix(OUT_FOLDER / f"stiffness_comparison_{quantity}_{variant}.png")


# Table of relative deviations (FEM vs experimental median), written as a
# CSV next to the figures and printed to the console.
DEVIATION_CSV_PATH = OUT_FOLDER / "stiffness_deviation.csv"

# Models excluded from the DEVIATION SUBPLOT ONLY in the "_trimmed" figure
# (matched against the "model" column, not the scenario number, so this
# list stays readable and does not break when scenario numbers are
# reassigned). The stiffness panels above and stiffness_deviation.csv are
# UNAFFECTED - they always show/contain every model.
#
# WHY: the treeqsm_pd002 runs were computed with a parameter setting known
# to be wrong (being fixed elsewhere - not diagnosed here). Its deviations
# run down to -2240%, which stretches the axis so far that every other
# model becomes unreadable. The "_all" figure still shows it (how the
# fault manifests is itself worth seeing); the "_trimmed" figure hides it
# from the subplot's y-axis scaling while keeping its legend entry, so it
# reads as "deliberately excluded", not "missing".
DEVIATION_EXCLUDE_MODELS = ["treeqsm_pd002"]

# None = every model is DRAWN. A list = only these models are drawn -
# matched against the "model" column (not the scenario number, so this
# stays readable if scenario numbers are ever reassigned), tolerant of the
# space/underscore inconsistency in the no-acel suffix (see
# normalize_model_name_for_matching()) but never rewriting a displayed
# name - that inconsistency is being fixed in ANSYS, not here.
#
# This is a FILTER ON WHAT IS DRAWN ONLY, independent of
# DEVIATION_EXCLUDE_MODELS above: fem_df and stiffness_deviation.csv always
# keep every row for every model, regardless of this setting.
#
# NOTE (Stage 7): this had been set to the two no_acel model names as a
# manual stand-in for "only show the gravity-off reference variant" - that
# is now GRAVITY_VARIANT's job (see below), applied structurally instead of
# by remembering to set this correctly. Reset to None so MODEL_INCLUDE goes
# back to its own job: narrowing by RECONSTRUCTION METHOD (e.g. AdTree
# only), not by gravity variant. Set it again only for that purpose.
MODEL_INCLUDE = None

# Appended to every figure filename whenever MODEL_INCLUDE is set (e.g.
# stiffness_comparison_rotation_all_noacel.png), so a filtered figure can
# never be confused with, or silently overwrite, the full one. main()
# refuses to run if MODEL_INCLUDE is set and this is empty.
MODEL_INCLUDE_SUFFIX = "no_acel"

# =========================================================================
# STALE-OUTPUT / PROVENANCE SAFEGUARDS (added Sep 2026)
# -------------------------------------------------------------------------
# WHY these exist: a filtered run (MODEL_INCLUDE set to the two "no acel"
# models) was left running unattended, its output PNGs sat in out/ next to
# the unfiltered ones, MODEL_INCLUDE was later reset to None, and every
# later re-run only refreshed the UNFILTERED files (writing a "_noacel"
# file is a no-op when MODEL_INCLUDE is None - see with_filter_suffix()).
# The stale "_noacel" files were then mistaken for current output because
# nothing on the image said when or how it was made. These three settings
# make that mistake impossible to repeat silently.
# =========================================================================

# True = stamp a small grey footer on every figure with: this script's
# filename, when the figure was generated, the newest modification time
# among the raw ANSYS files it depends on, the active GRAVITY_VARIANT/TREES/
# MODEL_INCLUDE filters, and how many trees/models were actually drawn -
# see add_provenance_footer(). Set to False only for a publication figure
# where the footer must not appear.
FIGURE_PROVENANCE_FOOTER = True

# True = before plotting, refuse to run if any raw ANSYS file this script
# ultimately depends on (out_static_bk*_scenario*.txt) is NEWER than
# FEM_CSV_PATH - i.e. compute_slopes.py/load_static_scenarios.py have not
# been re-run since new ANSYS results landed, so FEM_CSV_PATH is stale.
# False = print the same warning but keep going (useful when you know why
# it is stale and want the figures anyway). See check_input_freshness().
ABORT_ON_STALE_INPUT = True

# Folder holding the raw ANSYS output files, needed ONLY for the freshness
# check above (this script otherwise never reads raw ANSYS files directly -
# it reads FEM_CSV_PATH, which load_static_scenarios.py/compute_slopes.py
# already built from them). Same folder every other ansys_eval script reads.
ANSYS_RAW_FOLDER = str(paths.ANSYS_DIR)
ANSYS_RAW_PATTERN_FOR_FRESHNESS = "out_static_bk*_scenario*.txt"

# =========================================================================
# THE GRAVITY-VARIANT STRUCTURAL LOCK (added Stage 7)
# -------------------------------------------------------------------------
# WHY this exists: before this, "only show the gravity-off reference
# variant" was a CONVENTION - something to remember to set correctly via
# MODEL_INCLUDE, SCENARIOS, or DEVIATION_PANEL_VARIANT depending on which
# script you were looking at. Those disagreed with each other at least once
# (method_deviation_overview.py's MODEL_INCLUDE, set to the two no_acel
# names, intersected with base_models - "every NON-no_acel model" - to
# nothing, silently starving every panel).
#
# GRAVITY_VARIANT is applied ONCE, in main(), to fem_df_drawing - the
# dataframe every FIGURE reads - never to fem_df itself, which stays the
# complete, both-variants table that feeds stiffness_deviation.csv (see
# fem_df's own comment in main()). This mirrors MODEL_INCLUDE's existing
# guarantee exactly, and for the same reason: this script's CSV is a shared
# pipeline artifact (method_deviation_overview.py reads it), and it must
# stay ONE canonical, config-independent file - an earlier version of this
# stage's own work briefly filtered fem_df itself, which silently broke
# method_deviation_overview.py's self-weight panel (it needs BOTH gravity
# halves of every pair, and the "no_acel"-filtered CSV no longer had the
# gravity-on half at all).
# =========================================================================

# "no_acel"      = only the gravity-off runs (scenarios 3/4) - the reference
#                  variant, and the default.
# "with_gravity" = only the gravity-on runs (scenarios 1/2).
# "both"         = no gravity filtering at all (old, pre-Stage-7 behaviour).
GRAVITY_VARIANT = "no_acel"

# --- WHICH TREES TO PLOT -------------------------------------------------
# "all"  = every tree found in FEM_CSV_PATH (this is the old behaviour)
# or a list, mixed forms allowed:  [7, 12]  /  ["B21_S07", "S12"]  /  ["7"]
# NOTE: these are TREE NUMBERS, not the bk<N> index from the ANSYS filenames.
#       bk<N> is a ROW POSITION in the description file (bk3 = tree 7) - see
#       load_ansys_order_to_tree_id()'s docstring. normalize_tree_selector()
#       refuses a bare "bk3"-style entry outright rather than guessing.
TREES = "all"

# If True, the script raises an error (after printing the selection report)
# when a tree listed above has no data at all. Deliberate: a silently
# empty/narrower plot is what caused the Sep 2026 incident this option was
# added to catch. Set to False only if you knowingly want to ask for trees
# that may not exist yet.
#
# NOTE (Stage 8): this used to also gate a SCENARIOS constant. SCENARIOS
# was retired - it was fully redundant with GRAVITY_VARIANT (model family)
# x MODEL_INCLUDE (reconstruction method): a scenario number is just that
# pair, spelled as a number instead of by name - see Stage 7 Part A. Use
# GRAVITY_VARIANT for which gravity variant to draw.
STRICT_SELECTION = True

# How much horizontal space (out of 1.0 available per tree) the
# side-by-side boxes/points occupy within one tree's group.
GROUP_WIDTH = 0.6

# Fixed pastel colors for FEM models, assigned in alphabetical order of
# the model name. If there are more than 8 models, this list is cycled
# back to the start.
PASTEL_SCENARIO_COLORS = [
    "#94BBEA",  # pastel blue
    "#F5B399",  # pastel orange
    "#8DD7BC",  # pastel aqua
    "#F6D07F",  # pastel yellow
    "#F3BDD1",  # pastel magenta/pink
    "#7FC17F",  # pastel green
    "#A49CD3",  # pastel violet
    "#F1A4A3",  # pastel red
]

# --- Style config -----------------------------------------------------------
TITLE_FONTSIZE = 14
AXIS_LABEL_FONTSIZE = 12
TICK_LABEL_FONTSIZE = 7
LEGEND_FONTSIZE = 12
LEGEND_TITLE_FONTSIZE = 10
LINE_WIDTH = 1.5      # not really used here, but keep for consistency
MARKER_SIZE = 8        # FEM scenario point size
BOX_EDGE_LINEWIDTH = 1.0

# Legends, moved outside the axes. bbox_to_anchor is in axes fraction
# coordinates: LEGEND_BBOX_X > 1 sits just right of the axes edge,
# LEGEND_Y_TOP anchors the topmost stacked legend, and each further legend
# on the same axes is offset down by LEGEND_Y_STEP so they stack instead of
# overlapping. RIGHT_MARGIN is how much of the figure width (0-1) is left
# for the plot itself - the rest is blank margin reserved for the legends.
LEGEND_BBOX_X = 1.02
LEGEND_Y_TOP = 1.0
LEGEND_Y_STEP = 0.28
RIGHT_MARGIN = 0.80

# Deviation panel (replaces the old multi-line text footnote, which had
# grown to ~24 lines and ran over the axes): a grouped-bar subplot below
# the 3 stiffness panels, sharing their tree x-axis. One bar per model per
# FEM quantity, dodged with GROUP_WIDTH (reused from the boxplot dodging
# above) across models, plus a sub-offset for the 1-2 FEM columns within
# each model's slot (see draw_deviation_bars()).
DEVIATION_HEIGHT_RATIO = 0.4    # deviation panel height, as a fraction of the stiffness panels' height
DEVIATION_BAR_WIDTH = 0.09      # width of ONE FEM-column's bar within a model's dodge slot
DEVIATION_BASE_HATCH = None     # slope_rotation_base / slope_displacement: solid
DEVIATION_BASE2_HATCH = "///"   # slope_rotation_base2: hatched
DEVIATION_BAR_EDGECOLOR = "black"   # bar edge - hatch patterns need a visible edge colour to render
DEVIATION_BAR_EDGE_LINEWIDTH = 0.8
DEVIATION_ZERO_LINE_COLOR = "0.3"   # horizontal line at 0% deviation, so the sign is readable at a glance
DEVIATION_ZERO_LINE_WIDTH = 1.0

# Figure widths - see build_stiffness_figure(). The existing 3-panel
# figures (cam/nocam/combined) are 3 * STIFFNESS_PANEL_WIDTH_INCHES wide
# (18", unchanged from before this was a named constant); the simplified
# "_combined" figures (VARIANT 1: one panel only) get their own, wider
# single-panel width so 8+ trees plus the legend stack stay readable.
STIFFNESS_PANEL_WIDTH_INCHES = 6.0
STIFFNESS_COMBINED_WIDTH_INCHES = 10.0
# The full figures' RIGHT_MARGIN=0.80 reserves an ABSOLUTE legend width of
# 0.20 * 18" = 3.6" - plenty for the longest model name. Reusing 0.80 on a
# STIFFNESS_COMBINED_WIDTH_INCHES=10" figure would only reserve 2", too
# little (legend text was clipped at the figure edge) - so the combined
# figures get their own margin, sized for the same ~3.6" of legend room.
STIFFNESS_COMBINED_RIGHT_MARGIN = 0.62


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

# Every name assigned in the CONFIGURATION section above (between the
# "CONFIGURATION" and "CODE" headers), used by print_active_configuration()
# below (Stage 8 Part B3). TREES_SUFFIX/PROVENANCE_NEWEST_INPUT_MTIME/
# PROVENANCE_FILTERS_TEXT are deliberately left out - they are runtime
# state set once by main(), not something you edit, so their value here
# (before main() resolves anything) would just be a meaningless "" / None.
CONFIG_CONSTANT_NAMES = [
    "WORK_FOLDER", "DATA_FOLDER", "OUT_FOLDER", "EXPERIMENTAL_CSV_PATH", "FEM_CSV_PATH",
    "DAYS_TO_PROCESS", "MEASUREMENT", "TYPE_FILTER", "REDUCTION_LEVEL_FILTER", "TREE_FILTER",
    "TREE_ID_PREFIX", "ONLY_TREES_WITH_FEM", "DEVIATION_CSV_PATH", "DEVIATION_EXCLUDE_MODELS",
    "MODEL_INCLUDE", "MODEL_INCLUDE_SUFFIX", "FIGURE_PROVENANCE_FOOTER", "ABORT_ON_STALE_INPUT",
    "ANSYS_RAW_FOLDER", "ANSYS_RAW_PATTERN_FOR_FRESHNESS", "GRAVITY_VARIANT", "TREES",
    "STRICT_SELECTION", "GROUP_WIDTH", "PASTEL_SCENARIO_COLORS", "TITLE_FONTSIZE",
    "AXIS_LABEL_FONTSIZE", "TICK_LABEL_FONTSIZE", "LEGEND_FONTSIZE", "LEGEND_TITLE_FONTSIZE",
    "LINE_WIDTH", "MARKER_SIZE", "BOX_EDGE_LINEWIDTH", "LEGEND_BBOX_X", "LEGEND_Y_TOP",
    "LEGEND_Y_STEP", "RIGHT_MARGIN", "DEVIATION_HEIGHT_RATIO", "DEVIATION_BAR_WIDTH",
    "DEVIATION_BASE_HATCH", "DEVIATION_BASE2_HATCH", "DEVIATION_BAR_EDGECOLOR",
    "DEVIATION_BAR_EDGE_LINEWIDTH", "DEVIATION_ZERO_LINE_COLOR", "DEVIATION_ZERO_LINE_WIDTH",
    "STIFFNESS_PANEL_WIDTH_INCHES", "STIFFNESS_COMBINED_WIDTH_INCHES",
    "STIFFNESS_COMBINED_RIGHT_MARGIN",
]


def print_active_configuration(script_label: str, config_names: list) -> None:
    """
    Print the current value of every name in config_names, one per line -
    the exact, complete state of this script's CONFIGURATION section for
    this run (Stage 8 Part B3). Printed first, before any other output, so
    the console log alone says what produced a given figure - no need to
    re-open the script and check by eye.
    """
    print("\n" + "=" * 78)
    print(f"ACTIVE CONFIGURATION - {script_label}")
    print("=" * 78)
    for name in config_names:
        print(f"  {name} = {globals()[name]!r}")
    print("=" * 78)


def get_reduction_level(row) -> str:
    """
    Work out the reduction level for a single row of the table, from the
    three boolean reduction_20 / reduction_20+20 / reduction_40 columns
    (same logic as in filter_and_plot.py). Returns "20", "20+20", "40" or
    "none" (no reduction applied).
    """
    if row["reduction_20"]:
        return "20"
    elif row["reduction_20+20"]:
        return "20+20"
    elif row["reduction_40"]:
        return "40"
    else:
        return "none"


def add_reduction_level_column(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add a "reduction_level" column (see get_reduction_level above) if it
    is not already present in df.
    """
    if "reduction_level" in df.columns:
        return df
    df = df.copy()  # work on a copy so we never modify the original data
    df["reduction_level"] = df.apply(get_reduction_level, axis=1)
    return df


def filter_column(data: pd.DataFrame, column_name: str, value):
    """
    Generic filter helper (same convention as filter_and_plot.py):
      - value is None      -> no filtering, keep all rows
      - value is a list/tuple/set -> keep rows whose column is in it
      - otherwise           -> keep rows where the column equals value
    """
    if value is None:
        return data
    if isinstance(value, (list, tuple, set)):
        return data[data[column_name].isin(value)]
    return data[data[column_name] == value]


def tree_number_to_tree_id(tree_number: int) -> str:
    """
    Convert a plain tree number (e.g. 1, 7, 12) into the canonical project
    tree ID (e.g. "B21_S01", "B21_S07", "B21_S12"). :02d pads the number
    with a leading zero if it is a single digit.

    Copied from load_static_scenarios.py, not imported - these scripts are
    deliberately standalone with no shared imports.
    """
    return f"{TREE_ID_PREFIX}{tree_number:02d}"


def bk_id_to_canonical(bk_id: str) -> str:
    """
    Convert an experimental tree ID into the canonical project tree ID:
      "BK07" -> "B21_S07",  "BK12" -> "B21_S12"

    Only the digits are kept from the input, so the exact prefix in the
    experimental file does not matter. :02d pads a single digit with a
    leading zero.

    Why this is needed: the FEM table (slopes_combined.csv) uses the
    canonical form, the experimental table uses "BK..". Without this,
    every merge finds zero matches and the plots come out empty - and
    nothing errors, because "no matching rows" is not a crash.
    """
    digits_only = "".join(character for character in str(bk_id) if character.isdigit())
    return f"{TREE_ID_PREFIX}{int(digits_only):02d}"


def load_and_filter_experimental(csv_path) -> pd.DataFrame:
    """
    Load the experimental pulling-test table and filter it down to the
    rows configured above (days, measurement, type, reduction level,
    tree). Also adds "reduction_level" and "year" columns, and rewrites
    the "tree" column into the canonical B21_S.. form.
    """
    df = pd.read_csv(csv_path, low_memory=False)
    df = add_reduction_level_column(df)

    filtered = df[df["day"].isin(DAYS_TO_PROCESS)]
    filtered = filtered[filtered["measurement"] == MEASUREMENT]
    filtered = filtered[filtered["type"] == TYPE_FILTER]
    filtered = filter_column(filtered, "reduction_level", REDUCTION_LEVEL_FILTER)
    filtered = filter_column(filtered, "tree", TREE_FILTER)

    # The "year" column is used to group boxplots side by side within each
    # tree - two days in the same year would land in the same box.
    filtered = filtered.copy()
    filtered["year"] = pd.to_datetime(filtered["day"]).dt.year

    # Keep the original label in its own column, so nothing is lost, and
    # replace "tree" with the canonical form used everywhere else.
    filtered["tree_experimental_id"] = filtered["tree"]
    filtered["tree"] = filtered["tree"].apply(bk_id_to_canonical)

    print(f"Experimental data: {len(filtered)} rows after filtering")
    print(f"Experimental trees (renamed to canonical form): "
          f"{sorted(filtered['tree'].unique())}")
    return filtered


def pair_and_median(df: pd.DataFrame, col_a: str, col_b: str) -> pd.DataFrame:
    """
    Combine two per-sensor stiffness columns (e.g. slope_root_stiff_cam
    and slope_root_stiff_nocam) into one "combined" value per pull.

    Each physical pull is identified by (tree, day, measurement, pullNo).
    Within a pull, at most one row has a non-null value for col_a, and at
    most one row has a non-null value for col_b (different sensors write
    to different rows). This function collects whichever of the two
    values exist for a pull and takes their median:
      - both present  -> median = their average
      - only one present -> median = that single value
      - neither present  -> the pull is skipped entirely

    Returns a DataFrame with columns: tree, year, value - one row per pull
    that had at least one usable value.
    """
    grouped = df.groupby(["tree", "year", "day", "measurement", "pullNo"])

    combined_rows = []
    for group_key, group in grouped:
        tree, year = group_key[0], group_key[1]

        values_present = []
        value_a = group[col_a].dropna()
        if len(value_a) > 0:
            values_present.append(value_a.iloc[0])

        value_b = group[col_b].dropna()
        if len(value_b) > 0:
            values_present.append(value_b.iloc[0])

        if len(values_present) == 0:
            continue  # neither sensor recorded a value for this pull

        combined_rows.append({
            "tree": tree,
            "year": year,
            "value": float(np.median(values_present)),
        })

    return pd.DataFrame(combined_rows)


def load_fem_data(csv_path) -> pd.DataFrame:
    """
    Load the Stage 2 FEM results table (tree, scenario, model,
    slope_rotation_base, slope_rotation_base2, slope_displacement, ...).

    The "model" column carries the readable reconstruction name
    ("adtree_r10", "treeqsm_pd002", ...) and is what the plots use for
    colours and legends; "scenario" stays in the table as the raw ANSYS
    run number.
    """
    fem_df = pd.read_csv(csv_path)

    # Older Stage 2 outputs have no "model" column. Rather than failing,
    # fall back to the scenario number so the script still runs.
    if "model" not in fem_df.columns:
        print("WARNING: the FEM table has no 'model' column - falling back to "
              "scenario numbers. Re-run load_static_scenarios.py to get names.")
        fem_df["model"] = fem_df["scenario"].apply(lambda s: f"scenario{s}")

    print(f"FEM data: {len(fem_df)} row(s) ({fem_df['tree'].nunique()} tree(s), "
          f"models: {sorted(fem_df['model'].unique())})")
    return fem_df


def compute_dodge_offsets(n_groups: int, group_width: float) -> list:
    """
    Work out horizontal "dodge" offsets to place n_groups side-by-side
    sub-positions within one tree's x-axis slot, evenly spaced and
    centered on 0, spanning a total width of group_width.

    Example: compute_dodge_offsets(2, 0.6) -> [-0.3, 0.3]
    """
    if n_groups == 1:
        return [0.0]
    step = group_width / (n_groups - 1)
    return [-group_width / 2 + i * step for i in range(n_groups)]


def assign_model_colors(models: list) -> dict:
    """
    Assign a distinct, consistent color to each FEM model name, using the
    fixed PASTEL_SCENARIO_COLORS list (cycling with modulo if there are
    more models than colors), so the same model always gets the same
    color across all subplots and all figures.

    Sorting the names first means the colour of a model does not change
    when a new tree adds a model that was not there before.
    """
    n_colors = len(PASTEL_SCENARIO_COLORS)
    return {
        model: PASTEL_SCENARIO_COLORS[i % n_colors]
        for i, model in enumerate(sorted(models))
    }


def normalize_model_name_for_matching(name: str) -> str:
    """
    For COMPARISON only - never used to rewrite a displayed or written
    model name. Treats spaces and underscores as equivalent (collapsing
    runs of either into a single "_"), so a MODEL_INCLUDE entry spelled
    "adtree_r5mm_no_acel" still matches the data's actual
    "adtree_r5mm_no acel" (space) despite the inconsistent spelling in the
    source data - that inconsistency is being fixed in ANSYS, not here.
    """
    return re.sub(r"[ _]+", "_", str(name)).strip().lower()


# Matches either spelling of the no-acel suffix ("_no acel" with a space,
# "_no_acel" with an underscore) at the end of a model name. Copied from
# method_deviation_overview.py, not imported - these scripts are
# deliberately standalone with no shared imports.
NO_ACEL_PATTERN = re.compile(r"[ _]no[ _]acel$", re.IGNORECASE)


def is_no_acel_model(model_name: str) -> bool:
    """True for a gravity-off ("_no_acel"/"_no acel") model name."""
    return bool(NO_ACEL_PATTERN.search(str(model_name)))


def base_model_name(model_name: str) -> str:
    """Strip the no-acel suffix (either spelling), for matching only."""
    return NO_ACEL_PATTERN.sub("", str(model_name))


def filter_by_gravity_variant(df: pd.DataFrame, variant: str, model_col: str = "model") -> pd.DataFrame:
    """
    THE STRUCTURAL LOCK (see GRAVITY_VARIANT's own comment for why this
    exists): keep only the rows whose model matches the requested gravity
    variant. Called ONCE, in main(), immediately after fem_df is loaded -
    every figure and the CSV are built from what THIS function returns, so
    a row from the other variant cannot reach them, by construction rather
    than by remembering to filter correctly in N different places.

      "no_acel"      -> only _no_acel (gravity-off) models
      "with_gravity" -> only the base (gravity-on) models
      "both"         -> no filtering - returns df unchanged
    """
    if variant == "both":
        return df
    if variant == "no_acel":
        return df[df[model_col].apply(is_no_acel_model)]
    if variant == "with_gravity":
        return df[~df[model_col].apply(is_no_acel_model)]
    raise ValueError(f"Unknown GRAVITY_VARIANT: {variant!r} - use 'no_acel', 'with_gravity', or 'both'.")


def assert_gravity_variant(df: pd.DataFrame, expected_variant: str, function_name: str,
                           model_col: str = "model") -> None:
    """
    Loudly fail if df contains a model from the WRONG gravity variant - the
    check GRAVITY_VARIANT's structural lock is meant to make unnecessary in
    practice, kept as an assertion so a future code change that bypasses
    filter_by_gravity_variant() is caught immediately, at the one place a
    figure gets drawn, rather than discovered later by eye in a PNG.

    expected_variant="both" always passes (no restriction). Otherwise every
    model in df must agree with expected_variant's no_acel-ness; any that
    do not are named in the AssertionError, along with this function's own
    name and the variant that was expected.
    """
    if expected_variant == "both" or len(df) == 0:
        return
    models = df[model_col].unique()
    if expected_variant == "no_acel":
        offenders = [m for m in models if not is_no_acel_model(m)]
    elif expected_variant == "with_gravity":
        offenders = [m for m in models if is_no_acel_model(m)]
    else:
        raise ValueError(f"Unknown gravity variant to assert: {expected_variant!r}")
    assert not offenders, (
        f"{function_name}: expected ONLY '{expected_variant}' models but found "
        f"{offenders} - the gravity-variant structural lock (filter_by_gravity_variant, "
        f"applied once in main()) was bypassed somewhere upstream of this function."
    )


def resolve_model_filter(all_models: list, include_filter) -> list:
    """
    Resolve MODEL_INCLUDE (None, or a list of wanted model names) against
    the actual distinct model names in the data (all_models), using
    normalize_model_name_for_matching() for tolerant comparison. The
    returned list always uses the ORIGINAL spelling from all_models -
    nothing is rewritten.

    Prints every distinct model found, and which of them MODEL_INCLUDE
    matched. An include_filter entry matching nothing gets a WARNING (not
    an error), naming it.

    include_filter=None -> returns all_models unchanged (every model drawn).
    """
    print(f"Distinct model names found in the data: {sorted(all_models)}")
    if include_filter is None:
        print("MODEL_INCLUDE is None - every model is drawn.")
        return list(all_models)

    normalized_to_actual = {}
    for model in all_models:
        normalized_to_actual.setdefault(normalize_model_name_for_matching(model), []).append(model)

    matched = []
    for wanted in include_filter:
        key = normalize_model_name_for_matching(wanted)
        found = normalized_to_actual.get(key, [])
        if found:
            print(f"  MODEL_INCLUDE entry '{wanted}' matched: {found}")
            matched.extend(found)
        else:
            print(f"  WARNING: MODEL_INCLUDE entry '{wanted}' matched no model in the data.")

    matched = sorted(set(matched))
    print(f"Models drawn (MODEL_INCLUDE active): {matched}")
    return matched


# Rejects the ANSYS row-position spelling outright (lowercase "bk<N>",
# e.g. copy-pasted from "out_static_bk3_scenario1.txt") - see
# normalize_tree_selector()'s own docstring for why this must never be
# read as a tree number.
BK_INDEX_REJECT_PATTERN = re.compile(r"^bk\d+$")


def normalize_tree_selector(value) -> str:
    """
    Convert one TREES entry into the canonical tree ID "B21_Snn".

    Accepts: 7, "7", "07", "S07", "s07", "B21_S07", "BK07" (uppercase only).
    Rejects: "bk3"-style (lowercase "bk" + digits) with a clear error.

    WHY the BK/bk case matters (read before touching this function): ANSYS
    filenames use a LOWERCASE "bk<N>" for the tree's ROW POSITION in the
    description file (bk3 is position 3, i.e. whichever tree happens to sit
    at that row - see load_ansys_order_to_tree_id()'s own docstring). That
    number is NOT a tree number and must never be typed into TREES. An
    UPPERCASE "BK07" is accepted here only as another way of writing the
    TREE NUMBER 7 - deliberately spelled differently from the lowercase
    ANSYS row-position so a copy-pasted "bk3" from a filename is caught
    immediately instead of silently selecting the wrong tree.
    """
    text = str(value).strip()

    if BK_INDEX_REJECT_PATTERN.match(text):
        raise ValueError(
            f"TREES entry {value!r} looks like an ANSYS row-position index "
            f"(lowercase 'bk<N>', as it appears in a filename like "
            f"'out_static_bk3_scenario1.txt'). bk<N> is that tree's ROW "
            f"POSITION in the description file, NOT its tree number - see "
            f"load_ansys_order_to_tree_id()'s docstring. Use the tree's own "
            f"number instead, e.g. 7, 'S07', or 'B21_S07'."
        )

    for pattern in (r"^B21_S0*(\d+)$", r"^BK0*(\d+)$"):
        match = re.fullmatch(pattern, text)
        if match:
            return tree_number_to_tree_id(int(match.group(1)))

    match = re.fullmatch(r"^[sS]0*(\d+)$", text)
    if match:
        return tree_number_to_tree_id(int(match.group(1)))

    match = re.fullmatch(r"^\d+$", text)
    if match:
        return tree_number_to_tree_id(int(text))

    raise ValueError(
        f"TREES entry {value!r} is not a recognised tree selector. Use a "
        f"plain number (7), 'S07', 's07', 'B21_S07', or 'BK07'."
    )


def resolve_tree_selection(trees_config, available_trees: list) -> dict:
    """
    Resolve TREES ("all" or a list of mixed-form entries) against the trees
    actually present in the data (available_trees).

    Returns a dict:
      requested_display - what to print for "requested trees" (e.g. "all"
                           or the resolved canonical list)
      resolved          - sorted list of canonical tree IDs to actually use
                           (found in available_trees)
      missing           - sorted list of requested tree IDs NOT found in
                           available_trees (empty when trees_config="all")
    "all" never produces a MISSING entry by construction - it can only ask
    for what is already there.
    """
    if trees_config == "all":
        resolved = sorted(available_trees)
        return {"requested_display": "all", "resolved": resolved, "missing": [],
                "is_all": True}

    requested = sorted({normalize_tree_selector(v) for v in trees_config})
    available_set = set(available_trees)
    resolved = sorted(t for t in requested if t in available_set)
    missing = sorted(t for t in requested if t not in available_set)
    return {"requested_display": requested, "resolved": resolved, "missing": missing,
            "is_all": False}


def apply_tree_selection(df: pd.DataFrame, tree_sel: dict, tree_col: str = "tree") -> pd.DataFrame:
    """
    Filter df down to tree_sel["resolved"]. is_all=True keeps every row
    unchanged (no-op) - this is what keeps an "everything=all" run's data
    identical to before this feature existed.
    """
    if not tree_sel["is_all"]:
        df = df[df[tree_col].isin(tree_sel["resolved"])]
    return df


def build_filter_suffix(tree_sel: dict) -> str:
    """
    Build the deterministic filename-suffix fragment for an active TREES
    selection, e.g. tree_sel resolved to B21_S07/B21_S12 -> "trees-07-12".
    Returns "" when "all" (no-op - keeps default filenames unchanged by
    this particular filter). Numbers only (not the full "B21_Snn" string)
    to keep filenames short.

    NOTE (Stage 8): this used to also carry a SCENARIOS fragment ("sc-3-4")
    - removed along with SCENARIOS itself (see STRICT_SELECTION's own
    comment for why: fully redundant with GRAVITY_VARIANT x MODEL_INCLUDE).
    """
    if tree_sel["is_all"]:
        return ""
    numbers = sorted(t.split("_S")[-1] for t in tree_sel["resolved"])
    return "trees-" + "-".join(numbers)


def print_selection_report(script_label: str, tree_sel: dict) -> bool:
    """
    Print the SELECTION REPORT block (trees only - see STRICT_SELECTION's
    own comment for why the scenario section was retired in Stage 8) and
    return True if STRICT_SELECTION should raise for this script (i.e.
    something requested is genuinely missing from the data).
    """
    print(f"\nSELECTION REPORT - {script_label}")

    if tree_sel["is_all"]:
        print(f"  requested trees      : all  ->  resolved to {len(tree_sel['resolved'])}")
    else:
        print(f"  requested trees      : {tree_sel['requested_display']}")
    print(f"  found in data        : {len(tree_sel['resolved'])}")
    print(f"  plotted              : {len(tree_sel['resolved'])}")
    if tree_sel["missing"]:
        print(f"  requested but MISSING : {', '.join(tree_sel['missing'])}")
    else:
        print(f"  requested but MISSING : (none)")

    return bool(tree_sel["missing"])


def newest_mtime(folder, pattern: str):
    """
    Return (newest_mtime_float, newest_path) among files in folder matching
    pattern, or (None, None) if nothing matches. Used by
    check_input_freshness() and add_provenance_footer().
    """
    matches = sorted(Path(folder).glob(pattern))
    if not matches:
        return None, None
    newest_path = max(matches, key=lambda p: p.stat().st_mtime)
    return newest_path.stat().st_mtime, newest_path


def check_input_freshness(script_label: str, raw_folder, raw_pattern: str,
                          csv_paths_being_read: dict, abort_if_stale: bool) -> float:
    """
    Compare the newest modification time among raw ANSYS files (raw_folder/
    raw_pattern) against every CSV path in csv_paths_being_read ({label:
    path}). If any such CSV is OLDER than the newest raw file, print a
    clearly visible warning naming both files, and - if abort_if_stale -
    raise SystemExit instead of plotting from data that does not yet
    reflect the latest ANSYS results.

    Returns the newest raw-file mtime found (or None), for reuse by the
    provenance footer, so "newest input" always means the same thing on the
    figure as it did in this check.
    """
    newest_raw_mtime, newest_raw_path = newest_mtime(raw_folder, raw_pattern)
    if newest_raw_mtime is None:
        print(f"FRESHNESS CHECK ({script_label}): no files match '{raw_pattern}' in "
              f"{raw_folder} - skipping.")
        return None

    stale_found = False
    for label, csv_path in csv_paths_being_read.items():
        csv_path = Path(csv_path)
        if not csv_path.is_file():
            continue  # a missing file is reported elsewhere (FileNotFoundError on read)
        csv_mtime = csv_path.stat().st_mtime
        if csv_mtime < newest_raw_mtime:
            stale_found = True
            print(f"*** WARNING: STALE INPUT - {label} ({csv_path.name}, "
                  f"{datetime.fromtimestamp(csv_mtime):%Y-%m-%d %H:%M:%S}) is OLDER than "
                  f"{newest_raw_path.name} ({datetime.fromtimestamp(newest_raw_mtime):%Y-%m-%d %H:%M:%S}). "
                  f"The pipeline that builds {csv_path.name} has not been re-run since "
                  f"that raw ANSYS file last changed - this script would be plotting "
                  f"data that does not yet reflect it. ***")

    if stale_found and abort_if_stale:
        raise SystemExit(
            f"{script_label}: aborting because ABORT_ON_STALE_INPUT=True and at least one "
            f"input CSV is older than a raw ANSYS file it depends on (see WARNING(s) above). "
            f"Re-run the upstream pipeline (load_static_scenarios.py / compute_slopes.py) "
            f"first, or set ABORT_ON_STALE_INPUT=False to plot anyway."
        )
    return newest_raw_mtime


# Filter-suffix markers this family of scripts writes into filenames - used
# only to FIND and DESCRIBE a previously-written filtered output, never to
# read or write its contents. Extend this tuple if a new suffix style is
# ever introduced (the trees/scenarios suffix is matched separately, by
# TREES_SUFFIX_PATTERN, since it embeds tree/scenario numbers).
KNOWN_MODEL_SUFFIX_MARKERS = ("_noacel",)
TREES_SUFFIX_PATTERN = re.compile(r"_trees-[0-9-]+|_sc-[0-9-]+")


def find_filtered_output_files(out_folder, name_prefixes: tuple) -> list:
    """
    Return every FILE (not directory) in out_folder whose name starts with
    one of name_prefixes (this script's own output family) and carries a
    filter-suffix marker anywhere in its name (KNOWN_MODEL_SUFFIX_MARKERS,
    or a "_trees-.../_sc-..." tree/scenario suffix).
    """
    hits = []
    for path in sorted(Path(out_folder).iterdir()):
        if not path.is_file() or not path.name.startswith(name_prefixes):
            continue
        if any(marker in path.name for marker in KNOWN_MODEL_SUFFIX_MARKERS) \
                or TREES_SUFFIX_PATTERN.search(path.name):
            hits.append(path)
    return hits


def unfiltered_sibling_path(filtered_path) -> Path:
    """
    Strip every known filter-suffix marker from filtered_path's name to
    recover the UNFILTERED filename it was derived from, e.g.
    "stiffness_comparison_rotation_all_noacel.png" ->
    "stiffness_comparison_rotation_all.png". Used only to COMPARE
    timestamps in report_orphaned_filtered_outputs() - never to read or
    write that file automatically.
    """
    filtered_path = Path(filtered_path)
    stem = filtered_path.stem
    for marker in KNOWN_MODEL_SUFFIX_MARKERS:
        stem = stem.replace(marker, "")
    stem = TREES_SUFFIX_PATTERN.sub("", stem)
    return filtered_path.with_name(stem + filtered_path.suffix)


def report_orphaned_filtered_outputs(out_folder, name_prefixes: tuple, script_label: str) -> None:
    """
    Print (NEVER touch - no delete/move/rename) every existing filtered-
    output file in out_folder whose UNFILTERED sibling is newer than it -
    i.e. the unfiltered figure has been regenerated since this filtered one
    was last written, so the filtered one is now stale and, with nothing on
    either file saying which is which, could be mistaken for current. This
    is exactly the Sep 2026 incident that motivated this whole safeguard.
    """
    hits = find_filtered_output_files(out_folder, name_prefixes)
    if not hits:
        return
    print(f"\nORPHANED FILTERED OUTPUTS ({script_label}) - these were NOT (necessarily) "
          f"regenerated by this run and may be stale. NOTHING is touched - review them "
          f"by hand before using them:")
    for path in hits:
        sibling = unfiltered_sibling_path(path)
        mtime_text = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        if sibling.exists():
            sibling_mtime = sibling.stat().st_mtime
            if sibling_mtime > path.stat().st_mtime:
                sibling_text = datetime.fromtimestamp(sibling_mtime).strftime("%Y-%m-%d %H:%M")
                status = f"STALE - its unfiltered sibling ({sibling.name}) is newer ({sibling_text})"
            else:
                status = "unfiltered sibling is not newer than this file - may still be current"
        else:
            status = "no unfiltered sibling found to compare against"
        print(f"  {path}   {mtime_text}   {status}")


def add_provenance_footer(fig, script_filename: str, newest_input_mtime, active_filters_text: str,
                          n_trees: int, n_models: int) -> None:
    """
    Stamp a small grey monospace footer along the bottom edge of fig,
    below the axes, with: this script's filename, when the figure was
    generated (now), the newest modification time among the raw ANSYS
    files it was built from, the active filters, and how many trees/models
    were actually drawn - see FIGURE_PROVENANCE_FOOTER's own comment for
    why. A no-op when FIGURE_PROVENANCE_FOOTER is False.

    Stage 9 Part A1: also stashes the SAME provenance as fig.stage9_png_metadata,
    picked up by savefig_with_provenance() below when the figure is actually
    written to disk. The footer above is pixels only - reading it back needs
    OCR. This is a real PNG tEXt chunk, readable afterwards via PIL's
    Image.info without OCR or guessing a run's configuration from its
    filename - see check_png_provenance.py.
    """
    if not FIGURE_PROVENANCE_FOOTER:
        fig.stage9_png_metadata = {}
        return
    generated_text = datetime.now().strftime("%Y-%m-%d %H:%M")
    newest_input_text = (datetime.fromtimestamp(newest_input_mtime).strftime("%Y-%m-%d %H:%M")
                         if newest_input_mtime is not None else "unknown")
    line1 = f"{script_filename} | generated {generated_text} | newest input {newest_input_text}"
    line2 = f"{active_filters_text}  |  {n_trees} trees, {n_models} models"
    fig.text(0.01, -0.045, line1 + "\n" + line2, ha="left", va="bottom",
             fontsize=6, family="monospace", color="0.55")

    fig.stage9_png_metadata = {
        "Script": script_filename,
        "Generated": generated_text,
        "NewestInput": newest_input_text,
        "GravityVariant": str(GRAVITY_VARIANT),
        "ModelInclude": str(MODEL_INCLUDE),
        "Trees": str(TREES),
        "NTrees": str(n_trees),
        "NModels": str(n_models),
    }


def savefig_with_provenance(fig, path) -> None:
    """
    The ONLY way a figure should be written to disk in this script (Stage 9
    Part A1) - every build_*() function calls add_provenance_footer(fig,
    ...) before saving, which stashes fig.stage9_png_metadata; this passes
    that dict straight to savefig(metadata=...), embedding it as real PNG
    metadata alongside the existing visible footer. A PNG with no such
    metadata predates Stage 9 and its configuration cannot be verified this
    way - see check_png_provenance.py.
    """
    plt.savefig(path, dpi=150, bbox_inches="tight", metadata=getattr(fig, "stage9_png_metadata", {}))


def darker(hex_color: str, factor: float = 0.7):
    """
    Return a darker version of hex_color by scaling its RGB components by
    factor (0-1 - smaller means darker). Used to give each boxplot a subtle
    edge color that matches its own fill, instead of plain black.

    Example: darker("#94BBEA", 0.7) -> a slightly darker pastel blue.
    """
    r, g, b = mcolors.to_rgb(hex_color)
    return (r * factor, g * factor, b * factor)


def plot_experimental_boxplots(ax, df: pd.DataFrame, value_col: str, tree_order: list, year_order: list,
                               show_legend: bool = True):
    """
    Draw grey boxplots of df[value_col] onto ax: one x-position per tree
    (from tree_order), and within each tree's slot, one box per year (from
    year_order), placed side by side. Boxes are shaded from light to dark
    grey by year (earliest year = lightest), with a small legend mapping
    grey shade to year, drawn outside the axes (to the right).

    df must have columns: "tree", "year", value_col.

    show_legend=False skips drawing the legend - used when several subplots
    sit side by side (see build_rotation_figure/build_displacement_figure),
    so only the LAST subplot's external legend lands in the figure's own
    empty right margin instead of on top of a neighbouring subplot.
    """
    # Light background grid, drawn first (low zorder) so the boxplots and
    # FEM points are layered on top of it, not hidden behind it.
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=0)

    # Thin vertical separators between trees, at the midpoint between one
    # tree's slot and the next, so it is clear where each tree's data ends.
    for tree_index in range(len(tree_order) - 1):
        ax.axvline(tree_index + 0.5, color="#e1e0d9", linewidth=0.6, zorder=0)

    offsets = compute_dodge_offsets(len(year_order), GROUP_WIDTH)
    box_width = GROUP_WIDTH / len(year_order) * 0.8

    # Light-to-dark grey shades, one per year, avoiding pure white/black.
    grey_shades = {
        year: str(shade) for year, shade in
        zip(year_order, np.linspace(0.85, 0.45, len(year_order)))
    }

    for year, offset in zip(year_order, offsets):
        positions = []
        data_lists = []

        for tree_index, tree in enumerate(tree_order):
            values = df[(df["tree"] == tree) & (df["year"] == year)][value_col].dropna()
            if len(values) == 0:
                continue  # no data for this tree/year combination - leave the slot empty
            positions.append(tree_index + offset)
            data_lists.append(values.values)

        if len(data_lists) == 0:
            continue

        box_plot = ax.boxplot(
            data_lists, positions=positions, widths=box_width,
            patch_artist=True, manage_ticks=False, zorder=3,
        )
        # Each box's edge is a darker shade of its own fill color, instead
        # of plain black, for a softer look.
        box_edge_color = darker(grey_shades[year])
        for patch in box_plot["boxes"]:
            patch.set_facecolor(grey_shades[year])
            patch.set_edgecolor(box_edge_color)
            patch.set_linewidth(BOX_EDGE_LINEWIDTH)
        for median_line in box_plot["medians"]:
            median_line.set_color(box_edge_color)

    # Small legend mapping grey shade to year, moved outside the axes.
    # Kept as a separate "artist" so a later legend (e.g. from
    # overlay_fem_points) does not replace it.
    if show_legend:
        year_handles = [
            Patch(facecolor=grey_shades[year], edgecolor=darker(grey_shades[year]), label=str(year))
            for year in year_order
        ]
        year_legend = ax.legend(handles=year_handles, title="Year", loc="upper left",
                                bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP), frameon=False,
                                fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)
        ax.add_artist(year_legend)

    ax.set_xticks(range(len(tree_order)))
    ax.set_xticklabels(tree_order)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)


def overlay_fem_points(ax, fem_df: pd.DataFrame, value_specs: list, tree_order: list,
                       show_legend: bool = True, model_colors: dict = None):
    """
    Overlay FEM model results onto ax as colored points, on the same
    tree x-axis as the experimental boxplots.

    value_specs is a list of (column_name, marker_shape, label) tuples,
    e.g. [("slope_rotation_base", "o", "rotsum base"),
          ("slope_rotation_base2", "^", "rotsum base2")].
    Points are colored by reconstruction model (consistent across
    subplots via assign_model_colors), and their marker shape shows which
    FEM column they come from. Legends are drawn outside the axes (to the
    right), stacked below plot_experimental_boxplots()'s year legend.

    THE STRUCTURAL LOCK'S ASSERTION runs first (see assert_gravity_variant()):
    every row this function draws must already agree with GRAVITY_VARIANT.

    show_legend=False skips drawing the legend(s) - see
    plot_experimental_boxplots() for why.

    model_colors, when given, is used AS-IS instead of being recomputed
    from fem_df["model"] - this is what keeps a model's colour fixed
    across a MODEL_INCLUDE-filtered run and the full run (colours are
    assigned alphabetically over the FULL model set; recomputing them from
    an already-filtered fem_df would shift them). The legend still only
    lists models actually PRESENT in this fem_df, so a filtered-out
    model's colour is reserved but its legend entry does not appear.
    """
    assert_gravity_variant(fem_df, GRAVITY_VARIANT, "overlay_fem_points")

    if model_colors is None:
        model_colors = assign_model_colors(fem_df["model"].unique())

    for _, row in fem_df.iterrows():
        tree = row["tree"]
        if tree not in tree_order:
            print(f"WARNING: FEM tree '{tree}' (model {row['model']}) is not in the "
                  f"experimental tree list - skipping its points.")
            continue
        x = tree_order.index(tree)
        color = model_colors.get(row["model"])
        if color is None:
            continue  # not in the provided palette - should not normally happen

        for column_name, marker_shape, _ in value_specs:
            value = row[column_name]
            if pd.isna(value):
                continue
            ax.scatter(x, value, marker=marker_shape, color=color,
                       edgecolor="none", s=MARKER_SIZE ** 2, zorder=5)

    if not show_legend:
        return

    # Legend 1: color -> reconstruction model. Stacked below the year
    # legend added by plot_experimental_boxplots() on this same axes.
    # Restricted to models actually present in fem_df (not every key of
    # model_colors), so a MODEL_INCLUDE-filtered figure's legend lists only
    # what is actually drawn.
    present_models = [m for m in sorted(model_colors.keys()) if m in set(fem_df["model"].unique())]
    scenario_handles = [
        Line2D([], [], marker="o", linestyle="", color="black",
               markerfacecolor=model_colors[model], markersize=MARKER_SIZE, label=model)
        for model in present_models
    ]
    scenario_legend = ax.legend(handles=scenario_handles, title="Model",
                                 loc="upper left", bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - LEGEND_Y_STEP),
                                 frameon=False, fontsize=LEGEND_FONTSIZE,
                                 title_fontsize=LEGEND_TITLE_FONTSIZE)

    # Legend 2: marker shape -> which FEM column, only needed if there is
    # more than one column being compared. Stacked one slot further down.
    # Kept as an "artist" too (even though it is currently the last legend
    # call here) so a caller building the rotation figure's deviation-bar
    # hatch legend can stack a further legend after this one without it
    # being silently replaced.
    if len(value_specs) > 1:
        ax.add_artist(scenario_legend)  # keep the scenario legend visible
        marker_handles = [
            Line2D([], [], marker=marker_shape, linestyle="", color="black",
                   markerfacecolor="white", markersize=MARKER_SIZE, label=label)
            for _, marker_shape, label in value_specs
        ]
        marker_legend = ax.legend(handles=marker_handles, title="FEM value", loc="upper left",
                                  bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - 2 * LEGEND_Y_STEP), frameon=False,
                                  fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)
        ax.add_artist(marker_legend)


def compute_deviation_table(experimental: pd.DataFrame, combined: pd.DataFrame,
                            fem_df: pd.DataFrame, sensor_columns: dict,
                            fem_columns: list) -> pd.DataFrame:
    """
    Build a table of relative deviations: FEM value against the median of
    the experimental values, for every tree / model / FEM column, and for
    every experimental sensor variant.

    Arguments:
      experimental    - the filtered experimental table (one row per probe)
      combined        - output of pair_and_median(), columns tree/year/value
      fem_df          - the Stage 2 FEM table
      sensor_columns  - dict {label: column name in `experimental`}, e.g.
                        {"cam": "slope_root_stiff_cam", ...}. The label
                        "combined" is taken from `combined` instead.
      fem_columns     - list of FEM column names to compare, e.g.
                        ["slope_rotation_base", "slope_rotation_base2"]

    Relative deviation is defined as

        rel_dev_pct = 100 * (FEM - experimental_median) / experimental_median

    so a positive number means the model is STIFFER than the measurement
    and a negative number means it is softer. The experimental median is
    taken across all pulls and all years for that tree - the boxplots in
    the figure show that spread, this number does not.
    """
    rows = []

    for _, fem_row in fem_df.iterrows():
        tree = fem_row["tree"]

        for sensor_label, sensor_column in sensor_columns.items():
            # Collect the experimental values for this tree and sensor.
            if sensor_label == "combined":
                values = combined[combined["tree"] == tree]["value"].dropna()
            else:
                values = experimental[experimental["tree"] == tree][sensor_column].dropna()

            n_values = len(values)
            experimental_median = values.median() if n_values > 0 else float("nan")

            for fem_column in fem_columns:
                fem_value = fem_row[fem_column]

                # Guard against dividing by zero or by a missing median.
                if pd.isna(experimental_median) or experimental_median == 0 or pd.isna(fem_value):
                    rel_dev_pct = float("nan")
                else:
                    rel_dev_pct = 100.0 * (fem_value - experimental_median) / experimental_median

                rows.append({
                    "tree": tree,
                    "model": fem_row["model"],
                    "fem_column": fem_column,
                    "fem_value": fem_value,
                    "sensor": sensor_label,
                    "exp_median": experimental_median,
                    "exp_n": n_values,
                    "rel_dev_pct": rel_dev_pct,
                })

    table = pd.DataFrame(rows)
    if len(table) > 0:
        table = table.round({"fem_value": 4, "exp_median": 4, "rel_dev_pct": 1})
    return table


def draw_deviation_bars(ax_dev, deviation_table: pd.DataFrame, bar_specs: list, model_colors: dict,
                        tree_order: list, exclude_models: list, title_prefix: str,
                        visible_models: list = None) -> list:
    """
    Draw the deviation subplot: grouped bars of relative deviation (%),
    combined-sensor rows only (the same rows the old text footnote used),
    one bar per model per FEM quantity in bar_specs (a list of
    (fem_column, hatch, label) tuples - hatch=None means solid fill).

    Models are dodged with GROUP_WIDTH across ALL models in
    model_colors (or, if MODEL_INCLUDE narrowed what is drawn,
    visible_models - a more compact dodge, since the whole point of that
    filter is a less dense figure), regardless of exclude_models - so a
    model's x-slot is always in the same place whether or not it happens
    to be excluded in THIS call; an excluded model's slot is simply left
    EMPTY, not compressed away, so the gap itself stays legible (its
    legend entry, drawn by the caller, says it exists and is deliberately
    absent).

    visible_models=None means "no MODEL_INCLUDE filter active" - every
    model in model_colors is eligible (subject only to exclude_models, as
    before). A model outside visible_models is skipped WITHOUT being
    added to the returned excluded_rows list - that list is specifically
    for DEVIATION_EXCLUDE_MODELS (a different, independent reason for a
    gap) - but if the two filters COMBINED leave nothing at all drawn for
    a quantity that actually has data, a WARNING is printed naming it.

    Returns the list of rows (as dicts) that were excluded from drawing
    because their model is in exclude_models - the caller prints these so
    no value is ever dropped silently, even from the "trimmed" figure.
    """
    ax_dev.set_axisbelow(True)
    ax_dev.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=0)
    ax_dev.axhline(0.0, color=DEVIATION_ZERO_LINE_COLOR, linewidth=DEVIATION_ZERO_LINE_WIDTH, zorder=2)

    all_models = sorted(visible_models) if visible_models is not None else sorted(model_colors.keys())
    offsets = compute_dodge_offsets(len(all_models), GROUP_WIDTH)
    offset_by_model = dict(zip(all_models, offsets))

    quantity_names = [column for column, _, _ in bar_specs]
    n_quantities = len(bar_specs)
    # A single quantity (displacement) gets a fuller bar centered on the
    # model's slot instead of a half-width bar with an empty sibling slot.
    bar_width = DEVIATION_BAR_WIDTH if n_quantities > 1 else DEVIATION_BAR_WIDTH * 2
    sub_offsets = compute_dodge_offsets(n_quantities, bar_width * n_quantities)

    combined_rows = deviation_table[deviation_table["sensor"] == "combined"]
    tree_index_by_tree = {tree: i for i, tree in enumerate(tree_order)}

    models_with_data = sorted({
        row["model"] for _, row in combined_rows.iterrows()
        if row["tree"] in tree_index_by_tree and row["fem_column"] in quantity_names
        and pd.notna(row["rel_dev_pct"])
    })

    excluded_rows = []
    drawn_any = False
    for _, row in combined_rows.iterrows():
        if row["tree"] not in tree_index_by_tree or row["fem_column"] not in quantity_names:
            continue
        if visible_models is not None and row["model"] not in visible_models:
            continue   # not selected by MODEL_INCLUDE - a different, independent reason for a gap
        if row["model"] in exclude_models:
            excluded_rows.append(row.to_dict())
            continue   # gap left deliberately - see this function's own docstring
        if pd.isna(row["rel_dev_pct"]):
            continue
        drawn_any = True
        q_index = quantity_names.index(row["fem_column"])
        x = (tree_index_by_tree[row["tree"]] + offset_by_model.get(row["model"], 0.0)
             + sub_offsets[q_index])
        ax_dev.bar(x, row["rel_dev_pct"], width=bar_width, color=model_colors[row["model"]],
                   hatch=bar_specs[q_index][1], edgecolor=DEVIATION_BAR_EDGECOLOR,
                   linewidth=DEVIATION_BAR_EDGE_LINEWIDTH, zorder=3)

    if models_with_data and not drawn_any:
        print(f"  WARNING: MODEL_INCLUDE ({visible_models}) and DEVIATION_EXCLUDE_MODELS "
              f"({exclude_models}) between them remove EVERY model with data "
              f"({models_with_data}) from this deviation subplot - it will be empty.")

    ax_dev.set_xticks(range(len(tree_order)))
    ax_dev.set_xticklabels(tree_order)
    ax_dev.set_ylabel("Relative deviation (%)", fontsize=AXIS_LABEL_FONTSIZE)
    ax_dev.tick_params(labelsize=TICK_LABEL_FONTSIZE)

    title = title_prefix
    if exclude_models:
        title = f"{title_prefix} - {', '.join(sorted(exclude_models))} not shown"
    ax_dev.set_title(title, fontsize=TITLE_FONTSIZE)

    return excluded_rows


def count_fem_points(fem_df: pd.DataFrame, value_specs: list, tree_order: list) -> int:
    """
    Count exactly the points overlay_fem_points() would draw on ONE
    subplot for this fem_df/value_specs (same skip rules: tree not in
    tree_order, or a NaN value, are not drawn) - used only to print a
    verifiable point count, not to change what gets plotted.
    """
    count = 0
    for _, row in fem_df.iterrows():
        if row["tree"] not in tree_order:
            continue
        for column_name, _, _ in value_specs:
            if pd.notna(row[column_name]):
                count += 1
    return count


def save_stiffness_figure_variants(fig, ax_dev, deviation_table: pd.DataFrame, bar_specs: list,
                                   model_colors: dict, tree_order: list, quantity: str, title_prefix: str,
                                   right_margin: float = RIGHT_MARGIN, visible_models: list = None) -> None:
    """
    Render the deviation subplot twice into the SAME already-built figure
    (the 3 stiffness panels above are drawn ONCE by the caller and never
    touched here) and save two files:
        stiffness_comparison_<quantity>_all.png       every model shown
        stiffness_comparison_<quantity>_trimmed.png   DEVIATION_EXCLUDE_MODELS
                                                       excluded from the
                                                       SUBPLOT only
    Reusing the same fig/upper-axes object for both saves guarantees the
    stiffness panels are IDENTICAL between the two files BY CONSTRUCTION -
    not just "rebuilt the same way twice" - which is what this task's
    verification (d) asks to be confirmed.

    If DEVIATION_EXCLUDE_MODELS is empty, "trimmed" ends up identical to
    "all" - both are still written, so the file set is the same every run.

    Every excluded row is printed (tree/model/fem_column/rel_dev_pct) so a
    value never disappears from view without a trace, even though it is
    left out of the trimmed subplot's bars.
    """
    for variant, exclude_models in (("all", []), ("trimmed", DEVIATION_EXCLUDE_MODELS)):
        ax_dev.clear()
        excluded_rows = draw_deviation_bars(ax_dev, deviation_table, bar_specs, model_colors,
                                            tree_order, exclude_models, title_prefix,
                                            visible_models=visible_models)

        quantity_label = quantity.replace("_", " ")
        if exclude_models:
            print(f"\nExcluded from the {quantity_label} deviation subplot ('{variant}' figure, "
                  f"DEVIATION_EXCLUDE_MODELS={DEVIATION_EXCLUDE_MODELS}) - still in "
                  f"stiffness_deviation.csv and in the stiffness panels above, just not "
                  f"drawn as a bar here:")
            if excluded_rows:
                for row in excluded_rows:
                    value_text = ("no experimental data" if pd.isna(row["rel_dev_pct"])
                                  else f"{row['rel_dev_pct']:+.1f} %")
                    print(f"  {row['tree']} / {row['model']} / {row['fem_column']}: "
                          f"rel_dev_pct = {value_text}")
            else:
                print("  (none of the excluded model(s) appear in this quantity's data)")

        path = stiffness_plot_path(quantity, variant)
        plt.tight_layout()
        fig.subplots_adjust(right=right_margin)
        savefig_with_provenance(fig, path)
        print(f"{quantity_label.capitalize()} comparison plot ({variant}) saved to: {path}")


def build_stiffness_figure(quantity: str, panel_specs: list, fem_df: pd.DataFrame, value_specs: list,
                           bar_specs: list, model_colors: dict, tree_order: list, year_order: list,
                           ylabel: str, deviation_table: pd.DataFrame, figure_width: float,
                           title_prefix: str, print_label: str, right_margin: float = RIGHT_MARGIN,
                           visible_models: list = None) -> None:
    """
    Shared builder for every stiffness-comparison figure: the full 3-panel
    ones (cam/nocam/combined) AND the simplified combined-only ones
    (VARIANT 1) are built through this ONE function, parameterised by
    panel_specs/figure_width, rather than duplicated per variant.

    panel_specs: list of (title, df, value_col) tuples for the top row of
    stiffness subplots - 3 entries for the full figures, 1
    ("combined", combined_df, "value") for the simplified ones. Only the
    LAST panel gets legends (existing convention: drawing them on every
    panel would place one panel's "outside to the right" legend on top of
    its neighbour).

    deviation_table is passed in ALREADY COMPUTED, not recomputed here -
    the full and simplified figures for a given quantity are built from
    the exact SAME table, so their deviation numbers can never disagree
    with each other or with stiffness_deviation.csv. Its bars always come
    from the "combined" sensor rows regardless of which top panels are
    shown - draw_deviation_bars() (called via save_stiffness_figure_variants()
    below) filters on sensor=="combined" unconditionally, so the
    simplified figures' deviation subplot matches their single "combined"
    panel above it even though per-sensor rows are also present in the
    table that was passed in.

    Builds and saves TWO files via save_stiffness_figure_variants():
    stiffness_comparison_<quantity>_all.png / _trimmed.png.
    """
    n_panels = len(panel_specs)
    fig = plt.figure(figsize=(figure_width, 6 * (1.0 + DEVIATION_HEIGHT_RATIO)))
    gs = fig.add_gridspec(2, n_panels, height_ratios=[1.0, DEVIATION_HEIGHT_RATIO])
    axes = [fig.add_subplot(gs[0, 0])]
    for i in range(1, n_panels):
        axes.append(fig.add_subplot(gs[0, i], sharey=axes[0]))
    ax_dev = fig.add_subplot(gs[1, :], sharex=axes[0])

    last_index = n_panels - 1
    for i, (title, df, value_col) in enumerate(panel_specs):
        is_last = (i == last_index)
        plot_experimental_boxplots(axes[i], df, value_col, tree_order, year_order, show_legend=is_last)
        overlay_fem_points(axes[i], fem_df, value_specs, tree_order, show_legend=is_last, model_colors=model_colors)
        axes[i].set_title(title, fontsize=TITLE_FONTSIZE)

    # Deviation-bar fill legend (solid vs hatched = which FEM column) -
    # only needed when there is more than one FEM column to distinguish
    # (rotation: 2; displacement: 1) - shared with the upper panels
    # (attached to the LAST panel, drawn once, outside the axes), stacked
    # one slot below the "FEM value" legend overlay_fem_points() just
    # added there (which is add_artist'd so it survives this further
    # legend() call - see that function's docstring).
    if len(bar_specs) > 1:
        fill_handles = [
            Patch(facecolor="0.6", edgecolor=DEVIATION_BAR_EDGECOLOR, hatch=hatch, label=label)
            for _, hatch, label in bar_specs
        ]
        axes[last_index].legend(handles=fill_handles, title="Deviation bar fill", loc="upper left",
                                bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - 3 * LEGEND_Y_STEP), frameon=False,
                                fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    axes[0].set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)

    # Figure-level title stating which models are shown, whenever
    # MODEL_INCLUDE narrowed it - so a filtered figure is self-describing
    # even out of context (e.g. pasted alone into a write-up).
    if visible_models is not None:
        fig.suptitle(f"Models shown: {', '.join(visible_models)}", fontsize=TITLE_FONTSIZE)

    add_provenance_footer(fig, "plot_stiffness_comparison.py", PROVENANCE_NEWEST_INPUT_MTIME,
                          PROVENANCE_FILTERS_TEXT, n_trees=len(tree_order),
                          n_models=int(fem_df["model"].nunique()))

    panel_titles = "/".join(title for title, _, _ in panel_specs)
    n_points_per_panel = count_fem_points(fem_df, value_specs, tree_order)
    print(f"{print_label}: {n_points_per_panel} FEM point(s) per stiffness panel "
          f"({panel_titles} each) - identical in the 'all' and 'trimmed' files, "
          f"since both are saved from this same figure/axes without rebuilding them.")

    save_stiffness_figure_variants(fig, ax_dev, deviation_table, bar_specs, model_colors,
                                   tree_order, quantity, title_prefix, right_margin=right_margin,
                                   visible_models=visible_models)


def build_rotation_figure(experimental: pd.DataFrame, fem_df_full: pd.DataFrame, fem_df_drawing: pd.DataFrame,
                          tree_order: list, year_order: list, visible_models: list):
    """
    Build and save the rotation-stiffness comparison figures via
    build_stiffness_figure(): the full one (3 subplots - "cam", "nocam",
    "combined" - each with grey experimental boxplots and the FEM
    rotation-based results overlaid as colored points, PLUS a grouped-bar
    deviation subplot below) written as stiffness_comparison_rotation_all/
    _trimmed.png, AND a simplified "combined only" one (VARIANT 1 - the
    "cam"/"nocam" panels dropped, same deviation subplot) written as
    stiffness_comparison_rotation_combined_all/_trimmed.png - both pairs
    from the SAME deviation_table, computed once below.

    fem_df_full (every model) feeds compute_deviation_table() - the CSV is
    never affected by MODEL_INCLUDE. fem_df_drawing (MODEL_INCLUDE already
    applied, or identical to fem_df_full when that filter is off) is what
    actually gets drawn. visible_models is None (no filter) or the
    resolved list of models being drawn, passed straight to
    build_stiffness_figure() for its legend/title/dodge-width behaviour.
    """
    assert_gravity_variant(fem_df_full, "both", "build_rotation_figure (fem_df_full)")  # feeds the CSV - always both variants
    assert_gravity_variant(fem_df_drawing, GRAVITY_VARIANT, "build_rotation_figure (fem_df_drawing)")

    combined = pair_and_median(experimental, "slope_root_stiff_cam", "slope_root_stiff_nocam")

    value_specs = [
        ("slope_rotation_base", "o", "rotsum base"),
        ("slope_rotation_base2", "^", "rotsum base2"),
    ]
    bar_specs = [
        ("slope_rotation_base", DEVIATION_BASE_HATCH, "rotsum base"),
        ("slope_rotation_base2", DEVIATION_BASE2_HATCH, "rotsum base2"),
    ]
    model_colors = assign_model_colors(fem_df_full["model"].unique())
    ylabel = "Rotation stiffness (moment / rotation)"

    deviation_table = compute_deviation_table(
        experimental, combined, fem_df_full,
        sensor_columns={"cam": "slope_root_stiff_cam",
                        "nocam": "slope_root_stiff_nocam",
                        "combined": None},
        fem_columns=["slope_rotation_base", "slope_rotation_base2"],
    )

    build_stiffness_figure(
        "rotation",
        panel_specs=[("cam", experimental, "slope_root_stiff_cam"),
                     ("nocam", experimental, "slope_root_stiff_nocam"),
                     ("combined", combined, "value")],
        fem_df=fem_df_drawing, value_specs=value_specs, bar_specs=bar_specs, model_colors=model_colors,
        tree_order=tree_order, year_order=year_order, ylabel=ylabel, deviation_table=deviation_table,
        figure_width=3 * STIFFNESS_PANEL_WIDTH_INCHES,
        title_prefix="Relative deviation (combined sensor)", print_label="Rotation figure",
        visible_models=visible_models,
    )

    # VARIANT 1 - simplified "combined only" figure for the write-up: same
    # deviation_table (so its bars are guaranteed identical to the full
    # figure's), just the cam/nocam panels dropped.
    build_stiffness_figure(
        "rotation_combined",
        panel_specs=[("combined", combined, "value")],
        fem_df=fem_df_drawing, value_specs=value_specs, bar_specs=bar_specs, model_colors=model_colors,
        tree_order=tree_order, year_order=year_order, ylabel=ylabel, deviation_table=deviation_table,
        figure_width=STIFFNESS_COMBINED_WIDTH_INCHES, right_margin=STIFFNESS_COMBINED_RIGHT_MARGIN,
        title_prefix="Relative deviation (combined sensor)", print_label="Rotation figure (combined only)",
        visible_models=visible_models,
    )

    return combined, deviation_table


def build_displacement_figure(experimental: pd.DataFrame, fem_df_full: pd.DataFrame, fem_df_drawing: pd.DataFrame,
                              tree_order: list, year_order: list, visible_models: list):
    """
    Build and save the displacement-stiffness comparison figures via
    build_stiffness_figure(): the full one (3 subplots - "Pt3", "Pt4",
    "combined" - PLUS a grouped-bar deviation subplot below) written as
    stiffness_comparison_displacement_all/_trimmed.png, AND a simplified
    "combined only" one (VARIANT 1 - the "Pt3"/"Pt4" panels dropped, same
    deviation subplot) written as
    stiffness_comparison_displacement_combined_all/_trimmed.png - both
    pairs from the SAME deviation_table, computed once below.

    See build_rotation_figure()'s own docstring for fem_df_full vs.
    fem_df_drawing vs. visible_models.
    """
    assert_gravity_variant(fem_df_full, "both", "build_displacement_figure (fem_df_full)")  # feeds the CSV - always both variants
    assert_gravity_variant(fem_df_drawing, GRAVITY_VARIANT, "build_displacement_figure (fem_df_drawing)")

    combined = pair_and_median(experimental, "slope_Pt3", "slope_Pt4")

    value_specs = [("slope_displacement", "o", "displacement")]
    bar_specs = [("slope_displacement", DEVIATION_BASE_HATCH, "displacement")]
    model_colors = assign_model_colors(fem_df_full["model"].unique())
    ylabel = "Displacement stiffness (moment / displacement)"

    deviation_table = compute_deviation_table(
        experimental, combined, fem_df_full,
        sensor_columns={"Pt3": "slope_Pt3",
                        "Pt4": "slope_Pt4",
                        "combined": None},
        fem_columns=["slope_displacement"],
    )

    build_stiffness_figure(
        "displacement",
        panel_specs=[("Pt3", experimental, "slope_Pt3"),
                     ("Pt4", experimental, "slope_Pt4"),
                     ("combined", combined, "value")],
        fem_df=fem_df_drawing, value_specs=value_specs, bar_specs=bar_specs, model_colors=model_colors,
        tree_order=tree_order, year_order=year_order, ylabel=ylabel, deviation_table=deviation_table,
        figure_width=3 * STIFFNESS_PANEL_WIDTH_INCHES,
        title_prefix="Relative deviation (combined sensor)", print_label="Displacement figure",
        visible_models=visible_models,
    )

    # VARIANT 1 - simplified "combined only" figure for the write-up: same
    # deviation_table (so its bars are guaranteed identical to the full
    # figure's), just the Pt3/Pt4 panels dropped.
    build_stiffness_figure(
        "displacement_combined",
        panel_specs=[("combined", combined, "value")],
        fem_df=fem_df_drawing, value_specs=value_specs, bar_specs=bar_specs, model_colors=model_colors,
        tree_order=tree_order, year_order=year_order, ylabel=ylabel, deviation_table=deviation_table,
        figure_width=STIFFNESS_COMBINED_WIDTH_INCHES, right_margin=STIFFNESS_COMBINED_RIGHT_MARGIN,
        title_prefix="Relative deviation (combined sensor)", print_label="Displacement figure (combined only)",
        visible_models=visible_models,
    )

    return combined, deviation_table


def print_summary(experimental: pd.DataFrame, fem_df: pd.DataFrame):
    """
    Print, per tree/year, how many experimental pulls contributed data,
    and list the FEM tree/scenario combinations that were plotted.
    """
    print("\nExperimental pulls per tree/year (rows with at least one probe reading):")
    pull_counts = experimental.groupby(["tree", "year"])[["measurement"]].count()
    print(pull_counts.rename(columns={"measurement": "n_rows"}).to_string())

    print(f"\nFEM tree/model combinations plotted ({len(fem_df)} total):")
    print(fem_df[["tree", "scenario", "model"]].to_string(index=False))


def main():
    global TREES_SUFFIX, PROVENANCE_NEWEST_INPUT_MTIME, PROVENANCE_FILTERS_TEXT

    print_active_configuration("plot_stiffness_comparison.py", CONFIG_CONSTANT_NAMES)

    # --- Part 0: MODEL_INCLUDE / MODEL_INCLUDE_SUFFIX sanity check --------------
    # Writing filtered figures under the same filenames as the full ones
    # would silently and irreversibly overwrite them - refuse outright.
    if MODEL_INCLUDE is not None and not MODEL_INCLUDE_SUFFIX:
        print("MODEL_INCLUDE is set but MODEL_INCLUDE_SUFFIX is empty - refusing to run: "
              "writing filtered figures under the SAME filenames as the full ones would "
              "silently and irreversibly overwrite them. Set MODEL_INCLUDE_SUFFIX to a "
              "non-empty string and re-run.")
        return

    # --- Part 0b: freshness guard --------------------------------------------
    # Would have caught the Sep 2026 incident: raw ANSYS files for S04/S08/
    # S09/S10 finished AFTER slopes_combined.csv was last built from them.
    PROVENANCE_NEWEST_INPUT_MTIME = check_input_freshness(
        "plot_stiffness_comparison.py", ANSYS_RAW_FOLDER, ANSYS_RAW_PATTERN_FOR_FRESHNESS,
        {"FEM_CSV_PATH": FEM_CSV_PATH}, ABORT_ON_STALE_INPUT,
    )

    # --- Part 1: load and filter experimental data --------------------------
    experimental = load_and_filter_experimental(EXPERIMENTAL_CSV_PATH)

    if len(experimental) == 0:
        print("No experimental rows matched the configuration - check DAYS_TO_PROCESS, "
              "MEASUREMENT, TYPE_FILTER, REDUCTION_LEVEL_FILTER and TREE_FILTER.")
        return

    # --- Part 3: load FEM data ------------------------------------------------
    # UNFILTERED - always feeds stiffness_deviation.csv. Neither MODEL_INCLUDE
    # nor GRAVITY_VARIANT ever narrows fem_df itself - both are DRAW-ONLY
    # filters (see each one's own comment) - specifically so this script's
    # CSV output stays ONE complete, config-independent file that
    # method_deviation_overview.py can always read regardless of what THIS
    # script's own GRAVITY_VARIANT happens to be set to (an earlier version
    # of this stage's own work briefly broke that guarantee - see
    # GRAVITY_VARIANT's own comment for the postmortem).
    fem_df = load_fem_data(FEM_CSV_PATH)

    # --- TREES: resolve the selection against what is actually in fem_df,
    # print the SELECTION REPORT, then narrow fem_df to it. This runs
    # BEFORE MODEL_INCLUDE so a tree request is checked against the full
    # data, not an already model-narrowed slice. (SCENARIOS retired in
    # Stage 8 - see STRICT_SELECTION's own comment.) ---
    tree_sel = resolve_tree_selection(TREES, sorted(fem_df["tree"].unique()))
    any_missing = print_selection_report("plot_stiffness_comparison.py", tree_sel)
    if STRICT_SELECTION and any_missing:
        raise SystemExit(
            "plot_stiffness_comparison.py: STRICT_SELECTION=True and the SELECTION REPORT "
            "above lists a requested tree with no data - see 'requested but MISSING'. "
            "Set STRICT_SELECTION=False to plot the partial selection anyway."
        )
    fem_df = apply_tree_selection(fem_df, tree_sel)
    TREES_SUFFIX = build_filter_suffix(tree_sel)
    PROVENANCE_FILTERS_TEXT = (f"GRAVITY_VARIANT={GRAVITY_VARIANT}  MODEL_INCLUDE={MODEL_INCLUDE}  "
                               f"TREES={TREES}")

    # --- MODEL_INCLUDE: resolve which models are DRAWN (never affects fem_df
    # itself, which stays the full table for compute_deviation_table/CSV) ---
    all_models = sorted(fem_df["model"].unique())
    visible_models_list = resolve_model_filter(all_models, MODEL_INCLUDE)
    if MODEL_INCLUDE is not None and len(visible_models_list) == 0:
        print("MODEL_INCLUDE matched no model at all - nothing to draw. No figures written.")
        return
    fem_df_drawing = fem_df if MODEL_INCLUDE is None else fem_df[fem_df["model"].isin(visible_models_list)]
    visible_models = None if MODEL_INCLUDE is None else visible_models_list

    # --- THE GRAVITY-VARIANT STRUCTURAL LOCK - applied to DRAWING ONLY -----
    # Mirrors MODEL_INCLUDE exactly (same reasoning, same guarantee): this
    # narrows fem_df_drawing further, never fem_df/the CSV. See
    # GRAVITY_VARIANT's own comment for why this must never reach the CSV.
    fem_df_drawing_before_gravity = fem_df_drawing
    fem_df_drawing = filter_by_gravity_variant(fem_df_drawing, GRAVITY_VARIANT)
    print(f"GRAVITY_VARIANT={GRAVITY_VARIANT!r} - fem_df_drawing narrowed from "
          f"{len(fem_df_drawing_before_gravity)} to {len(fem_df_drawing)} row(s); "
          f"models drawn: {sorted(fem_df_drawing['model'].unique())}")

    # --- Shared axis categories -----------------------------------------------
    # With ONLY_TREES_WITH_FEM the x-axis is restricted to trees that have
    # at least one FEM result, so the figure is not mostly empty slots
    # while the ANSYS run plan is still incomplete. A tree that has FEM
    # results but none of them survive the MODEL_INCLUDE filter is ALSO
    # dropped from the drawing x-axis (reported separately below) - it is
    # never dropped from fem_df/the CSV.
    experimental_trees = set(experimental["tree"].unique())
    fem_trees_full = set(fem_df["tree"].unique())
    fem_trees_drawing = set(fem_df_drawing["tree"].unique())

    if ONLY_TREES_WITH_FEM:
        tree_order = sorted(experimental_trees & fem_trees_drawing)
        dropped = sorted(experimental_trees - fem_trees_full)
        if dropped:
            print(f"ONLY_TREES_WITH_FEM is on - {len(dropped)} experimental tree(s) "
                  f"with no FEM result at all are not drawn: {dropped}")
        dropped_by_filter = sorted((experimental_trees & fem_trees_full) - fem_trees_drawing)
        if dropped_by_filter:
            print(f"MODEL_INCLUDE dropped {len(dropped_by_filter)} tree(s) with no visible "
                  f"model left: {dropped_by_filter}")
    else:
        tree_order = sorted(experimental_trees | fem_trees_drawing)

    if len(tree_order) == 0:
        print("No tree appears in BOTH tables (after any MODEL_INCLUDE filtering) - "
              "nothing to plot. Check that the experimental IDs were converted to the "
              "canonical form, and that MODEL_INCLUDE matches real models.")
        return

    year_order = sorted(experimental["year"].unique())
    print(f"Trees on x-axis ({len(tree_order)}): {tree_order}")
    print(f"Years (boxplot groups) ({len(year_order)}): {year_order}")

    # --- Parts 4-5: build both figures ----------------------------------------
    _, rotation_deviation = build_rotation_figure(experimental, fem_df, fem_df_drawing, tree_order,
                                                  year_order, visible_models)
    _, displacement_deviation = build_displacement_figure(experimental, fem_df, fem_df_drawing, tree_order,
                                                           year_order, visible_models)

    # --- Part 6: deviation table, printed and saved ----------------------------
    deviation_table = pd.concat([rotation_deviation, displacement_deviation],
                                ignore_index=True)

    print("\nRelative deviation, FEM vs experimental median "
          "(positive = model stiffer than measurement):")
    print(deviation_table.to_string(index=False))

    deviation_csv_path = csv_filter_suffix(DEVIATION_CSV_PATH)
    deviation_table.to_csv(deviation_csv_path, index=False)
    print(f"\nDeviation table saved to: {deviation_csv_path}")

    # --- Part 7: console summary ----------------------------------------------
    print_summary(experimental, fem_df)

    # --- Part 8: orphaned filtered outputs (list only, never touched) --------
    report_orphaned_filtered_outputs(OUT_FOLDER, ("stiffness_comparison_", "stiffness_deviation"),
                                     "plot_stiffness_comparison.py")

    plt.show()


if __name__ == "__main__":
    main()
