"""
plot_modal_harmonic.py
------------------------
Produces, per tree:
  1. A modal-analysis figure showing which mode/frequency dominates the
     model's response in each direction (x, y, z), one bar per mode,
     colored by reconstruction model.
  2. A harmonic-response figure comparing the FEM frequency-response curve
     (colored, by model) against the experimentally measured fundamental
     frequency - both against a signed response (the default) and, when
     HARMONIC_USE_ABS is on, an additional absolute-value variant
     (harmonic_comparison_{tree}_{variant}_abs.png) where the resonance
     peak reads as a single clear maximum. Overlaid on a secondary y-axis:
     each mode's TIME_kg/TIME_incl_rots (Total translational Effective Mass,
     with/without rotational DOF) as points (no connecting line, colored by
     model, marker shape distinguishing the two columns) - so it is visible
     whether the harmonic peak coincides with the mode carrying the most
     effective mass. See load_modal_data()'s docstring for where these two
     columns come from (read from the END of the row, not positionally -
     the 19-name MODAL_COLUMN_NAMES list does not cover them) and
     print_time_mass_sanity_check()/print_frequency_range_check() for the
     printed checks run on them before any figure is drawn.
Plus one figure across all trees:
  3. A frequency-summary figure (frequency_summary_{elasto|median}.png):
     for every tree and model, the harmonic peak frequency and the modal
     frequency of the dominant EFF MASS direction, against the
     experimental frequency - with a printed and CSV-written table of
     relative deviations (frequency_deviation.csv).

Plus two simplified write-up variants, additions alongside the figures
above (never a replacement - the full ones keep being generated exactly
as before):
  - harmonic_simple_{tree}_{elasto|median}.png (VARIANT 2): a single-axis
    reduction of figure 2 - Top_ux only, always absolute value, no
    Half_ux/Quart_ux, no modal points, no effective-mass or mass-fraction
    axis. See plot_harmonic_simple_figure(), which shares its curve-
    drawing and x-range logic with plot_harmonic_figure() via
    plot_model_response_lines()/compute_harmonic_frequency_range() rather
    than reimplementing them.
  - frequency_summary_harmonic_{elasto|median}.png (VARIANT 3): figure 3
    with the modal_peak quantity left out everywhere (points, deviation
    bars, legend entries) - build_frequency_summary_figure()'s existing
    SUMMARY_SHOW_MODAL codepath, now exposed as a `show_modal` argument
    (default: the module-level SUMMARY_SHOW_MODAL, so every EXISTING
    caller is unaffected) and called here with show_modal=False
    explicitly. Its returned deviation table is discarded, not folded
    into frequency_deviation.csv - that CSV always contains both
    quantities exactly as before, for every variant.

In the modal and harmonic figures, the experimental fundamental frequency
is overlaid as a vertical band (25th-75th percentile) plus a median line -
grey, per the same grey=experiment / color=FEM convention as
plot_stiffness_comparison.py - drawn under the FEM data so it never hides
a curve or bar. Every legend is drawn outside the axes, to the right.

Also, report-only for now (no figure yet): mass-normalized effective-mass
quantities per mode (effmass_sum_kg, effmass_norm_kg, f_total, f_captured -
see add_mass_normalized_effmass_columns()), needed because the geometries
being compared differ in total mass by up to ~40%, which confounds "this
mode activates a lot of mass" with "this model is heavy". M_total_kg (the
model mass) is located in each modal file by find_total_mass_kg(), which
SEARCHES the first few lines for "TOTAL MASS of MODEL =" rather than
assuming a fixed line number. The per-mode table, including M_total_kg, is
written to modal_effmass.csv; six checks (extraction completeness,
per-direction coverage, top-5 modes by f_total, sum-vs-norm ranking
agreement, modes extracted, top-mode frequency vs. harmonic peak) are
printed - see print_modal_extraction_completeness() and its neighbours.

How to use it:
  - Edit the values in the "CONFIGURATION" section below.
  - Run the script from a terminal in VS Code:  python plot_modal_harmonic.py
  - Check the generated PNG files and the console summary.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folder holding the ANSYS output files (out_modal_*, out_harmonic_*, ...).
# This folder did not move - it is read from the same scripts\ansys folder
# every other script in this group reads it from.
INPUT_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\ansys"
MODAL_FILE_PATTERN = "out_modal_bk*_scenario*.txt"
HARMONIC_FILE_PATTERN = "out_harmonic_bk*_scenario*.txt"

# Name of the description file (one row per tree), inside INPUT_FOLDER. It is
# deliberately read from the SAME folder that ANSYS itself reads, so the two
# can never drift apart. Needed here ONLY to translate the position number in
# the filename into a real tree number - same approach as
# load_static_scenarios.py's load_ansys_order_to_tree_id().
GEOMETRY_FILE_NAME = "Popis_Babice_VSE_09012025_ans_lev6_0_azimcor.txt"
GEOMETRY_FILE_PATH = Path(INPUT_FOLDER) / GEOMETRY_FILE_NAME

# Canonical tree ID used across the whole reconstruction project:
# TREE_ID_PREFIX + zero-padded tree number, e.g. "B21_S" + "07" -> "B21_S07".
TREE_ID_PREFIX = "B21_S"

# Positions (0-based) in the description file. Position 0 is the row order,
# which equals the number after "bk" in the ANSYS filenames; position 2 is
# the tree's own number.
GEOMETRY_COLUMN_ANSYS_ORDER = 0
GEOMETRY_COLUMN_TREE_NUMBER = 2

# What each ANSYS scenario number actually is, under the settled model
# selection (handoff_2026-09-10_part2.md section 1):
#   1   adtree_r5mm     = AdTree r5mm, regression-perorder, -wlen,
#                         seg50-500-k80, calibrated against that tree's
#                         selected AdQSM variant
#   2   treeqsm_pd002   = TreeQSM pd2min 0.002, mo9, sr010, ri1, Simplified
#   999 adtree_r5_pilot = the cost-model pilot run on B21_S07, kept for
#                         reference (unchanged by the relabel)
# Scenario numbers 3 and 4 are deliberately left unassigned - no model is
# using them right now. ADD A LINE HERE whenever a new scenario is run -
# nothing else in the script needs changing. A scenario number that is
# missing from this dict still gets loaded, but is labelled "scenario<N>"
# and reported as a warning. Identical to the dict in
# load_static_scenarios.py.
SCENARIO_MODELS = {
        1: "adtree_r5mm",
        2: "treeqsm_pd002",
        3: "adtree_r5mm_no acel",
        4: "treeqsm_pd002_no_acel",
        999: "adtree_r5_pilot",   # high cylinder count run kept from the cost-model session
}

# Folder holding the input data and this script's output files.
DATA_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\data"
EXPERIMENTAL_CSV_PATH = str(Path(DATA_FOLDER) / "dynatree_data_merge_0926.csv")

DAYS_TO_PROCESS = ["2021-03-22", "2022-04-05"]   # configurable, same convention as before

# None = include ALL measurements (M01, M02, M03, M04...) that have a
# non-null first_frequency value. Unlike the static comparison (which used
# only M01), the fundamental frequency is typically only computed from the
# dynamic "pull and release" tests (M02+), so restricting to M01 here would
# likely give no data. Set to a specific measurement or list only if you
# want to narrow it further.
FREQUENCY_MEASUREMENT_FILTER = None
TYPE_FILTER = "normal"
REDUCTION_LEVEL_FILTER = "none"
TREE_FILTER = None

OUTPUT_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\out"

# Where the frequency-summary figure's deviation table is written (CHANGE 4).
# One combined file for both the elasto and median variants (see main()),
# with a "variant" column to tell the two apart.
FREQUENCY_DEVIATION_CSV_PATH = str(Path(OUTPUT_DIR) / "frequency_deviation.csv")

# Per-mode table with the mass-normalized effective-mass columns (see
# add_mass_normalized_effmass_columns()), one row per mode, including
# M_total_kg so the normalization can be rechecked later.
MODAL_EFFMASS_CSV_PATH = str(Path(OUTPUT_DIR) / "modal_effmass.csv")

# True  = also generate an absolute-value variant of every harmonic figure
#         (harmonic_comparison_{tree}_{variant}_abs.png), so the resonance
#         peak reads as a single clear maximum instead of a signed curve
#         that crosses zero.
# False = only generate the existing signed harmonic figures.
# The signed figures are ALWAYS generated regardless of this flag - this
# only controls whether the extra abs variant is ALSO written.
HARMONIC_USE_ABS = True

# True  = the frequency-summary figure (build_frequency_summary_figure())
#         shows the modal quantities (modal_peak points in the upper panel,
#         modal_peak bars in the deviation panel, and their legend entries)
#         alongside the harmonic ones.
# False = harmonic-only: every modal_peak point/bar/legend entry is left
#         out entirely - not drawn empty, not hidden behind a blank legend
#         swatch. frequency_deviation.csv is UNAFFECTED either way - it
#         always contains both quantities; this flag only controls what
#         build_frequency_summary_figure() draws.
# NOTE: whether the final figure keeps the modal quantities at all, or
# shows harmonic only, is not yet decided - this flag exists so that
# decision can be made later by flipping one value, not by deleting code.
SUMMARY_SHOW_MODAL = True

# None = every model is DRAWN. A list = only these models are drawn -
# matched against the "model" column (not the scenario number, so this
# stays readable if scenario numbers are ever reassigned), tolerant of the
# space/underscore inconsistency in the no-acel suffix (see
# normalize_model_name_for_matching()) but never rewriting a displayed
# name - that inconsistency is being fixed in ANSYS, not here.
#
# This is a FILTER ON WHAT IS DRAWN ONLY, independent of SUMMARY_SHOW_MODAL
# above: modal_df/harmonic_df/modal_effmass.csv/frequency_deviation.csv
# always keep every row for every model, regardless of this setting.
MODEL_INCLUDE = None

# Appended to every figure filename whenever MODEL_INCLUDE is set (e.g.
# modal_comparison_B21_S04_elasto_noacel.png), so a filtered figure can
# never be confused with, or silently overwrite, the full one. main()
# refuses to run if MODEL_INCLUDE is set and this is empty.
MODEL_INCLUDE_SUFFIX = "noacel"

# X-axis crop, in Hz, applied to every figure whose x-axis is frequency:
# both harmonic figures (signed and abs) and both modal figures. Set to
# None for auto (matplotlib's own autoscale - the previous behaviour).
# VIEW ONLY: this never filters/drops a data point from modal_df,
# harmonic_df, frequency_deviation.csv, or any printed number - it only
# changes what is visible on these specific figures. A point that falls
# outside this range is still in every table, just off the edge of the
# plot - see print_outside_axis_limits_report(), which is called in main()
# specifically to make that fact visible instead of silent.
# NOT applied to the frequency_summary figure - there, frequency is on the
# Y axis, not X, so this constant does not apply to it (see main()).
FREQUENCY_AXIS_LIMITS = (0.0,1)   # Hz, x axis; set to None for auto

# Same list as in plot_stiffness_comparison.py, for visual consistency.
PASTEL_SCENARIO_COLORS = [
    "#94BBEA", "#F5B399", "#8DD7BC", "#F6D07F",
    "#F3BDD1", "#7FC17F", "#A49CD3", "#F1A4A3",
]

# --- Style config -----------------------------------------------------------
TITLE_FONTSIZE = 16
AXIS_LABEL_FONTSIZE = 12
TICK_LABEL_FONTSIZE = 10
LEGEND_FONTSIZE = 10
LEGEND_TITLE_FONTSIZE = 10
LINE_WIDTH = 1.5          # FEM harmonic curves
BAR_EDGE_LINEWIDTH = 0.0  # modal EFF MASS bars (no edge by default)
BOX_EDGE_LINEWIDTH = 1.0  # experimental boxplots (still used by plot_stiffness_comparison.py's convention; kept for any future boxplot use)
MARKER_SIZE = 8            # FEM/summary point size

# Legends, moved outside the axes (CHANGE 2). bbox_to_anchor is in axes
# fraction coordinates: LEGEND_BBOX_X > 1 sits just right of the axes edge,
# LEGEND_Y_TOP anchors the topmost stacked legend, and each further legend
# on the same axes is offset down by LEGEND_Y_STEP so they stack instead of
# overlapping. RIGHT_MARGIN is how much of the figure width (0-1) is left
# for the plot itself - the rest is blank margin reserved for the legends.
LEGEND_BBOX_X = 1.02
LEGEND_Y_TOP = 1.0
LEGEND_Y_STEP = 0.32
RIGHT_MARGIN = 0.78

# Experimental frequency band (CHANGE 1): a shaded axvspan from the 25th to
# the 75th percentile, plus a solid axvline at the median, drawn under the
# FEM data (BAND_ZORDER < the FEM curves'/bars' zorder of 3).
BAND_ALPHA = 0.18          # shaded-band fill transparency
BAND_ZORDER = 1            # band fill zorder (below FEM data at zorder 3)
BAND_MIN_N = 3             # minimum measured values required to draw the shaded band; below this, only the median line is drawn
MEDIAN_LINE_WIDTH = 1.5    # width of the median axvline

# Frequency-summary figure (CHANGE 4).
SUMMARY_MARKER_HARMONIC = "o"       # harmonic peak frequency (max |Top_ux|)
SUMMARY_MARKER_MODAL = "^"          # modal frequency of the dominant EFF MASS direction
SUMMARY_MARKER_EXPERIMENTAL = "s"   # experimental median point
SUMMARY_EXPERIMENTAL_COLOR = "0.35"  # grey, per the grey=experiment convention
SUMMARY_ERRORBAR_LINEWIDTH = 1.5    # thin vertical line spanning the experimental 25th-75th percentile
SUMMARY_GROUP_WIDTH = 0.6           # horizontal space (of 1.0 per tree) the side-by-side model points occupy
# This figure's "Quantity" legend has longer labels ("Modal peak (max EFF
# MASS)") than the other figures' model-name-only legends, and with few
# trees the figure itself is narrow - RIGHT_MARGIN's fraction of a narrow
# figure is not enough absolute width for that text, so this figure gets
# its own, wider margin and a taller minimum width.
SUMMARY_MIN_WIDTH_INCHES = 11.0
SUMMARY_RIGHT_MARGIN = 0.68

# Deviation panel (CHANGE 2 - replaces the old multi-line text footnote,
# which ran to 24 lines and overran the axes): a grouped-bar subplot below
# the frequency-summary panel, sharing its x-axis (trees). One bar per
# model per quantity, dodged the same way as that panel's points
# (offset_by_model, reused) plus a sub-offset for the 1-2 quantities within
# each model's slot (sub_offsets, computed the same way via
# compute_dodge_offsets() so a single quantity - SUMMARY_SHOW_MODAL=False -
# centers on the model's slot instead of leaving a gap where its sibling
# bar would have been).
SUMMARY_DEVIATION_HEIGHT_RATIO = 0.4   # deviation panel height, as a fraction of the main panel's - named per this task's own instruction
SUMMARY_DEVIATION_BAR_WIDTH = 0.09     # width of ONE quantity's bar within a model's dodge slot
SUMMARY_HARMONIC_HATCH = "///"         # harmonic_peak bars: hatched (quantity distinguished by hatch, not colour - colour already means "model")
SUMMARY_MODAL_HATCH = None             # modal_peak bars: solid (no hatch)
SUMMARY_DEVIATION_BAR_EDGECOLOR = "black"  # bar edge - hatch patterns need a visible edge colour to render
SUMMARY_DEVIATION_BAR_EDGE_LINEWIDTH = 0.8
SUMMARY_ZERO_LINE_COLOR = "0.3"        # horizontal line at 0% deviation, so the sign is readable at a glance
SUMMARY_ZERO_LINE_WIDTH = 1.0

# Modal effective-mass overlay on the harmonic figure (secondary y-axis) -
# see plot_harmonic_figure()'s own docstring. Points only, no connecting
# line: colour = model (reuses model_colors, the same mapping as everywhere
# else), marker SHAPE distinguishes TIME_kg from TIME_incl_rots.
MODAL_MASS_MARKER_TIME_KG = "D"            # TIME_kg (Total Translational Effective Mass)
MODAL_MASS_MARKER_TIME_INCL_ROTS = "X"     # TIME_incl_rots (same, including rotational DOF)
MODAL_MASS_MARKER_SIZE = 70                # scatter marker area (points^2 - matplotlib's `s=`)
MODAL_MASS_LEGEND_MARKERSIZE = 8           # legend proxy marker size (points, not area)
MODAL_MASS_LEGEND_COLOR = "0.35"           # neutral grey for the shape-only legend (colour already means "model" elsewhere)
MODAL_MASS_ZORDER = 5                      # above the harmonic curve (zorder=3) so points are never hidden under it
MODAL_MASS_AXIS_LABEL = "Modal effective mass (kg)"
MODAL_MASS_AXIS_COLOR = "0.35"             # secondary y-axis spine/tick/label colour

# Mass-normalized effective-mass overlay (f_captured / f_total) on the
# harmonic figure - a THIRD y-axis, since these are dimensionless fractions
# (0-1) and cannot share the kg-scaled secondary axis above. Points only,
# same convention as the kg overlay: colour = model, marker SHAPE
# distinguishes f_captured from f_total. f_captured is the fraction of
# what the EXTRACTED modes captured, NOT of the whole model mass - see
# plot_harmonic_figure()'s own docstring for why f_total is plotted
# alongside it and the captured-mass percentage is annotated on the figure.
FRACTION_AXIS_OFFSET = 0.14                # how far right of the primary axis (axes-fraction units) this axis's spine sits, beyond the kg-axis's own spine at 1.0
FRACTION_AXIS_LABEL = "Fraction of mode-set mass (f_captured / f_total)"
FRACTION_AXIS_COLOR = "0.15"               # darker grey than MODAL_MASS_AXIS_COLOR (0.35), so the two secondary axes stay visually distinct
FRACTION_MARKER_F_CAPTURED = "P"           # f_captured: fraction of the EXTRACTED modes' mass
FRACTION_MARKER_F_TOTAL = "*"              # f_total: fraction of the WHOLE model mass
FRACTION_MARKER_SIZE = 70                  # scatter marker area (points^2), same convention as MODAL_MASS_MARKER_SIZE
FRACTION_LEGEND_MARKERSIZE = 8
FRACTION_LEGEND_COLOR = "0.15"             # neutral, matches FRACTION_AXIS_COLOR - colour already means "model" via the points themselves
FRACTION_ZORDER = 6                        # above the modal-mass kg points (zorder 5) so it is never hidden under them
FRACTION_ANNOTATION_FONTSIZE = 9

# plot_harmonic_figure()-only legend placement: pushed further out than
# LEGEND_BBOX_X/RIGHT_MARGIN (used by every OTHER figure, unchanged) because
# this figure's secondary y-axis now sits right where those legends used to
# start. Passed explicitly to add_experimental_frequency_band()'s bbox_x.
# Widened further (1.16 -> 1.34, 0.68 -> 0.60) to also clear the new third
# (fraction) axis added at FRACTION_AXIS_OFFSET beyond the kg axis.
HARMONIC_LEGEND_BBOX_X = 1.34
HARMONIC_RIGHT_MARGIN = 0.60

# Simplified single-axis harmonic figure (VARIANT 2: harmonic_simple_{tree}_
# {variant}.png) - Top_ux only, absolute values, no modal/mass-fraction
# axes, so it needs neither HARMONIC_LEGEND_BBOX_X above (sized for 3 extra
# stacked legends) nor a new figsize pair per subplot - it reuses the
# plain single-legend-stack LEGEND_BBOX_X used by plot_modal_figure().
# Model names here run long ("treeqsm_pd002_no_acel"), which the generic
# RIGHT_MARGIN=0.78 does not leave enough ABSOLUTE room for on a figure
# this narrow (legend text got clipped at the figure edge) - so this gets
# its own figsize/right-margin pair, wide enough for that text.
HARMONIC_SIMPLE_FIGSIZE = (10.0, 5.0)
HARMONIC_SIMPLE_RIGHT_MARGIN = 0.68

# Regular expressions used to pull tree/scenario numbers out of filenames.
MODAL_FILENAME_PATTERN = re.compile(r"out_modal_bk(\d+)_scenario(\d+)\.txt")
HARMONIC_FILENAME_PATTERN = re.compile(r"out_harmonic_bk(\d+)_scenario(\d+)\.txt")

# How many lines from the START of a modal file are searched for the
# "TOTAL MASS of MODEL = <value> <unit>" line (see find_total_mass_kg()).
# NOT assumed to be a fixed line number - the file format has already
# changed once this session (TIME_kg/TIME_incl_rots were added to the data
# rows after MODAL_COLUMN_NAMES was written), so a hardcoded line index for
# the mass line would break silently the next time a line is added/removed
# above it.
TOTAL_MASS_SEARCH_LINES = 5
TOTAL_MASS_PATTERN = re.compile(
    r"TOTAL\s+MASS\s+of\s+MODEL\s*=\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*([A-Za-z]*)",
    re.IGNORECASE,
)

# Fixed, positional column names for the modal data files (see PART 1 in
# the module docstring for why these can't be parsed from the header line).
MODAL_COLUMN_NAMES = [
    "FREQ", "PFACTx", "PFACTy", "PFACTz", "PFACTrotx", "PFACTroty", "PFACTrotz",
    "EFFMASSx_kg", "EFFMASSy_kg", "EFFMASSz_kg", "EFFMASSrotx_kgm2",
    "EFFMASSroty_kgm2", "EFFMASSrotz_kgm2", "TTME_kg", "TTME_incl_rots",
    "TTME_fsum_stem", "TTME_fsum_stem_incl_rots", "TTME_fsum_tree",
    "TTME_fsum_tree_incl_rots",
]

# TIME_kg / TIME_incl_rots: two columns the modal output files gained AFTER
# this positional list was written. NOT read positionally like the 19 names
# above - the file's own header repeats shorthand group names (PFACT/EFF
# MASS), so a column inserted anywhere in the MIDDLE would silently shift
# every name after it with no error. These two are instead read counting
# FROM THE END of the row (row[-N]), which stays correct even if another
# column is later inserted before them, at the front or in the middle - only
# a change in ROW LENGTH itself can move them, and that is exactly what the
# per-file value-count guard in load_modal_data() checks for and reports.
MODAL_COL_TIME_KG_FROM_END = 6
MODAL_COL_TIME_INCL_ROTS_FROM_END = 5

# Harmonic data column names, after the header's "TIME" is renamed to
# "frequency_hz" (that column is actually frequency-sweep steps in Hz).
HARMONIC_COLUMN_NAMES = [
    "frequency_hz", "Top_ux", "Top_uy", "Top_uz",
    "Half_ux", "Half_uy", "Half_uz", "Quart_ux", "Quart_uy", "Quart_uz",
]


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

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


# NOTE: no longer used by the loaders below (see load_ansys_order_to_tree_id()
# and tree_number_to_tree_id() instead) - kept only because nothing else in
# this task asked for it to be removed. This function embodies the bug this
# module used to have: it treats the number after "bk" in an ANSYS filename
# as the tree's own number, when it is actually the tree's ROW POSITION in
# the description file (see the module-level note near INPUT_FOLDER).
def bk_number_to_tree_id(bk_number: int) -> str:
    """
    Convert a plain tree number (e.g. 1, 4, 16) into the main dataset's
    tree ID format (e.g. "BK01", "BK04", "BK16") - :02d pads the number
    with a leading zero if it is a single digit.
    """
    return f"BK{bk_number:02d}"


def parse_czech_number(text: str) -> float:
    """
    Convert a number written in Czech format (comma as decimal separator,
    e.g. "6,43") into a Python float (6.43). The description file uses
    this format for every numeric value.

    Copied from load_static_scenarios.py, not imported - these scripts are
    deliberately standalone with no shared imports.
    """
    return float(text.replace(",", "."))


def tree_number_to_tree_id(tree_number: int) -> str:
    """
    Convert a plain tree number (e.g. 1, 7, 12) into the canonical project
    tree ID (e.g. "B21_S01", "B21_S07", "B21_S12"). :02d pads the number
    with a leading zero if it is a single digit.

    Copied from load_static_scenarios.py, not imported - these scripts are
    deliberately standalone with no shared imports.
    """
    return f"{TREE_ID_PREFIX}{tree_number:02d}"


def load_ansys_order_to_tree_id(geometry_file_path: str) -> dict:
    """
    Read the description file and return a dict that translates the number
    used in the ANSYS filenames (the tree's ROW POSITION, not its own
    number) into a canonical tree ID:

        {1: "B21_S01", 2: "B21_S04", 3: "B21_S07", ...}

    The file has 2 header/reference lines (column names, then column
    indices) which are skipped, followed by one tab-separated data row per
    tree. Numbers use a comma as the decimal separator, so they go through
    parse_czech_number() first.

    Why this exists: the ANSYS filenames carry the ROW POSITION, not the
    tree number. Without this lookup, "bk3" would be read as tree 3, which
    does not exist, and "bk8" would be read as tree 8, which does exist but
    is the wrong tree - a silent, invisible error.

    Copied from load_static_scenarios.py, not imported - these scripts are
    deliberately standalone with no shared imports.
    """
    path = Path(geometry_file_path)
    if not path.is_file():
        print(f"ERROR: description file not found: {path}")
        print("Text files present in that folder:")
        for candidate in sorted(path.parent.glob("*.txt")):
            print(f"  {candidate.name}")
        raise FileNotFoundError(path)

    with open(geometry_file_path, "r", encoding="utf-8") as f:
        lines = [line for line in f.read().splitlines() if line.strip() != ""]

    # Line 1 = column names, line 2 = column indices - both are only for
    # reference, so we skip them and start reading actual data from line 3.
    data_lines = lines[2:]

    order_to_tree_id = {}
    for line in data_lines:
        columns = line.split("\t")

        ansys_order = int(parse_czech_number(columns[GEOMETRY_COLUMN_ANSYS_ORDER]))
        tree_number = int(parse_czech_number(columns[GEOMETRY_COLUMN_TREE_NUMBER]))

        order_to_tree_id[ansys_order] = tree_number_to_tree_id(tree_number)

    return order_to_tree_id


def scenario_to_model_name(scenario_number: int) -> str:
    """
    Translate a scenario number into the reconstruction model name from
    SCENARIO_MODELS. An unlisted number is not an error - it gets a
    placeholder name so the row still appears in the output table, and the
    caller can print a warning.

    Copied from load_static_scenarios.py, not imported - these scripts are
    deliberately standalone with no shared imports.
    """
    return SCENARIO_MODELS.get(scenario_number, f"scenario{scenario_number}")


# WHY bk_id_to_canonical() exists: before this fix, the FEM side of this
# script produced tree IDs like "BK07" from a wrong assumption about tree
# numbering (treating the position after "bk" in an ANSYS filename as the
# tree's own number), and the experimental CSV's "tree" column independently
# also used "BK07" as its own label - a completely unrelated naming scheme
# that just happened to look the same. The two matched by coincidence, so
# the bug produced correct-looking plots and never raised an error. A file
# named bk8 would have been read as tree 8, while the tree actually at
# position 8 is tree 12 - silently pairing one tree's FEM model with a
# different tree's experimental measurements. Converting the experimental
# "tree" column to the canonical B21_S.. form (the same form the fixed FEM
# loaders now produce via load_ansys_order_to_tree_id()) makes the two sides
# use the same, correctly-derived tree identity instead of two independently
# wrong ones that coincidentally agreed.
def bk_id_to_canonical(bk_id: str) -> str:
    """
    Convert an experimental tree ID into the canonical project tree ID:
      "BK07" -> "B21_S07",  "BK12" -> "B21_S12"

    Only the digits are kept from the input, so the exact prefix in the
    experimental file does not matter. :02d pads a single digit with a
    leading zero.

    Copied from plot_stiffness_comparison.py, not imported - these scripts
    are deliberately standalone with no shared imports.
    """
    digits_only = "".join(character for character in str(bk_id) if character.isdigit())
    return f"{TREE_ID_PREFIX}{int(digits_only):02d}"


def darker(hex_color, factor: float = 0.7):
    """
    Return a darker version of hex_color by scaling its RGB components by
    factor (0-1 - smaller means darker). Used to give boxplots a subtle
    edge color that matches their own fill, instead of plain black.

    Example: darker("#94BBEA", 0.7) -> a slightly darker pastel blue.
    """
    r, g, b = mcolors.to_rgb(hex_color)
    return (r * factor, g * factor, b * factor)


def assign_scenario_colors(scenarios) -> dict:
    """
    Assign a distinct, consistent color to each FEM scenario number, using
    the fixed PASTEL_SCENARIO_COLORS list (cycling with modulo if there are
    more scenarios than colors), so the same scenario always gets the same
    color across all figures.
    """
    n_colors = len(PASTEL_SCENARIO_COLORS)
    return {
        scenario: PASTEL_SCENARIO_COLORS[i % n_colors]
        for i, scenario in enumerate(sorted(scenarios))
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


def add_model_include_suffix(path_str: str) -> str:
    """
    Insert MODEL_INCLUDE_SUFFIX before the file extension whenever
    MODEL_INCLUDE is set (e.g. ".../modal_comparison_B21_S04_elasto.png"
    -> "..._noacel.png") - see MODEL_INCLUDE's own comment for why. A no-op
    when MODEL_INCLUDE is None.
    """
    if MODEL_INCLUDE is None:
        return path_str
    root, ext = path_str.rsplit(".", 1)
    return f"{root}_{MODEL_INCLUDE_SUFFIX}.{ext}"


def add_light_grid(ax):
    """
    Draw the same light background grid used in plot_stiffness_comparison.py,
    placed behind the data (low zorder).
    """
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=0)


def get_year_grey_shades(years) -> dict:
    """
    Assign a light-to-dark grey shade to each year (earliest year =
    lightest), same convention as the year boxplots in
    plot_stiffness_comparison.py. Returns a dict {year: grey_shade_string}.
    """
    year_order = sorted(years)
    return {
        year: str(shade) for year, shade in
        zip(year_order, np.linspace(0.85, 0.45, len(year_order)))
    }


def add_experimental_frequency_band(ax, freq_df: pd.DataFrame, tree_id: str, year_colors: dict,
                                     show_legend: bool = True, legend_bbox_y: float = None,
                                     bbox_x: float = None):
    """
    Overlay the experimental fundamental-frequency distribution onto an
    existing ax whose x-axis is frequency, as a vertical band per year:
      - a light shaded band (axvspan) from the 25th to the 75th percentile
        of that year's measured first_frequency values
      - a solid vertical line (axvline) at the median
      - grey, per the existing convention that grey = experiment
      - drawn at BAND_ZORDER (below the FEM data's zorder of 3), so it
        never hides a curve or bar

    If a tree/year has fewer than BAND_MIN_N measured values, only the
    median line is drawn (no shaded band), and a note is printed explaining
    why - a band built from 1-2 points would misleadingly imply a spread
    that was never actually measured.

    show_legend controls whether this call also draws the "Experiment
    <year>" legend (stacked at legend_bbox_y, outside the axes) - callers
    that put several subplots side by side pass show_legend=False on every
    subplot except the last, to avoid one subplot's external legend landing
    on top of its neighbour (see plot_modal_figure()). When show_legend is
    True, the number of experimental values behind each year's band/line is
    also printed, so it is visible where a band is based on very few points.

    bbox_x overrides LEGEND_BBOX_X for this one legend - plot_harmonic_figure()
    passes HARMONIC_LEGEND_BBOX_X instead, since that figure's secondary
    (TIME_kg/TIME_incl_rots) y-axis now occupies the space right at
    LEGEND_BBOX_X and needs every stacked legend pushed further out to
    clear it. Every other caller leaves this at its default (None ->
    LEGEND_BBOX_X), unchanged.
    """
    if legend_bbox_y is None:
        legend_bbox_y = LEGEND_Y_TOP
    if bbox_x is None:
        bbox_x = LEGEND_BBOX_X

    tree_freq = freq_df[freq_df["tree"] == tree_id]
    years = sorted(tree_freq["year"].unique())

    if len(years) == 0:
        print(f"Note: no experimental frequency data for tree '{tree_id}' - "
              f"leaving this subplot without an experimental overlay.")
        return

    years_with_data = []
    for year in years:
        values = tree_freq[tree_freq["year"] == year]["first_frequency"].dropna()
        n = len(values)
        if n == 0:
            continue  # skip a year with no data
        years_with_data.append(year)

        color = year_colors[year]
        edge_color = darker(color)
        median = float(values.median())

        if n >= BAND_MIN_N:
            p25, p75 = float(values.quantile(0.25)), float(values.quantile(0.75))
            ax.axvspan(p25, p75, color=color, alpha=BAND_ALPHA, zorder=BAND_ZORDER)
        else:
            print(f"Note: tree '{tree_id}', year {year} has only {n} measured "
                  f"value(s) (fewer than BAND_MIN_N={BAND_MIN_N}) - the shaded "
                  f"band was skipped and only the median line is drawn.")

        ax.axvline(median, color=edge_color, linewidth=MEDIAN_LINE_WIDTH, zorder=BAND_ZORDER + 1)

        if show_legend:
            print(f"Experimental frequency band, tree '{tree_id}', year {year}: "
                  f"n={n} value(s), median={median:.3f} Hz.")

    if not show_legend or len(years_with_data) == 0:
        return

    # Small legend for the experimental band/line, separate from the
    # model-color legend (added by the caller). Moved outside the axes.
    year_handles = [
        Patch(facecolor=year_colors[year], edgecolor=darker(year_colors[year]),
              alpha=BAND_ALPHA, label=f"Experiment {year}")
        for year in years_with_data
    ]
    ax.legend(handles=year_handles, loc="upper left",
              bbox_to_anchor=(bbox_x, legend_bbox_y), frameon=False,
              fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)


def compute_dodge_offsets(n_groups: int, group_width: float) -> list:
    """
    Work out horizontal "dodge" offsets to place n_groups side-by-side
    sub-positions within one tree's x-axis slot, evenly spaced and centered
    on 0, spanning a total width of group_width.

    Example: compute_dodge_offsets(2, 0.6) -> [-0.3, 0.3]

    Copied from plot_stiffness_comparison.py, not imported - these scripts
    are deliberately standalone with no shared imports.
    """
    if n_groups == 1:
        return [0.0]
    step = group_width / (n_groups - 1)
    return [-group_width / 2 + i * step for i in range(n_groups)]


# --- PART 1: modal data ---------------------------------------------------

def find_total_mass_kg(lines: list, file_name: str):
    """
    Search the first TOTAL_MASS_SEARCH_LINES line(s) of a modal file's
    already-loaded lines for "TOTAL MASS of MODEL = <value> <unit>",
    wherever it actually sits - NOT assumed to be a fixed line number (see
    TOTAL_MASS_SEARCH_LINES's own comment for why).

    Prints the raw line it was found on, the parsed value, and the unit
    (or a note if no unit is present).

    Returns (total_mass_kg, mass_line_index):
      - mass_line_index is the 0-based position INTO `lines` where the
        line was found, so the caller can derive where the header/data
        rows start RELATIVE TO IT, instead of assuming a fixed offset -
        and so the mass line itself is never mistaken for a data row.
      - if the text is not found within the searched window, prints a
        WARNING naming the file and returns (None, None). The caller must
        leave the normalized columns blank for this file rather than
        substitute a value from another file or from the geometry.
    """
    for index, line in enumerate(lines[:TOTAL_MASS_SEARCH_LINES]):
        match = TOTAL_MASS_PATTERN.search(line)
        if match is None:
            continue
        value = float(match.group(1))
        unit = match.group(2) if match.group(2) else "(no unit given)"
        print(f"  '{file_name}' line {index + 1}: {line.strip()!r} -> "
              f"parsed mass = {value}, unit = {unit}")
        return value, index

    print(f"WARNING: '{file_name}' - 'TOTAL MASS of MODEL =' not found in the "
          f"first {TOTAL_MASS_SEARCH_LINES} line(s); leaving M_total_kg and every "
          f"mass-normalized column blank for this file.")
    return None, None


def load_modal_data(input_folder: str, file_pattern: str, order_to_tree_id: dict) -> pd.DataFrame:
    """
    Load every modal-analysis file matching file_pattern in input_folder
    into one combined DataFrame with columns: tree, ansys_no, scenario,
    model, source_file, M_total_kg, plus the 19 fixed MODAL_COLUMN_NAMES,
    plus TIME_kg and TIME_incl_rots (read from the END of the row - see
    MODAL_COL_TIME_KG_FROM_END/MODAL_COL_TIME_INCL_ROTS_FROM_END's own
    comment for why).

    order_to_tree_id translates the ANSYS filename's row POSITION into the
    canonical tree ID (see load_ansys_order_to_tree_id()) - a position that
    is not in that dict means the filename and the description file
    disagree, so the file is skipped rather than guessed.

    M_total_kg is the model mass located by find_total_mass_kg() -
    wherever it actually sits, not a fixed line number. The header line is
    assumed to immediately follow it (true of every file seen so far), so
    data rows start 2 lines after wherever the mass line was found. To
    make sure this never silently starts reading data one line too early
    or too late (e.g. by mistaking the mass line itself for a data row),
    this prints, per file, the data-row count under BOTH the old fixed
    "line 2" assumption and the new mass-line-relative one, and a final
    summary confirming whether every file's two counts agree.

    Also prints a per-file guard: the number of values found on each file's
    data rows, and a warning if that count is not the same across every
    file - a from-the-end index is only stable while the row length is
    stable, so a change in row length must be visible, not silently shift
    which values TIME_kg/TIME_incl_rots pick up.
    """
    matching_files = sorted(Path(input_folder).glob(file_pattern))
    print(f"Found {len(matching_files)} modal file(s) matching '{file_pattern}'")

    rows = []
    unknown_scenarios = set()
    row_length_by_file = {}   # filename -> the value-count of its data rows (see the guard below)
    data_row_counts_by_file = {}   # filename -> (old fixed-line-2 count, new mass-line-relative count)

    for file_path in matching_files:
        match = MODAL_FILENAME_PATTERN.match(file_path.name)
        if match is None:
            print(f"WARNING: '{file_path.name}' does not match the expected modal "
                  f"filename pattern - skipping this file.")
            continue

        ansys_order = int(match.group(1))
        scenario = int(match.group(2))

        tree = order_to_tree_id.get(ansys_order)
        if tree is None:
            print(f"WARNING: '{file_path.name}' refers to position bk{ansys_order}, "
                  f"which does not exist in the description file - skipping "
                  f"this file.")
            continue

        model = scenario_to_model_name(scenario)
        if scenario not in SCENARIO_MODELS:
            unknown_scenarios.add(scenario)

        with open(file_path, "r", encoding="utf-8") as f:
            lines = [line for line in f.read().splitlines() if line.strip() != ""]

        # Total mass: searched for, not assumed to be on a fixed line (see
        # find_total_mass_kg()). The header line is assumed to immediately
        # follow it - true of every file seen so far, and the only relative
        # position that has ever appeared - so data starts 2 lines after
        # wherever the mass line was actually found; the header itself
        # (several of whose names are ambiguous shorthand that repeats
        # between quantity groups) is skipped on purpose in favour of the
        # fixed MODAL_COLUMN_NAMES list.
        total_mass_kg, mass_line_index = find_total_mass_kg(lines, file_path.name)
        data_start_index = (mass_line_index + 2) if mass_line_index is not None else 2
        print(f"'{file_path.name}': tree={tree}, ansys_no={ansys_order}, "
              f"scenario={scenario}, model={model}, total mass = {total_mass_kg} kg")

        # Before/after data-row-count check: the OLD code always started
        # reading data at the fixed line 2 (0-based); the NEW code starts
        # at data_start_index, derived from wherever the mass line was
        # actually found. These must agree on every file seen so far (the
        # mass line has only ever appeared on line 0) - a mismatch would
        # mean the mass line, or something between it and the header, is
        # now being read as a data row (or vice versa).
        old_data_lines = lines[2:]
        new_data_lines = lines[data_start_index:]
        data_row_counts_by_file[file_path.name] = (len(old_data_lines), len(new_data_lines))
        print(f"  data rows: before (fixed line-2 start) = {len(old_data_lines)}, "
              f"after (mass-line-relative start) = {len(new_data_lines)} "
              f"{'(MATCH)' if len(old_data_lines) == len(new_data_lines) else '(MISMATCH - investigate!)'}")

        file_row_lengths = set()   # every value-count seen on this file's OWN data rows

        for row_number, line in enumerate(new_data_lines, start=data_start_index + 1):
            values = line.split()
            file_row_lengths.add(len(values))

            # Only the 19 NAMED columns require a minimum length - extra
            # trailing values (TIME_kg/TIME_incl_rots and anything else not
            # yet named) are fine and expected; a row is skipped only if it
            # is too SHORT to even fill the existing named columns, not
            # because it is longer than before (see this function's own
            # docstring / MODAL_COL_TIME_KG_FROM_END's comment for why
            # length is no longer required to match exactly).
            if len(values) < len(MODAL_COLUMN_NAMES):
                print(f"WARNING: '{file_path.name}' line {row_number} has "
                      f"{len(values)} value(s), expected at least {len(MODAL_COLUMN_NAMES)} "
                      f"- skipping this mode.")
                continue

            mode_values = dict(zip(MODAL_COLUMN_NAMES, [float(v) for v in values]))

            # TIME_kg / TIME_incl_rots, read from the end of the row (NOT
            # covered by the positional 19-name list above) - NaN if the
            # row is too short to contain them at all.
            mode_values["TIME_kg"] = (float(values[-MODAL_COL_TIME_KG_FROM_END])
                                       if len(values) >= MODAL_COL_TIME_KG_FROM_END else float("nan"))
            mode_values["TIME_incl_rots"] = (float(values[-MODAL_COL_TIME_INCL_ROTS_FROM_END])
                                              if len(values) >= MODAL_COL_TIME_INCL_ROTS_FROM_END else float("nan"))

            row = {"tree": tree, "ansys_no": ansys_order, "scenario": scenario,
                   "model": model, "source_file": file_path.name,
                   "M_total_kg": total_mass_kg if total_mass_kg is not None else float("nan")}
            row.update(mode_values)
            rows.append(row)

        if file_row_lengths:
            if len(file_row_lengths) > 1:
                print(f"WARNING: '{file_path.name}' has data rows of DIFFERENT lengths "
                      f"within the same file: {sorted(file_row_lengths)} - the from-the-end "
                      f"TIME_kg/TIME_incl_rots columns are not reliable for this file.")
            # One representative length per file for the cross-file guard
            # below - the max of what was seen (matches what a reader
            # checking "how many columns does this file have" would expect;
            # the within-file-inconsistency warning above already flagged
            # if that is not the whole story for this particular file).
            row_length_by_file[file_path.name] = max(file_row_lengths)

    if unknown_scenarios:
        print(f"WARNING: scenario number(s) {sorted(unknown_scenarios)} are not "
              f"listed in SCENARIO_MODELS (modal data). They were loaded with a "
              f"placeholder name - add them to SCENARIO_MODELS so the plots get "
              f"a proper label.")

    # --- STEP 1 guard: row length must be the same across every file -------
    # A from-the-end index (TIME_kg/TIME_incl_rots) is only stable while the
    # row length is stable - if one file has grown (or shrunk) relative to
    # the others, that must be visible here, not discovered later as a
    # silently wrong number on a plot.
    print(f"\nModal data row length per file (values per data row):")
    for name in sorted(row_length_by_file):
        print(f"  {name}: {row_length_by_file[name]} value(s)")
    distinct_lengths = sorted(set(row_length_by_file.values()))
    if len(distinct_lengths) > 1:
        by_length = {length: sorted(n for n, l in row_length_by_file.items() if l == length)
                     for length in distinct_lengths}
        print(f"WARNING: modal data row length is NOT consistent across files: {by_length}")
    elif distinct_lengths:
        n_values = distinct_lengths[0]
        n_unnamed = n_values - len(MODAL_COLUMN_NAMES)
        print(f"Every modal file has {n_values} value(s) per data row; MODAL_COLUMN_NAMES "
              f"names {len(MODAL_COLUMN_NAMES)} of them, leaving {n_unnamed} unnamed "
              f"(2 of which are read as TIME_kg/TIME_incl_rots via from-the-end indexing).")

    # --- mass-line search guard: confirm the search-based data start never
    # disagrees with the old fixed-line-2 assumption (see the per-file print
    # above for where this is computed) - printed again here as one final
    # summary so the answer to "did counts change?" does not require
    # scrolling back through every per-file line above.
    print(f"\nData-row count, before (fixed line-2 start) vs. after (mass-line-relative "
          f"start), per file:")
    mismatched_files = []
    for name in sorted(data_row_counts_by_file):
        old_count, new_count = data_row_counts_by_file[name]
        status = "MATCH" if old_count == new_count else "MISMATCH"
        if old_count != new_count:
            mismatched_files.append(name)
        print(f"  {name}: before={old_count}, after={new_count} ({status})")
    if mismatched_files:
        print(f"WARNING: {len(mismatched_files)} file(s) have a DIFFERENT data-row count "
              f"after the mass-line search change: {mismatched_files}")
    else:
        print(f"All {len(data_row_counts_by_file)} file(s): data-row counts are IDENTICAL "
              f"before and after the mass-line search change.")

    return pd.DataFrame(rows)


# --- PART 1b: mass-normalized effective mass -------------------------------
#
# The geometries being compared differ in total mass by up to ~40% - a raw
# effective-mass value confounds "this mode activates a lot of mass" with
# "this model is heavy". Normalizing by model mass separates the two. Both
# f_total and f_captured below are computed per (tree, model) group, never
# mixed across groups.

def add_mass_normalized_effmass_columns(modal_df: pd.DataFrame) -> pd.DataFrame:
    """
    Add, to a copy of modal_df, the per-mode columns:
      effmass_sum_kg  = EFFMASSx_kg + EFFMASSy_kg + EFFMASSz_kg
                        the physically additive total.
      effmass_norm_kg = sqrt(EFFMASSx_kg**2 + EFFMASSy_kg**2 + EFFMASSz_kg**2)
                        the vector norm - the same combination the APDL
                        macro's TTEM ranking already uses. Kept alongside
                        effmass_sum_kg (rather than replacing it) so it is
                        visible whether the two rank modes differently -
                        see print_ranking_agreement().
      f_total    = effmass_sum_kg / M_total_kg
                   the fraction of the WHOLE model mass this mode
                   activates.
      f_captured = effmass_sum_kg / (sum of effmass_sum_kg over every mode
                   of the SAME tree+model)
                   the fraction of what the EXTRACTED modes captured.
                   f_total and f_captured differ by however much mass the
                   extraction missed - see
                   print_modal_extraction_completeness().

    A file whose M_total_kg could not be found (find_total_mass_kg()
    returned None) gets NaN in f_total for every one of its modes, but
    still gets effmass_sum_kg/effmass_norm_kg/f_captured - those do not
    need the model mass.
    """
    df = modal_df.copy()  # work on a copy so we never modify the original data
    if len(df) == 0:
        return df

    df["effmass_sum_kg"] = df["EFFMASSx_kg"] + df["EFFMASSy_kg"] + df["EFFMASSz_kg"]
    df["effmass_norm_kg"] = np.sqrt(
        df["EFFMASSx_kg"] ** 2 + df["EFFMASSy_kg"] ** 2 + df["EFFMASSz_kg"] ** 2
    )

    df["f_total"] = df["effmass_sum_kg"] / df["M_total_kg"]

    group_captured_sum = df.groupby(["tree", "model"])["effmass_sum_kg"].transform("sum")
    df["f_captured"] = df["effmass_sum_kg"] / group_captured_sum

    return df


def print_modal_extraction_completeness(modal_df: pd.DataFrame) -> pd.DataFrame:
    """
    CHECK 1 - per tree and model: M_total_kg, the sum of effmass_sum_kg
    over every extracted mode, and their ratio (the captured mass
    fraction). A low ratio means the dominant mode may not be among the
    extracted modes at all.
    """
    print("\nModal extraction completeness (captured mass fraction), per tree/model:")
    if len(modal_df) == 0:
        print("  (no modal data)")
        return pd.DataFrame()
    rows = []
    for (tree, model), group in modal_df.groupby(["tree", "model"]):
        m_total = group["M_total_kg"].iloc[0]
        captured = group["effmass_sum_kg"].sum()
        ratio = (captured / m_total) if pd.notna(m_total) and m_total != 0 else float("nan")
        rows.append({
            "tree": tree, "model": model, "M_total_kg": m_total,
            "captured_effmass_sum_kg": captured, "captured_fraction": ratio,
        })
    result = pd.DataFrame(rows).sort_values(["tree", "model"]).reset_index(drop=True)
    if len(result) > 0:
        result = result.round({"M_total_kg": 2, "captured_effmass_sum_kg": 2, "captured_fraction": 4})
    print(result.to_string(index=False) if len(result) else "  (no modal data)")
    return result


def print_direction_coverage(modal_df: pd.DataFrame) -> pd.DataFrame:
    """
    CHECK 2 - per tree and model: the sum of EFFMASSx_kg, EFFMASSy_kg, and
    EFFMASSz_kg over every extracted mode, each as a percentage of
    M_total_kg. If the three differ a lot, the extraction covers the
    directions unevenly and mode selection is less reliable than it looks.
    """
    print("\nPer-direction effective-mass coverage (% of M_total_kg), per tree/model:")
    if len(modal_df) == 0:
        print("  (no modal data)")
        return pd.DataFrame()
    rows = []
    for (tree, model), group in modal_df.groupby(["tree", "model"]):
        m_total = group["M_total_kg"].iloc[0]
        row = {"tree": tree, "model": model, "M_total_kg": m_total}
        for column_name, direction in (("EFFMASSx_kg", "x"), ("EFFMASSy_kg", "y"), ("EFFMASSz_kg", "z")):
            direction_sum = group[column_name].sum()
            pct = (100.0 * direction_sum / m_total) if pd.notna(m_total) and m_total != 0 else float("nan")
            row[f"sum_EFFMASS{direction}_kg"] = direction_sum
            row[f"pct_{direction}"] = pct
        rows.append(row)
    result = pd.DataFrame(rows).sort_values(["tree", "model"]).reset_index(drop=True)
    if len(result) > 0:
        round_map = {"M_total_kg": 2}
        for direction in ("x", "y", "z"):
            round_map[f"sum_EFFMASS{direction}_kg"] = 2
            round_map[f"pct_{direction}"] = 2
        result = result.round(round_map)
    print(result.to_string(index=False) if len(result) else "  (no modal data)")
    return result


def print_top5_modes_by_f_total(modal_df: pd.DataFrame) -> pd.DataFrame:
    """
    CHECK 3 - per tree and model, the 5 modes with the largest f_total:
    frequency, EFFMASSx/y/z, effmass_sum_kg, effmass_norm_kg, f_total,
    f_captured.
    """
    print("\nTop 5 modes by f_total, per tree/model:")
    columns = ["tree", "model", "FREQ", "EFFMASSx_kg", "EFFMASSy_kg", "EFFMASSz_kg",
               "effmass_sum_kg", "effmass_norm_kg", "f_total", "f_captured"]
    if len(modal_df) == 0:
        print("  (no modal data)")
        return pd.DataFrame(columns=columns)
    top_groups = []
    for (tree, model), group in modal_df.groupby(["tree", "model"]):
        top_groups.append(group.nlargest(5, "f_total")[columns])
    result = (pd.concat(top_groups, ignore_index=True).sort_values(["tree", "model"], kind="stable")
              if top_groups else pd.DataFrame(columns=columns))
    if len(result) > 0:
        result = result.round({
            "FREQ": 4, "EFFMASSx_kg": 2, "EFFMASSy_kg": 2, "EFFMASSz_kg": 2,
            "effmass_sum_kg": 2, "effmass_norm_kg": 2, "f_total": 4, "f_captured": 4,
        })
    print(result.to_string(index=False) if len(result) else "  (no modal data)")
    return result


def print_ranking_agreement(modal_df: pd.DataFrame) -> pd.DataFrame:
    """
    CHECK 4 - per tree and model, whether the mode with the largest
    effmass_sum_kg is the SAME mode as the one with the largest
    effmass_norm_kg. Lists every case where the two rankings disagree.
    """
    print("\nTop-mode ranking agreement (sum-based vs. norm-based), per tree/model:")
    if len(modal_df) == 0:
        print("  (no modal data)")
        return pd.DataFrame()
    rows = []
    for (tree, model), group in modal_df.groupby(["tree", "model"]):
        top_by_sum = group.loc[group["effmass_sum_kg"].idxmax()]
        top_by_norm = group.loc[group["effmass_norm_kg"].idxmax()]
        agree = bool(top_by_sum["FREQ"] == top_by_norm["FREQ"])
        rows.append({
            "tree": tree, "model": model,
            "top_freq_by_sum": round(float(top_by_sum["FREQ"]), 4),
            "effmass_sum_kg_at_top_by_sum": round(float(top_by_sum["effmass_sum_kg"]), 2),
            "top_freq_by_norm": round(float(top_by_norm["FREQ"]), 4),
            "effmass_norm_kg_at_top_by_norm": round(float(top_by_norm["effmass_norm_kg"]), 2),
            "same_top_mode": agree,
        })
    result = pd.DataFrame(rows).sort_values(["tree", "model"]).reset_index(drop=True)
    print(result.to_string(index=False) if len(result) else "  (no modal data)")

    disagreements = result[~result["same_top_mode"]] if len(result) > 0 else result
    if len(disagreements) > 0:
        print(f"\n{len(disagreements)} tree/model case(s) where sum-based and norm-based "
              f"rankings pick a DIFFERENT top mode:")
        print(disagreements.to_string(index=False))
    else:
        print("\nEvery tree/model case: sum-based and norm-based rankings pick the SAME top mode.")

    return result


def print_n_modes_extracted(modal_df: pd.DataFrame) -> pd.DataFrame:
    """
    CHECK 5 - the number of modes extracted per tree and model.
    """
    print("\nNumber of modes extracted, per tree/model:")
    if len(modal_df) == 0:
        print("  (no modal data)")
        return pd.DataFrame()
    result = (modal_df.groupby(["tree", "model"]).size().rename("n_modes")
              .reset_index().sort_values(["tree", "model"]).reset_index(drop=True))
    print(result.to_string(index=False) if len(result) else "  (no modal data)")
    return result


def print_top_mode_vs_harmonic_peak(modal_df: pd.DataFrame, harmonic_df: pd.DataFrame) -> pd.DataFrame:
    """
    CHECK 6 - per tree and model: the frequency of the top mode under each
    ranking (sum-based and norm-based), next to the harmonic peak frequency
    for that tree/model (compute_harmonic_peak_frequency(), the same
    definition the frequency-summary figure and
    print_time_kg_peak_vs_harmonic_peak() already use), with the
    difference in percent of the harmonic peak.
    """
    print("\nTop mode frequency (sum-based and norm-based) vs. harmonic peak frequency, "
          "per tree/model:")
    if len(modal_df) == 0:
        print("  (no modal data)")
        return pd.DataFrame()
    rows = []
    for (tree, model), group in modal_df.groupby(["tree", "model"]):
        harmonic_peak = compute_harmonic_peak_frequency(harmonic_df, tree, model)
        freq_by_sum = float(group.loc[group["effmass_sum_kg"].idxmax(), "FREQ"])
        freq_by_norm = float(group.loc[group["effmass_norm_kg"].idxmax(), "FREQ"])

        if harmonic_peak is None or harmonic_peak == 0:
            diff_pct_sum = diff_pct_norm = float("nan")
        else:
            diff_pct_sum = 100.0 * (freq_by_sum - harmonic_peak) / harmonic_peak
            diff_pct_norm = 100.0 * (freq_by_norm - harmonic_peak) / harmonic_peak

        rows.append({
            "tree": tree, "model": model,
            "top_freq_by_sum_hz": round(freq_by_sum, 4),
            "top_freq_by_norm_hz": round(freq_by_norm, 4),
            "harmonic_peak_freq_hz": round(harmonic_peak, 4) if harmonic_peak is not None else None,
            "diff_pct_sum": round(diff_pct_sum, 2),
            "diff_pct_norm": round(diff_pct_norm, 2),
        })
    result = pd.DataFrame(rows).sort_values(["tree", "model"]).reset_index(drop=True)
    print(result.to_string(index=False) if len(result) else "  (no tree/model has both modal and harmonic data)")
    return result


# --- PART 2: harmonic data -------------------------------------------------

def load_harmonic_data(input_folder: str, file_pattern: str, order_to_tree_id: dict) -> pd.DataFrame:
    """
    Load every harmonic-response file matching file_pattern in
    input_folder into one combined DataFrame with columns: tree, ansys_no,
    scenario, model, source_file, frequency_hz, and the 9 displacement
    columns (NaN where a file does not provide a given column on a given
    row).

    order_to_tree_id translates the ANSYS filename's row POSITION into the
    canonical tree ID (see load_ansys_order_to_tree_id()) - a position that
    is not in that dict means the filename and the description file
    disagree, so the file is skipped rather than guessed.
    """
    matching_files = sorted(Path(input_folder).glob(file_pattern))
    print(f"Found {len(matching_files)} harmonic file(s) matching '{file_pattern}'")

    rows = []
    unknown_scenarios = set()

    for file_path in matching_files:
        match = HARMONIC_FILENAME_PATTERN.match(file_path.name)
        if match is None:
            print(f"WARNING: '{file_path.name}' does not match the expected harmonic "
                  f"filename pattern - skipping this file.")
            continue

        ansys_order = int(match.group(1))
        scenario = int(match.group(2))

        tree = order_to_tree_id.get(ansys_order)
        if tree is None:
            print(f"WARNING: '{file_path.name}' refers to position bk{ansys_order}, "
                  f"which does not exist in the description file - skipping "
                  f"this file.")
            continue

        model = scenario_to_model_name(scenario)
        if scenario not in SCENARIO_MODELS:
            unknown_scenarios.add(scenario)

        with open(file_path, "r", encoding="utf-8") as f:
            lines = [line for line in f.read().splitlines() if line.strip() != ""]

        # Line 1 is the header - we don't need to parse it since we already
        # know the fixed column order (HARMONIC_COLUMN_NAMES). Data starts
        # on line 2.
        for line in lines[1:]:
            values = [float(v) for v in line.split()]

            # A row may have fewer than 10 values if only some displacement
            # columns were populated for this run. Match whatever is present
            # to the FIRST N column names, and pad the rest with NaN.
            row_values = values + [float("nan")] * (len(HARMONIC_COLUMN_NAMES) - len(values))
            row = dict(zip(HARMONIC_COLUMN_NAMES, row_values))
            row["tree"] = tree
            row["ansys_no"] = ansys_order
            row["scenario"] = scenario
            row["model"] = model
            row["source_file"] = file_path.name
            rows.append(row)

    if unknown_scenarios:
        print(f"WARNING: scenario number(s) {sorted(unknown_scenarios)} are not "
              f"listed in SCENARIO_MODELS (harmonic data). They were loaded with a "
              f"placeholder name - add them to SCENARIO_MODELS so the plots get "
              f"a proper label.")

    harmonic_df = pd.DataFrame(rows)
    if len(harmonic_df) > 0:
        # Reorder columns to the order requested: tree, ansys_no, scenario,
        # model, source_file, then the rest.
        ordered_columns = ["tree", "ansys_no", "scenario", "model", "source_file"] + HARMONIC_COLUMN_NAMES
        harmonic_df = harmonic_df[ordered_columns]
    return harmonic_df


# --- PART 3: experimental frequency data -----------------------------------
#
# first_frequency is recorded per individual sensor/probe per pull (it is
# NOT pre-aggregated across sensors), and only for the dynamic "pull and
# release" measurements (M02-M05) - M01 rows always have first_frequency =
# NaN, since frequency cannot be estimated from a single static pull. That
# is expected, not a data problem.
#
# Because of the per-sensor duplication, there are two reasonable ways to
# turn this into one frequency value per pull, so we provide both:
#   - load_experimental_frequency_elasto():  use only the Elasto(90) probe
#     (the single sensor considered most reliable for this measurement).
#   - load_experimental_frequency_median():  pool whichever sensors are
#     present for a pull and take their median.

def _filter_experimental_base(csv_path: str, days: list, type_filter, reduction_level_filter, tree_filter, measurement_filter) -> pd.DataFrame:
    """
    Shared filtering step used by both frequency loaders below: load the
    CSV, convert the "tree" column to the canonical B21_S.. form (see
    bk_id_to_canonical()), add "reduction_level", and filter by
    day/type/reduction level/tree/measurement (None/value/list convention,
    same as filter_and_plot.py). Does NOT filter on probe or
    first_frequency - that part differs between the two loaders.
    """
    df = pd.read_csv(csv_path)

    # Keep the original experimental label, and convert "tree" to the same
    # canonical form the FEM loaders above now produce (see the WHY comment
    # above bk_id_to_canonical()).
    df["tree_experimental_id"] = df["tree"]
    df["tree"] = df["tree"].apply(bk_id_to_canonical)

    df = add_reduction_level_column(df)

    filtered = df[df["day"].isin(days)]
    filtered = filtered[filtered["type"] == type_filter]
    filtered = filter_column(filtered, "reduction_level", reduction_level_filter)
    filtered = filter_column(filtered, "tree", tree_filter)
    filtered = filter_column(filtered, "measurement", measurement_filter)
    return filtered


def load_experimental_frequency_elasto(csv_path: str, days: list, type_filter, reduction_level_filter, tree_filter, measurement_filter) -> pd.DataFrame:
    """
    Load experimental fundamental frequency readings from the Elasto(90)
    probe only. Returns a DataFrame with columns: tree, year,
    first_frequency - one row per Elasto(90) probe reading that has a
    valid (non-null) first_frequency.
    """
    filtered = _filter_experimental_base(csv_path, days, type_filter, reduction_level_filter, tree_filter, measurement_filter)
    filtered = filtered[filtered["probe"] == "Elasto(90)"]
    filtered = filtered[filtered["first_frequency"].notna()]

    filtered = filtered.copy()
    filtered["year"] = pd.to_datetime(filtered["day"]).dt.year

    print(f"Experimental frequency (Elasto(90) only): {len(filtered)} rows after filtering")
    return filtered[["tree", "year", "first_frequency"]]


def load_experimental_frequency_median(csv_path: str, days: list, type_filter, reduction_level_filter, tree_filter, measurement_filter) -> pd.DataFrame:
    """
    Load experimental fundamental frequency readings pooled across
    whatever sensors/probes recorded a value for each individual pull
    (Elasto(90), Pt3, Pt4, a02_z, a03_z, blueMaj, yellowMaj - whichever
    are non-null; not a fixed set).

    A "pull" is identified by (tree, day, measurement, pullNo). Within
    each pull, this takes the median of all non-null first_frequency
    values across its rows (i.e. across sensors), and skips the pull
    entirely if none of its sensors recorded a value.

    Returns a DataFrame with columns: tree, year, first_frequency - one
    row per pull that had at least one usable sensor reading.
    """
    filtered = _filter_experimental_base(csv_path, days, type_filter, reduction_level_filter, tree_filter, measurement_filter)

    filtered = filtered.copy()
    filtered["year"] = pd.to_datetime(filtered["day"]).dt.year

    grouped = filtered.groupby(["tree", "year", "day", "measurement", "pullNo"])

    rows = []
    for group_key, group in grouped:
        tree, year = group_key[0], group_key[1]
        values = group["first_frequency"].dropna()
        if len(values) == 0:
            continue  # no sensor recorded a frequency for this pull
        rows.append({"tree": tree, "year": year, "first_frequency": float(values.median())})

    result = pd.DataFrame(rows)
    print(f"Experimental frequency (median across sensors): {len(result)} pulls after filtering")
    return result


# --- PART 4: modal figure ---------------------------------------------------

def plot_modal_figure(tree: str, modal_df: pd.DataFrame, model_colors: dict, freq_df: pd.DataFrame, year_colors: dict, output_path: str):
    """
    Build and save the modal-analysis figure for one tree: 3 subplots
    ("EFF MASS x", "EFF MASS y", "EFF MASS z"), each showing one thin bar
    per mode at its frequency, height = effective mass in that direction,
    colored by reconstruction model, with the experimental fundamental
    frequency overlaid as a vertical band (see add_experimental_frequency_band
    above).

    These 3 subplots sit side by side, so only the LAST one gets its
    legends (model + experiment-year) - drawing them on every subplot would
    place subplot 0's and 1's "outside to the right" legends on top of
    their neighbouring subplots instead of in empty margin.
    """
    tree_modes = modal_df[modal_df["tree"] == tree]

    directions = [
        ("EFFMASSx_kg", "EFF MASS x"),
        ("EFFMASSy_kg", "EFF MASS y"),
        ("EFFMASSz_kg", "EFF MASS z"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Bar width is a small fraction of the frequency range, so bars stay
    # thin regardless of how wide the frequency axis ends up being.
    freq_range = tree_modes["FREQ"].max() - tree_modes["FREQ"].min()
    bar_width = max(freq_range * 0.01, 0.01)

    last_index = len(directions) - 1
    for i, (ax, (mass_column, title)) in enumerate(zip(axes, directions)):
        add_light_grid(ax)
        for model, group in tree_modes.groupby("model"):
            ax.bar(group["FREQ"], group[mass_column], width=bar_width,
                   color=model_colors[model], edgecolor="none",
                   linewidth=BAR_EDGE_LINEWIDTH, zorder=3)
        ax.set_title(title, fontsize=TITLE_FONTSIZE)
        ax.set_xlabel("Frequency (Hz)", fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylabel("Effective mass (kg)", fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
        # VIEW crop only (see FREQUENCY_AXIS_LIMITS's own comment) - every
        # bar is still drawn above, just not all necessarily visible.
        if FREQUENCY_AXIS_LIMITS is not None:
            ax.set_xlim(FREQUENCY_AXIS_LIMITS)

        is_last = (i == last_index)

        # Model legend, kept only on the last subplot - added as an
        # "artist" so it survives the experiment-year legend added next.
        # Restricted to models actually present in tree_modes (which is
        # already MODEL_INCLUDE-filtered by the caller when that filter is
        # active), not every key of model_colors.
        if is_last:
            present_models = [m for m in sorted(model_colors.keys()) if m in set(tree_modes["model"].unique())]
            model_handles = [
                Patch(facecolor=model_colors[model], edgecolor="none", label=model)
                for model in present_models
            ]
            model_legend = ax.legend(handles=model_handles, title="Model",
                                      loc="upper left", bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP),
                                      frameon=False, fontsize=LEGEND_FONTSIZE,
                                      title_fontsize=LEGEND_TITLE_FONTSIZE)
            ax.add_artist(model_legend)

        # Overlay the experimental fundamental frequency as a vertical band.
        add_experimental_frequency_band(ax, freq_df, tree, year_colors,
                                         show_legend=is_last,
                                         legend_bbox_y=LEGEND_Y_TOP - LEGEND_Y_STEP)

    title_suffix = f" - models shown: {', '.join(sorted(tree_modes['model'].unique()))}" if MODEL_INCLUDE is not None else ""
    fig.suptitle(f"Modal comparison - {tree}{title_suffix}")
    plt.tight_layout()
    fig.subplots_adjust(right=RIGHT_MARGIN)
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Modal figure saved to: {output_path}")


# --- PART 5: harmonic figure -------------------------------------------------

def compute_harmonic_frequency_range(tree_harmonic: pd.DataFrame):
    """
    Shared x-axis range logic for any harmonic figure - full or simplified
    - in this priority order:
      1. FREQUENCY_AXIS_LIMITS is set -> use it directly (the configured
         crop - VIEW ONLY, see that constant's own comment; nothing is
         filtered out of the underlying data).
      2. FREQUENCY_AXIS_LIMITS is None and there is harmonic data -> fall
         back to the harmonic curve's own frequency range - a point
         outside that range is then clipped at the plot's edge instead of
         silently expanding the whole x-axis to fit it.
      3. Neither is available -> (None, None) (matplotlib's own autoscale).
    Returns (freq_lo, freq_hi) - either can be None (case 3).
    """
    if FREQUENCY_AXIS_LIMITS is not None:
        return FREQUENCY_AXIS_LIMITS
    if len(tree_harmonic) > 0:
        return float(tree_harmonic["frequency_hz"].min()), float(tree_harmonic["frequency_hz"].max())
    return None, None


def plot_model_response_lines(ax, tree_harmonic: pd.DataFrame, column_name: str, model_colors: dict, use_abs: bool) -> None:
    """
    Draw one line per model of tree_harmonic[column_name] (or its absolute
    value if use_abs) vs frequency_hz, colored by model and labelled with
    the model name (so a plain ax.legend() picks the labels up) - the
    actual curve-drawing shared by every harmonic figure (the full 3-row
    one and the simplified single-quantity one), so a later style change
    to how these lines are drawn only has to be made here.
    """
    for model, group in tree_harmonic.groupby("model"):
        values = group[column_name].abs() if use_abs else group[column_name]
        if values.notna().sum() == 0:
            continue  # this model has no data for this column - skip its line
        ax.plot(group["frequency_hz"], values, color=model_colors[model],
                label=model, linewidth=LINE_WIDTH, zorder=3)


def add_model_color_legend(ax, bbox_x: float, bbox_y: float):
    """
    Small "Model" legend built from whatever labelled artists already sit
    on ax (e.g. the lines plot_model_response_lines() just drew) - shared
    by every figure that needs just this one legend, outside the axes.
    Does NOT add_artist() itself - the caller decides whether a further
    legend follows and needs this one protected.
    """
    return ax.legend(title="Model", loc="upper left", bbox_to_anchor=(bbox_x, bbox_y),
                     frameon=False, fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)


def plot_harmonic_simple_figure(tree: str, harmonic_df: pd.DataFrame, experimental_freq: pd.DataFrame,
                                model_colors: dict, year_colors: dict, output_path: str) -> None:
    """
    VARIANT 2 - simplified single-axis harmonic figure for the write-up:
    Top_ux ONLY, always as absolute values (Half_ux/Quart_ux dropped), with
    the experimental frequency band/median overlaid - no modal points, no
    effective-mass axis, no mass-fraction axis, a single plain y-axis.

    Reuses, rather than reimplements, the pieces plot_harmonic_figure()
    already needed for this: compute_harmonic_frequency_range() for the
    x-axis range (so FREQUENCY_AXIS_LIMITS applies here exactly as there),
    plot_model_response_lines() for the curves themselves, and
    add_experimental_frequency_band() for the band/median overlay.
    """
    tree_harmonic = harmonic_df[harmonic_df["tree"] == tree]
    freq_lo, freq_hi = compute_harmonic_frequency_range(tree_harmonic)

    fig, ax = plt.subplots(1, 1, figsize=HARMONIC_SIMPLE_FIGSIZE)
    add_light_grid(ax)

    plot_model_response_lines(ax, tree_harmonic, "Top_ux", model_colors, use_abs=True)

    ax.set_title("Top_ux (abs)", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Frequency (Hz)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("|Top_ux| (abs)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    if freq_lo is not None:
        ax.set_xlim(freq_lo, freq_hi)

    model_legend = add_model_color_legend(ax, LEGEND_BBOX_X, LEGEND_Y_TOP)
    ax.add_artist(model_legend)
    add_experimental_frequency_band(ax, experimental_freq, tree, year_colors, show_legend=True,
                                    legend_bbox_y=LEGEND_Y_TOP - LEGEND_Y_STEP, bbox_x=LEGEND_BBOX_X)

    title_suffix = f" - models shown: {', '.join(sorted(tree_harmonic['model'].unique()))}" if MODEL_INCLUDE is not None else ""
    fig.suptitle(f"Harmonic (simplified) - {tree}{title_suffix}")
    plt.tight_layout()
    fig.subplots_adjust(right=HARMONIC_SIMPLE_RIGHT_MARGIN)
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Simplified harmonic figure saved to: {output_path}")


def plot_harmonic_figure(tree: str, harmonic_df: pd.DataFrame, modal_df: pd.DataFrame, experimental_freq: pd.DataFrame, model_colors: dict, year_colors: dict, output_path: str, use_abs: bool = False):
    """
    Build and save the harmonic-response figure for one tree: 3 stacked
    line-plot rows ("Top_ux", "Half_ux", "Quart_ux") vs frequency, colored
    by reconstruction model, each with the experimental fundamental
    frequency overlaid as a vertical band (see add_experimental_frequency_band
    above), PLUS the tree's modal TIME_kg/TIME_incl_rots points on a
    secondary y-axis (points only, no line - see this module's own
    docstring) so it is visible whether the harmonic peak coincides with
    the mode carrying the most effective mass, PLUS f_captured/f_total on a
    THIRD y-axis (points only), since those are dimensionless fractions
    (0-1) and cannot share the kg-scaled secondary axis.

    f_captured is the fraction of the mass that the EXTRACTED modes
    captured, NOT of the model mass - on the current data the extraction
    captures anywhere from ~0.2% to ~18% of M_total depending on tree and
    model, so f_captured sums to 1 over a base that is itself tiny.
    Presented alone it would read as if a mode dominated the tree, when it
    may dominate almost nothing of the actual model. To keep that visible:
    f_total (fraction of the WHOLE model mass) is plotted alongside it with
    a different marker, and the captured-mass percentage for each model is
    annotated directly on the first row (not buried in a caption that could
    be read separately from the axis it explains).

    These 3 subplots are stacked vertically (not side by side), so each row
    has its own genuinely empty margin to its right - no collision, so each
    row keeps its own legend(s), just moved outside the axes.

    use_abs=True plots the absolute value of the response columns instead
    of the signed value (CHANGE 3) - the harmonic response is signed and
    crosses zero, which makes the resonance peak harder to read than a
    single clear maximum. The signed path (use_abs=False, the default)
    is unchanged and still available - this is an additional option, not
    a replacement.
    """
    tree_harmonic = harmonic_df[harmonic_df["tree"] == tree]
    tree_modal = modal_df[modal_df["tree"] == tree] if len(modal_df) > 0 else modal_df

    # X-axis limits, computed ONCE (used for both the primary curve axis and
    # the secondary modal-mass axis, which must share the same view) - see
    # compute_harmonic_frequency_range()'s own docstring for the 3-case
    # priority order. Whether anything actually falls outside the chosen
    # range is reported separately, never silently, by
    # print_frequency_range_check() (STEP 4) and
    # print_outside_axis_limits_report() - this function only draws.
    freq_lo, freq_hi = compute_harmonic_frequency_range(tree_harmonic)

    response_columns = [("Top_ux", "Top_ux"), ("Half_ux", "Half_ux"), ("Quart_ux", "Quart_ux")]

    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)

    for row_index, (ax, (column_name, title)) in enumerate(zip(axes, response_columns)):
        add_light_grid(ax)
        plot_model_response_lines(ax, tree_harmonic, column_name, model_colors, use_abs)
        ax.set_title(title + (" (abs)" if use_abs else ""), fontsize=TITLE_FONTSIZE)
        # Unambiguous y-axis label: the abs variant is explicitly marked so
        # the two variants (signed vs. absolute value) can never be confused
        # for one another when compared side by side.
        ax.set_ylabel(f"|{column_name}| (abs)" if use_abs else column_name,
                       fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)

        is_first = (row_index == 0)

        # Model legend, kept only on the first row (to avoid repeating it
        # 3 times) - added as an "artist" so it survives the
        # experiment-year legend added next.
        if is_first:
            model_legend = add_model_color_legend(ax, HARMONIC_LEGEND_BBOX_X, LEGEND_Y_TOP)
            ax.add_artist(model_legend)

        # Overlay the experimental fundamental frequency as a vertical band,
        # on the same axes as the FEM curves above. Stacked below the model
        # legend on the first row; at the top slot on the other two rows,
        # which have no model legend competing for that space. bbox_x is
        # pushed out (HARMONIC_LEGEND_BBOX_X, not the shared LEGEND_BBOX_X)
        # to clear the secondary axis added below.
        band_legend_y = (LEGEND_Y_TOP - LEGEND_Y_STEP) if is_first else LEGEND_Y_TOP
        add_experimental_frequency_band(ax, experimental_freq, tree, year_colors,
                                         show_legend=True, legend_bbox_y=band_legend_y,
                                         bbox_x=HARMONIC_LEGEND_BBOX_X)

        # --- STEP 2: modal TIME_kg / TIME_incl_rots overlay, secondary axis ---
        # Points only (no connecting line - these are discrete modes, not a
        # curve), drawn ABOVE the harmonic line (MODAL_MASS_ZORDER=5 >
        # zorder=3 above) so a point is never hidden under the curve.
        ax2 = ax.twinx()
        if freq_lo is not None:
            ax2.set_xlim(freq_lo, freq_hi)   # see the freq_lo/freq_hi comment above
        if len(tree_modal) > 0:
            for model, group in tree_modal.groupby("model"):
                color = model_colors.get(model)
                if color is None:
                    continue  # a modal-only model with no harmonic counterpart - not in model_colors
                ax2.scatter(group["FREQ"], group["TIME_kg"], color=color,
                            marker=MODAL_MASS_MARKER_TIME_KG, s=MODAL_MASS_MARKER_SIZE,
                            zorder=MODAL_MASS_ZORDER, edgecolors="none")
                ax2.scatter(group["FREQ"], group["TIME_incl_rots"], color=color,
                            marker=MODAL_MASS_MARKER_TIME_INCL_ROTS, s=MODAL_MASS_MARKER_SIZE,
                            zorder=MODAL_MASS_ZORDER, edgecolors="none")
        ax2.set_ylabel(MODAL_MASS_AXIS_LABEL, fontsize=AXIS_LABEL_FONTSIZE, color=MODAL_MASS_AXIS_COLOR)
        ax2.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE, labelcolor=MODAL_MASS_AXIS_COLOR)

        # Marker-shape legend for the two mass series, stacked one step
        # further down than the experiment-year legend - first row only,
        # same "avoid repeating 3 times" reasoning as the model legend.
        # Colour is neutral (MODAL_MASS_LEGEND_COLOR): colour already means
        # "model" via the legend above, shape alone distinguishes the two
        # mass series here.
        if is_first:
            mass_handles = [
                Line2D([0], [0], marker=MODAL_MASS_MARKER_TIME_KG, color="none",
                       markerfacecolor=MODAL_MASS_LEGEND_COLOR, markeredgecolor="none",
                       markersize=MODAL_MASS_LEGEND_MARKERSIZE, linestyle="None", label="TIME_kg"),
                Line2D([0], [0], marker=MODAL_MASS_MARKER_TIME_INCL_ROTS, color="none",
                       markerfacecolor=MODAL_MASS_LEGEND_COLOR, markeredgecolor="none",
                       markersize=MODAL_MASS_LEGEND_MARKERSIZE, linestyle="None", label="TIME_incl_rots"),
            ]
            ax2.legend(handles=mass_handles, title="Modal mass", loc="upper left",
                       bbox_to_anchor=(HARMONIC_LEGEND_BBOX_X, LEGEND_Y_TOP - 2 * LEGEND_Y_STEP),
                       frameon=False, fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

        # --- f_captured / f_total overlay, THIRD y-axis (CHANGE 1) -----------
        # A second twinx() on the same ax gives a new axes whose right spine
        # starts at the SAME position as ax2's (axes x=1.0) - moved further
        # right by FRACTION_AXIS_OFFSET so the two secondary axes don't sit
        # on top of one another. Scaled 0-1 (both f_captured and f_total are
        # fractions by construction) rather than autoscaled, so the axis
        # itself is a fixed, comparable 0-100% reference across every tree.
        ax3 = ax.twinx()
        ax3.spines["right"].set_position(("axes", 1.0 + FRACTION_AXIS_OFFSET))
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)   # keep it a pure overlay - no opaque background hiding ax/ax2
        ax3.set_ylim(0.0, 1.0)
        if freq_lo is not None:
            ax3.set_xlim(freq_lo, freq_hi)   # same shared view as ax2 above

        captured_pct_by_model = {}
        if len(tree_modal) > 0 and "f_captured" in tree_modal.columns:
            for model, group in tree_modal.groupby("model"):
                color = model_colors.get(model)
                if color is None:
                    continue  # a modal-only model with no harmonic counterpart - not in model_colors
                ax3.scatter(group["FREQ"], group["f_captured"], color=color,
                            marker=FRACTION_MARKER_F_CAPTURED, s=FRACTION_MARKER_SIZE,
                            zorder=FRACTION_ZORDER, edgecolors="none")
                ax3.scatter(group["FREQ"], group["f_total"], color=color,
                            marker=FRACTION_MARKER_F_TOTAL, s=FRACTION_MARKER_SIZE,
                            zorder=FRACTION_ZORDER, edgecolors="none")

                # Captured-mass percentage for THIS model - the number that
                # makes f_captured readable: it says what base the fractions
                # above are actually a share OF (see this function's own
                # docstring for why this cannot be left implicit).
                m_total = group["M_total_kg"].iloc[0] if "M_total_kg" in group.columns else float("nan")
                if pd.notna(m_total) and m_total != 0:
                    captured_pct_by_model[model] = 100.0 * group["effmass_sum_kg"].sum() / m_total

        ax3.set_ylabel(FRACTION_AXIS_LABEL, fontsize=AXIS_LABEL_FONTSIZE, color=FRACTION_AXIS_COLOR)
        ax3.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE, labelcolor=FRACTION_AXIS_COLOR)
        ax3.spines["right"].set_color(FRACTION_AXIS_COLOR)

        # Marker-shape legend for f_captured/f_total, stacked one step
        # further down than the modal-mass legend - first row only, same
        # "avoid repeating 3 times" reasoning as the other legends here.
        if is_first:
            fraction_handles = [
                Line2D([0], [0], marker=FRACTION_MARKER_F_CAPTURED, color="none",
                       markerfacecolor=FRACTION_LEGEND_COLOR, markeredgecolor="none",
                       markersize=FRACTION_LEGEND_MARKERSIZE, linestyle="None",
                       label="f_captured (of extracted modes)"),
                Line2D([0], [0], marker=FRACTION_MARKER_F_TOTAL, color="none",
                       markerfacecolor=FRACTION_LEGEND_COLOR, markeredgecolor="none",
                       markersize=FRACTION_LEGEND_MARKERSIZE, linestyle="None",
                       label="f_total (of M_total)"),
            ]
            ax3.legend(handles=fraction_handles, title="Mass fraction", loc="upper left",
                       bbox_to_anchor=(HARMONIC_LEGEND_BBOX_X, LEGEND_Y_TOP - 3 * LEGEND_Y_STEP),
                       frameon=False, fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

            # Captured-mass percentage, one line per model - see the
            # docstring's WHY: f_captured alone does not say what base it is
            # a share of, and that base is small (<20% of M_total here).
            if captured_pct_by_model:
                annotation_lines = [
                    f"captured {pct:.1f}% of M_total ({model})"
                    for model, pct in sorted(captured_pct_by_model.items())
                ]
                ax.text(0.01, 0.98, "\n".join(annotation_lines), transform=ax.transAxes,
                        ha="left", va="top", fontsize=FRACTION_ANNOTATION_FONTSIZE,
                        bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.85))

    axes[-1].set_xlabel("Frequency (Hz)", fontsize=AXIS_LABEL_FONTSIZE)

    title_suffix = f" - models shown: {', '.join(sorted(tree_harmonic['model'].unique()))}" if MODEL_INCLUDE is not None else ""
    fig.suptitle(f"Harmonic comparison - {tree}" + (" (abs)" if use_abs else "") + title_suffix)
    plt.tight_layout()
    fig.subplots_adjust(right=HARMONIC_RIGHT_MARGIN)
    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Harmonic figure saved to: {output_path}")


# --- PART 6: frequency summary figure ---------------------------------------

def compute_harmonic_peak_frequency(harmonic_df: pd.DataFrame, tree: str, model: str):
    """
    Return the frequency (Hz) at which the |Top_ux| harmonic response
    reaches its maximum for this tree/model, or None if there is no usable
    data (no rows, or Top_ux entirely missing).
    """
    subset = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)]
    values = subset["Top_ux"].abs()
    if values.notna().sum() == 0:
        return None
    idx = values.idxmax()
    return float(subset.loc[idx, "frequency_hz"])


def compute_modal_peak_frequency(modal_df: pd.DataFrame, tree: str, model: str):
    """
    Return (frequency_hz, direction) for the mode with the highest
    effective mass across EFFMASSx_kg, EFFMASSy_kg and EFFMASSz_kg, for
    this tree/model - direction is "x", "y" or "z", whichever direction's
    EFF MASS was largest for the winning mode. Returns (None, None) if
    there is no usable data.
    """
    subset = modal_df[(modal_df["tree"] == tree) & (modal_df["model"] == model)]
    if len(subset) == 0:
        return None, None

    direction_columns = [("EFFMASSx_kg", "x"), ("EFFMASSy_kg", "y"), ("EFFMASSz_kg", "z")]
    best_value, best_freq, best_direction = None, None, None
    for _, row in subset.iterrows():
        for column_name, direction in direction_columns:
            value = row[column_name]
            if pd.isna(value):
                continue
            if best_value is None or value > best_value:
                best_value = value
                best_freq = float(row["FREQ"])
                best_direction = direction

    return best_freq, best_direction


def compute_experimental_frequency_stats(freq_df: pd.DataFrame, tree: str):
    """
    Return (median, p25, p75, n) of first_frequency for this tree, pooled
    across all years, or (None, None, None, 0) if there is no data.
    """
    values = freq_df[freq_df["tree"] == tree]["first_frequency"].dropna()
    if len(values) == 0:
        return None, None, None, 0
    return (float(values.median()), float(values.quantile(0.25)),
            float(values.quantile(0.75)), len(values))


def print_time_mass_sanity_check(modal_df: pd.DataFrame) -> dict:
    """
    STEP 3 - unit sanity check on TIME_kg/TIME_incl_rots, PRINT ONLY: does
    not branch, does not adjust/swap/exclude anything based on what it
    finds (see this task's own instruction - both series are plotted on the
    shared kg axis regardless of the result printed here).

    Reports, across every loaded modal file: min/max of each column, the
    ratio TIME_incl_rots/TIME_kg (min/median/max), and the count of modes
    where TIME_incl_rots < TIME_kg (expected to be zero - rotational DOF
    add a contribution on top of the translational one, so the "incl rots"
    value cannot be smaller. If it is not zero, the offending files/modes
    are named, not silently accepted or removed).
    """
    print("\nTIME_kg / TIME_incl_rots unit sanity check (print only - see module docstring):")
    if len(modal_df) == 0 or "TIME_kg" not in modal_df.columns:
        print("  (no modal data loaded - nothing to check)")
        return {}

    time_kg = modal_df["TIME_kg"].dropna()
    time_incl = modal_df["TIME_incl_rots"].dropna()
    if len(time_kg) == 0 or len(time_incl) == 0:
        print("  (TIME_kg/TIME_incl_rots are entirely NaN - every row was too short to contain them)")
        return {}

    print(f"  TIME_kg        : min={time_kg.min():.4f} kg  max={time_kg.max():.4f} kg  (n={len(time_kg)})")
    print(f"  TIME_incl_rots : min={time_incl.min():.4f} kg  max={time_incl.max():.4f} kg  (n={len(time_incl)})")

    both = modal_df[["source_file", "tree", "model", "FREQ", "TIME_kg", "TIME_incl_rots"]].dropna(
        subset=["TIME_kg", "TIME_incl_rots"])
    ratio = both["TIME_incl_rots"] / both["TIME_kg"]
    print(f"  ratio TIME_incl_rots / TIME_kg : min={ratio.min():.4f}  median={ratio.median():.4f}  max={ratio.max():.4f}")

    below = both[both["TIME_incl_rots"] < both["TIME_kg"]]
    print(f"  modes where TIME_incl_rots < TIME_kg : {len(below)} (expected 0)")
    if len(below) > 0:
        print("  Rows where this happens - reported as-is, nothing adjusted/swapped/excluded:")
        print(below.to_string(index=False))

    return {
        "time_kg_min": float(time_kg.min()), "time_kg_max": float(time_kg.max()),
        "time_incl_rots_min": float(time_incl.min()), "time_incl_rots_max": float(time_incl.max()),
        "ratio_min": float(ratio.min()), "ratio_median": float(ratio.median()), "ratio_max": float(ratio.max()),
        "n_below": len(below),
    }


def print_frequency_range_check(modal_df: pd.DataFrame, harmonic_df: pd.DataFrame) -> pd.DataFrame:
    """
    STEP 4 - for every tree/model that has both modal and harmonic data,
    report the harmonic sweep's frequency range, the modal frequency range,
    and whether every modal frequency falls inside the sweep. A mode
    outside the swept band sits off the edge of the harmonic figure's
    secondary-axis overlay, which looks like missing data rather than what
    it is - this makes that visible instead of silent.
    """
    print("\nFrequency-range check (harmonic sweep vs. modal frequencies), per tree/model:")
    trees = sorted(set(modal_df["tree"].unique()) | set(harmonic_df["tree"].unique())) \
        if len(modal_df) > 0 or len(harmonic_df) > 0 else []

    rows = []
    for tree in trees:
        models = sorted(
            set(modal_df.loc[modal_df["tree"] == tree, "model"].unique() if len(modal_df) > 0 else [])
            | set(harmonic_df.loc[harmonic_df["tree"] == tree, "model"].unique() if len(harmonic_df) > 0 else [])
        )
        for model in models:
            h = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)]
            m = modal_df[(modal_df["tree"] == tree) & (modal_df["model"] == model)]
            if len(h) == 0 or len(m) == 0:
                continue
            sweep_lo, sweep_hi = float(h["frequency_hz"].min()), float(h["frequency_hz"].max())
            modal_freq = m["FREQ"].dropna()
            modal_lo, modal_hi = float(modal_freq.min()), float(modal_freq.max())
            inside = (modal_freq >= sweep_lo) & (modal_freq <= sweep_hi)
            rows.append({
                "tree": tree, "model": model,
                "sweep_lo_hz": round(sweep_lo, 4), "sweep_hi_hz": round(sweep_hi, 4),
                "modal_lo_hz": round(modal_lo, 4), "modal_hi_hz": round(modal_hi, 4),
                "all_modes_in_sweep": bool(inside.all()),
                "n_modes_outside": int((~inside).sum()),
            })

    result = pd.DataFrame(rows)
    print(result.to_string(index=False) if len(result) else "  (no tree/model has both modal and harmonic data)")
    return result


def print_time_kg_peak_vs_harmonic_peak(modal_df: pd.DataFrame, harmonic_df: pd.DataFrame) -> pd.DataFrame:
    """
    VERIFICATION (e) - for every tree/model, the frequency of the mode with
    the largest TIME_kg next to the harmonic peak frequency (from
    compute_harmonic_peak_frequency(), the same definition the frequency-
    summary figure already uses), and their difference in percent of the
    harmonic peak.
    """
    print("\nMode with the largest TIME_kg vs. harmonic peak frequency, per tree/model:")
    trees = sorted(set(modal_df["tree"].unique()) | set(harmonic_df["tree"].unique())) \
        if len(modal_df) > 0 or len(harmonic_df) > 0 else []

    rows = []
    for tree in trees:
        models = sorted(
            set(modal_df.loc[modal_df["tree"] == tree, "model"].unique() if len(modal_df) > 0 else [])
            | set(harmonic_df.loc[harmonic_df["tree"] == tree, "model"].unique() if len(harmonic_df) > 0 else [])
        )
        for model in models:
            m = modal_df[(modal_df["tree"] == tree) & (modal_df["model"] == model)]
            harmonic_peak = compute_harmonic_peak_frequency(harmonic_df, tree, model)
            time_kg = m["TIME_kg"].dropna()
            if len(time_kg) == 0 or harmonic_peak is None:
                continue
            modal_freq_at_max = float(m.loc[time_kg.idxmax(), "FREQ"])
            pct_diff = (100.0 * (modal_freq_at_max - harmonic_peak) / harmonic_peak
                        if harmonic_peak else float("nan"))
            rows.append({
                "tree": tree, "model": model,
                "modal_freq_at_max_TIME_kg_hz": round(modal_freq_at_max, 4),
                "harmonic_peak_freq_hz": round(harmonic_peak, 4),
                "diff_pct": round(pct_diff, 2),
            })

    result = pd.DataFrame(rows)
    print(result.to_string(index=False) if len(result) else "  (no tree/model has both modal TIME_kg and a harmonic peak)")
    return result


def print_outside_axis_limits_report(modal_df: pd.DataFrame, harmonic_df: pd.DataFrame,
                                      freq_elasto: pd.DataFrame, freq_median: pd.DataFrame):
    """
    Report, per tree and model, every plotted item whose frequency falls
    OUTSIDE FREQUENCY_AXIS_LIMITS: harmonic sweep points, modal modes, and
    experimental frequency values (both the elasto and median variants).

    FREQUENCY_AXIS_LIMITS crops the VIEW on the harmonic/modal figures only
    (see that constant's own comment) - nothing here is removed from
    modal_df, harmonic_df, frequency_deviation.csv, or any other printed
    number. This function exists so an item hidden by the axis crop is
    never mistaken for one that does not exist: it names every such item
    and gives its frequency, using the SAME full (uncropped) DataFrames
    every other check/figure in this script uses.

    Skipped (with a one-line note) if FREQUENCY_AXIS_LIMITS is None - there
    is no crop active, so nothing can be hidden by it.
    """
    print(f"\nItems outside the figure frequency-axis limits (FREQUENCY_AXIS_LIMITS) - "
          f"still in every table/CSV, just not visible on the harmonic/modal figures:")
    if FREQUENCY_AXIS_LIMITS is None:
        print("  FREQUENCY_AXIS_LIMITS is None - no view crop active, nothing can be hidden.")
        return

    lo, hi = FREQUENCY_AXIS_LIMITS
    print(f"  limits: {lo}-{hi} Hz")

    trees = sorted(
        set(modal_df["tree"].unique() if len(modal_df) > 0 else [])
        | set(harmonic_df["tree"].unique() if len(harmonic_df) > 0 else [])
    )
    any_found = False
    for tree in trees:
        models = sorted(
            set(modal_df.loc[modal_df["tree"] == tree, "model"].unique() if len(modal_df) > 0 else [])
            | set(harmonic_df.loc[harmonic_df["tree"] == tree, "model"].unique() if len(harmonic_df) > 0 else [])
        )
        for model in models:
            h = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)] if len(harmonic_df) > 0 else harmonic_df
            m = modal_df[(modal_df["tree"] == tree) & (modal_df["model"] == model)] if len(modal_df) > 0 else modal_df

            h_outside = h[(h["frequency_hz"] < lo) | (h["frequency_hz"] > hi)] if len(h) > 0 else h
            m_outside = m[(m["FREQ"] < lo) | (m["FREQ"] > hi)] if len(m) > 0 else m

            if len(h_outside) == 0 and len(m_outside) == 0:
                continue
            any_found = True
            print(f"  {tree} / {model}:")
            if len(h_outside) > 0:
                # The harmonic sweep is a dense, uniform frequency grid (a
                # continuous curve on the figure, visibly cut off at the
                # axis edge - not individually meaningful points), so this
                # is reported as count + range rather than all N values -
                # unlike the modal modes/experimental values below, which
                # ARE individually meaningful discrete points and are
                # listed in full.
                freqs = h_outside["frequency_hz"]
                print(f"    harmonic sweep points outside range: {len(h_outside)} of {len(h)} "
                      f"(range {float(freqs.min()):.4f}-{float(freqs.max()):.4f} Hz)")
            if len(m_outside) > 0:
                freqs = sorted(round(float(f), 4) for f in m_outside["FREQ"])
                print(f"    modal modes outside range: {len(m_outside)} of {len(m)} -> {freqs}")

    # Experimental frequency values, both variants - drawn on the SAME
    # x-axis as the harmonic/modal figures above, so subject to the same
    # crop. Reported per tree/year (experimental data carries no "model").
    for variant_label, freq_df in [("elasto", freq_elasto), ("median", freq_median)]:
        if len(freq_df) == 0:
            continue
        outside = freq_df[(freq_df["first_frequency"] < lo) | (freq_df["first_frequency"] > hi)]
        if len(outside) == 0:
            continue
        any_found = True
        print(f"  Experimental frequency ({variant_label} variant) outside range:")
        for _, row in outside.sort_values(["tree", "year"]).iterrows():
            print(f"    {row['tree']} / {row['year']}: {row['first_frequency']:.4f} Hz")

    if not any_found:
        print("  (nothing falls outside - every plotted item is visible within the current limits)")


def build_frequency_summary_figure(tree_order: list, modal_df: pd.DataFrame, harmonic_df: pd.DataFrame,
                                    model_colors: dict, freq_df: pd.DataFrame, variant_label: str,
                                    output_path: str, show_modal: bool = SUMMARY_SHOW_MODAL,
                                    visible_models: list = None) -> pd.DataFrame:
    """
    Build and save one summary figure with all trees side by side (same
    tree_order as the other figures): for each tree and model, two points -
    the harmonic peak frequency (max |Top_ux|) and the modal frequency of
    the dominant EFF MASS direction - colored by model, distinguished by
    marker shape, overlaid with the experimental frequency for that tree
    (a grey point at the median, with a thin vertical line spanning the
    25th-75th percentile).

    Prints, for every tree/model, the two FEM frequencies and which EFF
    MASS direction won the modal one, then the SAME relative-deviation
    numbers (FEM vs experimental median), using the formula:

        rel_dev_pct = 100 * (FEM - experimental_median) / experimental_median

    as a grouped-bar subplot BELOW the main panel (sharing its x-axis of
    trees) instead of the old multi-line text footnote, which ran to 24
    lines and overran the axes - see SUMMARY_DEVIATION_HEIGHT_RATIO's own
    comment. The bars are read straight off the same deviation_table this
    function already builds and returns for the CSV - the subplot is a
    display of those numbers, not a second calculation of them, so it can
    never disagree with FREQUENCY_DEVIATION_CSV_PATH.

    show_modal controls whether the modal_peak quantity (points in the
    upper panel, bars in the deviation panel, and their legend entries) is
    drawn at all - False leaves it out entirely rather than drawing it
    empty; frequency_deviation.csv still contains it either way (that
    table's own construction below does not depend on this argument) - and
    the deviation_rows this function builds and RETURNS also always
    contain both quantities regardless of show_modal, so a caller that
    does not want modal_peak rows duplicated into the CSV (see VARIANT 3 -
    frequency_summary_harmonic_{variant}.png in main()) must discard the
    return value rather than concatenate it in.

    Defaults to the module-level SUMMARY_SHOW_MODAL so every EXISTING
    caller (the plain frequency_summary_{variant}.png figures) is
    unaffected - only a caller that passes show_modal explicitly (the
    harmonic-only VARIANT 3 figures) gets different behaviour. The title
    gets a " - harmonic only" suffix whenever show_modal is False, so the
    two kinds of figure are never visually ambiguous.

    visible_models=None means "no MODEL_INCLUDE filter" - every model with
    data is drawn (subject only to show_modal, as before). When set, a
    model outside it is skipped when DRAWING (points, bars, legend
    entries) but its row is still added to deviation_rows/deviation_table
    UNCONDITIONALLY, same guarantee as show_modal above - MODEL_INCLUDE is
    a filter on what is drawn, never on what is computed or written to
    FREQUENCY_DEVIATION_CSV_PATH.

    Returns the deviation table (with a "variant" column identifying which
    experimental variant - elasto/median - it was computed against), for
    the caller to combine with the other variant's table and write to
    FREQUENCY_DEVIATION_CSV_PATH.
    """
    # Dodge geometry is computed over visible_models when MODEL_INCLUDE is
    # active (a more compact spread, matching that filter's whole point of
    # a less dense figure) - otherwise over every model in model_colors, as
    # before.
    all_models = sorted(visible_models) if visible_models is not None else sorted(model_colors.keys())
    offsets = compute_dodge_offsets(len(all_models), SUMMARY_GROUP_WIDTH)
    offset_by_model = dict(zip(all_models, offsets))

    base_height = 6.0
    fig, (ax, ax_dev) = plt.subplots(
        2, 1, figsize=(max(SUMMARY_MIN_WIDTH_INCHES, 1.8 * len(tree_order)),
                       base_height * (1.0 + SUMMARY_DEVIATION_HEIGHT_RATIO)),
        sharex=True, gridspec_kw={"height_ratios": [1.0, SUMMARY_DEVIATION_HEIGHT_RATIO]},
    )
    add_light_grid(ax)

    deviation_rows = []

    for tree_index, tree in enumerate(tree_order):
        exp_median, exp_p25, exp_p75, exp_n = compute_experimental_frequency_stats(freq_df, tree)

        if exp_n > 0:
            ax.vlines(tree_index, exp_p25, exp_p75, color=SUMMARY_EXPERIMENTAL_COLOR,
                      linewidth=SUMMARY_ERRORBAR_LINEWIDTH, zorder=2)
            ax.scatter(tree_index, exp_median, color=SUMMARY_EXPERIMENTAL_COLOR,
                       marker=SUMMARY_MARKER_EXPERIMENTAL, s=MARKER_SIZE ** 2, zorder=4)
        else:
            print(f"Note: no experimental frequency data for tree '{tree}' ({variant_label}) "
                  f"- no experimental point drawn for it in the summary figure.")

        tree_models = sorted(
            set(modal_df[modal_df["tree"] == tree]["model"].unique())
            | set(harmonic_df[harmonic_df["tree"] == tree]["model"].unique())
        )

        for model in tree_models:
            x = tree_index + offset_by_model.get(model, 0.0)
            color = model_colors[model]
            is_visible = visible_models is None or model in visible_models

            harmonic_freq = compute_harmonic_peak_frequency(harmonic_df, tree, model)
            modal_freq, winning_direction = compute_modal_peak_frequency(modal_df, tree, model)

            harmonic_text = f"{harmonic_freq:.4f} Hz (Top_ux)" if harmonic_freq is not None else "no data"
            modal_text = (f"{modal_freq:.4f} Hz (EFFMASS{winning_direction})"
                          if modal_freq is not None else "no data")
            exp_text = (f"{exp_median:.4f} Hz (n={exp_n})" if exp_n > 0 else "no data")
            print(f"{tree} / {model} ({variant_label}): harmonic peak = {harmonic_text}, "
                  f"modal peak = {modal_text}, experimental median = {exp_text}"
                  f"{'' if is_visible else '  [not drawn - excluded by MODEL_INCLUDE]'}")

            # MODEL_INCLUDE only gates DRAWING - deviation_rows below is
            # always built for every model, regardless of is_visible.
            if is_visible:
                if harmonic_freq is not None:
                    ax.scatter(x, harmonic_freq, color=color, marker=SUMMARY_MARKER_HARMONIC,
                               s=MARKER_SIZE ** 2, zorder=5)
                if show_modal and modal_freq is not None:
                    ax.scatter(x, modal_freq, color=color, marker=SUMMARY_MARKER_MODAL,
                               s=MARKER_SIZE ** 2, zorder=5)

            for quantity, fem_value in (("harmonic_peak", harmonic_freq), ("modal_peak", modal_freq)):
                if fem_value is None or exp_n == 0 or exp_median in (None, 0):
                    rel_dev_pct = float("nan")
                else:
                    rel_dev_pct = 100.0 * (fem_value - exp_median) / exp_median

                deviation_rows.append({
                    "tree": tree, "model": model, "quantity": quantity,
                    "fem_value_hz": fem_value, "exp_median_hz": exp_median,
                    "exp_n": exp_n, "rel_dev_pct": rel_dev_pct, "variant": variant_label,
                })

    ax.tick_params(labelbottom=False, labelsize=TICK_LABEL_FONTSIZE)  # x tick labels live on ax_dev below (shared x-axis)
    ax.set_ylabel("Frequency (Hz)", fontsize=AXIS_LABEL_FONTSIZE)
    title_suffix = "" if show_modal else " - harmonic only"
    if visible_models is not None:
        title_suffix += f" - models shown: {', '.join(visible_models)}"
    ax.set_title(f"Characteristic frequency summary ({variant_label}){title_suffix}", fontsize=TITLE_FONTSIZE)

    # Legend restricted to models actually drawn - a MODEL_INCLUDE-excluded
    # model's colour stays reserved (model_colors is unchanged) but its
    # legend entry does not appear here.
    drawn_models_with_data = {row["model"] for row in deviation_rows}
    present_models = [m for m in sorted(model_colors.keys())
                      if (visible_models is None or m in visible_models) and m in drawn_models_with_data]
    model_handles = [
        Line2D([], [], marker="o", linestyle="", color="black", markerfacecolor=model_colors[model],
               markersize=MARKER_SIZE, label=model)
        for model in present_models
    ]
    model_legend = ax.legend(handles=model_handles, title="Model", loc="upper left",
                              bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP), frameon=False,
                              fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)
    ax.add_artist(model_legend)

    # show_modal=False: the modal_peak marker/hatch entries are left OUT of
    # these lists entirely, not drawn disabled/greyed - an absent legend
    # entry, not an orphaned one pointing at nothing on the axes.
    quantity_handles = [
        Line2D([], [], marker=SUMMARY_MARKER_HARMONIC, linestyle="", color="black",
               markerfacecolor="white", markersize=MARKER_SIZE, label="Harmonic peak (Top_ux)"),
    ]
    if show_modal:
        quantity_handles.append(
            Line2D([], [], marker=SUMMARY_MARKER_MODAL, linestyle="", color="black",
                   markerfacecolor="white", markersize=MARKER_SIZE, label="Modal peak (max EFF MASS)")
        )
    quantity_handles.append(
        Line2D([], [], marker=SUMMARY_MARKER_EXPERIMENTAL, linestyle="", color=SUMMARY_EXPERIMENTAL_COLOR,
               markersize=MARKER_SIZE, label="Experiment (median, 25-75%)")
    )
    quantity_legend = ax.legend(handles=quantity_handles, title="Quantity", loc="upper left",
                                 bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - LEGEND_Y_STEP), frameon=False,
                                 fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)
    ax.add_artist(quantity_legend)

    deviation_table = pd.DataFrame(deviation_rows)
    if len(deviation_table) > 0:
        deviation_table = deviation_table.round(
            {"fem_value_hz": 4, "exp_median_hz": 4, "rel_dev_pct": 1})

    # --- Deviation panel (CHANGE 2) - grouped bars, read straight off
    # deviation_table (built above, unchanged by this display step) so the
    # bars can never disagree with FREQUENCY_DEVIATION_CSV_PATH. One bar
    # per model per quantity: solid = modal_peak, hatched = harmonic_peak
    # (SUMMARY_HARMONIC_HATCH/SUMMARY_MODAL_HATCH) - colour still means
    # "model", as everywhere else. ---
    add_light_grid(ax_dev)
    ax_dev.axhline(0.0, color=SUMMARY_ZERO_LINE_COLOR, linewidth=SUMMARY_ZERO_LINE_WIDTH, zorder=2)

    bar_quantities = [("harmonic_peak", SUMMARY_HARMONIC_HATCH, "Harmonic peak (hatched)")]
    if show_modal:
        bar_quantities.append(("modal_peak", SUMMARY_MODAL_HATCH, "Modal peak (solid)"))
    quantity_names = [q for q, _, _ in bar_quantities]
    # Sub-offsets within one model's dodge slot: 1 quantity -> centered on
    # the slot (no gap left where a sibling bar would have been when
    # show_modal is False); 2 quantities -> side by side.
    sub_offsets = compute_dodge_offsets(len(bar_quantities), SUMMARY_DEVIATION_BAR_WIDTH * len(bar_quantities))

    tree_index_by_tree = {tree: i for i, tree in enumerate(tree_order)}
    models_with_data = {row["model"] for row in deviation_rows
                        if row["quantity"] in quantity_names and pd.notna(row["rel_dev_pct"])}
    drawn_any_bar = False
    if len(deviation_table) > 0:
        for _, row in deviation_table.iterrows():
            if row["quantity"] not in quantity_names or pd.isna(row["rel_dev_pct"]):
                continue  # modal_peak rows still exist in deviation_table when SUMMARY_SHOW_MODAL is False (the CSV is unaffected) - just not drawn; a missing comparison is a gap, not a zero-height bar
            if visible_models is not None and row["model"] not in visible_models:
                continue  # excluded by MODEL_INCLUDE - a different, independent reason for a gap
            drawn_any_bar = True
            quantity_index = quantity_names.index(row["quantity"])
            x = (tree_index_by_tree[row["tree"]] + offset_by_model.get(row["model"], 0.0)
                 + sub_offsets[quantity_index])
            ax_dev.bar(x, row["rel_dev_pct"], width=SUMMARY_DEVIATION_BAR_WIDTH,
                       color=model_colors[row["model"]], hatch=bar_quantities[quantity_index][1],
                       edgecolor=SUMMARY_DEVIATION_BAR_EDGECOLOR,
                       linewidth=SUMMARY_DEVIATION_BAR_EDGE_LINEWIDTH, zorder=3)

    if models_with_data and not drawn_any_bar:
        print(f"  WARNING: MODEL_INCLUDE ({visible_models}) leaves no model with data "
              f"({sorted(models_with_data)}) drawn in this deviation subplot ({variant_label}) - "
              f"it will be empty.")

    ax_dev.set_xticks(range(len(tree_order)))
    ax_dev.set_xticklabels(tree_order)
    ax_dev.set_ylabel("Relative deviation (%)", fontsize=AXIS_LABEL_FONTSIZE)
    ax_dev.tick_params(labelsize=TICK_LABEL_FONTSIZE)

    # Fill/hatch legend for the deviation panel - shared with the upper
    # panel (attached to `ax`, drawn once, outside the axes), not a second
    # legend living on ax_dev. Absent (not disabled) when SUMMARY_SHOW_MODAL
    # is False, matching quantity_handles above.
    fill_handles = [
        Patch(facecolor="0.6", edgecolor=SUMMARY_DEVIATION_BAR_EDGECOLOR, hatch=hatch, label=label)
        for _, hatch, label in bar_quantities
    ]
    ax.legend(handles=fill_handles, title="Deviation bar fill", loc="upper left",
              bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - 2 * LEGEND_Y_STEP), frameon=False,
              fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    plt.tight_layout()
    fig.subplots_adjust(right=SUMMARY_RIGHT_MARGIN, hspace=0.08)

    plt.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"Frequency summary figure saved to: {output_path}")

    print(f"\nRelative deviation, FEM vs experimental median frequency ({variant_label}):")
    print(deviation_table.to_string(index=False) if len(deviation_table) > 0 else "  (no data)")

    return deviation_table


def main():
    # --- Part -1: MODEL_INCLUDE / MODEL_INCLUDE_SUFFIX sanity check -------------
    # Writing filtered figures under the same filenames as the full ones
    # would silently and irreversibly overwrite them - refuse outright.
    if MODEL_INCLUDE is not None and not MODEL_INCLUDE_SUFFIX:
        print("MODEL_INCLUDE is set but MODEL_INCLUDE_SUFFIX is empty - refusing to run: "
              "writing filtered figures under the SAME filenames as the full ones would "
              "silently and irreversibly overwrite them. Set MODEL_INCLUDE_SUFFIX to a "
              "non-empty string and re-run.")
        return

    # --- Part 0: build the position -> tree ID lookup ---------------------------
    print(f"Loading tree list from: {GEOMETRY_FILE_PATH}")
    order_to_tree_id = load_ansys_order_to_tree_id(GEOMETRY_FILE_PATH)
    print(f"Mapping built for {len(order_to_tree_id)} tree(s):")
    for ansys_order in sorted(order_to_tree_id):
        print(f"  bk{ansys_order}  ->  {order_to_tree_id[ansys_order]}")
    print()

    # --- Parts 1-3: load all data sources --------------------------------------
    modal_df = load_modal_data(INPUT_FOLDER, MODAL_FILE_PATTERN, order_to_tree_id)
    harmonic_df = load_harmonic_data(INPUT_FOLDER, HARMONIC_FILE_PATTERN, order_to_tree_id)

    # Two variants of the experimental frequency data (see PART 3 above for
    # why): Elasto(90) only, and median across whatever sensors recorded
    # each pull.
    freq_elasto = load_experimental_frequency_elasto(
        EXPERIMENTAL_CSV_PATH, DAYS_TO_PROCESS, TYPE_FILTER,
        REDUCTION_LEVEL_FILTER, TREE_FILTER, FREQUENCY_MEASUREMENT_FILTER,
    )
    freq_median = load_experimental_frequency_median(
        EXPERIMENTAL_CSV_PATH, DAYS_TO_PROCESS, TYPE_FILTER,
        REDUCTION_LEVEL_FILTER, TREE_FILTER, FREQUENCY_MEASUREMENT_FILTER,
    )

    if len(modal_df) == 0 and len(harmonic_df) == 0:
        print("No modal or harmonic files found - nothing to plot.")
        return

    # Mass-normalized effective-mass columns (effmass_sum_kg,
    # effmass_norm_kg, f_total, f_captured) - computed purely from columns
    # already loaded above (EFFMASSx/y/z_kg and M_total_kg), added here so
    # every check/figure below that uses modal_df sees them too.
    modal_df = add_mass_normalized_effmass_columns(modal_df)

    # --- STEP 3 / STEP 4 / VERIFICATION (e): print-only checks, run once on
    # the full combined data, before any figure is drawn. None of these
    # branch on their own result or change what gets plotted (see each
    # function's own docstring). ---
    print_time_mass_sanity_check(modal_df)
    print_frequency_range_check(modal_df, harmonic_df)
    print_time_kg_peak_vs_harmonic_peak(modal_df, harmonic_df)
    print_outside_axis_limits_report(modal_df, harmonic_df, freq_elasto, freq_median)

    # --- Mass-normalized effective-mass checks (report only, no plot yet) ---
    print_modal_extraction_completeness(modal_df)
    print_direction_coverage(modal_df)
    print_top5_modes_by_f_total(modal_df)
    print_ranking_agreement(modal_df)
    print_n_modes_extracted(modal_df)
    print_top_mode_vs_harmonic_peak(modal_df, harmonic_df)

    if len(modal_df) > 0:
        modal_df.to_csv(MODAL_EFFMASS_CSV_PATH, index=False)
        print(f"\nPer-mode effective-mass table saved to: {MODAL_EFFMASS_CSV_PATH}")
    else:
        print(f"\nNo modal data loaded - {MODAL_EFFMASS_CSV_PATH} was not written.")

    all_models = set(modal_df["model"].unique()) | set(harmonic_df["model"].unique())
    model_colors = assign_scenario_colors(all_models)   # ALWAYS the full set - keeps a model's colour fixed regardless of MODEL_INCLUDE

    # --- MODEL_INCLUDE: resolve which models are DRAWN (never affects
    # modal_df/harmonic_df themselves, which stay the full tables for
    # modal_effmass.csv/frequency_deviation.csv) ---
    visible_models_list = resolve_model_filter(sorted(all_models), MODEL_INCLUDE)
    if MODEL_INCLUDE is not None and len(visible_models_list) == 0:
        print("MODEL_INCLUDE matched no model at all - nothing to draw. No figures written.")
        return
    modal_df_drawing = modal_df if MODEL_INCLUDE is None else modal_df[modal_df["model"].isin(visible_models_list)]
    harmonic_df_drawing = harmonic_df if MODEL_INCLUDE is None else harmonic_df[harmonic_df["model"].isin(visible_models_list)]
    visible_models = None if MODEL_INCLUDE is None else visible_models_list

    # One consistent grey-shade mapping per year, shared across every
    # figure/tree so a given year always looks the same.
    all_years = set(freq_elasto["year"].unique()) | set(freq_median["year"].unique())
    year_colors = get_year_grey_shades(all_years)

    modal_trees = set(modal_df["tree"].unique()) if len(modal_df) > 0 else set()
    harmonic_trees = set(harmonic_df["tree"].unique()) if len(harmonic_df) > 0 else set()
    all_trees = sorted(modal_trees | harmonic_trees)

    # PER-TREE figures (modal_comparison_*/harmonic_comparison_*/
    # harmonic_simple_*) are drawn only for trees that still have at least
    # one VISIBLE model - a tree whose only models were excluded by
    # MODEL_INCLUDE gets no file at all for these (there would be nothing
    # of its own to show). This does NOT affect modal_effmass.csv (written
    # from the full modal_df above, before this filter existed) or which
    # trees/models compute_deviation_table-style rows exist for.
    modal_trees_drawing = set(modal_df_drawing["tree"].unique()) if len(modal_df_drawing) > 0 else set()
    harmonic_trees_drawing = set(harmonic_df_drawing["tree"].unique()) if len(harmonic_df_drawing) > 0 else set()
    dropped_by_filter = sorted((modal_trees | harmonic_trees) - (modal_trees_drawing | harmonic_trees_drawing))
    if dropped_by_filter:
        print(f"MODEL_INCLUDE dropped {len(dropped_by_filter)} tree(s) with no visible "
              f"model left (no per-tree figure written for them): {dropped_by_filter}")

    # CHANGE 5: this loop only ever visits trees that have FEM data (above).
    # Report, but do not change, which experimental trees get left out as a
    # result - matching what ONLY_TREES_WITH_FEM already prints in
    # plot_stiffness_comparison.py.
    experimental_trees = set(freq_elasto["tree"].unique()) | set(freq_median["tree"].unique())
    missing_trees = sorted(experimental_trees - set(all_trees))
    if missing_trees:
        print(f"{len(missing_trees)} experimental tree(s) have no FEM result "
              f"and are not drawn: {missing_trees}")

    # --- Parts 4-5: build the figures, in both "elasto" and "median" variants,
    # one pair per tree. Uses modal_df_drawing/harmonic_df_drawing (already
    # MODEL_INCLUDE-filtered, or identical to the full tables when that
    # filter is off) - modal_effmass.csv above was already written from the
    # FULL tables, so this filtering can never reach it. ---
    for tree in all_trees:
        if tree in modal_trees_drawing:
            plot_modal_figure(tree, modal_df_drawing, model_colors, freq_elasto, year_colors,
                               output_path=add_model_include_suffix(f"{OUTPUT_DIR}/modal_comparison_{tree}_elasto.png"))
            plot_modal_figure(tree, modal_df_drawing, model_colors, freq_median, year_colors,
                               output_path=add_model_include_suffix(f"{OUTPUT_DIR}/modal_comparison_{tree}_median.png"))
        if tree in harmonic_trees_drawing:
            plot_harmonic_figure(tree, harmonic_df_drawing, modal_df_drawing, freq_elasto, model_colors, year_colors,
                                  output_path=add_model_include_suffix(f"{OUTPUT_DIR}/harmonic_comparison_{tree}_elasto.png"))
            plot_harmonic_figure(tree, harmonic_df_drawing, modal_df_drawing, freq_median, model_colors, year_colors,
                                  output_path=add_model_include_suffix(f"{OUTPUT_DIR}/harmonic_comparison_{tree}_median.png"))
            if HARMONIC_USE_ABS:
                plot_harmonic_figure(tree, harmonic_df_drawing, modal_df_drawing, freq_elasto, model_colors, year_colors,
                                      output_path=add_model_include_suffix(f"{OUTPUT_DIR}/harmonic_comparison_{tree}_elasto_abs.png"),
                                      use_abs=True)
                plot_harmonic_figure(tree, harmonic_df_drawing, modal_df_drawing, freq_median, model_colors, year_colors,
                                      output_path=add_model_include_suffix(f"{OUTPUT_DIR}/harmonic_comparison_{tree}_median_abs.png"),
                                      use_abs=True)

            # VARIANT 2 - simplified single-axis harmonic figure for the
            # write-up: Top_ux only, absolute values, no modal/mass-fraction
            # axes (see plot_harmonic_simple_figure()'s own docstring). An
            # ADDITION alongside the full harmonic_comparison_* figures
            # above, produced in this same run - never a replacement.
            plot_harmonic_simple_figure(tree, harmonic_df_drawing, freq_elasto, model_colors, year_colors,
                                        output_path=add_model_include_suffix(f"{OUTPUT_DIR}/harmonic_simple_{tree}_elasto.png"))
            plot_harmonic_simple_figure(tree, harmonic_df_drawing, freq_median, model_colors, year_colors,
                                        output_path=add_model_include_suffix(f"{OUTPUT_DIR}/harmonic_simple_{tree}_median.png"))

    # --- Part 5b: frequency summary figures (CHANGE 4), one per experimental
    # variant. Both write into one combined deviation CSV below. tree_order
    # is deliberately the FULL all_trees, NOT MODEL_INCLUDE-filtered: this
    # is the one figure whose deviation_rows loop is also what feeds
    # FREQUENCY_DEVIATION_CSV_PATH, so trimming tree_order would silently
    # drop that tree's rows from the CSV too. visible_models instead gates
    # DRAWING per-model inside the function (see its own docstring) - a
    # tree with no visible model keeps its x-axis slot here (still showing
    # the grey experimental point, never a truly empty gap) rather than
    # being removed, which is the one place this script does not literally
    # drop the tree from the axis, precisely to protect the CSV.
    if dropped_by_filter:
        print(f"Note: frequency_summary figures keep all {len(all_trees)} tree(s) on the "
              f"x-axis (including {dropped_by_filter}) even though MODEL_INCLUDE leaves them "
              f"with no visible model - narrowing tree_order here would also drop those "
              f"trees' rows from frequency_deviation.csv, which must keep every row.")
    elasto_deviation = build_frequency_summary_figure(
        all_trees, modal_df, harmonic_df, model_colors, freq_elasto, "elasto",
        output_path=add_model_include_suffix(f"{OUTPUT_DIR}/frequency_summary_elasto.png"),
        visible_models=visible_models,
    )
    median_deviation = build_frequency_summary_figure(
        all_trees, modal_df, harmonic_df, model_colors, freq_median, "median",
        output_path=add_model_include_suffix(f"{OUTPUT_DIR}/frequency_summary_median.png"),
        visible_models=visible_models,
    )

    # VARIANT 3 - harmonic-only frequency summary for the write-up:
    # SUMMARY_SHOW_MODAL behaviour (already built and tested above) applied
    # PERMANENTLY via show_modal=False, regardless of the module-level
    # SUMMARY_SHOW_MODAL setting - an ADDITION alongside the plain
    # frequency_summary_{variant}.png figures just built, never a
    # replacement for them. The returned deviation tables are DISCARDED
    # (not concatenated below) - they contain the exact same rows
    # (deviation_rows is built unconditionally, see this function's own
    # docstring), so folding them in would duplicate every row in
    # frequency_deviation.csv. This call exists to render the figure, not
    # to add data to the CSV.
    _ = build_frequency_summary_figure(
        all_trees, modal_df, harmonic_df, model_colors, freq_elasto, "elasto",
        output_path=add_model_include_suffix(f"{OUTPUT_DIR}/frequency_summary_harmonic_elasto.png"),
        show_modal=False, visible_models=visible_models,
    )
    _ = build_frequency_summary_figure(
        all_trees, modal_df, harmonic_df, model_colors, freq_median, "median",
        output_path=add_model_include_suffix(f"{OUTPUT_DIR}/frequency_summary_harmonic_median.png"),
        show_modal=False, visible_models=visible_models,
    )

    frequency_deviation_table = pd.concat([elasto_deviation, median_deviation], ignore_index=True)
    frequency_deviation_table.to_csv(FREQUENCY_DEVIATION_CSV_PATH, index=False)
    print(f"\nFrequency deviation table saved to: {FREQUENCY_DEVIATION_CSV_PATH}")

    # --- Part 6: console summary ----------------------------------------------
    print(f"\nSummary: {len(modal_trees)} tree(s) with modal data "
          f"(models: {sorted(modal_df['model'].unique()) if len(modal_df) > 0 else []}), "
          f"{len(harmonic_trees)} tree(s) with harmonic data "
          f"(models: {sorted(harmonic_df['model'].unique()) if len(harmonic_df) > 0 else []}).")

    print("\nExperimental pulls per tree/year - Elasto(90) only:")
    if len(freq_elasto) > 0:
        print(freq_elasto.groupby(["tree", "year"])["first_frequency"].count()
              .rename("n_pulls").to_string())
    else:
        print("  (no data)")

    print("\nExperimental pulls per tree/year - median across sensors:")
    if len(freq_median) > 0:
        print(freq_median.groupby(["tree", "year"])["first_frequency"].count()
              .rename("n_pulls").to_string())
    else:
        print("  (no data)")


if __name__ == "__main__":
    main()
