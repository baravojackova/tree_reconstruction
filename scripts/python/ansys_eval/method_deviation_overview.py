"""
method_deviation_overview.py
------------------------------
Four write-up figures, deliberately split by whether a field reference
exists for the quantity being shown - mixing the two would suggest that
things like volume are validated when they are not (the beeches have no
destructive volume reference):

  WITH a field reference  : displacement stiffness, harmonic peak frequency,
                            trunk taper error         -> FIGURE A
  WITHOUT a field reference: volume (total/trunk/branch), model mass
                                                        -> FIGURE B
  Open-question regressions, using quantities from both groups:
    crown-mass difference vs frequency deviation       -> FIGURE C
    taper error vs stiffness deviation                 -> FIGURE D

Rotation stiffness is a special case even though it HAS a field reference:
the root-ball height was tuned on these same trees to match the static
pulling test, so its agreement measures the CALIBRATION, not the model. It
is shown in FIGURE A as a visually distinct, separately-labelled panel -
never presented alongside the other three as validation.

THE _no_acel MODELS (added after this script's first draft - see the
CONSISTENCY CHECK printed first thing in main()): scenario 3/4
("..._no acel" / "..._no_acel") are the SAME two reconstructions
(adtree_r5mm, treeqsm_pd002) with gravity/self-weight switched off, not
additional reconstruction variants. Two consequences, both load-bearing for
this script's structure:
  - the two spellings ("_no acel" with a space, "_no_acel" with an
    underscore) are NOT normalised anywhere a model name is DISPLAYED -
    only an internal matching helper (is_no_acel_model()/base_model_name())
    treats them the same, purely to pair a run with its twin. The
    consistency check prints exactly what is in the data.
  - FIGURE A gets a 5th panel, "self-weight contribution", built from the
    (model, no_acel-twin) pairs; every OTHER panel in every figure excludes
    the _no_acel runs entirely, since their geometry and mass are identical
    to their twins (confirmed empirically - see the consistency check) and
    would only add duplicate points on the volume/mass axes while carrying
    a different frequency.

FAMILY MAPPING (model name -> which volume/taper column applies), grounded
in exact numeric matches, not assumption - see check_family_mapping_evidence():
  "adtree_r5mm" (and its _no_acel twin) traces back to the calibrated AdTree
  geometry -> volume_summary_settled_*.csv column "adtree_calibrated" and
  trunk_taper_error_summary.csv's "adqsm_*" columns (the taper analysis was
  run on the AdQSM-derived cylinder model, one level below AdTree
  calibration).
  "treeqsm_pd002" (and its twin) -> volume column "treeqsm_simplified" and
  taper columns "treeqsm_*". Evidence: for B21_S12, adtree_calibrated trunk
  volume (0.777 m3) and treeqsm_simplified trunk volume (0.949 m3) match
  volume_results.csv's "AdTree calibrated r5mm (AdQSM 100)" (0.777740 m3)
  and "TreeQSM mine (man_pd04-002-05_mo9_sr010_ri1, Simplified (no
  islands))" (0.949298 m3) rows to 3 decimal places, and the treeqsm_tag
  "man_pd04-002-05_mo9_sr010_ri1" is exactly the geom_*.txt filename stem
  used for that ANSYS model.

NOMINAL WOOD DENSITY: 1184 kg/m3, read (not guessed) from
scripts/ansys/bk_kp.mac, line "densw=1184", applied via "mp,dens,1,densw"
(trunk material) and "mp,dens,2,densw" (branch material) - the actual
material density assigned to both wood materials in the FEM model.

INPUTS: located and reported (path found, or reported missing and that
figure/panel skipped - never substituted) at the top of main().

OUTPUTS: one PNG and one CSV per figure (A/B/C/D), all into ansys_eval/out/.
This script only READS its inputs - it never writes to any of them.

How to use it:
  - Edit the values in the "CONFIGURATION" section below.
  - Run the script from a terminal:  python method_deviation_overview.py
  - Check the four PNG/CSV pairs and the console output.
"""

import re
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

WORK_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval"
OUT_FOLDER = Path(WORK_FOLDER) / "out"
DATA_FOLDER = Path(WORK_FOLDER) / "data"
REPO_ROOT = Path(WORK_FOLDER).resolve().parents[2]
CSV_FOLDER = REPO_ROOT / "scripts" / "python" / "csv"

STIFFNESS_DEVIATION_CSV_PATH = OUT_FOLDER / "stiffness_deviation.csv"
FREQUENCY_DEVIATION_CSV_PATH = OUT_FOLDER / "frequency_deviation.csv"
MODAL_EFFMASS_CSV_PATH = OUT_FOLDER / "modal_effmass.csv"
EXPERIMENTAL_CSV_PATH = DATA_FOLDER / "dynatree_data_merge_0926.csv"
VOLUME_TOTAL_CSV_PATH = CSV_FOLDER / "volume_summary_settled_total.csv"
VOLUME_TRUNK_CSV_PATH = CSV_FOLDER / "volume_summary_settled_trunk.csv"
VOLUME_BRANCH_CSV_PATH = CSV_FOLDER / "volume_summary_settled_branch.csv"

# The taper CSV is located by SEARCHING, not assumed - see find_taper_csv().
TAPER_SEARCH_FOLDER = CSV_FOLDER
TAPER_CSV_NAME = "trunk_taper_error_summary.csv"

FIGURE_A_PNG_PATH = OUT_FOLDER / "method_deviation_vs_measurement.png"
FIGURE_A_CSV_PATH = OUT_FOLDER / "method_deviation_vs_measurement.csv"
FIGURE_B_PNG_PATH = OUT_FOLDER / "method_volume_disagreement.png"
FIGURE_B_CSV_PATH = OUT_FOLDER / "method_volume_disagreement.csv"
FIGURE_C_PNG_PATH = OUT_FOLDER / "frequency_vs_branch_volume.png"
FIGURE_C_CSV_PATH = OUT_FOLDER / "frequency_vs_branch_volume.csv"
FIGURE_D_PNG_PATH = OUT_FOLDER / "taper_vs_stiffness.png"
FIGURE_D_CSV_PATH = OUT_FOLDER / "taper_vs_stiffness.csv"

# None = every model is DRAWN. A list = only these models are drawn -
# matched against the "model" column (not the scenario number), tolerant
# of the space/underscore inconsistency in the no-acel suffix (see
# normalize_model_name_for_matching()) but never rewriting a displayed
# name - that inconsistency is being fixed in ANSYS, not here.
#
# This narrows base_models (see check_model_naming_consistency()) further -
# it can never remove a row from stiffness_deviation.csv/
# frequency_deviation.csv/modal_effmass.csv, which this script only reads.
# It CAN make FIGURE A's self-weight panel have nothing to pair, if it
# leaves no base model at all - that panel is then skipped with a printed
# reason, not drawn empty.
#
# NOTE (Stage 7): the comment above USED to say this "can never make a
# _no_acel model appear in FIGURE A's panels 1-4" - that stopped being true
# in Stage 5 (DEVIATION_PANEL_VARIANT) and is now GRAVITY_VARIANT's job
# entirely (see below), applied structurally. Setting MODEL_INCLUDE to the
# two no_acel names (as this had been left) INTERSECTS WITH base_models -
# "every NON-no_acel model" - to the EMPTY SET, silently starving every
# panel in the whole script. Reset to None; use GRAVITY_VARIANT for gravity
# selection and MODEL_INCLUDE only to narrow by RECONSTRUCTION METHOD.
MODEL_INCLUDE = None

# Appended to every figure filename (PNG and its companion CSV) whenever
# MODEL_INCLUDE is set (e.g. method_deviation_vs_measurement_noacel.png),
# so a filtered run can never be confused with, or silently overwrite, the
# full one. main() refuses to run if MODEL_INCLUDE is set and this is empty.
MODEL_INCLUDE_SUFFIX = "no_acel"

# =========================================================================
# STALE-OUTPUT / PROVENANCE SAFEGUARDS (added Sep 2026)
# -------------------------------------------------------------------------
# Same incident, same fix, as plot_stiffness_comparison.py's own copy of
# this block (see there for the full story).
# =========================================================================

# True = stamp a small grey footer on every figure - see
# add_provenance_footer(). False only for a publication figure.
FIGURE_PROVENANCE_FOOTER = True

# True = refuse to run if any raw ANSYS file is newer than
# STIFFNESS_DEVIATION_CSV_PATH / FREQUENCY_DEVIATION_CSV_PATH / MODAL_EFFMASS_CSV_PATH -
# i.e. plot_stiffness_comparison.py/plot_modal_harmonic.py have not been
# re-run since new ANSYS results landed, so this script's inputs are stale.
ABORT_ON_STALE_INPUT = True

# Folder holding the raw ANSYS output files, needed ONLY for the freshness
# check above - this script otherwise never reads raw ANSYS files directly.
ANSYS_RAW_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\ansys"
ANSYS_RAW_PATTERN_FOR_FRESHNESS = "out_*_bk*_scenario*.txt"

# What each ANSYS scenario number actually is - identical to the dict in
# load_static_scenarios.py, plot_modal_harmonic.py, and
# plot_frequency_spectra.py, copied (not imported) per the ansys_eval
# self-contained convention. Kept in sync across all four files by
# check_scenario_dicts_in_sync.py (Stage 8 Part C4) - run it after editing
# this dict. This copy is no longer used for row selection in this script
# (see STRICT_SELECTION's own comment for why the SCENARIOS constant that
# used to consult it was retired), only for tools/scripts that still want
# a scenario-number -> model-name mapping.
SCENARIO_MODELS = {
    1: "adtree_r5mm",
    2: "treeqsm_pd002",
    3: "adtree_r5mm_no_acel",
    4: "treeqsm_pd002_no_acel",
    999: "adtree_r5_pilot",
}

# --- WHICH TREES TO PLOT -------------------------------------------------
# "all"  = every tree found in the loaded deviation tables (old behaviour)
# or a list, mixed forms allowed:  [7, 12]  /  ["B21_S07", "S12"]  /  ["7"]
# NOTE: these are TREE NUMBERS, not the bk<N> index from the ANSYS filenames -
# see normalize_tree_selector()'s docstring.
TREES = "all"

# If True, raises an error (after printing the selection report) when a
# requested tree has no data at all.
#
# NOTE (Stage 8): this used to also gate a SCENARIOS constant, which chose
# which scenario numbers to REPORT on (it never filtered rows - see
# GRAVITY_VARIANT's own comment for why scenarios 3/4 need BOTH gravity
# halves together for the self-weight panel). SCENARIOS was retired -
# fully redundant with GRAVITY_VARIANT (model family) x MODEL_INCLUDE
# (reconstruction method), see Stage 7 Part A. FIGURES B/C/D now honour
# GRAVITY_VARIANT directly (Stage 8 Part A) - they are NOT "by design"
# gravity-on-only any more; only FIGURE A's self-weight panel still always
# needs both variants, by definition, regardless of GRAVITY_VARIANT.
STRICT_SELECTION = True

# Same filters plot_stiffness_comparison.py uses, needed here to recompute
# the measurement's OWN spread from the same population that pipeline's
# exp_median came from.
DAYS_TO_PROCESS = ["2021-03-22", "2022-04-05"]
MEASUREMENT = "M01"
TYPE_FILTER = "normal"
REDUCTION_LEVEL_FILTER = "none"
TREE_ID_PREFIX = "B21_S"

# M01 (a single static pull) never records first_frequency - only the
# dynamic "pull and release" tests (M02+) do (M01 rows have
# first_frequency = NaN structurally, not as missing data). Filtering the
# frequency measurement band to M01, as the stiffness quantities are,
# would silently zero out every tree's band. Same convention as
# plot_modal_harmonic.py's own FREQUENCY_MEASUREMENT_FILTER = None.
FREQUENCY_MEASUREMENT_FILTER = None

# Below this many pooled measured values, the 25th-75th percentile band is
# not trustworthy - only the zero (measurement-median) line is drawn, and
# why is printed (see compute_measurement_band()).
MIN_MEASURED_VALUES_FOR_BAND = 3

# Of the two rotation FEM columns available (slope_rotation_base / base2),
# FIGURE A's rotation panel shows only this one, to keep one point per
# tree/model - "base" (root-node rotation sum), not "base2".
ROTATION_FEM_COLUMN = "slope_rotation_base"

# frequency_deviation.csv has two experimental variants (elasto/median).
# Every figure here that needs ONE frequency-deviation value per tree/model
# uses this one - "median" (pooled across whatever sensors recorded each
# pull), the more robust of the two, matching this script's own measurement
# band (which also pools across sensors).
FREQUENCY_VARIANT_FOR_ANALYSIS = "median"

# trunk_taper_error_summary.csv has an "_all" and an "_above_base" variant
# per family. "_above_base" excludes the root-flare region (buttressing),
# which is usually excluded from taper fitting - used throughout.
TAPER_ERROR_COLUMN_SUFFIX = "err_above_base_pct"

# =========================================================================
# THE GRAVITY-VARIANT STRUCTURAL LOCK (added Stage 7)
# -------------------------------------------------------------------------
# GRAVITY_VARIANT is now the SINGLE SOURCE OF TRUTH for which gravity
# variant every figure in this script shows, applied ONCE in main() right
# after the input CSVs are loaded, to stiffness_df_variant/
# frequency_df_variant - the dataframes FIGURE A's panels 1-4 actually draw
# from. Panel 5 (self-weight contribution) is the ONE exception: it needs
# BOTH halves of every pair by definition, and reads the separately-named,
# never-filtered stiffness_df_all_variants/frequency_df_all_variants
# instead - see main()'s own comment for the full wiring, and
# assert_gravity_variant() for the check that fails loudly if a function
# ever sees the wrong one.
#
# WHY this replaced Stage 5's per-figure DEVIATION_PANEL_VARIANT: that
# constant only ever reached FIGURE A's 4 panels, and MODEL_INCLUDE was
# separately (and, for a while, WRONGLY - see MODEL_INCLUDE's own comment)
# being used as a stand-in gravity switch. Two knobs for one decision is
# exactly what broke: this script ran with base_models silently empty.
#
#   "no_acel"      = the gravity-off runs, scenarios 3 and 4 - the
#                    reference variant, and the new default.
#   "with_gravity" = the gravity-on runs, scenarios 1 and 2 (the OLD
#                    default - the 81 unsuffixed PNGs on disk mean this).
#   "both"         = both variants side by side in the same panel.
GRAVITY_VARIANT = "no_acel"

# --- DEVIATION_PANEL_VARIANT (Stage 5 name, kept as an ALIAS) -------------
# Stage 8 Part B2: this is no longer a second constant that has to AGREE
# with GRAVITY_VARIANT - it is not read anywhere in this file's actual
# logic any more (every function reads GRAVITY_VARIANT directly - see e.g.
# add_deviation_panel_variant_suffix()). The name is kept, set to
# GRAVITY_VARIANT's value, purely so any external code still written
# against the Stage 5 name (e.g. `import method_deviation_overview as m;
# m.DEVIATION_PANEL_VARIANT`) keeps working. Python has no live "constant
# alias" - this line copies GRAVITY_VARIANT's value ONCE, at import time -
# so if you ever need to change gravity variant PROGRAMMATICALLY after
# import (as opposed to editing this file and re-running, the normal way
# this script is used), set GRAVITY_VARIANT, not this name: only
# GRAVITY_VARIANT is actually consulted anywhere.
DEVIATION_PANEL_VARIANT = GRAVITY_VARIANT

# Model name (substring, case-insensitive) -> which volume_summary_settled_*
# column / trunk_taper_error_summary.csv family applies - see the module
# docstring's FAMILY MAPPING section for the numeric evidence.
MODEL_FAMILY_TO_VOLUME_COLUMN = {"adtree": "adtree_calibrated", "treeqsm": "treeqsm_simplified"}
MODEL_FAMILY_TO_TAPER_PREFIX = {"adtree": "adqsm", "treeqsm": "treeqsm"}

# kg/m3, read from scripts/ansys/bk_kp.mac ("densw=1184", mp,dens,1/2,densw)
# - the actual wood density assigned to both FEM materials. Reported
# alongside the computed M_total/volume ratio, never enforced or used to
# "correct" anything.
NOMINAL_WOOD_DENSITY_KG_M3 = 1184.0

# Fixed pastel colors, assigned in alphabetical order of the label -
# same convention (and same literal list) as every other ansys_eval script.
PASTEL_SCENARIO_COLORS = [
    "#94BBEA", "#F5B399", "#8DD7BC", "#F6D07F",
    "#F3BDD1", "#7FC17F", "#A49CD3", "#F1A4A3",
]

# --- Style config ------------------------------------------------------
TITLE_FONTSIZE = 13
SUPTITLE_FONTSIZE = 15
AXIS_LABEL_FONTSIZE = 10
TICK_LABEL_FONTSIZE = 9
LEGEND_FONTSIZE = 9
LEGEND_TITLE_FONTSIZE = 10
MARKER_SIZE = 9

LEGEND_BBOX_X = 1.02
LEGEND_Y_TOP = 1.0
LEGEND_Y_STEP = 0.28
RIGHT_MARGIN = 0.80

ZERO_LINE_COLOR = "0.3"
ZERO_LINE_WIDTH = 1.0
ONE_LINE_COLOR = "0.3"     # horizontal reference line at ratio=1.0 (FIGURE B)
ONE_LINE_WIDTH = 1.0

# FIGURE A's stiffness/frequency/rotation/self-weight panels mix small
# deviations (a few %, where the measurement band lives) with the
# known-bad treeqsm_pd002 parameter's extreme outliers (up to ~+11000%,
# see the module docstring). A linear axis puts every other point on top
# of the zero line; a plain log axis cannot show negative values at all.
# "symlog" keeps a linear region around zero (where the band comparison
# actually matters) and goes logarithmic beyond it - EVERY point stays
# visible and nothing is excluded or rescaled to "fix" the outlier.
DEVIATION_PANEL_YSCALE = "symlog"
DEVIATION_PANEL_LINTHRESH = 10.0   # +/-10%: linear inside, log outside

# Measurement-spread band (FIGURE A) - grey, per the grey=experiment
# convention used everywhere else in ansys_eval.
BAND_COLOR = "0.5"
BAND_ALPHA = 0.20
BAND_HALF_WIDTH = 0.32   # in tree-index units, i.e. band spans [tree-0.32, tree+0.32]

# Rotation panel (FIGURE A) - visually distinct background tint, since its
# agreement measures the calibration, not the model.
CALIBRATED_PANEL_BG_COLOR = "#FBEFD5"

# Horizontal per-model/method dodge, so two points at the same tree do not
# sit exactly on top of one another.
DODGE_WIDTH = 0.28

# Density reference line (FIGURE B)
DENSITY_LINE_COLOR = "0.3"
DENSITY_LINE_STYLE = "--"

# Regression lines (FIGURES C, D)
FIT_LINE_COLOR = "0.25"
FIT_LINE_WIDTH = 1.5
FIT_ANNOTATION_FONTSIZE = 9

# Per-point tree-number label next to each scatter point in FIGURE D (e.g.
# "07" for B21_S07) - small and grey so it reads as a point ID, not as data.
TAPER_POINT_LABEL_FONTSIZE = 7
TAPER_POINT_LABEL_COLOR = "0.35"


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
    "WORK_FOLDER", "OUT_FOLDER", "DATA_FOLDER", "REPO_ROOT", "CSV_FOLDER",
    "STIFFNESS_DEVIATION_CSV_PATH", "FREQUENCY_DEVIATION_CSV_PATH", "MODAL_EFFMASS_CSV_PATH",
    "EXPERIMENTAL_CSV_PATH", "VOLUME_TOTAL_CSV_PATH", "VOLUME_TRUNK_CSV_PATH",
    "VOLUME_BRANCH_CSV_PATH", "TAPER_SEARCH_FOLDER", "TAPER_CSV_NAME", "FIGURE_A_PNG_PATH",
    "FIGURE_A_CSV_PATH", "FIGURE_B_PNG_PATH", "FIGURE_B_CSV_PATH", "FIGURE_C_PNG_PATH",
    "FIGURE_C_CSV_PATH", "FIGURE_D_PNG_PATH", "FIGURE_D_CSV_PATH", "MODEL_INCLUDE",
    "MODEL_INCLUDE_SUFFIX", "FIGURE_PROVENANCE_FOOTER", "ABORT_ON_STALE_INPUT",
    "ANSYS_RAW_FOLDER", "ANSYS_RAW_PATTERN_FOR_FRESHNESS", "SCENARIO_MODELS", "TREES",
    "STRICT_SELECTION", "DAYS_TO_PROCESS", "MEASUREMENT", "TYPE_FILTER",
    "REDUCTION_LEVEL_FILTER", "TREE_ID_PREFIX", "FREQUENCY_MEASUREMENT_FILTER",
    "MIN_MEASURED_VALUES_FOR_BAND", "ROTATION_FEM_COLUMN", "FREQUENCY_VARIANT_FOR_ANALYSIS",
    "TAPER_ERROR_COLUMN_SUFFIX", "GRAVITY_VARIANT", "DEVIATION_PANEL_VARIANT",
    "MODEL_FAMILY_TO_VOLUME_COLUMN", "MODEL_FAMILY_TO_TAPER_PREFIX", "NOMINAL_WOOD_DENSITY_KG_M3",
    "PASTEL_SCENARIO_COLORS", "TITLE_FONTSIZE", "SUPTITLE_FONTSIZE", "AXIS_LABEL_FONTSIZE",
    "TICK_LABEL_FONTSIZE", "LEGEND_FONTSIZE", "LEGEND_TITLE_FONTSIZE", "MARKER_SIZE",
    "LEGEND_BBOX_X", "LEGEND_Y_TOP", "LEGEND_Y_STEP", "RIGHT_MARGIN", "ZERO_LINE_COLOR",
    "ZERO_LINE_WIDTH", "ONE_LINE_COLOR", "ONE_LINE_WIDTH", "DEVIATION_PANEL_YSCALE",
    "DEVIATION_PANEL_LINTHRESH", "BAND_COLOR", "BAND_ALPHA", "BAND_HALF_WIDTH",
    "CALIBRATED_PANEL_BG_COLOR", "DODGE_WIDTH", "DENSITY_LINE_COLOR", "DENSITY_LINE_STYLE",
    "FIT_LINE_COLOR", "FIT_LINE_WIDTH", "FIT_ANNOTATION_FONTSIZE",
    "TAPER_POINT_LABEL_FONTSIZE", "TAPER_POINT_LABEL_COLOR",
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


# --- generic helpers, copy-pasted from the other ansys_eval scripts,
# not imported - these scripts are deliberately standalone. -------------

def get_reduction_level(row) -> str:
    if row["reduction_20"]:
        return "20"
    elif row["reduction_20+20"]:
        return "20+20"
    elif row["reduction_40"]:
        return "40"
    else:
        return "none"


def add_reduction_level_column(df: pd.DataFrame) -> pd.DataFrame:
    if "reduction_level" in df.columns:
        return df
    df = df.copy()
    df["reduction_level"] = df.apply(get_reduction_level, axis=1)
    return df


def filter_column(data: pd.DataFrame, column_name: str, value):
    if value is None:
        return data
    if isinstance(value, (list, tuple, set)):
        return data[data[column_name].isin(value)]
    return data[data[column_name] == value]


def bk_id_to_canonical(bk_id: str) -> str:
    digits_only = "".join(character for character in str(bk_id) if character.isdigit())
    return f"{TREE_ID_PREFIX}{int(digits_only):02d}"


def tree_number_to_tree_id(tree_number: int) -> str:
    """
    Convert a plain tree number (e.g. 1, 7, 12) into the canonical project
    tree ID (e.g. "B21_S01", "B21_S07", "B21_S12"). Copied from
    load_static_scenarios.py, not imported.
    """
    return f"{TREE_ID_PREFIX}{tree_number:02d}"


def load_and_filter_experimental(csv_path, measurement_filter=MEASUREMENT, label: str = "") -> pd.DataFrame:
    """
    Load dynatree_data_merge_0926.csv and filter it down with the same
    days/type/reduction-level filters plot_stiffness_comparison.py uses,
    so the measurement spread computed here comes from the same population
    that pipeline's exp_median came from. measurement_filter defaults to
    "M01" (the stiffness pipeline's own filter) but is overridden to None
    for frequency - see FREQUENCY_MEASUREMENT_FILTER's own comment for why.
    """
    df = pd.read_csv(csv_path, low_memory=False)
    df = add_reduction_level_column(df)

    filtered = df[df["day"].isin(DAYS_TO_PROCESS)]
    filtered = filter_column(filtered, "measurement", measurement_filter)
    filtered = filtered[filtered["type"] == TYPE_FILTER]
    filtered = filter_column(filtered, "reduction_level", REDUCTION_LEVEL_FILTER)

    filtered = filtered.copy()
    filtered["tree"] = filtered["tree"].apply(bk_id_to_canonical)
    print(f"Experimental data (dynatree_data_merge_0926.csv{', ' + label if label else ''}): "
          f"{len(filtered)} row(s) after filtering (days={DAYS_TO_PROCESS}, "
          f"measurement={measurement_filter}, type={TYPE_FILTER}, "
          f"reduction_level={REDUCTION_LEVEL_FILTER})")
    return filtered


def assign_label_colors(labels) -> dict:
    """
    One fixed color per label (model OR reconstruction method - whichever
    this call is for), alphabetical assignment from PASTEL_SCENARIO_COLORS,
    same convention as every other ansys_eval script's assign_model_colors().
    """
    n_colors = len(PASTEL_SCENARIO_COLORS)
    return {
        label: PASTEL_SCENARIO_COLORS[i % n_colors]
        for i, label in enumerate(sorted(labels))
    }


def add_light_grid(ax) -> None:
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=0)


def compute_dodge_offsets(n_groups: int, group_width: float) -> list:
    if n_groups == 1:
        return [0.0]
    step = group_width / (n_groups - 1)
    return [-group_width / 2 + i * step for i in range(n_groups)]


def linear_fit(x: np.ndarray, y: np.ndarray) -> dict:
    """
    Ordinary least-squares fit y = slope*x + intercept, plus R-squared.
    Used as-is for FIGURE D (linear x), and on log10(x) for FIGURE C
    (so the fitted line is straight on that figure's log-x axis).
    """
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"slope": float(slope), "intercept": float(intercept), "r_squared": r_squared, "n": len(x)}


def fit_groups_for_gravity_variant(data: pd.DataFrame, x_col: str, y_col: str, variant: str,
                                   model_col: str = "model") -> list:
    """
    Return the list of (label, subset_df, fit_dict_or_None) regressions
    that should actually be DRAWN for FIGURES C/D, given GRAVITY_VARIANT.

    WHY this exists (Stage 8, Part A3): a model and its no_acel twin share
    an IDENTICAL x value here (taper error, phi - see D2's own comment)
    because gravity is a LOAD, not a geometry change. Pooling both halves
    into ONE fit would put two points at the same x for every tree,
    inflating n without adding an independent observation - a real
    statistical mistake, not just a display quirk.

      "no_acel" / "with_gravity" -> ONE fit, over all of `data` (already
        only one gravity variant, by construction - see main()).
      "both" -> TWO separate fits, one per gravity half, each with its own
        slope/intercept/R2/n - never pooled.

    x_col is whatever column the caller has already prepared (raw taper
    for FIGURE D, log10(phi) for FIGURE C - see linear_fit()'s own
    docstring for why C needs the log transform done first). A group with
    fewer than 2 points gets fit=None (printed/annotated as "no fit", not
    silently skipped from the figure - its POINTS are still drawn).
    """
    def _fit_or_none(subset: pd.DataFrame):
        if len(subset) < 2:
            return None
        return linear_fit(subset[x_col].to_numpy(), subset[y_col].to_numpy())

    if variant != "both":
        return [(None, data, _fit_or_none(data))]

    groups = []
    for label, is_no_acel in (("gravity ON", False), ("gravity OFF", True)):
        subset = data[data[model_col].apply(is_no_acel_model) == is_no_acel]
        groups.append((label, subset, _fit_or_none(subset)))
    return groups


# --- _no_acel matching (see module docstring) ---------------------------

NO_ACEL_PATTERN = re.compile(r"[ _]no[ _]acel$", re.IGNORECASE)


def is_no_acel_model(model_name: str) -> bool:
    return bool(NO_ACEL_PATTERN.search(str(model_name)))


def base_model_name(model_name: str) -> str:
    """
    Strip the no-acel suffix (either spelling) - for MATCHING a run to its
    twin only. Every printed/plotted label still uses the ORIGINAL string;
    this function's output is never displayed.
    """
    return NO_ACEL_PATTERN.sub("", str(model_name))


def normalize_model_name_for_matching(name: str) -> str:
    """
    For COMPARISON only (MODEL_INCLUDE matching) - never used to rewrite a
    displayed or written model name. Treats spaces and underscores as
    equivalent, same convention as is_no_acel_model()/base_model_name()
    above but generalised to the whole name, not just the no-acel suffix -
    so a MODEL_INCLUDE entry matches regardless of that inconsistency.
    """
    return re.sub(r"[ _]+", "_", str(name)).strip().lower()


def resolve_model_filter(all_models: list, include_filter) -> list:
    """
    Resolve MODEL_INCLUDE (None, or a list of wanted model names) against
    the actual distinct model names in the data (all_models, BASE and
    _no_acel alike), using normalize_model_name_for_matching() for
    tolerant comparison. The returned list always uses the ORIGINAL
    spelling from all_models - nothing is rewritten.

    Prints every distinct model found, and which of them MODEL_INCLUDE
    matched. An include_filter entry matching nothing gets a WARNING (not
    an error), naming it.

    include_filter=None -> returns all_models unchanged (every model drawn).
    """
    if include_filter is None:
        print("MODEL_INCLUDE is None - every model is eligible to be drawn.")
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


# Set once by main() (after TREES is resolved), read by
# add_model_include_suffix()/add_provenance_footer() below. Stays "" -
# meaning "no suffix" - whenever TREES is "all".
TREES_SUFFIX = ""
PROVENANCE_NEWEST_INPUT_MTIME = None
PROVENANCE_FILTERS_TEXT = ""


def add_model_include_suffix(path) -> str:
    """
    Insert MODEL_INCLUDE_SUFFIX and/or the active TREES suffix
    before the file extension whenever either is set (PNG or CSV alike -
    this script's own CSVs are companion "numbers behind the figure"
    outputs, not the authoritative stiffness_deviation.csv/
    frequency_deviation.csv/modal_effmass.csv it only reads, so they follow
    whatever this script draws, same as their figure). A no-op (path
    unchanged) when nothing is filtered.
    """
    path_str = str(path)
    parts = []
    if MODEL_INCLUDE is not None:
        parts.append(MODEL_INCLUDE_SUFFIX)
    if TREES_SUFFIX:
        parts.append(TREES_SUFFIX)
    if not parts:
        return path_str
    root, ext = path_str.rsplit(".", 1)
    return f"{root}_{'_'.join(parts)}.{ext}"


# Rejects the ANSYS row-position spelling outright (lowercase "bk<N>") -
# see normalize_tree_selector()'s own docstring. Copied from
# plot_stiffness_comparison.py.
BK_INDEX_REJECT_PATTERN = re.compile(r"^bk\d+$")


def normalize_tree_selector(value) -> str:
    """
    Convert one TREES entry into the canonical tree ID "B21_Snn". Accepts
    7, "7", "07", "S07", "s07", "B21_S07", "BK07" (uppercase only); rejects
    "bk3"-style (lowercase "bk" + digits). See
    plot_stiffness_comparison.py's copy of this function for the full WHY
    behind the BK/bk case distinction - copied here, not imported.
    """
    text = str(value).strip()

    if BK_INDEX_REJECT_PATTERN.match(text):
        raise ValueError(
            f"TREES entry {value!r} looks like an ANSYS row-position index "
            f"(lowercase 'bk<N>'). bk<N> is that tree's ROW POSITION in the "
            f"description file, NOT its tree number. Use the tree's own "
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
    actually present in the data. Copied from plot_stiffness_comparison.py.
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
    Filter df down to tree_sel["resolved"] - a no-op when is_all=True.
    """
    if tree_sel["is_all"]:
        return df
    return df[df[tree_col].isin(tree_sel["resolved"])]


def build_filter_suffix(tree_sel: dict) -> str:
    """
    Build the deterministic filename-suffix fragment for an active TREES
    selection. Copied from plot_stiffness_comparison.py.

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
    return True if STRICT_SELECTION should raise for this script.
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

    return bool(tree_sel["missing"]) or bool(missing_bits)


def newest_mtime(folder, pattern: str):
    """
    Return (newest_mtime_float, newest_path) among files in folder matching
    pattern, or (None, None). Copied from plot_stiffness_comparison.py.
    """
    matches = sorted(Path(folder).glob(pattern))
    if not matches:
        return None, None
    newest_path = max(matches, key=lambda p: p.stat().st_mtime)
    return newest_path.stat().st_mtime, newest_path


def check_input_freshness(script_label: str, raw_folder, raw_pattern: str,
                          csv_paths_being_read: dict, abort_if_stale: bool) -> float:
    """
    Compare the newest modification time among raw ANSYS files against
    every CSV in csv_paths_being_read; warn (and optionally abort) if any
    is older. Copied from plot_stiffness_comparison.py.
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
            continue
        csv_mtime = csv_path.stat().st_mtime
        if csv_mtime < newest_raw_mtime:
            stale_found = True
            print(f"*** WARNING: STALE INPUT - {label} ({csv_path.name}, "
                  f"{datetime.fromtimestamp(csv_mtime):%Y-%m-%d %H:%M:%S}) is OLDER than "
                  f"{newest_raw_path.name} ({datetime.fromtimestamp(newest_raw_mtime):%Y-%m-%d %H:%M:%S}). "
                  f"The upstream script that builds {csv_path.name} has not been re-run since "
                  f"that raw ANSYS file last changed. ***")

    if stale_found and abort_if_stale:
        raise SystemExit(
            f"{script_label}: aborting because ABORT_ON_STALE_INPUT=True and at least one "
            f"input CSV is older than a raw ANSYS file it depends on (see WARNING(s) above)."
        )
    return newest_raw_mtime


KNOWN_MODEL_SUFFIX_MARKERS = ("_noacel",)
TREES_SUFFIX_PATTERN = re.compile(r"_trees-[0-9-]+|_sc-[0-9-]+")


def find_filtered_output_files(out_folder, name_prefixes: tuple) -> list:
    """
    Return every FILE in out_folder whose name starts with one of
    name_prefixes and carries a filter-suffix marker. Copied from
    plot_stiffness_comparison.py.
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
    Strip every known filter-suffix marker to recover the UNFILTERED
    filename a filtered output was derived from. Copied from
    plot_stiffness_comparison.py.
    """
    filtered_path = Path(filtered_path)
    stem = filtered_path.stem
    for marker in KNOWN_MODEL_SUFFIX_MARKERS:
        stem = stem.replace(marker, "")
    stem = TREES_SUFFIX_PATTERN.sub("", stem)
    return filtered_path.with_name(stem + filtered_path.suffix)


def report_orphaned_filtered_outputs(out_folder, name_prefixes: tuple, script_label: str) -> None:
    """
    Print (NEVER touch) every filtered-output file whose unfiltered sibling
    is newer - i.e. it is stale. Copied from plot_stiffness_comparison.py.
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
    Stamp a small grey monospace footer along the bottom edge of fig.
    Copied from plot_stiffness_comparison.py. A no-op when
    FIGURE_PROVENANCE_FOOTER is False.

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
    Part A1) - every build_figure_*() function calls add_provenance_footer(
    fig, ...) before saving, which stashes fig.stage9_png_metadata; this
    passes that dict straight to savefig(metadata=...), embedding it as
    real PNG metadata alongside the existing visible footer. A PNG with no
    such metadata predates Stage 9 and its configuration cannot be verified
    this way - see check_png_provenance.py.
    """
    plt.savefig(path, dpi=150, bbox_inches="tight", metadata=getattr(fig, "stage9_png_metadata", {}))


def check_model_naming_consistency(model_sets: dict) -> tuple:
    """
    Print every distinct 'model' value found in each input CSV, then check
    whether the no-acel suffix is spelled inconsistently across models
    (space vs underscore) - reported, NOT normalised, so the naming can be
    fixed at the source.

    Returns (base_models, all_models): base_models is the sorted list of
    "base" (non-no_acel) model names found anywhere - what every
    panel/figure other than FIGURE A's self-weight panel plots by design;
    all_models is every distinct name found (base AND _no_acel), the
    universe MODEL_INCLUDE is matched against.
    """
    print("\nDistinct 'model' values found, per input file:")
    all_models = set()
    for label, models in model_sets.items():
        print(f"  {label}: {sorted(models)}")
        all_models |= set(models)

    no_acel_models = sorted(m for m in all_models if is_no_acel_model(m))
    print(f"\n'_no_acel' model names found across all files: {no_acel_models}")

    suffixes_found = sorted({NO_ACEL_PATTERN.search(m).group() for m in no_acel_models})
    print(f"Exact no-acel suffix spelling(s) found: {suffixes_found!r}")
    if len(suffixes_found) > 1:
        print(f"*** WARNING: the no-acel suffix is spelled INCONSISTENTLY across models "
              f"({suffixes_found!r}) - e.g. {no_acel_models}. This is a data-quality issue "
              f"in the upstream model-naming (the scenario -> model_name mapping); it is "
              f"NOT normalised here. Every literal string is still matched correctly for "
              f"pairing purposes (is_no_acel_model()/base_model_name() accept either "
              f"spelling), but the naming should be fixed at the source. ***")
    else:
        print("No inconsistency in the no-acel suffix spelling.")

    base_models = sorted(m for m in all_models if not is_no_acel_model(m))
    print(f"\nBase (gravity-on) models used for every panel except the self-weight one: {base_models}")
    return base_models, sorted(all_models)


def check_family_mapping_evidence() -> None:
    """
    Print the exact numeric evidence for MODEL_FAMILY_TO_VOLUME_COLUMN
    (see the module docstring) - run once, at startup, so the mapping is
    verifiable from this run's own console output, not just the docstring.
    """
    print("\nFamily-mapping evidence (model -> volume column), B21_S12 trunk volume:")
    try:
        trunk = pd.read_csv(VOLUME_TRUNK_CSV_PATH)
        row = trunk[trunk["tree"] == "B21_S12"].iloc[0]
        print(f"  adtree_calibrated = {row['adtree_calibrated']} m3 "
              f"(volume_results.csv 'AdTree calibrated r5mm (AdQSM 100)' trunk = 0.777740 m3)")
        print(f"  treeqsm_simplified = {row['treeqsm_simplified']} m3 "
              f"(volume_results.csv 'TreeQSM mine (man_pd04-002-05_mo9_sr010_ri1, "
              f"Simplified (no islands))' trunk = 0.949298 m3)")
    except Exception as exc:
        print(f"  could not verify (volume_summary_settled_trunk.csv not usable: {exc})")


def model_to_family(model_name: str):
    """'adtree' or 'treeqsm' (substring match on the model name), else None."""
    name_lower = model_name.lower()
    if "adtree" in name_lower:
        return "adtree"
    if "treeqsm" in name_lower:
        return "treeqsm"
    return None


def resolve_deviation_panel_models(base_models: list, all_models_found, variant: str) -> list:
    """
    Resolve DEVIATION_PANEL_VARIANT into the actual "model" names panels
    1-4 of FIGURE A should filter/loop on (see that constant's own
    CONFIGURATION comment for the three settings and why ALL FOUR panels
    always switch together).

      "with_gravity" -> base_models, unchanged (old behaviour).
      "no_acel"      -> each base model's no_acel TWIN, found via
                        is_no_acel_model()/base_model_name() so either
                        suffix spelling ("_no acel" / "_no_acel") pairs
                        correctly - never the base models themselves.
      "both"         -> base_models followed by their no_acel twins, so a
                        tree's slot draws both variants side by side.

    A base model with no no_acel twin actually present in the data is
    reported with a WARNING and left out of the "no_acel"/"both" result for
    that model - never silently guessed or substituted.
    """
    if variant == "with_gravity":
        return list(base_models)

    all_models_found = list(all_models_found)
    twins = []
    for base_model in base_models:
        twin = next((m for m in all_models_found
                    if is_no_acel_model(m) and base_model_name(m) == base_model), None)
        if twin is None:
            print(f"  WARNING: DEVIATION_PANEL_VARIANT={variant!r} but '{base_model}' has no "
                  f"no_acel twin in the loaded data - it is left out of panels 1-4 this run.")
            continue
        twins.append(twin)

    if variant == "no_acel":
        return twins
    if variant == "both":
        return list(base_models) + twins

    raise ValueError(
        f"Unknown DEVIATION_PANEL_VARIANT: {variant!r} - use 'with_gravity', 'no_acel', or 'both'."
    )


def gravity_variant_filename_suffix(variant: str) -> str:
    """
    Filename fragment for GRAVITY_VARIANT - see
    plot_stiffness_comparison.py's copy of this function for the full WHY
    ("with_gravity" is the one setting that keeps filenames unchanged).
    """
    return {"no_acel": "gravity-noacel", "both": "gravity-both"}.get(variant, "")


def resolve_input_path_for_gravity_variant(base_path, variant: str) -> Path:
    """
    stiffness_deviation.csv/frequency_deviation.csv/modal_effmass.csv are
    WRITTEN by plot_stiffness_comparison.py/plot_modal_harmonic.py under a
    GRAVITY_VARIANT-suffixed name whenever THEIR OWN GRAVITY_VARIANT is not
    "with_gravity" (e.g. "stiffness_deviation_gravity-noacel.csv" - see
    gravity_variant_filename_suffix() in those scripts, copied here as the
    same function above). This script must READ whichever file those
    scripts would have written under the SAME GRAVITY_VARIANT setting - the
    three scripts' GRAVITY_VARIANT is meant to be run in agreement - or
    every run here looks for a file that was never written under the new
    "no_acel" default, and the freshness guard (or a bare file-not-found)
    fires for the wrong reason.
    """
    base_path = Path(base_path)
    suffix = gravity_variant_filename_suffix(variant)
    if not suffix:
        return base_path
    return base_path.with_name(f"{base_path.stem}_{suffix}{base_path.suffix}")


def filter_by_gravity_variant(df: pd.DataFrame, variant: str, model_col: str = "model") -> pd.DataFrame:
    """
    THE STRUCTURAL LOCK (see GRAVITY_VARIANT's own comment). Keep only rows
    whose model matches the requested gravity variant. Called ONCE, in
    main(), to build stiffness_df_variant/frequency_df_variant - the
    dataframes FIGURE A's panels 1-4 draw from. Copied from
    plot_stiffness_comparison.py, not imported.
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
    Loudly fail if df contains a model from the WRONG gravity variant - see
    plot_stiffness_comparison.py's copy of this function for the full WHY.
    expected_variant="both" always passes - used here for the self-weight
    panel's inputs, which are SUPPOSED to contain both variants.
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


def deviation_panel_variant_title_suffix(variant: str) -> str:
    """
    Short, explicit text appended to every one of panels 1-4's titles so a
    reader looking only at the image (not the legend or this source file)
    can see which gravity variant is drawn - see DEVIATION_PANEL_VARIANT's
    own comment for why this must never be left to inference.

    Starts with its OWN "\\n" (not a " - " suffix on the existing first
    line): panels 1-4 sit in a 5-across, ~6"-wide row, and appending this
    text inline made the first title line long enough to visually collide
    with the neighbouring panel's title - a separate wrapped line reads
    cleanly at that width.
    """
    if variant == "no_acel":
        return "\ngravity OFF (scenarios 3, 4)"
    if variant == "both":
        return "\ngravity ON + OFF (scenarios 1/2 + 3/4)"
    return ""


def add_deviation_panel_variant_suffix(path) -> str:
    """
    Append a deterministic filename suffix for FIGURES A/B/C/D, whenever
    GRAVITY_VARIANT is not the shipped default "with_gravity" - so a
    gravity-off (or both-variant) figure can never overwrite the gravity-on
    file already on disk (81 existing unsuffixed PNGs all mean "gravity-on"
    - see GRAVITY_VARIANT's own comment). A no-op (path unchanged) for the
    default.

    Reads GRAVITY_VARIANT directly (Stage 8 Part B2) - NOT
    DEVIATION_PANEL_VARIANT, even though that name is kept as an alias (see
    its own comment) - so this always reflects whichever value was actually
    set, with no risk of the two drifting apart at runtime.
    """
    path_str = str(path)
    suffix = {"no_acel": "panels-noacel", "both": "panels-both"}.get(GRAVITY_VARIANT)
    if not suffix:
        return path_str
    root, ext = path_str.rsplit(".", 1)
    return f"{root}_{suffix}.{ext}"


def find_input(path: Path, label: str) -> bool:
    exists = path.is_file()
    status = "FOUND" if exists else "MISSING"
    print(f"  [{status}] {label}: {path}")
    return exists


def find_taper_csv(folder: Path) -> Path:
    candidates = sorted(folder.glob("*taper*.csv"))
    print(f"\nSearching {folder} for '*taper*.csv':")
    for c in candidates:
        print(f"  found: {c.name}")
    if not candidates:
        print("  (no candidates found)")
        return None

    chosen = folder / TAPER_CSV_NAME
    if chosen.is_file():
        print(f"Using: {chosen.name} - one row per tree, adqsm_err_*/treeqsm_err_* "
              f"columns, matching 'the trunk taper error summary CSV'.")
        return chosen
    print(f"WARNING: expected file '{TAPER_CSV_NAME}' not among the candidates above - "
          f"the taper panel (FIGURE A) and FIGURE D will be skipped.")
    return None


def compute_measurement_band(experimental: pd.DataFrame, tree: str, columns: list) -> dict:
    """
    Pool the raw per-sensor, per-pull values in `columns` for one tree
    (already filtered by the stiffness-pipeline filters) and return the
    measurement's own spread as a dict: n, median, p25_pct, p75_pct - the
    last two expressed as PERCENT OF THE MEDIAN
    (100*(p25-median)/median, 100*(p75-median)/median), i.e. directly on
    the same axis as rel_dev_pct = 100*(FEM-median)/median, so a model
    point at 0 sits exactly on the measurement median and a point inside
    [p25_pct, p75_pct] is not distinguishable from measurement noise.

    "Spread across sensors and pulls": every individual (sensor, pull)
    reading in `columns` is pooled into ONE array (e.g. tree's slope_Pt3
    values and slope_Pt4 values concatenated) - not first collapsed to a
    per-pull combined value the way the stiffness pipeline's own
    "combined" exp_median is.
    """
    tree_rows = experimental[experimental["tree"] == tree]
    pooled = pd.concat([tree_rows[c].dropna() for c in columns], ignore_index=True)
    n = len(pooled)
    if n == 0:
        return {"n": 0, "median": float("nan"), "p25_pct": float("nan"), "p75_pct": float("nan")}
    median = float(pooled.median())
    if median == 0:
        return {"n": n, "median": median, "p25_pct": float("nan"), "p75_pct": float("nan")}
    p25 = float(pooled.quantile(0.25))
    p75 = float(pooled.quantile(0.75))
    return {
        "n": n, "median": median,
        "p25_pct": 100.0 * (p25 - median) / median,
        "p75_pct": 100.0 * (p75 - median) / median,
    }


def print_measurement_bands(band_by_tree_and_quantity: dict) -> None:
    print("\nMeasurement band per tree and quantity (25th-75th percentile, as % of "
          "the measurement's own median; n = pooled sensor x pull readings):")
    for quantity, band_by_tree in band_by_tree_and_quantity.items():
        print(f"  {quantity}:")
        for tree, band in band_by_tree.items():
            if band["n"] == 0:
                print(f"    {tree}: n=0 - no measurement data, no reference line meaningful.")
            elif band["n"] < MIN_MEASURED_VALUES_FOR_BAND:
                print(f"    {tree}: n={band['n']} (< {MIN_MEASURED_VALUES_FOR_BAND}) - "
                      f"band skipped, only the median (zero) line is drawn.")
            else:
                print(f"    {tree}: n={band['n']}, median={band['median']:.4g}, "
                      f"band=[{band['p25_pct']:+.1f}%, {band['p75_pct']:+.1f}%]")


def compute_self_weight_rows(value_df: pd.DataFrame, tree_order: list, base_models: list,
                             value_col: str, quantity_label: str) -> list:
    """
    value_df: columns tree, model, <value_col> - the RAW FEM value (not a
    deviation), one row per (tree, model).

    For each tree and base (gravity-on) model, finds that model's row and
    its no-acel twin's row (matched via base_model_name(), so either
    spelling of the suffix pairs correctly) on the SAME tree only, and
    computes:
        100 * (with_gravity - no_acel) / no_acel
    A tree missing either run is printed as missing and skipped - never
    paired with a different tree's run.
    """
    rows = []
    for tree in tree_order:
        tree_df = value_df[value_df["tree"] == tree]
        for base_model in base_models:
            with_gravity = tree_df[tree_df["model"] == base_model]
            twin = tree_df[tree_df["model"].apply(
                lambda m: is_no_acel_model(m) and base_model_name(m) == base_model)]

            if len(with_gravity) == 0 or len(twin) == 0:
                missing = []
                if len(with_gravity) == 0:
                    missing.append(f"'{base_model}' (with gravity)")
                if len(twin) == 0:
                    missing.append(f"'{base_model}'-no_acel twin")
                print(f"  MISSING for {tree} / {quantity_label}: {' and '.join(missing)} - "
                      f"self-weight contribution not computed for this tree/model.")
                continue

            with_val = with_gravity[value_col].iloc[0]
            no_acel_val = twin[value_col].iloc[0]
            if pd.isna(with_val) or pd.isna(no_acel_val) or no_acel_val == 0:
                print(f"  {tree} / {base_model} / {quantity_label}: value(s) missing/zero - "
                      f"cannot compute contribution.")
                continue

            pct = 100.0 * (with_val - no_acel_val) / no_acel_val
            rows.append({"tree": tree, "model": base_model, "quantity": quantity_label,
                        "with_gravity_value": with_val, "no_acel_value": no_acel_val,
                        "self_weight_pct": pct})
            print(f"  {tree} / {base_model} / {quantity_label}: with-gravity={with_val:.4g}, "
                  f"no_acel={no_acel_val:.4g} -> self-weight contribution = {pct:+.2f}%")
    return rows


# --- FIGURE A ------------------------------------------------------------

def draw_deviation_panel(ax, panel_df: pd.DataFrame, tree_order: list, model_colors: dict,
                         band_by_tree: dict, offset_by_model: dict, title: str, ylabel: str,
                         draw_band: bool, inside_band_log: list, quantity_label: str,
                         use_symlog: bool = True) -> None:
    """
    One FIGURE A panel: grey measurement-spread band (if draw_band and
    enough measured values), zero line, one dodge-offset point per
    tree/model. Appends (tree, model, quantity) to inside_band_log for
    every model point that falls inside its tree's band.

    use_symlog=True (the default - every panel except taper, whose values
    are small and unsigned) applies DEVIATION_PANEL_YSCALE so the
    known-bad treeqsm_pd002 outliers (up to ~+11000%) do not squeeze every
    other point onto the zero line - see that constant's own comment.
    """
    add_light_grid(ax)
    ax.axhline(0.0, color=ZERO_LINE_COLOR, linewidth=ZERO_LINE_WIDTH, zorder=2)
    if use_symlog:
        ax.set_yscale(DEVIATION_PANEL_YSCALE, linthresh=DEVIATION_PANEL_LINTHRESH)

    if draw_band:
        for tree_index, tree in enumerate(tree_order):
            band = band_by_tree.get(tree)
            if band is None or band["n"] < MIN_MEASURED_VALUES_FOR_BAND or pd.isna(band["p25_pct"]):
                continue
            ax.fill_between([tree_index - BAND_HALF_WIDTH, tree_index + BAND_HALF_WIDTH],
                            band["p25_pct"], band["p75_pct"], color=BAND_COLOR,
                            alpha=BAND_ALPHA, zorder=1)

    for _, row in panel_df.iterrows():
        if row["tree"] not in tree_order:
            continue
        tree_index = tree_order.index(row["tree"])
        color = model_colors.get(row["model"])
        if color is None:
            continue
        x = tree_index + offset_by_model.get(row["model"], 0.0)
        ax.scatter(x, row["rel_dev_pct"], color=color, s=MARKER_SIZE ** 2, zorder=5, edgecolor="none")

        band = band_by_tree.get(row["tree"]) if draw_band else None
        if band and band["n"] >= MIN_MEASURED_VALUES_FOR_BAND and pd.notna(band["p25_pct"]):
            if band["p25_pct"] <= row["rel_dev_pct"] <= band["p75_pct"]:
                inside_band_log.append((row["tree"], row["model"], quantity_label))

    ax.set_xticks(range(len(tree_order)))
    ax.set_xticklabels(tree_order, rotation=45, ha="right")
    ax.set_title(title, fontsize=TITLE_FONTSIZE)
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)


def build_figure_a(stiffness_df_variant, frequency_df_variant, stiffness_df_all_variants,
                   frequency_df_all_variants, taper_df, experimental_df, experimental_freq_df,
                   base_models, tree_order, model_colors: dict) -> pd.DataFrame:
    # model_colors is passed in, computed ONCE in main() from every model
    # found in the loaded data - see that call's own comment for why this
    # is no longer assign_label_colors(base_models) here (a model's colour
    # used to shift when MODEL_INCLUDE narrowed base_models).
    #
    # THE STRUCTURAL LOCK'S ASSERTIONS (see assert_gravity_variant()): the
    # *_variant tables must contain ONLY GRAVITY_VARIANT's models (panels
    # 1-4); the *_all_variants tables are ALLOWED both (panel 5 needs both
    # halves of every pair by definition) - "both" here means "no
    # restriction", not "must contain both".
    assert_gravity_variant(stiffness_df_variant, GRAVITY_VARIANT, "build_figure_a (stiffness_df_variant)")
    assert_gravity_variant(frequency_df_variant, GRAVITY_VARIANT, "build_figure_a (frequency_df_variant)")
    assert_gravity_variant(stiffness_df_all_variants, "both", "build_figure_a (stiffness_df_all_variants)")
    assert_gravity_variant(frequency_df_all_variants, "both", "build_figure_a (frequency_df_all_variants)")

    # offset_by_model (base_models only) is kept for PANEL 5, which always
    # pairs a base model against its no_acel twin regardless of
    # GRAVITY_VARIANT - see that panel's own comment below.
    offset_by_model = dict(zip(sorted(base_models), compute_dodge_offsets(len(base_models), DODGE_WIDTH)))

    # --- panels 1-4's model set: read straight off the ALREADY gravity-
    # filtered tables passed in (see main()'s own comment for where
    # stiffness_df_variant/frequency_df_variant come from) - no filtering
    # happens here any more, this just names what is already true. ---
    panel_models = sorted(set(stiffness_df_variant["model"].unique()) | set(frequency_df_variant["model"].unique()))
    offset_by_panel_model = dict(zip(panel_models, compute_dodge_offsets(len(panel_models), DODGE_WIDTH)))
    title_suffix = deviation_panel_variant_title_suffix(GRAVITY_VARIANT)
    print(f"\nGRAVITY_VARIANT={GRAVITY_VARIANT!r} - panels 1-4 draw: {panel_models}")

    csv_rows = []

    fig, axes = plt.subplots(1, 5, figsize=(30, 6))

    # --- panel 1: displacement stiffness ---------------------------------
    disp_df = stiffness_df_variant[(stiffness_df_variant["fem_column"] == "slope_displacement")
                           & (stiffness_df_variant["sensor"] == "combined")]
    disp_bands = {tree: compute_measurement_band(experimental_df, tree, ["slope_Pt3", "slope_Pt4"])
                 for tree in tree_order}
    inside_band_log = []
    draw_deviation_panel(axes[0], disp_df, tree_order, model_colors, disp_bands, offset_by_panel_model,
                         f"Displacement stiffness deviation{title_suffix}\n(n={len(tree_order)} trees)",
                         "Relative deviation (%)", True, inside_band_log, "displacement_stiffness")
    for _, row in disp_df.iterrows():
        band = disp_bands.get(row["tree"], {})
        csv_rows.append({"panel": "displacement_stiffness", "tree": row["tree"], "model": row["model"],
                         "rel_dev_pct": row["rel_dev_pct"], "band_n": band.get("n"),
                         "band_p25_pct": band.get("p25_pct"), "band_p75_pct": band.get("p75_pct")})

    # --- panel 2: harmonic peak frequency ---------------------------------
    freq_df = frequency_df_variant[(frequency_df_variant["quantity"] == "harmonic_peak")
                           & (frequency_df_variant["variant"] == FREQUENCY_VARIANT_FOR_ANALYSIS)]
    freq_bands = {tree: compute_measurement_band(experimental_freq_df, tree, ["first_frequency"])
                 for tree in tree_order}
    draw_deviation_panel(axes[1], freq_df, tree_order, model_colors, freq_bands, offset_by_panel_model,
                         f"Harmonic peak frequency deviation{title_suffix}\n(n={len(tree_order)} trees)",
                         "Relative deviation (%)", True, inside_band_log, "harmonic_peak_frequency")
    for _, row in freq_df.iterrows():
        band = freq_bands.get(row["tree"], {})
        csv_rows.append({"panel": "harmonic_peak_frequency", "tree": row["tree"], "model": row["model"],
                         "rel_dev_pct": row["rel_dev_pct"], "band_n": band.get("n"),
                         "band_p25_pct": band.get("p25_pct"), "band_p75_pct": band.get("p75_pct")})

    # --- panel 3: trunk taper error (unsigned, no field-pull measurement
    # distribution exists for this quantity - no band) --------------------
    print("\nTaper panel: no per-pull/sensor measurement distribution exists in "
          "dynatree_data_merge_0926.csv for trunk taper (it is derived from field "
          "tape-measured diameters at fixed heights, not repeated dynamic pulls) - "
          "no measurement band drawn for this panel.")
    taper_rows = []
    if taper_df is not None:
        for model in panel_models:
            # base_model_name() first: a no_acel twin's OWN family lookup
            # must land on the SAME taper column as its base model (the
            # taper error comes from the cylinder geometry, which no_acel
            # shares with its base - see DEVIATION_PANEL_VARIANT's comment
            # and the data-integrity check this script prints at startup).
            # model_to_family()'s substring match would already resolve a
            # no_acel name correctly on its own (it still contains "adtree"/
            # "treeqsm"), but going through the base name first is the
            # explicit, unambiguous way to say "same family as the twin".
            family = model_to_family(base_model_name(model))
            column = f"{MODEL_FAMILY_TO_TAPER_PREFIX.get(family, family)}_{TAPER_ERROR_COLUMN_SUFFIX}"
            if family is None or column not in taper_df.columns:
                print(f"  WARNING: no taper family/column for model '{model}' - skipped in taper panel.")
                continue
            for _, trow in taper_df.iterrows():
                if trow["tree"] not in tree_order:
                    continue
                taper_rows.append({"tree": trow["tree"], "model": model, "rel_dev_pct": trow[column]})
    taper_panel_df = pd.DataFrame(taper_rows)
    draw_deviation_panel(axes[2], taper_panel_df, tree_order, model_colors, {}, offset_by_panel_model,
                         f"Trunk taper error ({TAPER_ERROR_COLUMN_SUFFIX}){title_suffix}\n(n={len(tree_order)} trees)",
                         "Taper error (%, unsigned)", False, inside_band_log, "trunk_taper_error",
                         use_symlog=False)
    for _, row in taper_panel_df.iterrows():
        csv_rows.append({"panel": "trunk_taper_error", "tree": row["tree"], "model": row["model"],
                         "rel_dev_pct": row["rel_dev_pct"], "band_n": None,
                         "band_p25_pct": None, "band_p75_pct": None})

    # --- panel 4: rotation stiffness - CALIBRATED, tinted, separately titled
    # CALIBRATED wording applies unconditionally, for every
    # DEVIATION_PANEL_VARIANT: the root-ball height was tuned once, on the
    # geometry, to match the static pulling test (see the module docstring)
    # - that tuning is shared by a model and its no_acel twin (gravity is a
    # LOAD, not a geometry change), and nothing in this pipeline re-tunes it
    # per scenario, so the "not independent validation" caveat is exactly as
    # true for the gravity-off runs as for the gravity-on ones.
    rot_df = stiffness_df_variant[(stiffness_df_variant["fem_column"] == ROTATION_FEM_COLUMN)
                          & (stiffness_df_variant["sensor"] == "combined")]
    rot_bands = {tree: compute_measurement_band(experimental_df, tree,
                                                ["slope_root_stiff_cam", "slope_root_stiff_nocam"])
                for tree in tree_order}
    axes[3].set_facecolor(CALIBRATED_PANEL_BG_COLOR)
    draw_deviation_panel(axes[3], rot_df, tree_order, model_colors, rot_bands, offset_by_panel_model,
                         f"Rotation stiffness ({ROTATION_FEM_COLUMN})\nCALIBRATED - not independent "
                         f"validation{title_suffix} (n={len(tree_order)} trees)",
                         "Relative deviation (%)", True, inside_band_log, "rotation_stiffness_calibrated")
    for _, row in rot_df.iterrows():
        band = rot_bands.get(row["tree"], {})
        csv_rows.append({"panel": "rotation_stiffness_calibrated", "tree": row["tree"], "model": row["model"],
                         "rel_dev_pct": row["rel_dev_pct"], "band_n": band.get("n"),
                         "band_p25_pct": band.get("p25_pct"), "band_p75_pct": band.get("p75_pct")})

    # --- panel 5: self-weight contribution (with-gravity vs no_acel twin) ---
    # Needs at least one BASE model to pair against its no_acel twin. main()
    # already refuses to run at all if MODEL_INCLUDE leaves base_models
    # empty (every OTHER panel/figure here needs a base model too), so this
    # is mainly a defensive check for any future caller of build_figure_a()
    # with an already-empty base_models - printed and skipped, not drawn
    # empty, exactly as asked.
    if not base_models:
        print("\nSelf-weight panel skipped: no BASE (gravity-on) model available to pair "
              "with its no_acel twin (MODEL_INCLUDE left none).")
        axes[4].set_axis_off()
        axes[4].text(0.5, 0.5, "no base model to pair\n(MODEL_INCLUDE)", ha="center", va="center",
                    transform=axes[4].transAxes, fontsize=TITLE_FONTSIZE)
    else:
        print("\nSelf-weight contribution (100 * (with_gravity - no_acel) / no_acel):")
        # Deliberately *_all_variants here, NOT *_variant - this panel needs
        # BOTH gravity halves of every pair, which is exactly what
        # *_variant (GRAVITY_VARIANT-filtered) no longer contains.
        disp_raw = stiffness_df_all_variants[(stiffness_df_all_variants["fem_column"] == "slope_displacement")
                                & (stiffness_df_all_variants["sensor"] == "combined")][["tree", "model", "fem_value"]]
        freq_raw = frequency_df_all_variants[(frequency_df_all_variants["quantity"] == "harmonic_peak")
                                & (frequency_df_all_variants["variant"] == FREQUENCY_VARIANT_FOR_ANALYSIS)][["tree", "model", "fem_value_hz"]]
        print(" displacement stiffness:")
        sw_disp = compute_self_weight_rows(disp_raw, tree_order, base_models, "fem_value", "displacement_stiffness")
        print(" harmonic peak frequency:")
        sw_freq = compute_self_weight_rows(freq_raw, tree_order, base_models, "fem_value_hz", "harmonic_peak_frequency")

        add_light_grid(axes[4])
        axes[4].axhline(0.0, color=ZERO_LINE_COLOR, linewidth=ZERO_LINE_WIDTH, zorder=2)
        axes[4].set_yscale(DEVIATION_PANEL_YSCALE, linthresh=DEVIATION_PANEL_LINTHRESH)
        sw_marker = {"harmonic_peak_frequency": "o", "displacement_stiffness": "^"}
        for row in sw_disp + sw_freq:
            tree_index = tree_order.index(row["tree"])
            x = tree_index + offset_by_model.get(row["model"], 0.0)
            axes[4].scatter(x, row["self_weight_pct"], color=model_colors[row["model"]],
                            marker=sw_marker[row["quantity"]], s=MARKER_SIZE ** 2, zorder=5, edgecolor="none")
            csv_rows.append({"panel": "self_weight_contribution", "tree": row["tree"], "model": row["model"],
                             "rel_dev_pct": row["self_weight_pct"], "band_n": None,
                             "band_p25_pct": None, "band_p75_pct": None})
        axes[4].set_xticks(range(len(tree_order)))
        axes[4].set_xticklabels(tree_order, rotation=45, ha="right")
        axes[4].set_title(f"Self-weight contribution\n(gravity vs. no-acel twin, n={len(tree_order)} trees)",
                          fontsize=TITLE_FONTSIZE)
        axes[4].set_ylabel("Self-weight contribution (%)", fontsize=AXIS_LABEL_FONTSIZE)
        axes[4].tick_params(labelsize=TICK_LABEL_FONTSIZE)

    # --- legends, outside the last panel -------------------------------
    # "Model" lists only what panels 1-4 actually draw (panel_models) - NOT
    # every key of model_colors, so under DEVIATION_PANEL_VARIANT="no_acel"
    # this legend shows exactly the two no_acel names, in their own colours
    # (pink/yellow), which is the at-a-glance confirmation the variant was
    # switched (see DEVIATION_PANEL_VARIANT's own comment).
    legend_step = 0
    model_handles = [Line2D([], [], marker="o", linestyle="", color="black", markerfacecolor=model_colors[m],
                            markersize=MARKER_SIZE, label=m) for m in sorted(panel_models)]
    model_legend = axes[4].legend(handles=model_handles, title="Model (panels 1-4)", loc="upper left",
                                  bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - legend_step * LEGEND_Y_STEP),
                                  frameon=False, fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)
    axes[4].add_artist(model_legend)
    legend_step += 1

    # Panel 5 always draws base_models (the gravity-on half of every pair),
    # regardless of DEVIATION_PANEL_VARIANT. When that set differs from what
    # panels 1-4 draw (variant "no_acel"), those models would otherwise have
    # NO legend entry anywhere on the figure - give them their own, labelled
    # to say where they appear. Under "with_gravity" or "both", panel_models
    # already includes every base model, so this is empty and skipped.
    panel5_only_models = sorted(set(base_models) - set(panel_models))
    if panel5_only_models:
        panel5_handles = [Line2D([], [], marker="o", linestyle="", color="black", markerfacecolor=model_colors[m],
                                 markersize=MARKER_SIZE, label=m) for m in panel5_only_models]
        panel5_legend = axes[4].legend(handles=panel5_handles, title="Model (panel 5 only)", loc="upper left",
                                       bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - legend_step * LEGEND_Y_STEP),
                                       frameon=False, fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)
        axes[4].add_artist(panel5_legend)
        legend_step += 1

    quantity_handles = [
        Line2D([], [], marker="o", linestyle="", color="black", markerfacecolor="white",
               markersize=MARKER_SIZE, label="harmonic peak frequency"),
        Line2D([], [], marker="^", linestyle="", color="black", markerfacecolor="white",
               markersize=MARKER_SIZE, label="displacement stiffness"),
    ]
    axes[4].legend(handles=quantity_handles, title="Self-weight quantity", loc="upper left",
                  bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP - legend_step * LEGEND_Y_STEP), frameon=False,
                  fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    # Panel 1's title grows an extra line whenever DEVIATION_PANEL_VARIANT
    # adds its own title line (see deviation_panel_variant_title_suffix()) -
    # push this legend up accordingly so it never sits on top of that line.
    band_legend_y = 1.18 if not title_suffix else 1.34
    band_handle = [Patch(facecolor=BAND_COLOR, alpha=BAND_ALPHA, edgecolor="none",
                         label="measurement 25-75% spread")]
    axes[0].legend(handles=band_handle, loc="upper left", bbox_to_anchor=(0.0, band_legend_y), frameon=False,
                  fontsize=LEGEND_FONTSIZE)

    model_suffix = f" - models shown: {', '.join(base_models)}" if MODEL_INCLUDE is not None else ""
    variant_suffix = f" | GRAVITY_VARIANT={GRAVITY_VARIANT}" if GRAVITY_VARIANT != "with_gravity" else ""
    fig.suptitle(f"Method deviation vs. measurement ({len(tree_order)} trees){model_suffix}{variant_suffix}",
                fontsize=SUPTITLE_FONTSIZE)
    # PROVENANCE_FILTERS_TEXT already states GRAVITY_VARIANT (set in main()).
    add_provenance_footer(fig, "method_deviation_overview.py", PROVENANCE_NEWEST_INPUT_MTIME,
                          PROVENANCE_FILTERS_TEXT, n_trees=len(tree_order), n_models=len(panel_models))
    plt.tight_layout()
    fig.subplots_adjust(right=RIGHT_MARGIN, top=0.82)
    figure_a_png_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_A_PNG_PATH))
    savefig_with_provenance(fig, figure_a_png_path)
    plt.close(fig)
    print(f"\nFigure A saved to: {figure_a_png_path}")

    print("\nModel points falling INSIDE the measurement band (indistinguishable from measurement):")
    if inside_band_log:
        for tree, model, quantity in inside_band_log:
            print(f"  {tree} / {model} / {quantity}")
    else:
        print("  (none)")

    figure_a_csv_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_A_CSV_PATH))
    csv_df = pd.DataFrame(csv_rows)
    csv_df.to_csv(figure_a_csv_path, index=False)
    print(f"Figure A data saved to: {figure_a_csv_path}")
    return csv_df


# --- FIGURE B --------------------------------------------------------------

def build_figure_b(volume_total, volume_trunk, volume_branch, modal_effmass_df, panel_models,
                   tree_order, model_colors: dict) -> pd.DataFrame:
    # model_colors passed in - see build_figure_a()'s own comment for why.
    # The TOP row (volume ratio) is keyed by reconstruction METHOD
    # (method_columns below), not ANSYS model, so it has no gravity axis to
    # switch at all - base_models/panel_models are irrelevant to it. The
    # BOTTOM row (model mass / density) is keyed by ANSYS model and draws
    # panel_models (base_models resolved through GRAVITY_VARIANT - see
    # build_figure_d()'s comment) - mass does not change with gravity
    # either, so the plotted VALUE is identical to the gravity-on run;
    # only the model label switches.
    method_columns = ["adqsm", "adtree_raw", "adtree_calibrated", "treeqsm_optimal", "treeqsm_simplified"]
    method_colors = assign_label_colors(method_columns)
    offset_by_method = dict(zip(sorted(method_columns), compute_dodge_offsets(len(method_columns), DODGE_WIDTH)))
    offset_by_model = dict(zip(sorted(panel_models), compute_dodge_offsets(len(panel_models), DODGE_WIDTH)))

    csv_rows = []
    fig = plt.figure(figsize=(22, 11))
    subfig_top, subfig_bottom = fig.subfigures(2, 1, height_ratios=[1.0, 1.0])
    axes_top = subfig_top.subplots(1, 3, sharey=True)
    axes_bottom = subfig_bottom.subplots(1, 2)

    print("\nFIGURE B - method disagreement (volume ratio to per-tree median across "
          f"{method_columns}), gravity-off (_no_acel) runs excluded (n/a here - this "
          "panel is by reconstruction METHOD, not ANSYS model):")
    for ax, (quantity_label, volume_df) in zip(
            axes_top, [("total", volume_total), ("trunk", volume_trunk), ("branch", volume_branch)]):
        add_light_grid(ax)
        ax.axhline(1.0, color=ONE_LINE_COLOR, linewidth=ONE_LINE_WIDTH, zorder=2)
        ax.set_yscale("log")
        for tree_index, tree in enumerate(tree_order):
            row = volume_df[volume_df["tree"] == tree]
            if len(row) == 0:
                continue
            row = row.iloc[0]
            values = {m: row[m] for m in method_columns if m in row and pd.notna(row[m])}
            if not values:
                continue
            median_value = float(np.median(list(values.values())))
            for method, value in values.items():
                ratio = value / median_value
                x = tree_index + offset_by_method.get(method, 0.0)
                ax.scatter(x, ratio, color=method_colors[method], s=MARKER_SIZE ** 2, zorder=5, edgecolor="none")
                csv_rows.append({"section": "volume_ratio", "quantity": quantity_label, "tree": tree,
                                 "method": method, "value_m3": value, "median_across_methods_m3": median_value,
                                 "ratio_to_median": ratio, "model": None, "M_total_kg": None,
                                 "volume_m3": None, "density_kg_m3": None})
        ax.set_xticks(range(len(tree_order)))
        ax.set_xticklabels(tree_order, rotation=45, ha="right")
        ax.set_title(f"{quantity_label} volume - ratio to per-tree median\n(n={len(tree_order)} trees)",
                    fontsize=TITLE_FONTSIZE)
        ax.set_ylabel("Ratio to per-tree median (log scale)", fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)

    method_handles = [Line2D([], [], marker="o", linestyle="", color="black", markerfacecolor=color,
                             markersize=MARKER_SIZE, label=method) for method, color in method_colors.items()]
    axes_top[-1].legend(handles=method_handles, title="Reconstruction method", loc="upper left",
                        bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP), frameon=False,
                        fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    # --- bottom row: model mass, and mass/volume density check -------------
    title_suffix = deviation_panel_variant_title_suffix(GRAVITY_VARIANT)
    print(f"\nFIGURE B - model mass and M_total/volume density (nominal wood density = "
          f"{NOMINAL_WOOD_DENSITY_KG_M3} kg/m3, from scripts/ansys/bk_kp.mac 'densw'), "
          f"GRAVITY_VARIANT={GRAVITY_VARIANT!r}:")
    ax_mass, ax_density = axes_bottom
    add_light_grid(ax_mass)
    add_light_grid(ax_density)
    ax_density.axhline(NOMINAL_WOOD_DENSITY_KG_M3, color=DENSITY_LINE_COLOR, linestyle=DENSITY_LINE_STYLE,
                       linewidth=ONE_LINE_WIDTH, zorder=2)

    for tree_index, tree in enumerate(tree_order):
        for model in panel_models:
            mass_rows = modal_effmass_df[(modal_effmass_df["tree"] == tree) & (modal_effmass_df["model"] == model)]
            if len(mass_rows) == 0:
                continue
            m_total = float(mass_rows["M_total_kg"].iloc[0])
            x = tree_index + offset_by_model.get(model, 0.0)
            ax_mass.scatter(x, m_total, color=model_colors[model], s=MARKER_SIZE ** 2, zorder=5, edgecolor="none")

            # base_model_name() first - see build_figure_d()'s own comment:
            # a no_acel model's volume family is its BASE model's family,
            # since volume is geometry, not load.
            family = model_to_family(base_model_name(model))
            volume_column = MODEL_FAMILY_TO_VOLUME_COLUMN.get(family)
            volume_row = volume_total[volume_total["tree"] == tree]
            volume_m3 = float(volume_row[volume_column].iloc[0]) if (len(volume_row) > 0 and volume_column) else float("nan")
            density = m_total / volume_m3 if volume_m3 and not pd.isna(volume_m3) else float("nan")
            if pd.notna(density):
                ax_density.scatter(x, density, color=model_colors[model], s=MARKER_SIZE ** 2, zorder=5, edgecolor="none")
            print(f"  {tree} / {model}: M_total = {m_total:.1f} kg, total volume "
                  f"({volume_column}) = {volume_m3:.4g} m3 -> density = {density:.1f} kg/m3 "
                  f"(nominal = {NOMINAL_WOOD_DENSITY_KG_M3:.0f} kg/m3)")
            csv_rows.append({"section": "mass_density", "quantity": None, "tree": tree, "method": None,
                             "value_m3": None, "median_across_methods_m3": None, "ratio_to_median": None,
                             "model": model, "M_total_kg": m_total, "volume_m3": volume_m3,
                             "density_kg_m3": density})

    ax_mass.set_xticks(range(len(tree_order)))
    ax_mass.set_xticklabels(tree_order, rotation=45, ha="right")
    ax_mass.set_title(f"Model mass{title_suffix} (n={len(tree_order)} trees)", fontsize=TITLE_FONTSIZE)
    ax_mass.set_ylabel("M_total (kg)", fontsize=AXIS_LABEL_FONTSIZE)
    ax_mass.tick_params(labelsize=TICK_LABEL_FONTSIZE)

    ax_density.set_xticks(range(len(tree_order)))
    ax_density.set_xticklabels(tree_order, rotation=45, ha="right")
    ax_density.set_title(f"M_total / total volume vs. nominal wood density\n(n={len(tree_order)} trees)",
                         fontsize=TITLE_FONTSIZE)
    ax_density.set_ylabel("Density (kg/m3)", fontsize=AXIS_LABEL_FONTSIZE)
    ax_density.tick_params(labelsize=TICK_LABEL_FONTSIZE)

    model_handles_2 = [Line2D([], [], marker="o", linestyle="", color="black", markerfacecolor=model_colors[m],
                              markersize=MARKER_SIZE, label=m) for m in panel_models]
    density_line_handle = Line2D([0], [0], color=DENSITY_LINE_COLOR, linestyle=DENSITY_LINE_STYLE,
                                 label=f"nominal ({NOMINAL_WOOD_DENSITY_KG_M3:.0f} kg/m3)")
    ax_density.legend(handles=model_handles_2 + [density_line_handle], title="Model", loc="upper left",
                      bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP), frameon=False,
                      fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    model_suffix = f" - models shown: {', '.join(panel_models)}" if MODEL_INCLUDE is not None else ""
    subfig_top.suptitle(f"Method disagreement on volume - no destructive reference exists "
                        f"(n={len(tree_order)} trees){model_suffix}", fontsize=TITLE_FONTSIZE + 1)
    subfig_bottom.suptitle(f"Model mass and implied density{title_suffix} - no destructive "
                           f"reference exists (n={len(tree_order)} trees){model_suffix}",
                           fontsize=TITLE_FONTSIZE + 1)
    subfig_top.subplots_adjust(top=0.78, right=0.86, bottom=0.22, wspace=0.3)
    subfig_bottom.subplots_adjust(top=0.78, right=0.86, bottom=0.22, wspace=0.3)

    add_provenance_footer(fig, "method_deviation_overview.py", PROVENANCE_NEWEST_INPUT_MTIME,
                          PROVENANCE_FILTERS_TEXT, n_trees=len(tree_order), n_models=len(panel_models))

    figure_b_png_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_B_PNG_PATH))
    savefig_with_provenance(fig, figure_b_png_path)
    plt.close(fig)
    print(f"\nFigure B saved to: {figure_b_png_path}")

    figure_b_csv_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_B_CSV_PATH))
    csv_df = pd.DataFrame(csv_rows)
    csv_df.to_csv(figure_b_csv_path, index=False)
    print(f"Figure B data saved to: {figure_b_csv_path}")
    return csv_df


# --- FIGURE C ----------------------------------------------------------

def build_figure_c(volume_branch, frequency_df, panel_models, tree_order, model_colors: dict) -> pd.DataFrame:
    # model_colors passed in - see build_figure_a()'s own comment for why.
    # panel_models is base_models resolved through GRAVITY_VARIANT - see
    # build_figure_d()'s own comment (identical reasoning here). Unlike
    # FIGURE D, phi itself needs NO base_model_name() resolution - it is
    # computed once per TREE from volume_branch's reconstruction-METHOD
    # columns, before any per-model loop, so it is already gravity- (and
    # method-column-) independent by construction.

    print("\nFIGURE C - phi (branch volume ratio) definition check:")
    print("  Task text says phi = AdTree / TreeQSM. Computed that way, phi ranges "
          "~0.055-0.66 (does not match the task's stated ~1.5-17.7).")
    print("  The INVERSE, phi = TreeQSM_simplified / AdTree_calibrated, ranges "
          "~1.5-18.1 - matching the stated range almost exactly. Using the inverse "
          "definition below; flagging the discrepancy rather than silently picking one.")

    rows = []
    for _, vrow in volume_branch.iterrows():
        tree = vrow["tree"]
        if tree not in tree_order:
            continue
        adtree_v = vrow.get("adtree_calibrated")
        treeqsm_v = vrow.get("treeqsm_simplified")
        if pd.isna(adtree_v) or pd.isna(treeqsm_v) or adtree_v == 0:
            continue
        phi = treeqsm_v / adtree_v

        for model in panel_models:
            freq_rows = frequency_df[(frequency_df["tree"] == tree) & (frequency_df["model"] == model)
                                     & (frequency_df["quantity"] == "harmonic_peak")
                                     & (frequency_df["variant"] == FREQUENCY_VARIANT_FOR_ANALYSIS)]
            if len(freq_rows) == 0 or pd.isna(freq_rows["rel_dev_pct"].iloc[0]):
                continue
            rows.append({"tree": tree, "model": model, "phi": phi,
                        "branch_volume_adtree_calibrated_m3": adtree_v,
                        "branch_volume_treeqsm_simplified_m3": treeqsm_v,
                        "freq_rel_dev_pct": float(freq_rows["rel_dev_pct"].iloc[0])})

    data = pd.DataFrame(rows)
    n = len(data)
    print(f"\nFigure C data points: n={n}")
    print(data.to_string(index=False) if n else "  (no data)")
    if n > 0:
        data["log_phi"] = np.log10(data["phi"].to_numpy())

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    add_light_grid(ax)
    ax.set_xscale("log")

    # --- fit(s): see fit_groups_for_gravity_variant()'s own docstring ---
    fit_groups = fit_groups_for_gravity_variant(data, "log_phi", "freq_rel_dev_pct",
                                                GRAVITY_VARIANT) if n > 0 else []
    fit_linestyles = {"gravity ON": "-", "gravity OFF": "--"}
    fit_rows_for_csv = []
    for i, (label, subset, fit) in enumerate(fit_groups):
        if fit is None:
            print(f"Fewer than 2 points{' (' + label + ')' if label else ''} - no regression fit computed.")
            continue
        x_line = np.linspace(subset["log_phi"].min(), subset["log_phi"].max(), 50)
        y_line = fit["slope"] * x_line + fit["intercept"]
        ax.plot(10 ** x_line, y_line, color=FIT_LINE_COLOR, linewidth=FIT_LINE_WIDTH, zorder=4,
               linestyle=fit_linestyles.get(label, "-"))
        label_text = f" ({label})" if label else ""
        print(f"Fit{label_text} (y = slope*log10(phi) + intercept): slope={fit['slope']:.4f}, "
              f"intercept={fit['intercept']:.4f}, R^2={fit['r_squared']:.4f}, n={fit['n']}")
        annotation = f"slope={fit['slope']:.3g}\nintercept={fit['intercept']:.3g}\nR2={fit['r_squared']:.3f}\nn={fit['n']}"
        if label:
            annotation = f"{label}\n{annotation}"
        ax.annotate(annotation, xy=(0.02, 0.02 + i * 0.34), xycoords="axes fraction",
                   fontsize=FIT_ANNOTATION_FONTSIZE, va="bottom", ha="left",
                   bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.85))
        for _, row in subset.iterrows():
            fit_rows_for_csv.append({"tree": row["tree"], "model": row["model"], "fit_group": label,
                                    "slope": fit["slope"], "intercept": fit["intercept"],
                                    "r_squared": fit["r_squared"], "n": fit["n"]})

    for model in panel_models:
        model_data = data[data["model"] == model]
        if len(model_data) == 0:
            continue
        ax.scatter(model_data["phi"], model_data["freq_rel_dev_pct"], color=model_colors[model],
                  s=MARKER_SIZE ** 2, zorder=5, edgecolor="none", label=model)

    ax.axhline(0.0, color=ZERO_LINE_COLOR, linewidth=ZERO_LINE_WIDTH, zorder=2)
    ax.set_xlabel("phi = branch volume (TreeQSM_simplified / AdTree_calibrated), log scale",
                 fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Harmonic peak frequency, relative deviation (%)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    model_suffix = f" - models shown: {', '.join(panel_models)}" if MODEL_INCLUDE is not None else ""
    title_suffix = deviation_panel_variant_title_suffix(GRAVITY_VARIANT)
    ax.set_title(f"Frequency deviation vs. crown-mass ratio phi{title_suffix} "
                f"(n={n} points){model_suffix}", fontsize=TITLE_FONTSIZE)

    ax.legend(title="Model", loc="upper left", bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP), frameon=False,
             fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    add_provenance_footer(fig, "method_deviation_overview.py", PROVENANCE_NEWEST_INPUT_MTIME,
                          PROVENANCE_FILTERS_TEXT, n_trees=len(tree_order), n_models=len(panel_models))

    plt.tight_layout()
    fig.subplots_adjust(right=RIGHT_MARGIN)
    figure_c_png_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_C_PNG_PATH))
    savefig_with_provenance(fig, figure_c_png_path)
    plt.close(fig)
    print(f"Figure C saved to: {figure_c_png_path}")

    figure_c_csv_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_C_CSV_PATH))
    fit_info = pd.DataFrame(fit_rows_for_csv) if fit_rows_for_csv else pd.DataFrame(
        columns=["tree", "model", "fit_group", "slope", "intercept", "r_squared", "n"])
    data = data.merge(fit_info, on=["tree", "model"], how="left") if len(fit_info) > 0 else data.assign(
        fit_group=None, slope=float("nan"), intercept=float("nan"), r_squared=float("nan"), n=0)
    data.to_csv(figure_c_csv_path, index=False)
    print(f"Figure C data saved to: {figure_c_csv_path}")
    return data


# --- FIGURE D ----------------------------------------------------------

def build_figure_d(taper_df, stiffness_df, panel_models, tree_order, model_colors: dict) -> pd.DataFrame:
    # model_colors passed in - see build_figure_a()'s own comment for why.
    # panel_models is base_models resolved through GRAVITY_VARIANT (see
    # main()) - e.g. the two no_acel names under the default, base_models
    # itself under "with_gravity", or both under "both".

    rows = []
    if taper_df is not None:
        for model in panel_models:
            # base_model_name() FIRST: a no_acel model's taper family is
            # its BASE model's family - taper comes from the cylinder
            # geometry, which a no_acel run shares with its twin (gravity
            # is a LOAD, not a geometry change) - same one-line resolution
            # panel 3 of FIGURE A already uses. The looked-up VALUE is
            # therefore identical to the gravity-on run; only the model
            # LABEL (and, via panel_models, which rows of stiffness_df get
            # pulled in below) changes - that is expected and correct, not
            # a bug to paper over.
            family = model_to_family(base_model_name(model))
            column = f"{MODEL_FAMILY_TO_TAPER_PREFIX.get(family, family)}_{TAPER_ERROR_COLUMN_SUFFIX}"
            if family is None or column not in taper_df.columns:
                print(f"  WARNING: no taper family/column for model '{model}' - skipped in FIGURE D.")
                continue
            for _, trow in taper_df.iterrows():
                tree = trow["tree"]
                if tree not in tree_order:
                    continue
                disp_rows = stiffness_df[(stiffness_df["tree"] == tree) & (stiffness_df["model"] == model)
                                         & (stiffness_df["fem_column"] == "slope_displacement")
                                         & (stiffness_df["sensor"] == "combined")]
                if len(disp_rows) == 0 or pd.isna(disp_rows["rel_dev_pct"].iloc[0]):
                    continue
                rows.append({"tree": tree, "model": model, "taper_family": family,
                            "taper_error_pct": trow[column],
                            "displacement_rel_dev_pct": float(disp_rows["rel_dev_pct"].iloc[0])})

    data = pd.DataFrame(rows)
    n = len(data)
    print(f"\nFigure D data points: n={n}")
    print(data.to_string(index=False) if n else "  (no data)")

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    add_light_grid(ax)
    # Same symlog rationale as FIGURE A (see DEVIATION_PANEL_YSCALE's own
    # comment) - the known-bad treeqsm_pd002 displacement outliers (up to
    # ~+7987%/-2240%) would otherwise flatten every other point onto y=0.
    ax.set_yscale(DEVIATION_PANEL_YSCALE, linthresh=DEVIATION_PANEL_LINTHRESH)

    # --- fit(s): see fit_groups_for_gravity_variant()'s own docstring for
    # why "both" gets TWO separate fits rather than one pooled fit -----
    fit_groups = fit_groups_for_gravity_variant(data, "taper_error_pct", "displacement_rel_dev_pct",
                                                GRAVITY_VARIANT) if n > 0 else []
    fit_linestyles = {"gravity ON": "-", "gravity OFF": "--"}
    fit_rows_for_csv = []
    for i, (label, subset, fit) in enumerate(fit_groups):
        if fit is None:
            print(f"Fewer than 2 points{' (' + label + ')' if label else ''} - no regression fit computed.")
            continue
        x_line = np.linspace(subset["taper_error_pct"].min(), subset["taper_error_pct"].max(), 50)
        y_line = fit["slope"] * x_line + fit["intercept"]
        ax.plot(x_line, y_line, color=FIT_LINE_COLOR, linewidth=FIT_LINE_WIDTH, zorder=4,
               linestyle=fit_linestyles.get(label, "-"))
        label_text = f" ({label})" if label else ""
        print(f"Fit{label_text} (y = slope*x + intercept): slope={fit['slope']:.4f}, "
              f"intercept={fit['intercept']:.4f}, R^2={fit['r_squared']:.4f}, n={fit['n']}")
        annotation = f"slope={fit['slope']:.3g}\nintercept={fit['intercept']:.3g}\nR2={fit['r_squared']:.3f}\nn={fit['n']}"
        if label:
            annotation = f"{label}\n{annotation}"
        ax.annotate(annotation, xy=(0.02, 0.02 + i * 0.34), xycoords="axes fraction",
                   fontsize=FIT_ANNOTATION_FONTSIZE, va="bottom", ha="left",
                   bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.85))
        for _, row in subset.iterrows():
            fit_rows_for_csv.append({"tree": row["tree"], "model": row["model"], "fit_group": label,
                                    "slope": fit["slope"], "intercept": fit["intercept"],
                                    "r_squared": fit["r_squared"], "n": fit["n"]})

    for model in panel_models:
        model_data = data[data["model"] == model]
        if len(model_data) == 0:
            continue
        ax.scatter(model_data["taper_error_pct"], model_data["displacement_rel_dev_pct"],
                  color=model_colors[model], s=MARKER_SIZE ** 2, zorder=5, edgecolor="none", label=model)
        # Tree-number label next to each point (e.g. "07" for B21_S07) - two
        # models per tree land at different (x, y), so the label just needs
        # a small constant offset, not per-tree collision avoidance.
        for _, row in model_data.iterrows():
            tree_label = row["tree"].split("_S")[-1]
            ax.annotate(tree_label, (row["taper_error_pct"], row["displacement_rel_dev_pct"]),
                       xytext=(4, 4), textcoords="offset points", fontsize=TAPER_POINT_LABEL_FONTSIZE,
                       color=TAPER_POINT_LABEL_COLOR, zorder=6)

    ax.axhline(0.0, color=ZERO_LINE_COLOR, linewidth=ZERO_LINE_WIDTH, zorder=2)
    ax.set_xlabel(f"Trunk taper error (%, {TAPER_ERROR_COLUMN_SUFFIX})", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Displacement stiffness, relative deviation (%)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    model_suffix = f" - models shown: {', '.join(panel_models)}" if MODEL_INCLUDE is not None else ""
    title_suffix = deviation_panel_variant_title_suffix(GRAVITY_VARIANT)
    ax.set_title(f"Displacement stiffness deviation vs. taper error{title_suffix} "
                f"(n={n} points){model_suffix}", fontsize=TITLE_FONTSIZE)

    ax.legend(title="Model", loc="upper left", bbox_to_anchor=(LEGEND_BBOX_X, LEGEND_Y_TOP), frameon=False,
             fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_TITLE_FONTSIZE)

    variant_filters_text = f"{PROVENANCE_FILTERS_TEXT}"
    add_provenance_footer(fig, "method_deviation_overview.py", PROVENANCE_NEWEST_INPUT_MTIME,
                          variant_filters_text, n_trees=len(tree_order), n_models=len(panel_models))

    plt.tight_layout()
    fig.subplots_adjust(right=RIGHT_MARGIN)
    figure_d_png_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_D_PNG_PATH))
    savefig_with_provenance(fig, figure_d_png_path)
    plt.close(fig)
    print(f"Figure D saved to: {figure_d_png_path}")

    figure_d_csv_path = add_deviation_panel_variant_suffix(add_model_include_suffix(FIGURE_D_CSV_PATH))
    fit_info = pd.DataFrame(fit_rows_for_csv) if fit_rows_for_csv else pd.DataFrame(
        columns=["tree", "model", "fit_group", "slope", "intercept", "r_squared", "n"])
    data = data.merge(fit_info, on=["tree", "model"], how="left") if len(fit_info) > 0 else data.assign(
        fit_group=None, slope=float("nan"), intercept=float("nan"), r_squared=float("nan"), n=0)
    data.to_csv(figure_d_csv_path, index=False)
    print(f"Figure D data saved to: {figure_d_csv_path}")
    return data


def main():
    global TREES_SUFFIX, PROVENANCE_NEWEST_INPUT_MTIME, PROVENANCE_FILTERS_TEXT

    print_active_configuration("method_deviation_overview.py", CONFIG_CONSTANT_NAMES)

    # --- MODEL_INCLUDE / MODEL_INCLUDE_SUFFIX sanity check ---------------------
    # Writing filtered figures under the same filenames as the full ones
    # would silently and irreversibly overwrite them - refuse outright.
    if MODEL_INCLUDE is not None and not MODEL_INCLUDE_SUFFIX:
        print("MODEL_INCLUDE is set but MODEL_INCLUDE_SUFFIX is empty - refusing to run: "
              "writing filtered figures under the SAME filenames as the full ones would "
              "silently and irreversibly overwrite them. Set MODEL_INCLUDE_SUFFIX to a "
              "non-empty string and re-run.")
        return

    # --- input paths are the bare, unsuffixed names, ALWAYS -----------------
    # stiffness_deviation.csv/frequency_deviation.csv/modal_effmass.csv are
    # written by plot_stiffness_comparison.py/plot_modal_harmonic.py from
    # their OWN "fem_df"/"modal_df"/"harmonic_df" - the CSV-feeding tables
    # those scripts deliberately never gravity-filter (GRAVITY_VARIANT there
    # is a draw-only filter, exactly like MODEL_INCLUDE - see GRAVITY_VARIANT's
    # own comment in either script for the postmortem of the earlier version
    # of this stage that briefly broke this). So these three names are
    # config-independent by construction and never need resolving against
    # THIS script's own GRAVITY_VARIANT - resolve_input_path_for_gravity_variant()
    # is kept, unused, only in case that guarantee is ever deliberately
    # revisited.
    stiffness_deviation_csv_path = STIFFNESS_DEVIATION_CSV_PATH
    frequency_deviation_csv_path = FREQUENCY_DEVIATION_CSV_PATH
    modal_effmass_csv_path = MODAL_EFFMASS_CSV_PATH

    # --- freshness guard --------------------------------------------------------
    PROVENANCE_NEWEST_INPUT_MTIME = check_input_freshness(
        "method_deviation_overview.py", ANSYS_RAW_FOLDER, ANSYS_RAW_PATTERN_FOR_FRESHNESS,
        {"STIFFNESS_DEVIATION_CSV_PATH": stiffness_deviation_csv_path,
         "FREQUENCY_DEVIATION_CSV_PATH": frequency_deviation_csv_path,
         "MODAL_EFFMASS_CSV_PATH": modal_effmass_csv_path},
        ABORT_ON_STALE_INPUT,
    )

    print("=" * 78)
    print("INPUTS")
    print("=" * 78)
    have_stiffness = find_input(stiffness_deviation_csv_path, "stiffness_deviation.csv")
    have_frequency = find_input(frequency_deviation_csv_path, "frequency_deviation.csv")
    have_modal_effmass = find_input(modal_effmass_csv_path, "modal_effmass.csv")
    have_experimental = find_input(EXPERIMENTAL_CSV_PATH, "dynatree_data_merge_0926.csv")
    have_volume_total = find_input(VOLUME_TOTAL_CSV_PATH, "volume_summary_settled_total.csv")
    have_volume_trunk = find_input(VOLUME_TRUNK_CSV_PATH, "volume_summary_settled_trunk.csv")
    have_volume_branch = find_input(VOLUME_BRANCH_CSV_PATH, "volume_summary_settled_branch.csv")
    taper_path = find_taper_csv(TAPER_SEARCH_FOLDER)

    stiffness_df = pd.read_csv(stiffness_deviation_csv_path) if have_stiffness else None
    frequency_df = pd.read_csv(frequency_deviation_csv_path) if have_frequency else None
    modal_effmass_df = pd.read_csv(modal_effmass_csv_path) if have_modal_effmass else None
    volume_total = pd.read_csv(VOLUME_TOTAL_CSV_PATH) if have_volume_total else None
    volume_trunk = pd.read_csv(VOLUME_TRUNK_CSV_PATH) if have_volume_trunk else None
    volume_branch = pd.read_csv(VOLUME_BRANCH_CSV_PATH) if have_volume_branch else None
    taper_df = pd.read_csv(taper_path) if taper_path is not None else None

    print("\n" + "=" * 78)
    print("CONSISTENCY CHECK - 'model' column across all input CSVs")
    print("=" * 78)
    model_sets = {}
    if stiffness_df is not None:
        model_sets["stiffness_deviation.csv"] = set(stiffness_df["model"].unique())
    if frequency_df is not None:
        model_sets["frequency_deviation.csv"] = set(frequency_df["model"].unique())
    if modal_effmass_df is not None:
        model_sets["modal_effmass.csv"] = set(modal_effmass_df["model"].unique())
    base_models, all_models_found = check_model_naming_consistency(model_sets) if model_sets else ([], [])

    # Colours assigned ONCE, from EVERY model found in the loaded data,
    # BEFORE MODEL_INCLUDE (or STRICT_SELECTION's tree/scenario narrowing)
    # touches base_models below - assign_label_colors() picks a colour by
    # alphabetical INDEX, so recomputing it from an already-narrowed list
    # (the old behaviour: each build_figure_*() called
    # assign_label_colors(base_models) itself) shifted a model's colour
    # between a filtered and an unfiltered run. Passed into every
    # build_figure_*() call below instead of being recomputed locally.
    model_colors = assign_label_colors(all_models_found)

    # --- TREES: resolve the selection against what is actually in the loaded
    # tables, print the SELECTION REPORT, then narrow by it. ---
    available_trees = sorted(set().union(
        stiffness_df["tree"].unique() if stiffness_df is not None else [],
        frequency_df["tree"].unique() if frequency_df is not None else [],
    ))
    tree_sel = resolve_tree_selection(TREES, available_trees)

    print("\n" + "=" * 78)
    print("TREES selection")
    print("=" * 78)
    any_missing = print_selection_report("method_deviation_overview.py", tree_sel)
    if STRICT_SELECTION and any_missing:
        raise SystemExit(
            "method_deviation_overview.py: STRICT_SELECTION=True and the SELECTION REPORT above "
            "lists a requested tree with no data - see 'requested but MISSING'. "
            "Set STRICT_SELECTION=False to plot the partial selection anyway."
        )
    stiffness_df = apply_tree_selection(stiffness_df, tree_sel) if stiffness_df is not None else None
    frequency_df = apply_tree_selection(frequency_df, tree_sel) if frequency_df is not None else None
    modal_effmass_df = apply_tree_selection(modal_effmass_df, tree_sel) if modal_effmass_df is not None else None
    if volume_total is not None:
        volume_total = apply_tree_selection(volume_total, tree_sel)
    if volume_trunk is not None:
        volume_trunk = apply_tree_selection(volume_trunk, tree_sel)
    if volume_branch is not None:
        volume_branch = apply_tree_selection(volume_branch, tree_sel)
    TREES_SUFFIX = build_filter_suffix(tree_sel)
    PROVENANCE_FILTERS_TEXT = (f"GRAVITY_VARIANT={GRAVITY_VARIANT}  MODEL_INCLUDE={MODEL_INCLUDE}  "
                               f"TREES={TREES}")

    # --- MODEL_INCLUDE: resolve which models are eligible to be DRAWN, on
    # top of the base_models/self-weight structure above (see MODEL_INCLUDE's
    # own comment - it can only narrow, never widen, what each figure was
    # already going to show; it never touches stiffness_deviation.csv/
    # frequency_deviation.csv/modal_effmass.csv, which this script only reads). ---
    print("\n" + "=" * 78)
    print("MODEL_INCLUDE resolution")
    print("=" * 78)
    visible_models_list = resolve_model_filter(all_models_found, MODEL_INCLUDE)
    if MODEL_INCLUDE is not None and len(visible_models_list) == 0:
        print("MODEL_INCLUDE matched no model at all - nothing to draw. No figures written.")
        return
    base_models_before_filter = base_models
    if MODEL_INCLUDE is not None:
        base_models = [m for m in base_models if m in visible_models_list]
        print(f"Base models after MODEL_INCLUDE (were {base_models_before_filter}): {base_models}")
        if not base_models:
            print("MODEL_INCLUDE leaves no BASE (gravity-on) model at all - every figure here "
                  "compares a base model against a measurement or another method, so there is "
                  "nothing left to draw. No figures written.")
            return

    check_family_mapping_evidence()

    if not have_stiffness or not have_frequency:
        print("\nstiffness_deviation.csv and/or frequency_deviation.csv missing - "
              "cannot build any figure. Stopping.")
        return

    # --- THE GRAVITY-VARIANT STRUCTURAL LOCK - Figure A panels 1-4 only -----
    # stiffness_df/frequency_df, as loaded and tree-selected above, are
    # renamed here to make explicit that they are the ONLY variables in
    # this script allowed to contain BOTH gravity variants - see
    # GRAVITY_VARIANT's own comment. FIGURE A's self-weight panel (panel 5)
    # reads these UNFILTERED tables, by definition (it needs both halves of
    # every pair). stiffness_df_variant/frequency_df_variant, derived from
    # them via panel_models (resolve_deviation_panel_models(), reused from
    # Stage 5 so MODEL_INCLUDE's narrowing of base_models is respected too),
    # are what FIGURE A's panels 1-4 actually draw from.
    #
    # FINDING (Stage 7 Part B2): MORE than the self-weight panel reads the
    # unfiltered tables - build_figure_b/c/d also do, unchanged from before
    # this stage. That is judged CORRECT, not a gap to fix: those three
    # figures' entire model axis is base_models (gravity-ON names) by
    # design, and their underlying quantities (volume, mass, taper) do not
    # have a separate gravity-off measurement to switch to - taper and mass
    # come from the cylinder geometry, which gravity does not change (see
    # Stage 5's data-integrity check: taper error is IDENTICAL between a
    # model and its no_acel twin, to the last digit). Passing the
    # gravity-filtered tables into figures B/C/D would only make their
    # existing base_models filter match nothing under GRAVITY_VARIANT=
    # "no_acel" - reported here rather than silently "fixed" into breaking.
    stiffness_df_all_variants = stiffness_df
    frequency_df_all_variants = frequency_df
    panel_models = resolve_deviation_panel_models(base_models, all_models_found, GRAVITY_VARIANT)
    stiffness_df_variant = stiffness_df_all_variants[stiffness_df_all_variants["model"].isin(panel_models)]
    frequency_df_variant = frequency_df_all_variants[frequency_df_all_variants["model"].isin(panel_models)]
    print(f"\nGRAVITY_VARIANT={GRAVITY_VARIANT!r} - FIGURE A panels 1-4 draw: {sorted(panel_models)} "
          f"({len(stiffness_df_variant)} of {len(stiffness_df_all_variants)} stiffness row(s), "
          f"{len(frequency_df_variant)} of {len(frequency_df_all_variants)} frequency row(s))")

    # Tree order: trees with at least one base-model row in BOTH deviation
    # tables - the population FIGURE A's main 4 panels can actually cover.
    # Uses base_models (not panel_models) deliberately - see
    # resolve_deviation_panel_models()'s docstring: a base model and its
    # no_acel twin cover the SAME trees in this dataset, so the tree axis
    # is unaffected by which variant is drawn; base_models is also what
    # FIGURES B/C/D need this same tree_order for.
    stiffness_trees = set(stiffness_df_all_variants[stiffness_df_all_variants["model"].isin(base_models)]["tree"].unique())
    frequency_trees = set(frequency_df_all_variants[frequency_df_all_variants["model"].isin(base_models)]["tree"].unique())
    tree_order = sorted(stiffness_trees & frequency_trees)
    print(f"\nTrees covered by FIGURE A (base models, both stiffness and frequency data): "
          f"{len(tree_order)} -> {tree_order}")

    print("\n" + "=" * 78)
    print("FIGURE A - deviation vs. measurement")
    print("=" * 78)
    if have_experimental:
        experimental_df = load_and_filter_experimental(EXPERIMENTAL_CSV_PATH, MEASUREMENT, "M01, for stiffness/rotation")
        experimental_freq_df = load_and_filter_experimental(
            EXPERIMENTAL_CSV_PATH, FREQUENCY_MEASUREMENT_FILTER, "all measurements, for frequency")
        band_by_tree_and_quantity = {
            "displacement_stiffness": {t: compute_measurement_band(experimental_df, t, ["slope_Pt3", "slope_Pt4"]) for t in tree_order},
            "harmonic_peak_frequency": {t: compute_measurement_band(experimental_freq_df, t, ["first_frequency"]) for t in tree_order},
            "rotation_stiffness": {t: compute_measurement_band(experimental_df, t, ["slope_root_stiff_cam", "slope_root_stiff_nocam"]) for t in tree_order},
        }
        print_measurement_bands(band_by_tree_and_quantity)
        build_figure_a(stiffness_df_variant, frequency_df_variant, stiffness_df_all_variants,
                      frequency_df_all_variants, taper_df, experimental_df, experimental_freq_df,
                      base_models, tree_order, model_colors)
    else:
        print("dynatree_data_merge_0926.csv missing - FIGURE A skipped entirely "
              "(its measurement band cannot be computed without it).")

    print("\n" + "=" * 78)
    print("FIGURE B - method volume disagreement")
    print("=" * 78)
    if have_volume_total and have_volume_trunk and have_volume_branch and have_modal_effmass:
        volume_tree_order = sorted(set(tree_order) & set(volume_total["tree"].unique())
                                   & set(modal_effmass_df["tree"].unique()))
        print(f"Trees covered by FIGURE B: {len(volume_tree_order)} -> {volume_tree_order}")
        build_figure_b(volume_total, volume_trunk, volume_branch, modal_effmass_df,
                      panel_models, volume_tree_order, model_colors)
    else:
        print("One or more volume_summary_settled_*.csv / modal_effmass.csv missing - FIGURE B skipped.")

    print("\n" + "=" * 78)
    print("FIGURE C - frequency deviation vs. branch volume ratio (phi)")
    print("=" * 78)
    if have_volume_branch:
        build_figure_c(volume_branch, frequency_df, panel_models, tree_order, model_colors)
    else:
        print("volume_summary_settled_branch.csv missing - FIGURE C skipped.")

    print("\n" + "=" * 78)
    print("FIGURE D - taper error vs. displacement stiffness deviation")
    print("=" * 78)
    if taper_df is not None:
        build_figure_d(taper_df, stiffness_df, panel_models, tree_order, model_colors)
    else:
        print("Taper CSV missing/not found - FIGURE D skipped.")

    # --- orphaned filtered outputs (list only, never touched) ------------------
    report_orphaned_filtered_outputs(
        OUT_FOLDER,
        ("method_deviation_vs_measurement", "method_volume_disagreement",
         "frequency_vs_branch_volume", "taper_vs_stiffness"),
        "method_deviation_overview.py",
    )

    print("\nDone.")


if __name__ == "__main__":
    main()
