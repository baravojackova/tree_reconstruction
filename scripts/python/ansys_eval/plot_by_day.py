"""
plot_by_day.py
---------------
Companion script to filter_and_plot.py.

This is a general-purpose "day vs. parameter" plotting tool: it plots one
or more chosen numeric parameters against the measurement day, on already
-filtered data. COLUMNS_TO_PLOT is not limited to force columns - it can
be any numeric column(s) in the table, e.g. a single parameter like "DBH",
or several shown side by side, like ["F_horizontal_reg_max",
"F_horizontal_total_max"].

What this script does, in plain terms:
  1. Loads the merged measurement table (dynatree_data_merge.csv).
  2. Narrows it down ("filters" it) to only the rows you are interested in,
     using the exact same filter logic as filter_and_plot.py.
  3. Draws the individual measurements of your chosen parameter(s) as
     points, grouped by measurement day, one subplot per parameter.
     This lets you see how a value changes from day to day and how much
     the individual pulls/measurements on the same day scatter.
  4. Saves the resulting figure as a PNG image.

How to use it:
  - Edit the values in the "CONFIGURATION" section below (that is the only
    part you normally need to touch).
  - Run the script from a terminal in VS Code:  python plot_by_day.py
  - The plot will be saved next to this script, and also opened on screen
    if your environment supports it.

You do not need to understand pandas or matplotlib in depth to use this -
just change the CONFIGURATION values and re-run the script. The comments
below explain what each part of the code is doing, in case you want to
learn more or adjust the logic later.
"""

import math
import random
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folders holding the input data and this script's output files.
DATA_FOLDER = Path(r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\data")
OUT_FOLDER = Path(r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\out")

# Path to the CSV file.
CSV_PATH = DATA_FOLDER / "dynatree_data_merge.csv"

# Quick reference: measurement days and whether the trees had leaves on
# that day. Handy when choosing DAY / YEAR below.
# 2021-03-22 - no leaves (early spring)
# 2021-06-29 - leaves
# 2021-08-03 - leaves
# 2022-04-05 - no leaves (early spring)
# 2022-08-16 - leaves
# 2023-07-17 - leaves
# 2024-01-16 - no leaves (winter/frost)
# 2024-04-10 - no leaves (early spring)
# 2024-09-02 - leaves
# 2025-04-01 - no leaves (early spring)

# --- Filters -------------------------------------------------------------
# For each filter, set a value to keep only matching rows, or set it to
# None to not filter on that column at all (i.e. keep all values).
#
# You can filter by a single value (e.g. TREE = "BK01") or a list of
# several values (e.g. TREE = ["BK01", "BK04"]).

TREE = None          # tree ID, e.g. "BK01", or a list like ["BK01", "BK04"], or None
YEAR = None            # measurement year, e.g. 2024, or a list like [2021, 2022], or None
DAY = ["2021-03-22", "2022-04-05"]     # exact measurement date as text, e.g. "2024-09-02", or None
TYPE = None           # measurement condition, e.g. "normal", "mokro", "mraz", or None
MEASUREMENT = "M01"     # repetition ID, e.g. "M01", or None
PROBE = None           # sensor name, e.g. "Elasto(90)", "Pt3", or None
PULLNO = None          # pull number (0, 1, 2, 3), or None
SIZE_CATEGORY = None   # tree size class, e.g. "kategorie menší (<45 cm)",
                       # "kategorie větší (>45 cm)", or a list of those, or None

# REDUCTION filters on the derived "reduction_level" column (see
# add_reduction_level_column() below), not on the raw reduction_* columns.
# Careful: None and "none" mean different things here!
#   REDUCTION = None      -> no filtering at all, keep every row
#   REDUCTION = "none"    -> keep only rows where NO reduction was applied
#   REDUCTION = "20"      -> keep only rows with the 20% reduction
#   REDUCTION = "20+20"   -> keep only rows with the 20+20% reduction
#   REDUCTION = "40"      -> keep only rows with the 40% reduction
REDUCTION = None

# --- What to plot ---------------------------------------------------------
# Names of the numeric columns to plot. Each column gets its own subplot.
COLUMNS_TO_PLOT = [
    "slope_Pt3",
    "slope_Pt4",
]

PLOT_TYPE = "points"
# Only "points" is implemented for now. More types (e.g. "box" for
# boxplots) will be added later.

# --- Subplot grid shape -----------------------------------------------------
# How the subplots (one per column in COLUMNS_TO_PLOT) are arranged in a
# grid. Leave both as None to get a sensible default (2 columns, as many
# rows as needed). Set one or both to an integer to control the layout
# yourself - see compute_grid_shape() below for the exact rules.
N_ROWS = None   # set to an integer to force the number of subplot rows
N_COLS = None   # set to an integer to force the number of subplot columns

# --- Point labels and jitter -----------------------------------------------
# Whether to write a small text label ("M1-1" style) next to every point.
SHOW_POINT_LABELS = True
LABEL_FONTSIZE = 8

# Small horizontal random offset ("jitter") so overlapping points/labels on
# the same day are easier to read. Set to 0 to plot points exactly stacked
# on top of each other with no horizontal offset.
JITTER_WIDTH = 0.08

# Where to save the resulting plot image.
OUTPUT_PLOT_PATH = None
# If None, the output filename is generated automatically from
# COLUMNS_TO_PLOT (see build_output_path() below). Set this to a
# specific path/filename to override that.


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def get_reduction_level(row) -> str:
    """
    Work out the reduction level for a single row of the table.

    The raw data has three separate True/False columns, one per reduction
    level: "reduction_20", "reduction_20+20" and "reduction_40". For any
    given row, at most one of them is True at a time - they all describe
    the same underlying category, just spread across three columns.

    This function looks at the three columns and turns them into a single
    text value that is easier to filter on:
      - "20"    if reduction_20 is True
      - "20+20" if reduction_20+20 is True
      - "40"    if reduction_40 is True
      - "none"  if none of the three are True (no reduction applied)
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
    Add a new "reduction_level" column to the table, computed row by row
    from the three boolean reduction_* columns (see get_reduction_level
    above). This makes it possible to filter on a single, simple value
    ("none", "20", "20+20" or "40") instead of juggling three separate
    True/False columns.
    """
    df = df.copy()  # work on a copy so we never modify the original data
    # apply(..., axis=1) runs get_reduction_level() once per row
    df["reduction_level"] = df.apply(get_reduction_level, axis=1)
    return df


def apply_filters(df: pd.DataFrame) -> pd.DataFrame:
    """
    Take the full table (df) and return only the rows that match the
    filters defined in the CONFIGURATION section above.

    The pattern used for every filter is the same:
      - if the filter variable is None, skip it (keep all rows)
      - if it is a single value, keep rows where the column equals it
      - if it is a list, keep rows where the column value is in that list
    """
    filtered = df.copy()  # work on a copy so we never modify the original data

    # Add the "reduction_level" column before filtering, so REDUCTION can
    # filter on it just like any other column (see add_reduction_level_column
    # above for how it is computed).
    filtered = add_reduction_level_column(filtered)

    def filter_column(data, column_name, value):
        if value is None:
            return data  # no filtering requested for this column
        if isinstance(value, (list, tuple, set)):
            return data[data[column_name].isin(value)]
        return data[data[column_name] == value]

    filtered = filter_column(filtered, "tree", TREE)
    filtered = filter_column(filtered, "day", DAY)
    filtered = filter_column(filtered, "type", TYPE)
    filtered = filter_column(filtered, "measurement", MEASUREMENT)
    filtered = filter_column(filtered, "probe", PROBE)
    filtered = filter_column(filtered, "pullNo", PULLNO)
    filtered = filter_column(filtered, "size_category", SIZE_CATEGORY)
    filtered = filter_column(filtered, "reduction_level", REDUCTION)

    # YEAR is a bit special: the "day" column stores full dates as text
    # (e.g. "2024-09-02"), so to filter by year we first convert that
    # column to real dates, then compare the year part.
    if YEAR is not None:
        day_as_date = pd.to_datetime(filtered["day"])
        if isinstance(YEAR, (list, tuple, set)):
            filtered = filtered[day_as_date.dt.year.isin(YEAR)]
        else:
            filtered = filtered[day_as_date.dt.year == YEAR]

    return filtered


def compute_grid_shape(n_plots: int, n_rows, n_cols) -> tuple[int, int]:
    """
    Work out how many subplot rows and columns to use.

    n_plots is how many subplots we actually need (= len(COLUMNS_TO_PLOT)).
    n_rows / n_cols come straight from the N_ROWS / N_COLS configuration
    variables, so each of them can be None or a fixed integer.

    Rules (checked in this order):
      - both None      -> use 2 columns, and as many rows as needed to fit
                           all n_plots subplots (rows = ceil(n_plots / 2))
      - only n_rows set -> compute n_cols = ceil(n_plots / n_rows)
      - only n_cols set -> compute n_rows = ceil(n_plots / n_cols)
      - both set        -> use them exactly as given
    """
    if n_rows is None and n_cols is None:
        n_cols = 2
        n_rows = math.ceil(n_plots / n_cols)
    elif n_rows is None:  # n_cols is set, n_rows is not
        n_rows = math.ceil(n_plots / n_cols)
    elif n_cols is None:  # n_rows is set, n_cols is not
        n_cols = math.ceil(n_plots / n_rows)
    # else: both are already set by the user, use them as-is

    return n_rows, n_cols


def format_point_label(measurement, pull_no) -> str:
    """
    Build a short text label for one measurement point, e.g. "M1-1".

    - "measurement" is a repetition ID like "M01", "M02", ... The leading
      zero is dropped, so "M01" -> "M1".
    - "pullNo" is the pull number (0, 1, 2, 3, ...). It is appended after
      a dash, as an integer, e.g. "-1". If pullNo is missing (NaN), no
      dash/number is added at all.

    Examples:
      format_point_label("M01", 1)   -> "M1-1"
      format_point_label("M02", 0)   -> "M2-0"
      format_point_label("M05", nan) -> "M5"
    """
    # int(measurement[1:]) turns "M01" into 1, then we put "M" back in
    # front without the leading zero.
    measurement_part = f"M{int(measurement[1:])}"

    if pd.isna(pull_no):
        return measurement_part

    return f"{measurement_part}-{int(pull_no)}"


def build_output_path(columns: list) -> str:
    """
    Automatically build an output PNG filename from the list of column
    names being plotted, so you don't have to update OUTPUT_PLOT_PATH by
    hand every time you change COLUMNS_TO_PLOT.

    Steps:
      1. "Sanitize" each column name by replacing any character that is
         not a letter, digit, underscore, plus sign, dot or hyphen with an
         underscore. This keeps the filename safe to use on Windows too.
      2. If there are 1-3 columns, join their sanitized names with "_" and
         add ".png", e.g. ["DBH", "height"] -> "DBH_height.png".
      3. If there are more than 3 columns, using all of their names would
         make a very long filename, so instead we just use the first
         column's name plus "_and_<n-1>_more.png", e.g. 4 columns starting
         with "F_horizontal_reg_max" -> "F_horizontal_reg_max_and_3_more.png".
    """
    # \w already matches letters, digits and underscore, so the pattern
    # below only needs to additionally allow "+", "." and "-".
    safe_characters = re.compile(r"[^\w+.-]")

    def sanitize(name: str) -> str:
        return safe_characters.sub("_", name)

    sanitized_columns = [sanitize(column) for column in columns]

    if len(sanitized_columns) <= 3:
        return "_".join(sanitized_columns) + ".png"

    n_more = len(sanitized_columns) - 1
    return f"{sanitized_columns[0]}_and_{n_more}_more.png"


def plot_column_by_day(ax, filtered: pd.DataFrame, column_name: str, days: list) -> tuple[int, int]:
    """
    Draw one subplot for a single force column, with one x-axis position
    per measurement day (the shared `days` list, so all subplots line up).

    Returns (n_plotted, n_skipped): how many points were plotted, and how
    many rows were skipped because the value was missing (NaN).
    """
    n_plotted = 0
    n_skipped = 0

    for _, row in filtered.iterrows():
        value = row[column_name]

        if pd.isna(value):
            n_skipped += 1
            continue  # nothing to plot for this row

        # Find where this row's day sits on the shared x-axis.
        day_index = days.index(row["day"])

        # Add a small random horizontal offset ("jitter") so points/labels
        # for the same day do not sit exactly on top of each other.
        if JITTER_WIDTH > 0:
            x = day_index + random.uniform(-JITTER_WIDTH, JITTER_WIDTH)
        else:
            x = day_index

        ax.scatter(x, value, color="tab:blue")

        if SHOW_POINT_LABELS:
            label = format_point_label(row["measurement"], row["pullNo"])
            ax.annotate(
                label,
                (x, value),
                fontsize=LABEL_FONTSIZE,
                xytext=(3, 3),
                textcoords="offset points",
            )

        n_plotted += 1

    ax.set_xticks(range(len(days)))
    ax.set_xticklabels(days, rotation=45, ha="right")
    ax.set_title(column_name)
    ax.set_ylabel(column_name)

    return n_plotted, n_skipped


def main():
    # --- Step 0: check the requested plot type is supported -----------------
    if PLOT_TYPE != "points":
        raise NotImplementedError(
            f"PLOT_TYPE '{PLOT_TYPE}' is not implemented yet - only 'points' "
            f"is currently supported."
        )

    # --- Step 1: load the data --------------------------------------------
    print(f"Loading data from: {CSV_PATH}")
    df = pd.read_csv(CSV_PATH)
    print(f"Full table: {len(df)} rows, {len(df.columns)} columns")

    # --- Step 2: apply the filters ------------------------------------------
    filtered = apply_filters(df)
    print(f"After filtering: {len(filtered)} rows match your criteria")

    if len(filtered) == 0:
        print("No rows matched your filters - check the CONFIGURATION values "
              "(e.g. is the tree name spelled correctly, does that year exist?).")
        return

    # --- Step 3: work out the shared list of days (x-axis categories) -------
    # sorted() works here because the dates are stored as "YYYY-MM-DD" text,
    # which sorts the same alphabetically as chronologically.
    days = sorted(filtered["day"].unique())
    print(f"Days present after filtering ({len(days)}): {days}")

    # --- Step 4: work out the subplot grid shape -----------------------------
    n_plots = len(COLUMNS_TO_PLOT)
    n_rows, n_cols = compute_grid_shape(n_plots, N_ROWS, N_COLS)
    print(f"Subplot grid: {n_rows} rows x {n_cols} columns for {n_plots} column(s) to plot")

    # --- Step 5: draw one subplot per column ---------------------------------
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), squeeze=False)
    # axes is always a 2D array now (thanks to squeeze=False), so flatten it
    # to a simple list we can index 0, 1, 2, ... in the same order as
    # COLUMNS_TO_PLOT.
    flat_axes = axes.flatten()

    for i, column_name in enumerate(COLUMNS_TO_PLOT):
        ax = flat_axes[i]

        if column_name not in filtered.columns:
            print(f"Column '{column_name}' does not exist in the table - skipping this subplot.")
            ax.set_visible(False)
            continue

        n_plotted, n_skipped = plot_column_by_day(ax, filtered, column_name, days)
        print(f"Column '{column_name}': {n_plotted} points plotted, "
              f"{n_skipped} skipped (missing values)")

    # Hide any leftover subplot cells that are not used (e.g. a 2x3 grid
    # with only 5 columns to plot leaves 1 empty cell).
    for j in range(n_plots, len(flat_axes)):
        flat_axes[j].set_visible(False)

    # --- Step 6: save and show the figure -------------------------------------
    # Work out the output filename: use OUTPUT_PLOT_PATH as-is if it was
    # set, otherwise build one automatically from COLUMNS_TO_PLOT.
    if OUTPUT_PLOT_PATH is None:
        output_path = OUT_FOLDER / build_output_path(COLUMNS_TO_PLOT)
    else:
        output_path = OUTPUT_PLOT_PATH

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Plot saved to: {output_path}")
    plt.show()


if __name__ == "__main__":
    main()
