"""
plot_by_tree.py
----------------
Companion script to filter_and_plot.py and plot_by_day.py.

This is a general-purpose "tree vs. parameter" plotting tool: it shows one
or more chosen numeric parameters grouped by tree on the x-axis, with the
different columns placed side by side (not overlapping) within each
tree's group, and every point colored by which measurement day it comes
from. This makes it easy to compare parameters across trees, and to see
whether a given day tends to sit higher or lower than others.

What this script does, in plain terms:
  1. Loads the merged measurement table (dynatree_data_merge.csv).
  2. Narrows it down ("filters" it) to only the rows you are interested in,
     using the exact same filter logic as filter_and_plot.py / plot_by_day.py.
  3. Draws the individual measurements of your chosen parameter(s) as
     points, grouped by tree, with one small "column" of points per
     parameter within each tree's group, colored by measurement day.
  4. Saves the resulting figure as a PNG image.

How to use it:
  - Edit the values in the "CONFIGURATION" section below (that is the only
    part you normally need to touch).
  - Run the script from a terminal in VS Code:  python plot_by_tree.py
  - The plot will be saved next to this script, and also opened on screen
    if your environment supports it.

You do not need to understand pandas or matplotlib in depth to use this -
just change the CONFIGURATION values and re-run the script. The comments
below explain what each part of the code is doing, in case you want to
learn more or adjust the logic later.
"""

import random
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

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

TREE = None            # tree ID, e.g. "BK01", or a list like ["BK01", "BK04"], or None
YEAR = None            # measurement year, e.g. 2024, or a list like [2021, 2022], or None
DAY = ["2021-03-22", "2022-04-05"]     # exact measurement date as text, e.g. "2024-09-02", or None
TYPE = 'normal'           # measurement condition, e.g. "normal", "mokro", "mraz", or None
MEASUREMENT = ["M02", "M03", "M04", "M05", "M06"]     # repetition ID, e.g. "M01", or None
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
COLUMNS_TO_PLOT = ["first_frequency"]
# List of numeric columns to compare side by side within each tree's
# group. Works with 1 or more columns.
# ["slope_root_stiff_cam", "slope_root_stiff_nocam"]

GROUP_WIDTH = 0.6
# How much horizontal space (out of 1.0 available per tree) the
# side-by-side columns occupy within one tree's group. Leaves a gap
# between neighboring trees' groups.

JITTER_WIDTH = 0.04
# Small random horizontal offset within each column's own sub-position,
# so points with similar values don't sit exactly on top of each other.
# Set to 0 to disable.

OUTPUT_PLOT_PATH = None
# If None, auto-generate the filename from COLUMNS_TO_PLOT (see
# build_output_path() below). Set to a specific path to override.


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


def compute_dodge_offsets(n_columns: int, group_width: float) -> list:
    """
    Work out the horizontal "dodge" offsets used to place n_columns
    side-by-side sub-positions within one tree's group, evenly spaced and
    centered on 0, spanning a total width of group_width.

    Example:
      compute_dodge_offsets(2, 0.6) -> [-0.15, 0.15]
      compute_dodge_offsets(1, 0.6) -> [0.0]
      compute_dodge_offsets(3, 0.6) -> [-0.2, 0.0, 0.2]
    """
    if n_columns == 1:
        return [0.0]

    # Evenly spaced positions from -group_width/2 to +group_width/2.
    step = group_width / (n_columns - 1)
    return [-group_width / 2 + i * step for i in range(n_columns)]


def assign_day_colors(days: list) -> dict:
    """
    Assign a distinct color to each measurement day, using matplotlib's
    "tab10" colormap (10 distinct colors). If there are more than 10 days,
    colors start repeating (cycling through the colormap with modulo).

    Returns a dict like {"2021-03-22": (r, g, b, a), "2022-04-05": (...), ...}
    """
    colormap = mpl.colormaps["tab10"]
    return {day: colormap(i % 10) for i, day in enumerate(days)}


def plot_tree_groups(ax, filtered: pd.DataFrame, columns: list, trees: list, day_colors: dict):
    """
    Draw the whole tree-grouped scatter plot onto `ax`.

    For each tree (x position i) and each column (dodge offset from
    compute_dodge_offsets), this plots every non-missing value as a point,
    colored by measurement day. It then labels each (tree, column)
    subgroup with the column name, placed just above its highest point.

    Returns nothing - it only draws onto `ax` and prints a per-subgroup
    "plotted vs skipped" console message.
    """
    offsets = compute_dodge_offsets(len(columns), GROUP_WIDTH)

    # --- Pass 1: plot all the points, and remember where each subgroup's
    # label should go (its x-center and its highest plotted y-value). ---
    subgroup_label_info = []  # list of (x_center, max_y, column_name) tuples

    for tree_index, tree in enumerate(trees):
        tree_rows = filtered[filtered["tree"] == tree]

        for column_name, offset in zip(columns, offsets):
            x_center = tree_index + offset
            n_plotted = 0
            n_skipped = 0
            max_y = None

            for _, row in tree_rows.iterrows():
                value = row[column_name]

                if pd.isna(value):
                    n_skipped += 1
                    continue  # nothing to plot for this row

                # Add a small random horizontal offset ("jitter") so points
                # with similar values do not sit exactly on top of each other.
                if JITTER_WIDTH > 0:
                    x = x_center + random.uniform(-JITTER_WIDTH, JITTER_WIDTH)
                else:
                    x = x_center

                color = day_colors[row["day"]]
                ax.scatter(x, value, color=color)

                if max_y is None or value > max_y:
                    max_y = value

                n_plotted += 1

            print(f"Tree '{tree}', column '{column_name}': {n_plotted} points plotted, "
                  f"{n_skipped} skipped (missing values)")

            if max_y is not None:
                subgroup_label_info.append((x_center, max_y, column_name))
            else:
                print(f"  -> no data for tree '{tree}', column '{column_name}' - skipping its label.")

    # --- Pass 2: now that every point is plotted, we know the y-axis
    # range, so we can offset the labels upward by a fixed fraction of it. ---
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    label_offset = 0.03 * y_range  # 3% of the y-axis range

    for x_center, max_y, column_name in subgroup_label_info:
        ax.text(x_center, max_y + label_offset, column_name,
                ha="center", va="bottom", fontsize=8)

    # --- x-axis: one tick per tree, labeled with the tree ID -----------------
    ax.set_xticks(range(len(trees)))
    ax.set_xticklabels(trees)

    # --- y-axis label: the column name if there is only one, else generic ---
    if len(columns) == 1:
        ax.set_ylabel(columns[0])
    else:
        ax.set_ylabel("value")

    ax.set_title("Distribution of " + ", ".join(columns) + " by tree")

    # --- Legend: one entry per day, in chronological order -------------------
    # We plot invisible ("proxy") scatter points just to generate legend
    # handles with the right color and label - they are not part of the data.
    sorted_days = sorted(day_colors.keys())
    legend_handles = [
        ax.scatter([], [], color=day_colors[day], label=day)
        for day in sorted_days
    ]
    ax.legend(handles=legend_handles, title="Measurement day",
              bbox_to_anchor=(1.02, 1), loc="upper left")


def main():
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

    # --- Step 3: work out the trees and days present after filtering --------
    trees = sorted(filtered["tree"].unique())
    days = sorted(filtered["day"].unique())
    print(f"Trees present after filtering ({len(trees)}): {trees}")
    print(f"Days present after filtering ({len(days)}): {days}")

    if len(trees) == 0:
        print("No trees found after filtering - nothing to plot.")
        return

    day_colors = assign_day_colors(days)

    # --- Step 4: draw the plot ------------------------------------------------
    fig, ax = plt.subplots(1, 1, figsize=(max(8, 1.5 * len(trees)), 6))
    plot_tree_groups(ax, filtered, COLUMNS_TO_PLOT, trees, day_colors)

    # --- Step 5: save and show the figure -------------------------------------
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
