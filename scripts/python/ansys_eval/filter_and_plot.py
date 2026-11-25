"""
filter_and_plot.py
-------------------
Step 1 of the DYNATREE data exploration workflow.

What this script does, in plain terms:
  1. Loads the merged measurement table (dynatree_data_merge.csv).
  2. Narrows it down ("filters" it) to only the rows you are interested in,
     e.g. one tree, one year, one measurement type.
  3. Prints how many rows matched, so you always know if a filter was too
     strict (0 rows) or did nothing (all rows).
  4. Draws a histogram (distribution plot) of one numeric column from the
     filtered data, and saves it as a PNG image.

How to use it:
  - Edit the values in the "CONFIGURATION" section below (that is the only
    part you normally need to touch).
  - Run the script from a terminal in VS Code:  python filter_and_plot.py
  - The plot will be saved next to this script, and also opened on screen
    if your environment supports it.

You do not need to understand pandas or matplotlib in depth to use this -
just change the CONFIGURATION values and re-run the script. The comments
below explain what each part of the code is doing, in case you want to
learn more or adjust the logic later.
"""

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

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

TREE = "BK01"          # tree ID, e.g. "BK01", or a list like ["BK01", "BK04"], or None
YEAR = [2021,2022]           # measurement year, e.g. 2024, or None for all years
DAY = ["2021-03-22", "2022-04-05"]        # exact measurement date as text, e.g. "2024-09-02", or None
TYPE = None            # measurement condition, e.g. "normal", "mokro", "mraz", or None
MEASUREMENT = None     # repetition ID, e.g. "M01", or None
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
# Name of the numeric column whose distribution (histogram) you want to see.
# Some example columns you might try: "Slope", "first_frequency", "DBH",
# "height", "fresh_mass".
COLUMN_TO_PLOT = "Slope"

# Number of bars ("bins") in the histogram. More bins = finer detail.
BINS = 30

# Where to save the resulting plot image.
OUTPUT_PLOT_PATH = OUT_FOLDER / "distribution_plot.png"


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

    # --- Step 3: check the column to plot exists and has data ---------------
    if COLUMN_TO_PLOT not in filtered.columns:
        print(f"Column '{COLUMN_TO_PLOT}' does not exist in the table. "
              f"Check the spelling against the column list in the CSV header.")
        return

    values = filtered[COLUMN_TO_PLOT].dropna()  # drop missing values before plotting
    n_missing = len(filtered) - len(values)
    print(f"Column '{COLUMN_TO_PLOT}': {len(values)} valid values, "
          f"{n_missing} missing (NaN) values excluded from the plot")

    if len(values) == 0:
        print(f"All values of '{COLUMN_TO_PLOT}' are missing for this filter - "
              f"nothing to plot.")
        return

    # --- Step 4: draw and save the histogram --------------------------------
    plt.figure(figsize=(8, 5))
    plt.hist(values, bins=BINS, edgecolor="black")
    plt.xlabel(COLUMN_TO_PLOT)
    plt.ylabel("Count")
    plt.title(f"Distribution of {COLUMN_TO_PLOT}\n(n = {len(values)} measurements)")
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT_PATH, dpi=150)
    print(f"Plot saved to: {OUTPUT_PLOT_PATH}")
    plt.show()


if __name__ == "__main__":
    main()
