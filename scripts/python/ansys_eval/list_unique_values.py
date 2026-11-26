"""
list_unique_values.py
-----------------------
Small helper script - separate from filter_and_plot.py on purpose, so you
can quickly check "what values exist in this column?" (e.g. which probe
names are available) without loading the full filtering/plotting logic.

How to use:
  - Edit COLUMN_NAME below to the column you want to inspect.
  - Run: python list_unique_values.py
"""

import pandas as pd
from pathlib import Path
import sys

# This script lives in scripts/python/ansys_eval, but paths.py lives one
# folder up (scripts/python). Python only looks for imports next to the
# script, so we first add the folder above to its search list.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

# =========================================================================
# CONFIGURATION
# =========================================================================

# Folder holding the input data.
DATA_FOLDER = paths.ANSYS_EVAL_DIR / "data"

CSV_PATH = DATA_FOLDER / "dynatree_data_merge.csv"

# Name of the column you want to inspect, e.g. "probe", "tree", "type",
# "measurement", "day", "size_category" ...
COLUMN_NAME = "day"

# If a column has more unique values than this, listing every single value
# is not useful (e.g. a continuous measurement column) - the script will
# print basic statistics instead. You can raise/lower this threshold.
MAX_VALUES_TO_LIST = 50


# =========================================================================
# CODE
# =========================================================================

def main():
    df = pd.read_csv(CSV_PATH)

    if COLUMN_NAME not in df.columns:
        print(f"Column '{COLUMN_NAME}' was not found in the table.")
        print("Available columns are:")
        for col in df.columns:
            print(f"  {col}")
        return

    column = df[COLUMN_NAME]
    n_unique = column.nunique(dropna=False)  # dropna=False also counts missing values as a category
    n_missing = column.isna().sum()

    print(f"Column: {COLUMN_NAME}")
    print(f"Number of unique values (including missing, if any): {n_unique}")
    print(f"Missing (NaN) values: {n_missing}")
    print()

    if n_unique > MAX_VALUES_TO_LIST:
        # Too many distinct values to list one by one - this is typically a
        # continuous numeric column, so summary statistics are more useful.
        print(f"This column has more than {MAX_VALUES_TO_LIST} unique values, "
              f"so it looks like a continuous/numeric column rather than a "
              f"category. Showing summary statistics instead:")
        print(column.describe())
    else:
        # Few distinct values - list each one with how many rows have it.
        # value_counts(dropna=False) also shows how many rows are NaN.
        counts = column.value_counts(dropna=False)
        print("Value -> number of rows:")
        for value, count in counts.items():
            print(f"  {value!r:30} -> {count}")


if __name__ == "__main__":
    main()
