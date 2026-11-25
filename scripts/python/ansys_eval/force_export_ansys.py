"""
force_export_ansys.py
----------------------
Produces a compact, one-row-per-tree force export table, meant to be used
as load input for ANSYS simulations.

Background on the data (so the logic below makes sense):
  - The table has one row per (day, tree, type, measurement, pullNo, probe)
    combination. For a single physical pull, the force/moment columns
    (F_horizontal_reg_max, F_vertical_reg_max, ...) are duplicated across
    every probe row of that pull - they describe the pull as a whole, not
    the individual probe.
  - However, F_horizontal_reg_max is NOT perfectly identical across all
    probes of the same pull: probes "Elasto(90)", "blueMaj", "yellowMaj"
    share one value, while "Pt3" and "Pt4" share a slightly different
    value (differing by about 1%), because the regression window (30-90%
    of the pull's max load) is chosen per-probe based on regression
    quality (R2). To avoid preferring one probe over another, this script
    averages across all probe rows of a pull before taking medians across
    pulls (see average_across_probes_per_pull() below).
  - The "probeCam" column (not "probe") identifies sensor role - "Inclino"
    vs "Pt" - because which physical probe (blueMaj/yellowMaj) plays
    "InclinoCamera" vs "InclinoNoCamera" flips between measurement days.

What this script does, in plain terms:
  1. Loads the merged measurement table and checks it is the extended
     (v2) dataset that has the force/moment columns.
  2. Filters to the measurement day(s), measurement ID and type you
     configure below.
  3. Double-checks (via assert) that no "reduction" (crown reduction) had
     been applied yet in this data - this script is only meant for the
     baseline, unreduced measurements.
  4. For each tree, computes median force values (in N), median raw
     Inclino/Pt sensor readings (converted to SI units - radians/meters),
     and its size category.
  5. Writes one compact tab-delimited .txt file, one row per tree, ready
     to be read by ANSYS.

How to use it:
  - Edit the values in the "CONFIGURATION" section below.
  - Run the script from a terminal in VS Code:  python force_export_ansys.py
  - Check the printed console messages and the generated .txt file before
    using it as ANSYS input.
"""

import math
from pathlib import Path

import pandas as pd

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folders holding the input data and this script's output files.
DATA_FOLDER = Path(r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\data")
OUT_FOLDER = Path(r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\out")

CSV_PATH = DATA_FOLDER / "dynatree_data_merge_0926.csv"

# List of measurement days to pool together into one combined median per
# tree. Can contain one day or several - all matching rows across all of
# these days are treated as one pool of measurements per tree.
DAYS_TO_PROCESS = ["2021-03-22", "2022-04-05"]

MEASUREMENT = "M01"    # M01 is the static pulling test.
TYPE = "normal"        # measurement condition to include.

FORCE_COLUMN_HORIZONTAL = "F_horizontal_reg_max"
FORCE_COLUMN_VERTICAL = "F_vertical_reg_max"

OUTPUT_PATH = OUT_FOLDER / "for_ces.txt"
DELIMITER = "\t"
DECIMAL_PLACES = 4     # number of decimal places for numeric output


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def load_and_validate_data(csv_path: str) -> pd.DataFrame:
    """
    Load the CSV file and check that it is the extended "v2" dataset that
    contains the force/moment columns (90 columns), not the older "v1"
    dataset (72 columns, no force/moment data). This script cannot produce
    a meaningful force export from the v1 dataset, so we fail loudly and
    early instead of silently producing an empty/wrong output file.
    """
    df = pd.read_csv(csv_path)

    if "F_horizontal_reg_max" not in df.columns:
        raise ValueError(
            f"This script requires the extended v2 dataset with "
            f"force/moment columns - got {len(df.columns)} columns, expected "
            f"90. Are you using the old dynatree_data_merge.csv (v1, 72 "
            f"columns)?"
        )

    return df


def filter_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only the rows for the configured days, measurement ID and type.
    This is the "raw" pool of rows (still one row per probe, per pull)
    that all the per-tree statistics below are computed from.
    """
    filtered = df[
        df["day"].isin(DAYS_TO_PROCESS)
        & (df["measurement"] == MEASUREMENT)
        & (df["type"] == TYPE)
    ]
    print(f"After filtering to day(s)={DAYS_TO_PROCESS}, measurement="
          f"{MEASUREMENT!r}, type={TYPE!r}: {len(filtered)} rows match")
    return filtered


def assert_no_reduction_applied(filtered: pd.DataFrame) -> None:
    """
    This script assumes the configured days are all "baseline" measurements,
    i.e. no crown reduction had been applied to any tree yet. Rather than
    silently filtering out reduced rows, we hard-fail if we find any - a
    reduced row sneaking into DAYS_TO_PROCESS would silently bias the
    medians below, which is worse than crashing here.
    """
    reduction_columns = ["reduction_20", "reduction_20+20", "reduction_40"]

    for reduction_column in reduction_columns:
        violating_rows = filtered[filtered[reduction_column]]
        if len(violating_rows) > 0:
            first_violation = violating_rows.iloc[0]
            raise AssertionError(
                f"Reduction flag '{reduction_column}' is True for tree "
                f"{first_violation['tree']!r} on day {first_violation['day']!r} "
                f"- this script assumes no reduction had been applied yet; "
                f"check DAYS_TO_PROCESS."
            )


def average_across_probes_per_pull(df: pd.DataFrame, force_column: str) -> pd.Series:
    """
    Collapse the several probe rows of a single physical pull into one
    averaged value, so that no single probe/sensor is preferred - every
    pull contributes exactly one value, blended equally across whichever
    probes happened to record it.

    A "pull" is identified by the combination (day, measurement, pullNo).
    Within each such group, this simply takes the mean of force_column
    across all its probe rows.

    Example: a pull recorded by 3 probes with force_column values
    10.0, 10.0, 9.9 -> averaged to (10.0 + 10.0 + 9.9) / 3 = 9.97

    Returns a pandas Series with one averaged value per pull (one entry
    per distinct (day, measurement, pullNo) group found in df).
    """
    grouped = df.groupby(["day", "measurement", "pullNo"])[force_column]
    return grouped.mean()


def get_tree_size_gt45(tree_rows: pd.DataFrame, tree: str) -> int:
    """
    Look up whether a tree belongs to the "větší (>45 cm)" (bigger) or
    "menší (<45 cm)" (smaller) size category. This is a per-tree constant,
    so we just read it off the first row for that tree.

    Returns 1 for the bigger category, 0 for the smaller category.
    """
    size_category = tree_rows["size_category"].iloc[0]

    if "větší" in size_category:
        return 1
    elif "menší" in size_category:
        return 0
    else:
        raise ValueError(
            f"Tree {tree!r} has an unrecognized size_category value: "
            f"{size_category!r}"
        )


def dedupe_pt_rows(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return one row per physical pull for the Pt3/Pt4 optical displacement
    sensors, removing duplicate readings.

    Unlike the Inclino columns, "Pt3_reg_max" and "Pt4_reg_max" are
    properties of the pull as a whole, not of a specific probe row: for a
    given pull, BOTH columns are already filled in on BOTH the Pt3-probe
    row and the Pt4-probe row. That means simply pooling all non-null
    Pt3_reg_max/Pt4_reg_max values from every row would count each real
    reading twice (once from the Pt3 row, once from the identical values
    on the Pt4 row).

    This function fixes that by:
      1. keeping only rows where probe is "Pt3" or "Pt4"
      2. dropping duplicate rows for the same pull - identified by
         (day, measurement, pullNo) - keeping just the first one, since
         both rows of a pull carry identical Pt3_reg_max/Pt4_reg_max
         values anyway

    Example: if pull (day="2021-03-22", measurement="M01", pullNo=0) has
    a Pt3 row and a Pt4 row, both with Pt3_reg_max=28.81 and
    Pt4_reg_max=27.47, this function keeps only one of those two rows.
    """
    pt_rows = df[df["probe"].isin(["Pt3", "Pt4"])]
    return pt_rows.drop_duplicates(subset=["day", "measurement", "pullNo"], keep="first")


def compute_tree_row(tree: str, tree_rows: pd.DataFrame) -> dict:
    """
    Compute the single output row (a dict of column name -> value) for one
    tree, from that tree's rows within the already-filtered data.
    """
    size_gt45 = get_tree_size_gt45(tree_rows, tree)

    # --- Overall horizontal/vertical force medians (in N): average across
    # probes per pull first (so no probe is preferred), then take the
    # median across the resulting per-pull values, then convert kN -> N. ---
    horiz_per_pull = average_across_probes_per_pull(tree_rows, FORCE_COLUMN_HORIZONTAL).dropna()
    vert_per_pull = average_across_probes_per_pull(tree_rows, FORCE_COLUMN_VERTICAL).dropna()

    f_horiz_med = horiz_per_pull.median() * 1000 if len(horiz_per_pull) > 0 else None
    f_vert_med = vert_per_pull.median() * 1000 if len(vert_per_pull) > 0 else None

    print(f"Tree '{tree}': {len(horiz_per_pull)} pulls contributed to F_horiz_med_N, "
          f"{len(vert_per_pull)} pulls contributed to F_vert_med_N")

    if f_horiz_med is None:
        print(f"  WARNING: tree '{tree}' has no valid data for F_horiz_med_N - leaving it blank.")
    if f_vert_med is None:
        print(f"  WARNING: tree '{tree}' has no valid data for F_vert_med_N - leaving it blank.")

    # --- Inclino median (in radians): blueMaj_reg_max and yellowMaj_reg_max
    # are mutually exclusive per row (at most one is non-null), so we can
    # simply pool both columns' non-null values with no deduplication. ---
    inclino_values = pd.concat([
        tree_rows["blueMaj_reg_max"].dropna(),
        tree_rows["yellowMaj_reg_max"].dropna(),
    ])
    inclino_median_rad = math.radians(inclino_values.median()) if len(inclino_values) > 0 else None

    print(f"Tree '{tree}': {len(inclino_values)} readings contributed to Inclino_reg_max_median_rad")

    if inclino_median_rad is None:
        print(f"  WARNING: tree '{tree}' has no valid data for Inclino_reg_max_median_rad - leaving it blank.")

    # --- Pt median (in meters): Pt3_reg_max/Pt4_reg_max are duplicated
    # across the Pt3 and Pt4 probe rows of the same pull, so we must
    # deduplicate by pull first (see dedupe_pt_rows above) before pooling. ---
    pt_pull_rows = dedupe_pt_rows(tree_rows)
    pt_values = pd.concat([
        pt_pull_rows["Pt3_reg_max"].dropna(),
        pt_pull_rows["Pt4_reg_max"].dropna(),
    ])
    pt_median_m = (pt_values.median() / 1000) if len(pt_values) > 0 else None

    print(f"Tree '{tree}': {len(pt_values)} readings contributed to Pt_reg_max_median_m")

    if pt_median_m is None:
        print(f"  WARNING: tree '{tree}' has no valid data for Pt_reg_max_median_m - leaving it blank.")

    return {
        "size_gt45": size_gt45,
        "tree": tree,
        "F_horiz_med_N": f_horiz_med,
        "F_vert_med_N": f_vert_med,
        "Inclino_reg_max_median_rad": inclino_median_rad,
        "Pt_reg_max_median_m": pt_median_m,
    }


def format_value(value) -> str:
    """
    Format one output value as a string: a fixed number of decimal places
    for numbers, or an empty string for missing (None) values.
    """
    if value is None:
        return ""
    return f"{value:.{DECIMAL_PLACES}f}"


def tree_id_to_number(tree_id: str) -> int:
    """
    Strip any non-digit characters from a tree ID string and return the
    remaining digits as an int.

    Example: tree_id_to_number("BK01") -> 1
             tree_id_to_number("BK16") -> 16
    """
    digits_only = "".join(character for character in tree_id if character.isdigit())
    return int(digits_only)


def write_output(rows: list, output_path: str) -> None:
    """
    Write the final export table as a DELIMITER-separated text file, using
    plain file writing (not pandas.to_csv) so the exact two-header-row
    structure is explicit and easy to follow:
      - row 1: short column names
      - row 2: column indices 0..6
      - then one row per tree: row index (starting at 1), then the 6
        computed values
    """
    column_names = [
        "index", "size_gt45", "tree", "F_horiz_med_N", "F_vert_med_N",
        "Inclino_reg_max_median_rad", "Pt_reg_max_median_m",
    ]
    column_indices = [str(i) for i in range(len(column_names))]

    with open(output_path, "w", encoding="utf-8") as output_file:
        output_file.write(DELIMITER.join(column_names) + "\n")
        output_file.write(DELIMITER.join(column_indices) + "\n")

        for row_index, row in enumerate(rows, start=1):
            line_values = [
                str(row_index),
                str(row["size_gt45"]),
                str(tree_id_to_number(row["tree"])),
                format_value(row["F_horiz_med_N"]),
                format_value(row["F_vert_med_N"]),
                format_value(row["Inclino_reg_max_median_rad"]),
                format_value(row["Pt_reg_max_median_m"]),
            ]
            output_file.write(DELIMITER.join(line_values) + "\n")


def main():
    # --- Step 1: load and validate the data ---------------------------------
    print(f"Loading data from: {CSV_PATH}")
    df = load_and_validate_data(CSV_PATH)
    print(f"Full table: {len(df)} rows, {len(df.columns)} columns")

    # --- Step 2: filter to the configured days/measurement/type -------------
    filtered = filter_rows(df)

    if len(filtered) == 0:
        print("No rows matched the configuration - check DAYS_TO_PROCESS, "
              "MEASUREMENT and TYPE.")
        return

    # --- Step 3: make sure no reduction had been applied yet ----------------
    assert_no_reduction_applied(filtered)

    # --- Step 4: determine the sorted list of trees (defines row order) -----
    trees = sorted(filtered["tree"].unique())
    print(f"Trees present after filtering ({len(trees)}): {trees}")

    # --- Steps 5-6: compute the per-tree statistics --------------------------
    rows = []
    for tree in trees:
        tree_rows = filtered[filtered["tree"] == tree]
        rows.append(compute_tree_row(tree, tree_rows))

    # --- Step 7: write the output file ---------------------------------------
    write_output(rows, OUTPUT_PATH)
    print(f"Export written to: {OUTPUT_PATH}")

    # --- Step 8: reminder about what this export does NOT cover -------------
    print("NOTE: columns Inclino_reg_max_median_rad / Pt_reg_max_median_m are "
          "raw sensor readings (inclination in radians, displacement in "
          "meters), not force values - see the module docstring for how "
          "each is computed.")


if __name__ == "__main__":
    main()
