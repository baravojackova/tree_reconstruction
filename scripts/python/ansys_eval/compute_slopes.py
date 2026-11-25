"""
compute_slopes.py
------------------
Stage 2 of the ANSYS post-processing workflow.

Combines the Stage 1 output table (static_scenarios_combined.csv, produced
by load_static_scenarios.py) with tree geometry data
(Popis_Babice_VSE_09012025_ans_lev6_0.txt) to compute lever arms, bending
moments, and stiffness slopes per tree/scenario, and saves the result as a
new CSV file.

What this script does, in plain terms:
  1. Loads the tree geometry file and, for each tree, works out two lever
     arms (vertical distances between anchorage height and two reference
     levels).
  2. Checks that the hardcoded column positions really point at the
     columns we think they do, and that the three heights are ordered
     sensibly for every tree - see verify_geometry_columns() below.
  3. Loads the Stage 1 table and merges in those lever arms by tree.
  4. Computes bending moments (force x lever arm) and stiffness slopes
     (moment / rotation or displacement).
  5. Saves the combined table as a CSV file for review.

WHY THE VERIFICATION STEP EXISTS: the columns in this file are read BY
POSITION, not by name, because the header names are not consistently
formatted. If the file is ever replaced by a different version with
different columns, reading by position does not fail - it just silently
returns the wrong numbers. The check in step 2 turns that silent failure
into a loud one.

How to use it:
  - Edit the values in the "CONFIGURATION" section below if needed.
  - Run the script from a terminal in VS Code:  python compute_slopes.py
  - Read the VERIFICATION block that gets printed first, then check the
    printed table and the generated CSV file.
"""

from pathlib import Path

import pandas as pd

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folder this script's own output, and Stage 1's output, live in.
OUT_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\out"
STAGE1_INPUT_PATH = Path(OUT_FOLDER) / "static_scenarios_combined.csv"   # output of load_static_scenarios.py

# The description file is read from the SAME folder that ANSYS reads, so
# the two copies can never drift apart. These two constants must match
# the ones in load_static_scenarios.py.
INPUT_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\ansys"
GEOMETRY_FILE_NAME = "Popis_Babice_VSE_09012025_ans_lev6_0_azimcor.txt"
GEOMETRY_FILE_PATH = Path(INPUT_FOLDER) / GEOMETRY_FILE_NAME

OUTPUT_PATH = Path(OUT_FOLDER) / "slopes_combined.csv"
DECIMAL_PLACES = 6

# Canonical tree ID prefix - must match load_static_scenarios.py.
TREE_ID_PREFIX = "B21_S"

# Positions (0-based) of the columns we need in the geometry file. We read
# by position rather than by name because the header names in that file
# are not all consistently formatted.
GEOMETRY_COLUMN_TREE_NUMBER = 2
GEOMETRY_COLUMN_HEIGHT_ANCHORAGE = 7
GEOMETRY_COLUMN_HEIGHT_LEVEL3 = 9
GEOMETRY_COLUMN_HEIGHT_LEVEL6 = 12

# Plausibility bounds used by the verification step. A height outside
# this range means we are almost certainly reading the wrong column
# (e.g. a diameter in cm, or an azimuth in degrees).
MIN_PLAUSIBLE_HEIGHT_M = 0.0
MAX_PLAUSIBLE_HEIGHT_M = 40.0


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def tree_number_to_tree_id(tree_number: int) -> str:
    """
    Convert a plain tree number (e.g. 1, 7, 12) into the canonical
    project tree ID (e.g. "B21_S01", "B21_S07", "B21_S12") - :02d pads
    the number with a leading zero if it is a single digit.
    """
    return f"{TREE_ID_PREFIX}{tree_number:02d}"


def parse_czech_number(text: str) -> float:
    """
    Convert a number written in Czech format (comma as decimal separator,
    e.g. "6,43") into a Python float (6.43). The geometry file uses this
    format for every numeric value.
    """
    return float(text.replace(",", "."))


def read_geometry_lines(geometry_file_path: str) -> list:
    """
    Read the geometry file and return its non-empty lines, so the header
    rows and the data rows can be handled separately by the callers.

    If the file is not there, the .txt files that ARE in that folder are
    listed, so the correct name can be copied out of the printed list
    rather than guessed.
    """
    path = Path(geometry_file_path)
    if not path.is_file():
        print(f"ERROR: geometry file not found: {path}")
        print("Text files present in that folder:")
        for candidate in sorted(path.parent.glob("*.txt")):
            print(f"  {candidate.name}")
        raise FileNotFoundError(path)

    with open(geometry_file_path, "r", encoding="utf-8") as f:
        return [line for line in f.read().splitlines() if line.strip() != ""]


def verify_geometry_columns(lines: list) -> None:
    """
    Print what actually sits at each configured column position, and
    check the values are plausible. This does not decide anything by
    itself - it makes a wrong position obvious instead of invisible.

    Three checks are run:
      1. Print the header name found at each configured position, so it
         can be read and confirmed by eye.
      2. Check every height is inside MIN/MAX_PLAUSIBLE_HEIGHT_M.
      3. Check the ordering anchorage > level3 > level6 holds for every
         tree - the lever arms are differences of these, so a wrong
         order produces negative lever arms and negative stiffness.
    """
    print("=" * 68)
    print("VERIFICATION - geometry file column positions")
    print("=" * 68)

    name_row = lines[0].split("\t")
    data_lines = lines[2:]

    positions = {
        "tree number": GEOMETRY_COLUMN_TREE_NUMBER,
        "height anchorage": GEOMETRY_COLUMN_HEIGHT_ANCHORAGE,
        "height level III": GEOMETRY_COLUMN_HEIGHT_LEVEL3,
        "height level VI": GEOMETRY_COLUMN_HEIGHT_LEVEL6,
    }

    print("What this script expects   | pos | header name found in the file")
    print("---------------------------+-----+------------------------------")
    for label, position in positions.items():
        found_name = name_row[position].strip() if position < len(name_row) else "<POSITION OUT OF RANGE>"
        print(f"{label:<26} | {position:>3} | {found_name}")
    print()

    # --- Checks 2 and 3, per tree ---------------------------------------
    problems = []
    for line in data_lines:
        columns = line.split("\t")
        tree_number = int(parse_czech_number(columns[GEOMETRY_COLUMN_TREE_NUMBER]))
        tree_id = tree_number_to_tree_id(tree_number)

        anchorage = parse_czech_number(columns[GEOMETRY_COLUMN_HEIGHT_ANCHORAGE])
        level3 = parse_czech_number(columns[GEOMETRY_COLUMN_HEIGHT_LEVEL3])
        level6 = parse_czech_number(columns[GEOMETRY_COLUMN_HEIGHT_LEVEL6])

        for name, value in (("anchorage", anchorage), ("level3", level3), ("level6", level6)):
            if not (MIN_PLAUSIBLE_HEIGHT_M <= value <= MAX_PLAUSIBLE_HEIGHT_M):
                problems.append(f"{tree_id}: {name} = {value} m is outside the "
                                f"plausible range {MIN_PLAUSIBLE_HEIGHT_M}-{MAX_PLAUSIBLE_HEIGHT_M} m")

        if not (anchorage > level3 >= level6):
            problems.append(f"{tree_id}: expected anchorage > level3 >= level6, "
                            f"got {anchorage} / {level3} / {level6}")

    if problems:
        print("PROBLEMS FOUND - do not use the results until these are resolved:")
        for problem in problems:
            print(f"  {problem}")
    else:
        print(f"All {len(data_lines)} tree(s) passed the plausibility and "
              f"ordering checks.")
    print()


def load_tree_geometry(lines: list) -> dict:
    """
    Turn the geometry file lines into a dict keyed by canonical tree ID
    (e.g. "B21_S07"), with each value being a dict of:
      height_anchorage_m, height_level3_m, height_level6_m

    The first two lines (column names, then column indices) are only for
    reference and are skipped. Numbers use a comma as the decimal
    separator, so they are run through parse_czech_number() before use.
    """
    data_lines = lines[2:]

    geometry_by_tree = {}
    for line in data_lines:
        columns = line.split("\t")

        tree_number = int(parse_czech_number(columns[GEOMETRY_COLUMN_TREE_NUMBER]))
        tree_id = tree_number_to_tree_id(tree_number)

        geometry_by_tree[tree_id] = {
            "height_anchorage_m": parse_czech_number(columns[GEOMETRY_COLUMN_HEIGHT_ANCHORAGE]),
            "height_level3_m": parse_czech_number(columns[GEOMETRY_COLUMN_HEIGHT_LEVEL3]),
            "height_level6_m": parse_czech_number(columns[GEOMETRY_COLUMN_HEIGHT_LEVEL6]),
        }

    return geometry_by_tree


def compute_lever_arms(geometry_by_tree: dict) -> dict:
    """
    Given the per-tree heights from load_tree_geometry(), compute the two
    lever arms (vertical distances, in meters) used later to turn forces
    into bending moments:
      lever_arm_rotation_m     = height_anchorage_m - height_level6_m
      lever_arm_displacement_m = height_anchorage_m - height_level3_m

    Returns a dict keyed by tree_id, with each value being a dict of
    lever_arm_rotation_m, lever_arm_displacement_m.
    """
    lever_arms_by_tree = {}
    for tree_id, heights in geometry_by_tree.items():
        lever_arms_by_tree[tree_id] = {
            "lever_arm_rotation_m": heights["height_anchorage_m"] - heights["height_level6_m"],
            "lever_arm_displacement_m": heights["height_anchorage_m"] - heights["height_level3_m"],
        }
    return lever_arms_by_tree


def add_lever_arm_columns(stage1: pd.DataFrame, lever_arms_by_tree: dict) -> pd.DataFrame:
    """
    Add lever_arm_rotation_m and lever_arm_displacement_m columns to the
    Stage 1 table, looked up per row by the "tree" column.

    If a tree has no matching entry in lever_arms_by_tree, a warning is
    printed (naming the tree and scenario) and NaN is left in those two
    columns for that row, instead of crashing.
    """
    stage1 = stage1.copy()  # work on a copy so we never modify the original data
    rotation_values = []
    displacement_values = []

    for _, row in stage1.iterrows():
        tree = row["tree"]
        lever_arms = lever_arms_by_tree.get(tree)

        if lever_arms is None:
            print(f"WARNING: no geometry data found for tree '{tree}' "
                  f"(scenario {row['scenario']}) - leaving lever arms/moments/slopes blank.")
            rotation_values.append(float("nan"))
            displacement_values.append(float("nan"))
        else:
            rotation_values.append(lever_arms["lever_arm_rotation_m"])
            displacement_values.append(lever_arms["lever_arm_displacement_m"])

    stage1["lever_arm_rotation_m"] = rotation_values
    stage1["lever_arm_displacement_m"] = displacement_values
    return stage1


def safe_divide(numerator: float, denominator: float, tree: str, scenario, column_name: str) -> float:
    """
    Divide numerator by denominator, but return NaN (with a printed
    warning naming the tree/scenario/column) instead of crashing when the
    denominator is zero (or missing).
    """
    if pd.isna(denominator) or denominator == 0:
        print(f"WARNING: '{column_name}' is zero/missing for tree '{tree}' "
              f"(scenario {scenario}) - cannot compute this slope, leaving it blank.")
        return float("nan")
    return numerator / denominator


def compute_moments_and_slopes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute bending moments (force x lever arm) and stiffness slopes
    (moment / rotation or displacement) for every row of df.

    Uses the already unit-converted Stage 1 columns directly - no unit
    reversal needed, since FORCEx_kN is already in kN, the ROT*_deg
    columns already in degrees, and DISPLx_level3_mm already in mm.
    """
    df = df.copy()  # work on a copy so we never modify the original data

    df["moment_rotation_kNm"] = df["FORCEx_kN"] * df["lever_arm_rotation_m"]
    df["moment_displacement_kNm"] = df["FORCEx_kN"] * df["lever_arm_displacement_m"]

    slope_rotation_base = []
    slope_rotation_base2 = []
    slope_displacement = []

    for _, row in df.iterrows():
        slope_rotation_base.append(
            safe_divide(row["moment_rotation_kNm"], row["ROTSUM_base_deg"],
                        row["tree"], row["scenario"], "ROTSUM_base_deg")
        )
        slope_rotation_base2.append(
            safe_divide(row["moment_rotation_kNm"], row["ROTSUM_base2_deg"],
                        row["tree"], row["scenario"], "ROTSUM_base2_deg")
        )
        slope_displacement.append(
            safe_divide(row["moment_displacement_kNm"], row["DISPLx_level3_mm"],
                        row["tree"], row["scenario"], "DISPLx_level3_mm")
        )

    df["slope_rotation_base"] = slope_rotation_base
    df["slope_rotation_base2"] = slope_rotation_base2
    df["slope_displacement"] = slope_displacement

    return df


def round_computed_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Round every newly computed column to DECIMAL_PLACES. NaN values stay
    NaN - rounding does not affect them.
    """
    df = df.copy()  # work on a copy so we never modify the original data
    computed_columns = [
        "lever_arm_rotation_m", "lever_arm_displacement_m",
        "moment_rotation_kNm", "moment_displacement_kNm",
        "slope_rotation_base", "slope_rotation_base2", "slope_displacement",
    ]
    for column_name in computed_columns:
        df[column_name] = df[column_name].round(DECIMAL_PLACES)
    return df


def main():
    # --- Part 1: load tree geometry, verify it, compute lever arms --------
    print(f"Loading tree geometry from: {GEOMETRY_FILE_PATH}")
    lines = read_geometry_lines(GEOMETRY_FILE_PATH)

    verify_geometry_columns(lines)

    geometry_by_tree = load_tree_geometry(lines)
    print(f"Loaded geometry for {len(geometry_by_tree)} tree(s): {sorted(geometry_by_tree.keys())}")

    lever_arms_by_tree = compute_lever_arms(geometry_by_tree)

    # --- Part 2: load Stage 1 data and merge in the lever arms --------------
    print(f"Loading Stage 1 data from: {STAGE1_INPUT_PATH}")
    stage1 = pd.read_csv(STAGE1_INPUT_PATH)
    print(f"Stage 1 table: {len(stage1)} row(s)")

    combined = add_lever_arm_columns(stage1, lever_arms_by_tree)

    # --- Part 3: compute moments and slopes ----------------------------------
    combined = compute_moments_and_slopes(combined)
    combined = round_computed_columns(combined)

    # --- Part 4: sort, order columns, print and save -------------------------
    combined = combined.sort_values(["tree", "scenario"]).reset_index(drop=True)

    output_columns = [
        "tree", "ansys_no", "scenario", "model", "source_file",
        "FORCEx_kN", "FORCEz_kN",
        "ROTSUM_base_deg", "ROTx_deg", "ROTy_deg", "ROTz_deg", "ROTSUM_base2_deg",
        "DISPLx_anchorage_mm", "DISPLx_level3_mm", "DISPLx_MAX_mm",
        "lever_arm_rotation_m", "lever_arm_displacement_m",
        "moment_rotation_kNm", "moment_displacement_kNm",
        "slope_rotation_base", "slope_rotation_base2", "slope_displacement",
    ]
    combined = combined[output_columns]

    print(f"\nCombined table ({len(combined)} row(s)):")
    print(combined.to_string(index=False))

    combined.to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved combined table to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
