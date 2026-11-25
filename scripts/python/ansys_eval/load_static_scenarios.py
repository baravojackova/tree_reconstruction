"""
load_static_scenarios.py
-------------------------
Stage 1 of the ANSYS post-processing workflow.

Loads all ANSYS static-analysis output files for multiple trees/scenarios
from one shared folder into a single combined table, converts the raw
units into more convenient ones, and saves the result as a CSV file for
manual review.

Each input file (e.g. "out_static_bk8_scenario1.txt") looks like this:

    FORCEx | FORCEz | ROTSUM base | ROTx | ...
    11.6036   4.4854   0.0000125   0.0000001738   ...

  - Line 1 lists the 10 column headers, separated by "|". Some headers
    contain spaces (e.g. "ROTSUM base"), so they must NOT be split on
    spaces - only on "|".
  - Line 2 has exactly one data row: the 10 numeric values, separated by
    plain whitespace (no "|" there).

TWO THINGS THAT ARE EASY TO GET WRONG (read before editing):

  1. The number after "bk" in the filename is the tree's POSITION in the
     description file (1, 2, 3, ...), NOT the tree's own number. E.g.
     "bk3" is position 3, which is tree 7, i.e. "B21_S07". This script
     therefore reads the description file, builds the mapping
     position -> tree number, and never guesses.

  2. The scenario number is just a run index. What it means (which
     reconstruction model was used) lives in SCENARIO_MODELS below and
     must be kept up to date by hand. An unknown scenario number does
     not crash the script - it gets a placeholder label and a warning.

How to use it:
  - Edit the values in the "CONFIGURATION" section below if needed.
  - Run the script from a terminal in VS Code:  python load_static_scenarios.py
  - Check the printed table and the generated CSV file before using it
    for anything further.
"""

import math
import re
from pathlib import Path

import pandas as pd

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folder holding the ANSYS output files AND the description file. The r""
# prefix means "raw string": it stops Python from treating backslashes in
# a Windows path as escape characters, so the path can be pasted exactly
# as Explorer shows it.
INPUT_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\ansys"

# Name of the description file (one row per tree), inside INPUT_FOLDER.
# It is deliberately read from the SAME folder that ANSYS itself reads,
# so the two can never drift apart. Needed here ONLY to translate the
# position number in the filename into a real tree number.
GEOMETRY_FILE_NAME = "Popis_Babice_VSE_09012025_ans_lev6_0_azimcor.txt"

# Full path, assembled from the two constants above. Path() / "name" is
# pathlib's way of joining a folder and a filename with the correct
# separator - it is the same as writing folder + "\" + name, but safe.
GEOMETRY_FILE_PATH = Path(INPUT_FOLDER) / GEOMETRY_FILE_NAME

FILE_PATTERN = "out_static_bk*_scenario*.txt"   # glob pattern - ONLY static files

# Folder this script's output CSV is written into.
OUTPUT_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\out"
OUTPUT_PATH = Path(OUTPUT_FOLDER) / "static_scenarios_combined.csv"
DECIMAL_PLACES = 6

# Canonical tree ID used across the whole reconstruction project:
# TREE_ID_PREFIX + zero-padded tree number, e.g. "B21_S" + "07" -> "B21_S07".
TREE_ID_PREFIX = "B21_S"

# Positions (0-based) in the description file. Position 0 is the row
# order, which equals the number after "bk" in the ANSYS filenames;
# position 2 is the tree's own number.
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
# and reported as a warning.
SCENARIO_MODELS = {
    1: "adtree_r5mm",
    2: "treeqsm_pd002",
    3: "adtree_r5mm_no acel",
    4: "treeqsm_pd002_no_acel",
    999: "adtree_r5_pilot",   # high cylinder count run kept from the cost-model session
}

# The 10 column names as they appear in the raw files, in their fixed
# order. Used only to sanity-check that a file has the columns we expect.
EXPECTED_COLUMN_NAMES = [
    "FORCEx", "FORCEz", "ROTSUM base", "ROTx", "ROTy", "ROTz",
    "ROTSUM base2", "DISPLx anchorage", "DISPLx level3", "DISPLx MAX",
]

# Regular expression used to pull the position number and scenario number
# out of a filename like "out_static_bk8_scenario1.txt".
FILENAME_PATTERN = re.compile(r"out_static_bk(\d+)_scenario(\d+)\.txt")


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def parse_czech_number(text: str) -> float:
    """
    Convert a number written in Czech format (comma as decimal separator,
    e.g. "6,43") into a Python float (6.43). The description file uses
    this format for every numeric value.
    """
    return float(text.replace(",", "."))


def tree_number_to_tree_id(tree_number: int) -> str:
    """
    Convert a plain tree number (e.g. 1, 7, 12) into the canonical
    project tree ID (e.g. "B21_S01", "B21_S07", "B21_S12").

    :02d pads the number with a leading zero if it is a single digit.
    """
    return f"{TREE_ID_PREFIX}{tree_number:02d}"


def load_ansys_order_to_tree_id(geometry_file_path: str) -> dict:
    """
    Read the description file and return a dict that translates the
    number used in the ANSYS filenames into a canonical tree ID:

        {1: "B21_S01", 2: "B21_S04", 3: "B21_S07", ...}

    The file has 2 header/reference lines (column names, then column
    indices) which are skipped, followed by one tab-separated data row
    per tree. Numbers use a comma as the decimal separator, so they go
    through parse_czech_number() first.

    Why this exists: the ANSYS filenames carry the ROW POSITION, not the
    tree number. Without this lookup, "bk3" would be read as tree 3,
    which does not exist, and "bk8" would be read as tree 8, which does
    exist but is the wrong tree - a silent, invisible error.
    """
    # If the filename is even slightly off, list what IS in the folder
    # instead of only reporting "file not found" - the correct name can
    # then be copied straight out of the printed list.
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
    placeholder name so the row still appears in the output table, and
    the caller can print a warning.
    """
    return SCENARIO_MODELS.get(scenario_number, f"scenario{scenario_number}")


def parse_static_file(file_path: Path) -> dict:
    """
    Read one ANSYS static output file and return its 10 column values as
    a dict {column_name: value}, with values converted to float.

    Returns None (and prints a warning) if the file does not contain
    exactly 10 header names and 10 numeric values, so a single malformed
    file does not crash the whole run.
    """
    # Read the file and drop any blank lines, so we always work with
    # "line 1 = headers, line 2 = values" even if there is a trailing
    # empty line at the end of the file.
    with open(file_path, "r", encoding="utf-8") as f:
        lines = [line for line in f.read().splitlines() if line.strip() != ""]

    if len(lines) < 2:
        print(f"WARNING: '{file_path.name}' does not have both a header line "
              f"and a data line - skipping this file.")
        return None

    header_line = lines[0]
    data_line = lines[1]

    # Headers are separated by "|" only - some of them contain spaces
    # (e.g. "ROTSUM base"), so splitting on whitespace would wrongly break
    # those apart into two columns.
    column_names = [name.strip() for name in header_line.split("|")]

    # The data line, on the other hand, has no "|" at all - the 10 values
    # are separated by plain whitespace.
    raw_values = data_line.split()

    if len(column_names) != 10 or len(raw_values) != 10:
        print(f"WARNING: '{file_path.name}' has {len(column_names)} column "
              f"name(s) and {len(raw_values)} value(s) - expected 10 of "
              f"each - skipping this file.")
        return None

    values = [float(v) for v in raw_values]

    # Zip the names and values together positionally (same order in both
    # lists) to build one dict per file, e.g. {"FORCEx": 11.6036, ...}.
    return dict(zip(column_names, values))


def extract_order_and_scenario(file_path: Path) -> tuple:
    """
    Extract the position number and scenario number from a filename like
    "out_static_bk8_scenario1.txt".

    Returns (ansys_order, scenario_number), e.g. (8, 1).
    Returns (None, None) if the filename does not match the expected
    pattern.

    Note what this function does NOT do: it does not convert the number
    into a tree ID. That conversion needs the description file and is
    done by the caller via the order_to_tree_id lookup.
    """
    match = FILENAME_PATTERN.match(file_path.name)
    if match is None:
        return None, None

    ansys_order = int(match.group(1))
    scenario_number = int(match.group(2))

    return ansys_order, scenario_number


def convert_units_static_row(raw_values: dict) -> dict:
    """
    Convert one file's raw column values (N, radians, meters) into the
    more convenient units used in the combined output table (kN, degrees,
    millimeters), and rename the columns to reflect the new units.

    Every converted value is rounded to DECIMAL_PLACES.
    """
    converted = {}

    # --- Force columns: N -> kN (divide by 1000) -----------------------
    converted["FORCEx_kN"] = round(raw_values["FORCEx"] / 1000, DECIMAL_PLACES)
    converted["FORCEz_kN"] = round(raw_values["FORCEz"] / 1000, DECIMAL_PLACES)

    # --- Rotation columns: radians -> degrees ---------------------------
    converted["ROTSUM_base_deg"] = round(math.degrees(raw_values["ROTSUM base"]), DECIMAL_PLACES)
    converted["ROTx_deg"] = round(math.degrees(raw_values["ROTx"]), DECIMAL_PLACES)
    converted["ROTy_deg"] = round(math.degrees(raw_values["ROTy"]), DECIMAL_PLACES)
    converted["ROTz_deg"] = round(math.degrees(raw_values["ROTz"]), DECIMAL_PLACES)
    converted["ROTSUM_base2_deg"] = round(math.degrees(raw_values["ROTSUM base2"]), DECIMAL_PLACES)

    # --- Displacement columns: meters -> millimeters (multiply by 1000) -
    converted["DISPLx_anchorage_mm"] = round(raw_values["DISPLx anchorage"] * 1000, DECIMAL_PLACES)
    converted["DISPLx_level3_mm"] = round(raw_values["DISPLx level3"] * 1000, DECIMAL_PLACES)
    converted["DISPLx_MAX_mm"] = round(raw_values["DISPLx MAX"] * 1000, DECIMAL_PLACES)

    return converted


def main():
    # --- Step 0: build the position -> tree ID lookup --------------------
    print(f"Loading tree list from: {GEOMETRY_FILE_PATH}")
    order_to_tree_id = load_ansys_order_to_tree_id(GEOMETRY_FILE_PATH)
    print(f"Mapping built for {len(order_to_tree_id)} tree(s):")
    for ansys_order in sorted(order_to_tree_id):
        print(f"  bk{ansys_order}  ->  {order_to_tree_id[ansys_order]}")
    print()

    # --- Step 1: find all matching files -------------------------------
    # Path.glob() finds every file in INPUT_FOLDER whose name matches
    # FILE_PATTERN.
    input_folder = Path(INPUT_FOLDER)
    matching_files = sorted(input_folder.glob(FILE_PATTERN))
    print(f"Found {len(matching_files)} file(s) matching '{FILE_PATTERN}' in '{INPUT_FOLDER}'")

    if len(matching_files) == 0:
        print("No files found - check INPUT_FOLDER and FILE_PATTERN.")
        return

    # --- Step 2: parse every file into one combined row -----------------
    rows = []
    unknown_scenarios = set()

    for file_path in matching_files:
        raw_values = parse_static_file(file_path)
        if raw_values is None:
            continue  # parse_static_file already printed a warning

        ansys_order, scenario = extract_order_and_scenario(file_path)
        if ansys_order is None:
            print(f"WARNING: '{file_path.name}' does not match the expected "
                  f"filename pattern 'out_static_bk<N>_scenario<M>.txt' - "
                  f"skipping this file.")
            continue

        # Translate the position number into a real tree ID. A position
        # that is not in the description file means the filename and the
        # description file disagree - refuse to guess, skip and report.
        tree = order_to_tree_id.get(ansys_order)
        if tree is None:
            print(f"WARNING: '{file_path.name}' refers to position bk{ansys_order}, "
                  f"which does not exist in the description file - skipping "
                  f"this file.")
            continue

        model = scenario_to_model_name(scenario)
        if scenario not in SCENARIO_MODELS:
            unknown_scenarios.add(scenario)

        converted_values = convert_units_static_row(raw_values)

        row = {
            "tree": tree,
            "ansys_no": ansys_order,
            "scenario": scenario,
            "model": model,
            "source_file": file_path.name,
        }
        row.update(converted_values)
        rows.append(row)

    if len(rows) == 0:
        print("No files could be parsed successfully - nothing to save.")
        return

    if unknown_scenarios:
        print()
        print(f"WARNING: scenario number(s) {sorted(unknown_scenarios)} are not "
              f"listed in SCENARIO_MODELS. They were loaded with a placeholder "
              f"name - add them to SCENARIO_MODELS so the plots get a proper "
              f"label.")

    # --- Step 3-4: build the combined table and sort it -----------------
    combined = pd.DataFrame(rows)
    combined = combined.sort_values(["tree", "scenario"]).reset_index(drop=True)

    # --- Step 5: show and save the result --------------------------------
    print(f"\nCombined table ({len(combined)} row(s)):")
    print(combined.to_string(index=False))

    print("\nCoverage (tree x model):")
    coverage = combined.pivot_table(index="tree", columns="model",
                                    values="FORCEx_kN", aggfunc="count",
                                    fill_value=0)
    print(coverage.to_string())

    combined.to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved combined table to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
