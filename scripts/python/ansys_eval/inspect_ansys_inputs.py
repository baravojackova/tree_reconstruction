"""
inspect_ansys_inputs.py
------------------------
READ-ONLY diagnostic. Stage 0 of the ANSYS-vs-experiment comparison.

This script writes nothing, creates nothing and deletes nothing. It only
reads two things and prints what it found, so that the column positions
and the tree numbering can be confirmed BEFORE any of the real scripts
(load_static_scenarios.py / compute_slopes.py) are pointed at the new
data.

What it does, in plain terms:
  Part A - looks into the ANSYS folder and lists every output file it
           finds, split into static / modal / harmonic, and prints a
           table of which ANSYS tree numbers and which scenario numbers
           actually exist on disk.
  Part B - opens the geometry description file and prints, for the first
           few rows, EVERY column together with its position number
           (0, 1, 2, ...). This is what lets us decide which position
           holds the anchorage height, level III and level VI.
  Part C - builds the mapping "ANSYS number -> tree number -> BK id"
           from column 0 and column 2, and prints it, so the mapping can
           be eyeballed against reality before it is used anywhere.

How to use it:
  - Edit the two paths in the CONFIGURATION section below.
  - Run it from a terminal in VS Code:  python inspect_ansys_inputs.py
  - Copy the whole printed output back into the chat.
"""

import re
from pathlib import Path
import sys

# This script lives in scripts/python/ansys_eval, but paths.py lives one
# folder up (scripts/python). Python only looks for imports next to the
# script, so we first add the folder above to its search list.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Folder that holds the ANSYS output files (out_static_*, out_modal_*, ...).
# The r"" prefix means "raw string" - it stops Python from treating the
# backslashes in a Windows path as escape characters, so the path can be
# pasted in exactly as Explorer shows it.
ANSYS_FOLDER = str(paths.ANSYS_DIR)

# Full path to the geometry / description file.
# NOTE: this is a DATA file - Git ignores it, so on a new computer (e.g. the
# Mac) it must be copied by hand into scripts/python/ before running.
GEOMETRY_FILE_PATH = str(paths.PYTHON_DIR / "Popis_Babice_VSE_09012025_ans_lev6_0.txt")

# How many data rows of the geometry file to print in full detail
# (column by column). 2 is usually enough to see the pattern.
N_DETAIL_ROWS = 2

# Column positions we BELIEVE hold the ANSYS number and the tree number.
# Part C prints the mapping built from these, so they can be verified.
# 0 = first column, 2 = third column (Python counts from zero).
COLUMN_ANSYS_ORDER = 0
COLUMN_TREE_NUMBER = 2


# =========================================================================
# CODE - nothing below needs editing
# =========================================================================

# A regular expression that recognises filenames of the form
#   out_static_bk1_scenario3.txt
# and captures three pieces: the analysis type, the bk number, the
# scenario number. Anything that does not match is reported separately
# rather than silently ignored.
FILENAME_PATTERN = re.compile(
    r"out_(static|modal|harmonic)_bk(\d+)_scenario(\d+)\.txt",
    re.IGNORECASE,
)


def parse_czech_number(text: str):
    """
    Convert a number written with a comma as the decimal separator
    ("6,43") into a float (6.43). Returns None if the text is not a
    number at all - this script must never crash on a text column.
    """
    try:
        return float(text.strip().replace(",", "."))
    except ValueError:
        return None


def inspect_ansys_folder(folder_path: str) -> None:
    """
    Part A: list the ANSYS output files and summarise which tree numbers
    and scenario numbers are present, per analysis type.
    """
    print("=" * 72)
    print("PART A - ANSYS output files")
    print("=" * 72)

    folder = Path(folder_path)
    if not folder.is_dir():
        print(f"ERROR: folder does not exist: {folder_path}")
        return

    all_txt_files = sorted(folder.glob("*.txt"))
    print(f"Folder: {folder_path}")
    print(f"Total .txt files found: {len(all_txt_files)}")
    print()

    # found[analysis_type] = set of (bk_number, scenario_number) pairs
    found = {"static": set(), "modal": set(), "harmonic": set()}
    unmatched = []

    for file_path in all_txt_files:
        match = FILENAME_PATTERN.match(file_path.name)
        if match is None:
            unmatched.append(file_path.name)
            continue
        analysis_type = match.group(1).lower()
        bk_number = int(match.group(2))
        scenario_number = int(match.group(3))
        found[analysis_type].add((bk_number, scenario_number))

    for analysis_type in ("static", "modal", "harmonic"):
        pairs = sorted(found[analysis_type])
        bk_numbers = sorted({bk for bk, _ in pairs})
        scenarios = sorted({sc for _, sc in pairs})
        print(f"--- {analysis_type} ---")
        print(f"  files            : {len(pairs)}")
        print(f"  ANSYS bk numbers : {bk_numbers}")
        print(f"  scenario numbers : {scenarios}")

        # A small presence grid: rows = bk number, columns = scenario.
        # "X" means the file exists, "." means it is missing. This makes
        # gaps in the run plan visible at a glance.
        if pairs:
            header = "      bk |" + "".join(f" s{sc}" for sc in scenarios)
            print(header)
            for bk in bk_numbers:
                cells = "".join("  X" if (bk, sc) in found[analysis_type] else "  ."
                                for sc in scenarios)
                print(f"      {bk:>2} |{cells}")
        print()

    if unmatched:
        print("Files that did NOT match the expected name pattern")
        print("(out_<type>_bk<N>_scenario<M>.txt) - listed, not touched:")
        for name in unmatched:
            print(f"  {name}")
        print()


def read_geometry_lines(geometry_file_path: str):
    """
    Read the geometry file and return its non-empty lines. Returns None
    if the file cannot be found, so the caller can report it cleanly.
    """
    path = Path(geometry_file_path)
    if not path.is_file():
        print(f"ERROR: geometry file does not exist: {geometry_file_path}")
        return None

    with open(path, "r", encoding="utf-8") as f:
        return [line for line in f.read().splitlines() if line.strip() != ""]


def inspect_geometry_columns(lines) -> None:
    """
    Part B: print the geometry file column by column, with the position
    number in front of every value.

    The file convention is: row 1 = column names, row 2 = column indices,
    row 3 onwards = one data row per tree. Columns are separated by TAB
    characters, so the line is split on "\t" and NOT on spaces (several
    column names contain spaces).
    """
    print("=" * 72)
    print("PART B - geometry file, column positions")
    print("=" * 72)

    name_row = lines[0].split("\t")
    index_row = lines[1].split("\t")
    data_rows = lines[2:]

    print(f"Columns in the name row : {len(name_row)}")
    print(f"Columns in the index row: {len(index_row)}")
    print(f"Data rows (trees)       : {len(data_rows)}")
    print()

    # Print name row and index row side by side, one column per output
    # line, prefixed by the Python position. This is the table we need in
    # order to fix the hardcoded positions in compute_slopes.py.
    print("pos | name (row 1)                          | index (row 2)")
    print("----+---------------------------------------+--------------")
    for position, name in enumerate(name_row):
        index_value = index_row[position] if position < len(index_row) else ""
        print(f"{position:>3} | {name.strip()[:37]:<37} | {index_value.strip()}")
    print()

    # Print the first few data rows in the same column-by-column form, so
    # the names above can be checked against real values (e.g. a height
    # that is 0.00 for most trees is level VI, the tree base).
    for row_number, line in enumerate(data_rows[:N_DETAIL_ROWS], start=1):
        columns = line.split("\t")
        print(f"--- data row {row_number} ({len(columns)} columns) ---")
        for position, value in enumerate(columns):
            name = name_row[position].strip()[:30] if position < len(name_row) else "?"
            print(f"  {position:>3} | {name:<30} | {value.strip()}")
        print()


def inspect_tree_mapping(lines) -> None:
    """
    Part C: build and print the mapping from the ANSYS order number
    (column COLUMN_ANSYS_ORDER) to the tree number (column
    COLUMN_TREE_NUMBER) and the BK id used in the experimental dataset.

    Nothing is saved - this is only here so the mapping can be checked
    by eye before any script relies on it.
    """
    print("=" * 72)
    print("PART C - ANSYS number -> tree number -> BK id")
    print("=" * 72)
    print(f"(built from column {COLUMN_ANSYS_ORDER} and column {COLUMN_TREE_NUMBER})")
    print()

    data_rows = lines[2:]

    print("ansys_no | tree_no | bk_id")
    print("---------+---------+------")
    for line in data_rows:
        columns = line.split("\t")

        ansys_value = parse_czech_number(columns[COLUMN_ANSYS_ORDER]) \
            if COLUMN_ANSYS_ORDER < len(columns) else None
        tree_value = parse_czech_number(columns[COLUMN_TREE_NUMBER]) \
            if COLUMN_TREE_NUMBER < len(columns) else None

        if ansys_value is None or tree_value is None:
            print(f"  ?      |   ?     |  ?     <- could not read this row: "
                  f"{columns[:4]}")
            continue

        ansys_number = int(ansys_value)
        tree_number = int(tree_value)
        # :02d pads a single digit with a leading zero, so 4 -> "BK04".
        bk_id = f"BK{tree_number:02d}"

        # Flag the rows where the two numbers differ - those are exactly
        # the trees that the old "bk<N> is the tree number" assumption
        # would have got wrong.
        flag = "" if ansys_number == tree_number else "   <- differs"
        print(f"   {ansys_number:>3}   |   {tree_number:>3}   | {bk_id}{flag}")
    print()


def main():
    inspect_ansys_folder(ANSYS_FOLDER)

    lines = read_geometry_lines(GEOMETRY_FILE_PATH)
    if lines is None:
        return
    if len(lines) < 3:
        print("ERROR: the geometry file has fewer than 3 non-empty lines - "
              "expected name row, index row and at least one data row.")
        return

    inspect_geometry_columns(lines)
    inspect_tree_mapping(lines)

    print("Done. Nothing was written or modified.")


if __name__ == "__main__":
    main()
