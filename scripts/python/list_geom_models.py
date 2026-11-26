r"""
list_geom_models.py
--------------------
Lists every geometry export file for one tree and how many cylinders it
contains, so the available models can be compared at a glance before
choosing what to run in ANSYS.

With respect to the geometry files it inspects (under ansys_geom\ and
scripts\matlab\), this script is READ-ONLY - it only ever opens them for
reading, never writes, moves or deletes anything there. Its only WRITE
side effect is its own CSV output (see WRITE_CSV below), which lands in a
completely separate folder (ansys_eval\out\), never among the files being
inspected.

It looks in two places, because the two reconstruction branches export to
different folders:
  - ansys_geom\<TREE>\   - AdTree / AdQSM exports, one subfolder per tree,
                           e.g. geom_B21_S21_r20mm_adqsm070_seg50-500-k80-wlen.txt
  - scripts\matlab\      - TreeQSM exports, all trees mixed in one folder,
                           e.g. geom_B21_S13_man_pd06-002-10_mo9_sr010_ri1_1.txt

WHY THE CYLINDER COUNT IS NOT SIMPLY THE LINE COUNT:
the file may or may not start with a header line. One extra line is
invisible among tens of thousands, but it shifts every count by one. This
script therefore tries to read the first line as numbers: if that fails,
the line is treated as a header, reported as such, and excluded from the
count. Blank lines are excluded too.

How to use it:
  - Set TREE_NAME below to the tree you want, then run:
        python list_geom_models.py
  - Or set SCAN_ALL_TREES = True to list every tree found.
  - The same table is also written as a CSV (see WRITE_CSV, OUT_FOLDER
    below), so the n_rad values can be reopened, sorted and pasted from
    later instead of scrolled back to in a terminal.
"""

import csv
from pathlib import Path

# paths.py (same folder as this script) knows where the project lives on
# THIS computer, so no folder path has to be typed in by hand here.
import paths

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Which tree to list. Only used when SCAN_ALL_TREES is False.
TREE_NAME = "B21_S04"

# True  = ignore TREE_NAME and list every tree found in both folders.
# False = list only TREE_NAME.
SCAN_ALL_TREES = True

# The r"" prefix means "raw string": it stops Python from treating
# backslashes in a Windows path as escape characters.
ANSYS_GEOM_ROOT = str(paths.PYTHON_DIR / "ansys_geom")
MATLAB_FOLDER = str(paths.MATLAB_DIR)

# Which files count as geometry exports.
FILE_GLOB = "geom_*.txt"

# Print the first data line of each file as well, to make it obvious what
# the columns look like. Set to False once you trust the format.
SHOW_FIRST_DATA_LINE = False

# Maximum width of the model-name column in the printed table. Names
# longer than this are still printed in full, the table just gets wider.
MIN_NAME_COLUMN_WIDTH = 40

# Folder the CSV is written into. Shared with the rest of the ANSYS-eval
# output, so the numbers sit next to everything else that gets typed into
# or compared against ANSYS.
OUT_FOLDER = str(paths.ANSYS_EVAL_DIR / "out")

# True  = also write the printed table to a CSV in OUT_FOLDER, so the
#         n_rad values can be reopened, sorted and pasted from later
#         instead of being scrolled back to in a terminal.
# False = print only, same behaviour as before the CSV was added.
WRITE_CSV = True


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def first_field_as_int(line: str):
    """
    Return the first whitespace-separated field of the line as an int, or
    None if it is not an integer. Used to read the row index that the
    geometry files carry in column 0.
    """
    fields = line.split()
    if len(fields) == 0:
        return None
    try:
        return int(fields[0])
    except ValueError:
        return None


def count_rows_for_ansys(file_path: Path) -> dict:
    """
    Read one geometry file and work out the row count that ANSYS uses,
    which is NOT simply the number of lines.

    The file starts with an index row "0 1 2 3 ... 10" that names the
    column positions. ANSYS does not count that row: the data rows are
    numbered from 1 in column 0, and the row count is the LAST index in
    column 0.

    Counting lines instead would be wrong by exactly one - and because
    the index row is made of numbers, a "is the first line numeric?"
    header test does not catch it either. So the count is taken straight
    from column 0, the same way ANSYS reads it.

    Returns a dict:
      n_rad          - the last index in column 0, i.e. what ANSYS uses
      n_data_lines   - how many data lines there actually are
      n_lines_total  - every line in the file
      has_index_row  - True if the first line starts with 0
      n_columns      - fields on the first data line
      index_gap      - True if the indices are not exactly 1..n_rad
      first_data_line
    """
    try:
        with open(file_path, "r", encoding="utf-8", errors="replace") as f:
            raw_lines = f.read().splitlines()
    except OSError as error:
        print(f"WARNING: could not read {file_path.name}: {error}")
        return None

    non_blank = [line for line in raw_lines if line.strip() != ""]

    if len(non_blank) == 0:
        return {"n_rad": 0, "n_data_lines": 0, "n_lines_total": len(raw_lines),
                "has_index_row": False, "n_columns": 0, "index_gap": False,
                "first_data_line": ""}

    # The index row is the one whose column 0 is 0. Everything after it
    # is data, numbered from 1.
    has_index_row = first_field_as_int(non_blank[0]) == 0
    data_lines = non_blank[1:] if has_index_row else non_blank

    indices = [first_field_as_int(line) for line in data_lines]
    valid_indices = [i for i in indices if i is not None]

    n_rad = max(valid_indices) if valid_indices else 0
    first_data_line = data_lines[0] if data_lines else ""
    n_columns = len(first_data_line.split()) if first_data_line else 0

    # If the indices are not exactly 1, 2, ... n_rad, then a row is
    # missing or duplicated. ANSYS would still read n_rad rows, but they
    # would not be the rows in this file - worth catching loudly.
    index_gap = sorted(valid_indices) != list(range(1, n_rad + 1))

    return {
        "n_rad": n_rad,
        "n_data_lines": len(data_lines),
        "n_lines_total": len(raw_lines),
        "has_index_row": has_index_row,
        "n_columns": n_columns,
        "index_gap": index_gap,
        "first_data_line": first_data_line,
    }


def model_label(file_name: str, tree_name: str) -> str:
    """
    Turn a full filename into a short model label by stripping the
    "geom_" prefix, the tree name and the ".txt" extension:

      "geom_B21_S21_r20mm_adqsm070_seg50-500-k80-wlen.txt"
        -> "r20mm_adqsm070_seg50-500-k80-wlen"

    If the expected prefix is not there, the filename is returned
    unchanged rather than being mangled.
    """
    label = file_name
    if label.endswith(".txt"):
        label = label[:-4]
    prefix = f"geom_{tree_name}_"
    if label.startswith(prefix):
        label = label[len(prefix):]
    elif label.startswith("geom_"):
        label = label[len("geom_"):]
    return label


def collect_files_for_tree(tree_name: str) -> list:
    """
    Find every geometry export belonging to one tree, in both folders.

    Returns a list of (source_label, file_path) tuples, where
    source_label is "ansys_geom" or "matlab", so the printed table can
    show which branch each model came from.

    The ansys_geom folder has one subfolder per tree, so the tree name
    selects the folder. The matlab folder mixes all trees together, so
    the tree name has to be matched inside the filename instead.
    """
    found = []

    # --- Branch 1: ansys_geom\<TREE>\ --------------------------------
    tree_folder = Path(ANSYS_GEOM_ROOT) / tree_name
    if tree_folder.is_dir():
        for file_path in sorted(tree_folder.glob(FILE_GLOB)):
            found.append(("ansys_geom", file_path))
    else:
        print(f"Note: no folder {tree_folder} - skipping the ansys_geom branch "
              f"for this tree.")

    # --- Branch 2: scripts\matlab\, filtered by tree name in filename -
    matlab_folder = Path(MATLAB_FOLDER)
    if matlab_folder.is_dir():
        for file_path in sorted(matlab_folder.glob(FILE_GLOB)):
            if tree_name in file_path.name:
                found.append(("matlab", file_path))
    else:
        print(f"Note: folder {matlab_folder} does not exist - skipping the "
              f"matlab branch.")

    return found


def find_all_tree_names() -> list:
    """
    Work out which trees have geometry exports at all, by looking at the
    subfolder names under ansys_geom and at the tree names embedded in
    the matlab filenames. Used only when SCAN_ALL_TREES is True.
    """
    tree_names = set()

    root = Path(ANSYS_GEOM_ROOT)
    if root.is_dir():
        for sub in root.iterdir():
            if sub.is_dir():
                tree_names.add(sub.name)

    matlab_folder = Path(MATLAB_FOLDER)
    if matlab_folder.is_dir():
        for file_path in matlab_folder.glob(FILE_GLOB):
            # Expected shape: geom_<TREE>_<rest>.txt, where <TREE> is two
            # underscore-separated parts, e.g. "B21_S13".
            parts = file_path.stem.split("_")
            if len(parts) >= 3 and parts[0] == "geom":
                tree_names.add(f"{parts[1]}_{parts[2]}")

    return sorted(tree_names)


def report_tree(tree_name: str) -> list:
    """
    Print the table of models and cylinder counts for one tree, and
    return the same rows as a list of dicts (one per geometry file, empty
    list if there was nothing to report), so the caller can also write
    them to a CSV without re-deriving anything - same numbers, same order.
    """
    print("=" * 78)
    print(f"TREE: {tree_name}")
    print("=" * 78)

    files = collect_files_for_tree(tree_name)
    if len(files) == 0:
        print("No geometry files found for this tree in either folder.")
        print()
        return []

    # Gather the numbers first, so the column widths can be sized to the
    # actual content rather than guessed.
    entries = []
    for source_label, file_path in files:
        info = count_rows_for_ansys(file_path)
        if info is None:
            continue
        entries.append({
            "tree": tree_name,
            "source": source_label,
            "path": file_path,
            "label": model_label(file_path.name, tree_name),
            **info,
        })

    if len(entries) == 0:
        print("No readable geometry files for this tree.")
        print()
        return []

    # Sort once, here, so the printed table and the returned (later:
    # CSV-written) rows are guaranteed to be in the same order.
    entries = sorted(entries, key=lambda e: e["n_rad"])

    name_width = max(MIN_NAME_COLUMN_WIDTH,
                     max(len(entry["label"]) for entry in entries))

    print(f"{'model':<{name_width}} | {'n_rad':>9} | {'lines':>9} | "
          f"{'cols':>4} | source")
    print("-" * name_width + "-+-" + "-" * 9 + "-+-" + "-" * 9 + "-+-"
          + "-" * 4 + "-+-------")

    for entry in entries:
        print(f"{entry['label']:<{name_width}} | {entry['n_rad']:>9,} | "
              f"{entry['n_lines_total']:>9,} | {entry['n_columns']:>4} | "
              f"{entry['source']}")

    total = sum(entry["n_rad"] for entry in entries)
    print("-" * (name_width + 40))
    print(f"{'TOTAL across ' + str(len(entries)) + ' model(s)':<{name_width}} | "
          f"{total:>9,} |")
    print()
    print("n_rad = last index in column 0, i.e. the value ANSYS needs.")
    print("lines = physical lines in the file, including the '0 1 2 ...' index row.")
    print()

    # --- Problems worth stopping for ------------------------------------
    for entry in entries:
        if not entry["has_index_row"]:
            print(f"WARNING: {entry['label']} does not start with an index row "
                  f"('0 1 2 ...'). Its first line was counted as data - check "
                  f"the file before using n_rad = {entry['n_rad']}.")
        if entry["index_gap"]:
            print(f"WARNING: {entry['label']} has gaps or duplicates in column 0 "
                  f"- the indices are not exactly 1..{entry['n_rad']} "
                  f"({entry['n_data_lines']} data lines for a last index of "
                  f"{entry['n_rad']}). ANSYS would read n_rad rows, but not the "
                  f"rows in this file.")

    column_counts = {entry["n_columns"] for entry in entries}
    if len(column_counts) > 1:
        print(f"WARNING: these files do not all have the same number of columns: "
              f"{sorted(column_counts)} - check before using them together.")

    if SHOW_FIRST_DATA_LINE:
        print("\nFirst data line of each file:")
        for entry in entries:
            print(f"  {entry['label']}")
            print(f"    {entry['first_data_line']}")
    print()

    return entries


# CSV columns, in this exact order - see write_geom_models_csv() below.
CSV_COLUMNS = ["tree", "model", "source", "n_rad", "n_lines_total", "n_columns",
               "has_index_row", "index_gap", "file_name"]


def write_geom_models_csv(entries: list, csv_path: Path) -> None:
    """
    Write the same rows already printed by report_tree() to a CSV at
    csv_path, one row per geometry file, in the same order they were
    printed in (entries is expected to already be sorted the way the
    caller wants the CSV sorted).

    Always overwrites csv_path rather than appending - appending would
    silently duplicate every row each time the script is re-run. This is
    the ONLY file this script ever writes; the geometry files themselves
    (under ansys_geom/ and scripts/matlab/) are only ever opened for
    reading, never touched.
    """
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(CSV_COLUMNS)
        for entry in entries:
            writer.writerow([
                entry["tree"],
                entry["label"],
                entry["source"],
                entry["n_rad"],
                entry["n_lines_total"],
                entry["n_columns"],
                entry["has_index_row"],
                entry["index_gap"],
                entry["path"].name,
            ])

    print(f"CSV written to: {csv_path}")


def main():
    if SCAN_ALL_TREES:
        tree_names = find_all_tree_names()
        print(f"SCAN_ALL_TREES is on - found {len(tree_names)} tree(s): {tree_names}\n")
    else:
        tree_names = [TREE_NAME]

    all_entries = []
    for tree_name in tree_names:
        all_entries.extend(report_tree(tree_name))

    if WRITE_CSV:
        if SCAN_ALL_TREES:
            csv_path = Path(OUT_FOLDER) / "geom_models_all.csv"
            # One file for every tree - sort by tree first, then n_rad,
            # so each tree's own rows stay grouped and ascending, same as
            # report_tree() already sorts within a single tree.
            csv_entries = sorted(all_entries, key=lambda e: (e["tree"], e["n_rad"]))
        else:
            csv_path = Path(OUT_FOLDER) / f"geom_models_{TREE_NAME}.csv"
            # Already sorted by n_rad, ascending, by report_tree().
            csv_entries = all_entries

        write_geom_models_csv(csv_entries, csv_path)

    print("Done. The geometry files themselves were only read, never written to.")


if __name__ == "__main__":
    main()
