# -*- coding: utf-8 -*-
# =====================================================================
#  Import MY OWN TreeQSM results (the MATLAB "volumes_<tree>_<run>.csv"
#  tables) into the shared master results table volume_results.csv.
# ---------------------------------------------------------------------
#  The MATLAB script run_treeqsm.m writes one table per run, with columns:
#     Tree, Run, Group, Attribute, N, Mean_m3, Std_m3, CV_pct
#  where Group is one of: All inputs / Optimal / Optimal (single) /
#  Estimated / Simplified, and Attribute is Total / Stem / Branches.
#
#  This script scans a folder for those files, picks ONE group per file
#  (by default "Estimated", which combines both runs and has the best
#  standard deviation), and writes one row per file into the shared
#  master table used by compare_volumes.py:
#     tree, method, total_m3, stem_m3, branch_m3, std_m3
#
#  It handles MANY TREES and MANY RUNS at once: every matching file
#  becomes its own row, labelled with its run tag, e.g.
#     IND01_054 | TreeQSM mine (v2, Estimated) | 3.97 | 1.57 | 2.40 | 0.19
#
#  Dependencies: none (standard library only).
# =====================================================================

import csv
import glob
import os

# =====================  PARAMETERS  ==================================
# Folder to scan for the MATLAB volume tables.
MATLAB_RESULTS_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\matlab"

# Which files to pick up (shell-style wildcard).
FILE_PATTERN = "volumes_*.csv"

# Which Group to import from each table. "Estimated" uses the optimal
# models plus the second run, so it has the most reliable mean and std.
# Other valid choices: "Optimal", "All inputs", "Optimal (single)", "Simplified", "Filtered <10cm".
IMPORT_GROUP = "Filtered <10cm"

# Shared master results table (read by compare_volumes.py).
RESULTS_CSV = "volume_results.csv"

# The MATLAB script may use a SHORT tree id (e.g. "IND01") while the Python
# scripts use the FULL id (e.g. "IND01_054"). Rows are only comparable when
# the ids match exactly, so map the short names to the full ones here.
# Any tree not listed is used unchanged.
TREE_ID_MAP = {
    "IND01": "IND01_054",
    # "IND07": "IND07_083",
}

# Reference heights [m] used for DBH (lower) and the taper metric (lower/
# upper). DBH is the stem diameter at TAPER_H_LOWER (1.3 m = breast height).
# Not used directly here (MATLAB doesn't export a diameter profile), but kept
# for consistency with the other scripts and to label what the optional
# height_<tree>_<run>.txt / dbh_<tree>_<run>.txt files are measured at.
TAPER_H_LOWER = 1.3    # lower reference height [m]
TAPER_H_UPPER = 10.0   # upper reference height [m]
# =====================================================================


def to_float(text):
    """Parse a number, or return None for blank/NaN/non-numeric cells."""
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    # MATLAB writes NaN for an undefined std (single-model groups)
    return None if value != value else value


def read_matlab_table(path):
    """Read one volumes_*.csv and return its rows as dicts."""
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def extract_group(rows, group):
    """From one MATLAB table, pull the Total/Stem/Branches values of one Group.
    Returns (tree, run, total, stem, branch, std_of_total) or None if missing."""
    picked = [r for r in rows if r["Group"].strip() == group]
    if not picked:
        return None

    tree = picked[0]["Tree"].strip()
    run = picked[0]["Run"].strip()

    # Build a lookup: Attribute -> (mean, std)
    by_attr = {}
    for r in picked:
        attr = r["Attribute"].strip()
        by_attr[attr] = (to_float(r["Mean_m3"]), to_float(r["Std_m3"]))

    total, std_total = by_attr.get("Total", (None, None))
    stem = by_attr.get("Stem", (None, None))[0]
    branch = by_attr.get("Branches", (None, None))[0]
    return tree, run, total, stem, branch, std_total


def upsert_result(csv_path, tree, method, total, stem, branch, std, dbh=None, height=None, taper=None,
                   trunk_len=None, branch_len=None):
    """Insert/update one (tree, method) row in the shared master results CSV.
    Re-running overwrites the previous row instead of duplicating it.
    Backward compatible: if csv_path still has an older/shorter header, its
    rows are read fine and rewritten with the new columns added as blank."""
    header = ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3",
              "dbh_m", "height_m", "taper_cm_per_m", "trunk_len_m", "branch_len_m"]

    def fmt(x):
        return "" if x is None else "%.6f" % x

    rows = []
    if os.path.exists(csv_path):
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    rows = [r for r in rows if not (r["tree"] == tree and r["method"] == method)]
    rows.append({"tree": tree, "method": method, "total_m3": fmt(total),
                 "stem_m3": fmt(stem), "branch_m3": fmt(branch), "std_m3": fmt(std),
                 "dbh_m": fmt(dbh), "height_m": fmt(height), "taper_cm_per_m": fmt(taper),
                 "trunk_len_m": fmt(trunk_len), "branch_len_m": fmt(branch_len)})

    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


def read_single_number(path):
    """Read a single float from a text file (e.g. 'dbh_<tree>_<run>.txt'),
    or return None if the file does not exist / is not a single number."""
    if not os.path.exists(path):
        return None
    with open(path, "r", encoding="utf-8") as f:
        text = f.read().strip()
    try:
        return float(text)
    except ValueError:
        return None


# =========================  RUN  =====================================
pattern = os.path.join(MATLAB_RESULTS_DIR, FILE_PATTERN)
files = sorted(glob.glob(pattern))

if not files:
    print("No files matching '%s' found in '%s'." % (FILE_PATTERN, MATLAB_RESULTS_DIR))
    raise SystemExit

print("Importing group '%s' from %d file(s):\n" % (IMPORT_GROUP, len(files)))
print("%-28s %-12s %10s %10s %10s %10s %8s %8s" %
      ("file", "tree", "total", "stem", "branch", "std", "dbh", "height"))
print("-" * 102)

imported = 0
for path in files:
    rows = read_matlab_table(path)
    found = extract_group(rows, IMPORT_GROUP)
    if found is None:
        print("%-28s  (no group '%s' in this file - skipped)"
              % (os.path.basename(path), IMPORT_GROUP))
        continue

    tree_raw, run, total, stem, branch, std = found
    tree = TREE_ID_MAP.get(tree_raw, tree_raw)      # map short id -> full id
    method = "TreeQSM mine (%s, %s)" % (run, IMPORT_GROUP)

    # Optional DBH/height/trunk-length/branch-length, read from
    # "dbh_<tree>_<run>.txt" / "height_<tree>_<run>.txt" /
    # "trunklen_<tree>_<run>.txt" / "branchlen_<tree>_<run>.txt" next to the
    # volumes table (single number each: metres for all four). The MATLAB
    # volumes_*.csv tables don't export a diameter profile, so taper cannot be
    # derived here and is always left as None.
    table_dir = os.path.dirname(path)
    dbh = read_single_number(os.path.join(table_dir, "dbh_%s_%s.txt" % (tree, run)))
    height = read_single_number(os.path.join(table_dir, "height_%s_%s.txt" % (tree, run)))
    taper = None
    trunk_len = read_single_number(os.path.join(table_dir, "trunklen_%s_%s.txt" % (tree, run)))
    branch_len = read_single_number(os.path.join(table_dir, "branchlen_%s_%s.txt" % (tree, run)))

    upsert_result(RESULTS_CSV, tree, method, total, stem, branch, std, dbh, height, taper,
                  trunk_len, branch_len)
    imported += 1

    def show(x):
        return "%10.4f" % x if x is not None else "         -"

    def show8(x):
        return "%8.3f" % x if x is not None else "       -"

    print("%-28s %-12s %s %s %s %s %s %s"
          % (os.path.basename(path)[:28], tree, show(total), show(stem),
             show(branch), show(std), show8(dbh), show8(height)))

print("-" * 102)
print("Imported %d row(s) into %s" % (imported, RESULTS_CSV))
print("Run compare_volumes.py to see them next to the other methods.")
