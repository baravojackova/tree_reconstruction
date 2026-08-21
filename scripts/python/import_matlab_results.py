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
IMPORT_GROUP = "Optimal"  # the de Tanago field crew physically only measured sections down to a 10 cm taper diameter (see AdQSM.pdf Appendix A) - this reference NEVER has a "none" (full/unfiltered) version, by methodology.

# Shared master results table (read by compare_volumes.py).
RESULTS_CSV = "volume_results.csv"

# The MATLAB script may use a SHORT tree id (e.g. "IND01") while the Python
# scripts use the FULL id (e.g. "IND01_054"). Rows are only comparable when
# the ids match exactly, so map the short names to the full ones here.
# Any tree not listed is used unchanged.
TREE_NAME_MAP = {
    "IND01": "IND01_054",
    "IND07": "IND07_083",
}

# Reference heights [m] used for DBH (lower) and the taper metric (lower/
# upper). DBH is the stem diameter at TAPER_H_LOWER (1.3 m = breast height).
# These two constants are NOT used directly in this script's own calculation
# (taper is computed by runsken.m in MATLAB, using these exact same two
# heights, and just read in below) - they're kept here to document what the
# height_<tree>_<run>.txt / dbh_<tree>_<run>.txt / taper_<tree>_<run>.txt
# files are measured at, and so this file's numbers stay comparable with
# every OTHER method in volume_results.csv, which all use these same two
# reference heights for their own taper_cm_per_m (see
# adtree_reconstruct_compare.py's raw_taper/cal_taper for the other side of
# that same convention).
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
                   trunk_len=None, branch_len=None, branch_filter="none"):
    """Insert/update one (tree, method) row in the shared master results CSV.
    Re-running overwrites the previous row instead of duplicating it.
    Backward compatible: if csv_path still has an older/shorter header, its
    rows are read fine and rewritten with the new columns added as blank.

    branch_filter is a plain string, NOT a number, so it is written as-is
    (no fmt()): "none" = full/unfiltered reconstruction (the default), "10cm"
    = trunk/branches restricted to diameter >= 10 cm (matching how the
    destructive reference was physically measured - see compare_volumes.py's
    header comment for why this distinction matters)."""
    header = ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3",
              "dbh_m", "height_m", "taper_cm_per_m", "trunk_len_m", "branch_len_m",
              "branch_filter"]

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
                 "trunk_len_m": fmt(trunk_len), "branch_len_m": fmt(branch_len),
                 "branch_filter": branch_filter})

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
    tree = TREE_NAME_MAP.get(tree_raw, tree_raw)      # map short id -> full id
    method = "TreeQSM mine (%s, %s)" % (run, IMPORT_GROUP)

    # Optional DBH/height/taper, read from "dbh_<tree>_<run>.txt" /
    # "height_<tree>_<run>.txt" / "taper_<tree>_<run>.txt" next to the
    # volumes table (single number each: metres for dbh/height, cm/m for
    # taper). All three are exported by runsken.m section 18 using the SAME
    # TAPER_H_LOWER/TAPER_H_UPPER reference heights defined above, which is
    # what makes taper here comparable to taper_cm_per_m for every other
    # method in volume_results.csv.
    #
    # taper has no "_filtered" counterpart (unlike trunk/branch LENGTH just
    # below): it's the diameter narrowing of the STEM between 1.3 m and
    # 10.0 m, and the stem at breast height is essentially always thicker
    # than the 10 cm cut-off anyway - same reason DBH/height don't get a
    # filtered variant either. So taper is read the same way regardless of
    # branch_filter (computed further down).
    table_dir = os.path.dirname(path)
    dbh = read_single_number(os.path.join(table_dir, "dbh_%s_%s.txt" % (tree, run)))
    height = read_single_number(os.path.join(table_dir, "height_%s_%s.txt" % (tree, run)))
    taper = read_single_number(os.path.join(table_dir, "taper_%s_%s.txt" % (tree, run)))

    # branch_filter: "10cm" for any group whose NAME says it's restricted to
    # branches >= 10 cm diameter (currently only IMPORT_GROUP = "Filtered
    # <10cm", from runsken.m's section 17b) - "none" for every other group
    # (Estimated/Optimal/Optimal (single)/All inputs/Simplified), which use
    # the full, unfiltered TreeQSM cylinder model. Computed BEFORE the
    # trunk/branch length read below, so that read can reuse this exact same
    # test to pick the right length file - see comment there.
    branch_filter = "10cm" if "Filtered" in IMPORT_GROUP else "none"

    # Trunk/branch length: runsken.m section 18 exports TWO separate pairs of
    # files - "trunklen_/branchlen_<tree>_<run>.txt" hold the UNFILTERED
    # model's length (correct for Optimal/Estimated/... groups), while
    # "trunklen_filtered_/branchlen_filtered_<tree>_<run>.txt" hold the
    # length AFTER the same 10 cm cut-off as the "Filtered <10cm" group's
    # volume. Reading the unfiltered file for a filtered group would silently
    # pair a filtered VOLUME with an unfiltered LENGTH - so which file to read
    # is decided by the SAME branch_filter test used just above, not a second,
    # differently-worded check.
    if branch_filter == "10cm":
        trunk_len_file = "trunklen_filtered_%s_%s.txt" % (tree, run)
        branch_len_file = "branchlen_filtered_%s_%s.txt" % (tree, run)
    else:
        trunk_len_file = "trunklen_%s_%s.txt" % (tree, run)
        branch_len_file = "branchlen_%s_%s.txt" % (tree, run)

    trunk_len = read_single_number(os.path.join(table_dir, trunk_len_file))
    branch_len = read_single_number(os.path.join(table_dir, branch_len_file))
    if branch_filter == "10cm" and (trunk_len is None or branch_len is None):
        # Don't crash - just flag it clearly so it's not silently wrong: an
        # older MATLAB run (before this fix) won't have written these two
        # files yet, and re-using the unfiltered length would be the exact
        # bug this change fixes. Re-run runsken.m for this tree/run to get them.
        print("  WARNING: missing filtered length file(s) for %s / %s (%s / %s) - "
              "trunk_len/branch_len left blank for this row. Re-run runsken.m "
              "for this tree/run to generate them."
              % (tree, run, trunk_len_file, branch_len_file))

    upsert_result(RESULTS_CSV, tree, method, total, stem, branch, std, dbh, height, taper,
                  trunk_len, branch_len, branch_filter=branch_filter)
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
