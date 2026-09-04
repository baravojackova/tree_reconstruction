# -*- coding: utf-8 -*-
# =====================================================================
#  Extract the DESTRUCTIVE REFERENCE VOLUME of one tree from the field
#  measurement table (e.g. IND_h_trees.txt from the de Tanago dataset).
# ---------------------------------------------------------------------
#  Each row in the source file is one measured tree SECTION. The columns
#  we care about are:
#     treeID    -> which tree the section belongs to (e.g. "IND01_054")
#     Fraction  -> what part of the tree it is: stem / branch / twig /
#                  leaf / stump / root / fruit
#     Volume    -> the section volume, in CUBIC CENTIMETRES (cm^3)
#
#  The script:
#     1) reads the file,
#     2) keeps only the rows of the tree you choose,
#     3) sums the Volume of the fractions you switch ON,
#     4) converts the sum to cubic metres (m^3) and prints a breakdown.
#
#  Dependencies: none (uses only Python's standard library).
# =====================================================================

import csv
import os

# =====================  PARAMETERS  ==================================
# Tree ID for THIS run. This is the ONLY thing you need to change to switch
# trees - it must match a value in the "treeID" column of the source file
# exactly (if it doesn't, the script prints the full list of available tree
# IDs so you can copy the right one), AND it builds DATA_DIR right below it
# automatically. NOTE: this only works when the data/<tree> folder is named
# EXACTLY like TREE_NAME (true for IND07_083, but NOT for IND01_054 - its
# folder is "IND01_54", a shorter form - if you switch back to that tree,
# set DATA_DIR directly instead of deriving it from TREE_NAME).
TREE_NAME = "IND01_054"  # e.g. "IND01_054" or "IND07_083" (must match the source file's treeID)

# Base folder holding every tree's data, one subfolder per tree (normally
# named exactly like TREE_NAME - see the note above for the one exception).
DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"

# Folder that contains the source measurement file.
DATA_DIR = os.path.join(DATA_ROOT, TREE_NAME)

# Name of the source measurement file (TAB-separated text), inside DATA_DIR.
SOURCE_FILE = os.path.join(DATA_DIR, "IND.h.trees.txt")

# Label used for this method in the shared results table (RESULTS_CSV).
METHOD_LABEL = "Reference (destructive)"

# Shared master results table (see compare_volumes.py). Each script upserts
# its own row(s) into this CSV so results from all methods live in one place.
RESULTS_CSV = "volume_results.csv"

# Which fractions to INCLUDE in the total. Set True to count that fraction,
# False to leave it out. "stem" + "branch" (+ "twig") is the woody structure
# that QSM methods (TreeQSM/AdQSM/AdTree) reconstruct; "leaf" is foliage (not
# wood), "stump"/"root" are at or below ground, "fruit" is fruits.
INCLUDE_FRACTIONS = {
    "stem":   True,     # main trunk sections
    "branch": True,     # branches
    "twig":   False,    # finest branches
    "stump":  False,    # stump (tree base)
    "leaf":   False,    # foliage (not wood)
    "root":   False,    # roots (below ground)
    "fruit":  False,    # fruits
}

# Unit of the "Volume" column in the source file. The de Tanago tables store
# section volume in cm^3. 1 m^3 = 1_000_000 cm^3, so the factor to m^3 is 1e-6.
# (If a file ever stored litres, use 1e-3; if already m^3, use 1.0.)
VOLUME_TO_M3 = 1e-6

# Unit of the "Db"/"De" columns (section-base/section-end diameter) in the
# source file. The de Tanago convention is centimetres; the raw Db of the
# first stem section is printed at runtime so you can verify this against
# the field notes. Change to 1e-3 if your file uses millimetres instead.
DIAMETER_TO_M = 1e-2

# Unit of the "L" column (section length). Same convention as the other
# linear measurements in the table (normally also centimetres).
LENGTH_TO_M = 1e-2

# Reference heights [m] used for DBH (lower) and the taper metric (lower/
# upper). DBH is the stem diameter at TAPER_H_LOWER (1.3 m = breast height).
TAPER_H_LOWER = 1.3    # lower reference height [m]
TAPER_H_UPPER = 10.0   # upper reference height [m]
# =====================================================================


def load_rows(path):
    """Read the TAB-separated file into a list of dict rows (one per section)."""
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def to_float(text):
    """Convert a cell to float, or return None for 'NA'/blank/non-numeric cells."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def stem_diameter_at_height(stem_sections, h):
    """Stem diameter [m] at height h [m] above the tree base, from the
    'stem' fraction's sections (in file order). Walks the sections,
    accumulating height as the running sum of L, and linearly interpolates
    between Db and De inside the section that spans h. Returns None if h
    falls outside the height actually covered by the stem sections (no
    extrapolation/guessing). Reused for both DBH and the taper metric.
    """
    running_h = 0.0
    for r in stem_sections:
        db, de, length = to_float(r["Db"]), to_float(r["De"]), to_float(r["L"])
        if db is None or de is None or length is None:
            continue
        length_m = length * LENGTH_TO_M
        h_start, h_end = running_h, running_h + length_m
        if h_start <= h <= h_end and length_m > 0:
            t = (h - h_start) / length_m
            diameter_m = (db + t * (de - db)) * DIAMETER_TO_M
            return diameter_m
        running_h = h_end
    return None


def upsert_result(csv_path, tree, method, total, trunk, branch, std, dbh=None, height=None, taper=None,
                   trunk_len=None, branch_len=None, branch_filter="none", n_cylinders=None,
                   mode=None, pd1=None, pd2min=None, pd2max=None, mincylrad=None,
                   simp_maxorder=None, simp_smallradii=None, simp_replaceiterations=None):
    """Insert/update one (tree, method) row in the shared master results CSV
    (see compare_volumes.py for its format). Reads csv_path if it exists
    (creating it with the header if not), removes any existing row with the
    same tree AND method, appends the new row, and writes the file back -
    so re-running a script overwrites its previous result instead of
    duplicating it. Numbers are formatted with 6 decimals; None -> "" (blank).
    Backward compatible: if csv_path still has an older/shorter header, its
    rows are read fine and rewritten with the new columns added as blank.

    branch_filter is a plain string, NOT a number, so it is written as-is
    (no fmt()): "none" = full/unfiltered reconstruction (the default), "10cm"
    = trunk/branches restricted to diameter >= 10 cm (matching how the
    destructive reference was physically measured - see compare_volumes.py's
    header comment for why this distinction matters).

    n_cylinders is left as None (-> blank cell) by every call in THIS file:
    the destructive field reference has no reconstructed cylinder model at
    all, so there is no cylinder count to report - blank is the correct
    value here, not 0 (0 would wrongly imply "a model with zero cylinders").

    mode/pd1/pd2min/pd2max/mincylrad/simp_maxorder/simp_smallradii/
    simp_replaceiterations: the ACTUAL TreeQSM reconstruction parameters
    used for a run (see runsken.m section 19's params_<tree>_<run>.csv) -
    all optional/None by default. The destructive field reference has no
    reconstruction run at all, so the one call in this file simply doesn't
    pass them and these columns stay blank for its row - purely additive,
    no behaviour change. mode is a plain string, written as-is like
    branch_filter (blank string, not the literal text "None", when not
    given)."""
    # n_cylinders is the LAST column - kept in sync with the header used by
    # every other upsert_result() copy (tree_geom_utils.py,
    # import_matlab_results.py, qsm_volume_mean.py) so they all write/read
    # the SAME shared volume_results.csv without column mismatches. mode/
    # pd1/.../simp_replaceiterations appended after n_cylinders, same
    # reason - kept in sync with those same three copies.
    header = ["tree", "method", "total_m3", "trunk_m3", "branch_m3", "std_m3",
              "dbh_m", "height_m", "taper_cm_per_m", "trunk_len_m", "branch_len_m",
              "branch_filter", "n_cylinders",
              "mode", "pd1_m", "pd2min_m", "pd2max_m", "mincylrad_m",
              "simp_maxorder", "simp_smallradii", "simp_replaceiterations",
              "adqsm_variant", "radius_threshold_mm", "seg_min_mm", "seg_max_mm", "seg_k_pct"]

    def fmt(x):
        return "" if x is None else "%.6f" % x

    def fmt_int(x):
        # Cylinder count is a whole number, not a measured float, so use
        # "%d" here instead of fmt()'s "%.6f" - still blank ("") when
        # n_cylinders is None, same convention as every other optional column.
        return "" if x is None else "%d" % int(x)

    rows = []
    if os.path.exists(csv_path):
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    rows = [r for r in rows if not (r["tree"] == tree and r["method"] == method)]
    rows.append({"tree": tree, "method": method, "total_m3": fmt(total),
                 "trunk_m3": fmt(trunk), "branch_m3": fmt(branch), "std_m3": fmt(std),
                 "dbh_m": fmt(dbh), "height_m": fmt(height), "taper_cm_per_m": fmt(taper),
                 "trunk_len_m": fmt(trunk_len), "branch_len_m": fmt(branch_len),
                 "branch_filter": branch_filter, "n_cylinders": fmt_int(n_cylinders),
                 "mode": mode or "", "pd1_m": fmt(pd1), "pd2min_m": fmt(pd2min),
                 "pd2max_m": fmt(pd2max), "mincylrad_m": fmt(mincylrad),
                 "simp_maxorder": fmt(simp_maxorder), "simp_smallradii": fmt(simp_smallradii),
                 "simp_replaceiterations": fmt(simp_replaceiterations),
                 "adqsm_variant": "", "radius_threshold_mm": "", "seg_min_mm": "",
                 "seg_max_mm": "", "seg_k_pct": ""})

    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


# =========================  RUN  =====================================
rows = load_rows(SOURCE_FILE)

# List of every tree present in the file (sorted), for validation/help.
all_trees = sorted({r["treeID"] for r in rows})

if TREE_NAME not in all_trees:
    print("Tree '%s' was not found in %s." % (TREE_NAME, SOURCE_FILE))
    print("Available tree IDs (%d):" % len(all_trees))
    for t in all_trees:
        print("   ", t)
    raise SystemExit  # stop cleanly so you can fix TREE_NAME and re-run

# Keep only this tree's sections.
tree_rows = [r for r in rows if r["treeID"] == TREE_NAME]

# Sum the section volumes per fraction (in the source unit, cm^3), and count
# how many sections each fraction has. Also sum the section LENGTHS ("L",
# same column already read by stem_diameter_at_height above) per fraction,
# the exact same way - this gives trunk_len_m/branch_len_m for the CSV below.
volume_by_fraction = {}
length_by_fraction = {}
count_by_fraction = {}
for r in tree_rows:
    frac = r["Fraction"]
    vol = to_float(r["Volume"])
    if vol is None:
        continue  # skip sections without a numeric volume (e.g. 'NA')
    volume_by_fraction[frac] = volume_by_fraction.get(frac, 0.0) + vol
    count_by_fraction[frac] = count_by_fraction.get(frac, 0) + 1
    length = to_float(r["L"])
    if length is not None:
        length_by_fraction[frac] = length_by_fraction.get(frac, 0.0) + length

# ---- print a per-fraction breakdown, marking which ones are included -------
print("Tree: %s   (source: %s)" % (TREE_NAME, SOURCE_FILE))
print()
print("%-8s %8s %14s   %s" % ("Fraction", "sections", "volume [m^3]", "included?"))
print("-" * 48)

selected_total_m3 = 0.0
# Show every fraction actually present for this tree, in a stable order.
fraction_order = ["stem", "branch", "twig", "stump", "leaf", "root", "fruit"]
present = [f for f in fraction_order if f in volume_by_fraction]
# also append any unexpected fraction names that are not in our known list
present += [f for f in volume_by_fraction if f not in fraction_order]

for frac in present:
    vol_m3 = volume_by_fraction[frac] * VOLUME_TO_M3
    include = INCLUDE_FRACTIONS.get(frac, False)  # unknown fractions default to OFF
    if include:
        selected_total_m3 += vol_m3
    print("%-8s %8d %14.4f   %s"
          % (frac, count_by_fraction[frac], vol_m3, "YES" if include else "no"))

print("-" * 48)
print("SELECTED TOTAL VOLUME: %.4f m^3  (%.1f L)"
      % (selected_total_m3, selected_total_m3 * 1000.0))

# For convenience, also show the full-tree total (all fractions, whatever is
# measured), so you can see how much the switched-off parts would add.
full_total_m3 = sum(volume_by_fraction.values()) * VOLUME_TO_M3
print("(all fractions combined would be: %.4f m^3)" % full_total_m3)

# ---- DBH / height / taper from the 'stem' fraction's sections --------------
stem_sections = [r for r in tree_rows if r["Fraction"] == "stem"]
if stem_sections:
    print("\nRaw Db of the first stem section: %r  (verify unit; DIAMETER_TO_M = %g)"
          % (stem_sections[0]["Db"], DIAMETER_TO_M))

d_lower = stem_diameter_at_height(stem_sections, TAPER_H_LOWER) if stem_sections else None
d_upper = stem_diameter_at_height(stem_sections, TAPER_H_UPPER) if stem_sections else None
dbh_m = d_lower
# The source table has no independent total-height measurement; do not guess it.
height_m = None
taper_cm_per_m = None
if d_lower is not None and d_upper is not None:
    taper_cm_per_m = (d_lower - d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)

print("DBH (stem diameter at %.1f m)   : %s"
      % (TAPER_H_LOWER, ("%.4f m" % dbh_m) if dbh_m is not None else "not resolved"))
print("Taper (%.1f-%.1f m)             : %s"
      % (TAPER_H_LOWER, TAPER_H_UPPER,
         ("%.2f cm/m" % taper_cm_per_m) if taper_cm_per_m is not None else "not resolved"))

# ---- upsert this result into the shared master results CSV -----------------
# "stem" here is the SOURCE FILE's own "Fraction" column value (and
# INCLUDE_FRACTIONS' matching key) - external data schema, not renamed. Our
# OWN output variable (what gets written to the trunk_m3 column) is named
# trunk_included/trunk_m3, consistent with the rest of the shared CSV.
trunk_included = INCLUDE_FRACTIONS.get("stem", False) and "stem" in volume_by_fraction
trunk_m3 = (volume_by_fraction["stem"] * VOLUME_TO_M3) if trunk_included else None
branch_m3 = (selected_total_m3 - trunk_m3) if trunk_m3 is not None else None

# trunk_len_m/branch_len_m: same "trunk included?" / "everything else selected"
# logic as trunk_m3/branch_m3 above, just applied to length_by_fraction instead
# of volume_by_fraction.
trunk_len_m = (length_by_fraction["stem"] * LENGTH_TO_M) if (trunk_included and "stem" in length_by_fraction) else None
selected_total_len_m = sum(length_by_fraction[f] for f in length_by_fraction
                            if INCLUDE_FRACTIONS.get(f, False)) * LENGTH_TO_M
branch_len_m = (selected_total_len_m - trunk_len_m) if trunk_len_m is not None else None

# branch_filter = "10cm": the de Tanago field crew physically only measured
# sections down to a 10 cm taper diameter (see AdQSM.pdf Appendix A) - this
# reference NEVER has a "none" (full/unfiltered) version, by methodology.
upsert_result(RESULTS_CSV, TREE_NAME, METHOD_LABEL, selected_total_m3, trunk_m3, branch_m3, None,
              dbh_m, height_m, taper_cm_per_m, trunk_len_m, branch_len_m, branch_filter="10cm")
