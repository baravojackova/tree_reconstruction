# -*- coding: utf-8 -*-
# =====================================================================
#  Average QSM volume of one tree from its 20 TreeQSM realization files.
# ---------------------------------------------------------------------
#  TreeQSM builds its cover sets from RANDOMLY chosen points, so every
#  reconstruction of the same tree with the same parameters gives a
#  slightly different model. The authors therefore reconstruct each tree
#  many times (here 20 files, one per realization) and report the MEAN
#  volume plus its spread (standard deviation).
#
#  Each file is a TreeQSM cylinder model: one row per cylinder. The
#  columns we use (0-based) are:
#     column 0  -> radius       [m]
#     column 1  -> length       [m]
#     column 11 -> BranchOrder  (0 = stem/trunk, >=1 = branches)
#  A cylinder's volume is pi * radius^2 * length, and the tree volume is
#  the sum over all cylinders. (radius & length are in metres, so the
#  volume comes out directly in cubic metres, m^3.)
#
#  The script reads every realization file of the chosen tree, computes
#  total / stem / branch volume for each, and prints the mean, standard
#  deviation and coefficient of variation across the realizations.
#
#  Dependencies: numpy   (install: pip install numpy)
# =====================================================================

import os
import re
import csv
import numpy as np

# =====================  PARAMETERS  ==================================
# Folder that contains the realization files. "." means "the folder this
# script is run from".
DATA_DIR = "C:\\Users\\Spravce\\Documents\\BARA\\01_Skeny_Babice\\skeny_tanago\\AdTree\\qsm_tango"

# Tree to report. This is the tree ID embedded in the file names, e.g. for
# "cyl_data_IND01_054.txt_0.5_0.55_5_0.025_0.075_3_4_1_t0.txt" ... "..._t19.txt"
# use "IND01_054". To report a different tree, just change this to its ID.
TREE_PREFIX = "IND01_054"

# Column indices inside each cylinder file (0-based), for the TreeQSM
# cylinder-data format. Change these only if your files use a different layout.
COL_RADIUS = 0
COL_LENGTH = 1
COL_ORDER = 11

# Where to write the CSV export (per-realization rows + summary rows).
# Set to None to skip CSV export.
CSV_PATH = os.path.join(DATA_DIR, "%s_volume_summary.csv" % TREE_PREFIX)

# Label used for this method in the shared results table (RESULTS_CSV).
METHOD_LABEL = "TreeQSM de Tanago (mean)"

# Shared master results table (see compare_volumes.py). Each script upserts
# its own row(s) into this CSV so results from all methods live in one place.
RESULTS_CSV = "volume_results.csv"
# =====================================================================


def find_realization_files(data_dir, prefix):
    """Return the tree's realization files, sorted by their trailing
    realization number. Matches names like
    'cyl_data_<prefix>.txt_<...params...>_t<n>.txt' (case-insensitive),
    e.g. 'cyl_data_IND01_054.txt_0.5_0.55_5_0.025_0.075_3_4_1_t1.txt'."""
    pattern = re.compile(
        r"^cyl_data_" + re.escape(prefix) + r"\.txt_.*_t(\d+)\.txt$",
        re.IGNORECASE,
    )
    hits = []
    for name in os.listdir(data_dir):
        m = pattern.match(name)
        if m:
            idx = int(m.group(1))                   # realization number (0-19)
            hits.append((idx, os.path.join(data_dir, name)))
    hits.sort(key=lambda t: t[0])                   # sort by realization number
    return hits


def cylinder_volumes(path):
    """Load one cylinder file and return (total, stem, branch) volume in m^3."""
    d = np.loadtxt(path)
    if d.ndim == 1:            # guard: a file with a single cylinder row
        d = d.reshape(1, -1)
    r = d[:, COL_RADIUS]
    length = d[:, COL_LENGTH]
    order = d[:, COL_ORDER].astype(int)
    vol = np.pi * r ** 2 * length          # volume of each cylinder [m^3]
    total = float(vol.sum())
    stem = float(vol[order == 0].sum())    # BranchOrder 0 = trunk/stem
    branch = float(vol[order >= 1].sum())
    return total, stem, branch


def summarize(values):
    """Return (mean, std, cv_percent) for a list of numbers. The standard
    deviation uses ddof=1 (sample standard deviation, as for repeated
    measurements). CV = std / mean * 100."""
    arr = np.asarray(values, dtype=float)
    mean = float(arr.mean())
    std = float(arr.std(ddof=1)) if len(arr) > 1 else 0.0
    cv = (std / mean * 100.0) if mean != 0 else 0.0
    return mean, std, cv


def upsert_result(csv_path, tree, method, total, stem, branch, std):
    """Insert/update one (tree, method) row in the shared master results CSV
    (see compare_volumes.py for its format). Reads csv_path if it exists
    (creating it with the header if not), removes any existing row with the
    same tree AND method, appends the new row, and writes the file back -
    so re-running a script overwrites its previous result instead of
    duplicating it. Numbers are formatted with 6 decimals; None -> "" (blank)."""
    header = ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3"]

    def fmt(x):
        return "" if x is None else "%.6f" % x

    rows = []
    if os.path.exists(csv_path):
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    rows = [r for r in rows if not (r["tree"] == tree and r["method"] == method)]
    rows.append({"tree": tree, "method": method, "total_m3": fmt(total),
                 "stem_m3": fmt(stem), "branch_m3": fmt(branch), "std_m3": fmt(std)})

    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


# =========================  RUN  =====================================
files = find_realization_files(DATA_DIR, TREE_PREFIX)

if not files:
    print("No files found for prefix '%s' in folder '%s'." % (TREE_PREFIX, DATA_DIR))
    # Help the user: show which tree prefixes ARE available in the folder.
    name_pattern = re.compile(r"^cyl_data_(.+)\.txt_.*_t\d+\.txt$", re.IGNORECASE)
    prefixes = sorted({name_pattern.match(f).group(1) for f in os.listdir(DATA_DIR)
                       if name_pattern.match(f)})
    if prefixes:
        print("Available tree prefixes in this folder:")
        for p in prefixes:
            print("   ", p)
    raise SystemExit

# Compute the three volumes for every realization.
totals, stems, branches = [], [], []
print("Tree: %s   (%d realization files)\n" % (TREE_PREFIX, len(files)))
print("%-24s %12s %12s %12s" % ("file", "total m^3", "stem m^3", "branch m^3"))
print("-" * 64)
for idx, path in files:
    total, stem, branch = cylinder_volumes(path)
    totals.append(total); stems.append(stem); branches.append(branch)
    print("%-24s %12.4f %12.4f %12.4f"
          % (os.path.basename(path), total, stem, branch))

# Summary statistics across the realizations.
print("-" * 64)
for label, vals in (("TOTAL", totals), ("STEM", stems), ("BRANCHES", branches)):
    mean, std, cv = summarize(vals)
    print("%-9s mean = %.4f m^3   std = %.4f m^3   CV = %.2f %%   (min %.4f, max %.4f)"
          % (label, mean, std, cv, min(vals), max(vals)))

mean_total = summarize(totals)[0]
print("\n=> The published QSM volume for this tree is the TOTAL mean above: "
      "%.4f m^3 (%.1f L)." % (mean_total, mean_total * 1000.0))

# ---- upsert this result into the shared master results CSV -----------------
std_total = summarize(totals)[1]
mean_stem = summarize(stems)[0]
mean_branch = summarize(branches)[0]
upsert_result(RESULTS_CSV, TREE_PREFIX, METHOD_LABEL,
              mean_total, mean_stem, mean_branch, std_total)

# Export the per-realization results and summary statistics to CSV.
if CSV_PATH:
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["tree", "file", "total_m3", "stem_m3", "branch_m3"])
        for (idx, path), total, stem, branch in zip(files, totals, stems, branches):
            writer.writerow([TREE_PREFIX, os.path.basename(path),
                              "%.6f" % total, "%.6f" % stem, "%.6f" % branch])
        writer.writerow([])
        writer.writerow(["stat", "total_m3", "stem_m3", "branch_m3"])
        for stat_name, stat_idx in (("mean", 0), ("std", 1), ("cv_percent", 2)):
            row = [stat_name]
            for vals in (totals, stems, branches):
                row.append("%.6f" % summarize(vals)[stat_idx])
            writer.writerow(row)
        writer.writerow(["min"] + ["%.6f" % min(vals) for vals in (totals, stems, branches)])
        writer.writerow(["max"] + ["%.6f" % max(vals) for vals in (totals, stems, branches)])
    print("\nCSV exported to: %s" % CSV_PATH)
