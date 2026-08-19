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
DATA_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data\IND01_54\qsm_tanago"

# Tree to report. This is the tree ID embedded in the file names, e.g. for
# "cyl_data_IND01_054.txt_0.5_0.55_5_0.025_0.075_3_4_1_t0.txt" ... "..._t19.txt"
# use "IND01_054". To report a different tree, just change this to its ID.
TREE_PREFIX = "IND01_054"

# Column indices inside each cylinder file (0-based), for the TreeQSM
# cylinder-data format. Change these only if your files use a different layout.
COL_RADIUS = 0
COL_LENGTH = 1
COL_START_Z = 4
COL_ORDER = 11

# Reference heights [m] used for DBH (lower) and the taper metric (lower/
# upper). DBH is the stem diameter at TAPER_H_LOWER (1.3 m = breast height).
TAPER_H_LOWER = 1.3    # lower reference height [m]
TAPER_H_UPPER = 10.0   # upper reference height [m]

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


def stem_diameter_at_height(order0_rows, base_z, h):
    """Diameter [m] of the stem (BranchOrder 0 cylinders) at height h [m]
    above the tree base (base_z). Each cylinder's approximate height span is
    [start_z - base_z, start_z - base_z + length]; the diameter of the span
    that contains h is returned, falling back to the nearest cylinder for
    small gaps between spans. Returns None if h is outside the height range
    actually covered by stem cylinders (no extrapolation/guessing). Reused
    for both DBH and the taper metric."""
    if len(order0_rows) == 0:
        return None
    spans = [(sz - base_z, sz - base_z + length, radius) for sz, length, radius in order0_rows]
    lo = min(s[0] for s in spans)
    hi = max(s[1] for s in spans)
    if h < lo or h > hi:
        return None
    for h_start, h_end, radius in spans:
        if h_start <= h <= h_end:
            return 2.0 * radius

    def dist(span):
        h_start, h_end, _ = span
        return h_start - h if h < h_start else h - h_end

    nearest = min(spans, key=dist)
    return 2.0 * nearest[2]


def cylinder_metrics(path):
    """Load one cylinder file and return a dict with total/stem/branch volume
    [m^3], DBH [m], tree height [m], and taper [cm/m] (dbh/taper may be None
    if they cannot be determined for this realization)."""
    d = np.loadtxt(path)
    if d.ndim == 1:            # guard: a file with a single cylinder row
        d = d.reshape(1, -1)
    r = d[:, COL_RADIUS]
    length = d[:, COL_LENGTH]
    start_z = d[:, COL_START_Z]
    order = d[:, COL_ORDER].astype(int)

    vol = np.pi * r ** 2 * length          # volume of each cylinder [m^3]
    total = float(vol.sum())
    stem = float(vol[order == 0].sum())    # BranchOrder 0 = trunk/stem
    branch = float(vol[order >= 1].sum())

    # trunk_len/branch_len: sum the SAME `length` array already loaded above,
    # split by BranchOrder the same way stem/branch volume are (just sum
    # length instead of pi*r^2*length).
    trunk_len = float(length[order == 0].sum())
    branch_len = float(length[order >= 1].sum())

    base_z = float(start_z.min())
    order0_rows = list(zip(start_z[order == 0], length[order == 0], r[order == 0]))
    dbh = stem_diameter_at_height(order0_rows, base_z, TAPER_H_LOWER)
    d_upper = stem_diameter_at_height(order0_rows, base_z, TAPER_H_UPPER)
    taper = ((dbh - d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
             if dbh is not None and d_upper is not None else None)

    top_idx = int(np.argmax(start_z))      # topmost cylinder by start_z
    height = float(start_z[top_idx] + length[top_idx] - base_z)

    return dict(total=total, stem=stem, branch=branch, dbh=dbh, height=height, taper=taper,
                trunk_len=trunk_len, branch_len=branch_len)


def mean_ignore_none(values):
    """Mean of the non-None values, or None if there are none."""
    vals = [v for v in values if v is not None]
    return (sum(vals) / len(vals)) if vals else None


def summarize(values):
    """Return (mean, std, cv_percent) for a list of numbers. The standard
    deviation uses ddof=1 (sample standard deviation, as for repeated
    measurements). CV = std / mean * 100."""
    arr = np.asarray(values, dtype=float)
    mean = float(arr.mean())
    std = float(arr.std(ddof=1)) if len(arr) > 1 else 0.0
    cv = (std / mean * 100.0) if mean != 0 else 0.0
    return mean, std, cv


def upsert_result(csv_path, tree, method, total, stem, branch, std, dbh=None, height=None, taper=None,
                   trunk_len=None, branch_len=None):
    """Insert/update one (tree, method) row in the shared master results CSV
    (see compare_volumes.py for its format). Reads csv_path if it exists
    (creating it with the header if not), removes any existing row with the
    same tree AND method, appends the new row, and writes the file back -
    so re-running a script overwrites its previous result instead of
    duplicating it. Numbers are formatted with 6 decimals; None -> "" (blank).
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

# Compute the three volumes (+ DBH/height/taper/trunk_len/branch_len) for every realization.
totals, stems, branches = [], [], []
dbhs, heights, tapers = [], [], []
trunk_lens, branch_lens = [], []
print("Tree: %s   (%d realization files)\n" % (TREE_PREFIX, len(files)))
print("%-24s %12s %12s %12s %8s %8s" % ("file", "total m^3", "stem m^3", "branch m^3", "dbh cm", "h m"))
print("-" * 78)
for idx, path in files:
    m = cylinder_metrics(path)
    totals.append(m["total"]); stems.append(m["stem"]); branches.append(m["branch"])
    dbhs.append(m["dbh"]); heights.append(m["height"]); tapers.append(m["taper"])
    trunk_lens.append(m["trunk_len"]); branch_lens.append(m["branch_len"])
    print("%-24s %12.4f %12.4f %12.4f %8s %8.2f"
          % (os.path.basename(path), m["total"], m["stem"], m["branch"],
             ("%.2f" % (m["dbh"] * 100.0)) if m["dbh"] is not None else "-", m["height"]))

# Summary statistics across the realizations.
print("-" * 78)
for label, vals in (("TOTAL", totals), ("STEM", stems), ("BRANCHES", branches)):
    mean, std, cv = summarize(vals)
    print("%-9s mean = %.4f m^3   std = %.4f m^3   CV = %.2f %%   (min %.4f, max %.4f)"
          % (label, mean, std, cv, min(vals), max(vals)))

mean_total = summarize(totals)[0]
print("\n=> The published QSM volume for this tree is the TOTAL mean above: "
      "%.4f m^3 (%.1f L)." % (mean_total, mean_total * 1000.0))

mean_dbh = mean_ignore_none(dbhs)
mean_height = mean_ignore_none(heights)
mean_taper = mean_ignore_none(tapers)
print("=> Mean DBH = %s   Mean height = %s   Mean taper = %s"
      % (("%.4f m" % mean_dbh) if mean_dbh is not None else "n/a",
         ("%.2f m" % mean_height) if mean_height is not None else "n/a",
         ("%.2f cm/m" % mean_taper) if mean_taper is not None else "n/a"))

# ---- upsert this result into the shared master results CSV -----------------
std_total = summarize(totals)[1]
mean_stem = summarize(stems)[0]
mean_branch = summarize(branches)[0]
# trunk_len/branch_len: mean across realizations, same averaging (mean_ignore_none)
# already used for dbh/height/taper above - every realization always HAS a
# trunk_len/branch_len (unlike dbh, which can be None), but mean_ignore_none
# works fine either way.
mean_trunk_len = mean_ignore_none(trunk_lens)
mean_branch_len = mean_ignore_none(branch_lens)
upsert_result(RESULTS_CSV, TREE_PREFIX, METHOD_LABEL,
              mean_total, mean_stem, mean_branch, std_total,
              mean_dbh, mean_height, mean_taper,
              mean_trunk_len, mean_branch_len)

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
