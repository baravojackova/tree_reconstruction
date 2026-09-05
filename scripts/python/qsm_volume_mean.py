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
#  total / trunk / branch volume for each, and prints the mean, standard
#  deviation and coefficient of variation across the realizations.
#
#  Dependencies: numpy   (install: pip install numpy)
# =====================================================================

import os
import re
import csv
import numpy as np

# =====================  PARAMETERS  ==================================
# Tree ID for THIS run - this is the ONLY thing you need to change to switch
# trees. It's the tree ID embedded in the realization file names, e.g. for
# "cyl_data_IND07_083.txt_0.5_0.55_5_0.025_0.075_3_4_1_t0.txt" ... "..._t19.txt"
# it's "IND07_083", AND it builds DATA_DIR right below it automatically.
# NOTE: this only works when the data/<tree> folder is named EXACTLY like
# TREE_NAME (true for IND07_083, but NOT for IND01_054 - its folder is
# "IND01_54", a shorter form - if you switch back to that tree, set
# DATA_DIR directly instead of deriving it from TREE_NAME).
TREE_NAME = "IND07_083"

# Base folder holding every tree's data, one subfolder per tree (normally
# named exactly like TREE_NAME - see the note above for the one exception).
DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"

# Folder that contains the realization files.
DATA_DIR = os.path.join(DATA_ROOT, TREE_NAME, "qsm_tanago")

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
CSV_PATH = os.path.join(DATA_DIR, "%s_volume_summary.csv" % TREE_NAME)

# Label used for this method in the shared results table (RESULTS_CSV).
METHOD_LABEL = "TreeQSM de Tanago (mean)"

# --- 10 cm diameter cut-off (second, filtered row) ----------------------
# The destructive field reference only ever measured trunk/branches down to
# a 10 cm taper diameter (de Tanago methodology, see AdQSM.pdf Appendix A -
# same cut-off already applied in runsken.m section 17b and to the
# "TreeQSM mine (*, Filtered<10cm)" rows). This script ALSO upserts a
# second row with exactly that cut-off applied to the published de Tanago
# realizations, so they can be compared to the reference on the same
# methodological footing. THIN_BRANCH_CUT_CM is the diameter threshold [cm];
# METHOD_LABEL_10CM is that second row's method name.
THIN_BRANCH_CUT_CM = 10.0
METHOD_LABEL_10CM = "TreeQSM de Tanago (mean, Filtered<10cm)"

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


def cylinder_metrics(path, cut_cm=THIN_BRANCH_CUT_CM):
    """Load one cylinder file and return a dict with total/trunk/branch volume
    [m^3], DBH [m], tree height [m], and taper [cm/m] (dbh/taper may be None
    if they cannot be determined for this realization) - PLUS the same three
    volumes computed a SECOND time, restricted to cylinders whose diameter is
    >= cut_cm (the "_10cm" keys below), using the exact same pi*r^2*length
    formula on a filtered subset of the SAME arrays, not a separate calculation."""
    d = np.loadtxt(path)
    if d.ndim == 1:            # guard: a file with a single cylinder row
        d = d.reshape(1, -1)
    r = d[:, COL_RADIUS]
    length = d[:, COL_LENGTH]
    start_z = d[:, COL_START_Z]
    order = d[:, COL_ORDER].astype(int)

    vol = np.pi * r ** 2 * length          # volume of each cylinder [m^3]
    total = float(vol.sum())
    trunk = float(vol[order == 0].sum())   # BranchOrder 0 = trunk
    branch = float(vol[order >= 1].sum())

    # --- 10 cm diameter cut-off variant (same principle as runsken.m 17b) ---
    # diameter_cm = 2 * radius[m] * 100; keep = diameter_cm >= cut_cm; then
    # total/trunk/branch volume again, but summed only over the kept cylinders.
    diameter_cm = 2.0 * r * 100.0
    keep = diameter_cm >= cut_cm
    n_total = len(r)
    n_kept = int(keep.sum())
    total_10cm = float(vol[keep].sum())
    trunk_10cm = float(vol[keep & (order == 0)].sum())
    branch_10cm = float(vol[keep & (order >= 1)].sum())

    # trunk_len/branch_len: sum the SAME `length` array already loaded above,
    # split by BranchOrder the same way trunk/branch volume are (just sum
    # length instead of pi*r^2*length).
    trunk_len = float(length[order == 0].sum())
    branch_len = float(length[order >= 1].sum())

    # trunk_len_10cm/branch_len_10cm: the filtered (>= cut_cm diameter)
    # counterpart of trunk_len/branch_len above - SAME `keep` mask already
    # computed for total_10cm/trunk_10cm/branch_10cm, just summing `length`
    # instead of `vol` over it (no new pass over the file, no new mask).
    trunk_len_10cm = float(length[keep & (order == 0)].sum())
    branch_len_10cm = float(length[keep & (order >= 1)].sum())

    base_z = float(start_z.min())
    order0_rows = list(zip(start_z[order == 0], length[order == 0], r[order == 0]))
    dbh = stem_diameter_at_height(order0_rows, base_z, TAPER_H_LOWER)
    d_upper = stem_diameter_at_height(order0_rows, base_z, TAPER_H_UPPER)
    taper = ((dbh - d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
             if dbh is not None and d_upper is not None else None)

    top_idx = int(np.argmax(start_z))      # topmost cylinder by start_z
    height = float(start_z[top_idx] + length[top_idx] - base_z)

    return dict(total=total, trunk=trunk, branch=branch, dbh=dbh, height=height, taper=taper,
                trunk_len=trunk_len, branch_len=branch_len,
                total_10cm=total_10cm, trunk_10cm=trunk_10cm, branch_10cm=branch_10cm,
                trunk_len_10cm=trunk_len_10cm, branch_len_10cm=branch_len_10cm,
                n_total=n_total, n_kept=n_kept)


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

    n_cylinders defaults to None (-> blank cell) so other/future callers in
    this file that don't have a cylinder count stay blank rather than
    wrongly writing 0 - but the two calls actually used below DO pass a
    real value (mean_n_total / mean_n_kept, already computed from each
    realization file's cylinder count).

    mode/pd1/pd2min/pd2max/mincylrad/simp_maxorder/simp_smallradii/
    simp_replaceiterations: the ACTUAL TreeQSM reconstruction parameters
    used for a run (see runsken.m section 19's params_<tree>_<run>.csv) -
    all optional/None by default. This file has no real values for them
    (it summarizes published/precomputed realizations, not a live
    runsken.m run), so both calls below simply don't pass them and these
    columns stay blank for this script's own rows - purely additive, no
    behaviour change. mode is a plain string, written as-is like
    branch_filter (blank string, not the literal text "None", when not
    given)."""
    # n_cylinders is the LAST column - kept in sync with the header used by
    # every other upsert_result() copy (tree_geom_utils.py,
    # import_matlab_results.py, reference_volume.py) so they all write/read
    # the SAME shared volume_results.csv without column mismatches. mode/
    # pd1/.../simp_replaceiterations appended after n_cylinders, same
    # reason - kept in sync with those same three copies.
    header = ["tree", "method", "total_m3", "trunk_m3", "branch_m3", "std_m3",
              "dbh_m", "height_m", "taper_cm_per_m", "trunk_len_m", "branch_len_m",
              "branch_filter", "n_cylinders",
              "mode", "pd1_m", "pd2min_m", "pd2max_m", "mincylrad_m",
              "simp_maxorder", "simp_smallradii", "simp_replaceiterations",
              "adqsm_variant", "radius_threshold_mm", "seg_min_mm", "seg_max_mm", "seg_k_pct",
              "calmethod"]

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
                 "seg_max_mm": "", "seg_k_pct": "", "calmethod": ""})

    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


# =========================  RUN  =====================================
files = find_realization_files(DATA_DIR, TREE_NAME)

if not files:
    print("No files found for prefix '%s' in folder '%s'." % (TREE_NAME, DATA_DIR))
    # Help the user: show which tree prefixes ARE available in the folder.
    name_pattern = re.compile(r"^cyl_data_(.+)\.txt_.*_t\d+\.txt$", re.IGNORECASE)
    prefixes = sorted({name_pattern.match(f).group(1) for f in os.listdir(DATA_DIR)
                       if name_pattern.match(f)})
    if prefixes:
        print("Available tree prefixes in this folder:")
        for p in prefixes:
            print("   ", p)
    raise SystemExit

# Compute the three volumes (+ DBH/height/taper/trunk_len/branch_len) for every realization,
# PLUS the 10cm-diameter-cutoff variant of total/trunk/branch (see cylinder_metrics()).
totals, trunks, branches = [], [], []
dbhs, heights, tapers = [], [], []
trunk_lens, branch_lens = [], []
totals_10cm, trunks_10cm, branches_10cm = [], [], []
trunk_lens_10cm, branch_lens_10cm = [], []
n_totals, n_kepts = [], []
print("Tree: %s   (%d realization files)\n" % (TREE_NAME, len(files)))
print("%-24s %12s %12s %12s %8s %8s" % ("file", "total m^3", "trunk m^3", "branch m^3", "dbh cm", "h m"))
print("-" * 78)
for idx, path in files:
    m = cylinder_metrics(path)
    totals.append(m["total"]); trunks.append(m["trunk"]); branches.append(m["branch"])
    dbhs.append(m["dbh"]); heights.append(m["height"]); tapers.append(m["taper"])
    trunk_lens.append(m["trunk_len"]); branch_lens.append(m["branch_len"])
    totals_10cm.append(m["total_10cm"]); trunks_10cm.append(m["trunk_10cm"]); branches_10cm.append(m["branch_10cm"])
    trunk_lens_10cm.append(m["trunk_len_10cm"]); branch_lens_10cm.append(m["branch_len_10cm"])
    n_totals.append(m["n_total"]); n_kepts.append(m["n_kept"])
    print("%-24s %12.4f %12.4f %12.4f %8s %8.2f"
          % (os.path.basename(path), m["total"], m["trunk"], m["branch"],
             ("%.2f" % (m["dbh"] * 100.0)) if m["dbh"] is not None else "-", m["height"]))

# Summary statistics across the realizations.
print("-" * 78)
for label, vals in (("TOTAL", totals), ("TRUNK", trunks), ("BRANCHES", branches)):
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

# mean_trunk/mean_branch/mean_trunk_len/mean_branch_len (needed by the 10cm
# report right below, and again later for the unfiltered upsert_result call).
mean_trunk = summarize(trunks)[0]
mean_branch = summarize(branches)[0]
mean_trunk_len = mean_ignore_none(trunk_lens)
mean_branch_len = mean_ignore_none(branch_lens)

# ---- 10 cm cut-off comparison, averaged over all realizations -------------
# Same style as "Cylinders, cut-off 10 cm" elsewhere in the project
# (runsken.m section 17b, report_thin_branch_volume() in tree_geom_utils.py):
# cylinder count kept + volume kept/removed, split into trunk/branch too.
mean_total_10cm = summarize(totals_10cm)[0]
mean_trunk_10cm = summarize(trunks_10cm)[0]
mean_branch_10cm = summarize(branches_10cm)[0]
# trunk_len_10cm/branch_len_10cm: same summarize() averaging as the volumes
# right above, just applied to the trunk_lens_10cm/branch_lens_10cm lists.
mean_trunk_len_10cm = summarize(trunk_lens_10cm)[0]
mean_branch_len_10cm = summarize(branch_lens_10cm)[0]
mean_n_total = sum(n_totals) / len(n_totals)
mean_n_kept = sum(n_kepts) / len(n_kepts)
vol_removed = mean_total - mean_total_10cm
branch_removed = mean_branch - mean_branch_10cm
trunk_removed = mean_trunk - mean_trunk_10cm

print("\n--- Cylinders, cut-off %.0f cm (TreeQSM de Tanago, mean over %d realizations) ---"
      % (THIN_BRANCH_CUT_CM, len(files)))
print("Cylinders total (mean)  : %.0f" % mean_n_total)
print("Cylinders kept  (mean)  : %.0f (%.1f %%)"
      % (mean_n_kept, (mean_n_kept / mean_n_total * 100.0) if mean_n_total else 0.0))
print("Volume total (mean)     : %.3f m3" % mean_total)
print("Volume kept  (mean)     : %.3f m3" % mean_total_10cm)
print("Volume removed (mean)   : %.3f m3 (%.1f %%)"
      % (vol_removed, (vol_removed / mean_total * 100.0) if mean_total else 0.0))
print("Trunk  volume kept (mean)  : %.3f m3  (removed %.3f m3, %.1f %%)"
      % (mean_trunk_10cm, trunk_removed, (trunk_removed / mean_trunk * 100.0) if mean_trunk else 0.0))
print("Branch volume kept (mean)  : %.3f m3  (removed %.3f m3, %.1f %%)"
      % (mean_branch_10cm, branch_removed, (branch_removed / mean_branch * 100.0) if mean_branch else 0.0))

# ---- same cut-off, but on LENGTH instead of volume (this turn's addition) --
trunk_len_removed = mean_trunk_len - mean_trunk_len_10cm
branch_len_removed = mean_branch_len - mean_branch_len_10cm
print("Trunk  length kept (mean)  : %.3f m    (removed %.3f m, %.1f %%)"
      % (mean_trunk_len_10cm, trunk_len_removed,
         (trunk_len_removed / mean_trunk_len * 100.0) if mean_trunk_len else 0.0))
print("Branch length kept (mean)  : %.3f m    (removed %.3f m, %.1f %%)"
      % (mean_branch_len_10cm, branch_len_removed,
         (branch_len_removed / mean_branch_len * 100.0) if mean_branch_len else 0.0))

# ---- upsert this result into the shared master results CSV -----------------
std_total = summarize(totals)[1]
# (mean_trunk/mean_branch/mean_trunk_len/mean_branch_len already computed
# above, right before the 10cm report)
# branch_filter = "none": these are the full TreeQSM cylinder realizations,
# not restricted to any diameter cut-off.
upsert_result(RESULTS_CSV, TREE_NAME, METHOD_LABEL,
              mean_total, mean_trunk, mean_branch, std_total,
              mean_dbh, mean_height, mean_taper,
              mean_trunk_len, mean_branch_len, branch_filter="none",
              # n_cylinders: mean_n_total (already computed above from
              # n_totals, one count per realization file) rounded to the
              # nearest whole number - a cylinder count can't be fractional,
              # same reasoning as runsken.m's round(mean(...)) for N_cylinders.
              n_cylinders=round(mean_n_total))

# ---- SECOND row: same realizations, but with the 10 cm cut-off applied ----
# total/trunk/branch = the FILTERED mean/std computed above (summarize(), the
# exact same function used for the unfiltered row - just applied to the
# totals_10cm/trunks_10cm/branches_10cm lists instead of totals/trunks/branches).
# trunk_len/branch_len = mean_trunk_len_10cm/mean_branch_len_10cm, same idea
# but for length (added this turn, so Trunk/Branch length charts and
# field_error_summary() also work for this row in the vs.-reference mode).
# dbh/height/taper are NOT recomputed: a 10 cm diameter cut-off barely touches
# the trunk near 1.3 m (it's always far thicker than 10 cm there) or the tree's
# overall height, so the unfiltered values above are reused as-is, per your
# instruction not to recompute something a cut-off this coarse wouldn't change.
std_total_10cm = summarize(totals_10cm)[1]
upsert_result(RESULTS_CSV, TREE_NAME, METHOD_LABEL_10CM,
              mean_total_10cm, mean_trunk_10cm, mean_branch_10cm, std_total_10cm,
              mean_dbh, mean_height, mean_taper,
              mean_trunk_len_10cm, mean_branch_len_10cm,
              branch_filter="10cm",
              # n_cylinders: mean_n_kept (already computed above from
              # n_kepts) - the count of cylinders that passed the >=10cm
              # filter, mirroring "n_cyl_kept" in tree_geom_utils.py's
              # report_thin_branch_volume() for the same filtered-row idea.
              n_cylinders=round(mean_n_kept))

# Export the per-realization results and summary statistics to CSV.
if CSV_PATH:
    with open(CSV_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["tree", "file", "total_m3", "trunk_m3", "branch_m3"])
        for (idx, path), total, trunk, branch in zip(files, totals, trunks, branches):
            writer.writerow([TREE_NAME, os.path.basename(path),
                              "%.6f" % total, "%.6f" % trunk, "%.6f" % branch])
        writer.writerow([])
        writer.writerow(["stat", "total_m3", "trunk_m3", "branch_m3"])
        for stat_name, stat_idx in (("mean", 0), ("std", 1), ("cv_percent", 2)):
            row = [stat_name]
            for vals in (totals, trunks, branches):
                row.append("%.6f" % summarize(vals)[stat_idx])
            writer.writerow(row)
        writer.writerow(["min"] + ["%.6f" % min(vals) for vals in (totals, trunks, branches)])
        writer.writerow(["max"] + ["%.6f" % max(vals) for vals in (totals, trunks, branches)])
    print("\nCSV exported to: %s" % CSV_PATH)
