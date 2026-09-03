# -*- coding: utf-8 -*-
# =====================================================================
#  Add a TreeGraph result to the shared master results table
#  (volume_results.csv, the same file compare_volumes.py reads).
# ---------------------------------------------------------------------
#  This script does NOT do any comparing/printing itself - it only reads
#  one TreeGraph .json (from tree2qsm.py) and upsert_result()s a row into
#  RESULTS_CSV, exactly like the AdTree/TreeQSM scripts already do. Once
#  the row is there, run compare_volumes.py as usual to see it side by
#  side with every other method/tree.
#
#  Re-running this script for the SAME tree+method OVERWRITES that row
#  (upsert_result()'s normal behaviour) instead of duplicating it - so if
#  you re-run tree2qsm.py with the same cluster_size, just re-run this too.
#
#  TAPER: computed the same way as the AdTree pipeline (stem diameter at
#  1.3 m vs. 10.0 m above the base, in cm/m), but from TreeGraph's own
#  "cyls" table via treegraph_stem_diameter_at_height() - see
#  tree_geom_utils.py.
#
#  Dependencies: pandas, numpy
# =====================================================================

from tree_geom_utils import load_treegraph_json, treegraph_stem_diameter_at_height, upsert_result

# =====================  PARAMETERS  ===================================
# Path to the TreeGraph output .json to add (from tree2qsm.py).
JSON_PATH = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\01_TreeGraph\treegraph_work\results\test_smoke\IND01_054-cs0.15-tipNone.json"

# The master results table (same file compare_volumes.py reads/writes).
RESULTS_CSV = "volume_results.csv"

# Method label for this variant in the CSV. Leave as None to build one
# automatically from the run's cluster_size, e.g. "TreeGraph (cs0.15)".
METHOD_LABEL = None

# Heights (m above the tree base) used for the taper metric, matching the
# AdTree pipeline's convention (see tree_geom_utils.stem_diameter_at_height).
TAPER_H_LOW = 1.3
TAPER_H_HIGH = 10.0

# "none" = full/unfiltered reconstruction (TreeGraph applies no diameter
# cut-off), matching the branch_filter convention in compare_volumes.py.
BRANCH_FILTER = "none"
# =====================================================================


cyls, tree_df, args, name = load_treegraph_json(JSON_PATH)
tree = tree_df.iloc[0]

if METHOD_LABEL is None:
    METHOD_LABEL = "TreeGraph (cs%s)" % args.get("cluster_size")

total_vol = float(tree["vol"])
trunk_vol = float(tree["trunk_vol"])
branch_vol = total_vol - trunk_vol
trunk_len = float(tree["trunk_length"])
branch_len = float(tree["length"]) - trunk_len
dbh = float(tree["DBH_from_qsm"])
height = float(tree["H_from_qsm"])

# TreeGraph's own base-correction already puts the tree base at the lowest
# point of the (corrected) skeleton, and DBH_from_qsm is measured at
# dbh_height (normally 1.3 m) above THAT base - so z_base here is the
# minimum start-height among branch_order==0 (trunk) cylinders, consistent
# with how TreeGraph itself measured DBH.
z_base = float(cyls[cyls.branch_order == 0].sz.min())

d_low = treegraph_stem_diameter_at_height(cyls, z_base, TAPER_H_LOW)
d_high = treegraph_stem_diameter_at_height(cyls, z_base, TAPER_H_HIGH)
taper = None
if d_low is not None and d_high is not None and (TAPER_H_HIGH - TAPER_H_LOW) > 0:
    taper = (d_low - d_high) * 100.0 / (TAPER_H_HIGH - TAPER_H_LOW)  # cm/m

upsert_result(RESULTS_CSV, tree=name, method=METHOD_LABEL,
              total=total_vol, trunk=trunk_vol, branch=branch_vol, std=None,
              dbh=dbh, height=height, taper=taper,
              trunk_len=trunk_len, branch_len=branch_len,
              branch_filter=BRANCH_FILTER, n_cylinders=len(cyls))

print("Added/updated row: tree=%s, method=%s" % (name, METHOD_LABEL))
print("  total_m3=%.4f  trunk_m3=%.4f  branch_m3=%.4f  n_cylinders=%d"
      % (total_vol, trunk_vol, branch_vol, len(cyls)))
print("  dbh_m=%.3f  height_m=%.2f  taper_cm_per_m=%s"
      % (dbh, height, ("%.2f" % taper) if taper is not None else "n/a"))
print("Now run compare_volumes.py to see it alongside the other methods.")
