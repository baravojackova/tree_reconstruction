# -*- coding: utf-8 -*-
# =====================================================================
#  Export a TreeGraph .json result to geom_*.txt for ANSYS (bk1.mac).
# ---------------------------------------------------------------------
#  Unlike the AdTree pipeline (adtree_reconstruct_compare.py + export_geom_
#  ansys.py), there is no calibration step here and no geometry to
#  recompute: TreeGraph's own "cyls" table already has start point, unit
#  axis, length and radius for every cylinder, so this script just loads
#  it and reformats it via write_geom_from_treegraph() (see
#  tree_geom_utils.py for what it does and why cylinder ids get renumbered).
#
#  Dependencies: pandas, numpy
# =====================================================================

import os

from tree_geom_utils import load_treegraph_json, write_geom_from_treegraph

# =====================  PARAMETERS  ===================================
# Path to the TreeGraph output .json to export (from tree2qsm.py).
JSON_PATH = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\01_TreeGraph\treegraph_work\results\test_smoke\IND01_054-cs0.15-tipNone.json"

# Where to write the geom.txt. Leave as None to write it next to JSON_PATH,
# named "geom_<treeid>-cs<cluster_size>.txt".
GEOM_PATH = None

# Recentre the model so the trunk base sits at x=0, y=0 (z is left as-is).
# Matches the convention used by the AdTree pipeline's write_geom() calls.
RECENTER_XY = True
# =====================================================================


cyls, tree, args, name = load_treegraph_json(JSON_PATH)

if GEOM_PATH is None:
    cs = args.get("cluster_size")
    GEOM_PATH = os.path.join(os.path.dirname(JSON_PATH), "geom_%s-cs%s.txt" % (name, cs))

write_geom_from_treegraph(GEOM_PATH, cyls, recenter_xy=RECENTER_XY)

print("Tree: %s" % name)
print("Cylinders written: %d" % len(cyls))
print("Wrote: %s" % GEOM_PATH)
