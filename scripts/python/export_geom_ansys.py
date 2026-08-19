# -*- coding: utf-8 -*-
# =====================================================================
#  Step 2 of the AdTree -> ANSYS pipeline: turn one "calib_*.npz" file
#  (written by adtree_reconstruct_compare.py) into the geom_*.txt that
#  ANSYS (buk.mac) actually reads.
# ---------------------------------------------------------------------
#  This script does NO calibration and NO comparison - it just loads the
#  geometry that adtree_reconstruct_compare.py already calibrated and
#  saved, and writes it out in the geom.txt format via write_geom() from
#  tree_geom_utils.py. That's the exact same function/format the old
#  single-file ply_to_geom.py used to call directly, so the output file
#  is identical to before - only WHEN it gets written has changed (now a
#  separate, fast step you can re-run any time without redoing calibration).
#
#  Dependencies: numpy   (install: pip install numpy)
# =====================================================================

import glob
import os

import numpy as np

from tree_geom_utils import write_geom

# =====================  PARAMETERS  ===================================
# Which .npz file (written by adtree_reconstruct_compare.py) to convert to
# geom_*.txt. Leave as None to just LIST the available files and stop (so
# you can see what's there before picking one) - set it to one of the
# printed names to actually export it, e.g.:
#   NPZ_FILE = "calib_IND01_054_r5mm.npz"
NPZ_FILE = "calib_IND01_054_r5mm.npz"

# Folder to look for "calib_*.npz" files in. "." means "the folder this
# script is run from" (matches how adtree_reconstruct_compare.py writes them).
NPZ_DIR = "."
# =====================================================================


def list_available_npz(npz_dir):
    """Print every "calib_*.npz" file found in npz_dir, so you can see what's
    available before setting NPZ_FILE above. Returns the sorted list of paths."""
    pattern = os.path.join(npz_dir, "calib_*.npz")
    paths = sorted(glob.glob(pattern))
    if not paths:
        print("No 'calib_*.npz' files found in '%s'." % npz_dir)
        print("Run adtree_reconstruct_compare.py first to create one.")
        return paths

    print("Available calibrated-geometry files in '%s':" % npz_dir)
    for path in paths:
        # Peek at each file's metadata (tree/variant/threshold) without
        # loading the (potentially large) xyz/cyl arrays into memory twice.
        with np.load(path, allow_pickle=False) as data:
            tree_name = str(data["tree_name"])
            variant_label = str(data["variant_label"])
            threshold_mm = round(float(data["threshold_m"]) * 1000)
            geom_filename = str(data["geom_filename"])
        variant_txt = (" variant=%s" % variant_label) if variant_label else ""
        print("  %-40s tree=%s%s threshold=%dmm -> would write '%s'"
              % (os.path.basename(path), tree_name, variant_txt, threshold_mm, geom_filename))
    return paths


def load_calibrated_geometry(path):
    """Load one "calib_*.npz" file back into the (xyz, cyl, root, recenter_xy)
    shape write_geom() expects. This is the exact inverse of the np.savez(...)
    call in adtree_reconstruct_compare.py - see the big comment there for what
    each field means and why it's enough to reconstruct the geometry losslessly."""
    with np.load(path, allow_pickle=False) as data:
        xyz = data["xyz"]                              # (N,3) float64, unchanged
        cyl_array = data["cyl"]                        # (n_cyl,4) float64: a,b,r,pid
        root = int(data["root"])                       # stored as a 0-d array -> plain int
        recenter_xy = bool(data["recenter_xy"])         # stored as a 0-d array -> plain bool
        geom_filename = str(data["geom_filename"])      # stored as a 0-d string array -> plain str

    # cyl_array's columns are all float64 (np.savez needs one dtype per array),
    # so a/b/pid (which are really integer node/cylinder indices) are cast
    # back to int here - write_geom() indexes xyz with them, which requires ints.
    cyl = [(int(a), int(b), float(r), int(pid)) for a, b, r, pid in cyl_array]
    return xyz, cyl, root, recenter_xy, geom_filename


# =========================  RUN  =====================================
available = list_available_npz(NPZ_DIR)

if NPZ_FILE is None:
    print("\nSet NPZ_FILE above to one of the files listed, then re-run this script.")
    raise SystemExit

npz_path = NPZ_FILE if os.path.dirname(NPZ_FILE) else os.path.join(NPZ_DIR, NPZ_FILE)
if not os.path.exists(npz_path):
    print("'%s' not found." % npz_path)
    raise SystemExit

xyz, cyl, root, recenter_xy, geom_filename = load_calibrated_geometry(npz_path)
print("Loaded %s: %d cylinders, root=%d, recenter_xy=%s" % (npz_path, len(cyl), root, recenter_xy))

write_geom(geom_filename, xyz, cyl, root, recenter_xy)
print("Wrote:", geom_filename)
