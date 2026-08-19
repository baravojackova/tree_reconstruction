# -*- coding: utf-8 -*-
# =====================================================================
#  Convert an AdTree skeleton (.ply) to the geom.txt format for ANSYS
# ---------------------------------------------------------------------
#  What the script does (step by step):
#   1) Reads a binary .ply (vertices with x,y,z,radius and edges).
#   2) Merges coincident vertices -> turns the disconnected skeleton into
#      a single connected tree with real branch points.
#   3) Roots the tree at the base (lowest z) and finds a unique parent
#      for every node (this removes any loops).
#   4) Smooths the centerline (Laplacian smoothing along branches) to
#      remove zig-zag noise from the AdTree skeleton.
#   5) Prunes thin twigs below a chosen radius threshold.
#   6) Resamples the dense skeleton into segments whose length adapts to
#      the local branch radius (keeps the curvature but drastically
#      reduces the element count).
#   7) Shifts the x,y coordinates to the origin (z is left untouched,
#      since it is the height above ground).
#   8) Optionally CALIBRATES the cylinder radii against AdQSM data (a taper
#      curve for the trunk, and median branch diameters per branch order for
#      the rest), replacing the AdTree radii before writing.
#   9) Writes geom.txt with columns: x y z length radius axis_x axis_y axis_z parent
#      -> exactly what your buk.mac reads.
#   10) Optionally shows/saves a 3D preview of the reduced model so you
#       can inspect it before importing into ANSYS.
#
#  Dependencies: numpy, scipy, matplotlib   (install: pip install numpy scipy matplotlib)
# =====================================================================

import os
import re
import csv

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components

# =====================  PARAMETERS  ===================================
# Directory holding this tree's source input files (the AdTree skeleton .ply,
# and the AdQSM taper/branch/params exports). Edit this to point at a
# different tree's folder; INPUT_PLY and the ADQSM_*_FILE paths below are all
# resolved relative to it.
AdQSM_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data\IND01_54\AdQSM\05"
AdTree_DIR = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data\IND01_54"

INPUT_PLY = os.path.join(AdTree_DIR, "IND01_054 - Cloud_skeleton.ply")   # input skeleton from AdTree

# Radius threshold (in METERS). You can give several values -> several variants.
# Example of a single variant:   RADIUS_THRESHOLDS = [0.010]
# Example of several variants:   RADIUS_THRESHOLDS = [0.010, 0.020, 0.030]
RADIUS_THRESHOLDS = [0.005]        # 0.030 m = 30 mm radius (60 mm diameter)

# Adaptive segment length used for resampling (in METERS). The target length at
# a given point is SEG_LEN_K * local_radius, clamped to [SEG_LEN_MIN, SEG_LEN_MAX]:
# thick branches (trunk) get long segments, thin twigs get short/fine ones.
# If SEG_LEN_MIN == SEG_LEN_MAX, this reduces to the old constant-length
# resampling (that fixed value, regardless of radius).
# Set SEG_LEN_MIN to 0 or None to disable resampling entirely (keep every point).
SEG_LEN_MIN = 0.01                # shortest allowed segment (m), for thin twigs
SEG_LEN_MAX = 0.3                # longest allowed segment (m), for the trunk
SEG_LEN_K = 0.5                   # target length = SEG_LEN_K * local_radius 

# Shortest permissible cylinder (in METERS). Shorter ones (e.g. at branch
# points) are not created, to avoid degenerate zero-length beams in ANSYS.
MIN_CYL_LEN = 0.005               # 0.1 mm

# Rounding used when merging coincident points (number of decimal places).
# 5 is a safe default.
MERGE_DECIMALS = 5

# Shift x,y to zero (removes a large georeferencing offset)? Recommended True.
RECENTER_XY = True

# --- Centerline smoothing (applied BEFORE pruning/resampling) --------
# Laplacian smoothing removes zig-zag noise from the AdTree skeleton. Each
# smoothing pass moves every free point along a branch toward the average of
# its two neighbours by SMOOTH_ALPHA. The root, junctions (degree >= 3), and
# branch tips (leaves) are never moved, so topology and branch endpoints stay
# fixed. Radii are not affected.
SMOOTH_ITERS = 5                  # number of smoothing passes (0 = off)
SMOOTH_ALPHA = 0.5                # 0..1 strength per pass

# --- Radius calibration against AdQSM (optional) ----------------------
# Replaces the AdTree skeleton radii with AdQSM-calibrated radii on the final
# (pruned + resampled) cylinders. The trunk (branch order 0) is taken from the
# AdQSM taper curve; every other branch order is scaled by a single factor
# (AdQSM median radius / AdTree median radius) computed for that order. Set
# CALIBRATE_RADII = False to skip this step entirely and use the raw AdTree
# radii, exactly like before this feature existed.
CALIBRATE_RADII = True

# AdQSM trunk taper curve: text file of "height[m] <TAB> diameter[m]" rows.
ADQSM_TAPER_FILE = os.path.join(AdQSM_DIR, "taper.txt")

# AdQSM per-branch table: text file with columns order, parent, diameter[m], ...
ADQSM_BRANCH_FILE = os.path.join(AdQSM_DIR, "BranchStructure.txt")

# AdQSM tree parameters (volumes in m^3), used only for the reference line in
# the volume comparison printout.
ADQSM_PARAMS_FILE = os.path.join(AdQSM_DIR, "TreesParams.txt")

# Measured trunk diameter at breast height (1.3 m), in METERS. If given, the
# taper curve is rescaled so its value at 1.3 m matches this measurement.
# Set to None to use the taper curve exactly as read from ADQSM_TAPER_FILE.
FIELD_DBH = None

# Tree ID used for this tree's rows in the shared results table (RESULTS_CSV).
TREE_NAME = "IND01_054"

# Shared master results table (see compare_volumes.py). When CALIBRATE_RADII
# is True, each generated threshold variant upserts its own row into this CSV
# so results from all methods live in one place.
RESULTS_CSV = "volume_results.csv"

# Reference heights [m] used for DBH (lower) and the taper metric (lower/
# upper). DBH is the stem diameter at TAPER_H_LOWER (1.3 m = breast height).
TAPER_H_LOWER = 1.3    # lower reference height [m]
TAPER_H_UPPER = 10.0   # upper reference height [m]
# =====================================================================

# --- Output file name(s) --------------------------------------------
# OUTPUT_NAME is used as-is when RADIUS_THRESHOLDS has exactly ONE value.
OUTPUT_NAME = "geom_r30_optim.txt"
# OUTPUT_PATTERN is used instead when RADIUS_THRESHOLDS has SEVERAL values,
# so the generated files don't overwrite each other. {r} is replaced by the
# radius threshold in mm (e.g. geom_r10.txt).
OUTPUT_PATTERN = "geom_r{r}.txt"
# =====================================================================

# =====================  VISUALIZATION  ================================
# Show an interactive 3D preview of the reduced beam model (one window per
# radius threshold) so you can check it before importing into ANSYS.
SHOW_PLOT = True

# Also save the preview as a PNG next to each output file (True/False).
# The PNG name is derived from the output file name (.txt -> .png).
SAVE_PLOT_PNG = False
# =====================================================================


def read_ply(path):
    """Read a binary .ply with vertices (x,y,z,radius) and edges (2 indices)."""
    raw = open(path, "rb").read()
    he = raw.find(b"end_header\n") + len(b"end_header\n")
    header = raw[:he].decode("ascii", "replace")
    body = raw[he:]

    nV = nE = None
    for line in header.splitlines():             #   find the vertex/edge counts in the header
        if line.startswith("element vertex"):    # the AdTree .ply uses "vertex" instead of "point"
            nV = int(line.split()[-1])
        elif line.startswith("element edge"):    # the AdTree .ply uses "edge" instead of "face"
            nE = int(line.split()[-1])
    if nV is None or nE is None:
        raise ValueError("Could not find the vertex/edge counts in the header.")

    # vertices: 4 float32 (x, y, z, radius)
    verts = np.frombuffer(body[:nV * 16], dtype="<f4").reshape(nV, 4)    # reshape vertices into (nV, 4) array
    xyz = verts[:, :3].astype(np.float64) 
    rad = verts[:, 3].astype(np.float64)

    # edges: each = uint32 count(=2) + 2x int32 index
    edges = np.frombuffer(body[nV * 16:], dtype="<u4").reshape(nE, 3)[:, 1:3].astype(np.int64)
    return xyz, rad, edges


def merge_vertices(xyz, rad, edges, decimals): 
    """Merge points at the same (rounded) location -> joins branches at forks."""
    key = np.round(xyz, decimals)
    uniq, inv = np.unique(key, axis=0, return_inverse=True) # find unique points and their indices
    U = len(uniq)                                           # number of unique points
    # average radius for each merged point
    rsum = np.zeros(U); cnt = np.zeros(U)   # put zero arrays on the same device as the input arrays
    np.add.at(rsum, inv, rad); np.add.at(cnt, inv, 1.0) # compute the average radius for each unique point
    rad_u = rsum / cnt
    # remap edges to the new indices, drop self-loops and duplicates
    e2 = inv[edges]              # inv maps old indices to new unique indices
    e2 = e2[e2[:, 0] != e2[:, 1]]   # drop self-loops (edges from a vertex to itself)
    e2 = np.unique(np.sort(e2, axis=1), axis=0)
    return uniq.astype(np.float64), rad_u, e2       # return the unique points, their averaged radii, and the new edges


def build_rooted_tree(nnodes, edges, xyz):
    """Root the graph at the lowest point and return the 'parent' array for each node."""
    root = int(np.argmin(xyz[:, 2]))          # base = lowest z
    A = csr_matrix((np.ones(len(edges) * 2),                        # build a sparse adjacency matrix for the graph
                    (np.concatenate([edges[:, 0], edges[:, 1]]),    # rows
                     np.concatenate([edges[:, 1], edges[:, 0]]))),  # columns
                   shape=(nnodes, nnodes))      
    # BFS tree from the root -> predecessors give a unique parent (also resolves loops)
    order, pred = breadth_first_order(A, root, directed=False, return_predecessors=True)
    return root, pred, order    # return the root index, parent array, and BFS order of nodes


def smooth_centerline(xyz, edges, root, iters, alpha):
    """Laplacian-smooth the vertex coordinates along degree-2 chains of the tree.

    Only "free" points (degree == 2, not the root) are moved, each pass pulling
    them toward the average of their two neighbours by `alpha`. Junctions
    (degree >= 3) and branch tips (degree <= 1, including the root) stay fixed,
    so topology and branch endpoints are unaffected. Radii are untouched.
    """
    if iters is None or iters <= 0 or alpha is None or alpha <= 0:
        return xyz

    n = len(xyz)    # number of nodes
    adj = [[] for _ in range(n)]    # build adjacency list
    for a, b in edges:  
        adj[a].append(b)
        adj[b].append(a)
    degree = np.array([len(a) for a in adj])    # compute the degree of each node

    movable = degree == 2   # only degree-2 nodes are movable; root is never moved
    movable[root] = False
    idx = np.where(movable)[0]  # get indices of movable nodes
    if len(idx) == 0:
        return xyz

    nb0 = np.array([adj[i][0] for i in idx]) # get the first neighbor of each movable node
    nb1 = np.array([adj[i][1] for i in idx]) # get the second neighbor of each movable node

    xyz = xyz.copy()
    for _ in range(iters):      # perform the smoothing iterations
        avg = 0.5 * (xyz[nb0] + xyz[nb1])   
        xyz[idx] = (1.0 - alpha) * xyz[idx] + alpha * avg   # move each movable node toward the average of its two neighbors
    return xyz


def local_seg_len(r, seg_min, seg_max, seg_k):  
    """Target resampling segment length at a point of local radius `r`."""
    if seg_min == seg_max:
        return seg_min
    return min(max(seg_k * r, seg_min), seg_max)    # clamp the segment length to [seg_min, seg_max]


def compute_branch_order(nnodes, kids, rad, bfs_order):
    """Branch order per node on a rooted tree: the root is order 0. At each
    node, the child with the LARGEST radius keeps the parent's order; every
    other child starts a new branch at order+1. `kids` is a list-of-lists of
    children (already restricted to whatever subtree is of interest, e.g. the
    pruned tree), and `bfs_order` must visit every node after its parent
    (as returned by build_rooted_tree) so orders propagate top-down."""
    node_order = np.zeros(nnodes, dtype=np.int64)   # initialize all nodes to order 0
    for n in bfs_order:  # iterate over nodes in BFS order (parent before children)
        c = kids[n]     # get the children of node n
        if not c:
            continue
        thick = c[int(np.argmax([rad[k] for k in c]))]  # find the child with the largest radius (the "thick" child)
        for k in c:
            node_order[k] = node_order[n] if k == thick else node_order[n] + 1  # assign the same order to the thick child, and order+1 to the others
    return node_order


def parse_adqsm_taper_file(path):
    """Parse an AdQSM taper.txt file: rows of 'height[m] <TAB> diameter[m]'.
    Robust to non-UTF8 (Chinese) headers, blank lines, and several
    concatenated blocks (the file may repeat the same export multiple times)
    - only the FIRST block of numeric rows is used, per AdQSM convention."""
    heights, diameters = [], []
    collecting = False
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                try:
                    h, d = float(parts[0]), float(parts[1])
                except ValueError:
                    if collecting:
                        break
                    continue
                heights.append(h)
                diameters.append(d)
                collecting = True
            elif collecting:
                break
    if not heights:
        raise ValueError("No numeric height/diameter rows found in %s" % path)
    idx = np.argsort(heights)
    return np.asarray(heights)[idx], np.asarray(diameters)[idx]


def make_trunk_radius_func(taper_heights, taper_diameters, field_dbh=None):
    """Build a trunk_radius(height) -> radius function from a taper curve,
    via linear interpolation (clamped to the end values outside the measured
    range). If `field_dbh` (a measured trunk diameter at 1.3 m) is given, the
    whole curve is rescaled so its value at 1.3 m matches it."""
    radii = taper_diameters / 2.0
    valid = radii > 0                       # drop invalid taper rows (e.g. diameter 0 at the top)
    taper_heights = taper_heights[valid]
    radii = radii[valid]
    if field_dbh is not None:
        r_at_dbh = float(np.interp(1.3, taper_heights, radii))
        if r_at_dbh > 0:
            radii = radii * ((field_dbh / 2.0) / r_at_dbh)

    def trunk_radius(height):
        return float(np.interp(height, taper_heights, radii))

    return trunk_radius


def parse_adqsm_branch_file(path):
    """Parse an AdQSM BranchStructure.txt file and return the MEDIAN AdQSM
    radius per branch order. Only lines that split by TAB into >=7 fields
    whose first field is an integer are used (column 0 = order, column 2 =
    diameter [m]); the file may contain several concatenated runs, all
    matching rows are used."""
    diam_by_order = {}
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 7:
                continue
            try:
                order = int(parts[0])
                diameter = float(parts[2])
            except ValueError:
                continue
            diam_by_order.setdefault(order, []).append(diameter / 2.0)
    return {o: float(np.median(r)) for o, r in diam_by_order.items()}


def parse_adqsm_params_file(path):
    """Parse an AdQSM TreesParams.txt file and return a dict shaped like the
    ones `volume_stats` produces (so it can be printed by `print_volume_stats`),
    built from the tree-level volumes/lengths AdQSM itself reports.

    The data is a single line of TAB-separated "Key: value" tokens. The file
    may contain several concatenated blocks - only the FIRST block containing
    a "TreeVolume:" token is used. Returns None if the file is missing or no
    TreeVolume is found (never raises)."""
    if not os.path.exists(path):
        return None

    with open(path, "r", encoding="latin-1") as f:
        lines = f.readlines()

    def extract(key, text):
        m = re.search(key + r":\s*([-\d.eE+]+)", text)
        return float(m.group(1)) if m else None

    for line in lines:
        if "TreeVolume:" not in line:
            continue
        tree_vol = extract("TreeVolume", line)
        if tree_vol is None:
            continue
        trunk_vol = extract("TrunkVolume", line) or 0.0
        branch_vol = extract("BranchVolume", line) or 0.0
        trunk_len = extract("TrunkLength", line) or 0.0
        branch_len = extract("BranchLength", line) or 0.0
        branches_num = extract("BranchesNum", line) or 0.0
        tree_height = extract("TreeHeight", line)
        return dict(n=int(branches_num), total_len=trunk_len + branch_len, total_vol=tree_vol,
                    trunk_n=0, trunk_len=trunk_len, trunk_vol=trunk_vol,
                    branch_n=0, branch_len=branch_len, branch_vol=branch_vol,
                    height=tree_height)

    return None


def calibrate_cylinder_radii(xyz, cyl, cyl_order, trunk_radius_func, adqsm_median_by_order):
    """Compute AdQSM-calibrated radii for the final cylinders.

    Order 0 (trunk) cylinders get their radius from the taper curve, evaluated
    at the cylinder's midpoint height. Every other order is scaled by a single
    factor = AdQSM_median_radius(order) / AdTree_median_radius(order), where
    the AdTree median is computed over these same cylinders; orders absent
    from the AdQSM table get factor 1 (unchanged). Returns (new_radii, factors)
    where new_radii is an array parallel to `cyl` and factors maps order -> factor.
    """
    orders = np.asarray(cyl_order)
    orig_r = np.array([c[2] for c in cyl])
    new_r = orig_r.copy()

    factors = {}
    for o in sorted(set(orders.tolist()) - {0}):
        mask = orders == o
        adtree_med = float(np.median(orig_r[mask]))
        adqsm_med = adqsm_median_by_order.get(o)
        factor = (adqsm_med / adtree_med) if (adqsm_med is not None and adtree_med > 0) else 1.0
        factors[o] = factor
        new_r[mask] = orig_r[mask] * factor

    z_base = float(xyz[:, 2].min())   # tree base; AdQSM taper heights are measured from here
    trunk_idx = np.nonzero(orders == 0)[0]
    for i in trunk_idx:
        a, b = cyl[i][0], cyl[i][1]
        h = 0.5 * (xyz[a, 2] + xyz[b, 2]) - z_base   # height ABOVE the tree base
        new_r[i] = trunk_radius_func(h)

    new_r = np.maximum(new_r, 1e-6)   # never allow radius <= 0
    return new_r, factors


def volume_stats(lengths, radii, order_arr):
    """Cylinder-volume summary (pi * r^2 * length), split into trunk (order 0)
    vs. branches (order >= 1)."""
    volumes = np.pi * radii ** 2 * lengths
    n = len(lengths)
    total_len, total_vol = float(lengths.sum()), float(volumes.sum())
    trunk_mask = order_arr == 0
    t_n = int(trunk_mask.sum())
    t_len, t_vol = float(lengths[trunk_mask].sum()), float(volumes[trunk_mask].sum())
    return dict(n=n, total_len=total_len, total_vol=total_vol,
                trunk_n=t_n, trunk_len=t_len, trunk_vol=t_vol,
                branch_n=n - t_n, branch_len=total_len - t_len, branch_vol=total_vol - t_vol)


def upsert_result(csv_path, tree, method, total, stem, branch, std, dbh=None, height=None, taper=None):
    """Insert/update one (tree, method) row in the shared master results CSV
    (see compare_volumes.py for its format). Reads csv_path if it exists
    (creating it with the header if not), removes any existing row with the
    same tree AND method, appends the new row, and writes the file back -
    so re-running a script overwrites its previous result instead of
    duplicating it. Numbers are formatted with 6 decimals; None -> "" (blank).
    Backward compatible: if csv_path still has the old 6-column header, its
    rows are read fine and rewritten with the new columns added as blank."""
    header = ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3",
              "dbh_m", "height_m", "taper_cm_per_m"]

    def fmt(x):
        return "" if x is None else "%.6f" % x

    rows = []
    if os.path.exists(csv_path):
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    rows = [r for r in rows if not (r["tree"] == tree and r["method"] == method)]
    rows.append({"tree": tree, "method": method, "total_m3": fmt(total),
                 "stem_m3": fmt(stem), "branch_m3": fmt(branch), "std_m3": fmt(std),
                 "dbh_m": fmt(dbh), "height_m": fmt(height), "taper_cm_per_m": fmt(taper)})

    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)


def stem_diameter_at_height(xyz, cyl, cyl_order, z_base, h):
    """Diameter [m] of the trunk (branch order 0 cylinders) at height h [m]
    above the tree base z_base. Uses the cylinder whose height span
    [min(za,zb), max(za,zb)] contains h, falling back to the nearest
    cylinder for small gaps between spans. Returns None if h is outside the
    height range actually covered by trunk cylinders (no extrapolation/
    guessing). Reused for both DBH and the taper metric."""
    spans = []
    for (a, b, r, pid), o in zip(cyl, cyl_order):
        if o != 0:
            continue
        za, zb = xyz[a, 2] - z_base, xyz[b, 2] - z_base
        spans.append((min(za, zb), max(za, zb), r))
    if not spans:
        return None
    lo = min(s[0] for s in spans)
    hi = max(s[1] for s in spans)
    if h < lo or h > hi:
        return None
    for h_start, h_end, r in spans:
        if h_start <= h <= h_end:
            return 2.0 * r

    def dist(span):
        h_start, h_end, _ = span
        return h_start - h if h < h_start else h - h_end

    nearest = min(spans, key=dist)
    return 2.0 * nearest[2]


def _fmt_dbh(x):
    return "%.4f m" % x if x is not None else "n/a"


def _fmt_taper(x):
    return "%.2f cm/m" % x if x is not None else "n/a"


def print_volume_stats(label, s):
    print("    %-26s: %6d cyl, length %8.1f m, volume %7.3f m^3 (%8.0f L)  "
          "[trunk %6.3f m^3 | branches %6.3f m^3]"
          % (label, s["n"], s["total_len"], s["total_vol"], s["total_vol"] * 1000.0,
             s["trunk_vol"], s["branch_vol"]))


def raw_skeleton_stats(xyz_raw, rad, edges):
    """Volume/length stats for the full, un-pruned, un-resampled AdTree
    skeleton (as a single merged tree, before smoothing), split into trunk
    (order 0) vs. branches with the same branch-order rule used for the final
    cylinders. This is the baseline the calibration is checked against."""
    nnodes = len(xyz_raw)
    root, parent, bfs = build_rooted_tree(nnodes, edges, xyz_raw)
    children = [[] for _ in range(nnodes)]
    for v in bfs:
        p = parent[v]
        if p >= 0:
            children[p].append(v)
    node_order = compute_branch_order(nnodes, children, rad, bfs)

    child_idx = np.nonzero(parent >= 0)[0]
    parent_idx = parent[child_idx]
    lengths = np.linalg.norm(xyz_raw[child_idx] - xyz_raw[parent_idx], axis=1)
    radii = 0.5 * (rad[child_idx] + rad[parent_idx])
    orders = node_order[child_idx]
    return volume_stats(lengths, radii, orders)


def convert(xyz, rad, edges, thr, seg_len_min, seg_len_max, seg_len_k):
    """Main conversion for a single radius threshold 'thr'. Returns
    (root, cyl, cyl_order): the cylinder list and, parallel to it, the branch
    order of each cylinder's end node (see compute_branch_order)."""
    root, parent, order = build_rooted_tree(len(xyz), edges, xyz)

    # list of children for each node (from the parent pointers)
    children = [[] for _ in range(len(xyz))]
    for v in order:
        p = parent[v]
        if p >= 0:
            children[p].append(v)

    # --- pruning: from the root down, only continue into nodes with radius >= thr
    keep = np.zeros(len(xyz), dtype=bool)
    stack = [root]; keep[root] = True
    while stack:
        n = stack.pop()
        for c in children[n]:
            if rad[c] >= thr:
                keep[c] = True
                stack.append(c)

    # remaining children (filtered through keep)
    kids = [[c for c in children[n] if keep[c]] for n in range(len(xyz))]

    # branch order per node on this pruned tree (root = 0, thickest child keeps
    # the parent's order, other children start order+1) - used below to tag
    # each cylinder, and later for AdQSM radius calibration.
    node_order = compute_branch_order(len(xyz), kids, rad, order)

    # --- resampling + cylinder creation (DFS from the root)
    # rule: a node is "kept as a model point" if it is a branch point (>=2 children),
    #       a branch tip (0 children), or we have already travelled far enough since
    #       the last kept point. The distance threshold is adaptive: it is derived
    #       from the LOCAL radius at the node being visited (thick branches get
    #       longer segments, thin twigs get finer ones).
    cyl = []                    # (start_idx, end_idx, radius, parent_cyl_id)
    cyl_order = []               # branch order of each cylinder's end node, parallel to cyl
    end_cyl = {}                # end_cyl[node] = id of the cylinder ending at this node
    no_resample = seg_len_min is None or seg_len_min <= 0

    # the stack carries: (node, anchor=last kept ancestor, distance travelled since anchor)
    st = [(c, root, float(np.linalg.norm(xyz[c] - xyz[root]))) for c in kids[root]]
    while st:
        n, anchor, dist = st.pop()
        kc = len(kids[n])
        if no_resample:
            retain = True
        else:
            target = local_seg_len(rad[n], seg_len_min, seg_len_max, seg_len_k)
            retain = (kc != 1) or (dist >= target)
        chord = float(np.linalg.norm(xyz[n] - xyz[anchor]))   # actual cylinder length
        if retain and chord >= MIN_CYL_LEN:
            # create cylinder  anchor -> n
            pid = end_cyl.get(anchor, 0)                 # 0 = root (no parent)
            cid = len(cyl) + 1                            # 1-based id
            r = 0.5 * (rad[anchor] + rad[n])
            cyl.append((anchor, n, r, pid))
            cyl_order.append(int(node_order[n]))
            end_cyl[n] = cid
            new_anchor, base = n, 0.0
        elif retain:
            # the branch point practically coincides with the anchor -> don't create a
            # zero-length cylinder, but attach children to the cylinder that ended at anchor
            end_cyl[n] = end_cyl.get(anchor, 0)
            new_anchor, base = n, 0.0
        else:
            new_anchor, base = anchor, dist              # skip this point
        for c in kids[n]:
            st.append((c, new_anchor, base + float(np.linalg.norm(xyz[c] - xyz[n]))))

    return root, cyl, cyl_order


def write_geom(path, xyz, cyl, root, recenter_xy):
    """Write the cylinders to geom.txt in the format read by buk.mac."""
    off = np.array([xyz[root, 0], xyz[root, 1], 0.0]) if recenter_xy else np.zeros(3)
    with open(path, "w") as f:
        # header row with column indices (0 = corner, then 1..9) - required by *TREAD
        f.write("0\t1\t2\t3\t4\t5\t6\t7\t8\t9\n")
        for k, (a, b, r, pid) in enumerate(cyl, start=1):
            s = xyz[a] - off                       # start point
            v = xyz[b] - xyz[a]                     # direction vector
            L = float(np.linalg.norm(v))
            ax = v / L if L > 0 else np.array([0.0, 0.0, 1.0])
            f.write("%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%d\n"
                    % (k, s[0], s[1], s[2], L, r, ax[0], ax[1], ax[2], pid))


def cylinder_metrics(xyz, cyl):
    """Return (lengths, radii) arrays for the given cylinders, computed the exact
    same way as the values written to geom.txt, so downstream metrics (volume,
    total length) match the model ANSYS will actually see."""
    if not cyl:
        return np.zeros(0), np.zeros(0)
    lengths = np.array([float(np.linalg.norm(xyz[b] - xyz[a])) for a, b, r, pid in cyl])
    radii = np.array([r for a, b, r, pid in cyl])
    return lengths, radii


def find_trunk_cylinders(cyl):
    """Identify the cylinders forming the main trunk: starting from the root,
    always follow the thickest child cylinder, stopping at the first real
    bifurcation (a cylinder with 2+ child cylinders) or at a branch tip.
    Returns a list of 1-based cylinder ids."""
    kids_of = {}
    for k, (a, b, r, pid) in enumerate(cyl, start=1):
        kids_of.setdefault(pid, []).append(k)

    trunk = []
    options = kids_of.get(0, [])
    while options:
        cur = max(options, key=lambda k: cyl[k - 1][2])   # thickest by radius
        trunk.append(cur)
        options = kids_of.get(cur, [])
        if len(options) >= 2:
            break
    return trunk


def report_volume(xyz, cyl, thr):
    """Print a volume/length verification summary for one radius-threshold
    variant, computed from the exact cylinders written to geom.txt (constant-
    radius cylinder volume = pi * radius^2 * length per segment)."""
    r_mm = int(round(thr * 1000))
    if not cyl:
        print("threshold %d mm : no cylinders" % r_mm)
        return

    lengths, radii = cylinder_metrics(xyz, cyl)
    volumes = np.pi * radii ** 2 * lengths
    total_len = float(lengths.sum())
    total_vol = float(volumes.sum())

    print("threshold %d mm : %d cylinders, length %.1f m, volume %.3f m^3 (%.0f L)"
          % (r_mm, len(cyl), total_len, total_vol, total_vol * 1000.0))

    # Branch order is not tracked by the model, so at minimum split trunk vs. branches.
    trunk_ids = find_trunk_cylinders(cyl)
    if trunk_ids:
        trunk_mask = np.zeros(len(cyl), dtype=bool)
        trunk_mask[np.array(trunk_ids) - 1] = True
        t_len, t_vol = float(lengths[trunk_mask].sum()), float(volumes[trunk_mask].sum())
        b_len, b_vol = total_len - t_len, total_vol - t_vol
        print("    trunk    : %d cylinders, length %.1f m, volume %.3f m^3 (%.0f L)"
              % (int(trunk_mask.sum()), t_len, t_vol, t_vol * 1000.0))
        print("    branches : %d cylinders, length %.1f m, volume %.3f m^3 (%.0f L)"
              % (len(cyl) - int(trunk_mask.sum()), b_len, b_vol, b_vol * 1000.0))
    else:
        print("    (no distinguishable trunk chain - reporting overall total only)")


def plot_model(xyz, cyl, root, recenter_xy, thr, out_path, show, save_png):
    """Draw the reduced beam model in 3D (fast, via Line3DCollection) for a quick
    visual check before importing into ANSYS."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    off = np.array([xyz[root, 0], xyz[root, 1], 0.0]) if recenter_xy else np.zeros(3)

    segments = np.empty((len(cyl), 2, 3))
    radii = np.empty(len(cyl))
    for i, (a, b, r, pid) in enumerate(cyl):
        segments[i, 0] = xyz[a] - off
        segments[i, 1] = xyz[b] - off
        radii[i] = r

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    lc = Line3DCollection(segments, array=radii, cmap="viridis", linewidths=1.2)
    ax.add_collection3d(lc)

    pts = segments.reshape(-1, 3)
    mins, maxs = pts.min(axis=0), pts.max(axis=0)
    ax.set_xlim(mins[0], maxs[0])
    ax.set_ylim(mins[1], maxs[1])
    ax.set_zlim(mins[2], maxs[2])
    ax.set_box_aspect(maxs - mins)     # true/equal aspect ratio from the data extents

    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title("Radius threshold: %d mm  |  cylinders: %d"
                 % (int(round(thr * 1000)), len(cyl)))

    cbar = fig.colorbar(lc, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label("Radius [m]")

    fig.tight_layout()

    if save_png:
        png_path = os.path.splitext(out_path)[0] + ".png"
        fig.savefig(png_path, dpi=200)
        print("  Saved plot:", png_path)

    if show:
        plt.show()
    else:
        plt.close(fig)


# =========================  RUN  =====================================
print("Reading:", INPUT_PLY)
xyz, rad, edges = read_ply(INPUT_PLY)
print("  vertices: %d, edges: %d" % (len(xyz), len(edges)))

xyz, rad, edges = merge_vertices(xyz, rad, edges, MERGE_DECIMALS)
ncomp, _ = connected_components(
    csr_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
               shape=(len(xyz), len(xyz))), directed=False)
print("  after merging: %d points, %d edges, connected components: %d" % (len(xyz), len(edges), ncomp))
xyz_raw = xyz.copy()   # AdTree geometry/radii before smoothing, kept as the calibration baseline

smooth_root = int(np.argmin(xyz[:, 2]))
xyz = smooth_centerline(xyz, edges, smooth_root, SMOOTH_ITERS, SMOOTH_ALPHA)
print("  centerline smoothing: %d passes, alpha=%.2f" % (SMOOTH_ITERS, SMOOTH_ALPHA))

if CALIBRATE_RADII:
    print("\nLoading AdQSM calibration data...")
    taper_heights, taper_diameters = parse_adqsm_taper_file(ADQSM_TAPER_FILE)
    trunk_radius_func = make_trunk_radius_func(taper_heights, taper_diameters, FIELD_DBH)
    adqsm_median_by_order = parse_adqsm_branch_file(ADQSM_BRANCH_FILE)
    print("  %s: %d height/diameter rows, %.1f-%.1f m"
          % (ADQSM_TAPER_FILE, len(taper_heights), taper_heights.min(), taper_heights.max()))
    if FIELD_DBH is not None:
        print("  taper curve rescaled so radius at 1.3 m = FIELD_DBH/2 = %.4f m" % (FIELD_DBH / 2.0))
    print("  %s: AdQSM median radius by order: %s"
          % (ADQSM_BRANCH_FILE,
             ", ".join("%d=%.4f m" % (o, r) for o, r in sorted(adqsm_median_by_order.items()))))

    raw_stats = raw_skeleton_stats(xyz_raw, rad, edges)
    print("  raw PLY skeleton baseline (all tree edges, AdTree radii):")
    print_volume_stats("raw skeleton (AdTree)", raw_stats)

print("\n%-12s %-12s %-12s %-12s" % ("threshold", "cylinders", "length [m]", "file"))
print("Volume verification below is computed from the exact cylinders written to "
      "each geom file (pi * radius^2 * length per cylinder).\n")
z_base = float(xyz[:, 2].min())   # tree base; DBH/height/taper are measured from here
multiple_thresholds = len(RADIUS_THRESHOLDS) > 1
for thr in RADIUS_THRESHOLDS:
    root, cyl, cyl_order = convert(xyz, rad, edges, thr, SEG_LEN_MIN, SEG_LEN_MAX, SEG_LEN_K)
    out = OUTPUT_PATTERN.format(r=int(round(thr * 1000))) if multiple_thresholds else OUTPUT_NAME

    # Height of the pruned model: z-range of the nodes actually used by these
    # cylinders. Unaffected by radius calibration (geometry doesn't change).
    node_ids = sorted({idx for a, b, r, pid in cyl for idx in (a, b)})
    height_m = float(xyz[node_ids, 2].max() - xyz[node_ids, 2].min()) if node_ids else None

    if CALIBRATE_RADII:
        orig_lengths, orig_radii = cylinder_metrics(xyz, cyl)
        orig_stats = volume_stats(orig_lengths, orig_radii, np.asarray(cyl_order))
        # DBH/taper of the UNCALIBRATED (raw AdTree) trunk, before radii are replaced.
        raw_dbh = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_LOWER)
        raw_d_upper = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_UPPER)
        raw_taper = ((raw_dbh - raw_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                     if raw_dbh is not None and raw_d_upper is not None else None)

        new_radii, factors = calibrate_cylinder_radii(
            xyz, cyl, cyl_order, trunk_radius_func, adqsm_median_by_order)
        cyl = [(a, b, float(new_radii[i]), pid) for i, (a, b, r, pid) in enumerate(cyl)]

    write_geom(out, xyz, cyl, root, RECENTER_XY)
    total_len = sum(float(np.linalg.norm(xyz[b] - xyz[a])) for a, b, _, _ in cyl)
    print("%-12s %-12d %-12.1f %-12s" % ("%d mm" % round(thr * 1000), len(cyl), total_len, out))
    report_volume(xyz, cyl, thr)   # uses the (possibly calibrated) radii above

    if CALIBRATE_RADII:
        print("  AdQSM calibration factors (order: AdQSM_median_radius / AdTree_median_radius):")
        for o in sorted(factors):
            print("    order %d : factor = %.3f" % (o, factors[o]))
        cal_lengths, cal_radii = cylinder_metrics(xyz, cyl)
        cal_stats = volume_stats(cal_lengths, cal_radii, np.asarray(cyl_order))

        # DBH/taper of the CALIBRATED trunk, using the now-replaced radii.
        cal_dbh = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_LOWER)
        cal_d_upper = stem_diameter_at_height(xyz, cyl, cyl_order, z_base, TAPER_H_UPPER)
        cal_taper = ((cal_dbh - cal_d_upper) * 100.0 / (TAPER_H_UPPER - TAPER_H_LOWER)
                     if cal_dbh is not None and cal_d_upper is not None else None)

        # ---- upsert both the uncalibrated and calibrated rows for this threshold ----
        upsert_result(RESULTS_CSV, TREE_NAME, "AdTree raw r%dmm" % round(thr * 1000),
                      orig_stats["total_vol"], orig_stats["trunk_vol"], orig_stats["branch_vol"], None,
                      raw_dbh, height_m, raw_taper)
        upsert_result(RESULTS_CSV, TREE_NAME, "AdTree calibrated r%dmm" % round(thr * 1000),
                      cal_stats["total_vol"], cal_stats["trunk_vol"], cal_stats["branch_vol"], None,
                      cal_dbh, height_m, cal_taper)

        print("  DBH (at %.1f m)   : raw AdTree = %s   |   calibrated = %s"
              % (TAPER_H_LOWER, _fmt_dbh(raw_dbh), _fmt_dbh(cal_dbh)))
        print("  Taper (%.1f-%.1f m): raw AdTree = %s   |   calibrated = %s"
              % (TAPER_H_LOWER, TAPER_H_UPPER, _fmt_taper(raw_taper), _fmt_taper(cal_taper)))
        print("  Height (pruned model): %s" % (("%.2f m" % height_m) if height_m is not None else "n/a"))

        print("  Volume comparison (a) raw skeleton vs. (b) processed/AdTree vs. (c) processed/calibrated:")
        print_volume_stats("(a) raw skeleton (AdTree)", raw_stats)
        print_volume_stats("(b) processed, AdTree radii", orig_stats)
        print_volume_stats("(c) processed, CALIBRATED", cal_stats)
        adqsm_ref = parse_adqsm_params_file(ADQSM_PARAMS_FILE)
        if adqsm_ref is not None:
            print_volume_stats("(d) AdQSM reference (TreesParams)", adqsm_ref)
            if adqsm_ref.get("height") is not None:
                print("      AdQSM TreeHeight: %.2f m" % adqsm_ref["height"])
    print()

    if SHOW_PLOT or SAVE_PLOT_PNG:
        plot_model(xyz, cyl, root, RECENTER_XY, thr, out, SHOW_PLOT, SAVE_PLOT_PNG)
