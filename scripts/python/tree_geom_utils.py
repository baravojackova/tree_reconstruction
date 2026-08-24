# -*- coding: utf-8 -*-
# =====================================================================
#  Shared helper functions for the AdTree -> ANSYS geometry pipeline.
# ---------------------------------------------------------------------
#  This file has NO "RUN" section and does nothing by itself - it only
#  DEFINES functions, to be imported by the scripts that actually do work:
#
#    - adtree_reconstruct_compare.py : reads the .ply skeleton, calibrates
#      it against AdQSM, prints/saves comparison stats, and saves the
#      final calibrated geometry to a "calib_*.npz" file per threshold.
#    - export_geom_ansys.py : loads one "calib_*.npz" file and writes the
#      final geom_*.txt that ANSYS (buk.mac) actually reads.
#
#  Every function here is "pure" in the sense that it only depends on its
#  own arguments (no module-level PARAMETERS block, no printing tied to a
#  particular run) - that's what makes it safe to import from two different
#  scripts without dragging along unrelated state.
#
#  Dependencies: numpy, scipy   (matplotlib is imported lazily, only inside
#  plot_model(), so importing this file doesn't require it unless you
#  actually call that one function).
# =====================================================================

import os
import re
import csv

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import breadth_first_order


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
    - only the FIRST block of numeric rows is used, per AdQSM convention.

    Isolated "spike" rows (a single row with an implausibly large diameter
    compared to both its neighbours) are dropped before returning - see
    _reject_taper_spikes() below for why this is necessary and how it
    decides what counts as a spike."""
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
    heights = np.asarray(heights)[idx]
    diameters = np.asarray(diameters)[idx]
    return _reject_taper_spikes(heights, diameters, path)


def _reject_taper_spikes(heights, diameters, path, factor=2.0):
    """Drop isolated "spike" rows from an AdQSM taper curve before it's used
    for radius calibration.

    WHY THIS EXISTS: a tree trunk's diameter should change fairly smoothly
    with height - it should never jump to several times the diameter of the
    rows immediately above AND below it. In practice, AdQSM's own taper.txt
    export CAN contain an isolated garbage row (observed for real, in
    data/IND07_083/05/taper.txt: height 26.6 m reports diameter 4.43 m,
    sandwiched between 0.55 m and 0.49 m at the neighbouring heights -
    clearly an export/fitting artifact, not a real 4+ metre-thick trunk).

    If a spike like that is left in, make_trunk_radius_func()'s linear
    interpolation draws a straight line up to it and back down, so ANY
    AdTree trunk cylinder whose height happens to fall near the spike gets
    assigned a hugely inflated radius - which can multiply the CALIBRATED
    stem volume several-fold, even though DBH and taper_cm_per_m (measured
    at 1.3 m / 10.0 m, usually nowhere near a spike further up the trunk)
    come out looking completely normal. That mismatch - "calibrated stem
    volume way too high, but DBH/taper look fine" - is exactly the symptom
    this function prevents.

    An interior row i (never the very first or last row - there's no
    "neighbour on both sides" to compare those against) is dropped only if
    its diameter is more than `factor` times BOTH of its immediate
    neighbours' diameters, i.e. it stands out as a spike relative to what's
    on EITHER side, not merely part of an ordinary gentle taper (a real
    trunk can widen slightly lower down before tapering - e.g. a root
    flare - so this check deliberately only fires on an isolated,
    much-larger-than-both-neighbours point, never on a normal small
    increase). Never silent: prints exactly which row(s) were dropped,
    since losing a row changes the calibration result and that should be
    visible, not a silent correction.
    """
    keep = np.ones(len(diameters), dtype=bool)
    for i in range(1, len(diameters) - 1):
        left, right = diameters[i - 1], diameters[i + 1]
        neighbor_max = max(left, right)
        if neighbor_max > 0 and diameters[i] > factor * neighbor_max:
            keep[i] = False
            print("  WARNING: dropping implausible taper.txt row at height %.2f m "
                  "(diameter %.4f m vs. neighbouring rows %.4f m / %.4f m) from %s "
                  "- looks like an AdQSM export artifact, not real trunk data."
                  % (heights[i], diameters[i], left, right, path))
    return heights[keep], diameters[keep]


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


def _read_adqsm_branch_header(path):
    """Return the BranchStructure.txt header row (list of column names) if
    the file has one (a line whose first TAB-separated field is literally
    "order"), else None. Used by the two diagnostic functions below so they
    both agree on what each column means, instead of guessing twice."""
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if parts and parts[0] == "order":
                return parts
    return None


def _find_adqsm_column(header_cols, keywords):
    """Return the index of the first header column whose name contains ANY
    of `keywords` (case-insensitive), or None if there is no header or no
    match. E.g. _find_adqsm_column(header, ["length"]) -> index of the
    "length(m)" column, whatever its exact spelling."""
    if not header_cols:
        return None
    for i, name in enumerate(header_cols):
        if any(k in name.lower() for k in keywords):
            return i
    return None


def print_adqsm_branch_file_sample(path, n=10):
    """DIAGNOSTIC: print the first `n` valid data rows of an AdQSM
    BranchStructure.txt, split by TAB, one line per column, showing the
    column's index AND its name (read from the file's own header row, if
    present) AND its value. Run this before trusting any calculation based
    on the file, so you can see its real layout instead of assuming it.
    """
    header_cols = _read_adqsm_branch_header(path)
    if header_cols:
        print("  Header row found: %s" % "\t".join(header_cols))
    else:
        print("  No header row found ('order\\t...') - columns are unlabeled below.")

    shown = 0
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if not parts or parts[0] == "order" or not parts[0]:
                continue           # skip the header row and blank lines
            try:
                int(parts[0])       # data rows start with an integer branch order
            except ValueError:
                continue
            print("  --- row %d (%d columns) ---" % (shown, len(parts)))
            for i, value in enumerate(parts):
                col_name = header_cols[i] if header_cols and i < len(header_cols) else "col%d" % i
                print("    [%d] %-16s = %s" % (i, col_name, value))
            shown += 1
            if shown >= n:
                break


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
        # DBH (trunk diameter at breast height), exported by AdQSM as "TreeDBH" in cm.
        tree_dbh = extract("TreeDBH", line)/100
        return dict(n=int(branches_num), total_len=trunk_len + branch_len, total_vol=tree_vol,
                    trunk_n=0, trunk_len=trunk_len, trunk_vol=trunk_vol,
                    branch_n=0, branch_len=branch_len, branch_vol=branch_vol,
                    height=tree_height, dbh=tree_dbh)

    return None


def report_adqsm_thin_branch(path, cut_cm=10.0, params_file=None):
    """AdQSM equivalent of report_thin_branch_volume() (see below), built
    directly from BranchStructure.txt.

    IMPORTANT: this file gives only ONE diameter per row (no proximal/distal
    pair), so AdQSM.pdf's own Smalian's-formula volume can't be reproduced
    here. If a "length" column exists, this function instead approximates
    each row's volume as a CONSTANT-diameter cylinder (pi*(d/2)^2*length) -
    a rough estimate that OVER-counts volume for branches that taper a lot
    along their length (thin, elongated ones especially). To let you judge
    how large that bias is, the cylinder-approximation totals are printed
    next to AdQSM's own official TrunkVolume/BranchVolume from
    TreesParams.txt (via params_file), if given.

    The file's own "volume(...)" column is intentionally NOT used: its sum
    is off by several orders of magnitude from TreesParams.txt's official
    BranchVolume/TrunkVolume (not a simple unit-conversion factor), so it
    does not appear to be a trustworthy per-branch wood volume.

    If NO length column can be found at all, falls back to a proxy: just
    the COUNT and PERCENTAGE of branches thinner than cut_cm - clearly
    labelled as a proxy, since no volume can honestly be computed from
    diameter alone.
    """
    header_cols = _read_adqsm_branch_header(path)
    idx_diam = _find_adqsm_column(header_cols, ["diameter"])
    if idx_diam is None:
        idx_diam = 2   # fallback: column 2 is "diameter(m)" in this dataset
    idx_length = _find_adqsm_column(header_cols, ["length"])

    orders, diam_m, length_m = [], [], []
    with open(path, "r", encoding="latin-1") as f:
        for line in f:
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 7:
                continue
            try:
                order = int(parts[0])
                diameter = float(parts[idx_diam])
            except (ValueError, IndexError):
                continue
            orders.append(order)
            diam_m.append(diameter)
            if idx_length is not None:
                try:
                    length_m.append(float(parts[idx_length]))
                except (ValueError, IndexError):
                    length_m.append(None)
            else:
                length_m.append(None)

    orders = np.asarray(orders)
    diam_m = np.asarray(diam_m)
    diam_cm = diam_m * 100.0
    is_stem = orders == 0
    is_branch = orders >= 1
    keep = diam_cm >= cut_cm       # True = at/above the cut-off (kept by the reference)

    print("\n--- AdQSM BranchStructure.txt, cut-off %.0f cm diameter ---" % cut_cm)
    have_length = idx_length is not None and all(v is not None for v in length_m)

    if have_length:
        length_arr = np.asarray(length_m)
        volumes = np.pi * (diam_m / 2.0) ** 2 * length_arr   # constant-diameter cylinder approx.
        diam_name = header_cols[idx_diam] if header_cols else "diameter"
        length_name = header_cols[idx_length] if header_cols else "length"
        print("  Using columns: diameter='%s' (col %d), length='%s' (col %d)"
              % (diam_name, idx_diam, length_name, idx_length))
        print("  NOTE: single-diameter cylinder approximation (no proximal/distal pair")
        print("  available for Smalian's formula) - this OVER-estimates volume for")
        print("  tapering branches. See the official-totals comparison below.")

        stem_vol = float(volumes[is_stem].sum())
        branch_vol = float(volumes[is_branch].sum())
        stem_vol_kept = float(volumes[is_stem & keep].sum())
        branch_vol_kept = float(volumes[is_branch & keep].sum())

        for label, v_total, v_kept in (("Stem", stem_vol, stem_vol_kept),
                                        ("Branch", branch_vol, branch_vol_kept)):
            v_removed = v_total - v_kept
            pct_removed = (v_removed / v_total * 100.0) if v_total else 0.0
            print("  %-6s: total %.3f m3 (cylinder approx) | kept %.3f m3 | "
                  "removed %.3f m3 (%.1f %%)" % (label, v_total, v_kept, v_removed, pct_removed))

        if params_file is not None:
            official = parse_adqsm_params_file(params_file)
            if official is not None:
                print("  For comparison, AdQSM's OWN official totals (TreesParams.txt):")
                print("    TrunkVolume  = %7.3f m3   (cylinder approx above: %7.3f m3, %.1fx)"
                      % (official["trunk_vol"], stem_vol,
                         (stem_vol / official["trunk_vol"]) if official["trunk_vol"] else float("nan")))
                print("    BranchVolume = %7.3f m3   (cylinder approx above: %7.3f m3, %.1fx)"
                      % (official["branch_vol"], branch_vol,
                         (branch_vol / official["branch_vol"]) if official["branch_vol"] else float("nan")))
    else:
        # No usable length column - deliberately do NOT invent a volume;
        # report only what the file actually supports: branch count/share.
        print("  No usable length column found in BranchStructure.txt - volume can't")
        print("  be computed. PROXY METRIC ONLY (branch count/share, NOT volume):")
        for label, mask in (("Stem", is_stem), ("Branch", is_branch)):
            n_total = int(mask.sum())
            n_removed = int((mask & ~keep).sum())
            pct_removed = (n_removed / n_total * 100.0) if n_total else 0.0
            print("    %-6s: %5d total, %5d thinner than %.0f cm (%.1f %%)"
                  % (label, n_total, n_removed, cut_cm, pct_removed))


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


def report_thin_branch_volume(lengths, radii, order_arr, cut_cm=10.0, source_label="AdTree calibrated"):
    """Diagnostic: how much cylinder volume sits in cylinders THINNER than
    `cut_cm` diameter, split into stem (order 0) and branches (order >= 1).
    Mirrors "Cylinders, cut-off 10 cm" in runsken.m (section 17b) so the two
    printouts are directly comparable: the de Tanago field reference only
    measured branches down to a 10 cm taper diameter, so any volume below
    that cut-off could never have been checked against it anyway.

    lengths/radii/order_arr are arrays parallel to the final cylinders (same
    ones cylinder_metrics()/volume_stats() use). Returns a dict with the
    kept (>= cut_cm) volumes, so the RUN section can optionally write them
    into RESULTS_CSV as a second, filtered row.

    `source_label` names WHICH cylinder set this call is reporting on, in
    the printed header only (e.g. "AdTree raw" vs. "AdTree calibrated") -
    this function is called TWICE per threshold in adtree_reconstruct_compare.py
    (once on the raw/uncalibrated cylinders, once on the calibrated ones), so
    without this label the two printouts would be visually indistinguishable
    in the console output. Defaults to "AdTree calibrated" to match every
    call site that existed before this parameter was added.
    """
    lengths = np.asarray(lengths)
    radii = np.asarray(radii)
    order_arr = np.asarray(order_arr)

    volumes = np.pi * radii ** 2 * lengths      # per-cylinder volume [m^3]
    diam_cm = 2.0 * radii * 100.0               # per-cylinder diameter [cm]
    keep = diam_cm >= cut_cm                    # True = at/above the cut-off (kept)

    is_stem = order_arr == 0
    is_branch = order_arr >= 1

    n_total = len(radii)
    n_kept = int(keep.sum())
    vol_total = float(volumes.sum())
    vol_kept = float(volumes[keep].sum())
    vol_removed = vol_total - vol_kept

    print("\n--- Cylinders, cut-off %.0f cm (%s) ---" % (cut_cm, source_label))
    print("Cylinders total  : %d" % n_total)
    print("Cylinders kept   : %d (%.1f %%)" % (n_kept, (n_kept / n_total * 100.0) if n_total else 0.0))
    print("Volume total     : %.3f m3" % vol_total)
    print("Volume kept      : %.3f m3" % vol_kept)
    print("Volume removed   : %.3f m3 (%.1f %%)"
          % (vol_removed, (vol_removed / vol_total * 100.0) if vol_total else 0.0))

    # n_cyl_kept: count of cylinders that passed the diameter filter (Task B) -
    # reuses n_kept, which is already int(keep.sum()) computed above, so the
    # count and the vol_kept/len_kept values below all come from the exact
    # same "keep" mask instead of being recomputed separately.
    result = dict(total_vol=vol_total, total_vol_kept=vol_kept, n_cyl_kept=n_kept)
    for label, key, mask in (("Stem", "trunk", is_stem), ("Branch", "branch", is_branch)):
        v_total = float(volumes[mask].sum())
        v_kept = float(volumes[mask & keep].sum())
        v_removed = v_total - v_kept
        pct_removed = (v_removed / v_total * 100.0) if v_total else 0.0
        print("%-6s volume kept  : %.3f m3  (removed %.3f m3, %.1f %%)"
              % (label, v_kept, v_removed, pct_removed))
        result["%s_vol" % key] = v_total
        result["%s_vol_kept" % key] = v_kept

        # Same idea as the volume block above, but summing cylinder LENGTH
        # instead of volume. `lengths` is already an input parameter of this
        # function (used to build `volumes` above), so no extra data is
        # needed - we just also sum it directly, using the same `mask`
        # (stem vs. branch) and the same `keep` (>= cut_cm diameter) filter.
        len_total = float(lengths[mask].sum())      # total length of this group (stem or branch), before filtering
        len_kept = float(lengths[mask & keep].sum())  # length remaining after removing thin (< cut_cm) cylinders
        len_removed = len_total - len_kept
        pct_len_removed = (len_removed / len_total * 100.0) if len_total else 0.0
        print("%-6s length kept  : %.3f m   (removed %.3f m, %.1f %%)"
              % (label, len_kept, len_removed, pct_len_removed))
        result["%s_len" % key] = len_total
        result["%s_len_kept" % key] = len_kept
    return result


def upsert_result(csv_path, tree, method, total, stem, branch, std, dbh=None, height=None, taper=None,
                   trunk_len=None, branch_len=None, branch_filter="none", n_cylinders=None):
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
    header comment for why this distinction matters)."""
    # n_cylinders is the LAST column (Task A), added after branch_filter so
    # every existing column keeps its position - old rows/readers relying
    # on column position elsewhere are unaffected.
    header = ["tree", "method", "total_m3", "stem_m3", "branch_m3", "std_m3",
              "dbh_m", "height_m", "taper_cm_per_m", "trunk_len_m", "branch_len_m",
              "branch_filter", "n_cylinders"]

    def fmt(x):
        return "" if x is None else "%.6f" % x

    def fmt_int(x):
        # Cylinder count is always a whole number (not a measured float
        # like the other columns), so use "%d" instead of fmt()'s "%.6f" -
        # still blank ("") when n_cylinders is None, same missing-value
        # convention as every other optional column here.
        return "" if x is None else "%d" % int(x)

    rows = []
    if os.path.exists(csv_path):
        with open(csv_path, "r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))

    rows = [r for r in rows if not (r["tree"] == tree and r["method"] == method)]
    rows.append({"tree": tree, "method": method, "total_m3": fmt(total),
                 "stem_m3": fmt(stem), "branch_m3": fmt(branch), "std_m3": fmt(std),
                 "dbh_m": fmt(dbh), "height_m": fmt(height), "taper_cm_per_m": fmt(taper),
                 "trunk_len_m": fmt(trunk_len), "branch_len_m": fmt(branch_len),
                 "branch_filter": branch_filter, "n_cylinders": fmt_int(n_cylinders)})

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


def convert(xyz, rad, edges, thr, seg_len_min, seg_len_max, seg_len_k, min_cyl_len):
    """Main conversion for a single radius threshold 'thr'. Returns
    (root, cyl, cyl_order): the cylinder list and, parallel to it, the branch
    order of each cylinder's end node (see compute_branch_order).

    `min_cyl_len` is MIN_CYL_LEN from the calling script's PARAMETERS block
    (passed in explicitly, since this module has no parameters of its own)."""
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
        if retain and chord >= min_cyl_len:
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


def write_geom(path, xyz, cyl, root, recenter_xy, cyl_order):
    """Write the cylinders to geom.txt in the format read by buk.mac.

    NEW PARAMETER cyl_order: a list/array of ints, one branch-order value per
    cylinder in `cyl`, IN THE SAME ORDER as `cyl` (cyl_order[i] belongs to
    cyl[i]). This is the same list produced by compute_branch_order() /
    convert() elsewhere in this file. WHY it's added: the MATLAB script
    myfun.m (function result_ansys) writes an 11th column with the branch
    order of each cylinder, and geom_*.txt needs to match that format so
    downstream ANSYS/MATLAB tooling that expects 11 columns keeps working.
    """
    # WHY this check: cyl_order must line up 1-to-1 with cyl (same length,
    # same order). If some caller passes a mismatched list (e.g. from a
    # different run, or forgets to update it after filtering `cyl`), the
    # bug would otherwise be silent - each cylinder would silently get the
    # WRONG branch order instead of an obvious crash. Failing loudly here
    # makes that mistake easy to spot immediately.
    if len(cyl_order) != len(cyl):
        raise ValueError(
            "write_geom: cyl_order has %d entries but cyl has %d - they must "
            "be the same length and in the same order." % (len(cyl_order), len(cyl)))

    off = np.array([xyz[root, 0], xyz[root, 1], 0.0]) if recenter_xy else np.zeros(3)
    with open(path, "w") as f:
        # header row with column indices (0 = corner, then 1..9, now also 10
        # for the new branch-order column) - required by *TREAD
        f.write("0\t1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n")
        for k, (a, b, r, pid) in enumerate(cyl, start=1):
            s = xyz[a] - off                       # start point
            v = xyz[b] - xyz[a]                     # direction vector
            L = float(np.linalg.norm(v))
            ax = v / L if L > 0 else np.array([0.0, 0.0, 1.0])
            # k is 1-based (enumerate(..., start=1)) but cyl_order is a plain
            # 0-based Python list, so the cylinder written on loop iteration
            # k corresponds to cyl_order[k - 1], NOT cyl_order[k].
            order_k = cyl_order[k - 1]
            f.write("%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%d\t%d\n"
                    % (k, s[0], s[1], s[2], L, r, ax[0], ax[1], ax[2], pid, order_k))


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
