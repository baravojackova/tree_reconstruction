# -*- coding: utf-8 -*-
# =====================================================================
#  Single shared source of truth for every output directory this project
#  writes into (charts, CSVs, ANSYS geometry, calibration caches).
#
#  WHY THIS FILE EXISTS: scratch_output_path_audit.md (a prior session's
#  read-only audit) found ensure_plots_dir() independently re-implemented
#  in THREE different files (plot_volumes.py, calmethod_decision_summary.py,
#  reconstruction_method_decision_summary.py), each with its own local
#  PLOTS_DIR = "plots" - all three happened to agree today, but nothing
#  enforced that, so a future path change made in only one of them would
#  silently desync the other two. Centralising the directory NAMES and the
#  ensure_*() helpers here means there is exactly one place to change any
#  of them.
#
#  WHY THIS FILE IMPORTS NOTHING FROM THE REST OF THE PROJECT: every other
#  module in scripts/python may end up needing a path from here (plotting
#  scripts, the AdTree pipeline, the ANSYS export step, CSV-writing
#  scripts...) - if this file imported, say, compare_volumes.py (which
#  itself might one day want a path from here), that would create a
#  circular import the moment compare_volumes.py tried to import
#  something from this file. Keeping this module's only dependency on the
#  standard library (os) means it can be imported from literally anywhere
#  in this project with zero risk of an import cycle, now or later.
#
#  Every directory name below is relative, matching this project's
#  existing convention throughout scripts/python: every script assumes it
#  is LAUNCHED with scripts/python as the current working directory (see
#  scratch_output_path_audit.md's own opening note), so a bare relative
#  name here resolves under scripts/python exactly like the pre-existing
#  "plots"/"npz" literals already did.
# =====================================================================

import os

# ---- Directory NAMES (the literal strings) ---------------------------
# PLOTS_DIR: every PNG chart this project produces, of any kind, lives
# somewhere under here - either directly in a per-tree or all-trees
# subfolder (see ensure_tree_plots_dir()/ensure_all_plots_dir() below),
# never loose in PLOTS_DIR itself (that bare-root placement is exactly
# the "90 misplaced files" bug this reorganisation exists to fix).
PLOTS_DIR = "plots"

# CSV_DIR: every CSV this project WRITES as an output (decision summaries,
# coverage reports, cross-tree sensitivity tables, ...) lives here, FLAT -
# no per-tree subfolders. The tree name is already embedded in every
# per-tree CSV's own filename (e.g. "coverage_B21_S01.csv"), so a nested
# csv/<tree>/ layout would be redundant rather than clarifying - this is
# the "PNG goes under plots/, CSV goes under csv/" rule applied literally.
# volume_results.csv itself is the one exception: it is SHARED STATE (the
# project's master input/output table, not a derived report) and stays in
# scripts/python's own root, untouched by this constant.
CSV_DIR = "csv"

# NPZ_DIR: calibrated-geometry cache files (calib_*.npz) written by
# adtree_reconstruct_compare.py and read back by export_geom_ansys.py.
# Unchanged by this reorganisation - listed here only so every output
# location is named in one place, per this file's own purpose.
NPZ_DIR = "npz"

# ANSYS_GEOM_DIR: the actual geom_*.txt ANSYS INPUT geometry text files
# written by export_geom_ansys.py's write_geom() call. Before this file
# existed, that call wrote to a bare filename with NO directory component
# at all, landing wherever the process's current working directory
# happened to be - see ensure_ansys_geom_dir() below for the per-tree
# fix. Deliberately a SEPARATE top-level directory from PLOTS_DIR: this
# is ANSYS's own input geometry (a .txt file ANSYS *tread*s), not a chart,
# so it does not belong under plots/ even though its filename also
# happens to contain the substring "geom" (see the DIFFERENT, PNG-based
# "geom_*.png" diagnostic charts under plots/<tree>/geom/, produced by an
# unrelated function - keep these two "geom" things apart; they are not
# the same artifact).
ANSYS_GEOM_DIR = "ansys_geom"

# ALL_TREES_DIRNAME: the subfolder name (under PLOTS_DIR) for charts that
# POOL data across multiple trees (e.g. total_volume_by_tree.png,
# error_boxplot_*.png) - as opposed to a chart whose content concerns
# exactly one tree, which goes under plots/<that tree's name>/ instead.
ALL_TREES_DIRNAME = "all"


# ---- ensure_*() helpers -------------------------------------------------
# Every one of these ONLY creates a directory (os.makedirs(..., exist_ok=True)
# - safe to call every time, never errors if the directory already exists,
# and never touches/removes anything already inside it) and returns its
# path. None of them ever deletes or empties anything - see this project's
# own safety rules (do not clear plots/ or any subfolder, ever, for any
# reason) for why that matters here specifically.

def ensure_plots_dir():
    """Create PLOTS_DIR (the top-level "plots" folder) if it doesn't exist
    yet, and return its path. Mostly useful as the parent a per-tree/
    all-trees subfolder is created under - most callers want
    ensure_tree_plots_dir()/ensure_all_plots_dir() instead of writing
    directly into this bare directory (see PLOTS_DIR's own comment above
    for why nothing should land loose here any more)."""
    os.makedirs(PLOTS_DIR, exist_ok=True)
    return PLOTS_DIR


def ensure_tree_plots_dir(tree):
    """Create plots/<tree>/ if it doesn't exist yet (creating PLOTS_DIR
    itself too, if needed), and return its path. Use for any chart whose
    content concerns exactly one tree."""
    path = os.path.join(ensure_plots_dir(), tree)
    os.makedirs(path, exist_ok=True)
    return path


def ensure_tree_geom_dir(tree):
    """Create plots/<tree>/geom/ if it doesn't exist yet, and return its
    path. For the tree-SHAPE render PNGs written by
    tree_geom_utils.plot_model() specifically - NOT for the
    adtree_adqsm_radius_regression_perorder_*.png diagnostic charts (a
    different artifact from a different function, plot_radius_regression_
    per_order(), which stays directly under plots/<tree>/ - see paths.py's
    own ANSYS_GEOM_DIR comment for the same "two different things share
    the substring 'geom'" caveat)."""
    path = os.path.join(ensure_tree_plots_dir(tree), "geom")
    os.makedirs(path, exist_ok=True)
    return path


def ensure_all_plots_dir():
    """Create plots/all/ if it doesn't exist yet, and return its path.
    Use for any chart that pools/aggregates across more than one tree
    (e.g. total_volume_by_tree.png, error_boxplot_*.png,
    calmethod_decision_*.png) - as opposed to a chart scoped to one tree,
    which goes under plots/<tree>/ instead."""
    path = os.path.join(ensure_plots_dir(), ALL_TREES_DIRNAME)
    os.makedirs(path, exist_ok=True)
    return path


def ensure_csv_dir():
    """Create CSV_DIR (the flat "csv" folder) if it doesn't exist yet, and
    return its path. Use for every CSV this project WRITES as an output
    (not volume_results.csv itself - see CSV_DIR's own comment above)."""
    os.makedirs(CSV_DIR, exist_ok=True)
    return CSV_DIR


def ensure_ansys_geom_dir(tree):
    """Create ansys_geom/<tree>/ if it doesn't exist yet, and return its
    path. For the actual geom_*.txt ANSYS input geometry text file
    written by export_geom_ansys.py - see ANSYS_GEOM_DIR's own comment
    above for why this is a separate top-level directory from plots/."""
    path = os.path.join(ANSYS_GEOM_DIR, tree)
    os.makedirs(path, exist_ok=True)
    return path
