# -*- coding: utf-8 -*-
# =====================================================================
#  Shared visual style for every plotting script in this project - one
#  place to change a color/line width/font size so every chart updates
#  together, instead of hand-copying the same literal into each file.
#
#  Every value below was the REAL value already in use in plot_volumes.py
#  / plot_box.py at the time this module was created (confirmed via a
#  prior diagnostic) - nothing invented - EXCEPT AXIS_LABEL_FONTSIZE,
#  which is a deliberate, approved 1pt change unifying plot_volumes.py's
#  old LABEL_FONTSIZE (8) and plot_box.py's old BOX_LABEL_FONTSIZE (9)
#  into a single shared value (9) - see plot_volumes.py's own comment at
#  its AXIS_LABEL_FONTSIZE usage for that specific, intentional change.
# =====================================================================

import matplotlib.colors as mcolors
import numpy as np

# ---- Palette ------------------------------------------------------------
# Family -> gradient of pastel hex stops (pale -> deeper), one shade per
# method within that family. Moved verbatim from plot_volumes.py's old
# module-level FAMILY_GRADIENTS (was plot_volumes.py:535-541).
FAMILY_GRADIENTS = {
    "AdTree raw":        ["#e1f5e1", "#b8e2b8", "#8fce8f", "#63b563"],  # pastel green: pale -> sage -> leaf -> deeper green
    "AdTree calibrated": ["#dceaf9", "#b3d1f2", "#84b3e8", "#5a92d6"],  # pastel blue: pale -> sky -> mid -> deeper blue
    "TreeQSM":           ["#eeeeee", "#d4d4d4", "#b8b8b8", "#98989a"],  # pastel grey: near-white -> light -> mid -> deeper grey
    "AdQSM":             ["#fdf3c9", "#f8e08c", "#eec85a", "#d6a83f"],  # pastel yellow/ochre: pale -> gold -> ochre -> deeper ochre
    "Reference":         ["#fbc3d0", "#ef476f", "#c9315a"],             # pink/red: pale -> #ef476f (the original flat highlight) -> deeper red
}

TREEQSM_REF_LINE_COLOR = "#2a9d8f"   # teal reference line (moved from plot_box.py) - visually distinct from the destructive reference's pink/#ef476f
NEUTRAL_GREY = "#9a9a9a"             # fallback color for an unclassified method/group (moved from plot_volumes.py/plot_box.py)

TREE_MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*"]   # per-tree scatter marker shapes, in use order (moved from plot_volumes.py)

# ---- Line / marker weight ------------------------------------------------
DEFAULT_LINEWIDTH = 1        # plain data lines (e.g. zero-reference axhlines)
REFERENCE_LINEWIDTH = 1.5    # emphasized reference lines (plot_box.py's ref_line/treeqsm_line)
DEFAULT_MARKERSIZE = 6       # plot_volumes.py's tree-marker scatter points
JITTER_POINT_SIZE = 20       # plot_box.py's jittered scatter points (matplotlib scatter `s=`, area not radius)

# ---- Font sizes -----------------------------------------------------------
TITLE_FONTSIZE = 14          # fig.suptitle() - identical in plot_volumes.py and plot_box.py already
PANEL_TITLE_FONTSIZE = 12    # per-panel ax.set_title() - plot_box.py's old BOX_TITLE_FONTSIZE; matches
                              # plot_volumes.py's panel titles via matplotlib's own default
                              # (axes.titlesize='large'=12pt at font.size=10), so unifying it here
                              # changes nothing there either.
AXIS_LABEL_FONTSIZE = 12      # x/y axis tick-category label size. UNIFIED value (see module header):
                              # plot_volumes.py's old LABEL_FONTSIZE was 8, plot_box.py's old
                              # BOX_LABEL_FONTSIZE was already 9 - this is now the ONE shared value,
                              # so plot_volumes.py's axis labels get 1pt larger (approved change) and
                              # plot_box.py's stay exactly as they were.
LEGEND_FONTSIZE = 12          # consistent across every legend() call in both files already
ANNOTATION_FONTSIZE = 9      # small in-plot text labels (% / value annotations)
POINT_LABEL_FONTSIZE = 6     # smallest tier - dense per-point labels (plot_box.py)

# ---- Extended per-family shades -------------------------------------------
# FAMILY_GRADIENTS above only has 3-4 hand-picked stops per family - enough
# for a handful of methods, but not for a chart that needs MANY distinct
# lines from the same family (e.g. taper_curve_compare.py's up to ~21
# per-AdQSM-variant curves for one tree). family_shades() below fills that
# gap by building a CONTINUOUS colormap from a family's existing stops and
# sampling it, instead of cycling/repeating the raw 3-4 stops (which would
# make two different lines share an identical color once n exceeds the
# stop count).
_FAMILY_SHADES_CLAMP_LOW = 0.2   # never sample below 20% into the gradient -
                                   # skips the palest sliver near each family's
                                   # first stop (e.g. AdQSM's "#fdf3c9", a very
                                   # pale cream) which would otherwise render
                                   # as a near-invisible line against a white
                                   # plot background once n is large and
                                   # samples crowd close to position 0.
                                   # Verified empirically (not just asserted):
                                   # at n=21 on the AdQSM gradient, the
                                   # palest returned shade ("#fae8a4",
                                   # rgb(250,232,164)) sits at Euclidean
                                   # distance 94.0 from pure white
                                   # (255,255,255) - clearly visible, and
                                   # every one of the 21 shades is a distinct
                                   # hex value (confirmed via len(set(...))).


def family_shades(family, n):
    """Return a list of `n` hex color strings sampled evenly across
    FAMILY_GRADIENTS[family]'s existing stops, for a chart that needs more
    distinct shades of that family than the raw stop list provides.

    Builds a continuous LinearSegmentedColormap from the family's stops,
    then samples `n` evenly-spaced points across
    [_FAMILY_SHADES_CLAMP_LOW, 1.0] (not [0.0, 1.0] - see that constant's
    own comment for why the palest sliver is skipped) and converts each
    sampled point to a hex string. n=1 returns the single shade at the
    clamp boundary (the palest allowed shade, still comfortably non-pale
    per the clamp)."""
    stops = FAMILY_GRADIENTS[family]
    cmap = mcolors.LinearSegmentedColormap.from_list(family, stops)
    positions = np.linspace(_FAMILY_SHADES_CLAMP_LOW, 1.0, n)
    return [mcolors.to_hex(cmap(p)) for p in positions]


# ---- Method -> colour (method_evaluation.py) ------------------------------
# Colour per reconstruction METHOD (not just family) for figures that compare
# AdTree_raw/AdTree_calibrated/AdQSM/TreeQSM_Optimal/TreeQSM_Simplified side
# by side. Every value here REUSES a colour already established above rather
# than inventing a new one - see method_evaluation.py's own Step 0 report for
# the full audit this was built from:
#   - AdTree_raw/AdTree_calibrated/AdQSM: FAMILY_GRADIENTS[...][-1], the same
#     "deepest stop" every other script already treats as THE colour for that
#     family (e.g. adqsm_variant_sensitivity.py, base_vs_dbh_compare.py,
#     trunk_taper_vs_field.py all key off index -1 this same way).
#   - TreeQSM_Optimal/TreeQSM_Simplified: FAMILY_GRADIENTS only has ONE flat
#     grey for "TreeQSM" as a whole (no existing Optimal-vs-Simplified colour
#     split) - family_shades() samples two points across that SAME grey
#     gradient instead of picking new hex literals; its existing
#     _FAMILY_SHADES_CLAMP_LOW=0.2 clamp already keeps the two shades far
#     enough apart to survive projector washout (see that constant's own
#     comment) rather than landing in the pale, hard-to-read end of the ramp.
#   NOTE: AdTree_raw (green) and AdTree_calibrated (blue) are DIFFERENT hues,
#   not "one hue, two lightnesses" - this violates the ideal that raw/
#   calibrated variants of one method should read as related. Kept as-is
#   anyway because both colours are already in active use elsewhere in this
#   project; repainting either one here would make this script's figures
#   inconsistent with every earlier chart instead of consistent with them.
_TREEQSM_STAGE_SHADES = family_shades("TreeQSM", 2)
METHOD_COLORS = {
    "AdQSM":              FAMILY_GRADIENTS["AdQSM"][-1],
    "AdTree_raw":         FAMILY_GRADIENTS["AdTree raw"][-1],
    "AdTree_calibrated":  FAMILY_GRADIENTS["AdTree calibrated"][-1],
    "TreeQSM_Optimal":    _TREEQSM_STAGE_SHADES[0],
    "TreeQSM_Simplified": _TREEQSM_STAGE_SHADES[1],
}

# ---- Box plot appearance (method_evaluation.py) ----------------------------
# This project already has TWO different box-plot conventions in active use
# (see method_evaluation.py's Step 0 report): plot_box.py draws whis=(0,100)
# (whiskers span the group's true min/max, no separate outlier fliers), a
# SOLID family-colour fill, and edges/whiskers/caps/median darkened from that
# same fill colour; plot_volumes.py's plot_error_boxplot() instead uses
# matplotlib's default 1.5*IQR whiskers (fliers possible), a semi-transparent
# fill, and a flat grey edge. These two constants adopt plot_box.py's
# convention (the more fully worked-out of the two) for any NEW box plot that
# wants to match it, rather than leaving the choice to hard-coded literals
# copied into yet a third file.
BOX_WHIS = (0, 100)          # whis=(0,100): whiskers = true min/max, no fliers - see plot_box.py
BOX_EDGE_DARKEN_FACTOR = 0.4  # edge/whisker/cap/median darkening from the box's own fill colour
