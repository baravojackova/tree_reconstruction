# -*- coding: utf-8 -*-
# =====================================================================
#  Draw and save charts that summarize volume_results.csv (the shared
#  master results table produced by ply_to_geom.py, qsm_volume_mean.py,
#  reference_volume.py, ... and printed in detail by compare_volumes.py).
# ---------------------------------------------------------------------
#  Like compare_volumes.py, this script works in TWO separate comparison
#  MODES (see that file's header comment for the full "why"):
#    - branch_filter == "10cm": vs. the destructive field reference
#      (REFERENCE_METHOD) - the fair, apples-to-apples accuracy check.
#    - branch_filter == "none": methods vs. EACH OTHER, using AdQSM
#      (REFERENCE_METHOD_NONE) as a common yardstick, since the destructive
#      reference can never appear in this mode.
#  Charts that compare "vs. a reference" are drawn ONCE PER MODE (so you get
#  two separate PNGs, one per mode, never mixed together in one chart).
#
#  This script makes the following PNG charts every time you run it:
#
#   a) total_volume_by_tree.png
#        Bar chart of total_m3 per method, with one GROUP of bars per tree.
#        The reference method's bar is drawn in a different colour and
#        labelled "(reference)" in the legend, so it's easy to spot. Only
#        drawn for the "10cm" mode (see plot_total_volume_by_tree()'s
#        comment for why). Always drawn, no matter how many trees are in the CSV.
#
#   b) error_boxplot_10cm.png / error_boxplot_none.png
#        Box plot of the percentage error of each method's total_m3 vs. that
#        mode's reference, one box per method, built from the SAME trees
#        used for that method's box (percent error computed the same way as
#        pct_diff() in compare_volumes.py - imported from there, not
#        re-implemented). Needs at least 2 trees (in THAT mode) to say
#        anything meaningful (comparing across trees is the whole point) -
#        with fewer, the RUN section below skips that mode's chart and
#        prints why, instead of drawing something misleading.
#
#   c) error_metrics_bar_10cm.png / error_metrics_bar_none.png
#        Bar chart of Bias / MAE / RMSE per method, using the exact same
#        calculation as compare_volumes.py's error_metrics() table
#        (via compute_error_metrics(), imported from compare_volumes.py -
#        so the numbers here can never drift out of sync with the printed
#        table). Same "needs >=2 trees (in that mode)" rule as (b).
#
#   d) tree_overview_<tree>_10cm.png / tree_overview_<tree>_none.png (one
#      pair PER tree in the CSV)
#        A single figure with a 2x3 grid of bar charts for ONE tree AND ONE
#        mode: total volume, DBH, height, taper, trunk length, branch length
#        - one bar per method in each subplot. Unlike (b)/(c), this makes
#        sense even with just ONE tree in the CSV (it doesn't compare across
#        trees, only across methods for the SAME tree), so both mode's PNGs
#        are always drawn for every tree found in the CSV.
#
#  All PNGs are written into a "plots" subfolder next to this script
#  (created automatically if it doesn't exist yet). The path of each saved
#  file is printed to the console.
#
#  Colour scheme: every chart above shares ONE {method: colour} mapping per
#  mode (see build_method_color_map()), so the same method is always the
#  same colour across every chart in that mode - the reference method always
#  gets a fixed highlight colour, every other method is spread across a
#  smooth green -> gray -> yellow -> orange gradient.
#
#  Dependencies: matplotlib   (install: pip install matplotlib)
#  This script also IMPORTS a few things from compare_volumes.py, which
#  must live in the same folder:
#     RESULTS_CSV, REFERENCE_METHOD, REFERENCE_METHOD_NONE, load_results,
#     pct_diff, compute_error_metrics, filter_by_branch_filter
# =====================================================================

import os
import random   # used in plot_error_boxplot() to jitter the individual-tree scatter points sideways
import re        # used by shorten_method_label() to parse the CSV's full method strings
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors     # used in _darken_color() and build_method_color_map()
import matplotlib.patches as mpatches   # used to build the shared legend in plot_tree_overview()

from compare_volumes import (
    RESULTS_CSV,
    REFERENCE_METHOD,
    REFERENCE_METHOD_NONE,
    load_results,
    pct_diff,
    compute_error_metrics,
    filter_by_branch_filter,
)

# =====================  PARAMETERS  ==================================
# Folder (relative to this script's working directory) where the PNG
# charts are saved. Created automatically if it doesn't exist.
PLOTS_DIR = "plots"
# =====================================================================


def ensure_plots_dir():
    """Create PLOTS_DIR if it doesn't exist yet, and return its path."""
    if not os.path.isdir(PLOTS_DIR):
        os.makedirs(PLOTS_DIR)   # makedirs (not mkdir) also creates parent folders if needed
    return PLOTS_DIR


def save_and_report(fig, filename):
    """Save one matplotlib figure into PLOTS_DIR, close it (frees memory),
    and print the full path so you know where to look for it."""
    out_path = os.path.join(ensure_plots_dir(), filename)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", out_path)


# ----------------------------------------------------------------------
# Short display labels, used by EVERY chart in this file wherever a
# method's name is actually DRAWN on the figure (x-tick labels, legend
# text, titles that name a specific method). This is DISPLAY-ONLY - see
# shorten_method_label()'s own docstring - volume_results.csv, color_map
# dict keys, and every lookup/filter/comparison elsewhere in this file
# keep using the FULL, untouched method string throughout.
# ----------------------------------------------------------------------

# TreeQSM's "group" qualifier (the second, comma-separated part of
# "TreeQSM mine (<run>, <group>)") doesn't shorten by one consistent rule -
# "Optimal (single)" drops its parenthetical entirely, while "Simplified
# (no islands)" keeps it (as an underscore-joined suffix, so it stays
# visually distinguishable from plain "Simplified" rather than collapsing
# onto it). Hard-coded per the groups actually seen in volume_results.csv
# today; extend this table if a new TreeQSM group string appears.
_TREEQSM_GROUP_SHORT = {
    "Optimal (single)": "Optimal",
    "Optimal": "Optimal",   # current volume_results.csv uses this form (no "(single)" suffix)
    "Simplified": "Simplified",
    "Simplified (no islands)": "Simplified_noislands",
    "Filtered <10cm": "Filtered10cm",
}

# Known "[calmethod=...]" tag VALUES -> short suffix word (used by
# shorten_method_label() below). calref=... tags don't need this table -
# their value (e.g. "min5mm") is already short enough to use as-is.
_CALMETHOD_SHORT = {
    "regression-perorder": "regperorder",
    "regression": "regression",
    "interpolation-perorder": "interpperorder",
}


def _calibration_tag_suffix(tag):
    """Turn a "calref=X" or "calmethod=X" tag (captured WITHOUT the
    surrounding brackets) into a short suffix to append to an
    "AT_Calib_..." label - e.g. "calref=min5mm" -> "_calref-min5mm",
    "calmethod=regression-perorder" -> "_regperorder". Returns "" if `tag`
    is None (no tag present in this method string at all - the old,
    pre-multi-method rows had none)."""
    if tag is None:
        return ""
    if tag.startswith("calref="):
        return "_calref-%s" % tag[len("calref="):]
    if tag.startswith("calmethod="):
        value = tag[len("calmethod="):]
        short = _CALMETHOD_SHORT.get(value, value.replace("-", ""))
        # General rule (not a one-off replace): whenever "reg" and
        # "perorder" appear together in the short form (e.g. "regperorder"
        # from the dict above, or "regressionfooperorder" from a future
        # calmethod value not yet added to the dict), drop the "perorder"
        # part - it's redundant with per-order regression already being
        # implied by context, and keeps the label short without needing a
        # new dict entry for every future reg*-perorder variant.
        short = re.sub(r'reg\w*?perorder', 'reg', short, count=1)
        return "_%s" % short
    return ""   # unrecognized tag shape - shouldn't happen, but never crash over it


def shorten_method_label(method):
    """Map a full volume_results.csv method string to a short display
    label, for legends/titles/axis labels ONLY - never used to look up or
    filter rows (matching against the CSV must always use the full,
    untouched method string; color_map keys, reference_method comparisons,
    etc. all keep using `method` as-is, never this function's output).

    Handles the method-string "families" as they actually appear in
    volume_results.csv today (see each branch below for the exact
    pattern); anything that doesn't match a known pattern is returned
    UNCHANGED - a long label drawn once is far less harmful than silently
    mislabelling or dropping a method from a chart.
    """
    # "AdTree calibrated r{N}mm (AdQSM {variant}) [{tag}]? (>=10cm only)?_seg{min}-{max}-k{k}"
    # -> "AT_Calib_{N}_{variant}_{min}/{max}/{k}" (+ a short tag suffix, if
    # a "[calref=...]"/"[calmethod=...]" tag is present). The "(>=10cm
    # only)" marker is dropped entirely (not folded into the short name) -
    # branch_filter is already shown elsewhere (e.g. tree_overview's own
    # suptitle), so repeating it here would be redundant.
    m = re.match(
        r'^AdTree calibrated r(\d+)mm \(AdQSM ([^)]+)\)'
        r'(?: \[(calref=[^\]]+|calmethod=[^\]]+)\])?'
        r'(?: \(>=\d+cm only\))?'
        r'_seg(\d+)-(\d+)-k(\d+)$',
        method)
    if m:
        n_mm, variant, tag, seg_min, seg_max, seg_k = m.groups()
        return "AT_Calib_%s_%s_%s/%s/%s" % (n_mm, variant, seg_min, seg_max, seg_k) \
            + _calibration_tag_suffix(tag)

    # "AdTree raw r{N}mm (>=10cm only)?_seg{min}-{max}-k{k}"
    # -> "AT_Raw_{N}_{min}/{max}/{k}" - no variant (raw isn't calibrated
    # against any particular AdQSM variant).
    m = re.match(r'^AdTree raw r(\d+)mm(?: \(>=\d+cm only\))?_seg(\d+)-(\d+)-k(\d+)$', method)
    if m:
        n_mm, seg_min, seg_max, seg_k = m.groups()
        return "AT_Raw_%s_%s/%s/%s" % (n_mm, seg_min, seg_max, seg_k)

    # "TreeQSM mine ({run}, {group})" -> "TQ_{run}_{group_short}". `group`
    # can itself contain parentheses (e.g. "Optimal (single)"), so this is
    # parsed by stripping the outer "TreeQSM mine (...)" wrapper and
    # splitting on the FIRST ", " rather than with one regex.
    if method.startswith("TreeQSM mine (") and method.endswith(")"):
        inner = method[len("TreeQSM mine ("):-1]
        if ", " in inner:
            run, group = inner.split(", ", 1)
            group_short = _TREEQSM_GROUP_SHORT.get(group)
            if group_short is not None:
                return "TQ_%s_%s" % (run, group_short)

    # "TreeQSM de Tanago (mean)" / "TreeQSM de Tanago (mean, Filtered<10cm)"
    # -> "TQ_ref" / "TQ_ref_Filtered10cm" (exact matches - only these two
    # variants exist today).
    if method == "TreeQSM de Tanago (mean)":
        return "TQ_ref"
    if method == "TreeQSM de Tanago (mean, Filtered<10cm)":
        return "TQ_ref_Filtered10cm"

    # "AdQSM (TreesParams) (AdQSM {variant})" -> "AQ_Params_{variant}"
    m = re.match(r'^AdQSM \(TreesParams\) \(AdQSM ([^)]+)\)$', method)
    if m:
        return "AQ_Params_%s" % m.group(1)

    # "AdQSM (BranchStructure, cyl. approx.) (AdQSM {variant})"
    # -> "AQ_BranchApprox_{variant}"
    m = re.match(r'^AdQSM \(BranchStructure, cyl\. approx\.\) \(AdQSM ([^)]+)\)$', method)
    if m:
        return "AQ_BranchApprox_%s" % m.group(1)

    # "Reference (destructive)" -> "Ref"
    if method == "Reference (destructive)":
        return "Ref"

    return method   # unrecognized shape - show it in full rather than guess


# ----------------------------------------------------------------------
# Shared colour scheme, used by EVERY chart in this file (plot_tree_overview,
# plot_total_volume_by_tree, plot_error_boxplot, plot_error_metrics_bar).
#
# WHY a single shared function: before this change, each chart picked its
# own colours independently (plot_tree_overview built its own ad-hoc
# palette, plot_total_volume_by_tree hard-coded "tab:orange" for the
# reference and left everything else to matplotlib's default cycle,
# plot_error_boxplot/plot_error_metrics_bar didn't colour by method at all).
# That meant the SAME method (e.g. "AdQSM (TreesParams)") could show up in a
# different colour in each chart, making it harder to visually track one
# method across the whole set of PNGs. This function is now the ONE place
# that decides "method -> colour", called once per branch_filter mode in the
# RUN section below and threaded into every chart that needs it, so the
# mapping is guaranteed identical everywhere it's used.
# ----------------------------------------------------------------------
def build_method_color_map(methods, reference_method):
    """Return a {method: color} dict for the given list of methods.

    Non-reference methods are split into FOUR groups by the FULL,
    untouched method string's prefix - deliberately independent of
    shorten_method_label() (this function never shortens a label, it only
    classifies the original string), so a future change to the short-label
    rules can never accidentally change which colour a method gets:
      - "AdTree raw ..."        -> pastel GREEN gradient
      - "AdTree calibrated ..." -> pastel BLUE gradient
      - "TreeQSM ..."           -> pastel GREY gradient (covers both
        "TreeQSM mine (...)" and "TreeQSM de Tanago (...)")
      - "AdQSM ..."             -> pastel YELLOW/OCHRE gradient (covers
        both "AdQSM (TreesParams) ..." and "AdQSM (BranchStructure...) ...")
    Each group is spread evenly across ITS OWN gradient's 0..1 range
    independently (same `t = i/(n-1)` logic as before, just computed once
    PER GROUP, same as the previous 2-group scheme) - so adding/removing a
    method in one group never shifts another group's shades.

    `reference_method` (the method this chart is comparing everything else
    against - REFERENCE_METHOD for the destructive-reference mode, or
    REFERENCE_METHOD_NONE for the AdQSM-as-yardstick mode) is checked
    against the four family prefixes FIRST, before any other
    classification:
      - If reference_method's prefix does NOT match any of the four
        families (e.g. REFERENCE_METHOD, "Reference (destructive)"), it
        is excluded from the four groups and instead gets a single fixed
        PINK/RED highlight colour (not a gradient, same "one flat
        highlight colour" treatment the old "coral" reference had - just
        a pink/red hue now), so it can never accidentally land on the
        same shade as one of the gradient-coloured methods.
      - If reference_method's prefix DOES match one of the four families
        (e.g. REFERENCE_METHOD_NONE, which literally IS an "AdQSM ..."
        method string, in the "none" mode), it is NOT isolated - it is
        classified into that family's group like any other member (sorted
        FIRST within that group) and shares that family's gradient with
        its siblings, so the "none" mode shows every AdQSM variant
        together in one gradient family instead of splitting the
        reference out into its own colour while its siblings stay
        yellow.

    Any non-reference method whose prefix matches NONE of "AdTree raw",
    "AdTree calibrated", "TreeQSM", "AdQSM" is unexpected (every method
    string produced by this project's pipeline today starts with one of
    those) - rather than crashing, it's printed as a clear warning and
    coloured a flat neutral grey (a different, flatter grey than the
    TreeQSM pastel-grey GRADIENT, so an unclassified method is still
    visually distinguishable from an actual TreeQSM one).

    Also prints, for every method in `methods`, which group (green/blue/
    grey/yellow/pink) it was classified into - a diagnostic so the colour
    assignment can be reviewed at a glance, same spirit as
    shorten_method_label()'s own diagnostic table printed in the RUN
    section below.
    """
    # Hand-picked, pastel (light, low-saturation) stops per group - four
    # stops each, same style as this file's previous hand-built gradients,
    # so a dozen+ methods in one group still stays readable instead of
    # clashing the way matplotlib's default saturated cycle would.
    green_gradient = mcolors.LinearSegmentedColormap.from_list(
        "adtree_raw_gradient",
        ["#e1f5e1", "#b8e2b8", "#8fce8f", "#63b563"],  # pastel green: pale -> sage -> leaf -> deeper green
    )
    blue_gradient = mcolors.LinearSegmentedColormap.from_list(
        "adtree_calib_gradient",
        ["#dceaf9", "#b3d1f2", "#84b3e8", "#5a92d6"],  # pastel blue: pale -> sky -> mid -> deeper blue
    )
    grey_gradient = mcolors.LinearSegmentedColormap.from_list(
        "treeqsm_gradient",
        ["#eeeeee", "#d4d4d4", "#b8b8b8", "#98989a"],  # pastel grey: near-white -> light -> mid -> deeper grey
    )
    yellow_gradient = mcolors.LinearSegmentedColormap.from_list(
        "adqsm_gradient",
        ["#fdf3c9", "#f8e08c", "#eec85a", "#d6a83f"],  # pastel yellow/ochre: pale -> gold -> ochre -> deeper ochre
    )

    family_prefixes = ("AdTree raw", "AdTree calibrated", "TreeQSM", "AdQSM")
    reference_is_family_member = any(
        reference_method.startswith(p) for p in family_prefixes) if reference_method else False
    print("  build_method_color_map(): reference '%s' %s a known family prefix - %s"
          % (reference_method,
             "matches" if reference_is_family_member else "does NOT match",
             "included within its family gradient (sorted first in that group), not isolated"
             if reference_is_family_member else "isolated with its own fixed highlight colour"))

    non_ref_methods = [m for m in methods if m != reference_method]

    # Classify BEFORE assigning colors, so each group's `t = i/(n-1)`
    # spread is computed over ONLY that group's own count - a method being
    # added/removed in one group must never shift another group's shades.
    # When the reference matches a family, it's put in the pool FIRST so it
    # lands at the start of its family's group list (see docstring above).
    classify_pool = ([reference_method] + non_ref_methods) if reference_is_family_member else non_ref_methods
    groups = {"AdTree raw": [], "AdTree calibrated": [], "TreeQSM": [], "AdQSM": []}
    unclassified = []
    for m in classify_pool:
        if m.startswith("AdTree raw"):
            groups["AdTree raw"].append(m)
        elif m.startswith("AdTree calibrated"):
            groups["AdTree calibrated"].append(m)
        elif m.startswith("TreeQSM"):
            groups["TreeQSM"].append(m)
        elif m.startswith("AdQSM"):
            groups["AdQSM"].append(m)
        else:
            unclassified.append(m)

    color_of = {}
    group_of_method = {}   # for the diagnostic print below only
    for group_name, group_gradient, color_name in (
            ("AdTree raw", green_gradient, "green"),
            ("AdTree calibrated", blue_gradient, "blue"),
            ("TreeQSM", grey_gradient, "grey"),
            ("AdQSM", yellow_gradient, "yellow")):
        group_methods = groups[group_name]
        n = len(group_methods)
        for i, m in enumerate(group_methods):
            # Spread methods evenly across the gradient's full 0..1 range. With
            # only one method in this group, t=0.5 (the middle of the gradient)
            # is used instead of dividing by (n - 1) = 0.
            t = (i / (n - 1)) if n > 1 else 0.5
            color_of[m] = group_gradient(t)
            group_of_method[m] = color_name

    if unclassified:
        print("WARNING: build_method_color_map(): %d method(s) matched none of "
              "'AdTree raw'/'AdTree calibrated'/'TreeQSM'/'AdQSM' and are not the "
              "reference method - falling back to flat neutral grey: %s"
              % (len(unclassified), unclassified))
        for m in unclassified:
            color_of[m] = "#9a9a9a"   # flat neutral grey - distinct from the TreeQSM grey GRADIENT
            group_of_method[m] = "grey (unclassified fallback)"

    if reference_method in methods and not reference_is_family_member:
        color_of[reference_method] = "#ef476f"   # fixed pink/red highlight (was "coral")
        group_of_method[reference_method] = "pink (reference)"

    print("  build_method_color_map(): colour-group assignment for this mode "
          "(reference='%s'):" % reference_method)
    for m in methods:
        print("    %-75s -> %s" % (m, group_of_method.get(m, "?")))

    return color_of


def order_with_reference_first(methods, reference_method):
    """Return `methods` with reference_method moved to the front (if
    present), keeping the relative order of everyone else unchanged.
    Used everywhere a method list drives x-axis/box/legend order, so
    the reference is always the first bar/box in every chart - easier
    to anchor on visually than wherever it happened to first appear
    in the CSV."""
    if reference_method in methods:
        return [reference_method] + [m for m in methods if m != reference_method]
    return methods


def group_order_methods(methods, reference_method):
    """Order `methods` into a fixed 4-group sequence, so chart x-axis/box/
    legend order is consistent across every chart regardless of which order
    methods happen to first appear in the CSV:

      1. reference_method's own family group (if reference_method's prefix
         matches one of "AdQSM"/"TreeQSM"/"AdTree raw"/"AdTree calibrated" -
         e.g. the "none" mode, where the reference IS an AdQSM variant):
         reference_method FIRST, then its remaining family siblings in
         their existing relative order. If reference_method matches NO
         family (e.g. "Reference (destructive)" in "10cm" mode), this
         group is just reference_method alone, first.
      2. AdTree raw   (remaining methods with this prefix, existing
         relative order, skipping anything already placed in group 1)
      3. AdTree calibrated (same pattern)
      4. TreeQSM      (same pattern)

    Any method matching none of the four family prefixes and not already
    placed as the reference is appended at the end, in its original
    relative order, with a printed warning - should not normally happen,
    matching this file's existing defensive-fallback style (see
    build_method_color_map()'s "unclassified" handling).
    """
    family_prefixes = ("AdTree raw", "AdTree calibrated", "TreeQSM", "AdQSM")
    reference_family = None
    if reference_method is not None:
        for p in family_prefixes:
            if reference_method.startswith(p):
                reference_family = p
                break

    placed = set()
    ordered = []

    # Group 1: reference_method's own family (or just itself, if no family match).
    if reference_method in methods:
        ordered.append(reference_method)
        placed.add(reference_method)
    if reference_family is not None:
        for m in methods:
            if m.startswith(reference_family) and m not in placed:
                ordered.append(m)
                placed.add(m)

    # Groups 2-4: the three remaining fixed families, in this fixed order,
    # skipping reference_family (already handled as group 1 above) and
    # anything already placed.
    for family in ("AdTree raw", "AdTree calibrated", "TreeQSM"):
        if family == reference_family:
            continue
        for m in methods:
            if m.startswith(family) and m not in placed:
                ordered.append(m)
                placed.add(m)

    # Anything left over matched no known family prefix - append as-is.
    leftover = [m for m in methods if m not in placed]
    if leftover:
        print("WARNING: group_order_methods(): %d method(s) matched none of "
              "'AdTree raw'/'AdTree calibrated'/'TreeQSM'/'AdQSM' - appending at "
              "the end in original order: %s" % (len(leftover), leftover))
        ordered.extend(leftover)

    return ordered


def _darken_color(color, factor=0.6):
    """Return a DARKER version of `color` (same hue, just scaled toward
    black) - used ONLY by plot_error_boxplot()'s per-tree scatter dots.

    WHY darker-same-color instead of a flat gray for the dots: the whole
    point of build_method_color_map() is "one consistent colour per method,
    everywhere in this file" - using plain gray dots for every method would
    throw that identity away right where it matters most (the individual
    observations you're overlaying so you can see per-tree spread, not just
    the summary box). A darker shade of the SAME colour keeps "which method
    is this" recognisable at a glance, while still reading as clearly
    different from the (lighter, semi-transparent) box fill underneath it.

    `factor` (0..1) is how much of the original colour to keep - 0.6 means
    "60% of the original brightness", which is dark enough to stand out
    against the box fill without going all the way to black (which would
    make every method's dots look identical again, defeating the purpose).
    """
    r, g, b = mcolors.to_rgb(color)   # to_rgb() accepts hex strings, named colours, AND RGBA tuples alike
    return (r * factor, g * factor, b * factor)


# ----------------------------------------------------------------------
# a) Bar chart: total_m3 per method, grouped by tree.
# ----------------------------------------------------------------------
def plot_total_volume_by_tree(rows, color_map):
    # Only branch_filter == "10cm" rows: this chart draws the REFERENCE
    # method's bar highlighted, so it's implicitly a "vs. reference" chart -
    # mixing in "none" (full/unfiltered) rows here would compare some
    # methods' full reconstructions against a reference that only ever
    # measured >= 10 cm, exactly the unfair comparison branch_filter exists
    # to prevent (see compare_volumes.py's header comment). Filtering here
    # too (not just in the RUN section below) means this function stays
    # correct even if called directly with unfiltered rows.
    rows = filter_by_branch_filter(rows, "10cm")

    trees = sorted({r["tree"] for r in rows})
    # dict.fromkeys() keeps the methods in the order they first appear in
    # the CSV (a plain set() would print them in a random order every run).
    # group_order_methods() then puts REFERENCE_METHOD (the destructive
    # reference - this chart is always "10cm" mode) first, followed by the
    # fixed AdTree raw / AdTree calibrated / TreeQSM group order, so bars
    # are always in the same order across every chart.
    methods = group_order_methods(list(dict.fromkeys(r["method"] for r in rows)), REFERENCE_METHOD)

    # Quick lookup: (tree, method) -> total_m3 (may be missing -> None).
    total_of = {(r["tree"], r["method"]): r["total"] for r in rows}

    n_methods = len(methods)
    # Each tree gets a group of bars 0.8 units wide, split evenly between
    # the methods, with a small gap (the leftover 0.2) to the next tree's group.
    bar_width = 0.8 / max(n_methods, 1)

    fig, ax = plt.subplots(figsize=(max(6, 1.6 * len(trees) * n_methods), 6))

    for i, method in enumerate(methods):
        heights = [total_of.get((t, method)) or 0.0 for t in trees]  # missing -> 0 (drawn as no bar)
        # x position of this method's bar within each tree's group
        x = [tree_idx + i * bar_width for tree_idx in range(len(trees))]
        is_ref = (method == REFERENCE_METHOD)
        # color_map (built once in the RUN section via build_method_color_map(),
        # shared across every chart) replaces the old hard-coded "tab:orange" -
        # this guarantees the reference bar here uses the EXACT same highlight
        # colour as every other chart's reference bar/box.
        ax.bar(x, heights, width=bar_width,
               color=color_map.get(method),
               label=shorten_method_label(method) + (" (reference)" if is_ref else ""))

    # Put the x-axis tick for each tree in the middle of its group of bars.
    group_centers = [tree_idx + bar_width * (n_methods - 1) / 2.0 for tree_idx in range(len(trees))]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(trees)
    ax.set_xlabel("Tree")
    ax.set_ylabel("Total volume [m^3]")
    ax.set_title("Total volume by method, per tree")
    # Legend outside the plot area (to the right) so it doesn't cover bars.
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    save_and_report(fig, "total_volume_by_tree.png")


# ----------------------------------------------------------------------
# d) One tree's overview: 2x3 grid of bar charts (total volume, DBH,
#    height, taper, trunk length, branch length), one bar per method, for
#    a SINGLE tree AND a SINGLE branch_filter ("none" or "10cm"). Makes
#    sense even with only 1 tree in the CSV (unlike the boxplot/RMSE
#    charts, which compare ACROSS trees).
#
#    branch_filter MUST be passed explicitly (no default) - drawing "none"
#    (full reconstruction) and "10cm" (diameter-restricted) rows in the SAME
#    bar chart would silently mix two different methodologies (very
#    different branch_len scale especially - see compare_volumes.py's
#    header comment for why they're kept apart everywhere else too).
# ----------------------------------------------------------------------
def plot_tree_overview(rows, tree, branch_filter, color_map):
    # Filter to THIS tree AND THIS branch_filter before anything else, so
    # every method/color/field computed below only ever sees rows from one
    # consistent methodology.
    tree_rows = [r for r in rows if r["tree"] == tree and r["branch_filter"] == branch_filter]
    if not tree_rows:
        print("No rows for tree '%s' with branch_filter='%s' - skipping this overview."
              % (tree, branch_filter))
        return

    # WHICH method is "the reference" depends on branch_filter, same as
    # everywhere else in this file: the destructive field reference
    # (REFERENCE_METHOD) only ever has "10cm" rows, so it can never be found
    # in "none"-mode data - AdQSM (REFERENCE_METHOD_NONE) plays that role
    # there instead. Before this fix, this line was hard-coded to
    # REFERENCE_METHOD, which meant ref_row was ALWAYS None in "none" mode
    # (since that method never appears there) - so the percent-difference
    # annotations further down were silently never drawn for "none" mode
    # charts, even though AdQSM WAS present and perfectly usable as a
    # yardstick. Computing the right reference per-mode here fixes that.
    # (Moved above the `methods` list below so group_order_methods()
    # can use it right away - it only depends on branch_filter, not on the
    # tree's actual rows, so computing it first is safe.)
    reference_method = REFERENCE_METHOD if branch_filter == "10cm" else REFERENCE_METHOD_NONE

    # Methods in the order they first appear for THIS tree (dict.fromkeys()
    # trick again, see plot_total_volume_by_tree), then reference_method
    # moved to the front (Task 1) - this order drives BOTH the per-subplot
    # bar order AND the shared legend order below, so the reference is
    # always the first bar/legend entry in every subplot of this figure.
    # The actual COLOUR per method still comes from `color_map` (built once
    # in the RUN section via build_method_color_map(), shared across every
    # chart in this file), not from an ad-hoc palette built locally here.
    methods = group_order_methods(
        list(dict.fromkeys(r["method"] for r in tree_rows)), reference_method)
    # One row dict per method, for quick lookups below (assumes at most one
    # row per (tree, method) pair, which is how upsert_result() keeps the CSV).
    row_of = {r["method"]: r for r in tree_rows}

    ref_row = row_of.get(reference_method)   # None if this tree has no row for that reference

    # color_of is just an alias into the shared color_map here (rather than
    # `color_map` directly) so the rest of this function's code below didn't
    # need to change when this was refactored to take color_map as a parameter.
    color_of = color_map

    # Which field (as returned by load_results()) goes in which subplot, and
    # what to title that subplot. This list is the ONLY thing you touch to
    # add/remove/reorder subplots - the loop below draws one subplot per
    # entry, so adding a field is a one-line change here, not a copy-pasted
    # block of plotting code. Just keep the subplot COUNT matching the grid
    # shape passed to plt.subplots() right below.
    #
    # ALL 8 fields load_results() provides are now shown (previously 6 - the
    # "stem"/"branch" volume panels below are NEW, added because you asked to
    # see per-method stem/branch volume side by side with everything else,
    # not just the combined total). Order groups the three VOLUME fields
    # first (total, stem, branch), then the rest in the same order as before.
    fields = [
        ("total",      "Total volume [m^3]"),
        ("stem",       "Stem volume [m^3]"),
        ("branch",     "Branch volume [m^3]"),
        ("dbh",        "DBH [m]"),
        ("height",     "Height [m]"),
        ("taper",      "Taper [cm/m]"),
        ("trunk_len",  "Trunk length [m]"),
        ("branch_len", "Branch length [m]"),
        # n_cylinders: was in volume_results.csv/load_results() already
        # (Task A), just never shown in any chart - added here as a 9th
        # panel so you can see, per method, how many cylinders its
        # reconstruction used (methods with no count, e.g. the destructive
        # reference, are simply skipped in this panel like any other
        # missing value - see the "present_methods" filter below).
        ("n_cylinders", "Number of cylinders"),
    ]

    # Unit shown next to the absolute-difference annotation below, matching
    # each field's own subplot_title unit above exactly (m^3 -> "m3" without
    # the caret, since it's plain annotation text, not an axis label; "" for
    # n_cylinders, which has no unit - see the annotation loop's "n_cylinders"
    # special case, which formats it as a plain signed integer instead).
    FIELD_UNITS = {
        "total": "m3", "stem": "m3", "branch": "m3",
        "dbh": "m", "height": "m",
        "taper": "cm/m",
        "trunk_len": "m", "branch_len": "m",
        "n_cylinders": "",
    }

    # GRID SIZE CHOICE: 3x3 (9 slots), matches the 9 fields above exactly -
    # no empty panel to explain away. (Previously 2x4/8 fields fit its own
    # grid perfectly too, for the same reason - adding n_cylinders as a 9th
    # field made 3x3 the new perfect fit instead of leaving an empty slot
    # in a 2x5 grid.) figsize scaled to keep each panel roughly the same
    # size as before (21/4 =~5.25 wide, 9/2 =~4.5 tall per panel before ->
    # 3*5.25 =~16 wide, 3*4.5 =~13.5 tall now). Height bumped again to
    # 3*5.3 =~16 tall (was 13.5) - taller panels give the rotated 90 deg
    # deviation annotations (see the annotation loop below) more vertical
    # room, on top of the ylim headroom fix, instead of relying on ylim
    # headroom alone to keep them from crossing into the panel above.
    fig, axes = plt.subplots(3, 3, figsize=(16, 16))

    for ax, (field_key, subplot_title) in zip(axes.flat, fields):
        # Skip methods with no value (None) for THIS field entirely, instead
        # of trying to draw a bar for them (which would crash matplotlib).
        present_methods = [m for m in methods if row_of[m][field_key] is not None]

        # Task 3: trunk_len/branch_len ONLY - drop "AdTree raw" variants from
        # THIS panel. Radius calibration only rescales/replaces cylinder
        # RADII, it never changes cylinder length or count (see
        # adtree_reconstruct_compare.py's comment on this same fact), so
        # "AdTree raw ..." and "AdTree calibrated ..." always draw IDENTICAL
        # bars here - showing both is pure visual clutter on these two
        # panels specifically (every other panel, including volume/DBH/
        # n_cylinders, is left completely untouched, since those genuinely
        # DO differ between raw and calibrated).
        if field_key in ("trunk_len", "branch_len"):
            present_methods = [m for m in present_methods if not m.startswith("AdTree raw")]

        values = [row_of[m][field_key] for m in present_methods]
        colors = [color_of[m] for m in present_methods]

        # Spacing factor > 1 widens the gap between bars (and thus between
        # their labels) without changing bar width itself - needed because
        # with many methods the rotated labels below start overlapping at
        # spacing=1 (bars packed edge-to-edge).
        bar_spacing = 1.6
        bar_width = 0.7
        x_positions = [i * bar_spacing for i in range(len(present_methods))]
        ax.bar(x_positions, values, width=bar_width, color=colors)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([shorten_method_label(m) for m in present_methods],
                           rotation=40, ha="right", fontsize=8)
        ax.set_title(subplot_title)

        # Small FYI note (not the heavy AdQSM data-quality warning further
        # below - this isn't a data-quality issue, just an explanation of why
        # fewer bars appear here than in the other panels) explaining the
        # omission above, so it isn't mistaken for missing/bad data.
        if field_key in ("trunk_len", "branch_len"):
            ax.text(0.98, 0.98,
                    "AdTree raw omitted - identical to calibrated\n(calibration only rescales radius)",
                    ha="right", va="top", fontsize=6, color="#666666", transform=ax.transAxes)

        if not present_methods:
            # Nothing to plot for this field at all (e.g. no method has a
            # taper value yet) - say so instead of leaving a blank mystery panel.
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            continue

        # Percent-difference label above every non-reference bar, using the
        # same pct_diff() calculation compare_volumes.py uses for its table.
        # Compared against `reference_method` (REFERENCE_METHOD for "10cm",
        # REFERENCE_METHOD_NONE for "none" - see the comment where that
        # variable is set above), NOT a hard-coded REFERENCE_METHOD, so this
        # annotation now actually appears in "none"-mode charts too.
        #
        # Also shows the SIGNED ABSOLUTE difference (value - ref_value), on
        # a second line, in this field's own unit (matching whatever unit
        # subplot_title already shows for it above - see FIELD_UNITS below).
        # WHY both: a percentage alone can make a small absolute difference
        # on a small reference value look dramatic, while an equally-small
        # absolute difference on a large reference value looks negligible -
        # showing both numbers lets the viewer judge practical relevance
        # themselves instead of only seeing one side of that picture.
        # Extra headroom above the tallest bar, set BEFORE the annotation
        # loop below (not after) - the two-line, rotated 90 deg deviation
        # annotations need vertical room reserved above the bars they sit
        # on, so the axis already has that room when the text is placed,
        # instead of relying on a post-hoc ylim bump to keep it from
        # running into subplot_title/the neighbouring panel. Based on
        # `values` (the actual bar heights) rather than whatever ylim
        # ax.bar() happened to autoscale to, so it's not sensitive to
        # matplotlib's own default margin. Bumped from the previous 18% to
        # 35% - 18% still left a few labels crossing into the panel above
        # at small fontsize/rotation=90 in practice.
        finite_values = [v for v in values if v is not None]
        if finite_values and max(finite_values) > 0:
            ax.set_ylim(top=max(finite_values) * 1.35)

        ref_value = ref_row[field_key] if ref_row is not None else None
        unit = FIELD_UNITS.get(field_key, "")
        for xi, m in zip(x_positions, present_methods):
            if m == reference_method:
                continue
            value = row_of[m][field_key]
            d_pct = pct_diff(value, ref_value)
            if d_pct is None:   # no reference value to compare against - leave blank
                continue
            d_abs = value - ref_value   # ref_value is guaranteed non-None here (d_pct was not None)
            if field_key == "n_cylinders":
                abs_text = "%+d" % round(d_abs)   # a cylinder COUNT - never fractional
            else:
                abs_text = ("%+.2f %s" % (d_abs, unit)).strip()
            ax.annotate("%+.0f%%\n(%s)" % (d_pct, abs_text),
                        xy=(xi, value),
                        xytext=(0, 3), textcoords="offset points",   # 3 points above the bar top
                        ha="center", va="bottom", fontsize=5, linespacing=1.15, rotation=90)

    # Subtitle spells out which methodology this figure shows, so it's clear
    # at a glance even without reading the filename.
    filter_label = ("full reconstruction, branch_filter='none'" if branch_filter == "none"
                     else "diameter >= 10 cm only, branch_filter='10cm'")
    fig.suptitle("Tree overview: %s  (%s)" % (tree, filter_label), fontsize=14)

    # ONLY for the "10cm" mode: a clearly-visible warning that AdQSM's
    # numbers in this filtered subset may not be trustworthy. WHY: AdQSM's
    # BranchStructure.txt has no reliable per-branch volume of its own for
    # this - its "volume(...)" column is off by orders of magnitude from
    # AdQSM's own official totals (see report_adqsm_thin_branch() in
    # tree_geom_utils.py), so any ">=10cm" AdQSM volume has to be
    # approximated as simple constant-radius cylinders, which measurably
    # over-estimates volume for tapering branches. This warning does NOT
    # apply to "none" mode, since that mode uses AdQSM's own official,
    # un-filtered TreesParams.txt totals directly - no cylinder
    # approximation involved there. Drawn as fig.text() (a separate, coloured
    # line - NOT folded into the small suptitle above) with a light
    # background box, specifically so it can't be mistaken for a routine
    # subtitle and skimmed past.
    if branch_filter == "10cm":
        fig.text(0.5, 0.955,
                 "NOTE: AdQSM values in this >=10cm subset are approximate - "
                 "AdQSM has no reliable per-branch volume source for this cut-off "
                 "(see BranchStructure.txt volume(...) column discussion).",
                 ha="center", va="top", fontsize=9, color="firebrick", fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff3cd", edgecolor="firebrick"))

    # ONE legend for the whole figure (not one per subplot, which would just
    # repeat the same method names 4 times) - built from coloured squares
    # ("patches") rather than real bar handles, since not every method has a
    # bar in every subplot. "(reference)" is tagged onto whichever method is
    # THIS mode's reference_method (see where that's computed above), not a
    # hard-coded REFERENCE_METHOD, so AdQSM correctly gets tagged in "none" mode.
    legend_handles = [
        mpatches.Patch(color=color_of[m],
                       label=shorten_method_label(m) + (" (reference)" if m == reference_method else ""))
        for m in methods
    ]
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=min(3, len(methods)), fontsize=8, bbox_to_anchor=(0.5, -0.02))

    # rect leaves room at the top for suptitle (+ the AdQSM warning line, when
    # present, in "10cm" mode) and at the bottom for the legend.
    top_margin = 0.90 if branch_filter == "10cm" else 0.96
    fig.tight_layout(rect=(0, 0.06, 1, top_margin))
    save_and_report(fig, "tree_overview_%s_%s.png" % (tree, branch_filter))


# ----------------------------------------------------------------------
# b) Box plot: percentage error vs. reference, per method, across trees.
#
#    branch_filter/reference_method are now REQUIRED parameters (instead of
#    the old hard-coded "10cm"/REFERENCE_METHOD) so this SAME function can
#    draw BOTH comparison modes from compare_volumes.py:
#      - branch_filter="10cm",  reference_method=REFERENCE_METHOD      (vs. the destructive reference)
#      - branch_filter="none",  reference_method=REFERENCE_METHOD_NONE (vs. AdQSM, methods-vs-each-other)
#    The output filename includes branch_filter (see save_and_report() call
#    below) so the two modes' PNGs never overwrite each other.
# ----------------------------------------------------------------------
def plot_error_boxplot(rows, branch_filter, reference_method, color_map):
    # Restrict to the requested branch_filter mode - same reasoning as
    # plot_total_volume_by_tree above: mixing "10cm" and "none" rows here
    # would compute a percent error against a reference method that isn't
    # even the right one for half the rows.
    rows = filter_by_branch_filter(rows, branch_filter)

    trees = sorted({r["tree"] for r in rows})
    # reference_method is already excluded here (this chart never draws a
    # box for it, only the axhline(0.0) below stands in for "perfect match
    # with the reference"), so group_order_methods()'s "group 1" is a no-op
    # in practice (reference_method isn't in the list to move) - applied
    # anyway for consistency with every other method-list build in this
    # file, and in case that exclusion ever changes.
    methods = group_order_methods(
        [m for m in dict.fromkeys(r["method"] for r in rows) if m != reference_method],
        reference_method)
    total_of = {(r["tree"], r["method"]): r["total"] for r in rows}

    if len(trees) < 3:
        print("NOTE: only %d tree(s) currently in %s (branch_filter='%s') - a box plot "
              "only becomes meaningful with 3+ trees (with fewer, each 'box' is really "
              "just 1-2 points). Drawing it anyway so you can see the layout."
              % (len(trees), RESULTS_CSV, branch_filter))

    data = []     # one list of % errors per method (only methods with >=1 value)
    labels = []
    for m in methods:
        errors_pct = []
        for t in trees:
            d = pct_diff(total_of.get((t, m)), total_of.get((t, reference_method)))
            if d is not None:
                errors_pct.append(d)
        if errors_pct:
            data.append(errors_pct)
            labels.append(m)

    if not data:
        print("No method has both a total_m3 value and a reference value "
              "(branch_filter='%s') - skipping box plot." % branch_filter)
        return

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(labels)), 6))
    # patch_artist=True turns the boxes into fillable patches (by default
    # boxplot() draws unfilled outlines only) so each box's facecolor can be
    # set from color_map below - this is what makes THIS chart's per-method
    # colours match every other chart's, instead of every box being the
    # same default matplotlib blue.
    # tick_labels gets the SHORTENED display text - `labels` itself stays
    # full-string throughout (it's still used below for color_map lookups
    # keyed by the full method string, e.g. the zip()s further down).
    bplot = ax.boxplot(data, tick_labels=[shorten_method_label(m) for m in labels], patch_artist=True)
    for patch, m in zip(bplot["boxes"], labels):
        # facecolor gets an ALPHA of 0.55 (via to_rgba, not patch.set_alpha())
        # specifically so only the FILL becomes semi-transparent - this is
        # what lets the individual-tree scatter dots (added further below)
        # show through the box clearly instead of being hidden underneath a
        # fully opaque one. Using to_rgba(..., alpha=...) instead of
        # patch.set_alpha() keeps the outline's own colour/opacity
        # independent of this, so setting the outline colour next isn't
        # also accidentally faded by the same alpha.
        patch.set_facecolor(mcolors.to_rgba(color_map.get(m, "#999999"), alpha=0.55))
        # Outline was matplotlib's default (black) - lightened to a soft gray
        # so it reads as a subtle boundary rather than a heavy border now
        # that the box interior is also busy with overlaid scatter points.
        patch.set_edgecolor("#888888")
    ax.axhline(0.0, color="gray", linestyle="--", linewidth=1)  # 0% error = perfect match

    # ---- overlay each tree's individual % error as a small dot ----------
    # WHY: with only 1-2 trees right now, the box itself is a degenerate
    # (near-meaningless) summary - showing the actual observations makes the
    # real data visible underneath the statistic, and stays useful later too
    # once more trees are added (you can see spread AND the box summary at
    # the same time, instead of choosing one or the other).
    #
    # boxplot() places method i's box at x = i + 1 (1-based, left to right in
    # `labels` order) - `data`/`labels` are already in that same order (built
    # together in the loop above), so zip()-ing them here lines up each
    # method's dots with its own box automatically.
    jitter_width = 0.08   # small horizontal spread, in the same x units as the boxes (box width = 0.5 by default)
    for i, (m, errors_pct) in enumerate(zip(labels, data)):
        # random.uniform() jitters each point sideways by a small random
        # amount so points with the same (or very close) y-value don't all
        # stack up in one indistinguishable vertical line - purely a visual
        # spread, it does NOT change any actual data value being plotted.
        x_jittered = [(i + 1) + random.uniform(-jitter_width, jitter_width) for _ in errors_pct]
        dot_color = _darken_color(color_map.get(m, "#999999"), factor=0.6)
        ax.plot(x_jittered, errors_pct, "o", color=dot_color, markersize=5,
                markeredgewidth=0, alpha=0.9, zorder=3)   # zorder=3: draw dots ON TOP of the (semi-transparent) boxes
    ax.set_ylabel("Error vs. reference [%]")
    ax.set_title("Total-volume error distribution by method (across trees)\n"
                  "vs. '%s'  (branch_filter='%s')" % (shorten_method_label(reference_method), branch_filter))
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    # Filename includes branch_filter so the "10cm" and "none" versions of
    # this chart are two separate files, e.g. error_boxplot_10cm.png /
    # error_boxplot_none.png, instead of the second run overwriting the first.
    save_and_report(fig, "error_boxplot_%s.png" % branch_filter)


# ----------------------------------------------------------------------
# c) Bar chart: Bias / MAE / RMSE per method (same numbers as
#    compare_volumes.py's printed "Error metrics" table).
#
#    branch_filter/reference_method - same reasoning as plot_error_boxplot()
#    above: this one function now draws both comparison modes.
# ----------------------------------------------------------------------
def plot_error_metrics_bar(rows, branch_filter, reference_method, color_map):
    # Restrict to the requested branch_filter mode - see plot_error_boxplot().
    rows = filter_by_branch_filter(rows, branch_filter)

    metrics = compute_error_metrics(rows, reference_method)   # reuses compare_volumes.py's own calculation
    if not metrics:
        print("No method could be compared to '%s' (branch_filter='%s') - "
              "skipping error-metrics chart." % (reference_method, branch_filter))
        return

    # compute_error_metrics() already excludes reference_method from its
    # results (same reasoning as plot_error_boxplot() above), so
    # group_order_methods()'s "group 1" is a no-op here too (reference_method
    # isn't in the list to move) - applied anyway for consistency with every
    # other method-list build in this file.
    methods = group_order_methods([m["method"] for m in metrics], reference_method)
    bias = [m["bias"] for m in metrics]
    mae = [m["mae"] for m in metrics]
    rmse = [m["rmse"] for m in metrics]

    x = list(range(len(methods)))
    width = 0.25   # 3 bars per method (bias, mae, rmse), each this wide

    fig, ax = plt.subplots(figsize=(max(6, 1.4 * len(methods)), 6))
    # NOTE on colour here: this chart groups THREE bars per METHOD (one each
    # for Bias/MAE/RMSE), so a single "one colour per method" mapping
    # (color_map) doesn't fit the same way it does in the other three charts
    # (there, each bar/box IS one method). Instead, the three metrics
    # themselves are drawn in a fixed mint/turquoise/teal triplet (a flat,
    # method-independent style choice for this chart specifically - not tied
    # to build_method_color_map()'s groups), and each method's x-axis tick
    # label is tinted with its color_map colour below - that's how this
    # chart still shows "this method = this colour", consistent with the rest.
    ax.bar([xi - width for xi in x], bias, width=width, label="Bias", color="#8fe0cf")   # mint
    ax.bar(x,                        mae,  width=width, label="MAE",  color="#4fbfae")   # turquoise
    ax.bar([xi + width for xi in x], rmse, width=width, label="RMSE", color="#2f8f8a")   # deep teal
    ax.axhline(0.0, color="gray", linestyle="-", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([shorten_method_label(m) for m in methods], rotation=30, ha="right")
    # Tint each x-axis tick label with that method's shared color_map colour,
    # so this chart still visually agrees with the other three on "which
    # colour is which method", even though its BARS are coloured by metric
    # (Bias/MAE/RMSE) rather than by method - see the comment above. Zips
    # against the FULL-string `methods` list (not the shortened tick text),
    # since color_map is keyed by the full method string.
    for tick_label, m in zip(ax.get_xticklabels(), methods):
        tick_label.set_color(color_map.get(m, "#333333"))
    ax.set_ylabel("Volume error [m^3]")
    ax.set_title("Error metrics vs. '%s'  (branch_filter='%s')" % (shorten_method_label(reference_method), branch_filter))
    ax.legend()
    fig.tight_layout()
    # Filename includes branch_filter, same reason as plot_error_boxplot() above.
    save_and_report(fig, "error_metrics_bar_%s.png" % branch_filter)


# =========================  RUN  =====================================
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        raise SystemExit(
            "'%s' not found - run compare_volumes.py (or one of the volume "
            "scripts, e.g. ply_to_geom.py) first so it gets created." % RESULTS_CSV
        )

    rows = load_results(RESULTS_CSV)
    all_trees = sorted({r["tree"] for r in rows})

    # ---- DIAGNOSTIC: full method string -> short display label ------------
    # Printed BEFORE shorten_method_label() gets wired into any chart below,
    # so the mapping can be reviewed against real CSV data before trusting
    # it - see that function's own docstring for the mapping rules. Every
    # DISTINCT method string currently in volume_results.csv (both modes
    # combined), sorted alphabetically.
    print("=" * 90)
    print("DIAGNOSTIC: method name -> short label mapping (review before trusting on charts)")
    print("=" * 90)
    for m in sorted({r["method"] for r in rows}):
        print("  %-75s -> %s" % (m, shorten_method_label(m)))
    print()

    # Pre-filtered once here too (mirrors compare_volumes.py's RUN section),
    # and passed into the "vs. reference" chart functions below - they also
    # filter internally (see each function's comment), so this is a
    # belt-and-suspenders double-filter: harmless (filtering "10cm" rows by
    # "10cm" again is a no-op) and keeps the RUN section's intent explicit.
    rows_10cm = filter_by_branch_filter(rows, "10cm")
    rows_none = filter_by_branch_filter(rows, "none")

    # Build the TWO shared colour maps used by every chart below - one per
    # branch_filter mode, since the "reference" method (and therefore which
    # method gets the highlight colour, and which methods share the
    # gradient) differs between modes: REFERENCE_METHOD (the destructive
    # reference) for "10cm", REFERENCE_METHOD_NONE (AdQSM) for "none". Built
    # from the FULL set of methods present in each mode (not per-chart
    # subsets), so a method's colour is guaranteed identical across ALL
    # charts that draw it in the same mode (plot_total_volume_by_tree,
    # plot_tree_overview, plot_error_boxplot, plot_error_metrics_bar all
    # receive the SAME dict for a given mode, instead of each recomputing
    # its own local mapping that could drift out of sync with the others).
    methods_10cm = list(dict.fromkeys(r["method"] for r in rows_10cm))
    methods_none = list(dict.fromkeys(r["method"] for r in rows_none))
    color_map_10cm = build_method_color_map(methods_10cm, REFERENCE_METHOD)
    color_map_none = build_method_color_map(methods_none, REFERENCE_METHOD_NONE)

    # Always makes sense, regardless of how many trees are in the CSV.
    plot_total_volume_by_tree(rows_10cm, color_map_10cm)

    # TWO overview PNGs per tree currently in the CSV - one per branch_filter
    # value, so "10cm" (vs.-reference) and "none" (full reconstruction) rows
    # are never drawn together in the same bar chart (see plot_tree_overview()'s
    # docstring for why mixing them would be misleading). If a tree has no
    # rows for one of the two filters, plot_tree_overview() just prints a
    # skip message for that one and moves on - harmless.
    for tree in all_trees:
        plot_tree_overview(rows, tree, "10cm", color_map_10cm)
        plot_tree_overview(rows, tree, "none", color_map_none)

    # The boxplot and RMSE/Bias/MAE charts compare methods ACROSS trees vs. a
    # reference, so what matters for EACH mode is how many trees have a row
    # in THAT mode (a tree could exist in the CSV with rows in only one of
    # the two modes) - not the raw tree count, and NOT shared between modes:
    # a tree with 2+ "10cm" rows but only 1 "none" row should still get the
    # "10cm" charts, just not the "none" ones (and vice versa).
    trees_10cm = sorted({r["tree"] for r in rows_10cm})
    trees_none = sorted({r["tree"] for r in rows_none})

    if len(trees_10cm) >= 2:
        plot_error_boxplot(rows_10cm, "10cm", REFERENCE_METHOD, color_map_10cm)
        plot_error_metrics_bar(rows_10cm, "10cm", REFERENCE_METHOD, color_map_10cm)
    else:
        print("Only %d tree(s) with branch_filter='10cm' rows in %s - skipping "
              "error_boxplot_10cm.png and error_metrics_bar_10cm.png (both compare "
              "methods ACROSS trees vs. the reference, so they need at least 2 such "
              "trees to be meaningful). Add more trees to the CSV and re-run to get them."
              % (len(trees_10cm), RESULTS_CSV))

    if len(trees_none) >= 2:
        plot_error_boxplot(rows_none, "none", REFERENCE_METHOD_NONE, color_map_none)
        plot_error_metrics_bar(rows_none, "none", REFERENCE_METHOD_NONE, color_map_none)
    else:
        print("Only %d tree(s) with branch_filter='none' rows in %s - skipping "
              "error_boxplot_none.png and error_metrics_bar_none.png (both compare "
              "methods ACROSS trees vs. AdQSM, so they need at least 2 such trees "
              "to be meaningful). Add more trees to the CSV and re-run to get them."
              % (len(trees_none), RESULTS_CSV))
