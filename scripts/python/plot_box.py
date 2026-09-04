# -*- coding: utf-8 -*-
# =====================================================================
#  Box plots comparing GROUPS of related method variants (not individual
#  methods) for one tree, across the same 9 metric fields shown in
#  plot_volumes.py's plot_tree_overview() panels.
# ---------------------------------------------------------------------
#  WHY this is a separate script from plot_volumes.py: plot_tree_overview()
#  draws one BAR per individual method (e.g. every single RADIUS_THRESHOLDS
#  x SEG_VARIANT_SUFFIX combination gets its own bar) - useful to see every
#  variant, but with many variants the bars/labels get crowded and it's hard
#  to see "how much does AdTree raw vary as a FAMILY" at a glance. This
#  script instead groups several related method-variant rows together (via
#  GROUP_RULES below) into ONE box per group, so the box's spread
#  (min/max/median/quartiles) shows the variation BETWEEN those variants'
#  values - not repeated measurements of a single method.
#
#  GROUP_RULES is a plain, user-editable list of (regex, group label) pairs,
#  matched against shorten_method_label(method) (the SAME shortened string
#  plot_tree_overview() already shows on its x-axis - not the raw,
#  untouched volume_results.csv "method" string), first match wins. Edit it
#  any time you want to regroup the same data differently (by segmentation
#  method, by radius threshold, by TreeQSM manual vs. auto, etc.) - nothing
#  else in this script needs to change.
#
#  Reused from elsewhere, not duplicated:
#    - compare_volumes.py: RESULTS_CSV, REFERENCE_METHOD, load_results(),
#      to_float() - the shared CSV loader, so this file can never drift out
#      of sync with how every other script reads volume_results.csv.
#    - plot_volumes.py: OVERVIEW_NCOLS (the same fixed-column-count
#      per-metric-field grid plot_tree_overview() uses - panel WIDTH/HEIGHT
#      and label styling are now this file's OWN BOX_* parameters instead,
#      see the PARAMETERS block below, since box plots need different
#      width/label tuning than the per-method overview bar charts),
#      shorten_method_label() (method -> short display string),
#      FAMILY_GRADIENTS/classify_family() (the shared "method family ->
#      colour" scheme - see plot_volumes.py's "Shared colour scheme"
#      section for the full rationale), and ensure_plots_dir()/PLOTS_DIR
#      (shared output folder).
#
#  Colour scheme note: since this chart colours GROUPS (not individual
#  methods), a group's colour is derived from classify_family() applied to
#  its FIRST member's raw, untouched CSV method string (by construction,
#  every member of one group should share a family - GROUP_RULES exists
#  precisely to keep families from mixing within a group; if a group's
#  members DO span more than one family, a warning names the mismatch and
#  the first member's family is used anyway, rather than crashing). Groups
#  within the SAME family are then spread across that family's existing
#  gradient by index, same t = i/(n-1) logic build_method_color_map() uses
#  for individual methods in plot_volumes.py.
#
#  Dependencies: matplotlib (install: pip install matplotlib). This script
#  also IMPORTS from compare_volumes.py and plot_volumes.py, which must
#  live in the same folder.
# =====================================================================

import math
import os
import re

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from compare_volumes import RESULTS_CSV, REFERENCE_METHOD, load_results, to_float
from plot_volumes import (
    OVERVIEW_NCOLS,
    shorten_method_label, FAMILY_GRADIENTS, classify_family,
    ensure_plots_dir,
    treeqsm_pd_token,
)

# to_float is imported for parity with the other scripts' import lists (see
# compare_volumes.py's load_results(), which already applies it to every
# numeric cell) - not called directly here, since load_results() already
# returns parsed floats/None for every field this script reads.
_ = to_float

# =====================  PARAMETERS  ==================================
# One tree at a time for now. The plotting logic below is organized as
# build_boxplot_figure(rows, tree, branch_filter) taking tree/branch_filter
# as plain arguments (same shape as plot_volumes.py's plot_tree_overview()),
# so a future side-by-side multi-tree view is a matter of looping that call
# once per tree, not restructuring anything in this file.
SELECT_TREE = "IND01_054"

# "none" = full reconstruction, "10cm" = >=10cm-only comparison - switch
# this to re-run for the other branch_filter variant.
BRANCH_FILTER = "10cm"

# Draw a horizontal dashed line at the destructive-reference value on every
# panel where a reference row/value exists for that field. Silently skipped
# (no warning) wherever it doesn't - that's expected for many fields (e.g.
# "none" mode never has a REFERENCE_METHOD row at all), not a data problem.
SHOW_REFERENCE_LINE = True

# Same treatment as SHOW_REFERENCE_LINE/REFERENCE_METHOD above, but for the
# TreeQSM published mean from qsm_volume_mean.py - which method string
# shows up depends on branch_filter (mutually exclusive, so at most ONE of
# these two is ever present for a given BRANCH_FILTER): "TreeQSM de Tanago
# (mean)" for "none", "TreeQSM de Tanago (mean, Filtered<10cm)" for "10cm".
SHOW_TREEQSM_REF_LINE = True
TREEQSM_REF_METHODS = ["TreeQSM de Tanago (mean)",
                        "TreeQSM de Tanago (mean, Filtered<10cm)"]
TREEQSM_REF_LINE_COLOR = "#2a9d8f"   # teal - visually distinct from the destructive reference's pink/#ef476f

SHOW_PLOT = False
SAVE_PLOT_PNG = True   # always save, consistent with this project's convention elsewhere

# ---- panel sizing / label styling - INDEPENDENT from plot_volumes.py's
# PANEL_WIDTH/PANEL_HEIGHT/LABEL_FONTSIZE/LABEL_ROTATION/BOTTOM_MARGIN -----
# Box plots need different width tuning than plot_tree_overview()'s
# per-method bar panels (far fewer x-axis entries per panel here - one per
# GROUP, not one per individual method variant - so the same wide panels
# felt oversized). Tune these freely; they no longer affect, or get
# affected by, plot_volumes.py's own panel styling.
BOX_PANEL_WIDTH = 5     # inches per panel
BOX_PANEL_HEIGHT = 6  # inches per panel
BOX_LABEL_FONTSIZE = 9  # x-axis group-label font size
BOX_LABEL_ROTATION = 30  # x-axis group-label rotation angle (degrees)
BOX_BOTTOM_MARGIN = 0.15 # fraction of figure height reserved for x-axis labels
BOX_TITLE_FONTSIZE = 12 # per-panel title font size

# Which GROUP a method belongs to, decided by matching shorten_method_label
# (method) - the SAME shortened string plot_tree_overview() shows on its
# x-axis, not the raw CSV method string - against these regexes IN ORDER;
# the FIRST pattern that matches wins. A method matching NO pattern is
# skipped from the plot entirely (with one warning line printed per skipped
# method) rather than dumped into a catch-all "(ungrouped)" box, since that
# would mix unrelated method families together and defeat the point of
# family-consistent colouring. Edit this list any time you want to regroup
# the same data differently - nothing else in this script needs to change.
#   (r"TreeQSM mine.*Optimal", "TreeQSM manual/Optimal"),
#   (r"TreeQSM.*[Aa]uto",      "TreeQSM auto"),
#   (r"^AQ_Params_",           "AdQSM"
GROUP_RULES = [
    # (regex matched against shorten_method_label(method), group label)
    (r"^AT_Raw_",              "AdTree raw"),
    (r"^AT_Calib_\d+_04_",     "AdTree calib (seg 04)"),
    (r"^AT_Calib_\d+_05_",     "AdTree calib (seg 05)"),
    (r"^AT_Calib_\d+_06_",     "AdTree calib (seg 06)"),
    #(r"^AT_Calib_\d+_07_",     "AdTree calib (seg 07)"),
    #(r"^AT_Calib_\d+_08_",     "AdTree calib (seg 08)"),
    #(r"^AT_Calib_\d+_09_",     "AdTree calib (seg 09)"),
    #(r"^AT_Calib_\d+_10",     "AdTree calib (seg 10)"),
    #(r"TreeQSM mine.*Optimal", "TreeQSM manual/Optimal"),
    #(r"TreeQSM mine.*Simplified (no islands)", "TreeQSM manual/Simplified"),
    #(r"TreeQSM mine.*Filtered <10cm", "TreeQSM manual/Filtered")
]

# ----------------------------------------------------------------------
# STRUCTURED TreeQSM grouping - a SEPARATE mechanism from GROUP_RULES
# above, parallel to it, for "TreeQSM mine (...)" rows specifically (any
# row whose raw method string starts with "TreeQSM mine (" - the same
# check plot_volumes.py's shorten_method_label() already uses to recognize
# this row shape). GROUP_RULES itself is untouched and keeps handling
# every AdTree/AdQSM row exactly as before; both mechanisms' resulting
# groups are merged into ONE combined list that drives the same box-plot
# panels, TreeQSM groups appended after the GROUP_RULES ones.
#
# WHY a separate, non-regex mechanism: TreeQSM rows written via the
# params_*.csv pipeline (see runsken.m section 19 / compare_volumes.py's
# load_results()) carry their ACTUAL reconstruction parameters as real
# structured columns (mode, pd1, pd2min, pd2max, simp_maxorder,
# simp_smallradii, simp_replaceiterations) - grouping/filtering on those
# directly is simpler and less fragile than regex-parsing a text label.
# Older "TreeQSM mine (v1aut, ...)"/"(v1man, ...)" rows predate that
# pipeline (mode == "" for them) and are excluded from THIS mechanism with
# one batched warning - same as today's silent exclusion by GROUP_RULES,
# not a regression.
TREEQSM_STAGE_FILTER = "Filtered <10cm"   # one of "Optimal", "Simplified", "Simplified (no islands)",
                                            # "Filtered <10cm", or None (show all stages mixed - not
                                            # recommended, but not blocked either)

TREEQSM_VARY_BY = ["ri"]        # short param name(s) (see TREEQSM_PARAM_SHORT_NAMES below) that become
                                  # the per-point spread INSIDE each box; every other short name becomes
                                  # part of the box-defining group key

TREEQSM_PARAM_FILTERS = {}      # optional: {short_name: value} - restrict to rows matching these exact
                                  # values (numeric values matched with a +/-1e-4 tolerance, to absorb
                                  # float round-tripping through the CSV's "%.6f" text - not to blur
                                  # genuinely different settings, which in practice differ by >=0.01)
                                  # only applied if this dict is non-empty. e.g. {"pd1": 0.08}

TREEQSM_PARAM_SHORT_NAMES = {   # short name -> actual load_results() row dict key, used by both
                                  # TREEQSM_VARY_BY and TREEQSM_PARAM_FILTERS above (and to build box/
                                  # point labels - see treeqsm_param_tokens() below)
    "mode": "mode", "pd1": "pd1", "pd2min": "pd2min", "pd2max": "pd2max",
    "maxorder": "simp_maxorder", "sr": "simp_smallradii", "ri": "simp_replaceiterations",
}
# ----------------------------------------------------------------------



# Label each jittered point with the specific varying parameter that
# distinguishes it from its group-mates (e.g. the radius threshold for
# AdTree groups, the AdQSM reconstruction variant number for AdQSM) -
# purely cosmetic, so a group/member with no entry (or a member whose
# shorten_method_label(method) doesn't match its group's regex) simply
# gets no point label, with no warning.
SHOW_POINT_LABELS = True    # master on/off switch
POINT_LABEL_FONTSIZE = 6    # smaller than BOX_LABEL_FONTSIZE - many dense small labels next to each other

JITTER_RANGE = 0.08            # half-width of horizontal point jitter (data units)
JITTER_POINT_SIZE = 20       # scatter marker size (matplotlib's s=)
JITTER_POINT_ALPHA = 0.8     # 0 = fully transparent, 1 = fully opaque
JITTER_DARKEN_FACTOR = 0.4     # 0 = same colour as box fill, 1 = black - controls jitter POINT colour,
                                # independently from LABEL_DARKEN_FACTOR (point labels) and the box edge's
                                # own darken_color() call (still its default factor, unparametrized for now)
LABEL_VERTICAL_STAGGER = 6    # extra points of vertical offset per point index within a group - see the
                                # annotate() call below for why (horizontal jitter alone doesn't separate
                                # labels for points whose actual VALUES are near-identical, e.g. Trunk
                                # volume's AdTree raw group)

BOX_TOP_MARGIN = 0.15   # fraction of each panel's own data range added above its highest whisker/point/
                         # label, so the topmost whisker/point/label never sits right against (or gets
                         # clipped by) that panel's title - applied PER PANEL (each of the 9 subplots gets
                         # its own independent headroom, based on ITS OWN autoscaled range), not one shared
                         # ylim across all panels.

LABEL_DARKEN_FACTOR = 1   # 0 = same colour as box fill, 1 = black - controls POINT-LABEL TEXT colour
                             # only, independently from the box edge/whisker/cap/median/jitter-point colour
                             # (still darken_color()'s own default factor=0.4, unrelated to this constant -
                             # see the annotate() call below).

# {group_label: (regex, format)} - regex matched against the SAME
# shorten_method_label(method) string GROUP_RULES already used to place
# that row into this group (NOT re-matched independently against some
# other candidate group), with ONE capture group pinpointing the value to
# show; format ("%s"-style) wraps that capture for display.
POINT_LABEL_PATTERNS = {
    "AdTree raw":              (r"^AT_Raw_(\d+)_",      "r%s"),   # "AT_Raw_5_..." -> "r5"
    "AdTree calib (seg 04)":   (r"^AT_Calib_(\d+)_04_", "r%s"),   # "AT_Calib_5_04_..." -> "r5"
    "AdTree calib (seg 05)":   (r"^AT_Calib_(\d+)_05_", "r%s"),
    "AdTree calib (seg 06)":   (r"^AT_Calib_(\d+)_06_", "r%s"),
    "AdQSM":                   (r"^AQ_Params_(\d+)$",   "v%s"),   # "AQ_Params_04" -> "v04"
    # TreeQSM entries intentionally omitted - add once real TreeQSM rows
    # exist in volume_results.csv and shorten_method_label()'s actual
    # output shape can be verified, rather than shipping an unverified guess.
}
# =====================================================================


# ----------------------------------------------------------------------
# The 9 metric fields shown - SAME set, same order, same display labels/
# units as plot_volumes.py's plot_tree_overview() (see its own `fields`/
# FIELD_UNITS, which aren't module-level there so they're restated here,
# not imported - kept manually in sync since both describe the exact same
# load_results() row shape).
# ----------------------------------------------------------------------
FIELDS = [
    ("total",       "Total volume [m^3]"),
    ("trunk",       "Trunk volume [m^3]"),
    ("branch",      "Branch volume [m^3]"),
    ("dbh",         "DBH [m]"),
    ("height",      "Height [m]"),
    ("taper",       "Taper [cm/m]"),
    ("trunk_len",   "Trunk length [m]"),
    ("branch_len",  "Branch length [m]"),
    ("n_cylinders", "Number of cylinders"),
]


def assign_groups(rows):
    """Bucket non-reference rows into groups via GROUP_RULES.

    Returns (groups, group_order):
      groups       {group_label: [row, ...]} - rows in their original order
      group_order  [group_label, ...] - GROUP_RULES' own order, restricted
                    to labels that actually got at least one row, so the
                    x-axis/box order is deterministic and author-controlled
                    (via GROUP_RULES' own ordering) rather than depending on
                    whatever order rows happen to appear in the CSV.

    Any row whose shorten_method_label(method) matches NO pattern in
    GROUP_RULES is printed as a warning (one line each) and excluded.
    """
    groups = {}
    skipped = []
    for r in rows:
        label = shorten_method_label(r["method"])
        group_label = None
        for pattern, candidate_label in GROUP_RULES:
            if re.search(pattern, label):
                group_label = candidate_label
                break
        if group_label is None:
            skipped.append(r["method"])
            continue
        groups.setdefault(group_label, []).append(r)

    if skipped:
        print("WARNING: plot_box.py: %d method(s) matched no GROUP_RULES pattern - "
              "excluded from the plot:" % len(skipped))
        for m in skipped:
            print("    %s  (shortened: %s)" % (m, shorten_method_label(m)))

    group_order = [label for _, label in GROUP_RULES if label in groups]
    return groups, group_order


def parse_treeqsm_method(method):
    """Split a "TreeQSM mine (run, stage)" method string into (run, stage),
    or return None if `method` isn't shaped that way - the SAME parsing
    plot_volumes.py's shorten_method_label() uses to recognize this row
    shape (see its own "TreeQSM mine (" branch)."""
    if not (method.startswith("TreeQSM mine (") and method.endswith(")")):
        return None
    inner = method[len("TreeQSM mine ("):-1]
    if ", " not in inner:
        return None
    run, stage = inner.split(", ", 1)
    return run, stage


def format_treeqsm_param_token(short_name, value):
    """Format ONE structured TreeQSM param into its compact label token
    (e.g. "ri" + 0 -> "r0"), using format B's conventions (m/s/r instead
    of the old mo/sr/ri) - shorten_method_label()'s TreeQSM branch and
    plot_box.py's box/point labels stay consistent with each other this way."""
    if short_name == "mode":
        return str(value)[:3]   # "manual" -> "man", "auto" -> "aut"
    if short_name == "maxorder":
        return "m%d" % int(value)
    if short_name == "sr":
        return "s%d" % round(value * 1000)
    if short_name == "ri":
        return "r%d" % int(value)
    if short_name in ("pd1", "pd2min", "pd2max"):
        return "%s%d" % (short_name, round(value * 100))
    return "%s%s" % (short_name, value)   # fallback for any future short name added later


def treeqsm_param_tokens(row, short_names):
    """Build the list of compact tokens for `short_names` (an ordered
    subset of TREEQSM_PARAM_SHORT_NAMES' keys) from `row`'s actual values,
    IN `short_names`' OWN ORDER (e.g. "mode" first if it's first in
    `short_names`) - so the caller's chosen order (TREEQSM_PARAM_SHORT_NAMES'
    declared order, by default) is what ends up in the label, matching
    shorten_method_label()'s "mode, pd, m, s, r" format.

    pd1/pd2min/pd2max are combined into ONE "p#-#-#" token, via the SAME
    treeqsm_pd_token() shorten_method_label() uses (format B), when all
    three are present together in `short_names` - the common case (they're
    always fixed or all varied together in practice) - emitted at the
    position of the FIRST of the three in `short_names`. If only SOME of
    the three are present (an unusual TREEQSM_VARY_BY split), each present
    one falls back to its own individual "pd#"-style token instead, since
    there's no established combined format for a partial trio.
    """
    pd_trio_present = all(k in short_names for k in ("pd1", "pd2min", "pd2max"))
    tokens = []
    pd_emitted = False
    for short in short_names:
        if pd_trio_present and short in ("pd1", "pd2min", "pd2max"):
            if pd_emitted:
                continue   # already emitted the combined token at the trio's first position
            pd1 = row[TREEQSM_PARAM_SHORT_NAMES["pd1"]]
            pd2min = row[TREEQSM_PARAM_SHORT_NAMES["pd2min"]]
            pd2max = row[TREEQSM_PARAM_SHORT_NAMES["pd2max"]]
            tokens.append(treeqsm_pd_token(pd1, pd2min, pd2max))
            pd_emitted = True
            continue
        value = row[TREEQSM_PARAM_SHORT_NAMES[short]]
        if value is None:
            continue
        tokens.append(format_treeqsm_param_token(short, value))
    return tokens


def assign_treeqsm_groups(rows):
    """Bucket "TreeQSM mine (...)" rows into groups via the STRUCTURED
    TREEQSM_STAGE_FILTER/TREEQSM_VARY_BY/TREEQSM_PARAM_FILTERS mechanism -
    parallel to (and independent of) assign_groups()/GROUP_RULES above.

    Returns (groups, group_order) in the SAME shape as assign_groups(), so
    the two can be merged into one combined groups dict / group_order list
    by the caller: {group_label: [row, ...]} / [group_label, ...], the
    latter in first-CSV-appearance order (there's no author-curated rule
    list here the way GROUP_RULES has one - group labels are discovered
    dynamically from whatever fixed-param combinations actually appear).

    Rows are excluded (no error) when: their stage doesn't match
    TREEQSM_STAGE_FILTER (silent - that's the filter's whole point), they
    don't satisfy TREEQSM_PARAM_FILTERS (also silent, same reason), or
    they have no structured params at all (mode == "", an older pre-
    params_*.csv row) - THAT case is reported as one batched warning,
    since it's a data-shape issue rather than an intentional filter.
    """
    fixed_short_names = [s for s in TREEQSM_PARAM_SHORT_NAMES if s not in TREEQSM_VARY_BY]

    rows_by_key = {}
    key_order = []
    skipped_no_params = []
    for r in rows:
        parsed = parse_treeqsm_method(r["method"])
        if parsed is None:
            continue   # not actually "TreeQSM mine (...)" shaped - shouldn't happen given the caller's own routing check
        _run, stage = parsed

        if TREEQSM_STAGE_FILTER is not None and stage != TREEQSM_STAGE_FILTER:
            continue

        if not r["mode"]:
            skipped_no_params.append(r["method"])
            continue

        if TREEQSM_PARAM_FILTERS:
            matched = True
            for short, target in TREEQSM_PARAM_FILTERS.items():
                value = r[TREEQSM_PARAM_SHORT_NAMES[short]]
                if isinstance(target, str):
                    if value != target:
                        matched = False
                        break
                else:
                    if value is None or abs(value - target) >= 1e-4:
                        matched = False
                        break
            if not matched:
                continue

        fixed_key = tuple(r[TREEQSM_PARAM_SHORT_NAMES[s]] for s in fixed_short_names)
        if fixed_key not in rows_by_key:
            rows_by_key[fixed_key] = []
            key_order.append(fixed_key)
        rows_by_key[fixed_key].append(r)

    if skipped_no_params:
        print("WARNING: plot_box.py: %d TreeQSM row(s) have no structured reconstruction "
              "parameters (older import, from before params_*.csv existed) - excluded from "
              "TreeQSM grouping: %s" % (len(skipped_no_params), skipped_no_params))

    # Turn the (fixed_key -> rows) buckets into the {label: rows}/[label,...]
    # shape assign_groups() already produces, using each group's formatted
    # token string (treeqsm_param_tokens(), joined with spaces) as its
    # actual label from here on.
    groups = {}
    group_order = []
    for fixed_key in key_order:
        members = rows_by_key[fixed_key]
        # "TQ_" prefix + underscore-joined tokens, matching
        # shorten_method_label()'s TreeQSM format (e.g. "TQ_man_p7-2-7_m8_s5") -
        # fixed_short_names is TREEQSM_PARAM_SHORT_NAMES' own declared order
        # (mode, pd1/pd2min/pd2max, maxorder, sr) minus whatever's in
        # TREEQSM_VARY_BY, so "mode" naturally comes first.
        label = "TQ_" + "_".join(treeqsm_param_tokens(members[0], fixed_short_names))
        groups[label] = members
        group_order.append(label)
    return groups, group_order


def classify_groups(groups, group_order):
    """Return {group_label: family_name} - one family per group, derived
    from classify_family() on the group's FIRST member's raw, untouched CSV
    method string. If any OTHER member of that group classifies into a
    DIFFERENT family (a GROUP_RULES authoring mistake - members of one
    group are supposed to share a family), a warning names the mismatch;
    the first member's family is used regardless, rather than crashing."""
    family_of_group = {}
    for label in group_order:
        members = groups[label]
        first_family = classify_family(members[0]["method"], REFERENCE_METHOD)
        for m in members[1:]:
            m_family = classify_family(m["method"], REFERENCE_METHOD)
            if m_family != first_family:
                print("WARNING: plot_box.py: group '%s' mixes families - "
                      "'%s' classified as %s, but '%s' classified as %s. "
                      "Using %s (the first member's family) for this "
                      "group's colour - check GROUP_RULES if this is "
                      "unintended." % (label, members[0]["method"], first_family,
                                       m["method"], m_family, first_family))
                break
        family_of_group[label] = first_family
    return family_of_group


def build_group_color_map(group_order, family_of_group):
    """Return {group_label: color}, spreading groups within the SAME family
    across that family's FAMILY_GRADIENTS gradient by index (same t =
    i/(n-1) logic build_method_color_map() in plot_volumes.py uses for
    individual methods) - so groups sharing a family still get visually
    distinct, but clearly related, shades."""
    family_gradient = {
        family: mcolors.LinearSegmentedColormap.from_list(
            "%s_gradient" % family.lower().replace(" ", "_"), stops)
        for family, stops in FAMILY_GRADIENTS.items()
    }

    groups_by_family = {family: [] for family in FAMILY_GRADIENTS}
    unclassified = []
    for label in group_order:
        family = family_of_group[label]
        if family is None:
            unclassified.append(label)
        else:
            groups_by_family[family].append(label)

    color_of_group = {}
    for family, labels in groups_by_family.items():
        n = len(labels)
        for i, label in enumerate(labels):
            t = (i / (n - 1)) if n > 1 else 0.5
            color_of_group[label] = family_gradient[family](t)

    if unclassified:
        print("WARNING: plot_box.py: %d group(s) classified into NO known family - "
              "falling back to flat neutral grey: %s" % (len(unclassified), unclassified))
        for label in unclassified:
            color_of_group[label] = "#9a9a9a"

    return color_of_group


def darken_color(color, factor=0.4):
    """Return a DARKER version of `color` (same hue, scaled toward black) -
    used for each box's edge/whisker/cap/median colour AND its scattered
    data-point colour, so both read clearly against that box's own
    (lighter) fill colour instead of a flat black/gray that would look
    identical across every box regardless of its family.

    Plain RGB-channel scaling (multiply each of r/g/b by `1 - factor`), not
    an HSV round-trip - simpler, and for a "make it darker without shifting
    hue" goal the two give visually indistinguishable results; RGB scaling
    is also what this project's plot_volumes.py already uses for the exact
    same purpose (its own _darken_color(), used for plot_error_boxplot()'s
    per-tree scatter dots) - kept consistent with that rather than
    introducing a second, differently-behaved darkening method.

    `factor` (0..1) is how much to darken TOWARD black - 0 means "same
    colour as the fill, unchanged", 1 means "black", 0.4 (this function's
    own default, and every DARKEN_FACTOR parameter's own default) is dark
    enough to stand out against the box fill without going all the way to
    black. (Fixed from an earlier version of this function where `factor`
    meant the OPPOSITE - "brightness kept", so factor=1 was a no-op and
    factor=0 was black - which silently contradicted every DARKEN_FACTOR
    parameter's own documented contract; every existing default below was
    rescaled at the same time so nothing's on-screen appearance changed.)
    """
    r, g, b = mcolors.to_rgb(color)   # to_rgb() accepts hex strings, named colours, AND RGBA tuples alike
    keep = 1 - factor
    return (r * keep, g * keep, b * keep)


def build_boxplot_figure(rows, tree, branch_filter):
    """Build and save/show ONE figure: a 2-column grid of box-plot panels
    (one per FIELDS entry), each panel showing one box per GROUP_RULES
    group, for `tree`/`branch_filter`. Same (rows, tree, branch_filter)
    argument shape as plot_volumes.py's plot_tree_overview(), so a future
    multi-tree view is just a loop around this call."""
    tree_rows = [r for r in rows if r["tree"] == tree and r["branch_filter"] == branch_filter]
    if not tree_rows:
        print("No rows for tree '%s' with branch_filter='%s' - skipping this boxplot."
              % (tree, branch_filter))
        return

    # Reference row set aside (for the horizontal reference lines) and
    # excluded from grouping/boxing - it's a single destructive-reference
    # measurement, not a group of variants to build a box from.
    ref_row = next((r for r in tree_rows if r["method"] == REFERENCE_METHOD), None)

    # TreeQSM published-mean row, same treatment as ref_row above - also
    # set aside and excluded from grouping/boxing (it's a single mean
    # value, not a group of variants). At most one of TREEQSM_REF_METHODS
    # is ever present for a given branch_filter (mutually exclusive), so
    # this is a single lookup, not a list.
    treeqsm_ref_row = next((r for r in tree_rows if r["method"] in TREEQSM_REF_METHODS), None)

    non_ref_rows = [r for r in tree_rows if r is not ref_row and r is not treeqsm_ref_row]

    # Route each row to EITHER the existing regex-based GROUP_RULES
    # mechanism (AdTree/AdQSM/anything else) OR the new structured
    # TreeQSM mechanism (assign_treeqsm_groups()) - same
    # "TreeQSM mine (" prefix check parse_treeqsm_method() itself uses, so
    # a row can never accidentally go through both. Both mechanisms return
    # the SAME {label: rows}/[label,...] shape, so their results merge into
    # ONE combined groups/group_order that drives the box-drawing loop
    # below - AdTree/AdQSM groups first, TreeQSM groups appended after
    # (per Bara's request that they appear side by side, not as two
    # separate charts). group_kind tracks which mechanism produced each
    # label, so the per-point-label code further down knows which of the
    # two label-building strategies to use for it.
    group_rules_rows = [r for r in non_ref_rows if not r["method"].startswith("TreeQSM mine (")]
    treeqsm_candidate_rows = [r for r in non_ref_rows if r["method"].startswith("TreeQSM mine (")]

    groups, group_order = assign_groups(group_rules_rows)
    treeqsm_groups, treeqsm_group_order = assign_treeqsm_groups(treeqsm_candidate_rows)

    group_kind = {label: "rules" for label in group_order}
    group_kind.update({label: "treeqsm" for label in treeqsm_group_order})

    groups = {**groups, **treeqsm_groups}
    group_order = group_order + treeqsm_group_order

    if not group_order:
        print("No groups matched any GROUP_RULES pattern, and no TreeQSM row matched "
              "TREEQSM_STAGE_FILTER/TREEQSM_PARAM_FILTERS, for tree '%s' branch_filter='%s' - "
              "skipping this boxplot." % (tree, branch_filter))
        return

    family_of_group = classify_groups(groups, group_order)
    color_of_group = build_group_color_map(group_order, family_of_group)

    # Reference line colour: the exact "#ef476f" hex FAMILY_GRADIENTS'
    # "Reference" family is built around (its middle stop - see
    # plot_volumes.py's FAMILY_GRADIENTS definition), so the reference line
    # here matches the pink reference highlight used everywhere else.
    reference_line_color = FAMILY_GRADIENTS["Reference"][1]

    n_fields = len(FIELDS)
    n_rows = math.ceil(n_fields / OVERVIEW_NCOLS)
    fig, axes = plt.subplots(n_rows, OVERVIEW_NCOLS,
                              figsize=(OVERVIEW_NCOLS * BOX_PANEL_WIDTH, n_rows * BOX_PANEL_HEIGHT))
    for ax in axes.flat[n_fields:]:
        ax.axis("off")

    for ax, (field_key, subplot_title) in zip(axes.flat, FIELDS):
        box_data = []
        box_labels = []
        box_colors = []
        box_methods = []   # per-point raw CSV method strings, same shape/order as box_data - for GROUP_RULES-origin point labels below
        box_rows = []      # per-point full row dicts, same shape/order as box_data - for TreeQSM-origin point labels below (needs actual param values, not just the method string)
        for label in group_order:
            present = [r for r in groups[label] if r[field_key] is not None]
            if not present:
                print("WARNING: plot_box.py: group '%s' has zero valid values for "
                      "field '%s' (tree '%s', branch_filter='%s') - no box drawn "
                      "for it in this panel." % (label, field_key, tree, branch_filter))
                continue
            box_data.append([r[field_key] for r in present])
            box_labels.append(label)
            box_colors.append(color_of_group[label])
            box_methods.append([r["method"] for r in present])
            box_rows.append(present)

        ax.set_title(subplot_title, fontsize=BOX_TITLE_FONTSIZE)
        ax.grid(False)

        if not box_data:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
            continue

        # whis=(0, 100): whiskers extend all the way to each group's actual
        # min/max, instead of the default 1.5*IQR cutoff - so the box+
        # whiskers alone show the TRUE range, with no separately-plotted
        # "outlier" dots (there technically can't be any past whis=100).
        bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True, whis=(0, 100))

        # Per-box styling: fill = that group's family colour (as before);
        # edge/whiskers/caps/median = a DARKER shade of that SAME colour
        # (darken_color(), above) instead of matplotlib's default black -
        # so a box's outline/median stay visually tied to its own family
        # colour rather than looking identical across every box.
        edge_colors = [darken_color(c) for c in box_colors]
        for patch, fill, edge in zip(bp["boxes"], box_colors, edge_colors):
            patch.set_facecolor(fill)
            patch.set_edgecolor(edge)
        for i, edge in enumerate(edge_colors):
            bp["whiskers"][2 * i].set_color(edge)
            bp["whiskers"][2 * i + 1].set_color(edge)
            bp["caps"][2 * i].set_color(edge)
            bp["caps"][2 * i + 1].set_color(edge)
            bp["medians"][i].set_color(edge)

        # Individual data points overlaid on top of each box, jittered
        # sideways (uniform +/-JITTER_RANGE around the box's own x
        # position, which is i+1 - matplotlib's default boxplot() x-tick
        # position for the (i+1)'th box, since `positions=` isn't passed
        # above) so points from the same group don't all stack in one
        # vertical line. zorder=3 keeps them drawn on TOP of the box
        # (patch_artist boxes default to a lower zorder) - same darker
        # per-group colour as the box's own edge, so a point still reads
        # as "which group" at a glance, rather than flat black/gray that
        # would look the same for every group.
        rng = np.random.default_rng()
        for i, (values, edge, label, methods, rows_for_box) in enumerate(
                zip(box_data, edge_colors, box_labels, box_methods, box_rows)):
            jitter = rng.uniform(-JITTER_RANGE, JITTER_RANGE, size=len(values))
            x_positions = np.full(len(values), i + 1) + jitter
            # Jitter point colour: darken_color() applied to this group's
            # own FILL colour (box_colors[i]) at JITTER_DARKEN_FACTOR,
            # computed explicitly here rather than reusing `edge` (the box
            # edge/whisker/cap/median colour, still darken_color()'s
            # default factor) - independent from both that AND
            # LABEL_DARKEN_FACTOR (point labels' own colour, below).
            jitter_color = darken_color(box_colors[i], factor=JITTER_DARKEN_FACTOR)
            ax.scatter(x_positions, values,
                       s=JITTER_POINT_SIZE, alpha=JITTER_POINT_ALPHA,
                       color=jitter_color, edgecolors="none", zorder=3)

            # Point labels (SHOW_POINT_LABELS): the specific varying
            # parameter that distinguishes this point from its group-mates
            # (e.g. "r5", "v04") - purely cosmetic, so a group with no
            # POINT_LABEL_PATTERNS entry, or a member whose
            # shorten_method_label(method) doesn't match that entry's
            # regex, simply gets no label here (no warning - this isn't a
            # data-integrity check like GROUP_RULES itself).
            #
            # xytext's vertical component is staggered by point INDEX j
            # within this group's own point list (values/methods, iterated
            # in the SAME order box_data/box_methods were built in above -
            # plain list order, never re-sorted by value - so the same
            # method always lands at the same stagger position run-to-run)
            # rather than a flat (3, 3) for every point: horizontal jitter
            # alone doesn't separate labels when the points' actual VALUES
            # are near-identical (e.g. Trunk volume's AdTree raw group),
            # since the label offset is relative to each point's own (x, y)
            # - a vertical stagger makes those labels fan out into a small
            # staircase instead of landing on top of each other. With many
            # points in one group (more than ~6-8) this could push labels
            # quite far vertically - acceptable for today's group sizes
            # (4-7 points), revisit if larger groups make this unreadable.
            # Label TEXT colour: darken_color() applied to this group's own
            # FILL colour (box_colors[i], NOT `edge`) at LABEL_DARKEN_FACTOR
            # - independent from `edge` (the box edge/whisker/cap/median/
            # jitter-point colour, darken_color()'s own default factor)
            # so label darkness can be tuned on its own without touching
            # any of those other elements.
            label_color = darken_color(box_colors[i], factor=LABEL_DARKEN_FACTOR)

            if SHOW_POINT_LABELS and group_kind.get(label) == "treeqsm":
                # TreeQSM-origin group: label built directly from
                # TREEQSM_VARY_BY's actual value(s) on this point's own row
                # (treeqsm_param_tokens(), same formatting as the box
                # label) - no regex needed, we already have the real value.
                for j, (x, y, row) in enumerate(zip(x_positions, values, rows_for_box)):
                    point_label = "/".join(treeqsm_param_tokens(row, TREEQSM_VARY_BY))
                    if point_label:
                        ax.annotate(point_label, xy=(x, y),
                                    xytext=(3, 3 + LABEL_VERTICAL_STAGGER * j),
                                    textcoords="offset points",
                                    fontsize=POINT_LABEL_FONTSIZE, ha="left", va="bottom",
                                    color=label_color, zorder=4)
            elif SHOW_POINT_LABELS and label in POINT_LABEL_PATTERNS:
                pattern, fmt = POINT_LABEL_PATTERNS[label]
                for j, (x, y, method) in enumerate(zip(x_positions, values, methods)):
                    m = re.search(pattern, shorten_method_label(method))
                    if m:
                        ax.annotate(fmt % m.group(1), xy=(x, y),
                                    xytext=(3, 3 + LABEL_VERTICAL_STAGGER * j),
                                    textcoords="offset points",
                                    fontsize=POINT_LABEL_FONTSIZE, ha="left", va="bottom",
                                    color=label_color, zorder=4)

        ax.set_xticklabels(box_labels, rotation=BOX_LABEL_ROTATION, ha="right", fontsize=BOX_LABEL_FONTSIZE)

        # Capture axhline()'s own Line2D return value for each line actually
        # drawn, so the legend below can be built from real handles instead
        # of re-deriving color/label - and so the legend is only added when
        # at least one of the two lines exists for THIS field (skipped
        # entirely otherwise, rather than showing an empty/misleading box).
        ref_line = None
        if SHOW_REFERENCE_LINE and ref_row is not None and ref_row[field_key] is not None:
            ref_line = ax.axhline(ref_row[field_key], linestyle="--", color=reference_line_color,
                                   linewidth=1.5, label="Reference (destructive)")

        treeqsm_line = None
        if SHOW_TREEQSM_REF_LINE and treeqsm_ref_row is not None and treeqsm_ref_row[field_key] is not None:
            treeqsm_line = ax.axhline(treeqsm_ref_row[field_key], linestyle="--", color=TREEQSM_REF_LINE_COLOR,
                                       linewidth=1.5, label=treeqsm_ref_row["method"])

        ref_line_handles = [h for h in (ref_line, treeqsm_line) if h is not None]
        if ref_line_handles:
            # loc="best" (not a fixed corner like "upper right") - matplotlib
            # picks whichever corner overlaps the LEAST with what's already
            # drawn (boxes/whiskers/points), since those are already on the
            # axes by this point - safer than hardcoding a corner on a panel
            # with many groups (boxes can fill the full width, including the
            # right edge). Check a busy panel (e.g. "Total volume") for
            # placement that still looks awkward despite this.
            ax.legend(handles=ref_line_handles, fontsize=BOX_LABEL_FONTSIZE, loc="best")

        # Top headroom (BOX_TOP_MARGIN): applied AFTER everything for this
        # panel is drawn (boxes, whiskers, jittered points, reference line),
        # so ax.get_ylim() below reflects the PANEL'S OWN fully-autoscaled
        # data range. NOTE: matplotlib's autoscaling does NOT factor in
        # ax.annotate() text extents (only actual plotted data - boxes/
        # whiskers/scatter points - drives autoscale), so this margin is a
        # data-range-based approximation, same spirit as plot_volumes.py's
        # own TOP_MARGIN_BASE for its bar-chart annotations - not a literal
        # "measure the label's pixel height" calculation. Computed fresh
        # per panel (this is still inside the per-field
        # `for ax, ... in zip(axes.flat, FIELDS):` loop), so each of the 9
        # subplots gets its own independent headroom based on its own data
        # range, not one shared ylim copied across panels.
        y_min, y_max = ax.get_ylim()
        ax.set_ylim(y_min, y_max + BOX_TOP_MARGIN * (y_max - y_min))

    filter_label = ("full reconstruction, branch_filter='none'" if branch_filter == "none"
                     else "diameter >= 10 cm only, branch_filter='10cm'")
    fig.suptitle("Group comparison: %s  (%s)" % (tree, filter_label), fontsize=14)

    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    fig.subplots_adjust(bottom=BOX_BOTTOM_MARGIN)

    if SAVE_PLOT_PNG:
        out_path = os.path.join(ensure_plots_dir(), "%s_boxplot_%s.png" % (tree, branch_filter))
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        print("Saved:", out_path)

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close(fig)


# =========================  RUN  =====================================
if __name__ == "__main__":
    all_rows = load_results(RESULTS_CSV)
    build_boxplot_figure(all_rows, SELECT_TREE, BRANCH_FILTER)
