# -*- coding: utf-8 -*-
# =====================================================================
#  Compare tree-volume results from different methods against a reference.
# ---------------------------------------------------------------------
#  All results live in ONE simple "master" CSV (RESULTS_CSV below), one row
#  per (tree, method). Columns:
#
#     tree, method, total_m3, trunk_m3, branch_m3, std_m3, dbh_m, height_m,
#     taper_cm_per_m, trunk_len_m, branch_len_m, branch_filter
#
#  - tree           : tree ID, e.g. "IND01_054"
#  - method         : a free-text label, e.g. "AdQSM", "TreeQSM mine v2"
#  - total_m3       : total wood volume [m^3]
#  - trunk_m3       : trunk volume [m^3]   (may be blank)
#  - branch_m3      : branch volume [m^3]        (may be blank)
#  - std_m3         : standard deviation of the total, if known (may be blank)
#  - dbh_m          : stem diameter at 1.3 m above the tree base [m] (may be blank)
#  - height_m       : total tree height [m]                          (may be blank)
#  - taper_cm_per_m : stem taper between two reference heights [cm/m] (may be blank)
#  - trunk_len_m    : total trunk/stem centerline length [m]          (may be blank)
#  - branch_len_m   : total branch centerline length [m]              (may be blank)
#  - branch_filter  : "none" (full/unfiltered reconstruction) or "10cm"
#                     (trunk/branches restricted to diameter >= 10 cm).
#                     Blank/missing (older rows) is treated as "none".
#
#  trunk_len_m/branch_len_m exist to tell apart TWO different reasons a
#  method could report less volume than another: either it reconstructs a
#  shorter/less-complete branch structure (visible as a length difference),
#  or it reconstructs the SAME length but with systematically different
#  radii (length matches, volume doesn't) - see field_error_summary() calls
#  for "Trunk length"/"Branch length" at the bottom of this file.
#
#  WHY branch_filter EXISTS - the destructive field reference physically
#  never measured branches thinner than a 10 cm taper diameter (de Tanago
#  methodology, see AdQSM.pdf Appendix A). That makes a "vs. reference"
#  comparison using a method's FULL (unfiltered) reconstruction methodologically
#  unfair to that method - it's being penalised for wood the reference never
#  even tried to measure. So this script prints TWO separate sections:
#
#    A) "vs. reference" - only branch_filter == "10cm" rows (the reference
#       always has one; this is the fair, apples-to-apples accuracy check).
#    B) "methods vs. each other" - only branch_filter == "none" rows (full
#       reconstructions, no reference present at all) - here you're comparing
#       what each method reconstructs INCLUDING thin branches, purely against
#       each other, not against any "truth".
#
#  This script:
#     1) reads that CSV,
#     2) for the chosen tree (or ALL trees) prints every method side by side,
#        with its difference from the REFERENCE method (absolute and %),
#     3) across all trees, computes error metrics of every method vs the
#        reference: Bias, MAE, RMSE and CV-RMSE.
#
#  To ADD ANOTHER TREE later: just add more rows to the CSV (same columns)
#  and re-run. Nothing else changes.
#
#  If RESULTS_CSV does not exist yet, this script writes a starter file
#  pre-filled with the IND01_054 numbers we already have, so you can run it
#  immediately and then edit/extend it.
#
#  Dependencies: none (standard library only).
# =====================================================================

import csv
import os
import math

# =====================  PARAMETERS  ==================================
# The master results table. Edit this file to add trees/methods.
RESULTS_CSV = "volume_results.csv"

# Which tree to show in detail: a specific ID like "IND01_054", or "ALL".
SELECT_TREE = "ALL"

# The method used as the "truth" that everything is compared against. Must
# match a value in the 'method' column exactly.
REFERENCE_METHOD = "Reference (destructive)"

# The "reference" for MODE B (branch_filter == "none", methods compared to
# EACH OTHER - see the header comment above). The destructive field
# reference can never appear in that mode (it only ever has a "10cm" row),
# so there is no absolute "truth" to compare against there - instead we pick
# AdQSM's own reported numbers as a common yardstick, purely so the
# Bias/MAE/RMSE/field-error tables have SOME baseline to express every other
# method's difference against. This is NOT an accuracy claim about AdQSM -
# it's just "how far is each method from AdQSM", nothing more.
REFERENCE_METHOD_NONE = "AdQSM (TreesParams) (AdQSM 05)"
# =====================================================================


# Starter content written only if RESULTS_CSV does not exist yet. The two
# length columns (trunk_len_m, branch_len_m) are left blank here since this
# starter data predates them - real runs of the other scripts fill them in.
# branch_filter is filled in per the same rule the real scripts use: "10cm"
# for the destructive reference (it can only ever be that), "none" for
# everything else here (none of these starter rows are a "Filtered"/
# "(>=10cm only)" variant).
# n_cylinders (Task A) added as the LAST column, same as the real
# upsert_result() writers - starter rows predate it too, so it's left
# blank ("") here, same pattern already used for trunk_len_m/branch_len_m.
STARTER_ROWS = [
    ["tree", "method", "total_m3", "trunk_m3", "branch_m3", "std_m3",
     "dbh_m", "height_m", "taper_cm_per_m", "trunk_len_m", "branch_len_m", "branch_filter", "n_cylinders"],
    ["IND01_054", "Reference (destructive)", "1.7169", "1.2557", "0.4612", "", "", "", "", "", "", "10cm", ""],
    ["IND01_054", "AdQSM (TreesParams)",     "1.6905", "1.1302", "0.5603", "", "", "", "", "", "", "none", ""],
    ["IND01_054", "AdTree calibrated",       "1.9240", "1.3020", "0.6220", "", "", "", "", "", "", "none", ""],
    ["IND01_054", "TreeQSM de Tanago (mean 20)", "3.0838", "1.8383", "1.2455", "0.3113", "", "", "", "", "", "none", ""],
    ["IND01_054", "TreeQSM mine v1",         "3.6806", "1.5317", "2.1488", "0.1915", "", "", "", "", "", "none", ""],
    ["IND01_054", "TreeQSM mine v2",         "3.9744", "1.5737", "2.4007", "0.1923", "", "", "", "", "", "none", ""],
]


def to_float(text):
    """Parse a number, or return None for blank/non-numeric cells."""
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def load_results(path):
    """Read the master CSV into a list of dict rows with parsed numbers."""
    rows = []
    with open(path, "r", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "tree": r["tree"].strip(),
                "method": r["method"].strip(),
                "total": to_float(r.get("total_m3")),
                "trunk": to_float(r.get("trunk_m3")),
                "branch": to_float(r.get("branch_m3")),
                "std": to_float(r.get("std_m3")),
                "dbh": to_float(r.get("dbh_m")),
                "height": to_float(r.get("height_m")),
                "taper": to_float(r.get("taper_cm_per_m")),
                "trunk_len": to_float(r.get("trunk_len_m")),
                "branch_len": to_float(r.get("branch_len_m")),
                # branch_filter: blank cell OR the column missing entirely
                # (an older row, from before this column existed) both fall
                # back to "none" - the "or" chain handles both None and "".
                "branch_filter": (r.get("branch_filter") or "").strip() or "none",
                # n_cylinders (Task A): same to_float() helper as every other
                # numeric column - already returns None for blank/missing
                # cells (e.g. rows from before this column existed).
                "n_cylinders": to_float(r.get("n_cylinders")),
                # Reconstruction parameters (mode/PD1/PD2Min/PD2Max/
                # MinCylRad/simp_*) - see runsken.m section 19's
                # params_<tree>_<run>.csv and import_matlab_results.py's
                # read_params_file()/upsert_result(). "mode" follows
                # branch_filter's own blank-or-missing-column convention
                # (falls back to "" here, not "none" - there's no
                # meaningful default mode the way "none" is a meaningful
                # default branch_filter); every numeric one uses the same
                # to_float()-returns-None-if-missing pattern as every other
                # optional column above (AdTree rows, and any row from
                # before this column existed, simply get None here).
                "mode": (r.get("mode") or "").strip(),
                "pd1": to_float(r.get("pd1_m")),
                "pd2min": to_float(r.get("pd2min_m")),
                "pd2max": to_float(r.get("pd2max_m")),
                "mincylrad": to_float(r.get("mincylrad_m")),
                "simp_maxorder": to_float(r.get("simp_maxorder")),
                "simp_smallradii": to_float(r.get("simp_smallradii")),
                "simp_replaceiterations": to_float(r.get("simp_replaceiterations")),
            })
    return rows


def fmt(x, width=10, dec=4):
    """Format a number or a blank cell to a fixed width."""
    if x is None:
        return " " * (width - 1) + "-"
    return ("%*.*f" % (width, dec, x))


def pct_diff(value, ref_value):
    """Percent difference of value vs. ref_value, or None if either is
    missing (blank cell in the printed table)."""
    if value is None or ref_value is None or ref_value == 0:
        return None
    return (value - ref_value) / ref_value * 100.0


def fmt_pct(x, width=8):
    """Format a percent difference, or a blank cell if None."""
    if x is None:
        return " " * (width - 1) + "-"
    return "%+*.1f%%" % (width - 1, x)


def print_tree_block(tree, rows, no_reference_note=None):
    """Print the side-by-side method comparison for a single tree.

    `no_reference_note`, if given, REPLACES the generic "(no reference
    method...)" line printed when this tree has no REFERENCE_METHOD row in
    `rows` - used by the "methods vs. each other" section in the RUN part
    below, where that's expected (not a data problem), so a clearer,
    section-specific note is shown instead. Purely a wording change - the
    comparison math above is identical either way."""
    tree_rows = [r for r in rows if r["tree"] == tree]
    ref = next((r for r in tree_rows if r["method"] == REFERENCE_METHOD), None)

    print("=" * 118)
    print("Tree: %s" % tree)
    print("-" * 118)
    print("%-28s %10s %10s %10s %10s   %8s %8s   %8s %8s   %8s %8s" %
          ("method", "total", "trunk", "branch", "d(total)",
           "DBH[m]", "d(DBH)", "H[m]", "d(H)", "taper", "d(taper)"))
    print("-" * 118)

    # reference first (if present), then the rest in file order
    ordered = ([ref] if ref else []) + [r for r in tree_rows if r is not ref]
    for r in ordered:
        if r is None:
            continue
        if ref and ref["total"] is not None and r["total"] is not None:
            d_abs = r["total"] - ref["total"]
            d_pct = d_abs / ref["total"] * 100.0 if ref["total"] else 0.0
            if r is ref:
                dcol = "reference"
            else:
                dcol = "%+.3f (%+.0f%%)" % (d_abs, d_pct)
        else:
            dcol = "-"

        dbh_pct = None if r is ref else pct_diff(r["dbh"], ref["dbh"] if ref else None)
        height_pct = None if r is ref else pct_diff(r["height"], ref["height"] if ref else None)
        taper_pct = None if r is ref else pct_diff(r["taper"], ref["taper"] if ref else None)

        print("%-28s %s %s %s   %-16s %s %s   %s %s   %s %s" %
              (r["method"][:28], fmt(r["total"]), fmt(r["trunk"]),
               fmt(r["branch"]), dcol,
               fmt(r["dbh"], width=8, dec=3), fmt_pct(dbh_pct),
               fmt(r["height"], width=8, dec=2), fmt_pct(height_pct),
               fmt(r["taper"], width=8, dec=2), fmt_pct(taper_pct)))
    if ref is None:
        if no_reference_note is not None:
            print(no_reference_note)
        else:
            print("(no reference method '%s' for this tree - showing raw values only)"
                  % REFERENCE_METHOD)
    print()


def compute_error_metrics(rows, reference_method):
    """Same Bias/MAE/RMSE/CV-RMSE calculation as error_metrics() below, but
    RETURNS the numbers instead of printing them (one dict per method), so
    other scripts (e.g. plot_volumes.py) can reuse this exact calculation
    instead of duplicating the math and risking the two going out of sync.
    Only trees where BOTH that method and the reference exist are used.

    `reference_method` is now a REQUIRED parameter (previously this function
    always used the global REFERENCE_METHOD) so the SAME calculation can
    serve both comparison modes further down this file: MODE A passes
    REFERENCE_METHOD (the destructive field reference, branch_filter=="10cm"),
    MODE B passes REFERENCE_METHOD_NONE (AdQSM, branch_filter=="none") -
    without this parameter the two modes would need two near-duplicate
    copies of this function, one per hard-coded reference."""
    trees = sorted({r["tree"] for r in rows})
    methods = [m for m in dict.fromkeys(r["method"] for r in rows)
               if m != reference_method]

    # quick lookup: (tree, method) -> total
    total_of = {(r["tree"], r["method"]): r["total"] for r in rows}

    results = []
    for m in methods:
        pairs = []
        for t in trees:
            est = total_of.get((t, m))
            ref = total_of.get((t, reference_method))
            if est is not None and ref is not None:
                pairs.append((est, ref))
        if not pairs:
            continue
        n = len(pairs)
        errs = [e - r for e, r in pairs]
        refs = [r for _, r in pairs]
        bias = sum(errs) / n
        mae = sum(abs(e) for e in errs) / n
        rmse = math.sqrt(sum(e * e for e in errs) / n)
        cv_rmse = rmse / (sum(refs) / n) * 100.0 if refs else 0.0
        results.append(dict(method=m, n=n, bias=bias, mae=mae, rmse=rmse, cv_rmse=cv_rmse))
    return results


def error_metrics(rows, reference_method):
    """Across all trees, compute (via compute_error_metrics) and PRINT
    Bias/MAE/RMSE/CV-RMSE of each method vs `reference_method`.

    `reference_method` is REQUIRED (see compute_error_metrics() above for
    why) - it also appears literally in the printed header line, so MODE B's
    call (reference_method=REFERENCE_METHOD_NONE) automatically prints
    "Error metrics vs 'AdQSM (TreesParams)'" instead of the destructive
    reference's name, with no separate wording needed."""
    results = compute_error_metrics(rows, reference_method)

    print("=" * 78)
    print("Error metrics vs '%s'  (total volume, across trees)" % reference_method)
    print("-" * 78)
    print("%-28s %5s %9s %9s %9s %9s" %
          ("method", "n", "Bias", "MAE", "RMSE", "CV-RMSE%"))
    print("-" * 78)

    for r in results:
        print("%-28s %5d %9.3f %9.3f %9.3f %9.1f" %
              (r["method"][:35], r["n"], r["bias"], r["mae"], r["rmse"], r["cv_rmse"]))
    print("(Bias/MAE/RMSE in m^3. n = number of trees compared. With one tree,")
    print(" MAE = RMSE = |error| and CV-RMSE is that error relative to the reference.)")
    print()


def field_error_summary(rows, key, label, reference_method):
    """Print one short line per method: its mean percent error vs.
    `reference_method` for a single field (dbh/height/taper/trunk_len/
    branch_len), across all trees where both values are known. Kept OUT of
    the volume RMSE metrics block above by design - these fields are not
    volumes.

    `reference_method` is REQUIRED for the same reason as in
    compute_error_metrics()/error_metrics() above: this one function now
    serves both MODE A (vs. the destructive reference) and MODE B (vs.
    AdQSM) - the caller decides which by passing the right reference_method,
    and the printed header line ("... error vs '<reference_method>':")
    reflects whichever one was actually used."""
    trees = sorted({r["tree"] for r in rows})
    methods = [m for m in dict.fromkeys(r["method"] for r in rows) if m != reference_method]
    value_of = {(r["tree"], r["method"]): r[key] for r in rows}

    print("%s error vs '%s':" % (label, reference_method))
    printed_any = False
    for m in methods:
        diffs = []
        for t in trees:
            d = pct_diff(value_of.get((t, m)), value_of.get((t, reference_method)))
            if d is not None:
                diffs.append(d)
        if not diffs:
            continue
        printed_any = True
        print("  %-28s %+6.1f%%  (n=%d)" % (m[:28], sum(diffs) / len(diffs), len(diffs)))
    if not printed_any:
        print("  (no method has both a value and a reference value for this field)")
    print()


def filter_by_branch_filter(rows, value):
    """Keep only the rows whose branch_filter equals `value` ("none" or
    "10cm"). This is the ONLY new filtering step for the two-mode RUN
    section below - it happens BEFORE rows are handed to print_tree_block/
    error_metrics/field_error_summary, so none of that existing math changes,
    it just runs on a smaller (pre-filtered) list of rows."""
    return [r for r in rows if r["branch_filter"] == value]


# =========================  RUN  =====================================
# Guarded by __name__ == "__main__" so this file can also be IMPORTED (e.g.
# by plot_volumes.py, to reuse load_results/pct_diff/compute_error_metrics)
# without re-running all these prints / writing the starter CSV as a side effect.
if __name__ == "__main__":
    if not os.path.exists(RESULTS_CSV):
        with open(RESULTS_CSV, "w", encoding="utf-8", newline="") as f:
            csv.writer(f).writerows(STARTER_ROWS)
        print("No results file found - wrote a starter '%s' with the IND01_054 data.\n"
              % RESULTS_CSV)

    rows = load_results(RESULTS_CSV)

    # ====================================================================
    # MODE A: comparison AGAINST THE REFERENCE - only branch_filter == "10cm"
    # rows (the destructive reference only ever has a "10cm" row, so this is
    # the fair, apples-to-apples accuracy check - see the header comment for why).
    # ====================================================================
    print("#" * 118)
    print("=== COMPARISON AGAINST REFERENCE (branches/trunks >= 10 cm) ===")
    print("#" * 118)
    print()

    rows_10cm = filter_by_branch_filter(rows, "10cm")
    trees_10cm = sorted({r["tree"] for r in rows_10cm})

    if SELECT_TREE.upper() == "ALL":
        for t in trees_10cm:
            print_tree_block(t, rows_10cm)
    elif SELECT_TREE in trees_10cm:
        print_tree_block(SELECT_TREE, rows_10cm)
    else:
        print("Tree '%s' not found among branch_filter='10cm' rows. Available: %s"
              % (SELECT_TREE, ", ".join(trees_10cm)))

    # Error metrics make sense whenever a reference exists (works for 1 or many trees).
    # reference_method=REFERENCE_METHOD: this is MODE A, so the yardstick is
    # the destructive field reference (the fair, "same 10cm cut-off" comparison).
    error_metrics(rows_10cm, REFERENCE_METHOD)

    # DBH/height/taper are not volumes, so their error is reported separately
    # instead of folding them into the RMSE metrics block above.
    field_error_summary(rows_10cm, "dbh", "DBH", REFERENCE_METHOD)
    field_error_summary(rows_10cm, "height", "Height", REFERENCE_METHOD)
    field_error_summary(rows_10cm, "taper", "Taper", REFERENCE_METHOD)

    # Trunk/branch LENGTH (not volume) - same "one line per method" summary
    # as DBH/height/taper above, not squeezed into print_tree_block()'s
    # already-wide per-tree table (that table is 118 characters wide with 5
    # metric pairs already; two more pairs would make it wrap/unreadable in
    # a normal terminal). Lets you tell apart "shorter/less-complete branch
    # structure" (length is off) from "same length, different radii" (length
    # matches, volume doesn't) - see the header comment at the top of this file.
    field_error_summary(rows_10cm, "trunk_len", "Trunk length", REFERENCE_METHOD)
    field_error_summary(rows_10cm, "branch_len", "Branch length", REFERENCE_METHOD)

    # ====================================================================
    # MODE B: methods COMPARED TO EACH OTHER - only branch_filter == "none"
    # rows (full/unfiltered reconstructions, thin branches included). The
    # destructive reference NEVER appears here (it methodologically can't -
    # it only ever has a "10cm" row), so print_tree_block always finds no
    # reference for this mode; a custom note explains that's expected, not
    # a data problem.
    #
    # error_metrics/field_error_summary are NOW ALSO printed here (this used
    # to be skipped entirely, since the destructive reference can't appear in
    # this mode) - but using REFERENCE_METHOD_NONE (AdQSM) as the yardstick
    # instead of the destructive reference. This is NOT an accuracy
    # evaluation (AdQSM is just another reconstruction, not ground truth) -
    # it only tells you how far each OTHER method's full/unfiltered
    # reconstruction sits from AdQSM's, which is still a useful cross-check
    # even without a "true" reference in this mode.
    # ====================================================================
    print()
    print("#" * 118)
    print("=== COMPARISON OF METHODS (full reconstruction, including thin branches) ===")
    print("#" * 118)
    print()

    rows_none = filter_by_branch_filter(rows, "none")
    trees_none = sorted({r["tree"] for r in rows_none})
    no_ref_note = ("(No destructive-reference row here BY DESIGN - the reference only ever\n"
                   " measured branches >= 10 cm, see the '10cm' section above. This is a pure\n"
                   " method-vs-method comparison, not an accuracy evaluation.)")

    if SELECT_TREE.upper() == "ALL":
        for t in trees_none:
            print_tree_block(t, rows_none, no_reference_note=no_ref_note)
    elif SELECT_TREE in trees_none:
        print_tree_block(SELECT_TREE, rows_none, no_reference_note=no_ref_note)
    else:
        print("Tree '%s' not found among branch_filter='none' rows. Available: %s"
              % (SELECT_TREE, ", ".join(trees_none)))

    # Same error-metrics/field-error blocks as MODE A above, just pointed at
    # REFERENCE_METHOD_NONE (AdQSM) instead of REFERENCE_METHOD (the
    # destructive reference) - reusing the exact same functions, so the
    # calculation can never silently drift between the two modes.
    error_metrics(rows_none, REFERENCE_METHOD_NONE)
    field_error_summary(rows_none, "dbh", "DBH", REFERENCE_METHOD_NONE)
    field_error_summary(rows_none, "height", "Height", REFERENCE_METHOD_NONE)
    field_error_summary(rows_none, "taper", "Taper", REFERENCE_METHOD_NONE)
    field_error_summary(rows_none, "trunk_len", "Trunk length", REFERENCE_METHOD_NONE)
    field_error_summary(rows_none, "branch_len", "Branch length", REFERENCE_METHOD_NONE)
