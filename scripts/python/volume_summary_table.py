r"""
volume_summary_table.py
------------------------
READ-ONLY with respect to volume_results.csv. Builds one settled summary
table PER VALUE COLUMN (total/trunk/branch m3): one row per tree, one column
per reconstruction method, using the SETTLED model configuration decided in
handoff_2026-09-10_part2.md section 1.

The settled configuration is:

    AdQSM         the selected variant for that tree (NOT the 999 median)
    AdTree raw    r5mm, -wlen, seg50-500-k80
    AdTree calib  r5mm, -wlen, seg50-500-k80, regression-perorder,
                  calibrated against the selected AdQSM variant
    TreeQSM       pd2min = 0.002, mo9, sr010, ri1
                  (both the Optimal and the Simplified stage are reported)

WHY THERE IS NO MEDIAN HERE:
an earlier version of this table took a median across the SEG_LEN grid for
AdTree and across the manual parameter sweep for TreeQSM. Under the settled
configuration every filter resolves to exactly ONE row per tree, so there is
nothing to take a median of. The script checks this and complains loudly if a
filter ever returns zero rows or more than one - that would mean the settled
configuration no longer matches what is in the CSV.

WHY THE TreeQSM TAG IS A HARDCODED TABLE:
the tag encodes per-tree reconstruction parameters (pd04, pd045, pd05, ...)
which differ from tree to tree and cannot be derived from anything else in the
CSV. It is also NOT unique across trees - pd05-002-09 belongs to S08, S10 and
S14 alike - so the tree must be carried alongside it. The table below is
copied from the handoff and must be updated by hand if the selection changes.

WHY THE AdQSM VARIANT IS NOT HARDCODED:
each tree has exactly two AdQSM rows: the 999 synthetic median variant and one
real selected variant. "Not 999" therefore identifies the selection uniquely,
and the script derives it rather than relying on a second hand-maintained list
that could drift out of step with the first.

WHY total/trunk/branch SHARE ONE SELECTION PASS (not three):
total_m3, trunk_m3 and branch_m3 are three columns of the SAME settled row in
volume_results.csv - the settled AdTree/TreeQSM/AdQSM filters do not depend on
which of the three is being read. So the matching row is looked up ONCE per
(tree, method) via exactly_one_row(), and all three value columns are read off
that SAME row. Filtering separately per value column would let a bug silently
pick a DIFFERENT row for trunk than for total; reading one shared row makes
that impossible by construction (see exactly_one_row()'s own docstring).

How to use it:
  - Check the paths in CONFIGURATION, then run:
        python volume_summary_table.py
  - Read the printed tables and the checks above/below them before using the
    three output CSVs.
"""

from pathlib import Path

import pandas as pd

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# The r"" prefix means "raw string": it stops Python from treating the
# backslashes in a Windows path as escape characters.
WORK_FOLDER = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python"

INPUT_CSV = Path(WORK_FOLDER) / "volume_results.csv"
CSV_OUT_DIR = Path(WORK_FOLDER) / "csv"

# Which volume columns to summarise - one output CSV per entry, in this
# order. "total_m3" is the whole tree; "trunk_m3"/"branch_m3" split it.
# This is the ONE place these three names are written down - everything
# below reads this list, never the column names themselves.
VALUE_COLUMNS = ["total_m3", "trunk_m3", "branch_m3"]

# Consistency-check tolerance (Step 6): trunk_m3 + branch_m3 vs total_m3,
# as a percentage of total_m3. Cells over this get flagged, not explained.
CONSISTENCY_TOLERANCE_PCT = 0.5

DECIMAL_PLACES = 3

# --- Settled AdTree configuration ----------------------------------------
ADTREE_RADIUS_THRESHOLD_MM = 5
ADTREE_SEG_MIN_MM = 50
ADTREE_SEG_MAX_MM = 500
ADTREE_SEG_K_PCT = 80
ADTREE_WEIGHTING_SUFFIX = "-wlen"      # length-weighted radius statistics
ADTREE_CALMETHOD = "regression-perorder"

# --- Settled TreeQSM configuration ---------------------------------------
# The per-tree parameter tag, from handoff_2026-09-10_part2.md section 1.
# The suffix below is appended to each tag to form the full parameter set.
TREEQSM_TAG_SUFFIX = "_mo9_sr010_ri1"

TREEQSM_TAG_BY_TREE = {
    "B21_S12": "man_pd04-002-05",
    "B21_S16": "man_pd045-002-05",
    "B21_S09": "man_pd05-002-08",
    "B21_S10": "man_pd05-002-09",
    "B21_S08": "man_pd05-002-09",
    "B21_S14": "man_pd05-002-09",
    "B21_S04": "man_pd05-002-10",
    "B21_S11": "man_pd06-002-11",
    "B21_S13": "man_pd06-002-10",
    "B21_S21": "man_pd066-002-099",
    "B21_S01": "man_pd07-002-10",
    "B21_S07": "man_pd08-002-14",
}

# The methods, in the fixed column order every output table uses (also the
# order the consistency check walks in Step 6).
METHOD_NAMES = ["adqsm", "adtree_raw", "adtree_calibrated", "treeqsm_optimal", "treeqsm_simplified"]


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def output_csv_path(value_column: str) -> Path:
    """
    csv/volume_summary_settled_<short>.csv, where <short> is `value_column`
    with its "_m3" suffix stripped (total_m3 -> total, trunk_m3 -> trunk,
    branch_m3 -> branch). Derived from VALUE_COLUMNS itself, not a second
    hardcoded name list - see this module's own docstring for why the three
    names are written down in exactly one place.
    """
    short = value_column[:-len("_m3")] if value_column.endswith("_m3") else value_column
    return CSV_OUT_DIR / f"volume_summary_settled_{short}.csv"


def selected_adqsm_variant_by_tree(df: pd.DataFrame) -> dict:
    """
    Return {tree: variant string}, e.g. {"B21_S01": "080", ...}.

    Every tree has exactly two AdQSM rows: variant 999, which is the
    synthetic median built from the 999 variant grid, and one real
    selected variant. The selected one is simply "the one that is not
    999", so it does not need a second hand-maintained lookup table.

    The variant number lives inside the method string, e.g.
    "AdQSM (TreesParams) (AdQSM 080)", so it is pulled out with a regular
    expression rather than read from the adqsm_variant column, which is
    empty on these rows.
    """
    adqsm_rows = df[df["method"].str.startswith("AdQSM (TreesParams)")].copy()
    adqsm_rows["variant"] = adqsm_rows["method"].str.extract(r"AdQSM (\d+)\)")

    selected = {}
    for tree, group in adqsm_rows.groupby("tree"):
        non_median = sorted(v for v in group["variant"].dropna().unique() if v != "999")
        if len(non_median) == 1:
            selected[tree] = non_median[0]
        else:
            print(f"WARNING: tree {tree} has {len(non_median)} non-999 AdQSM "
                  f"variant(s) {non_median} - expected exactly 1. This tree will "
                  f"be left blank in the AdQSM columns.")
    return selected


def exactly_one_row(matched_rows: pd.DataFrame, tree: str, label: str):
    """
    Return the single row in `matched_rows` (as a pandas Series), or None
    with a printed warning if there is not exactly one.

    Selected ONCE per (tree, method) - not once per value column - so
    total_m3/trunk_m3/branch_m3 for a given cell all come from the SAME
    matched row (see this module's own docstring, "WHY total/trunk/branch
    SHARE ONE SELECTION PASS").

    Under the settled configuration each filter must select one model, so
    zero rows means the settled combination is missing from the CSV, and
    more than one means the filter is not specific enough. Both are
    reported rather than silently averaged away - averaging would hide the
    fact that the configuration and the data no longer agree.
    """
    if len(matched_rows) == 1:
        return matched_rows.iloc[0]
    if len(matched_rows) == 0:
        print(f"WARNING: {tree} / {label}: no row matches the settled "
              f"configuration - leaving blank.")
    else:
        preview = matched_rows[VALUE_COLUMNS].round(3).to_dict("records")
        print(f"WARNING: {tree} / {label}: {len(matched_rows)} rows match the settled "
              f"configuration, expected 1 - leaving blank. Values: {preview}")
    return None


def values_from_row(row, tree: str, label: str) -> dict:
    """
    Read every column in VALUE_COLUMNS off one matched row (from
    exactly_one_row()), each rounded to DECIMAL_PLACES. If `row` is None
    (selection failed), every value column comes back None - a failed
    selection leaves ALL THREE tables blank for that cell, not just one,
    since there is no valid row to read any of them from.
    """
    if row is None:
        return {col: None for col in VALUE_COLUMNS}
    return {col: round(float(row[col]), DECIMAL_PLACES) for col in VALUE_COLUMNS}


def build_summaries(df: pd.DataFrame) -> dict:
    """
    Build one settled summary table per entry in VALUE_COLUMNS.

    Returns {value_column: DataFrame}, each DataFrame with the same column
    layout as before: tree, adqsm_variant, treeqsm_tag, adqsm, adtree_raw,
    adtree_calibrated, treeqsm_optimal, treeqsm_simplified.

    The matching row for each (tree, method) is found exactly ONCE (via
    exactly_one_row()) and then every value column is read off that SAME
    row (via values_from_row()) - see this module's own docstring.
    """
    selected_variants = selected_adqsm_variant_by_tree(df)
    rows_by_value_column = {col: [] for col in VALUE_COLUMNS}

    for tree in sorted(TREEQSM_TAG_BY_TREE):
        tree_rows = df[df["tree"] == tree]
        if len(tree_rows) == 0:
            print(f"WARNING: tree {tree} is in TREEQSM_TAG_BY_TREE but has no rows "
                  f"in the CSV at all - skipping it.")
            continue

        method = tree_rows["method"]
        variant = selected_variants.get(tree)

        # --- AdTree geometry filter, shared by raw and calibrated --------
        # These four settings together define the settled AdTree model.
        adtree_settings = (
            (tree_rows["radius_threshold_mm"] == ADTREE_RADIUS_THRESHOLD_MM)
            & (tree_rows["seg_min_mm"] == ADTREE_SEG_MIN_MM)
            & (tree_rows["seg_max_mm"] == ADTREE_SEG_MAX_MM)
            & (tree_rows["seg_k_pct"] == ADTREE_SEG_K_PCT)
        )
        is_weighted = method.str.endswith(ADTREE_WEIGHTING_SUFFIX)

        # --- AdQSM: the selected variant ---------------------------------
        if variant is None:
            adqsm_row = None
        else:
            adqsm_row = exactly_one_row(
                tree_rows[method.str.startswith("AdQSM (TreesParams)")
                          & method.str.contains(f"AdQSM {variant})", regex=False)],
                tree, f"AdQSM {variant}")

        # --- AdTree raw ---------------------------------------------------
        raw_row = exactly_one_row(
            tree_rows[method.str.startswith("AdTree raw") & is_weighted & adtree_settings],
            tree, "AdTree raw")

        # --- AdTree calibrated --------------------------------------------
        # Calibrated rows additionally carry the calibration method and the
        # AdQSM variant they were calibrated against.
        if variant is None:
            calib_row = None
        else:
            calib_row = exactly_one_row(
                tree_rows[method.str.startswith("AdTree calibrated")
                          & is_weighted
                          & adtree_settings
                          & (tree_rows["calmethod"] == ADTREE_CALMETHOD)
                          & (tree_rows["adqsm_variant"] == int(variant))],
                tree, "AdTree calibrated")

        # --- TreeQSM, both stages of the same parameter set ----------------
        parameter_set = TREEQSM_TAG_BY_TREE[tree] + TREEQSM_TAG_SUFFIX
        optimal_row = exactly_one_row(
            tree_rows[method.str.contains(f"({parameter_set}, Optimal)", regex=False)],
            tree, f"TreeQSM Optimal {parameter_set}")
        simplified_row = exactly_one_row(
            tree_rows[method.str.contains(f"({parameter_set}, Simplified", regex=False)],
            tree, f"TreeQSM Simplified {parameter_set}")

        matched_row_by_method = {
            "adqsm": adqsm_row,
            "adtree_raw": raw_row,
            "adtree_calibrated": calib_row,
            "treeqsm_optimal": optimal_row,
            "treeqsm_simplified": simplified_row,
        }

        # One values_from_row() call per method (NOT per value column) -
        # each call reads all of VALUE_COLUMNS off that method's ONE
        # matched row at once.
        values_by_method = {
            method_name: values_from_row(matched_row, tree, method_name)
            for method_name, matched_row in matched_row_by_method.items()
        }

        for value_col in VALUE_COLUMNS:
            row_dict = {"tree": tree, "adqsm_variant": variant, "treeqsm_tag": parameter_set}
            for method_name in METHOD_NAMES:
                row_dict[method_name] = values_by_method[method_name][value_col]
            rows_by_value_column[value_col].append(row_dict)

    return {col: pd.DataFrame(rows) for col, rows in rows_by_value_column.items()}


def print_and_report_missing(value_column: str, summary: pd.DataFrame):
    """Print one summary table under a heading, then report how many
    cells (if any) could not be filled - same check as before, run once
    per value column."""
    print(f"\nSettled volume summary ({value_column}), m3:")
    print(summary.to_string(index=False))

    n_missing = summary[METHOD_NAMES].isna().sum().sum()
    if n_missing:
        print(f"\n{n_missing} cell(s) could not be filled - see the warnings above.")
    else:
        print(f"\nAll {len(summary)} tree(s) resolved to exactly one model per method.")


def print_consistency_check(summaries: dict):
    """
    Step 6: for every tree and method, compare trunk_m3 + branch_m3 against
    total_m3 - printed as an absolute difference and as a percentage of
    total_m3, for every cell, then a separate flagged list of every cell
    over CONSISTENCY_TOLERANCE_PCT. Cells where any of the three values is
    missing (a warning already fired above) are skipped and counted, not
    silently treated as zero.
    """
    total_df = summaries["total_m3"].set_index("tree")
    trunk_df = summaries["trunk_m3"].set_index("tree")
    branch_df = summaries["branch_m3"].set_index("tree")

    print("\n" + "=" * 100)
    print(f"CONSISTENCY CHECK: trunk_m3 + branch_m3 vs total_m3 (tolerance {CONSISTENCY_TOLERANCE_PCT} %)")
    print("=" * 100)

    check_rows = []
    skipped = []
    for tree in total_df.index:
        for method in METHOD_NAMES:
            total_v = total_df.loc[tree, method]
            trunk_v = trunk_df.loc[tree, method]
            branch_v = branch_df.loc[tree, method]
            if pd.isna(total_v) or pd.isna(trunk_v) or pd.isna(branch_v):
                skipped.append((tree, method))
                continue
            sum_parts = trunk_v + branch_v
            abs_diff = abs(sum_parts - total_v)
            pct_diff = (100.0 * abs_diff / total_v) if total_v else float("nan")
            check_rows.append({
                "tree": tree, "method": method,
                "trunk_m3": trunk_v, "branch_m3": branch_v,
                "trunk_plus_branch": round(sum_parts, DECIMAL_PLACES),
                "total_m3": total_v,
                "abs_diff": round(abs_diff, DECIMAL_PLACES),
                "pct_of_total": round(pct_diff, 3),
            })

    check_df = pd.DataFrame(check_rows)
    print(check_df.to_string(index=False))

    if skipped:
        print(f"\n{len(skipped)} tree/method cell(s) skipped (missing trunk/branch/total value "
              f"- see warnings above): {skipped}")

    flagged = check_df[check_df["pct_of_total"] > CONSISTENCY_TOLERANCE_PCT]
    print()
    if len(flagged):
        print(f"{len(flagged)} cell(s) exceed {CONSISTENCY_TOLERANCE_PCT} %:")
        print(flagged.to_string(index=False))
    else:
        print(f"No cell exceeds {CONSISTENCY_TOLERANCE_PCT} %.")


def main():
    print(f"Loading: {INPUT_CSV}")
    if not INPUT_CSV.is_file():
        print(f"ERROR: file not found: {INPUT_CSV}")
        return
    # low_memory=False stops pandas guessing column types from the first
    # chunk only, which produces mixed-type warnings on this file.
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    print(f"{len(df)} row(s), {df['tree'].nunique()} tree(s)\n")

    summaries = build_summaries(df)

    for value_column in VALUE_COLUMNS:
        print_and_report_missing(value_column, summaries[value_column])

    print_consistency_check(summaries)

    print()
    CSV_OUT_DIR.mkdir(parents=True, exist_ok=True)
    for value_column in VALUE_COLUMNS:
        out_path = output_csv_path(value_column)
        summaries[value_column].to_csv(out_path, index=False)
        print(f"Saved to: {out_path}")


if __name__ == "__main__":
    main()
