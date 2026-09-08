# -*- coding: utf-8 -*-
# =====================================================================
#  Build ONE new synthetic "median" AdQSM variant folder for a tree, from
#  a configurable list of already-reconstructed source variants (e.g.
#  B21_S01's HS=0.4..1.0 sweep: "040".."100"), by taking the MEDIAN of
#  each source variant's taper.txt/BranchStructure.txt/TreesParams.txt
#  values, row-by-row (per height / per cylinder / per field).
#
#  WHY row-by-row (not per-order) median is safe here: two prior read-
#  only diagnostics confirmed, for B21_S01's 7-variant sweep specifically:
#    - taper.txt: all 7 variants share an IDENTICAL height grid (28
#      points, 0..27.6 m), no spike rejected in any of them.
#    - BranchStructure.txt: `order`, `parent`, `length(m)`, `height(m)`,
#      `angle(deg)`, `azimuth(deg)`, `zenith(deg)` are BYTE-IDENTICAL,
#      row-by-row, across all 21,263 rows, comparing the most extreme
#      outlier variant against two normal-cluster variants - i.e. row
#      index IS a stable, direct reference to "the same physical
#      cylinder" across variants; only diameter(m)/volume(L)/area(m^2)
#      (all diameter-derived) actually vary.
#  This script does NOT re-run those diagnostics for other trees - it
#  re-checks the SAME two invariants (identical height grid; identical
#  row count) at runtime for whatever tree/variant list it's given, and
#  raises a clear SystemExit if either fails, rather than silently
#  assuming they hold (see build_median_taper()/build_median_branch_
#  structure() below).
#
#  Only CREATES data/<TREE_NAME>/<NEW_VARIANT_NAME>/ with 3 new files
#  inside. Never touches any existing file/folder - refuses to run at
#  all (SystemExit) if that folder already exists, rather than silently
#  overwriting it.
#
#  Dependencies: numpy (install: pip install numpy)
# =====================================================================

import os

import numpy as np

# Reuse (do not re-implement): the same taper.txt parser (with its own
# spike-rejection) and TreesParams.txt parser already used everywhere
# else in this codebase that reads AdQSM output.
from tree_geom_utils import parse_adqsm_taper_file

# =====================  PARAMETERS  ===================================
TREE_NAME = "B21_S04"

# Source AdQSM variant folder names to median together - a list, not a
# range, so it works for any subset/tree without editing the logic
# below. No variant name is hard-coded anywhere past this block.
SOURCE_VARIANTS = ["040", "050", "060", "070", "080", "090", "100"]

# New synthetic variant's folder name. "999" is a plain \d+-matching
# string (confirmed against every AdQSM-variant-number regex used
# elsewhere in this codebase - e.g. reconstruction_method_decision_
# summary.py's _adqsm_variant_of() `r"AdQSM (\d+)\)$"`,
# compare_volumes.py's resolve_reference_method_none() `r"...(\d+)\)$"`,
# taper_curve_compare.py's discover_variants() `sorted(labels, key=int)`
# - all parse "999" the same as any other all-digit folder name) and
# confirmed to NOT collide with any existing variant folder for ANY tree
# currently under DATA_ROOT (checked data/*/  - only B21_S01 uses 3-digit
# folders, "040".."100"; no tree has a "999" folder today).
NEW_VARIANT_NAME = "999"

DATA_ROOT = r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"
# =====================================================================


def _variant_dir(tree_name, variant):
    return os.path.join(DATA_ROOT, tree_name, variant)


# ---------------------------------------------------------------------
# taper.txt
# ---------------------------------------------------------------------

def _taper_data_block(path):
    """Replicate parse_adqsm_taper_file()'s own line-classification state
    machine (tree_geom_utils.py) just far enough to find WHICH raw lines
    are the numeric height/diameter data block in `path` - returns
    (lines, start, stop) where lines[start:stop] is that block (stop
    exclusive), so the surrounding header/footer lines (the "DataName:"
    line, the Chinese column-header line, the two Chinese summary footer
    lines) can be preserved byte-for-byte by the caller. Does NOT do
    spike-rejection or return parsed values itself -
    parse_adqsm_taper_file() (called separately, in build_median_taper()
    below) remains the single source of truth for the actual height/
    diameter numbers used to compute the median; this function is only
    used to locate where in the copied template to splice the new block."""
    with open(path, "r", encoding="latin-1", newline="") as f:
        lines = f.readlines()
    collecting = False
    start = stop = None
    for i, line in enumerate(lines):
        parts = line.strip().split("\t")
        is_data_row = False
        if len(parts) == 2:
            try:
                float(parts[0])
                float(parts[1])
                is_data_row = True
            except ValueError:
                is_data_row = False
        if is_data_row:
            if not collecting:
                start = i
                collecting = True
        elif collecting:
            stop = i
            break
    if collecting and stop is None:
        stop = len(lines)
    if start is None:
        raise ValueError("No numeric height/diameter rows found in %s" % path)
    return lines, start, stop


def build_median_taper(tree_name, source_variants, out_dir):
    """Write out_dir/taper.txt: same header/footer lines as the FIRST
    source variant's taper.txt (copied verbatim), same height grid, but
    diameter at each height replaced by the MEDIAN diameter across all
    source variants at that height.

    Re-verifies (does not assume) that every source variant's height
    grid, as returned by parse_adqsm_taper_file() (i.e. AFTER its own
    spike-rejection), is identical to the first variant's - a
    row-by-row median is only valid when that holds. Raises SystemExit
    with the offending variant name if it doesn't, rather than silently
    guessing/interpolating."""
    per_variant = {}
    for v in source_variants:
        path = os.path.join(_variant_dir(tree_name, v), "taper.txt")
        h, d = parse_adqsm_taper_file(path)
        per_variant[v] = (h, d)

    ref_v = source_variants[0]
    ref_h, _ref_d = per_variant[ref_v]
    for v in source_variants[1:]:
        h, _d = per_variant[v]
        if len(h) != len(ref_h) or not np.allclose(h, ref_h):
            raise SystemExit(
                "adqsm_build_median_variant.py: taper.txt height grid for variant %r does not "
                "match variant %r for tree %r (lengths %d vs %d) - a row-by-row median is only "
                "valid when every source variant shares an identical height grid. Re-run the "
                "height-grid diagnostic for THIS tree/variant list before proceeding."
                % (v, ref_v, tree_name, len(h), len(ref_h)))

    median_diam = np.median(np.stack([per_variant[v][1] for v in source_variants], axis=0), axis=0)

    ref_path = os.path.join(_variant_dir(tree_name, ref_v), "taper.txt")
    lines, start, stop = _taper_data_block(ref_path)

    new_rows = ["%g\t%.6g\r\n" % (h, d) for h, d in zip(ref_h, median_diam)]
    out_lines = lines[:start] + new_rows + lines[stop:]

    out_path = os.path.join(out_dir, "taper.txt")
    with open(out_path, "w", encoding="latin-1", newline="") as f:
        f.writelines(out_lines)

    return ref_h, median_diam, per_variant, out_path


# ---------------------------------------------------------------------
# BranchStructure.txt
# ---------------------------------------------------------------------

def _is_branch_data_row(parts):
    """Same detection tree_geom_utils.py's parse_adqsm_branch_file_raw()
    uses: >=7 tab fields, first field parses as an integer (branch
    order)."""
    if len(parts) < 7:
        return False
    try:
        int(parts[0])
    except ValueError:
        return False
    return True


def build_median_branch_structure(tree_name, source_variants, out_dir):
    """Write out_dir/BranchStructure.txt: every line copied verbatim from
    the FIRST source variant's file EXCEPT each data row's diameter(m)
    column (index 2), which is replaced by the MEDIAN of that SAME ROW
    INDEX's diameter(m) across all source variants.

    order/parent/volume(L)/area(m^2)/length(m)/height(m)/angle(deg)/
    azimuth(deg)/zenith(deg) are left exactly as they were in the
    reference variant's file - in particular volume(L) and area(m^2) are
    NOT recomputed from the new median diameter, so they go stale
    relative to it. This is deliberate, not an oversight:
    parse_adqsm_branch_file()/parse_adqsm_branch_file_raw() (the only
    code in this project that reads BranchStructure.txt downstream) only
    ever read `order` (column 0) and `diameter(m)` (column 2) - confirmed
    in a prior diagnostic - so volume(L)/area(m^2) have no downstream
    effect on anything, and recomputing them would require guessing
    AdQSM's own internal cylinder-volume/area formula without having
    verified it.

    Re-verifies (does not assume) that every source variant has the SAME
    number of data rows as the first - the row-by-row correspondence
    this whole approach depends on was only actually confirmed, in a
    prior diagnostic, for B21_S01's specific 7-variant set."""
    ref_v = source_variants[0]
    ref_path = os.path.join(_variant_dir(tree_name, ref_v), "BranchStructure.txt")
    with open(ref_path, "r", encoding="latin-1", newline="") as f:
        ref_lines = f.readlines()

    per_variant_diam = {}
    n_rows = None
    for v in source_variants:
        path = os.path.join(_variant_dir(tree_name, v), "BranchStructure.txt")
        diam = []
        with open(path, "r", encoding="latin-1", newline="") as f:
            for line in f:
                parts = line.rstrip("\r\n").split("\t")
                if not _is_branch_data_row(parts):
                    continue
                diam.append(float(parts[2]))
        per_variant_diam[v] = diam
        if n_rows is None:
            n_rows = len(diam)
        elif len(diam) != n_rows:
            raise SystemExit(
                "adqsm_build_median_variant.py: BranchStructure.txt row count for variant %r "
                "(%d) does not match variant %r (%d) for tree %r - row-by-row median needs the "
                "SAME cylinder count in every source variant. Re-run the row-correspondence "
                "diagnostic for THIS tree/variant list before proceeding."
                % (v, len(diam), ref_v, n_rows, tree_name))

    median_diam = np.median(np.array([per_variant_diam[v] for v in source_variants]), axis=0)

    out_lines = []
    row_i = 0
    for line in ref_lines:
        parts = line.rstrip("\r\n").split("\t")
        if _is_branch_data_row(parts):
            parts[2] = "%.6g" % median_diam[row_i]
            out_lines.append("\t".join(parts) + "\r\n")
            row_i += 1
        else:
            out_lines.append(line)   # header row / blank separator lines - copied verbatim
    assert row_i == n_rows, "internal error: rewrote %d rows, expected %d" % (row_i, n_rows)

    out_path = os.path.join(out_dir, "BranchStructure.txt")
    with open(out_path, "w", encoding="latin-1", newline="") as f:
        f.writelines(out_lines)

    return ref_lines, per_variant_diam, median_diam, out_path


# ---------------------------------------------------------------------
# TreesParams.txt
# ---------------------------------------------------------------------

def _parse_all_params_fields(path):
    """Parse EVERY numeric 'Key: value' token from a TreesParams.txt
    file's single data line (the one containing 'TreeVolume:'), as an
    ordered {key: float} dict, in the file's own field order.

    parse_adqsm_params_file() (tree_geom_utils.py) is NOT used here: it
    returns a fixed, renamed SUBSET of fields (n/total_len/total_vol/
    trunk_vol/branch_vol/height/dbh/...), shaped for volume_stats()/
    print_volume_stats() - it does not expose CBH, TrunkSurfaceArea,
    BranchesSurfaceArea, TotalArea, CrownWidth, or any of the
    CrownArea*/CrownVolume* fields, and this script needs EVERY field in
    the file median-ed, not that subset (per the task's own requirement).
    This generalizes the exact same "Key: value" regex idiom that
    function's own extract() helper uses to every token actually present
    in the line, instead of a hard-coded field list - so it keeps working
    unchanged if AdQSM's TreesParams.txt export ever gains or loses a
    field, and it does not hard-code TreeHeight/TreeVolume/etc. by name
    anywhere in this function."""
    with open(path, "r", encoding="latin-1", newline="") as f:
        lines = f.readlines()
    for line in lines:
        if "TreeVolume:" not in line:
            continue
        fields = {}
        for tok in line.strip().split("\t"):
            if ":" not in tok:
                continue
            key, _, value = tok.partition(":")
            key = key.strip()
            try:
                fields[key] = float(value.strip())
            except ValueError:
                continue
        if fields:
            return fields
    return None


def build_median_params(tree_name, source_variants, out_dir):
    """Write out_dir/TreesParams.txt: same layout as the source files
    (the "DataName:" line, then one tab-separated "Key: value" data
    line), with EVERY numeric field independently set to its own
    per-field median across the source variants - not a single "closest"
    real variant's row copied wholesale, so this file can't accidentally
    duplicate whatever variant gets separately picked as the "best"
    reconstruction elsewhere in the workflow.

    Because each field is medianed independently (its median can come
    from a DIFFERENT source variant per field), TreeVolume will in
    general NOT exactly equal TrunkVolume + BranchVolume in the output -
    this is expected, not a bug (see run()'s printed check below). It is
    harmless because TreesParams.txt only ever feeds a comparison row in
    volume_results.csv (via report_adqsm_thin_branch()/
    parse_adqsm_params_file() elsewhere) - the actual AdTree calibration
    reads taper.txt and BranchStructure.txt (properly row-by-row medianed
    above), never TreesParams.txt's own volume fields.

    Fields already identical across all source variants (e.g. TreeHeight,
    BranchesNum, CBH, CrownWidth-ish geometry fields - confirmed
    identical for B21_S01's 7-variant set in a prior diagnostic) trivially
    median to that same value - no special-casing needed here.

    Re-verifies (does not assume) that every source variant reports the
    SAME set of field names - raises SystemExit naming the mismatch if
    not, since a per-field median across mismatched field sets isn't
    well-defined."""
    per_variant_fields = {}
    for v in source_variants:
        path = os.path.join(_variant_dir(tree_name, v), "TreesParams.txt")
        fields = _parse_all_params_fields(path)
        if fields is None:
            raise SystemExit("adqsm_build_median_variant.py: no 'TreeVolume:' data line found in %s" % path)
        per_variant_fields[v] = fields

    ref_v = source_variants[0]
    field_order = list(per_variant_fields[ref_v].keys())
    for v in source_variants[1:]:
        missing = [k for k in field_order if k not in per_variant_fields[v]]
        extra = [k for k in per_variant_fields[v] if k not in field_order]
        if missing or extra:
            raise SystemExit(
                "adqsm_build_median_variant.py: TreesParams.txt field set for variant %r differs "
                "from variant %r for tree %r - missing=%s extra=%s. Every source variant must "
                "report the SAME set of fields for a per-field median to be well-defined."
                % (v, ref_v, tree_name, missing, extra))

    median_fields = {k: float(np.median([per_variant_fields[v][k] for v in source_variants]))
                      for k in field_order}

    ref_path = os.path.join(_variant_dir(tree_name, ref_v), "TreesParams.txt")
    with open(ref_path, "r", encoding="latin-1", newline="") as f:
        ref_lines = f.readlines()

    out_lines = []
    for line in ref_lines:
        if "TreeVolume:" not in line:
            out_lines.append(line)   # "DataName:" line, trailing blank line - copied verbatim
            continue
        new_tokens = []
        for tok in line.rstrip("\r\n").split("\t"):
            if ":" not in tok:
                new_tokens.append(tok)   # a stray trailing empty token (the line's own trailing tab)
                continue
            key, _, _value = tok.partition(":")
            key = key.strip()
            if key in median_fields:
                # AdQSM's own field formatting isn't documented/verified, so
                # match the ~6-significant-figure precision already visible
                # in the source files ("%.6g") rather than guessing at a
                # fixed decimal count.
                new_tokens.append("%s: %.6g" % (key, median_fields[key]))
            else:
                new_tokens.append(tok)
        out_lines.append("\t".join(new_tokens) + "\r\n")

    out_path = os.path.join(out_dir, "TreesParams.txt")
    with open(out_path, "w", encoding="latin-1", newline="") as f:
        f.writelines(out_lines)

    return per_variant_fields, median_fields, out_path


def run():
    out_dir = _variant_dir(TREE_NAME, NEW_VARIANT_NAME)
    if os.path.exists(out_dir):
        raise SystemExit(
            "adqsm_build_median_variant.py: %r already exists - refusing to overwrite an "
            "existing folder (this script only ever CREATES a new variant folder). Remove it "
            "yourself first if you really want to regenerate it." % out_dir)
    os.makedirs(out_dir)   # the ONLY folder this script ever creates - nothing else is touched

    print("=" * 90)
    print("Tree: %s   source variants: %s   -> new variant: %s"
          % (TREE_NAME, SOURCE_VARIANTS, NEW_VARIANT_NAME))
    print("New folder:", out_dir)
    print("-" * 90)

    # ---- taper.txt ----------------------------------------------------
    ref_h, median_diam, taper_per_variant, taper_path = build_median_taper(TREE_NAME, SOURCE_VARIANTS, out_dir)
    print("Saved:", taper_path)
    i_1_3 = int(np.argmin(np.abs(ref_h - 1.3)))
    print("  diameter @ height=%.4g m (source values vs. median):" % ref_h[i_1_3])
    for v in SOURCE_VARIANTS:
        print("    source %s: %.6f m" % (v, taper_per_variant[v][1][i_1_3]))
    print("    MEDIAN  : %.6f m" % median_diam[i_1_3])
    print()

    # ---- BranchStructure.txt -------------------------------------------
    _ref_bs_lines, branch_per_variant, median_branch_diam, branch_path = build_median_branch_structure(
        TREE_NAME, SOURCE_VARIANTS, out_dir)
    print("Saved:", branch_path)
    n_rows = len(median_branch_diam)
    sample_idx = sorted(set([0, n_rows // 4, n_rows // 2, n_rows - 1]))
    print("  sample rows (trunk + 3 spread across the row range), source diameters(m) -> median:")
    for idx in sample_idx:
        old_vals = [branch_per_variant[v][idx] for v in SOURCE_VARIANTS]
        print("    row %d: %s -> median=%.6f m"
              % (idx, ["%.6f" % x for x in old_vals], median_branch_diam[idx]))
    print()

    # ---- TreesParams.txt ------------------------------------------------
    per_variant_fields, median_fields, params_path = build_median_params(TREE_NAME, SOURCE_VARIANTS, out_dir)
    print("Saved:", params_path)
    header = "  %-22s" % "field" + "".join("%14s" % v for v in SOURCE_VARIANTS) + "%14s" % "MEDIAN"
    print(header)
    for key in median_fields:
        row = "".join("%14.6g" % per_variant_fields[v][key] for v in SOURCE_VARIANTS)
        print("  %-22s%s%14.6g" % (key, row, median_fields[key]))

    tv = median_fields.get("TreeVolume")
    trv = median_fields.get("TrunkVolume")
    brv = median_fields.get("BranchVolume")
    print()
    if tv is not None and trv is not None and brv is not None:
        exact = abs(tv - (trv + brv)) < 1e-9
        print("  TreeVolume == TrunkVolume + BranchVolume for the new median file: %s "
              "(TreeVolume=%.6g, TrunkVolume+BranchVolume=%.6g, diff=%.6g)" % (exact, tv, trv + brv, tv - (trv + brv)))
        print("  EXPECTED to be False: each field's median can come from a different source "
              "variant (see build_median_params()'s own docstring above) - harmless, since "
              "TreesParams.txt only feeds a comparison row, not the actual calibration.")
    else:
        print("  WARNING: TreeVolume/TrunkVolume/BranchVolume not all found in the median fields - "
              "cannot check TreeVolume == TrunkVolume + BranchVolume.")
    print("=" * 90)


# =========================  RUN  =====================================
if __name__ == "__main__":
    run()
