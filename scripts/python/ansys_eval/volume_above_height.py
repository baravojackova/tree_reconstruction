"""
volume_above_height.py
-----------------------
Standalone report: how much of a tree's cylinder-model volume sits above a
given height.

WHY THIS MATTERS: mass high in the crown dominates the first bending mode,
mass near the base barely contributes. The split height is therefore a
PARAMETER, not a fixed number - the answer changes completely depending on
where the plane is put. On B21_S12 TreeQSM, splitting at the crown midpoint
gives 59.5% of total volume above; splitting at half the TREE height gives
98.5% of CROWN volume above. Both are correct, they answer different
questions. This script always reports which definition (SPLIT_MODE +
SPLIT_VALUE) produced a given number, and the actual height in meters that
definition resolved to.

GEOM FILE FORMAT (established against volume_results.csv - do not re-derive,
do not guess):
    line 1          "0 1 2 ... 10", a column-index row. Skipped.
    column 0        row index, 1..n
    columns 1,2,3   x, y, z of the cylinder START point
    column 4        LENGTH (m)
    column 5        RADIUS (m)
    columns 6,7,8   unit direction vector
    column 9        parent index
    column 10       branch order, 0 = trunk
Cylinder volume = pi * radius^2 * length. This layout is the same for the
AdTree exports and the TreeQSM exports.

Heights used throughout this script (tree_height_m, crown_bottom_m,
crown_top_m, split_height_m, and the split itself) are measured ABOVE THE
TREE'S OWN BASE (the lowest point of any cylinder in the file), not as the
raw z-coordinate - the geom files do not all start exactly at z = 0 (e.g.
the very first cylinder in the B21_S12 file starts at z = -0.0126 m), and
measuring from the file's own minimum keeps "1.3 m" meaning the same thing
in every file regardless of that offset.

ANTI-DELETION / READ-ONLY: the geom_*.txt files are inputs and are never
written to. This script produces exactly one output file:
    ansys_eval/out/volume_above_height.csv
Nothing else on disk is touched, moved, or regenerated.

How to use it:
  - Edit the CONFIGURATION block below.
  - Run it from a terminal:  python volume_above_height.py
  - Read the CHECKS that get printed first, then the results table.
"""

import csv
import math
import re
from pathlib import Path

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

REPO_ROOT = Path(__file__).resolve().parents[3]

# Folders to search for geom_*.txt files. Searched recursively, but any
# path with an "archive" folder in it is skipped - those are superseded
# copies kept for history, not current models.
GEOM_FOLDERS = [
    REPO_ROOT / "scripts" / "python" / "ansys_geom",
    REPO_ROOT / "scripts" / "matlab",
]

# None = every tree. A single tree ID ("B21_S12"), or a list of tree IDs,
# restricts the run to those trees only.
TREE_FILTER = None

# None = every model. A substring that must appear in the filename (e.g.
# "r5mm" or "man_pd") restricts the run to matching models only.
MODEL_FILTER = None

# How the split height is defined. One of:
#   "absolute"      - SPLIT_VALUE is a height in meters above the tree base
#   "tree_fraction"  - SPLIT_VALUE is a fraction (0-1) of total tree height
#   "crown_fraction" - SPLIT_VALUE is a fraction (0-1) of the crown extent,
#                      measured from the lowest branch-order>=CROWN_MIN_ORDER
#                      cylinder to the top of the tree
SPLIT_MODE = "tree_fraction"

# The number(s) that go with SPLIT_MODE. A single float, or a list of
# floats to compare several split heights in one run without re-running.
SPLIT_VALUE = [0.5]

# How a cylinder that straddles the split plane is counted:
#   "midpoint" - the cylinder belongs wholly to whichever side its
#                MIDPOINT height falls on (default)
#   "split"    - the cylinder's volume is divided in proportion to how much
#                of its length lies on each side of the plane
STRADDLE_HANDLING = "midpoint"

# Branch order at which the "crown" is considered to start. Order 0 is
# always the trunk.
CROWN_MIN_ORDER = 1

OUTPUT_CSV = REPO_ROOT / "scripts" / "python" / "ansys_eval" / "out" / "volume_above_height.csv"

# Smoke-test window (meters above tree base) used to sanity-check that
# columns 4 (length) and 5 (radius) have not been swapped: the diameter
# averaged over trunk cylinders in this band should look like a real DBH
# (tens of centimeters), not like a swapped-column artifact (~1.8 m).
DBH_CHECK_MIN_M = 1.0
DBH_CHECK_MAX_M = 1.6
DBH_SWAP_WARNING_THRESHOLD_M = 1.0


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

FILENAME_PATTERN = re.compile(r"^geom_([^_]+_[^_]+)_(.+)\.txt$", re.IGNORECASE)


def normalize_filter(value):
    """None stays None; a single value is wrapped in a list; a list is kept."""
    if value is None:
        return None
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def find_geom_files(folders):
    """
    Recursively find geom_*.txt under each folder, skipping anything with
    an "archive" path component. Returns a sorted list of Path objects.
    """
    found = []
    for folder in folders:
        folder = Path(folder)
        if not folder.is_dir():
            print(f"WARNING: geom folder does not exist, skipping: {folder}")
            continue
        for path in folder.rglob("geom_*.txt"):
            if any(part.lower() == "archive" for part in path.parts):
                continue
            found.append(path)
    return sorted(set(found))


def parse_tree_and_model(filename: str):
    """
    Split "geom_<TREE>_<model>.txt" into (tree, model). TREE is taken as
    the first two underscore-separated tokens after "geom_" (this project's
    tree IDs are of the form "B21_S12"); everything after that, with the
    .txt stripped, is the model name.
    Returns (None, None) if the filename does not match this pattern.
    """
    match = FILENAME_PATTERN.match(filename)
    if match is None:
        return None, None
    return match.group(1), match.group(2)


class GeomFileError(Exception):
    pass


def read_cylinders(path: Path):
    """
    Read one geom_*.txt file and return a list of cylinder dicts plus the
    checks needed to confirm the column layout held for this specific file.

    Each cylinder dict holds: row_index, x, y, z (start point, raw
    z-coordinate, not yet base-relative), length, radius, dz (unit
    direction z-component), parent, order, volume.

    Raises GeomFileError with a human-readable message if the file does
    not look like the expected format at all (wrong column count) - a
    missing/duplicate row index is a WARNING, not a hard failure, since
    the rest of the file may still be usable.
    """
    with open(path, "r", encoding="utf-8") as f:
        lines = [line for line in f.read().splitlines() if line.strip() != ""]

    if not lines:
        raise GeomFileError("file is empty")

    header_tokens = lines[0].split()
    warnings = []
    if header_tokens != [str(i) for i in range(len(header_tokens))]:
        warnings.append(
            f"first line is not a clean '0 1 2 ...' index row: {lines[0]!r}"
        )

    cylinders = []
    bad_column_count_rows = 0
    seen_indices = []
    for line in lines[1:]:
        columns = line.split()
        if len(columns) != 11:
            bad_column_count_rows += 1
            continue
        values = [float(c) for c in columns]
        row_index = int(values[0])
        seen_indices.append(row_index)
        x, y, z = values[1], values[2], values[3]
        length = values[4]
        radius = values[5]
        dz = values[7]
        parent = int(values[9])
        order = int(values[10])
        volume = math.pi * radius * radius * length
        cylinders.append({
            "row_index": row_index,
            "x": x, "y": y, "z": z,
            "length": length, "radius": radius, "dz": dz,
            "parent": parent, "order": order, "volume": volume,
        })

    if bad_column_count_rows:
        warnings.append(f"{bad_column_count_rows} data row(s) did not have 11 columns")

    expected_indices = list(range(1, len(cylinders) + bad_column_count_rows + 1))
    if seen_indices != expected_indices[:len(seen_indices)] and bad_column_count_rows == 0:
        if seen_indices != list(range(1, len(seen_indices) + 1)):
            warnings.append(
                f"column-0 row indices do not run 1..{len(seen_indices)} without gaps"
            )

    if not cylinders:
        raise GeomFileError("no usable data rows found")

    return cylinders, warnings


def cylinder_heights(cyl, base_z):
    """
    Return (bottom_height, top_height) for one cylinder, in meters above
    base_z. The end point is start + length * direction; since the
    direction vector is a unit vector, length * dz is the z-extent.
    """
    start_h = cyl["z"] - base_z
    end_h = start_h + cyl["length"] * cyl["dz"]
    return (start_h, end_h) if start_h <= end_h else (end_h, start_h)


def resolve_split_height(split_mode, split_value, tree_height_m, crown_bottom_m, crown_top_m):
    if split_mode == "absolute":
        return float(split_value)
    if split_mode == "tree_fraction":
        return float(split_value) * tree_height_m
    if split_mode == "crown_fraction":
        return crown_bottom_m + float(split_value) * (crown_top_m - crown_bottom_m)
    raise ValueError(f"unknown SPLIT_MODE: {split_mode!r}")


def split_volume_above_below(cylinders, base_z, split_height_m, straddle_handling):
    """
    Given cylinders (with raw z already read) and a split height in meters
    above base_z, return:
      volume_above, volume_below,
      crown_volume_above (restricted to order >= CROWN_MIN_ORDER),
      n_straddling, straddling_volume_fraction (of TOTAL volume)

    "midpoint": a cylinder counts wholly above/below based on the height
    of its own midpoint.
    "split": a cylinder crossing the plane has its volume divided in
    proportion to the length on each side (linear along the cylinder).
    """
    volume_above = 0.0
    volume_below = 0.0
    crown_volume_above = 0.0
    total_volume = 0.0
    straddling_volume = 0.0
    n_straddling = 0

    for cyl in cylinders:
        bottom_h, top_h = cylinder_heights(cyl, base_z)
        vol = cyl["volume"]
        total_volume += vol
        is_crown = cyl["order"] >= CROWN_MIN_ORDER

        crosses = bottom_h < split_height_m < top_h
        if crosses:
            n_straddling += 1
            straddling_volume += vol

        if straddle_handling == "split":
            extent = top_h - bottom_h
            if extent <= 0:
                above_fraction = 1.0 if bottom_h >= split_height_m else 0.0
            elif top_h <= split_height_m:
                above_fraction = 0.0
            elif bottom_h >= split_height_m:
                above_fraction = 1.0
            else:
                above_fraction = (top_h - split_height_m) / extent
            vol_above = vol * above_fraction
            vol_below = vol - vol_above
        else:  # "midpoint"
            midpoint_h = (bottom_h + top_h) / 2.0
            if midpoint_h > split_height_m:
                vol_above, vol_below = vol, 0.0
            else:
                vol_above, vol_below = 0.0, vol

        volume_above += vol_above
        volume_below += vol_below
        if is_crown:
            crown_volume_above += vol_above

    straddling_fraction = (straddling_volume / total_volume) if total_volume > 0 else 0.0
    return volume_above, volume_below, crown_volume_above, n_straddling, straddling_fraction


def analyze_file(path: Path, tree: str, model: str, split_values):
    cylinders, format_warnings = read_cylinders(path)
    for msg in format_warnings:
        print(f"  WARNING [{path.name}]: {msg}")

    n_cylinders = len(cylinders)

    all_bottoms = [cyl["z"] for cyl in cylinders]
    all_bottoms += [cyl["z"] + cyl["length"] * cyl["dz"] for cyl in cylinders]
    base_z = min(all_bottoms)

    heights = [cylinder_heights(cyl, base_z) for cyl in cylinders]
    tree_height_m = max(top for _, top in heights)

    crown_heights = [
        cylinder_heights(cyl, base_z)
        for cyl in cylinders
        if cyl["order"] >= CROWN_MIN_ORDER
    ]
    if crown_heights:
        crown_bottom_m = min(bottom for bottom, _ in crown_heights)
        crown_top_m = max(top for _, top in crown_heights)
    else:
        crown_bottom_m = float("nan")
        crown_top_m = float("nan")
        print(f"  WARNING [{path.name}]: no cylinders with order >= {CROWN_MIN_ORDER} - no crown found")

    volume_total_m3 = sum(cyl["volume"] for cyl in cylinders)
    volume_trunk_m3 = sum(cyl["volume"] for cyl in cylinders if cyl["order"] == 0)
    volume_crown_m3 = sum(cyl["volume"] for cyl in cylinders if cyl["order"] >= CROWN_MIN_ORDER)

    # --- smoke test: trunk diameter between DBH_CHECK_MIN_M and DBH_CHECK_MAX_M ---
    dbh_diameters = []
    for cyl, (bottom_h, top_h) in zip(cylinders, heights):
        if cyl["order"] == 0 and DBH_CHECK_MIN_M <= bottom_h <= DBH_CHECK_MAX_M:
            dbh_diameters.append(2.0 * cyl["radius"])
    dbh_estimate_m = sum(dbh_diameters) / len(dbh_diameters) if dbh_diameters else float("nan")

    if dbh_diameters and dbh_estimate_m > DBH_SWAP_WARNING_THRESHOLD_M:
        print(
            f"  *** LOUD WARNING [{path.name}]: trunk diameter at "
            f"{DBH_CHECK_MIN_M}-{DBH_CHECK_MAX_M} m = {dbh_estimate_m:.3f} m - "
            f"not physical, columns 4 and 5 (length/radius) may be swapped! ***"
        )
    elif not dbh_diameters:
        print(f"  WARNING [{path.name}]: no trunk cylinders found in the "
              f"{DBH_CHECK_MIN_M}-{DBH_CHECK_MAX_M} m band - cannot run DBH smoke test")
    else:
        print(f"  smoke test [{path.name}]: trunk diameter at "
              f"{DBH_CHECK_MIN_M}-{DBH_CHECK_MAX_M} m = {dbh_estimate_m:.3f} m (looks physical)")

    rows = []
    for split_value in split_values:
        split_height_m = resolve_split_height(
            SPLIT_MODE, split_value, tree_height_m, crown_bottom_m, crown_top_m
        )
        (volume_above, volume_below, crown_volume_above,
         n_straddling, straddling_fraction) = split_volume_above_below(
            cylinders, base_z, split_height_m, STRADDLE_HANDLING
        )

        if STRADDLE_HANDLING == "midpoint" and n_straddling > 0:
            print(
                f"  straddle check [{path.name}, split={split_height_m:.3f} m]: "
                f"{n_straddling} cylinder(s) cross the plane, holding "
                f"{straddling_fraction * 100:.2f}% of total volume "
                f"(assigned whole-cylinder by midpoint rule)"
            )

        pct_above_total = (volume_above / volume_total_m3 * 100.0) if volume_total_m3 > 0 else float("nan")
        pct_above_crown = (crown_volume_above / volume_crown_m3 * 100.0) if volume_crown_m3 > 0 else float("nan")

        rows.append({
            "tree": tree,
            "model": model,
            "source_folder": str(path.parent),
            "n_cylinders": n_cylinders,
            "split_mode": SPLIT_MODE,
            "split_value": split_value,
            "split_height_m": round(split_height_m, 4),
            "tree_height_m": round(tree_height_m, 4),
            "crown_bottom_m": round(crown_bottom_m, 4) if crown_heights else None,
            "crown_top_m": round(crown_top_m, 4) if crown_heights else None,
            "volume_total_m3": round(volume_total_m3, 4),
            "volume_trunk_m3": round(volume_trunk_m3, 4),
            "volume_crown_m3": round(volume_crown_m3, 4),
            "volume_above_m3": round(volume_above, 4),
            "volume_below_m3": round(volume_below, 4),
            "pct_above_total": round(pct_above_total, 2),
            "crown_volume_above_m3": round(crown_volume_above, 4),
            "pct_above_crown": round(pct_above_crown, 2),
        })

    return rows


def print_table(rows):
    if not rows:
        print("(no rows)")
        return
    columns = list(rows[0].keys())
    widths = {c: max(len(c), max(len(str(r[c])) for r in rows)) for c in columns}
    header = " | ".join(c.ljust(widths[c]) for c in columns)
    print(header)
    print("-" * len(header))
    for r in rows:
        print(" | ".join(str(r[c]).ljust(widths[c]) for c in columns))


def main():
    split_values = normalize_filter(SPLIT_VALUE)
    tree_filter = normalize_filter(TREE_FILTER)

    print("=" * 78)
    print("CONFIGURATION")
    print("=" * 78)
    print(f"GEOM_FOLDERS      : {[str(f) for f in GEOM_FOLDERS]}")
    print(f"TREE_FILTER       : {tree_filter}")
    print(f"MODEL_FILTER      : {MODEL_FILTER}")
    print(f"SPLIT_MODE        : {SPLIT_MODE}")
    print(f"SPLIT_VALUE       : {split_values}")
    print(f"STRADDLE_HANDLING : {STRADDLE_HANDLING}")
    print(f"CROWN_MIN_ORDER   : {CROWN_MIN_ORDER}")
    print()

    all_files = find_geom_files(GEOM_FOLDERS)
    print(f"Found {len(all_files)} geom_*.txt file(s) before filename/filter checks.")

    selected = []
    for path in all_files:
        tree, model = parse_tree_and_model(path.name)
        if tree is None:
            print(f"  WARNING: filename does not match 'geom_<TREE>_<model>.txt', skipping: {path}")
            continue
        if tree_filter is not None and tree not in tree_filter:
            continue
        if MODEL_FILTER is not None and MODEL_FILTER not in model:
            continue
        selected.append((path, tree, model))

    print(f"{len(selected)} file(s) selected after TREE_FILTER / MODEL_FILTER.")
    print()

    print("=" * 78)
    print("CHECKS (per file)")
    print("=" * 78)

    all_rows = []
    for path, tree, model in selected:
        print(f"[{tree} / {model}]  {path}")
        try:
            rows = analyze_file(path, tree, model, split_values)
        except GeomFileError as exc:
            print(f"  ERROR: {exc} - skipping this file")
            continue
        all_rows.extend(rows)
        print()

    print("=" * 78)
    print("RESULTS")
    print("=" * 78)
    print_table(all_rows)

    if not all_rows:
        print("\nNo rows produced - nothing written.")
        return

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    if OUTPUT_CSV.exists():
        print(f"\nSTOPPING: output file already exists, not overwriting: {OUTPUT_CSV}")
        return

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nSaved {len(all_rows)} row(s) to: {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
