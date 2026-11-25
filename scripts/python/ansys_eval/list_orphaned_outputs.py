"""
list_orphaned_outputs.py
-------------------------
Standalone helper: lists filtered (MODEL_INCLUDE / TREES / SCENARIOS -
suffixed) output files in ansys_eval/out/ whose UNFILTERED sibling is newer
than they are - i.e. they were left behind by an older, differently-
configured run and may now be stale.

WHY this exists: on 11 Sep 2026, a MODEL_INCLUDE-filtered run of
plot_stiffness_comparison.py / plot_modal_harmonic.py /
method_deviation_overview.py wrote "_noacel"-suffixed PNGs into out/.
MODEL_INCLUDE was later reset to None and every later run only refreshed
the UNFILTERED files (writing a "_noacel" file is a no-op when
MODEL_INCLUDE is None - see with_filter_suffix()/add_model_include_suffix()
in the three scripts). The stale "_noacel" files sat there with nothing on
either file saying which was current, and were mistaken for current output.

This script does ONLY the listing - the same check each of the three
scripts already runs on itself at the end of its own main() (via
report_orphaned_filtered_outputs()) - so it can be run on its own, without
regenerating any figure. It NEVER deletes, moves, or renames anything:
printing the list is the entire feature. Deciding what to do with a listed
file is left to you.

How to use it:
  - Run the script from a terminal:  python list_orphaned_outputs.py
  - Read the printed list. Nothing is written or changed by running this.
"""

import re
from pathlib import Path
from datetime import datetime

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

OUT_FOLDER = Path(r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\scripts\python\ansys_eval\out")

# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

# Same markers/pattern as the copy of this logic inside each of the three
# plotting scripts (report_orphaned_filtered_outputs() and its neighbours) -
# copied here, not imported, per the ansys_eval self-contained convention.
KNOWN_MODEL_SUFFIX_MARKERS = ("_noacel",)
TREES_SCENARIOS_SUFFIX_PATTERN = re.compile(r"_trees-[0-9-]+|_sc-[0-9-]+")


def find_filtered_output_files(out_folder: Path) -> list:
    """
    Return every FILE (not directory) in out_folder that carries a filter-
    suffix marker anywhere in its name (KNOWN_MODEL_SUFFIX_MARKERS, or a
    "_trees-.../_sc-..." tree/scenario suffix). Unlike the per-script copy
    of this function, this standalone tool has no single "name_prefixes"
    family to restrict itself to - it scans the WHOLE out/ folder, since it
    is meant to be run independently of any one script.
    """
    hits = []
    for path in sorted(out_folder.iterdir()):
        if not path.is_file():
            continue
        if any(marker in path.name for marker in KNOWN_MODEL_SUFFIX_MARKERS) \
                or TREES_SCENARIOS_SUFFIX_PATTERN.search(path.name):
            hits.append(path)
    return hits


def unfiltered_sibling_path(filtered_path: Path) -> Path:
    """
    Strip every known filter-suffix marker from filtered_path's name to
    recover the UNFILTERED filename it was derived from, e.g.
    "stiffness_comparison_rotation_all_noacel.png" ->
    "stiffness_comparison_rotation_all.png". Used only to COMPARE
    timestamps below - never to read or write that file automatically.
    """
    stem = filtered_path.stem
    for marker in KNOWN_MODEL_SUFFIX_MARKERS:
        stem = stem.replace(marker, "")
    stem = TREES_SCENARIOS_SUFFIX_PATTERN.sub("", stem)
    return filtered_path.with_name(stem + filtered_path.suffix)


def main():
    if not OUT_FOLDER.is_dir():
        print(f"ERROR: OUT_FOLDER does not exist: {OUT_FOLDER}")
        return

    hits = find_filtered_output_files(OUT_FOLDER)
    if not hits:
        print(f"No filtered-output files (carrying '_noacel' or a '_trees-.../_sc-...' "
              f"suffix) found in {OUT_FOLDER}.")
        return

    print(f"Filtered-output files found in {OUT_FOLDER} ({len(hits)}):")
    print("NOTHING is touched by this script - it only lists and compares timestamps.\n")

    n_stale = 0
    for path in hits:
        sibling = unfiltered_sibling_path(path)
        mtime_text = datetime.fromtimestamp(path.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        if sibling.exists():
            sibling_mtime = sibling.stat().st_mtime
            if sibling_mtime > path.stat().st_mtime:
                n_stale += 1
                sibling_text = datetime.fromtimestamp(sibling_mtime).strftime("%Y-%m-%d %H:%M:%S")
                status = f"STALE - its unfiltered sibling ({sibling.name}) is newer ({sibling_text})"
            else:
                status = "unfiltered sibling is not newer than this file - may still be current"
        else:
            status = "no unfiltered sibling found to compare against"
        print(f"  {path.name}")
        print(f"    modified : {mtime_text}")
        print(f"    status   : {status}")

    print(f"\n{n_stale} of {len(hits)} filtered-output file(s) are STALE (unfiltered sibling is newer). "
          f"Review them by hand - this script changes nothing.")


if __name__ == "__main__":
    main()
