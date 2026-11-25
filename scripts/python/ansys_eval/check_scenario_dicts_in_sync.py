"""
check_scenario_dicts_in_sync.py
--------------------------------
Read-only check for Stage 8 Part C4. SCENARIO_MODELS (which ANSYS scenario
number maps to which readable model name) is copy-pasted - not imported -
into FOUR ansys_eval scripts, per this project's "ansys_eval stays
self-contained" convention (see each script's own copy for why: no
imports between these scripts). Copy-pasted code drifts silently unless
something checks it, so this script:

  1. Reads the four source files listed in SOURCE_FILES below.
  2. Extracts each file's SCENARIO_MODELS dict (via Python's own parser -
     ast.literal_eval - not a hand-rolled regex, so it can't be fooled by
     a stray "}" inside a comment or string).
  3. Compares all four dicts.
  4. Prints what it found, then exits 0 if they all agree, or a NON-ZERO
     exit code with a clear message naming exactly which file(s) disagree
     and how, if they don't.

This script only READS the four files - it never writes to them. If it
reports a mismatch, fix the disagreeing file(s) by hand and re-run it.

How to use it:
  python check_scenario_dicts_in_sync.py
"""

import ast
import re
import sys
from pathlib import Path

# The four files that each keep their own copy of SCENARIO_MODELS. Kept in
# one place here (rather than re-discovering them by searching the folder)
# so an new fifth copy in some future script is invisible to this checker
# until it is deliberately added below - a silent gap, not a crash, so if
# you add a SCENARIO_MODELS copy to a new script, add it to this list too.
SOURCE_FILES = [
    "load_static_scenarios.py",
    "plot_modal_harmonic.py",
    "method_deviation_overview.py",
    "plot_frequency_spectra.py",
]


def extract_scenario_models(file_path: Path) -> dict:
    """
    Find the "SCENARIO_MODELS = { ... }" assignment in file_path's source
    text and parse it into a real Python dict via ast.literal_eval (which
    tolerates the trailing "# comment" some copies have on the 999 line -
    Python's own tokenizer strips comments, a hand-rolled regex would not
    reliably). Raises SystemExit with a clear message if the file has no
    such assignment at all (a copy went missing, not just out of sync).
    """
    source = file_path.read_text()
    match = re.search(r"^SCENARIO_MODELS\s*=\s*\{", source, re.MULTILINE)
    if match is None:
        raise SystemExit(
            f"check_scenario_dicts_in_sync.py: {file_path.name} has no top-level "
            f"'SCENARIO_MODELS = {{' assignment at all - has it been renamed or removed?"
        )

    # Walk forward from the opening "{" counting brace depth, so the exact
    # end of the dict literal is found regardless of what is inside it.
    start = match.end() - 1  # index of the opening "{"
    depth = 0
    end = None
    for i in range(start, len(source)):
        if source[i] == "{":
            depth += 1
        elif source[i] == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end is None:
        raise SystemExit(
            f"check_scenario_dicts_in_sync.py: {file_path.name}'s SCENARIO_MODELS "
            f"dict literal never closes (unbalanced braces) - cannot parse it."
        )

    dict_text = source[start:end]
    return ast.literal_eval(dict_text)


def main():
    work_folder = Path(__file__).resolve().parent
    print("=" * 78)
    print("check_scenario_dicts_in_sync.py")
    print("=" * 78)

    dicts_by_file = {}
    for name in SOURCE_FILES:
        path = work_folder / name
        if not path.exists():
            raise SystemExit(f"check_scenario_dicts_in_sync.py: expected file not found: {path}")
        dicts_by_file[name] = extract_scenario_models(path)
        print(f"\n{name}:")
        for number in sorted(dicts_by_file[name]):
            print(f"  {number}: {dicts_by_file[name][number]!r}")

    reference_name = SOURCE_FILES[0]
    reference_dict = dicts_by_file[reference_name]

    mismatches = []
    for name in SOURCE_FILES[1:]:
        if dicts_by_file[name] != reference_dict:
            mismatches.append(name)

    print("\n" + "=" * 78)
    if not mismatches:
        print(f"IN SYNC - all {len(SOURCE_FILES)} SCENARIO_MODELS copies are identical.")
        print("=" * 78)
        return

    print(f"OUT OF SYNC - the following file(s) disagree with {reference_name}:")
    for name in mismatches:
        this_dict = dicts_by_file[name]
        all_numbers = sorted(set(reference_dict) | set(this_dict))
        for number in all_numbers:
            ref_value = reference_dict.get(number, "<missing>")
            this_value = this_dict.get(number, "<missing>")
            if ref_value != this_value:
                print(f"  {name}: scenario {number} = {this_value!r}, "
                      f"but {reference_name} has {ref_value!r}")
    print("=" * 78)
    sys.exit(1)


if __name__ == "__main__":
    main()
