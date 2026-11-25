"""
check_png_provenance.py
-------------------------
Read-only audit for Stage 9 Part A3.

Every current ansys_eval script embeds its run's configuration as real PNG
metadata now (see savefig_with_provenance() in each script - Stage 9 Part
A1), not only as pixels in the visible footer. The Stage 8 D1 audit had to
guess a PNG's GRAVITY_VARIANT/MODEL_INCLUDE from its FILENAME, because the
footer was pixels only (no OCR available) - that guess could never catch a
file whose name disagrees with its actual configuration, which is exactly
the failure this script exists to catch.

This script reads the real metadata back (via PIL, no OCR needed) for
every PNG in out/ and lists:
  1. every PNG with NO Stage 9 metadata at all - it predates Stage 9, and
     its configuration cannot be verified this way;
  2. every PNG whose STORED GravityVariant or ModelInclude differs from
     the CURRENT default of the script named in its own "Script" field -
     the current default is read straight from that script's own source
     (not hardcoded here a second time), so this stays correct if a
     script's default ever changes.

This is a listing tool only - it changes nothing.

How to use it:
  python check_png_provenance.py
"""

import re
from pathlib import Path

from PIL import Image

WORK_FOLDER = Path(__file__).resolve().parent
OUT_FOLDER = WORK_FOLDER / "out"

# Every script that writes PNGs into out/ with Stage 9 metadata - matched
# against a PNG's own "Script" field to find which file's CURRENT default
# it should be compared against.
SCRIPT_FILES = [
    "plot_stiffness_comparison.py",
    "plot_modal_harmonic.py",
    "method_deviation_overview.py",
    "plot_frequency_spectra.py",
]


def read_current_default(script_name: str, constant_name: str):
    """
    Read script_name's own source for a "CONSTANT_NAME = <value>" line (the
    CONFIGURATION-section convention every ansys_eval script follows) and
    return the value as plain text (quotes stripped), or None if that
    script has no such constant at all - e.g. plot_frequency_spectra.py has
    no MODEL_INCLUDE knob, so callers must treat None as "do not compare",
    not as "the default is None".
    """
    source = (WORK_FOLDER / script_name).read_text()
    match = re.search(rf"^{constant_name} = (.+)$", source, re.MULTILINE)
    if match is None:
        return None
    return match.group(1).strip().strip('"').strip("'")


def main():
    print("=" * 78)
    print("check_png_provenance.py")
    print("=" * 78)

    current_gravity_variant = {name: read_current_default(name, "GRAVITY_VARIANT") for name in SCRIPT_FILES}
    current_model_include = {name: read_current_default(name, "MODEL_INCLUDE") for name in SCRIPT_FILES}
    print("\nCurrent defaults read from each script's own source:")
    for name in SCRIPT_FILES:
        print(f"  {name}: GRAVITY_VARIANT={current_gravity_variant[name]!r}  "
              f"MODEL_INCLUDE={current_model_include[name]!r}")

    png_paths = sorted(OUT_FOLDER.glob("*.png"))
    no_metadata = []
    non_default = []

    for path in png_paths:
        with Image.open(path) as img:
            info = dict(img.info)

        script = info.get("Script")
        if script is None:
            no_metadata.append(path.name)
            continue

        reasons = []
        expected_gravity = current_gravity_variant.get(script)
        stored_gravity = info.get("GravityVariant")
        if expected_gravity is not None and stored_gravity != expected_gravity:
            reasons.append(f"GravityVariant={stored_gravity!r} (current default for {script}: {expected_gravity!r})")

        expected_model_include = current_model_include.get(script)
        stored_model_include = info.get("ModelInclude")
        if expected_model_include is not None and stored_model_include != expected_model_include:
            reasons.append(f"ModelInclude={stored_model_include!r} "
                          f"(current default for {script}: {expected_model_include!r})")

        if reasons:
            non_default.append((path.name, script, reasons))

    print(f"\nTotal PNGs in out/: {len(png_paths)}")

    print(f"\n--- No Stage 9 metadata at all (predates Stage 9): {len(no_metadata)} ---")
    for name in no_metadata:
        print(f"  {name}")

    print(f"\n--- Non-default GravityVariant/ModelInclude per stored metadata: {len(non_default)} ---")
    for name, script, reasons in non_default:
        print(f"  {name} ({script}): {'; '.join(reasons)}")


if __name__ == "__main__":
    main()
