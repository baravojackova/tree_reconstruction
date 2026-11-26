"""
plot_frequency_spectra.py
---------------------------
Frequency-spectrum figures for the two reconstruction models actually
exported to ANSYS: adtree_r5mm and treeqsm_pd002 (there is no "AdQSM" model
in the ANSYS results - AdTree's geometry traces back to an AdQSM-derived
taper curve, which is a modelling detail upstream of ANSYS, not a third
model here - see main()'s own model-count check, which stops the script if
it ever finds anything other than exactly these two "_no_acel" models).

Reads the harmonic sweep files (out_harmonic_bk*_scenario*.txt, ~100
frequency points per tree/model) that plot_modal_harmonic.py also reads,
but this script's whole job is the shape of that sweep curve - not a single
peak number - which did not fit into a script already at 2980 lines.

Produces:
  1. frequency_spectra_combined.png   - ALL trees, both models, ONE axes
  2. frequency_spectra_grid.png       - one small panel per tree, both
                                        models, shared axes (see main()'s
                                        own report for which of these two is
                                        actually readable)
  3. frequency_spectra_auc.png + .csv - area under each tree/model's curve
                                        (trapezoid rule) and the relative
                                        deviation between the two models
  4. A validity check (printed only, no figure) that MUST be read before
     trusting figure 3 at all - see check_spectrum_validity()'s own
     docstring for why an integral over a truncated, undecayed peak is
     meaningless.

Same GRAVITY_VARIANT structural lock as the other three ansys_eval scripts
(default "no_acel" - the reference variant), same provenance footer, same
freshness guard, same self-contained convention (every helper below is
copy-pasted from plot_modal_harmonic.py, not imported).

How to use it:
  - Edit the values in the "CONFIGURATION" section below.
  - Run the script from a terminal:  python plot_frequency_spectra.py
  - Read the VALIDITY CHECK printed first - it tells you whether the AUC
    figure/CSV can be trusted for each tree, before you look at them.
"""

import re
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import sys

# This script lives in scripts/python/ansys_eval, but paths.py lives one
# folder up (scripts/python). Python only looks for imports next to the
# script, so we first add the folder above to its search list.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import paths

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

INPUT_FOLDER = str(paths.ANSYS_DIR)
HARMONIC_FILE_PATTERN = "out_harmonic_bk*_scenario*.txt"

GEOMETRY_FILE_NAME = "Popis_Babice_VSE_09012025_ans_lev6_0_azimcor.txt"
GEOMETRY_FILE_PATH = Path(INPUT_FOLDER) / GEOMETRY_FILE_NAME
TREE_ID_PREFIX = "B21_S"
GEOMETRY_COLUMN_ANSYS_ORDER = 0
GEOMETRY_COLUMN_TREE_NUMBER = 2

# Identical to the dict in load_static_scenarios.py/plot_modal_harmonic.py/
# method_deviation_overview.py - copied, not imported. Kept in sync by
# check_scenario_dicts_in_sync.py (Stage 8 Part C4) - run it after editing
# this dict.
SCENARIO_MODELS = {
    1: "adtree_r5mm",
    2: "treeqsm_pd002",
    3: "adtree_r5mm_no_acel",
    4: "treeqsm_pd002_no_acel",
    999: "adtree_r5_pilot",
}

OUTPUT_DIR = paths.ANSYS_EVAL_DIR / "out"
SPECTRA_COMBINED_PNG_PATH = OUTPUT_DIR / "frequency_spectra_combined.png"
SPECTRA_GRID_PNG_PATH = OUTPUT_DIR / "frequency_spectra_grid.png"
SPECTRA_AUC_PNG_PATH = OUTPUT_DIR / "frequency_spectra_auc.png"
SPECTRA_AUC_CSV_PATH = OUTPUT_DIR / "frequency_spectra_auc.csv"

# The harmonic sweep column used throughout this pipeline as "the"
# frequency response - the same one compute_harmonic_peak_frequency() in
# plot_modal_harmonic.py uses for the harmonic peak (max |Top_ux|). Top_ux
# is the tree-top x-displacement response; the "_hz" column is the swept
# frequency itself, not a quantity choice.
RESPONSE_COLUMN = "Top_ux"

# --- THE GRAVITY-VARIANT STRUCTURAL LOCK (same as the other 3 scripts) ---
# "no_acel" = the gravity-off runs (scenarios 3/4) - the reference variant,
#             and the default. "with_gravity" = scenarios 1/2. "both" = no
#             filtering. Applied to DRAWING only, never to the loaded
#             harmonic_df itself, so this script's own CSV (if it ever
#             wrote one from harmonic_df directly - it currently does not)
#             would stay complete regardless - see the other 3 scripts'
#             own GRAVITY_VARIANT comment for the postmortem of getting
#             this backwards.
GRAVITY_VARIANT = "no_acel"

# True = stamp a small grey provenance footer on every figure.
FIGURE_PROVENANCE_FOOTER = True

# True = refuse to run if any raw harmonic file is newer than... nothing,
# actually - this script reads out_harmonic_bk* directly on every run, so
# there is no intermediate CSV to go stale (same situation as
# plot_modal_harmonic.py's own ABORT_ON_STALE_INPUT comment). Kept for
# interface consistency; check_input_freshness() reports "nothing to
# check" for this script and always returns the newest raw mtime, which
# still feeds the provenance footer.
ABORT_ON_STALE_INPUT = True

# --- SPECTRUM FIGURE OPTIONS (C1) -----------------------------------------
# True = divide each tree/model's curve by its OWN peak amplitude before
# plotting, so curve SHAPES can be compared when magnitudes differ by
# orders of magnitude across trees. False = plot raw magnitudes. Only
# affects SPECTRUM_COMBINED_STYLE="all_curves" and the per-tree grid figure
# - "summary" ALWAYS normalizes each tree first (see its own comment), so
# this flag has no effect on the combined figure while "summary" is active.
SPECTRUM_NORMALIZE = False

# True = log-scale the y-axis - amplitudes span several orders of
# magnitude across trees (see the validity check below), so try both and
# see which is legible; False = linear y-axis.
SPECTRUM_LOG_Y = False

# --- COMBINED SPECTRUM FIGURE STYLE (Stage 9 Part C) ----------------------
# "summary" (default) = one MEDIAN curve per model, each tree normalized to
#                        its own peak first (so one large tree cannot set
#                        the median), with a shaded 25th-75th percentile
#                        band showing the spread across trees. Answers "do
#                        the two reconstructions produce the same overall
#                        response shape?" - the question the combined
#                        figure was for.
# "all_curves"          = the original layout: every tree's raw (or
#                          SPECTRUM_NORMALIZE'd) curve drawn on the same
#                          axes - kept available, not removed, since
#                          Stage 8's own verdict was that neither
#                          SPECTRUM_NORMALIZE nor SPECTRUM_LOG_Y actually
#                          solves 22 overlapping lines on one axes; the
#                          per-tree grid figure (build_grid_spectrum_figure)
#                          remains the figure for seeing every curve.
SPECTRUM_COMBINED_STYLE = "summary"

# --- AREA-UNDER-CURVE OPTIONS (C2) ----------------------------------------
# Which model's AUC is the NUMERATOR's "other side" - i.e. the deviation
# formula is:
#   100 * (AUC[AUC_DEVIATION_MODEL_B] - AUC[AUC_DEVIATION_MODEL_A]) / AUC[AUC_DEVIATION_MODEL_A]
# Flip these two to flip the sign convention WITHOUT touching the body of
# compute_auc_deviation() - see that function's own docstring. Positive
# means MODEL_B's total response magnitude is larger.
AUC_DEVIATION_MODEL_A = "adtree_r5mm_no_acel"
AUC_DEVIATION_MODEL_B = "treeqsm_pd002_no_acel"

# --- TEXT SIZE -------------------------------------------------------------
# Minimum readable size for every title/axis-label/tick-label/legend across
# all three figures in this script.
TITLE_FONTSIZE = 15
AXIS_LABEL_FONTSIZE = 15
TICK_LABEL_FONTSIZE = 15
LEGEND_FONTSIZE = 15

# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

# Every name assigned in the CONFIGURATION section above (between the
# "CONFIGURATION" and "CODE" headers), used by print_active_configuration()
# below (Stage 8 Part B3).
CONFIG_CONSTANT_NAMES = [
    "INPUT_FOLDER", "HARMONIC_FILE_PATTERN", "GEOMETRY_FILE_NAME", "GEOMETRY_FILE_PATH",
    "TREE_ID_PREFIX", "GEOMETRY_COLUMN_ANSYS_ORDER", "GEOMETRY_COLUMN_TREE_NUMBER",
    "SCENARIO_MODELS", "OUTPUT_DIR", "SPECTRA_COMBINED_PNG_PATH", "SPECTRA_GRID_PNG_PATH",
    "SPECTRA_AUC_PNG_PATH", "SPECTRA_AUC_CSV_PATH", "RESPONSE_COLUMN", "GRAVITY_VARIANT",
    "FIGURE_PROVENANCE_FOOTER", "ABORT_ON_STALE_INPUT", "SPECTRUM_NORMALIZE", "SPECTRUM_LOG_Y",
    "SPECTRUM_COMBINED_STYLE", "AUC_DEVIATION_MODEL_A", "AUC_DEVIATION_MODEL_B",
    "TITLE_FONTSIZE", "AXIS_LABEL_FONTSIZE", "TICK_LABEL_FONTSIZE", "LEGEND_FONTSIZE",
]


def print_active_configuration(script_label: str, config_names: list) -> None:
    """
    Print the current value of every name in config_names, one per line -
    the exact, complete state of this script's CONFIGURATION section for
    this run (Stage 8 Part B3). Printed first, before any other output, so
    the console log alone says what produced a given figure - no need to
    re-open the script and check by eye.
    """
    print("\n" + "=" * 78)
    print(f"ACTIVE CONFIGURATION - {script_label}")
    print("=" * 78)
    for name in config_names:
        print(f"  {name} = {globals()[name]!r}")
    print("=" * 78)


PASTEL_SCENARIO_COLORS = [
    "#7E9FC7", "#D09882", "#78B7A0", "#D1B16C",
    "#CFA1B2", "#6CA46C", "#8B85B3", "#CD8B8B",
]

NO_ACEL_PATTERN = re.compile(r"[ _]no[ _]acel$", re.IGNORECASE)


def is_no_acel_model(model_name: str) -> bool:
    return bool(NO_ACEL_PATTERN.search(str(model_name)))


def base_model_name(model_name: str) -> str:
    return NO_ACEL_PATTERN.sub("", str(model_name))


def parse_czech_number(text: str) -> float:
    return float(text.replace(",", "."))


def tree_number_to_tree_id(tree_number: int) -> str:
    return f"{TREE_ID_PREFIX}{tree_number:02d}"


def load_ansys_order_to_tree_id(geometry_file_path) -> dict:
    """Position -> canonical tree ID. Copied from load_static_scenarios.py."""
    path = Path(geometry_file_path)
    if not path.is_file():
        print(f"ERROR: description file not found: {path}")
        raise FileNotFoundError(path)
    with open(geometry_file_path, "r", encoding="utf-8") as f:
        lines = [line for line in f.read().splitlines() if line.strip() != ""]
    data_lines = lines[2:]
    order_to_tree_id = {}
    for line in data_lines:
        columns = line.split("\t")
        ansys_order = int(parse_czech_number(columns[GEOMETRY_COLUMN_ANSYS_ORDER]))
        tree_number = int(parse_czech_number(columns[GEOMETRY_COLUMN_TREE_NUMBER]))
        order_to_tree_id[ansys_order] = tree_number_to_tree_id(tree_number)
    return order_to_tree_id


def scenario_to_model_name(scenario_number: int) -> str:
    return SCENARIO_MODELS.get(scenario_number, f"scenario{scenario_number}")


HARMONIC_FILENAME_PATTERN = re.compile(r"out_harmonic_bk(\d+)_scenario(\d+)\.txt")
HARMONIC_COLUMN_NAMES = [
    "frequency_hz", "Top_ux", "Top_uy", "Top_uz",
    "Half_ux", "Half_uy", "Half_uz", "Quart_ux", "Quart_uy", "Quart_uz",
]


def load_harmonic_data(input_folder, file_pattern: str, order_to_tree_id: dict) -> pd.DataFrame:
    """
    Load every harmonic sweep file into one combined DataFrame with
    columns tree/ansys_no/scenario/model/source_file plus
    HARMONIC_COLUMN_NAMES. Simplified copy of plot_modal_harmonic.py's own
    load_harmonic_data() (that version also reads the mass-line search
    used by the modal loader - not needed here, harmonic files have no
    mass line).
    """
    matching_files = sorted(Path(input_folder).glob(file_pattern))
    print(f"Found {len(matching_files)} harmonic file(s) matching '{file_pattern}'")
    rows = []
    for file_path in matching_files:
        match = HARMONIC_FILENAME_PATTERN.match(file_path.name)
        if match is None:
            print(f"WARNING: '{file_path.name}' does not match the expected harmonic "
                  f"filename pattern - skipping.")
            continue
        ansys_order, scenario = int(match.group(1)), int(match.group(2))
        tree = order_to_tree_id.get(ansys_order)
        if tree is None:
            print(f"WARNING: '{file_path.name}' refers to position bk{ansys_order}, which "
                  f"does not exist in the description file - skipping.")
            continue
        model = scenario_to_model_name(scenario)

        with open(file_path, "r", encoding="utf-8") as f:
            lines = [line for line in f.read().splitlines() if line.strip() != ""]
        # Line 1 is the header, data starts on line 2 - matches
        # plot_modal_harmonic.py's own load_harmonic_data() exactly,
        # including the tolerance below (copied from there, not imported).
        for line in lines[1:]:
            values = [float(v) for v in line.split()]
            # A row may have FEWER than 10 values if only some displacement
            # columns were populated for this run (seen in practice, not a
            # hypothetical - RESPONSE_COLUMN="Top_ux" is always present,
            # which is all this script actually uses). Pad the rest with
            # NaN rather than discarding the row - discarding was tried
            # first and silently produced an EMPTY dataframe for files that
            # only ever populate 3 of the 10 columns.
            row_values = values + [float("nan")] * (len(HARMONIC_COLUMN_NAMES) - len(values))
            row = dict(zip(HARMONIC_COLUMN_NAMES, row_values))
            row.update({"tree": tree, "ansys_no": ansys_order, "scenario": scenario,
                       "model": model, "source_file": file_path.name})
            rows.append(row)
    df = pd.DataFrame(rows)
    if len(df) > 0:
        df = df.sort_values(["tree", "model", "frequency_hz"]).reset_index(drop=True)
    return df


def filter_by_gravity_variant(df: pd.DataFrame, variant: str, model_col: str = "model") -> pd.DataFrame:
    if variant == "both":
        return df
    if variant == "no_acel":
        return df[df[model_col].apply(is_no_acel_model)]
    if variant == "with_gravity":
        return df[~df[model_col].apply(is_no_acel_model)]
    raise ValueError(f"Unknown GRAVITY_VARIANT: {variant!r} - use 'no_acel', 'with_gravity', or 'both'.")


def assert_gravity_variant(df: pd.DataFrame, expected_variant: str, function_name: str,
                           model_col: str = "model") -> None:
    if expected_variant == "both" or len(df) == 0:
        return
    models = df[model_col].unique()
    if expected_variant == "no_acel":
        offenders = [m for m in models if not is_no_acel_model(m)]
    elif expected_variant == "with_gravity":
        offenders = [m for m in models if is_no_acel_model(m)]
    else:
        raise ValueError(f"Unknown gravity variant to assert: {expected_variant!r}")
    assert not offenders, (
        f"{function_name}: expected ONLY '{expected_variant}' models but found "
        f"{offenders} - the gravity-variant structural lock was bypassed upstream."
    )


def newest_mtime(folder, pattern: str):
    matches = sorted(Path(folder).glob(pattern))
    if not matches:
        return None, None
    newest_path = max(matches, key=lambda p: p.stat().st_mtime)
    return newest_path.stat().st_mtime, newest_path


def check_input_freshness(script_label: str, raw_folder, raw_pattern: str,
                          csv_paths_being_read: dict, abort_if_stale: bool) -> float:
    newest_raw_mtime, newest_raw_path = newest_mtime(raw_folder, raw_pattern)
    if newest_raw_mtime is None:
        print(f"FRESHNESS CHECK ({script_label}): no files match '{raw_pattern}' in {raw_folder}.")
        return None
    if not csv_paths_being_read:
        print(f"FRESHNESS CHECK ({script_label}): reads {raw_pattern} directly every run - "
              f"nothing to check.")
        return newest_raw_mtime
    stale_found = False
    for label, csv_path in csv_paths_being_read.items():
        csv_path = Path(csv_path)
        if not csv_path.is_file():
            continue
        if csv_path.stat().st_mtime < newest_raw_mtime:
            stale_found = True
            print(f"*** WARNING: STALE INPUT - {label} is older than {newest_raw_path.name}. ***")
    if stale_found and abort_if_stale:
        raise SystemExit(f"{script_label}: aborting because ABORT_ON_STALE_INPUT=True.")
    return newest_raw_mtime


def add_provenance_footer(fig, script_filename: str, newest_input_mtime, active_filters_text: str,
                          n_trees: int, n_models: int) -> None:
    """
    Stamp a small grey monospace footer along the bottom edge of fig. A
    no-op when FIGURE_PROVENANCE_FOOTER is False.

    Stage 9 Part A1: also stashes the SAME provenance as fig.stage9_png_metadata,
    picked up by savefig_with_provenance() below when the figure is actually
    written to disk. The footer above is pixels only - reading it back needs
    OCR. This is a real PNG tEXt chunk, readable afterwards via PIL's
    Image.info without OCR or guessing a run's configuration from its
    filename - see check_png_provenance.py.
    """
    if not FIGURE_PROVENANCE_FOOTER:
        fig.stage9_png_metadata = {}
        return
    generated_text = datetime.now().strftime("%Y-%m-%d %H:%M")
    newest_input_text = (datetime.fromtimestamp(newest_input_mtime).strftime("%Y-%m-%d %H:%M")
                         if newest_input_mtime is not None else "unknown")
    line1 = f"{script_filename} | generated {generated_text} | newest input {newest_input_text}"
    line2 = f"{active_filters_text}  |  {n_trees} trees, {n_models} models"
    fig.text(0.01, -0.045, line1 + "\n" + line2, ha="left", va="bottom",
             fontsize=6, family="monospace", color="0.55")

    fig.stage9_png_metadata = {
        "Script": script_filename,
        "Generated": generated_text,
        "NewestInput": newest_input_text,
        "GravityVariant": str(GRAVITY_VARIANT),
        "ModelInclude": "N/A (this script has no MODEL_INCLUDE knob)",
        "Trees": "N/A (this script has no TREES knob - always draws every tree found)",
        "NTrees": str(n_trees),
        "NModels": str(n_models),
    }


def savefig_with_provenance(fig, path) -> None:
    """
    The ONLY way a figure should be written to disk in this script (Stage 9
    Part A1) - every build_*_figure() function calls add_provenance_footer(
    fig, ...) before saving, which stashes fig.stage9_png_metadata; this
    passes that dict straight to savefig(metadata=...), embedding it as
    real PNG metadata alongside the existing visible footer. A PNG with no
    such metadata predates Stage 9 and its configuration cannot be verified
    this way - see check_png_provenance.py.
    """
    plt.savefig(path, dpi=150, bbox_inches="tight", metadata=getattr(fig, "stage9_png_metadata", {}))


def assign_model_colors(models) -> dict:
    """
    Same alphabetical-assignment convention as every other ansys_eval
    script's assign_model_colors()/assign_label_colors() - assigned from
    the FULL model set (all 4), not just the two drawn under
    GRAVITY_VARIANT="no_acel", so a model keeps the SAME colour it has in
    every other figure in this pipeline (pink for adtree_r5mm_no_acel,
    yellow for treeqsm_pd002_no_acel - the same two colours Stage 5's
    method_deviation_overview.py figures use).
    """
    n_colors = len(PASTEL_SCENARIO_COLORS)
    return {m: PASTEL_SCENARIO_COLORS[i % n_colors] for i, m in enumerate(sorted(models))}


# --- C3: validity check ----------------------------------------------------

def check_spectrum_validity(harmonic_df: pd.DataFrame, tree_order: list, models: list,
                            decay_ratio_threshold: float = 0.1,
                            edge_fraction_threshold: float = 0.95) -> pd.DataFrame:
    """
    THE CHECK THAT DECIDES WHETHER THE AUC FIGURE IS VALID AT ALL. An
    integral over a curve whose peak sits at (or past) the edge of the
    swept range is meaningless - it would report "total response" for a
    response that was actually cut off. This is the SAME failure mode
    already found once in this project (the 20-mode "max EFF MASS" pick
    landing on mode 19/20 - i.e. at the edge of the extracted window, in a
    different part of this pipeline) - this function exists so it is not
    repeated silently here.

    For every tree/model, prints: the peak frequency, the amplitude AT the
    peak, the amplitude at the HIGHEST swept frequency, and their ratio
    (amplitude_at_edge / amplitude_at_peak - small is good, i.e. the curve
    decayed by the edge of the sweep).

    A curve is flagged INVALID for the AUC comparison if EITHER:
      - its peak sits at or past edge_fraction_threshold (95%) of the way
        through the swept frequency range (peak too close to the edge to
        trust that the true peak was captured at all), OR
      - the edge/peak amplitude ratio exceeds decay_ratio_threshold (10%)
        (the response has not meaningfully decayed by the top of the
        sweep, so the integral is missing an unknown amount of area past
        the edge).

    Returns a DataFrame (tree, model, peak_freq_hz, amplitude_at_peak,
    amplitude_at_edge, edge_to_peak_ratio, peak_position_fraction, valid)
    for the caller to gate compute_auc() on.
    """
    print("\n" + "=" * 78)
    print("C3 - SPECTRUM VALIDITY CHECK (read this before trusting the AUC figure/CSV)")
    print("=" * 78)
    rows = []
    for tree in tree_order:
        for model in models:
            curve = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)]
            curve = curve.sort_values("frequency_hz")
            if len(curve) == 0:
                continue
            amplitude = curve[RESPONSE_COLUMN].abs()
            freq = curve["frequency_hz"]
            peak_idx = amplitude.idxmax()
            peak_freq = float(freq.loc[peak_idx])
            amplitude_at_peak = float(amplitude.loc[peak_idx])
            amplitude_at_edge = float(amplitude.iloc[-1])
            edge_to_peak_ratio = (amplitude_at_edge / amplitude_at_peak) if amplitude_at_peak > 0 else float("nan")
            freq_lo, freq_hi = float(freq.min()), float(freq.max())
            peak_position_fraction = (peak_freq - freq_lo) / (freq_hi - freq_lo) if freq_hi > freq_lo else 1.0

            peak_near_edge = peak_position_fraction >= edge_fraction_threshold
            not_decayed = not pd.isna(edge_to_peak_ratio) and edge_to_peak_ratio > decay_ratio_threshold
            valid = not (peak_near_edge or not_decayed)

            rows.append({
                "tree": tree, "model": model, "peak_freq_hz": round(peak_freq, 4),
                "amplitude_at_peak": amplitude_at_peak, "amplitude_at_edge": amplitude_at_edge,
                "edge_to_peak_ratio": round(edge_to_peak_ratio, 4) if not pd.isna(edge_to_peak_ratio) else None,
                "peak_position_fraction": round(peak_position_fraction, 3),
                "valid": valid,
            })
            flag = "" if valid else "  <-- INVALID for AUC: "
            if not valid:
                reasons = []
                if peak_near_edge:
                    reasons.append(f"peak at {peak_position_fraction:.0%} of swept range (>= {edge_fraction_threshold:.0%})")
                if not_decayed:
                    reasons.append(f"edge/peak amplitude ratio {edge_to_peak_ratio:.3f} (> {decay_ratio_threshold})")
                flag += "; ".join(reasons)
            print(f"  {tree} / {model}: peak={peak_freq:.4f} Hz (amp={amplitude_at_peak:.4g}), "
                  f"edge amp={amplitude_at_edge:.4g}, edge/peak ratio={edge_to_peak_ratio:.4g}{flag}")

    validity_df = pd.DataFrame(rows)
    n_invalid = int((~validity_df["valid"]).sum()) if len(validity_df) > 0 else 0
    print(f"\n{n_invalid} of {len(validity_df)} curve(s) flagged INVALID for the AUC comparison.")
    return validity_df


# --- C2: area under the curve and model-vs-model deviation -----------------

def compute_auc(harmonic_df: pd.DataFrame, tree_order: list, models: list) -> pd.DataFrame:
    """
    Trapezoid-rule integral (numpy.trapezoid) of |RESPONSE_COLUMN| over
    frequency, for every tree/model. VERIFIES the two models' frequency
    grids are identical before integrating (they should be - both are the
    same ANSYS harmonic sweep setup) - if they are not, resamples the
    narrower one onto the wider one's common overlapping range via linear
    interpolation and prints that this happened, rather than integrating
    two different ranges silently.
    """
    print("\nFrequency-grid check (must be identical, or this function resamples and says so):")
    rows = []
    for tree in tree_order:
        grids = {}
        for model in models:
            curve = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)].sort_values("frequency_hz")
            if len(curve) == 0:
                continue
            grids[model] = curve
        if len(grids) < 2:
            continue

        freq_arrays = {m: c["frequency_hz"].to_numpy() for m, c in grids.items()}
        reference_model = models[0]
        grids_identical = all(
            len(freq_arrays[m]) == len(freq_arrays[reference_model])
            and np.allclose(freq_arrays[m], freq_arrays[reference_model])
            for m in freq_arrays
        )
        print(f"  {tree}: grids identical across models = {grids_identical}")

        if grids_identical:
            common_freq = freq_arrays[reference_model]
            amplitudes = {m: grids[m][RESPONSE_COLUMN].abs().to_numpy() for m in grids}
        else:
            lo = max(a.min() for a in freq_arrays.values())
            hi = min(a.max() for a in freq_arrays.values())
            common_freq = np.linspace(lo, hi, 200)
            print(f"    RESAMPLED onto a common {lo:.4f}-{hi:.4f} Hz grid (200 points, linear "
                  f"interpolation) because the raw grids differed.")
            amplitudes = {
                m: np.interp(common_freq, freq_arrays[m], grids[m][RESPONSE_COLUMN].abs().to_numpy())
                for m in grids
            }

        for model, amp in amplitudes.items():
            auc = float(np.trapezoid(amp, common_freq))
            rows.append({"tree": tree, "model": model, "auc": auc,
                        "freq_lo_hz": float(common_freq.min()), "freq_hi_hz": float(common_freq.max()),
                        "n_points": len(common_freq), "grids_identical": grids_identical})
    return pd.DataFrame(rows)


def compute_auc_deviation(auc_df: pd.DataFrame, model_a: str, model_b: str) -> pd.DataFrame:
    """
    rel_dev_pct = 100 * (AUC[model_b] - AUC[model_a]) / AUC[model_a]

    model_a/model_b come from AUC_DEVIATION_MODEL_A/_B (CONFIGURATION) - flip
    those two constants to flip the sign convention; this function's body
    never needs editing for that.
    """
    rows = []
    for tree in sorted(auc_df["tree"].unique()):
        tree_auc = auc_df[auc_df["tree"] == tree].set_index("model")["auc"]
        if model_a not in tree_auc.index or model_b not in tree_auc.index:
            print(f"  MISSING for {tree}: cannot compute deviation without both '{model_a}' and '{model_b}'.")
            continue
        auc_a, auc_b = tree_auc[model_a], tree_auc[model_b]
        rel_dev_pct = 100.0 * (auc_b - auc_a) / auc_a if auc_a != 0 else float("nan")
        rows.append({"tree": tree, "model_a": model_a, "auc_a": auc_a,
                    "model_b": model_b, "auc_b": auc_b, "rel_dev_pct": rel_dev_pct})
    return pd.DataFrame(rows)


# --- C1: the two spectrum figures -------------------------------------------

def build_combined_spectrum_figure(harmonic_df: pd.DataFrame, tree_order: list, models: list,
                                   model_colors: dict, newest_input_mtime, filters_text: str) -> None:
    """
    ONE axes, both models. SPECTRUM_COMBINED_STYLE picks the layout - see
    its own comment. Stage 8's verdict on the old "all_curves" layout was
    that SPECTRUM_NORMALIZE helps a little (every curve becomes visible)
    but SPECTRUM_LOG_Y does not, and neither solves 22 overlapping lines on
    one axes - so "summary" (median + spread per model) is now the
    default; "all_curves" is kept, not removed, behind the config switch.
    """
    assert_gravity_variant(harmonic_df, GRAVITY_VARIANT, "build_combined_spectrum_figure")
    fig, ax = plt.subplots(1, 1, figsize=(12, 7))
    ax.set_axisbelow(True)
    ax.grid(color="#e1e0d9", linewidth=0.6, zorder=0)

    if SPECTRUM_COMBINED_STYLE == "all_curves":
        for tree in tree_order:
            for model in models:
                curve = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)].sort_values("frequency_hz")
                if len(curve) == 0:
                    continue
                amplitude = curve[RESPONSE_COLUMN].abs()
                if SPECTRUM_NORMALIZE and amplitude.max() > 0:
                    amplitude = amplitude / amplitude.max()
                ax.plot(curve["frequency_hz"], amplitude, color=model_colors[model], alpha=0.45, linewidth=1.2)
        if SPECTRUM_LOG_Y:
            ax.set_yscale("log")
        ylabel = f"|{RESPONSE_COLUMN}|" + (" (normalized to each curve's own peak)" if SPECTRUM_NORMALIZE else " (raw amplitude)")
        title = f"Frequency spectra - all {len(tree_order)} trees, {len(models)} models (combined axes)"

    else:  # "summary" (default) - see SPECTRUM_COMBINED_STYLE's own comment
        for model in models:
            # Each tree's curve normalized to ITS OWN peak first (always -
            # not gated by SPECTRUM_NORMALIZE), otherwise the tree with the
            # single largest raw peak would dominate the across-tree median
            # by sheer magnitude rather than shape.
            per_tree_curves = []
            for tree in tree_order:
                curve = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)].sort_values("frequency_hz")
                if len(curve) == 0:
                    continue
                amplitude = curve[RESPONSE_COLUMN].abs()
                peak = amplitude.max()
                if peak > 0:
                    amplitude = amplitude / peak
                per_tree_curves.append(pd.Series(amplitude.to_numpy(),
                                                 index=curve["frequency_hz"].to_numpy(), name=tree))
            if not per_tree_curves:
                continue
            # Align every tree's curve on the shared frequency grid (a tree
            # missing a point contributes NaN there, not a crash), then take
            # the median/25th/75th percentile ACROSS TREES at each
            # frequency, ignoring any NaN.
            aligned = pd.concat(per_tree_curves, axis=1)
            median_curve = aligned.median(axis=1, skipna=True)
            band_low = aligned.quantile(0.25, axis=1)
            band_high = aligned.quantile(0.75, axis=1)
            ax.fill_between(aligned.index, band_low, band_high, color=model_colors[model],
                           alpha=0.25, linewidth=0)
            ax.plot(aligned.index, median_curve, color=model_colors[model], linewidth=2.2)
        if SPECTRUM_LOG_Y:
            ax.set_yscale("log")
        # Short on purpose - at AXIS_LABEL_FONTSIZE (15) the longer original
        # wording ("...each tree normalized to its own peak, then median
        # across trees") was taller than the axes and overflowed past the
        # top of the figure. The full explanation is still in the title.
        ylabel = f"|{RESPONSE_COLUMN}| (normalized, median across trees)"
        title = (f"Frequency spectra summary - median per model across n={len(tree_order)} trees "
                f"(shaded band: 25th-75th percentile)")

    ax.set_xlabel("Frequency (Hz)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE)
    ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
    legend_handles = [Line2D([0], [0], color=model_colors[m], linewidth=2, label=m) for m in models]
    ax.legend(handles=legend_handles, title="Model", loc="upper right", frameon=False,
             fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_FONTSIZE)

    add_provenance_footer(fig, "plot_frequency_spectra.py", newest_input_mtime, filters_text,
                          n_trees=len(tree_order), n_models=len(models))
    plt.tight_layout()
    savefig_with_provenance(fig, SPECTRA_COMBINED_PNG_PATH)
    plt.close(fig)
    print(f"\nCombined spectrum figure ({SPECTRUM_COMBINED_STYLE}) saved to: {SPECTRA_COMBINED_PNG_PATH}")


def build_grid_spectrum_figure(harmonic_df: pd.DataFrame, tree_order: list, models: list,
                               model_colors: dict, newest_input_mtime, filters_text: str) -> None:
    """One small panel per tree, both models, shared axes - usually the readable layout."""
    assert_gravity_variant(harmonic_df, GRAVITY_VARIANT, "build_grid_spectrum_figure")
    n_trees = len(tree_order)
    n_cols = 3
    n_rows = -(-n_trees // n_cols)  # ceiling division
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.2 * n_rows),
                             sharex=True, sharey=not SPECTRUM_NORMALIZE)
    axes = np.atleast_1d(axes).flatten()

    for i, tree in enumerate(tree_order):
        ax = axes[i]
        ax.set_axisbelow(True)
        ax.grid(color="#e1e0d9", linewidth=0.6, zorder=0)
        for model in models:
            curve = harmonic_df[(harmonic_df["tree"] == tree) & (harmonic_df["model"] == model)].sort_values("frequency_hz")
            if len(curve) == 0:
                continue
            amplitude = curve[RESPONSE_COLUMN].abs()
            if SPECTRUM_NORMALIZE and amplitude.max() > 0:
                amplitude = amplitude / amplitude.max()
            ax.plot(curve["frequency_hz"], amplitude, color=model_colors[model], linewidth=1.3)
        if SPECTRUM_LOG_Y:
            ax.set_yscale("log")
        ax.set_title(tree, fontsize=TITLE_FONTSIZE)
        ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)

    for j in range(n_trees, len(axes)):
        axes[j].set_axis_off()

    fig.supxlabel("Frequency (Hz)", fontsize=AXIS_LABEL_FONTSIZE)
    ylabel = f"|{RESPONSE_COLUMN}|" + (" (normalized)" if SPECTRUM_NORMALIZE else "")
    fig.supylabel(ylabel, fontsize=AXIS_LABEL_FONTSIZE)
    legend_handles = [Line2D([0], [0], color=model_colors[m], linewidth=2, label=m) for m in models]
    # Placed inside the first UNUSED grid slot (there is always at least
    # one whenever n_trees is not an exact multiple of n_cols) rather than
    # as a figure-level legend at a guessed coordinate - a fixed coordinate
    # collided with the last real panel's title when the grid had no spare
    # row (e.g. exactly 3, 6, or 9 trees).
    if n_trees < len(axes):
        axes[n_trees].legend(handles=legend_handles, title="Model", loc="center", frameon=False,
                             fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_FONTSIZE)
    else:
        fig.legend(handles=legend_handles, title="Model", loc="upper right", frameon=False,
                  bbox_to_anchor=(0.99, 0.99), fontsize=LEGEND_FONTSIZE, title_fontsize=LEGEND_FONTSIZE)
    fig.suptitle(f"Frequency spectra - one panel per tree (n={n_trees})", fontsize=TITLE_FONTSIZE)

    add_provenance_footer(fig, "plot_frequency_spectra.py", newest_input_mtime, filters_text,
                          n_trees=n_trees, n_models=len(models))
    plt.tight_layout()
    savefig_with_provenance(fig, SPECTRA_GRID_PNG_PATH)
    plt.close(fig)
    print(f"Grid spectrum figure saved to: {SPECTRA_GRID_PNG_PATH}")


def build_auc_figure(deviation_df: pd.DataFrame, tree_order: list, newest_input_mtime, filters_text: str) -> None:
    """One bar per tree: relative deviation in area-under-curve between the two models."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e1e0d9", linewidth=0.6, zorder=0)
    ax.axhline(0.0, color="0.3", linewidth=1.0, zorder=2)

    plotted_trees = [t for t in tree_order if t in set(deviation_df["tree"])]
    values = [deviation_df.set_index("tree").loc[t, "rel_dev_pct"] for t in plotted_trees]
    ax.bar(range(len(plotted_trees)), values, color="#7E9FC7", edgecolor="black", linewidth=0.8, zorder=3)
    ax.set_xticks(range(len(plotted_trees)))
    ax.set_xticklabels(plotted_trees, rotation=45, ha="right", fontsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE)

    model_a = AUC_DEVIATION_MODEL_A
    model_b = AUC_DEVIATION_MODEL_B
    # Stage 8 Part C1: the AUC formula used to live in the (rotated) y-axis
    # label, where it was long enough to poke up past the top of the axes
    # and overlap the title text above it. Moved into the title block
    # instead (plain horizontal text has much more width to work with),
    # leaving the y-axis label short enough to always fit within the axes.
    ax.set_ylabel("Relative deviation (%)", fontsize=AXIS_LABEL_FONTSIZE)
    # This whole block is a multi-line PARAGRAPH (title + formula + caveat),
    # not a short label - TITLE_FONTSIZE (15) here forced the figure far
    # wider than its own bar chart (bbox_inches="tight" expanding the
    # canvas to fit one very long wrapped line). Kept smaller than
    # TITLE_FONTSIZE on purpose, still up from the original 8.
    ax.set_title(f"Area under the harmonic response curve - deviation between models\n"
                f"Relative deviation = 100 x (AUC[{model_b}] - AUC[{model_a}]) / AUC[{model_a}]\n"
                f"CAVEAT: area measures total response magnitude across the swept band, "
                f"NOT agreement in peak frequency - two curves can share an area with "
                f"peaks at different frequencies. See the printed C3 validity check before "
                f"trusting this figure for any tree it flags INVALID.", fontsize=11)

    add_provenance_footer(fig, "plot_frequency_spectra.py", newest_input_mtime, filters_text,
                          n_trees=len(plotted_trees), n_models=2)
    plt.tight_layout()
    savefig_with_provenance(fig, SPECTRA_AUC_PNG_PATH)
    plt.close(fig)
    print(f"AUC deviation figure saved to: {SPECTRA_AUC_PNG_PATH}")


def main():
    print("=" * 78)
    print("plot_frequency_spectra.py")
    print("=" * 78)

    print_active_configuration("plot_frequency_spectra.py", CONFIG_CONSTANT_NAMES)

    newest_input_mtime = check_input_freshness(
        "plot_frequency_spectra.py", INPUT_FOLDER, "out_harmonic_bk*_scenario*.txt", {}, ABORT_ON_STALE_INPUT,
    )

    order_to_tree_id = load_ansys_order_to_tree_id(GEOMETRY_FILE_PATH)
    harmonic_df = load_harmonic_data(INPUT_FOLDER, HARMONIC_FILE_PATTERN, order_to_tree_id)
    if len(harmonic_df) == 0:
        print("No harmonic files found - nothing to plot.")
        return

    all_models_found = sorted(harmonic_df["model"].unique())
    model_colors = assign_model_colors(all_models_found)  # from ALL models, so colours match every other figure

    # --- naming check, per the task's own explicit caution -----------------
    no_acel_models = sorted(m for m in all_models_found if is_no_acel_model(m))
    print(f"\nDistinct models found in the harmonic data: {all_models_found}")
    print(f"'_no_acel' models found: {no_acel_models}")
    if len(no_acel_models) != 2:
        raise SystemExit(
            f"plot_frequency_spectra.py: expected EXACTLY 2 '_no_acel' models in the harmonic "
            f"data, found {len(no_acel_models)}: {no_acel_models}. Stopping rather than "
            f"guessing which two are meant - see this script's module docstring."
        )

    harmonic_df_drawing = filter_by_gravity_variant(harmonic_df, GRAVITY_VARIANT)
    models = sorted(harmonic_df_drawing["model"].unique())
    tree_order = sorted(harmonic_df_drawing["tree"].unique())
    print(f"GRAVITY_VARIANT={GRAVITY_VARIANT!r} - drawing models: {models}, "
          f"{len(tree_order)} tree(s): {tree_order}")

    filters_text = f"GRAVITY_VARIANT={GRAVITY_VARIANT}  SPECTRUM_NORMALIZE={SPECTRUM_NORMALIZE}  SPECTRUM_LOG_Y={SPECTRUM_LOG_Y}"

    # --- C1 ---
    build_combined_spectrum_figure(harmonic_df_drawing, tree_order, models, model_colors,
                                   newest_input_mtime, filters_text)
    build_grid_spectrum_figure(harmonic_df_drawing, tree_order, models, model_colors,
                               newest_input_mtime, filters_text)

    # --- C3 (must run, and be read, before C2) ---
    validity_df = check_spectrum_validity(harmonic_df_drawing, tree_order, models)

    # --- C2 ---
    print("\n" + "=" * 78)
    print("C2 - AREA UNDER THE HARMONIC RESPONSE CURVE")
    print("=" * 78)
    auc_df = compute_auc(harmonic_df_drawing, tree_order, models)
    deviation_df = compute_auc_deviation(auc_df, AUC_DEVIATION_MODEL_A, AUC_DEVIATION_MODEL_B)

    print(f"\nRaw AUC values (units: |{RESPONSE_COLUMN}| x Hz - trapezoid integral of the "
          f"absolute harmonic response amplitude over the swept frequency range):")
    print(auc_df.to_string(index=False))
    print(f"\nDeviation table (100 x (AUC[{AUC_DEVIATION_MODEL_B}] - AUC[{AUC_DEVIATION_MODEL_A}]) "
          f"/ AUC[{AUC_DEVIATION_MODEL_A}]):")
    print(deviation_df.to_string(index=False))

    invalid_trees = set(validity_df.loc[~validity_df["valid"], "tree"]) if len(validity_df) > 0 else set()
    if invalid_trees:
        print(f"\n*** {len(invalid_trees)} tree(s) flagged INVALID by the C3 check above - "
              f"their AUC deviation numbers are printed and saved anyway (never silently "
              f"dropped), but should NOT be treated as a trustworthy area comparison: "
              f"{sorted(invalid_trees)} ***")

    build_auc_figure(deviation_df, tree_order, newest_input_mtime, filters_text)

    merged = auc_df.merge(validity_df[["tree", "model", "valid", "peak_freq_hz", "edge_to_peak_ratio"]],
                          on=["tree", "model"], how="left")
    merged.to_csv(SPECTRA_AUC_CSV_PATH, index=False)
    print(f"\nAUC + validity CSV saved to: {SPECTRA_AUC_CSV_PATH}")


if __name__ == "__main__":
    main()
