"""
make_boxplot.py
----------------
Box plot of the signed deviation from measurement, across all trees, for the
two gravity-off FEM models.

Reads the two deviation tables produced by the plotting pipeline and writes a
PNG and an SVG. Nothing else is written; the input CSVs are never modified.

Three quantities are shown side by side:
  - displacement stiffness  (slope_displacement, combined sensor)
  - rotation stiffness      (slope_rotation_base, combined sensor)
  - harmonic peak frequency (harmonic_peak, elasto variant)

The deviation is SIGNED, not absolute: a box below zero means the model is
consistently softer / slower than the measurement, which is a different
statement from "the model is off by X %".
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # no display available; write straight to file

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Path(__file__).parent is the folder this script sits in, so the paths hold
# wherever the terminal happens to be - and they survive a move to another
# machine without editing. The CSVs and figures live in ansys_eval/out/.
OUT_FOLDER = Path(__file__).parent / "out"

STIFFNESS_CSV = OUT_FOLDER / "stiffness_deviation.csv"
FREQUENCY_CSV = OUT_FOLDER / "frequency_deviation.csv"

OUTPUT_PNG = OUT_FOLDER / "model_comparison_boxplot.png"
OUTPUT_SVG = OUT_FOLDER / "model_comparison_boxplot.svg"

SENSOR = "combined"             # which experimental sensor variant to use
FREQUENCY_VARIANT = "elasto"    # which experimental frequency variant to use
ONLY_NO_ACEL = True             # True = gravity-off runs only

# Pastel fills with a darker edge from the same colour family, following the
# project convention (never a plain black edge).
MODEL_STYLE = {
    "adtree":  {"label": "adtree_r5mm_no_acel",   "fill": "#F5B399", "edge": "#993C1D"},
    "treeqsm": {"label": "treeqsm_pd002_no_acel", "fill": "#F6D07F", "edge": "#854F0B"},
}

BOX_WIDTH = 0.28
POINT_SIZE = 26
POINT_ALPHA = 0.75
FIGSIZE = (12, 5.5)

# Match method_deviation_overview.py: symlog y axis, one panel per quantity.
# Linear axis puts the +186 % outlier at the top and squeezes everything else
# into the bottom third; symlog keeps the band around zero readable.
USE_SYMLOG = True
SYMLOG_LINTHRESH = 10.0
DPI = 200

TITLE_FONTSIZE = 13
AXIS_LABEL_FONTSIZE = 11
TICK_LABEL_FONTSIZE = 10
LEGEND_FONTSIZE = 10


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def load_groups() -> list:
    """
    Return a list of (panel_label, {model_key: array of deviations}).

    One entry per quantity. Each inner dict holds the per-tree deviation
    values for that quantity and model, so the box and the individual
    points come from exactly the same numbers.
    """
    for path in (STIFFNESS_CSV, FREQUENCY_CSV):
        if not Path(path).is_file():
            print(f"ERROR: input not found: {path}")
            print("CSV files present in that folder:")
            for candidate in sorted(Path(path).parent.glob("*.csv")):
                print(f"  {candidate.name}")
            raise SystemExit(1)

    stiffness = pd.read_csv(STIFFNESS_CSV)
    frequency = pd.read_csv(FREQUENCY_CSV)

    if ONLY_NO_ACEL:
        stiffness = stiffness[stiffness["model"].str.contains("no_acel")]
        frequency = frequency[frequency["model"].str.contains("no_acel")]

    panels = []

    for column, label in [("slope_displacement", "Displacement stiffness"),
                          ("slope_rotation_base", "Rotation stiffness")]:
        subset = stiffness[(stiffness["fem_column"] == column)
                           & (stiffness["sensor"] == SENSOR)]
        panels.append((label, {
            key: subset[subset["model"].str.contains(key)]["rel_dev_pct"].values
            for key in MODEL_STYLE
        }))

    subset = frequency[(frequency["quantity"] == "harmonic_peak")
                       & (frequency["variant"] == FREQUENCY_VARIANT)]
    panels.append(("Harmonic peak frequency", {
        key: subset[subset["model"].str.contains(key)]["rel_dev_pct"].values
        for key in MODEL_STYLE
    }))

    return panels


def draw() -> None:
    panels = load_groups()

    fig, axes = plt.subplots(1, len(panels), figsize=FIGSIZE, sharey=True)

    n_trees = 0
    for ax, (panel_label, by_model) in zip(axes, panels):
        ax.axhline(0, color="#898781", linewidth=1.2, linestyle="--", zorder=1)

        for model_index, (key, style) in enumerate(MODEL_STYLE.items()):
            values = by_model[key]
            n_trees = max(n_trees, len(values))
            x = model_index

            box = ax.boxplot(values, positions=[x], widths=BOX_WIDTH,
                             patch_artist=True, showfliers=False, zorder=2)
            box["boxes"][0].set(facecolor=style["fill"], edgecolor=style["edge"],
                                linewidth=1.4, alpha=0.85)
            for part in ("whiskers", "caps"):
                for line in box[part]:
                    line.set(color=style["edge"], linewidth=1.4)
            box["medians"][0].set(color=style["edge"], linewidth=2.2)

            rng = np.random.default_rng(0)
            jitter = rng.uniform(-0.04, 0.04, size=len(values))
            ax.scatter(np.full(len(values), x + 0.22) + jitter, values,
                       s=POINT_SIZE, color=style["edge"], alpha=POINT_ALPHA,
                       edgecolors="none", zorder=3)

        if USE_SYMLOG:
            # symlog is linear within +-SYMLOG_LINTHRESH and logarithmic
            # outside it, so both the band near zero and the outliers stay
            # readable on one axis.
            ax.set_yscale("symlog", linthresh=SYMLOG_LINTHRESH)

        ax.set_xlim(-0.5, 1.6)
        ax.set_xticks([])
        ax.set_title(panel_label, fontsize=AXIS_LABEL_FONTSIZE)
        ax.grid(axis="y", color="#e1e0d9", linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#c3c2b7")

    axes[0].set_ylabel("deviation from measurement (%)\nbelow zero = model softer / slower",
                       fontsize=AXIS_LABEL_FONTSIZE)
    axes[0].tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE)

    fig.suptitle(f"FEM vs measurement, {n_trees} trees, gravity off",
                 fontsize=TITLE_FONTSIZE, y=1.0)

    handles = [Patch(facecolor=s_["fill"], edgecolor=s_["edge"], label=s_["label"])
               for s_ in MODEL_STYLE.values()]
    axes[-1].legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=False, fontsize=LEGEND_FONTSIZE)

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTPUT_SVG, bbox_inches="tight")
    print(f"Saved: {OUTPUT_PNG}")
    print(f"Saved: {OUTPUT_SVG}")

    for panel_label, by_model in panels:
        print(f"\n{panel_label}")
        for key, style in MODEL_STYLE.items():
            v = by_model[key]
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            print(f"  {style['label']:24} n={len(v):2}  min {v.min():7.1f}  "
                  f"q1 {q1:7.1f}  median {med:7.1f}  q3 {q3:7.1f}  max {v.max():7.1f}")


if __name__ == "__main__":
    draw()
