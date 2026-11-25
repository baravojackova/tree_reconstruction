"""
make_size_category_plot.py
---------------------------
Displacement stiffness deviation against tree size, with the dataset's own
size categories marked.

Why this plot: the deviation from measurement is not random across trees - it
tracks tree size. Large trees sit close to the measurement, small trees are
consistently too soft. One tree (B21_S11) breaks the pattern, and it is also
the tree with by far the worst trunk taper error, which is shown as a ring
around its point.

Reads taper_vs_stiffness_panels-noacel.csv (produced by
method_deviation_overview.py) and writes a PNG and an SVG. Nothing else is
written and the input is never modified.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# =========================================================================
# CONFIGURATION - edit these values, then run the script
# =========================================================================

# Path(__file__).parent is the folder this script sits in, so the paths hold
# wherever the terminal happens to be - and they survive a move to another
# machine without editing. The CSVs and figures live in ansys_eval/out/.
OUT_FOLDER = Path(__file__).parent / "out"

INPUT_CSV = OUT_FOLDER / "taper_vs_stiffness_panels-noacel.csv"
OUTPUT_PNG = OUT_FOLDER / "stiffness_by_size_category.png"
OUTPUT_SVG = OUT_FOLDER / "stiffness_by_size_category.svg"

# Field DBH per tree, in cm, from the run-plan table in the handoff.
DBH_CM = {
    "B21_S12": 30, "B21_S16": 33, "B21_S09": 36, "B21_S10": 37,
    "B21_S08": 39, "B21_S14": 41, "B21_S04": 43, "B21_S11": 46,
    "B21_S13": 49, "B21_S21": 49, "B21_S01": 56, "B21_S07": 64,
}

# The dataset's own size split, in cm. Trees at or above this go in the
# larger category.
SIZE_SPLIT_CM = 45

# Trees to ring as taper outliers, with the reason shown in the legend.
TAPER_OUTLIER_THRESHOLD_PCT = 15.0

MODEL_STYLE = {
    "adtree":  {"label": "adtree_r5mm_no_acel",   "fill": "#F5B399", "edge": "#993C1D"},
    "treeqsm": {"label": "treeqsm_pd002_no_acel", "fill": "#F6D07F", "edge": "#854F0B"},
}

BAND_SMALL_COLOR = "#F1EFE8"
BAND_LARGE_COLOR = "#E4E2D8"

POINT_SIZE = 90
FIGSIZE = (10, 6.5)
DPI = 200
LABEL_FONTSIZE = 8
AXIS_LABEL_FONTSIZE = 11
TITLE_FONTSIZE = 13
LEGEND_FONTSIZE = 9

SHOW_TREND_LINE = True


# =========================================================================
# CODE - you normally do not need to change anything below this line
# =========================================================================

def load() -> pd.DataFrame:
    """
    Load the panel CSV and reshape it to one row per tree, with a column
    per model for both the taper error and the stiffness deviation.
    """
    if not Path(INPUT_CSV).is_file():
        print(f"ERROR: input not found: {INPUT_CSV}")
        print("CSV files present in that folder:")
        for candidate in sorted(Path(INPUT_CSV).parent.glob("*.csv")):
            print(f"  {candidate.name}")
        raise SystemExit(1)

    df = pd.read_csv(INPUT_CSV)
    wide = df.pivot_table(index="tree", columns="taper_family",
                          values=["taper_error_pct", "displacement_rel_dev_pct"])
    wide.columns = ["dev_adtree", "dev_treeqsm", "taper_adtree", "taper_treeqsm"]
    wide["dbh_cm"] = [DBH_CM[t] for t in wide.index]
    return wide.reset_index()


def draw() -> None:
    data = load()

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # --- Size-category bands, drawn first so everything sits on top -----
    x_min, x_max = 27, 68
    ax.axvspan(x_min, SIZE_SPLIT_CM, color=BAND_SMALL_COLOR, zorder=0)
    ax.axvspan(SIZE_SPLIT_CM, x_max, color=BAND_LARGE_COLOR, zorder=0)
    ax.axvline(SIZE_SPLIT_CM, color="#c3c2b7", linewidth=1.2, zorder=1)
    ax.axhline(0, color="#898781", linewidth=1.2, linestyle="--", zorder=1)

    ax.text((x_min + SIZE_SPLIT_CM) / 2, 22, f"smaller category (< {SIZE_SPLIT_CM} cm)",
            ha="center", fontsize=9, color="#898781")
    ax.text((SIZE_SPLIT_CM + x_max) / 2, 22, f"larger category (> {SIZE_SPLIT_CM} cm)",
            ha="center", fontsize=9, color="#898781")

    for key, style in MODEL_STYLE.items():
        x = data["dbh_cm"].values
        y = data["dev_" + key].values

        ax.scatter(x, y, s=POINT_SIZE, facecolor=style["fill"],
                   edgecolor=style["edge"], linewidth=1.2, zorder=3)

        # Ring the trees whose taper error marks them out as badly
        # reconstructed, so a point that breaks the trend is explained on
        # the figure rather than in a caption.
        outliers = data[data["taper_" + key] > TAPER_OUTLIER_THRESHOLD_PCT]
        if len(outliers):
            ax.scatter(outliers["dbh_cm"], outliers["dev_" + key],
                       s=POINT_SIZE * 3.2, facecolor="none", edgecolor="#993C1D",
                       linewidth=1.6, linestyle=(0, (3, 2)), zorder=4)

        if SHOW_TREND_LINE:
            slope, intercept = np.polyfit(x, y, 1)
            r = np.corrcoef(x, y)[0, 1]
            xs = np.array([x_min + 1, x_max - 1])
            ax.plot(xs, slope * xs + intercept, color=style["edge"],
                    linewidth=1.6, alpha=0.65, zorder=2)
            ax.text(x_max - 1.5, slope * (x_max - 1) + intercept + 2.5,
                    f"R² = {r*r:.2f}", color=style["edge"], fontsize=9,
                    ha="right", va="bottom")

    # Tree labels - one under EACH model's point (not just the lower one),
    # so a treeqsm point can be identified even when it sits above its
    # adtree twin.
    for _, row in data.iterrows():
        tree_label = row["tree"].replace("B21_", "")
        for key in MODEL_STYLE:
            ax.annotate(tree_label, (row["dbh_cm"], row["dev_" + key]),
                        textcoords="offset points", xytext=(0, -14), ha="center",
                        fontsize=LABEL_FONTSIZE, color=MODEL_STYLE[key]["edge"])

    ax.set_xlim(x_min, x_max)
    ax.set_xlabel("field DBH (cm)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("displacement stiffness deviation (%)\nbelow zero = model softer than measured",
                  fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title("Stiffness deviation tracks tree size (12 trees, gravity off)",
                 fontsize=TITLE_FONTSIZE, pad=12)

    ax.grid(axis="y", color="#e1e0d9", linewidth=0.8)
    ax.set_axisbelow(False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c3c2b7")

    handles = [Patch(facecolor=s["fill"], edgecolor=s["edge"], label=s["label"])
               for s in MODEL_STYLE.values()]
    handles.append(Line2D([], [], marker="o", linestyle="", markerfacecolor="none",
                          markeredgecolor="#993C1D", markersize=12,
                          label=f"taper error > {TAPER_OUTLIER_THRESHOLD_PCT:.0f} %"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              frameon=False, fontsize=LEGEND_FONTSIZE)

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=DPI, bbox_inches="tight")
    fig.savefig(OUTPUT_SVG, bbox_inches="tight")
    print(f"Saved: {OUTPUT_PNG}\nSaved: {OUTPUT_SVG}\n")

    # Group means, so the figure can be checked against numbers.
    for key, style in MODEL_STYLE.items():
        small = data[data["dbh_cm"] < SIZE_SPLIT_CM]["dev_" + key]
        large = data[data["dbh_cm"] >= SIZE_SPLIT_CM]["dev_" + key]
        r = np.corrcoef(data["dbh_cm"], data["dev_" + key])[0, 1]
        print(f"{style['label']:24} smaller n={len(small):2} median {small.median():7.1f} | "
              f"larger n={len(large):2} median {large.median():7.1f} | r={r:.3f} R2={r*r:.3f}")


if __name__ == "__main__":
    draw()
