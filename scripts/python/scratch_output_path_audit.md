
# Output Path Audit — scripts/python

READ-ONLY static analysis, prepared as input for a future output-path reorganisation.
No files were run, imported, moved, deleted, or edited to produce this report — everything
below comes from reading source files and listing the filesystem. This document does not
propose any new layout; it only records the current state.

Working-directory convention observed throughout the codebase: scripts are written assuming
`scripts/python` is the process's current working directory. A bare relative literal such as
`"volume_results.csv"` or `"plots"` therefore resolves under `scripts/python`, not the repo root,
*only if the script is actually launched from that directory*.

---

## STEP 1 — Inventory

### 1a. `scripts/python` root — `.py` files and one-line purpose (from each file's own header comment)

| File | Purpose (paraphrased from the file's own header) |
|---|---|
| `adqsm_build_median_variant.py` | Build ONE new synthetic "median" AdQSM variant folder for a tree from a configurable list of already-reconstructed source variants, by taking the median of each source variant's taper.txt/BranchStructure.txt/TreesParams.txt row-by-row. Only creates `data/<TREE>/<NEW_VARIANT>/` with 3 new files; refuses to run (SystemExit) if that folder already exists. |
| `adqsm_variant_sensitivity.py` | Reproduce/save the AdQSM-variant-sensitivity analysis: for a dense grid of AdQSM segmentation-setting variants on one tree, overlay AdQSM's own total volume against AdTree's calibrated total volume on one log-scale chart, saved as CSV/PNG for later review. |
| `adtree_reconstruct_compare.py` | Calibrate an AdTree skeleton (.ply) against AdQSM, compare it to other reconstruction methods, and save the calibrated geometry. Step 1 of a 2-step pipeline (see `export_geom_ansys.py` for step 2). Saves final calibrated geometry to `calib_<tree>_r<mm>mm<variant>.npz` per threshold/variant; no longer writes `geom_*.txt` itself. |
| `calmethod_decision_summary.py` | Decide between AdTree's two radius-calibration methods (`calmethod="min5mm"` vs `"regression-perorder"`) as the primary calibration going forward, using an objective, decomposed (total/trunk/branch, separately) comparison. Produces OUTPUT_CSV, three PNG charts, and a printed recommendation. |
| `check_adqsm_thin_branch_length.py` | Standalone diagnostic: how much of AdQSM's branch length/volume sits in very thin branches, binned by diameter. Explicitly states it does NOT touch `volume_results.csv` or any other shared pipeline file — reads BranchStructure.txt/TreesParams.txt only and prints/plots a report. |
| `compare_volumes.py` | Compares tree-volume results from different reconstruction methods against a reference. Reads the "master" CSV `volume_results.csv`, prints side-by-side comparison plus Bias/MAE/RMSE/CV-RMSE metrics. |
| `export_geom_ansys.py` | Step 2 of the AdTree→ANSYS pipeline: converts one `calib_*.npz` file (written by `adtree_reconstruct_compare.py`) into the `geom_*.txt` file ANSYS (`buk.mac`) reads, via `write_geom()` from `tree_geom_utils.py`. Does no calibration/comparison itself. |
| `import_matlab_results.py` | Imports the author's own TreeQSM MATLAB results (`volumes_<tree>_<run>.csv`, written by `run_treeqsm.m`) into the shared master results table `volume_results.csv`. |
| `parameter_sensitivity.py` | Parameter-sensitivity analysis of how AdTree's or TreeQSM's input parameters affect reconstructed outputs and which combination lands closest to reference — three chart types plus pooled multi-tree variants, tornado diagrams, ANOVA/eta-squared chart. |
| `plot_box.py` | Box plots comparing GROUPS of related method variants (not individual methods) for one tree, across the same 9 metric fields `plot_volumes.py`'s tree-overview shows; grouping configured via `GROUP_RULES`. |
| `plot_param_sweep.py` | Plot how TreeQSM's reconstruction metrics respond to ONE swept parameter (e.g. `pd2min`) across the small number of distinct reconstructions actually run for one tree — deliberately separate from `plot_box.py` because n is too small for box plots/correlations. |
| `plot_style.py` | Shared visual style for every plotting script — one place to change a color/line width/font size so every chart updates together. |
| `plot_volumes.py` | Draw and save charts summarizing `volume_results.csv` (the shared master results table). All PNGs written into a `plots` subfolder next to this script (created automatically), path of each saved file printed to console. |
| `qsm_volume_mean.py` | Average QSM volume of one tree from its 20 TreeQSM realization files — reads every realization, computes total/trunk/branch volume per realization, prints mean/std/CV. |
| `reconstruction_method_decision_summary.py` | Compare reconstruction APPROACHES (AdTree raw/calibrated, TreeQSM manual/auto, TreeQSM de Tanago, other AdQSM variants) against the appropriate reference in 3 separate blocks. Deliberately NOT a calibration-method decision (see `calmethod_decision_summary.py`). |
| `reference_volume.py` | Extract the destructive reference volume of one tree from the field measurement table (e.g. `IND_h_trees.txt`, de Tanago dataset). No dependencies beyond the standard library. |
| `scratch_calib_sensitivity.py` | READ-ONLY diagnostic: reproduces `adtree_reconstruct_compare.py`'s pipeline through the per-order regression fit and calref=min5mm scalar-factor computation for TWO SEG_LEN settings. States explicitly: "Writes NOTHING to disk ... Console output only." |
| `scratch_order_alignment.py` | READ-ONLY diagnostic: tests whether AdQSM's branch order corresponds to the same physical category as AdTree's, or lags by a constant shift k. Header claims "Writes NOTHING to disk" — **this claim is contradicted by the actual code** (see Step 2). |
| `scratch_skeleton_diag.py` | READ-ONLY diagnostic: reproduces `adtree_reconstruct_compare.py`'s input-preparation stage (read_ply → merge_vertices → smooth_centerline), stopping before `convert()`, to measure skeleton edge-length spacing and raw uncalibrated trunk radius. States "Writes NOTHING to disk." Confirmed true. |
| `taper_curve_compare.py` | Diagnostic: overlay every available AdQSM variant's `taper.txt` curve for a given tree on one chart, for visual/numeric inspection of DBH/trunk-volume differences between variants. Variant discovery is automatic (scans the tree's data folder). |
| `tree_coverage_check.py` | Pure read/print-only diagnostic (per header): for tree X, which parameter combinations already have rows in `volume_results.csv`. Header states "This script never writes anything" — **this claim is contradicted by the actual code** (see Step 2). |
| `tree_geom_utils.py` | Shared helper functions for the AdTree→ANSYS geometry pipeline. Explicitly has NO "RUN" section — only defines functions, imported by `adtree_reconstruct_compare.py` and `export_geom_ansys.py` (and, per its TreeGraph-support section, also by the TreeGraph pipeline scripts). |
| `treegraph_add_to_comparison.py` | Add a TreeGraph result to the shared master results table `volume_results.csv` (same file `compare_volumes.py` reads), via `upsert_result()`. |
| `treegraph_export_geom.py` | Export a TreeGraph `.json` result to `geom_*.txt` for ANSYS (`bk1.mac`), via `write_geom_from_treegraph()`. Unlike the AdTree pipeline there is no calibration step. |
| `visualize_pointcloud_region.py` | (Header in Czech.) Visualization of a single tree's point cloud, focused on an arbitrary selected region (e.g. a suspicious point cluster causing implausibly large radii in AdTree/AdQSM/TreeQSM reconstruction). Saves three images into `OUTPUT_DIR`. |

### 1b. Subfolders of `scripts/python` (one level, contents only)

| Subfolder | Contents (one level) |
|---|---|
| `__pycache__/` | Compiled bytecode cache (`.pyc`), not a data output. |
| `archiv_bak_2026-09-06/` | Old numbered `.py.bak`/`.bak2`...`.bak8` snapshots of `adtree_reconstruct_compare.py`, `compare_volumes.py`, `parameter_sensitivity.py`, `plot_volumes.py`, `reconstruction_method_decision_summary.py`, plus an un-suffixed stray copy of `tree_coverage_check.py`. Manual backups of source code, not outputs. |
| `archive_plots_IND/` | Archived PNG plots for the IND* trees, both loose PNGs at its root and per-tree subfolders `IND01_054/`, `IND03_088/`, `IND07_083/`, `IND_all/`. |
| `archive_vol_res_IND/` | Archived snapshots of `volume_results*.csv` for the IND* trees (6 CSVs directly inside), plus subfolders `IND_supporting_res/` and `IND_sweep_branChOrders/`. |
| `backup_BK/` | Manual backup snapshots: `B21_S01/` (with `adqs_var_testing/`, `final_B21_S01/`, `solving_calib_issue/` subfolders plus 2 loose CSVs), `B21_S04/` (a `plots/` subfolder plus 2 loose CSVs incl. a full `volume_results.csv` snapshot), and a top-level `plots/` (90 archived PNGs matching the same `adtree_adqsm_radius_regression_perorder_*` pattern seen live in `plots/`). |
| `npz/` | Flat store of `calib_*.npz` cache files (see Step 7c: 360 files, ~4.48 GiB). |
| `plots/` | Live plot output tree — per-tree subfolders plus some loose files (see Step 7a/7b). |

### 1c. Whole-project search for `.m`, `.mac`, `.bat`, `.ps1` files

No `.bat` or `.ps1` files exist anywhere in the repository. Total non-Python script files: **93** (87 `.m`, 6 `.mac`).

**`scripts/ansys/` (.mac, 6 files):**
| File | Purpose (one line, if header states it) |
|---|---|
| `bk1.mac` | ANSYS beam model of the beech tree from the Babice scans, dynamic analysis, units [m,N,Pa,kg/m3]. |
| `bk1_qsm.mac` | Same as `bk1.mac` but reads a QSM-based geometry input. |
| `koef_mat.mac` | No header comment; appears to be a raw coefficient/index data table. |
| `stat.mac` | Static solution and post-processing for the "buk babice" model. |
| `trans.mac` | Transient solution and post-processing for the "strom_dyna" (tree dynamics) model. |
| `vysledkytrans.mac` | Writes transient-analysis results (displacements over time) to a text file. |

**`scripts/matlab/` and subfolders (87 `.m` files)** — full list retained in Step 6 discussion below for the files that matter to Python-output coordination; the complete inventory spans `scripts/matlab/` (top level: `create_input.m`, `estimate_precision.m`, `make_models.m`, `make_models_parallel.m`, `myfun.m`, `runsken.m`, `select_optimum.m`, `treeqsm.m`), `least_squares_fitting/` (11 files), `main_steps/` (10 files), `plotting/` (18 files), `tools/` (33 files), `triangulation/` (5 files). Of all 87, only `runsken.m` (top level) references any Python-output filename pattern — see Step 6.

---

## STEP 2 — Every write operation (scripts/python root `.py` files only)

Legend: "→" concrete resolved path relative to `scripts/python` where staticaly knowable; "(abs)" = absolute path; "(dyn)" = depends on a runtime/loop variable, general form shown.

### `adqsm_build_median_variant.py`
| Line | Expression | Resolves to |
|---|---|---|
| 156 | `open(out_path, "w", ...)`, `out_path=os.path.join(out_dir,"taper.txt")` | `<DATA_ROOT>\B21_S04\999\taper.txt` (abs, outside scripts/python) |
| 244 | `open(out_path, "w", ...)`, `out_path=os.path.join(out_dir,"BranchStructure.txt")` | `<DATA_ROOT>\B21_S04\999\BranchStructure.txt` (abs) |
| 370 | `open(out_path, "w", ...)`, `out_path=os.path.join(out_dir,"TreesParams.txt")` | `<DATA_ROOT>\B21_S04\999\TreesParams.txt` (abs) |
| 383 | `os.makedirs(out_dir)` | `<DATA_ROOT>\B21_S04\999\` — "the ONLY folder this script ever creates" per its own comment |

### `adqsm_variant_sensitivity.py`
| Line | Expression | Resolves to |
|---|---|---|
| 326 | `os.makedirs(out_dir, exist_ok=True)`, `out_dir=os.path.join(ensure_plots_dir(), tree_name)` | `plots\B21_S04\` (current `TREE_NAME`) |
| 328 | `fig.savefig(png_path, ...)`, `png_path=os.path.join(out_dir, "adqsm_variant_sensitivity_%s.png"%tree_name)` | `plots\B21_S04\adqsm_variant_sensitivity_B21_S04.png` |
| 345 | `open(csv_path, "w", ...)` + `csv.DictWriter`, `csv_path=os.path.join(out_dir, "adqsm_variant_sensitivity_%s.csv"%tree_name)` | `plots\B21_S04\adqsm_variant_sensitivity_B21_S04.csv` |

### `adtree_reconstruct_compare.py`
| Line | Expression | Resolves to |
|---|---|---|
| 402 | `os.makedirs(NPZ_DIR, exist_ok=True)` | `npz\` (`NPZ_DIR="npz"`, line 332) |
| 403 | `os.makedirs(FIGURES_DIR, exist_ok=True)` | `plots\B21_S04\` (`FIGURES_DIR=os.path.join("plots",TREE_NAME)`, line 333) |
| 971 | `np.savez(npz_name, ...)`, built lines 968-969 from `NPZ_DIR`, `TREE_NAME`, loop vars `thr`∈`RADIUS_THRESHOLDS`, `variant_suffix` (from `ADQSM_VARIANTS`), `SEG_VARIANT_SUFFIX` (from `SEG_LEN_MIN_LIST`/`SEG_LEN_K_LIST`) | (dyn) e.g. `npz\calib_B21_S04_r5mm_adqsm090_seg35-500-k50.npz` |
| 513, 1005, 1049, 1081, 1094, 1114, 1127 | `upsert_result(RESULTS_CSV, TREE_NAME, ...)` — write itself happens inside `tree_geom_utils.upsert_result()` | `volume_results.csv` (`RESULTS_CSV="volume_results.csv"`, line 323) |
| 1160 | `plot_model(..., os.path.join(FIGURES_DIR, out), ...)` — actual `fig.savefig` inside `tree_geom_utils.plot_model()`; `out` built line 823 from `TREE_NAME`/`thr`/`variant_suffix`/`SEG_VARIANT_SUFFIX` | (dyn) e.g. `plots\B21_S04\geom_B21_S04_r5mm_adqsm090_seg35-500-k50.png` — **note: this is the PNG twin of the geom filename, not the geom text file itself (see Step 5/6)** |

### `calmethod_decision_summary.py`
| Line | Expression | Resolves to |
|---|---|---|
| 166 | `os.makedirs(PLOTS_DIR)` (inside `ensure_plots_dir()`) | `plots\` (`PLOTS_DIR="plots"`, line 63) |
| 222 | `open(OUTPUT_CSV,"w",...)` + `csv.DictWriter` | `calmethod_decision_summary.csv` (line 62) |
| 288 | `fig.savefig(out_path,...)`, `out_path=os.path.join(ensure_plots_dir(), "calmethod_decision_%s.png"%field)`, `field`∈`FIELDS=["total","trunk","branch"]` | `plots\calmethod_decision_total.png`, `plots\calmethod_decision_trunk.png`, `plots\calmethod_decision_branch.png` |

### `check_adqsm_thin_branch_length.py`
| Line | Expression | Resolves to |
|---|---|---|
| 256 | `os.makedirs(PLOTS_DIR)` | `plots\` (`PLOTS_DIR="plots"`, line 72) |
| 273 | `fig.savefig(plot_path,...)`, `plot_path=os.path.join(PLOTS_DIR, PLOT_FILENAME)` | `plots\adqsm_thin_branch_length_share.png` (`PLOT_FILENAME`, line 73) |

### `compare_volumes.py`
| Line | Expression | Resolves to |
|---|---|---|
| 433-434 | `open(RESULTS_CSV,"w",...)` + `csv.writer(f).writerows(STARTER_ROWS)`, only if file doesn't exist yet (`__main__` guard) | `volume_results.csv` (line 71) |

### `export_geom_ansys.py`
No direct write call in the file. Line 114 `write_geom(geom_filename, ...)` triggers the actual `open(path,"w")` inside `tree_geom_utils.write_geom()`. `geom_filename` is loaded verbatim from the `.npz` metadata (`geom_filename=str(data["geom_filename"])`), which was stored as a **bare filename with no directory prefix** (`out`, see `adtree_reconstruct_compare.py` line 977: `geom_filename=np.array(out)`). Concretely: `write_geom()` writes to whatever bare filename string was captured, e.g. `geom_B21_S04_r5mm_adqsm090_seg35-500-k50.txt`, landing in the **current working directory at the time `export_geom_ansys.py` is run** (i.e. `scripts/python\` if launched conventionally) — **not** into `plots/<tree>/`.

### `import_matlab_results.py`
| Line | Expression | Resolves to |
|---|---|---|
| 189-190 | `open(csv_path,"w",...)` + `csv.DictWriter`, inside `upsert_result()`, called at line 351 with `RESULTS_CSV` | `volume_results.csv` (`RESULTS_CSV`, line 43) |

### `parameter_sensitivity.py`
| Line | Expression | Resolves to |
|---|---|---|
| 273 | `os.makedirs(out_dir, exist_ok=True)` inside `_output_dir(tree)`, `out_dir=os.path.join(ensure_plots_dir(), tree)` | `plots\<tree>\` (dyn, `tree`∈`REFERENCE_TREES`) |
| 291 | `fig.savefig(out_path,...)` inside `_save(...)`; per-tree: `plots/<tree>/sensitivity_<family>_<branch_filter>_<metric>_<chart>.png`; pooled (`tree is None`): `plots/sensitivity_multitree_<family>_<branch_filter>_<metric>_<chart>.png` | (dyn, many combinations) |

### `plot_box.py`
| Line | Expression | Resolves to |
|---|---|---|
| 1172 | `fig.savefig(out_path,...)`, `out_path=os.path.join(ensure_plots_dir(), "%s_boxplot_%s.png"%(tree,branch_filter))`; current call site line 1184 uses `SELECT_TREE="B21_S01"` (line 105), `BRANCH_FILTER="none"` (line 109) | `plots\B21_S01_boxplot_none.png` — **note: this file lands directly in `plots/`, not in a per-tree subfolder, despite being per-tree data** |

### `plot_param_sweep.py`
| Line | Expression | Resolves to |
|---|---|---|
| 450 | `ensure_plots_dir()` (imported from `plot_volumes.py`) | ensures `plots\` exists |
| 452 | `fig.savefig(out_path,...)`, `out_path=os.path.join(PLOTS_DIR, "param_sweep_%s_%s.png"%(SELECT_TREE,SWEEP_PARAM))` | `plots\param_sweep_B21_S01_pd2min.png` (current `SELECT_TREE`/`SWEEP_PARAM`) — **same loose-in-plots/ pattern as plot_box.py** |

### `plot_style.py`
No write operations — pure style constants and a color-sampling helper.

### `plot_volumes.py`
| Line | Expression | Resolves to |
|---|---|---|
| 279 | `os.makedirs(PLOTS_DIR)` (inside `ensure_plots_dir()`) | `plots\` (`PLOTS_DIR="plots"`, line 109) |
| 287 | `fig.savefig(out_path, dpi=PLOT_DPI, ...)` inside `save_and_report(fig, filename)`; call sites: line 883 `"total_volume_by_tree.png"`; line 1238 `"tree_overview_%s_%s%s.png"%(tree,branch_filter,PNG_FILENAME_SUFFIX)`; line 1400 `"error_boxplot_%s%s.png"%(branch_filter,suffix)`; line 1483 `"error_metrics_bar_%s%s.png"%(branch_filter,suffix)` | `plots\total_volume_by_tree.png`, `plots\tree_overview_<tree>_<branch_filter>.png`, `plots\error_boxplot_<branch_filter>[_<field>].png`, `plots\error_metrics_bar_<branch_filter>[_<field>].png` — **all loose directly in `plots/`, including the per-tree `tree_overview_*` file** |

### `qsm_volume_mean.py`
| Line | Expression | Resolves to |
|---|---|---|
| 296 | `open(csv_path,"w",...)` + `csv.DictWriter` inside `upsert_result()`, called lines 413/435 with `RESULTS_CSV` | `volume_results.csv` (`RESULTS_CSV`, line 83) |
| 448-461 | `open(CSV_PATH,"w",...)` + `csv.writer` + multiple `writerow()` | `<DATA_ROOT>\IND07_083\qsm_tanago\IND07_083_volume_summary.csv` (abs; `CSV_PATH` line 64, `DATA_DIR` line 48, `DATA_ROOT` line 45) |

### `reconstruction_method_decision_summary.py`
| Line | Expression | Resolves to |
|---|---|---|
| 132 | `os.makedirs(PLOTS_DIR)` (locally re-defined `ensure_plots_dir()`, not imported) | `plots\` (`PLOTS_DIR="plots"`, line 82) |
| 433-437 | `open(OUTPUT_CSV,"w",...)` + `csv.DictWriter` | `reconstruction_method_decision_summary.csv` (line 80) |
| 447-455 | `open(BEST_PICKS_CSV,"w",...)` + `csv.DictWriter` | `reconstruction_method_best_picks.csv` (line 81) |
| 490 | `fig.savefig(out_path,...)`, `out_path=os.path.join(ensure_plots_dir(), "reconstruction_decision_%s_%s.png"%(block["name"],metric))` | `plots\reconstruction_decision_<block>_<metric>.png` (12 possible files: 3 blocks × 4 metrics) |
| 637 | same pattern, `"...__bestof_pooled_%s.png"` | `plots\reconstruction_decision_<block>__bestof_pooled_<metric>.png` |
| 760 | same pattern, `"...__bestof_pertree_%s.png"` | `plots\reconstruction_decision_<block>__bestof_pertree_<metric>.png` |

### `reference_volume.py`
| Line | Expression | Resolves to |
|---|---|---|
| 210-211 | `open(csv_path,"w",...)` + `csv.DictWriter` inside `upsert_result()`, called line 321 with `RESULTS_CSV` | `volume_results.csv` (`RESULTS_CSV`, line 51) |

### `scratch_calib_sensitivity.py`
No writes — confirmed by inspection, matches its header's explicit claim.

### `scratch_order_alignment.py`
| Line | Expression | Resolves to |
|---|---|---|
| 130-131 | `open(path,"a",...)` + `csv.writer`, `path=RESULTS_SUMMARY_PATH` (lines 117-119) | An absolute path in the **Claude session's own scratchpad temp directory** (`...\AppData\Local\Temp\claude\...\scratchpad\order_alignment_summary.csv`) — outside the repo entirely. **The file's own header claims "Writes NOTHING to disk" — this is contradicted by the actual code.** |

### `scratch_skeleton_diag.py`
No writes — confirmed by inspection, matches its header's explicit claim.

### `taper_curve_compare.py`
| Line | Expression | Resolves to |
|---|---|---|
| 227 | `os.makedirs(out_dir, exist_ok=True)`, `out_dir=os.path.join("plots", tree_name)` (bare literal, inlined, not a constant) | `plots\B21_S04\` (dyn per `tree_name`∈`TREES_TO_RUN`) |
| 229 | `fig.savefig(out_path,...)`, `out_path=os.path.join(out_dir, "taper_curve_compare_%s.png"%tree_name)` | `plots\B21_S04\taper_curve_compare_B21_S04.png` |
| 286-289 | `open(SUMMARY_CSV_PATH,"w",...)` + `csv.DictWriter`, guarded by `if all_rows:` | `taper_curve_compare_summary.csv` (`SUMMARY_CSV_PATH`, line 80) |

### `tree_coverage_check.py`
| Line | Expression | Resolves to |
|---|---|---|
| 289-292 | `open(out_path,"w",...)` + `csv.DictWriter` inside `write_coverage_csv()`, called line 311 with `out_name` computed line 310 from `SELECT_TREE` | `coverage_all_trees.csv` (default, `SELECT_TREE="ALL"`) or `coverage_<tree>.csv` — **the file's own header claims "This script never writes anything" — contradicted by the actual code.** |

### `tree_geom_utils.py` (shared library — imported by other scripts, no RUN section of its own)
| Line | Function | Expression | Resolves to |
|---|---|---|---|
| 941 | `plot_radius_regression_per_order(..., plots_dir="plots")` | `os.makedirs(plots_dir)` | `plots` (default, caller-overridable) |
| 986 | same | `fig.savefig(plot_path,...)`, built lines 984-985 from `plots_dir`/`tree_name`/`variant_tag`/`filename_suffix` | (dyn) `plots/adtree_adqsm_radius_regression_perorder_<tree>_AdQSM<variant><suffix>.png` |
| 1184-1187 | `upsert_result(csv_path, ...)` | `open(csv_path,"w",...)` + `csv.DictWriter` | caller-supplied `csv_path` — in every observed call site this is `RESULTS_CSV`="volume_results.csv" |
| 1366 | `write_geom(path, ...)` | `open(path,"w")` | caller-supplied `path` — bare filename with no directory (see `export_geom_ansys.py` above) |
| 1500 | `plot_model(..., out_path, ..., save_png)` | `fig.savefig(png_path,...)`, `png_path=os.path.splitext(out_path)[0]+".png"` | derived from caller-supplied `out_path` |
| 1618 | `write_geom_from_treegraph(path, ...)` | `open(path,"w")` | caller-supplied `path` — in `treegraph_export_geom.py`'s call this resolves to an **absolute** path next to the TreeGraph JSON input, outside `scripts/python` entirely |

### `treegraph_add_to_comparison.py`
| Line | Expression | Resolves to |
|---|---|---|
| 75-79 | `upsert_result(RESULTS_CSV, ...)` (write itself inside `tree_geom_utils.upsert_result`) | `volume_results.csv` (`RESULTS_CSV`, line 31) |

### `treegraph_export_geom.py`
| Line | Expression | Resolves to |
|---|---|---|
| 39 | `write_geom_from_treegraph(GEOM_PATH, ...)`, `GEOM_PATH` computed lines 35-37 as `os.path.join(os.path.dirname(JSON_PATH), "geom_%s-cs%s.txt"%(name,cs))` | **Absolute** path: same directory as the (hardcoded, absolute) `JSON_PATH` input, i.e. under `C:\...\01_TreeGraph\treegraph_work\results\test_smoke\` — entirely outside `tree_reconstruction` |

### `visualize_pointcloud_region.py`
| Line | Expression | Resolves to |
|---|---|---|
| 79 | `fig.savefig(OUTPUT_DIR + "/1_cely_strom_zboku.png", ...)` | `1_cely_strom_zboku.png` (`OUTPUT_DIR="."`, line 42) |
| 103 | `fig.savefig(OUTPUT_DIR + "/2_pohled_shora.png", ...)` | `2_pohled_shora.png` |
| 128 | `fig.savefig(OUTPUT_DIR + "/3_detail_zboku.png", ...)` | `3_detail_zboku.png` — all three land directly in `scripts/python\` (or wherever CWD is), not under `plots/` at all |

**Total distinct write-operation call sites found across the 25 root `.py` files: 40** (counting each `open(...,"w"/"a")`/`savefig`/`csv.DictWriter`/`makedirs` call site once; several of these are the *same* logical write reached from multiple call sites, e.g. `upsert_result(RESULTS_CSV,...)` is called from 6 different files but the actual `open()` lives once, in `tree_geom_utils.py`).

---

## STEP 3 — Path constants

| Constant | File | Line | Value | Abs/Rel |
|---|---|---|---|---|
| `RESULTS_CSV` | `compare_volumes.py` | 71 | `"volume_results.csv"` | Relative |
| `RESULTS_CSV` | `adtree_reconstruct_compare.py` | 323 | `"volume_results.csv"` | Relative |
| `RESULTS_CSV` | `adqsm_variant_sensitivity.py` | 93 | `"volume_results.csv"` | Relative |
| `RESULTS_CSV` | `import_matlab_results.py` | 43 | `"volume_results.csv"` | Relative |
| `RESULTS_CSV` | `qsm_volume_mean.py` | 83 | `"volume_results.csv"` | Relative |
| `RESULTS_CSV` | `reference_volume.py` | 51 | `"volume_results.csv"` | Relative |
| `RESULTS_CSV` | `treegraph_add_to_comparison.py` | 31 | `"volume_results.csv"` | Relative |
| `NPZ_DIR` | `adtree_reconstruct_compare.py` | 332 | `"npz"` | Relative |
| `NPZ_DIR` | `export_geom_ansys.py` | 36 | `"npz"` | Relative |
| `NPZ_FILE` | `export_geom_ansys.py` | 31 | `"calib_B21_S01_r15mm_adqsm999_seg35-500-k50.npz"` | Relative |
| `FIGURES_DIR` | `adtree_reconstruct_compare.py` | 333 | `os.path.join("plots", TREE_NAME)` | Relative (derived) |
| `PLOTS_DIR` | `plot_volumes.py` | 109 | `"plots"` | Relative |
| `PLOTS_DIR` | `calmethod_decision_summary.py` | 63 | `"plots"` | Relative |
| `PLOTS_DIR` | `check_adqsm_thin_branch_length.py` | 72 | `"plots"` | Relative |
| `PLOTS_DIR` | `reconstruction_method_decision_summary.py` | 82 | `"plots"` | Relative |
| `OUTPUT_CSV` | `calmethod_decision_summary.py` | 62 | `"calmethod_decision_summary.csv"` | Relative |
| `OUTPUT_CSV` | `reconstruction_method_decision_summary.py` | 80 | `"reconstruction_method_decision_summary.csv"` | Relative |
| `BEST_PICKS_CSV` | `reconstruction_method_decision_summary.py` | 81 | `"reconstruction_method_best_picks.csv"` | Relative |
| `SUMMARY_CSV_PATH` | `taper_curve_compare.py` | 80 | `"taper_curve_compare_summary.csv"` | Relative |
| `TAPER_SUMMARY_CSV` | `adqsm_variant_sensitivity.py` | 94 | `"taper_curve_compare_summary.csv"` | Relative (same file as above, separately named constant) |
| `CSV_PATH` | `qsm_volume_mean.py` | 64 | `os.path.join(DATA_DIR, "%s_volume_summary.csv"%TREE_NAME)` | **Absolute** (derived from absolute `DATA_ROOT`) |
| `DATA_ROOT` | `adqsm_build_median_variant.py` | 64 | `r"C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\data"` | Absolute |
| `DATA_ROOT` | `adtree_reconstruct_compare.py` | 85 | same literal | Absolute |
| `DATA_ROOT` | `check_adqsm_thin_branch_length.py` | 50 | same literal | Absolute |
| `DATA_ROOT` | `qsm_volume_mean.py` | 45 | same literal | Absolute |
| `DATA_ROOT` | `reference_volume.py` | 38 | same literal | Absolute |
| `DATA_ROOT` | `scratch_calib_sensitivity.py` | 65 | same literal | Absolute |
| `DATA_ROOT` | `scratch_order_alignment.py` | 86 | same literal | Absolute |
| `DATA_ROOT` | `scratch_skeleton_diag.py` | 42 | same literal | Absolute |
| `DATA_ROOT` | `taper_curve_compare.py` | 63 | same literal | Absolute |
| `MATLAB_RESULTS_DIR` | `import_matlab_results.py` | 30 | `r"...\tree_reconstruction\scripts\matlab"` | Absolute |
| `JSON_PATH` | `treegraph_add_to_comparison.py` | 28 | `r"C:\...\01_TreeGraph\treegraph_work\results\test_smoke\IND01_054-cs0.15-tipNone.json"` | Absolute (outside repo) |
| `JSON_PATH` | `treegraph_export_geom.py` | 21 | same literal | Absolute (outside repo) |
| `GEOM_PATH` | `treegraph_export_geom.py` | 25 | `None` initially, computed absolute (see Step 2) | Absolute |
| `OUTPUT_DIR` | `visualize_pointcloud_region.py` | 42 | `"."` | Relative |
| `RESULTS_SUMMARY_PATH` | `scratch_order_alignment.py` | 117-119 | absolute path into the Claude session's own scratchpad temp dir | Absolute, outside repo |

`BEST_PICKS_CSV` and `SUMMARY_CSV_PATH`/`TAPER_SUMMARY_CSV` are all bare relative filenames, i.e. no directory component; they resolve to `scripts/python/` itself when the script's CWD assumption holds.

---

## STEP 4 — Path helper functions

### Directory-ensuring helpers (`ensure_plots_dir()` pattern)

| Defined in | Line | Behavior |
|---|---|---|
| `plot_volumes.py` | 276-280 | `if not os.path.isdir(PLOTS_DIR): os.makedirs(PLOTS_DIR)`; returns `PLOTS_DIR` ("plots"). **This is the canonical/original definition.** |
| `calmethod_decision_summary.py` | 164-167 | Identical logic, its own local `PLOTS_DIR="plots"` — **independently re-defined, not imported.** |
| `reconstruction_method_decision_summary.py` | 130-133 | Identical logic, its own local `PLOTS_DIR="plots"` — **independently re-defined, not imported.** |

**Importers of `plot_volumes.ensure_plots_dir`:**
- `plot_param_sweep.py` line 35: `from plot_volumes import ensure_plots_dir, PLOTS_DIR, treeqsm_pd_token, _pd_field_token`; called line 450.
- `plot_box.py` line 65-68: `from plot_volumes import (... ensure_plots_dir ...)`; called line 1171.
- `adqsm_variant_sensitivity.py` line 65: `from plot_volumes import ensure_plots_dir`; called line 325.
- `parameter_sensitivity.py` line 43: `from plot_volumes import classify_family, ensure_plots_dir, shorten_method_label`; called lines 272, 290.

So there are **three separate definitions** of the same `ensure_plots_dir()` logic (in `plot_volumes.py`, `calmethod_decision_summary.py`, `reconstruction_method_decision_summary.py`) and **four importers** of the `plot_volumes.py` original. All five resolve to the same literal `"plots"` today, but the two independent re-definitions would not automatically follow a change made only in `plot_volumes.py`.

### Other helper functions building/ensuring an output path

| Function | Defined in | Line | Builds |
|---|---|---|---|
| `_variant_dir(tree_name, variant)` | `adqsm_build_median_variant.py` | 68-69 | `os.path.join(DATA_ROOT, tree_name, variant)` — absolute; used both for reading source variants and (via module constants) for the one output folder it creates |
| `_output_dir(tree)` | `parameter_sensitivity.py` | 268-274 | `os.path.join(ensure_plots_dir(), tree)`, then `os.makedirs(out_dir, exist_ok=True)` |
| `_save(fig, family, branch_filter, metric, chart, tree=None)` | `parameter_sensitivity.py` | 277-293 | Assembles the final PNG path (via `_output_dir()` when per-tree) and performs the `fig.savefig()` |
| `save_and_report(fig, filename)` | `plot_volumes.py` | 283-289 | `os.path.join(ensure_plots_dir(), filename)`, then saves — the shared save routine reused throughout that file |
| `plot_radius_regression_per_order(..., plots_dir="plots")` | `tree_geom_utils.py` | 899-988 | Builds `plot_path` inline and calls `os.makedirs(plots_dir)` (line 941) directly in the function body — not delegated to a separate helper |
| `write_coverage_csv(csv_rows, out_path)` | `tree_coverage_check.py` | 283-293 | Takes `out_path` as a plain parameter; does not build or ensure a directory itself |

### Scripts that build output paths INLINE (no shared helper call) instead of using any of the above

- `adqsm_build_median_variant.py` — inline `os.path.join(DATA_ROOT, ...)` via its own private `_variant_dir()` (not shared with any other file).
- `adtree_reconstruct_compare.py` — inline `os.path.join(NPZ_DIR, ...)`, `os.path.join(FIGURES_DIR, out)`, `os.makedirs(NPZ_DIR/FIGURES_DIR, exist_ok=True)` directly at module level (not organized into any function at all — this is a flat top-to-bottom script).
- `export_geom_ansys.py` — inline `os.path.join(NPZ_DIR, NPZ_FILE)`; no makedirs call anywhere (assumes `npz/` already exists).
- `qsm_volume_mean.py` — inline `os.path.join(DATA_DIR, ...)`; no makedirs call anywhere.
- `reference_volume.py`, `import_matlab_results.py`, `treegraph_add_to_comparison.py` — no local directory-building logic; delegate the CSV write entirely to `upsert_result()` (which itself takes a fully-formed path and never calls makedirs).
- `taper_curve_compare.py` — inline `os.path.join("plots", tree_name)` and `os.makedirs(out_dir, exist_ok=True)` directly inside `plot_taper_variants()`, using the bare literal `"plots"` rather than any `PLOTS_DIR` constant or `ensure_plots_dir()` call.
- `treegraph_export_geom.py` — inline `os.path.join(os.path.dirname(JSON_PATH), ...)` at module level.
- `visualize_pointcloud_region.py` — inline string concatenation `OUTPUT_DIR + "/filename.png"`.

**Summary:** directory-creation for `plots/` is centralized around `ensure_plots_dir()`, but that function itself exists in **three independent copies**. Everything else (the `data/` tree, `npz/`, per-script CSVs, the TreeGraph absolute paths) is built ad hoc, inline, per-script, with no shared helper at all.

---

## STEP 5 — Read/write dependency graph

| Produced file / pattern | Written by (file:line) | Read by (file:line, or "none found") |
|---|---|---|
| `volume_results.csv` | `compare_volumes.py:433` (starter only); `import_matlab_results.py:189` (via `upsert_result`); `qsm_volume_mean.py:296`; `reference_volume.py:210`; `treegraph_add_to_comparison.py:75`; `adtree_reconstruct_compare.py:513,1005,1049,1081,1094,1114,1127`; underlying `open()` in all cases lives in `tree_geom_utils.py:1184` (`upsert_result`) | `compare_volumes.py:438` (`load_results`); `plot_volumes.py:1494`; `plot_box.py:1183`; `plot_param_sweep.py:270`; `parameter_sensitivity.py:1089`; `calmethod_decision_summary.py:378`; `reconstruction_method_decision_summary.py:823`; `tree_coverage_check.py:301`; `adqsm_variant_sensitivity.py:224` |
| `taper_curve_compare_summary.csv` | `taper_curve_compare.py:286` (`SUMMARY_CSV_PATH`) | `adqsm_variant_sensitivity.py:232` (`read_taper_summary(TAPER_SUMMARY_CSV,...)`, its own separately-named constant pointing at the same literal filename) |
| `calib_*.npz` (in `npz/`) | `adtree_reconstruct_compare.py:971` | `export_geom_ansys.py:111` (`load_calibrated_geometry`, via `list_available_npz`/manual `NPZ_FILE` selection) |
| `geom_*.txt` (actual ANSYS text geometry) | `export_geom_ansys.py:114` (`write_geom`, bare filename, lands in CWD); `treegraph_export_geom.py:39` (`write_geom_from_treegraph`, absolute path outside repo) | **None found within `scripts/python`.** Only consumed by `.mac`/`.m` files outside Python (see Step 6) — and even there, the reference is indirect (see Step 6 detail: no script reads the Python-produced file by that exact resolved path). |
| `geom_*_r*mm*.png` (diagnostic radius-regression plots, filename shares the `geom_` prefix with the text file above but is a different artifact) | `adtree_reconstruct_compare.py:1160` → `tree_geom_utils.py:986` (`plot_radius_regression_per_order`) and `tree_geom_utils.py:1500` (`plot_model`, via `os.path.splitext(out_path)[0]+".png"`) | none found (visual/diagnostic only) |
| `calmethod_decision_summary.csv` | `calmethod_decision_summary.py:222` | none found |
| `reconstruction_method_decision_summary.csv` | `reconstruction_method_decision_summary.py:433` | none found |
| `reconstruction_method_best_picks.csv` | `reconstruction_method_decision_summary.py:447` | none found |
| `coverage_all_trees.csv` / `coverage_<tree>.csv` | `tree_coverage_check.py:289` | none found |
| `<tree>_volume_summary.csv` (absolute, under `data/<tree>/qsm_tanago/`) | `qsm_volume_mean.py:448` | none found |
| `plots/*.png` general (tree-overview, error-boxplot, error-metrics-bar, total-volume-by-tree, sensitivity charts, param-sweep, boxplots, calmethod/reconstruction-decision charts, taper-curve/variant-sensitivity charts) | `plot_volumes.py`, `parameter_sensitivity.py`, `plot_box.py`, `plot_param_sweep.py`, `calmethod_decision_summary.py`, `reconstruction_method_decision_summary.py`, `taper_curve_compare.py`, `adqsm_variant_sensitivity.py`, `tree_geom_utils.py` (see Step 2 for exact lines) | none found — every PNG is a terminal visual artifact, never re-read by another Python script |
| `order_alignment_summary.csv` (scratchpad, outside repo) | `scratch_order_alignment.py:130` | none found — not part of the repo/pipeline at all |

**Flags per the task instructions:**
- **Written by MORE than one script:** `volume_results.csv` (7 writer call sites across 6 files, funneled through one shared `upsert_result()`), `taper_curve_compare_summary.csv`/`TAPER_SUMMARY_CSV` is written by one script but *named* by two constants in two files (not a multi-writer case, but a naming-duplication risk), `ensure_plots_dir()`/`PLOTS_DIR="plots"` is independently defined identically in 3 files (not a data file, but the same *directory* is created redundantly by non-shared code).
- **Read by a script OTHER than the one that wrote it:** `volume_results.csv` (written by 6 scripts, read by 9 different scripts — the single biggest coordination surface in the project) and `taper_curve_compare_summary.csv` (written by `taper_curve_compare.py`, read by `adqsm_variant_sensitivity.py`).
- **Produced-but-never-read files:** the vast majority of PNG/CSV outputs (`calmethod_decision_summary.csv`, `reconstruction_method_decision_summary.csv`, `reconstruction_method_best_picks.csv`, `coverage_*.csv`, all PNGs) are terminal — nothing downstream in the Python codebase reads them back.

---

## STEP 6 — External (non-Python) consumers

Full sweep of all 87 `.m` and 6 `.mac` files for `geom_*.txt`, `volume_results.csv`, `plots/`, `npz/`.

**No `.m`/`.mac` file references `plots/` or `npz/` as a directory path anywhere** (the word "plots" appears only as plain English verb usage, e.g. "plots the QSM").

**`volume_results.csv`** appears in exactly one place, as a comment only:
- `scripts/matlab/runsken.m:1378` — descriptive comment noting a taper value should be "directly comparable to theirs in volume_results.csv." Not an actual read/write reference.

**`geom_*.txt`:**

| File | Line | Reference | Type |
|---|---|---|---|
| `scripts/ansys/bk1_qsm.mac` | 27 | `*tread,geom,geom_1,txt` → `geom_1.txt` | Hardcoded, relative (ANSYS reads from its own CWD at execution time, not set within the `.mac` file) |
| `scripts/ansys/bk1.mac` | 17 | `*tread,geom,geom_IND01_054_v1auto_1,txt` → `geom_IND01_054_v1auto_1.txt` | Hardcoded, relative |
| `scripts/matlab/runsken.m` | 76, 284 | `'geom_*.txt'` glob pattern, used via `fullfile(search_dir, pattern)` where `search_dir` derives from `clean_search_dirs = {data_dir, fullfile(data_dir,'results')}` | Built from a variable; `data_dir` itself is hardcoded absolute at line 16 (`'C:\...\tree_reconstruction\scripts\matlab'`) — **points at the MATLAB folder, not `scripts/python`** |
| `scripts/matlab/runsken.m` | 590 | `export_prefix = ['geom_' tree_id '_' run_tag '_'];` | Built from variables `tree_id`/`run_tag`; written relative to `data_dir` (hardcoded absolute, line 16) after `cd(data_dir)` at line 594 |
| `scripts/matlab/runsken.m` | 1636, 1642-1643 | `ansys_export_prefix = ['geom_' tree_id '_' ansys_tag '_'];` then `file_name = sprintf('%s%d.txt', ansys_export_prefix, i);` then `writematrix(geom_table, file_name, ...)` | Built from variables; same `data_dir` base |

**Conclusion:** MATLAB's own `geom_*.txt` export/cleanup logic in `runsken.m` operates entirely within `scripts/matlab/` (and its `results` subfolder) — a **separate namespace** from the Python pipeline's `geom_*.txt` (written by `export_geom_ansys.py`/`treegraph_export_geom.py` to a bare CWD-relative filename, or an absolute TreeGraph-adjacent path, respectively — see Step 2/5). The two ANSYS `.mac` files that actually `*tread` a `geom_*.txt` file (`bk1.mac`, `bk1_qsm.mac`) use hardcoded relative filenames with no directory logic at all, meaning whoever runs ANSYS is responsible for having placed the correct `geom_*.txt` file in ANSYS's working directory by hand — there is no automated or scripted hand-off from any Python or MATLAB script into the exact filename/location these two `.mac` files expect.

---

## STEP 7 — Actual filesystem state

### 7a. Non-`.py` files directly in `scripts/python/`

| File | Size | Modified |
|---|---|---|
| `coverage_B21_S04.csv` | 1,337 B | 2026-09-08 11:39:19 |
| `coverage_all_trees.csv` | 2,610 B | 2026-09-08 11:40:57 |
| `qsm_run.log` | 3,757 B | 2026-08-24 15:55:22 (notably older than everything else) |
| `taper_curve_compare_summary.csv` | 714 B | 2026-09-08 10:43:30 |
| `volume_results.csv` | 208,282 B | 2026-09-08 11:35:28 |

### 7b. `scripts/python/plots/` — two levels deep

Two per-tree subfolders exist, both flat (no further subfolders):

**`plots/B21_S01/` — 216 files**, all prefixed `geom_B21_S01_...` (PNG twins of the radius-regression diagnostic, not the ANSYS text geometry — see Step 5/6). Breakdown: 4 radius variants (r5mm/r10mm/r15mm/r20mm, 54 each) × 3 AdQSM versions (adqsm080/090/999, 72 each) × 6 seg/k combos (`seg{35,50,65}-500-k{50,80,100}`) × 2 (`-wlen` suffix present/absent) = 216.

**`plots/B21_S04/` — 145 files**: 144 files prefixed `geom_B21_S04_...` (4 radius variants × 2 AdQSM versions [adqsm090/999 only — no adqsm080 for this tree] × 6 seg/k combos × 2 wlen-suffix = 144), plus 1 file `taper_curve_compare_B21_S04.png` (182,223 B, 2026-09-08 10:43).

### 7c. Files sitting DIRECTLY in `plots/` (not in any per-tree subfolder)

**90 files**, all prefixed `adtree_adqsm_radius_regression_perorder_` — **these are the files the user suspects are misplaced.** Split: 54 for `B21_S01` (AdQSM080/090/999 × 6 seg-k combos × 2 wlen variants), 36 for `B21_S04` (AdQSM090/999 × 6 seg-k combos × 2 wlen variants). Total combined size ≈ 14,914,236 B (≈14.2 MiB); individual files ≈152,000–177,000 B; all modified 2026-09-08 between 09:52 and 11:06. This is exactly the same logical output type as the "geom_*" PNGs inside the per-tree subfolders (both produced by `tree_geom_utils.plot_radius_regression_per_order()`), but landing loose at the `plots/` root instead of under a per-tree subfolder — i.e., **the same function's output is inconsistently placed depending on caller/parameters** (see Step 8).

### 7d. `scripts/python/npz/`

Flat, no subfolders. **360 files, total 4,807,226,368 bytes (≈4,584.5 MiB / ≈4.48 GiB).**

### 7e. Other output-holding folders under `scripts/python/`

No folders exist beyond `plots/`, `npz/`, `__pycache__/`, and the four archive/backup folders already inventoried in Step 1b (`archiv_bak_2026-09-06/`, `archive_plots_IND/`, `archive_vol_res_IND/`, `backup_BK/`). Nothing further to report.

---

## STEP 8 — Classification

| Pattern | Classification | Evidence |
|---|---|---|
| `plots/<TREE>/geom_<TREE>_r<mm>mm_adqsm<variant>_seg<...>[-wlen].png` | **PER-TREE** (already under `plots/<TREE>/`) | `<TREE>` is a literal path segment of the containing folder AND embedded in the filename itself (`tree_name` parameter of `plot_radius_regression_per_order`/`plot_model`) — unambiguous. |
| `plots/adtree_adqsm_radius_regression_perorder_<TREE>_AdQSM<variant>_seg<...>[-wlen].png` (loose in `plots/` root) | **PER-TREE**, currently misfiled at `plots/` root instead of `plots/<TREE>/` | Same producing function (`tree_geom_utils.plot_radius_regression_per_order`) as the row above; `<TREE>` is embedded in the filename (`B21_S01`/`B21_S04`) even though the containing folder is not per-tree — tree identity is present in the name, just not reflected in the folder location. |
| `plots/<TREE>_boxplot_<branch_filter>.png` (from `plot_box.py`) | **PER-TREE**, loose in `plots/` root | `<TREE>` (`SELECT_TREE`) is the leading filename segment; not placed in a per-tree subfolder. |
| `plots/tree_overview_<TREE>_<branch_filter>.png` (from `plot_volumes.py`) | **PER-TREE**, loose in `plots/` root | `<TREE>` (loop variable from `TREES_TO_PLOT`) is embedded in the filename; not placed in a per-tree subfolder. |
| `plots/<TREE>/sensitivity_<family>_<branch_filter>_<metric>_<chart>.png` (from `parameter_sensitivity.py`, per-tree call) | **PER-TREE** (already under `plots/<TREE>/`) | `<TREE>` is both the containing folder (via `_output_dir(tree)`) and not otherwise in the filename — folder is the sole tree-identifying element here. |
| `plots/sensitivity_multitree_<family>_<branch_filter>_<metric>_<chart>.png` | **ALL-TREES** | Explicitly the `tree is None` branch of `_save()`; filename literally contains `multitree`; no single tree's data. |
| `plots/param_sweep_<TREE>_<param>.png` | **PER-TREE**, loose in `plots/` root | `<TREE>` (`SELECT_TREE`) is the second filename segment. |
| `plots/<TREE>/adqsm_variant_sensitivity_<TREE>.{png,csv}` | **PER-TREE** (already under `plots/<TREE>/`) | `<TREE>` is both the folder and repeated in the filename. |
| `plots/<TREE>/taper_curve_compare_<TREE>.png` | **PER-TREE** (already under `plots/<TREE>/`) | Same pattern as above. |
| `plots/total_volume_by_tree.png` | **ALL-TREES** | Filename explicitly says "by_tree" (aggregates every tree into one chart); no single-tree scoping anywhere in the code path (`save_and_report(fig, "total_volume_by_tree.png")`, no tree parameter). |
| `plots/error_boxplot_<branch_filter>[_<field>].png`, `plots/error_metrics_bar_<branch_filter>[_<field>].png` | **ALL-TREES** | Built from `rows` spanning `load_results(RESULTS_CSV)` filtered only by `branch_filter`/`field`, never by a single tree — these are project-wide error-metric summaries across every tree in `volume_results.csv`. |
| `plots/calmethod_decision_<field>.png` | **ALL-TREES** | `field`∈{total,trunk,branch} only; no tree dimension in the function signature or filename — `plot_calmethod_decision_chart(summary_rows, field)` pools across whatever trees are present in `RESULTS_CSV`. |
| `plots/reconstruction_decision_<block>_<metric>.png` and its `__bestof_pooled_`/`__bestof_pertree_` variants | **ALL-TREES**, EXCEPT the `__bestof_pertree_` variant is ambiguous | The base and `__bestof_pooled_` variants pool across trees (no tree in filename/signature). The `__bestof_pertree_` variant's *content* likely varies per tree internally (per its name), but the **output filename itself still has no tree segment** — it is one chart covering all trees' per-tree best picks together, so classified ALL-TREES by the filename-based rule used throughout this audit, with the caveat noted. |
| `calmethod_decision_summary.csv`, `reconstruction_method_decision_summary.csv`, `reconstruction_method_best_picks.csv`, `coverage_all_trees.csv` | **ALL-TREES** (or **SHARED-STATE**-adjacent) | All are decision/summary tables spanning multiple trees' rows (from `RESULTS_CSV`), living at `scripts/python/` root, not under `plots/` at all — flagged as ambiguous between ALL-TREES-summary and SHARED-STATE, since unlike `volume_results.csv` these are themselves *derived/output* artifacts (not the master input table), yet like `volume_results.csv` they are project-wide CSVs sitting outside `plots/`. Not proposing a category resolution — flagging the ambiguity as instructed. |
| `coverage_<TREE>.csv` | **PER-TREE** | `<TREE>` (`SELECT_TREE`, non-"ALL") is the filename segment; only produced when a single tree is selected. |
| `<TREE>_volume_summary.csv` (absolute, under `data/<TREE>/qsm_tanago/`) | **PER-TREE**, but **not under `scripts/python` at all** | `<TREE>` is both the containing folder segment (`data/<TREE>/`) and the filename prefix; however this file lives entirely outside `scripts/python`'s output trees (`plots/`, `npz/`) — flagged as out-of-scope for any `plots/`-rooted reorganisation. |
| `geom_*.txt` (actual ANSYS geometry text, from `export_geom_ansys.py`/`treegraph_export_geom.py`) | **GEOMETRY**, but currently NOT under any `plots/<TREE>/geom/`-style location | Filename contains "geom" and `<TREE>` is embedded in the AdTree-pipeline variant (`geom_%s_r%dmm...`); however the AdTree-pipeline version resolves to a **bare CWD-relative filename** (no directory at all) and the TreeGraph-pipeline version resolves to an **absolute path entirely outside the repo** — neither currently lands anywhere near `plots/`. |
| `plots/<TREE>/geom_*.png` and `plots/geom_*...perorder...png` (the diagnostic PNGs) | **GEOMETRY** (per the literal "geom" substring rule) — see caveat | Filename contains "geom_"; `<TREE>` present. Caveat: these PNGs are a **different artifact** from the actual `geom_*.txt` ANSYS text file above (a radius-regression diagnostic chart, not tree geometry itself) — flagged as a naming ambiguity: the classification rule "filename contains 'geom'" catches both a true geometry export and an unrelated diagnostic plot that happens to share the naming convention. |
| `volume_results.csv` | **SHARED-STATE** | Per the task's explicit instruction — the project-wide master results table, written by 6 different scripts and read by 9 different scripts (Step 5); not a plot output. Classified as SHARED-STATE only, no move proposed. |
| `taper_curve_compare_summary.csv` | **SHARED-STATE**-like (ambiguous) | Written by one script, read by another (`adqsm_variant_sensitivity.py`) — functionally a small shared intermediate table, similar in role to `volume_results.csv` but scoped to one diagnostic. Not confidently ALL-TREES nor a pure per-run artifact; flagged as ambiguous rather than guessed. |
| `calib_*.npz` (in `npz/`) | **CACHE** | Explicit per the task's definition; filename embeds `<TREE>` (`calib_%s_r%dmm...`% TREE_NAME) plus radius/variant/seg parameters — these are intermediate calibrated-geometry caches consumed only by `export_geom_ansys.py`, never a terminal deliverable. |
| `order_alignment_summary.csv` (scratchpad, outside repo) | **Out of scope** | Lives in the Claude session's own temp scratchpad directory, not in the project at all — not classified under any of the five categories since it isn't a project output. |

---

## Risk flags

In priority order, per the task's instructions:

1. **Outputs whose path change would require coordinated edits in more than one script.**
   - `volume_results.csv`: the constant `RESULTS_CSV = "volume_results.csv"` is independently re-declared as a **literal string** in 7 different files (`compare_volumes.py:71`, `adtree_reconstruct_compare.py:323`, `adqsm_variant_sensitivity.py:93`, `import_matlab_results.py:43`, `qsm_volume_mean.py:83`, `reference_volume.py:51`, `treegraph_add_to_comparison.py:31`) rather than imported from one place — 4 further files import `RESULTS_CSV` from `compare_volumes.py` (`plot_box.py`, `plot_param_sweep.py`, `parameter_sensitivity.py`, `tree_coverage_check.py`), and `plot_volumes.py` and `reconstruction_method_decision_summary.py` also import it from `compare_volumes.py`. A path change to this file requires touching at least the 7 files with independent literals, plus verifying every importer still resolves correctly.
   - `PLOTS_DIR`/`ensure_plots_dir()`: independently defined in 3 files (`plot_volumes.py`, `calmethod_decision_summary.py`, `reconstruction_method_decision_summary.py`) and imported by 4 more — a path change here requires editing 3 separate function bodies, not just one shared helper.
   - `taper_curve_compare_summary.csv`: named by two different constants in two different files (`SUMMARY_CSV_PATH` in `taper_curve_compare.py`, `TAPER_SUMMARY_CSV` in `adqsm_variant_sensitivity.py`) that must be kept in sync by hand — no shared import exists between them.

2. **Outputs consumed by `.m` or `.mac` files.**
   - `geom_*.txt` is the only such case (Step 6). The two ANSYS macros that actually read it (`bk1.mac`, `bk1_qsm.mac`) use hardcoded relative filenames with zero directory logic, meaning any path reorganisation of Python's `geom_*.txt` output must also account for a currently-entirely-manual hand-off into ANSYS's working directory — there is no script anywhere that automatically stages the file where these `.mac` files expect it. `runsken.m`'s own `geom_*.txt` glob/export logic is a parallel, separate namespace rooted at `scripts/matlab/` and does not consume the Python-produced files at all — worth distinguishing so a Python-side path change is not mistaken for affecting the MATLAB pipeline.

3. **Bare relative paths that would silently write to the wrong place if launched from a different working directory than `scripts/python`.**
   - Every `RESULTS_CSV="volume_results.csv"`, `PLOTS_DIR="plots"`/`NPZ_DIR="npz"`, `OUTPUT_CSV`/`BEST_PICKS_CSV`/`SUMMARY_CSV_PATH`/`coverage_*.csv` constant listed in Step 3 is a bare relative literal with this exposure.
   - Most acute case: `export_geom_ansys.py`'s `write_geom(geom_filename, ...)` writes to a **bare filename with no directory prefix at all** (not even `"plots"` or `"."`) — it lands directly in whatever the current process's CWD happens to be at the moment the script runs, with no fallback and no `os.makedirs` call anywhere in that file.
   - `visualize_pointcloud_region.py`'s `OUTPUT_DIR = "."` has the same exposure — its three PNGs go straight into the launch-time CWD.
   - `taper_curve_compare.py`'s `plots/<tree_name>` in `plot_taper_variants()` is built from the bare literal `"plots"` inline (not even its own module constant), independent of `PLOTS_DIR` elsewhere.

4. **Two different scripts writing the same filename.**
   - No case was found where two *different literal filenames* collide outright. The closest cases are (a) `volume_results.csv`, which is written by 6 different scripts but this is the file's documented, intentional shared-upsert design, not an accidental collision; and (b) the `geom_` filename-prefix ambiguity noted in Step 8 — the AdTree-pipeline's actual ANSYS-input `geom_*.txt` text file and the diagnostic radius-regression `geom_*...perorder...png` chart share the same "geom_" naming convention and similar tree/threshold/variant tokens in their names, despite being unrelated artifact types written by different code paths (`tree_geom_utils.write_geom()` vs. `tree_geom_utils.plot_radius_regression_per_order()`/`plot_model()`) — a filename-pattern collision risk even though no single concrete filename was found written twice.

---

*End of audit. No files other than this report were created; no existing file was modified, moved, deleted, or run to produce it.*
