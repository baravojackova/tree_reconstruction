# Tree Reconstruction & Volume Comparison Pipeline — Guide

> Working notes for the project `tree_reconstruction`. The project can sit
> anywhere on disk (Windows or Mac): every path is computed from the place
> where `scripts/python/paths.py` lives, see §2.
>
> This file is meant to be edited over time (by you, or by Claude Code on
> your instructions) as the pipeline evolves. Keep it in the project root
> so both you and any AI assistant working on the code can see the current
> state at a glance.

---

## 1. What this project is for

You are reconstructing the 3D branch/trunk structure of individual trees
from terrestrial laser-scanning (TLS) point clouds, using **several
independent methods**, so that:

1. **The volume/shape estimates can be cross-checked** against each other
   and against a **destructive reference measurement** (a tree that was
   physically cut apart and measured section by section).
2. **The best-suited reconstruction can be exported into ANSYS** (via a
   `.mac` macro, `buk.mac`) as a beam/cylinder model for structural
   analysis.

The methods currently in play:

| Method | What it produces | Where it runs | Status in this pipeline |
|---|---|---|---|
| **AdTree** | A connected centerline skeleton (`.ply`) with a radius per vertex | External tool (AdTree software), output read by Python | Fully scripted (Python) |
| **AdQSM** | A cylinder-based QSM: `taper.txt`, `BranchStructure.txt`, `TreesParams.txt`, plus a `.obj` mesh | External tool (AdQSM), its text exports are only *read* here | Used as calibration reference + its own result rows |
| **TreeQSM** | A cylinder QSM built from many random realizations | MATLAB (`runsken.m`) | Fully scripted (MATLAB) + imported into Python |
| **TreeGraph** | A skeleton graph with per-node radius estimates, aimed at *topology* rather than just volume | Python (your own pipeline, in active development) | Newest method; not yet wired into `compare_volumes.py`/`volume_results.csv` — see §7 |
| **Destructive reference** | Ground-truth section volumes from a felled tree | Field measurement, read from a table | Read-only ground truth |

Every method that produces a volume estimate writes (or is imported into)
**one shared table**, `volume_results.csv`, so that `compare_volumes.py`
and `plot_volumes.py` can compare all of them side by side without you
having to re-type any numbers.

---

## 2. Folder / naming conventions you need to know

- **`scripts/python/paths.py` is the single place that knows where things
  are.** It computes every location from where `paths.py` itself lives
  (`PROJECT_ROOT`, `DATA_DIR`, `PLOTS_DIR`, `CSV_DIR`, `NPZ_DIR`,
  `ANSYS_GEOM_DIR`, `VOLUME_RESULTS_CSV`, ...), so **no path has to be
  edited when the project moves** to another folder or computer, and the
  scripts can be started from any folder (VS Code, terminal, project root).
- **Paths that differ between computers** (scan folders, the TreeGraph
  results folder outside the project, the MATLAB data folder) are in
  `local_paths.json` in the project root. That file is **not** in Git (every
  computer has its own); the template `local_paths.example.json` is. On a
  new computer: copy the template to `local_paths.json` and fill in the
  paths (details in `README.md`). Python reads it with
  `paths.get_local_path("<key>")`, MATLAB (`runsken.m`) reads the key
  `matlab_data_dir`.
- **Tree data: `data/<series folder>/<tree ID>/`.** `data/` is split by
  *series* (one plot / scan year / state); each tree has its own folder:

  ```
  data/
    2021_bezlisty/             series B21  (leafless, 2021)
      B21_S01/  B21_S04/  ...
    2024_bezlisty_po_redukci/  series B24  (leafless, 2024, after thinning)
      B24_S01/  ...
    IND_tanago2017/            series IND  (González de Tánago et al. 2017)
      IND01_054/  IND03_088/  IND07_083/
  ```

  Tree ID = `<SERIES>_Sxx` (B = Babice plot, 21/24 = scan year, S01 = tree
  number). The series is the text before `_S<number>` (`B21_S10_NOBUTT` is
  still series `B21`). The table `paths.SERIES_DATA_DIRS` says which folder
  belongs to which series, and **scripts never build a tree path by hand**:
  they call `paths.tree_data_dir(TREE_NAME)` (returns
  `data/<series folder>/<TREE_NAME>`). A tree from a series that is not in
  the table stops with a readable error. Adding a series = one line in
  `SERIES_DATA_DIRS`.
- Inside a tree's folder you will typically find:
  - the AdTree skeleton `<TREE_NAME>_noplate_clean_skeleton.ply` (B21; the
    IND trees use e.g. `IND07_083 - Cloud_nobutt_skeleton.ply`).
  - one or more AdQSM run subfolders (e.g. `05`, `060`), each containing
    `taper.txt`, `BranchStructure.txt`, `TreesParams.txt`.
  - `IND.h.trees.txt` — the destructive field-measurement table (shared
    across trees, filtered by `treeID` column).
  - a `qsm_tanago` subfolder with 20 TreeQSM "de Tanago" realization
    files (published reference reconstructions, not yours).
- The Python scripts can be run **from any folder**: `volume_results.csv`,
  `plots/`, `csv/`, `npz/`, `ansys_geom/` and `data/` are all resolved by
  `paths.py` to absolute paths. (One exception is a data file that Git
  ignores: `Popis_Babice_VSE_09012025_ans_lev6_0.txt` in `scripts/python/`
  has to be copied to a new computer by hand.)
- Output files are tagged with the settings that produced them so that
  re-running a script with different parameters **never silently
  overwrites** a previous result — e.g.
  `calib_IND01_054_r5mm_seg100-500-k50.npz` encodes the tree, the radius
  pruning threshold (5 mm) and the resampling settings
  (`seg100-500-k50` = `SEG_LEN_MIN=100mm, SEG_LEN_MAX=500mm, SEG_LEN_K=0.50`).

### 2.1 Jak uložit výsledky AdQSM pro nový strom

AdQSM se spouští ručně mimo Python, skripty výsledky jen **čtou**. Aby je
našly, ulož je takto:

1. **Série musí být v `paths.SERIES_DATA_DIRS`** (`scripts/python/paths.py`).
   Pro novou sérii přidej jeden řádek, např.
   `"B23_po_ol": DATA_DIR / "2023_olisteny_po_redukci",`. Klíč je text
   před `_S<číslo>` v názvu stromu. Bez toho skript skončí chybou, která
   sérii jmenuje.
2. **Založ složku stromu** `data/<složka série>/<ID stromu>/`, např.
   `data/2024_bezlisty_po_redukci/B24_S07/`. ID stromu je vždy
   `<SÉRIE>_Sxx`. Číslo stromu se musí shodovat s ostatními sériemi
   (z `B21_S07` se pro SAG stane `BK07`).
3. **Jeden běh AdQSM = jedna podsložka varianty** uvnitř složky stromu,
   např. `data/2024_bezlisty_po_redukci/B24_S07/060/`. Název varianty je
   libovolný, ale ze samých číslic (skripty je řadí číselně a rozeznají
   podle toho, že obsahují `taper.txt`). U B21 jsou `040`…`100`, u IND
   `04`…`10`. Složka `999` je vyhrazená pro **syntetický medián** variant,
   který vytvoří `adqsm_build_median_variant.py` — ručně do ní nic neukládej.
4. **Do podsložky varianty zkopíruj výstup AdQSM** (celou jeho složku).
   Skripty čtou tyto tři soubory (názvy přesně, včetně velkých písmen):

   | Soubor | K čemu |
   |---|---|
   | `taper.txt` | průměr kmene po výšce (taper, kalibrace průměrů) |
   | `BranchStructure.txt` | tabulka větví (řád, průměr, délka, objem, …) |
   | `TreesParams.txt` | celkové parametry stromu (objem, délky, výška, …) |

   Další soubory z AdQSM (`BranchOrderInfo.txt`, `TreesParams.csv`) se
   zkopírují také, neškodí.
5. Výsledná struktura jednoho stromu:

   ```
   data/2024_bezlisty_po_redukci/B24_S07/
     060/  taper.txt  BranchStructure.txt  TreesParams.txt  ...
     070/  ...
   ```

   Volitelně vedle toho leží vstup (mračno) a skelet AdTree
   `B24_S07_noplate_clean_skeleton.ply`, jen pokud poběží AdTree skripty.
6. Jméno stromu pak zadáš jen jako `TREE_NAME = "B24_S07"` na začátku
   skriptu. Cestu ke složce dopočítá `paths.tree_data_dir()`. Ověření:
   `python -c "import paths; print(paths.tree_data_dir('B24_S07'))"` ve
   `scripts/python` vypíše očekávanou cestu.

---

## 3. The shared results table: `volume_results.csv`

Every comparison ultimately reads/writes this one CSV. Columns:

| Column | Meaning |
|---|---|
| `tree` | Tree ID, e.g. `IND01_054` |
| `method` | Free-text label, e.g. `"AdQSM (TreesParams)"`, `"AdTree calibrated r5mm"` |
| `total_m3`, `stem_m3`, `branch_m3` | Wood volume [m³], total / trunk / branches |
| `std_m3` | Standard deviation of `total_m3` across repeated realizations, if known |
| `dbh_m` | Stem diameter at 1.3 m above ground (breast height) |
| `height_m` | Total tree height |
| `taper_cm_per_m` | Stem narrowing between 1.3 m and 10.0 m, in cm per metre |
| `trunk_len_m`, `branch_len_m` | Centerline length of trunk / branches |
| `branch_filter` | `"none"` = full reconstruction; `"10cm"` = trunk/branches restricted to ≥10 cm diameter |
| `n_cylinders` | How many cylinders make up that row's model |

**Why `branch_filter` exists:** the destructive field crew only ever
measured down to a 10 cm taper diameter — they never cut and weighed the
finest twigs. Comparing a method's *full* reconstruction (which does
include those twigs) against the reference would unfairly penalize it for
wood the reference never tried to capture. So every reconstruction method
can optionally write **a second, `"10cm"`-filtered row**, and
`compare_volumes.py` keeps the two comparisons (`"vs. reference"` and
`"methods vs. each other"`) strictly separate.

A row is uniquely identified by `(tree, method)`. Every writer script
**overwrites** the existing row for that pair instead of duplicating it,
so re-running a script after a parameter tweak is always safe — as long
as you also change something in the method *name* (a threshold, a
resampling suffix, an AdQSM variant tag, …) if you want to keep the old
result *and* the new one side by side.

---

## 4. Method 1 — AdTree pipeline (Python, 2 steps)

### Step 1: `adtree_reconstruct_compare.py`

**Purpose:** turn the raw AdTree skeleton into a clean, pruned,
resampled cylinder model; optionally calibrate its radii against AdQSM;
compute/print/record volume, DBH, height and taper; save the result.

**Required input files** (all read automatically, built from `TREE_NAME`):
| File | From | Used for |
|---|---|---|
| `<TREE_NAME> - Cloud_skeleton.ply` | AdTree software output | The raw skeleton itself (vertices + radius + edges) |
| `taper.txt` (in the AdQSM run folder) | AdQSM software output | Trunk taper curve, used to calibrate trunk radii |
| `BranchStructure.txt` (same folder) | AdQSM software output | Per-branch diameter/length table, used for median-radius-by-order calibration of non-trunk branches |
| `TreesParams.txt` (same folder) | AdQSM software output | Whole-tree totals, printed as an independent reference row |

**What it does, step by step** (matches the numbered comment block at the
top of the script):
1. Reads the binary `.ply` (vertices with x, y, z, radius; edges).
2. Merges coincident vertices, turning the raw (often disconnected)
   skeleton into one connected tree with real branch points.
3. Roots the tree at its lowest point and derives a unique parent per
   node (removes any loops).
4. **Smooths the centerline** (Laplacian smoothing) to remove zig-zag
   scanning noise — the root, junctions, and branch tips are never moved.
5. **Prunes** any branch whose radius is below `RADIUS_THRESHOLDS`
   (the trunk is never pruned, even if thin).
6. **Resamples** the dense skeleton into segments whose length adapts to
   local radius — long segments on the thick trunk, short ones on thin
   twigs.
7. Shifts x, y to the origin (z is left alone — it is height above
   ground).
8. **Optionally calibrates** the cylinder radii against AdQSM (trunk from
   the taper curve, other branch orders scaled by a single factor per
   order).
9. Prints/records volume, DBH, height, taper, trunk/branch length into
   `volume_results.csv`.
10. Saves the finished geometry to a `calib_*.npz` file (it does **not**
    write the final ANSYS text file itself any more — see Step 2 below).
11. Optionally shows/saves a 3D preview.

**Key parameters and what to set them to**

| Parameter | What it controls | Current value | Guidance |
|---|---|---|---|
| `TREE_NAME` | Which tree; drives every other path | `"IND01_054"` | The only thing you normally change to switch trees |
| `RADIUS_THRESHOLDS` | Prune any non-trunk branch thinner than this (metres). Give a *list* to produce several variants in one run. | `[0.005]` (5 mm) | Try a few values (e.g. `[0.005, 0.010, 0.020]`) and compare — this changes **topology** (which branches exist at all), not just detail |
| `SEG_LEN_MIN` / `SEG_LEN_MAX` / `SEG_LEN_K` | Adaptive resampling: target segment length = `SEG_LEN_K * local_radius`, clamped to `[SEG_LEN_MIN, SEG_LEN_MAX]` | `0.1 / 0.5 / 0.5` m | This is a *mesh-density* choice (like element size in ANSYS), it should **not** change the underlying volume — only how many cylinders represent it. Smaller `SEG_LEN_MIN` = more, shorter elements = smoother curvature but more cylinders for ANSYS to solve |
| `MIN_CYL_LEN` | Any cylinder shorter than this after resampling gets merged into its parent | `0.005` m | Safety net against zero-length beams at junctions; leave as-is unless ANSYS complains about degenerate elements |
| `MERGE_DECIMALS` | Rounding used to detect "same point" when merging vertices | `5` | Safe default, rarely needs changing |
| `RECENTER_XY` | Shift x/y to remove the scanner's georeferencing offset | `True` | Keep `True` unless you have a specific reason to keep absolute coordinates |
| `SMOOTH_ITERS` / `SMOOTH_ALPHA` | Laplacian smoothing passes / strength (0–1) | `5 / 0.5` | More passes = smoother but flatter centerline; radii are never touched by this step |
| `CALIBRATE_RADII` | Whether to replace AdTree's own radii with AdQSM-calibrated ones | `True` | Keep `True` unless you specifically want to see raw (heavily overestimated) AdTree radii — see §6.3 finding |
| `FIELD_DBH` | If you have an independently measured DBH, rescale the taper curve to match it exactly at 1.3 m | `None` | Set this only if you trust a field caliper measurement more than the AdQSM taper curve |
| `THIN_BRANCH_CUT_CM` | Diameter cut-off for the fair "vs. destructive reference" comparison | `10.0` cm | Matches the de Tanago field methodology (see AdQSM.pdf, Appendix A) — do not change unless your reference used a different cut-off |
| `TAPER_H_LOWER` / `TAPER_H_UPPER` | Heights (m) at which DBH / taper-upper diameter are measured | `1.3 / 10.0` | 1.3 m is the standard forestry "breast height"; keep these consistent across *every* method or the taper numbers become incomparable |

**Output:** `calib_<TREE_NAME>_r<mm>mm<adqsm-variant><seg-suffix>.npz` — the
calibrated geometry, ready for Step 2. Also several rows appended/updated
in `volume_results.csv` (raw AdTree, calibrated AdTree, optionally a
10 cm-filtered variant of each).

### Step 2: `export_geom_ansys.py`

**Purpose:** turn one already-calibrated `.npz` (from Step 1) into the
`geom_*.txt` file your ANSYS macro (`buk.mac`) actually reads. This step
does **no** calibration or comparison — it is a pure, fast format
conversion, so you can re-export to ANSYS any time without redoing Step 1.

**How to use it:**
1. Run it once with `NPZ_FILE = None` — it will just print every
   `calib_*.npz` file it finds in the current folder, so you can see what
   is available.
2. Copy one of the printed filenames into `NPZ_FILE`, re-run.
3. It writes `geom_*.txt` with 11 columns (x, y, z, length, radius,
   axis_x, axis_y, axis_z, parent, …, **branch order** as the 11th
   column).

**Required input:** the `calib_*.npz` file produced by Step 1 (nothing
else — everything needed is already stored inside it).

---

## 5. Method 2 — TreeQSM pipeline (MATLAB, `runsken.m`)

This runs entirely in MATLAB, cell by cell (`Ctrl+Enter` per `%%`
section). It is **not** part of the Python codebase, but its results feed
into it via `import_matlab_results.py` (§6.2).

**Required input file:** a plain-text point cloud, `<tree_id>.txt`
(three columns X Y Z, one tree only — remove other trees/ground/leaves as
much as possible beforehand, since TreeQSM interprets every point as
wood).

**What it does, section by section (high level):**
1. Sets `src_dir` (TreeQSM MATLAB code) and `data_dir` (your point clouds
   + where results are written).
2. *(Optional, guarded by two switches + a dry-run)* — cleans out all
   previous outputs for a completely fresh start.
3. Loads the point cloud, shifts it so it starts at the origin, plots it
   for a visual sanity check, and saves it as a `.mat` file.
4. **Defines input parameters** — either automatically
   (`define_input(P,1,1,1)`) or manually (see the parameter table below).
   If you define mannual run nPD1, nPD2min, nPD2 max = 1, 
   For automatic set manual_patchdiam = false!
5. Runs a **first, coarse batch** of `n_models_first` reconstructions
   (default 5) — TreeQSM's cover-set generation is randomized, so several
   models are always built even for one parameter combination.
6. **Selects the optimal model(s)** from that batch (`select_optimum`,
   based on point-to-cylinder distance).
7. Runs a **second, larger batch** (`n_models_opt`, default 25) using the
   optimal parameters, to get a reliable mean volume and standard
   deviation.
8. Exports volume/DBH/height/taper/trunk-and-branch-length tables and
   (optionally) a `geom_*.txt` for ANSYS, plus a `Filtered <10cm` group
   matching the destructive reference's methodology.
9. Detects and (optionally) removes **disconnected "branch islands"** —
   see the key finding in §6, item 4.9 — cylinders that TreeQSM could not
   reconnect to the main structure, which otherwise break the ANSYS
   import.

**Key parameters and recommended ranges**

TreeQSM's own manual and the AdQSM paper (which tuned TreeQSM on similar
tropical trees) both point at the same handful of parameters:

| Parameter | What it controls | Guidance from the papers |
|---|---|---|
| `PatchDiam1` | Size of the *first* (coarse, uniform) cover set — used mainly to strip non-tree points and get an initial branch layout | TreeQSM manual: **5–15 cm is a good range**; the exact value matters little for the final result. AdQSM paper (Table 4, tuned per study site): Peruvian 0.12 m, **Indonesian 0.08 m**, Guyanese 0.10 m |
| `PatchDiam2Min` | Minimum patch size of the *second* cover, used at branch tips | AdQSM paper: Peruvian 0.04 m, **Indonesian 0.02 m**, Guyanese 0.01 m. Smaller = more detail captured but risk of a branch being wrongly split into several |
| `PatchDiam2Max` | Maximum patch size of the *second* cover, used at the stem base | AdQSM paper: **0.07 m** for all three sites (this is the most consistent of the three) |
| `MinCylRad` | Smallest cylinder radius TreeQSM will fit | `0.0025` m in your current setup — a physically-motivated floor, not usually swept |
| `nPD1`, `nPD2Min`, `nPD2Max` | How many candidate values to test per parameter (grid search) | Set to `1` for a manual/fixed value (your current setup); increase to try several values and let `select_optimum` pick the best — TreeQSM's own manual notes this "usually only 3 parameters" (`PatchDiam2Max`, `PatchDiam2Min`, and cylinder length `lcyl`) are worth grid-searching |
| `n_models_first` / `n_models_opt` | Reconstructions per parameter combination (coarse / refined) | `5` / `25` — more models = tighter standard deviation but longer runtime; TreeQSM's cover-set generation is intentionally randomized so this is not optional |
| `TAPER_H_LOWER` / `TAPER_H_UPPER` | Same convention as AdTree — 1.3 m / 10.0 m | **Must match every other method** for taper numbers to be comparable |

> Rule of thumb from the TreeQSM manual: **smaller patch size → more
> detail but a near-quadratic increase in compute time**, and very small
> patches can wrongly split a branch that isn't fully covered by points.
> Bigger patches are faster and merge small branches together, and may
> make the fitted cylinder at a branch's start too large (it "eats"
> points from the child branch).

---

## 6. Comparing results

### 6.1 `qsm_volume_mean.py` — averaging TreeQSM's published (de Tanago) realizations

Reads 20 pre-existing published TreeQSM reconstruction files for a tree
(not your own MATLAB run) and reports mean/std/CV of stem, branch and
total volume across the realizations, then upserts one row into
`volume_results.csv` under the label `"TreeQSM de Tanago (mean)"`.

**Required input:** the `qsm_tanago` folder of 20 cylinder-model text
files for that tree (columns: radius, length, …, branch order — see the
script header for exact column indices).

### 6.2 `import_matlab_results.py` — bringing your MATLAB TreeQSM runs in

Scans a folder for `volumes_*.csv` tables (written by `runsken.m`),
picks one `Group` per file (default `"Optimal (single)"`), and
upserts a row per file into `volume_results.csv`, along with DBH,
height, taper, and trunk/branch length read from matching sidecar `.txt`
files in the same folder.

**Required input:** the `volumes_<tree>_<run>.csv` table plus its sidecar
`dbh_/height_/taper_/trunklen_/branchlen_<tree>_<run>.txt` files, all
written by `runsken.m`.

**Parameter to know:** `IMPORT_GROUP` — which of TreeQSM's result
groups to import (`"Optimal"`, `"Optimal (single)"`, `"Estimated"`,
`"All inputs"`, `"Simplified"`, `"Filtered <10cm"`, …). `"Estimated"`
combines both MATLAB runs and gives the best standard deviation;
`"Filtered <10cm"` is the one that gets `branch_filter="10cm"` for the
fair reference comparison.

### 6.3 `reference_volume.py` — the destructive ground truth

Reads the field measurement table (`IND.h.trees.txt`, de Tanago dataset
format), keeps only the chosen tree's rows, sums section volumes for the
fractions you switch on (`stem`, `branch`, `twig`, `leaf`, `stump`,
`root`, `fruit` — by default only `stem` + `branch`), and computes
DBH/taper directly from the measured stem sections. Always writes with
`branch_filter="10cm"`, because the field crew never measured anything
thinner.

**Required input:** `IND.h.trees.txt` in the tree's data folder.

### 6.4 `compare_volumes.py` — the actual comparison

Reads `volume_results.csv` and prints two clearly separated sections:

- **Section A — "vs. reference"**: only `branch_filter == "10cm"` rows.
  This is the fair, apples-to-apples accuracy check against the
  destructive measurement (`REFERENCE_METHOD = "Reference (destructive)"`).
- **Section B — "methods vs. each other"**: only `branch_filter == "none"`
  rows (full reconstructions, thin branches included). The destructive
  reference can never appear here, so AdQSM's own reported numbers
  (`REFERENCE_METHOD_NONE = "AdQSM (TreesParams)"`) are used as a common
  yardstick instead — **this is not an accuracy claim about AdQSM**, just
  a way to give every method something to be compared against.

For each section it prints a per-tree table, then Bias/MAE/RMSE/CV-RMSE
across all trees, plus separate one-line-per-method summaries for
DBH/height/taper/trunk-length/branch-length error (kept out of the volume
metrics since they are not volumes).

**Parameter to know:** `SELECT_TREE` — a specific tree ID, or `"ALL"`.

### 6.5 `plot_volumes.py`

Draws the same comparisons as charts, saved as PNGs into a `plots/`
subfolder: total volume by tree, error box plots, Bias/MAE/RMSE bar
charts, and a per-tree 2×3 overview grid (volume, DBH, height, taper,
trunk length, branch length). Everything is drawn twice — once for the
`"10cm"` mode, once for the `"none"` mode — using one shared
`{method: colour}` mapping per mode so the same method always gets the
same colour.

---

## 7. Where TreeGraph fits in

TreeGraph is your **actively developed, primary** reconstruction pipeline
(`build_graph.py`, `estimate_radius.py`, `taper.py`), but it currently
lives **outside** the `volume_results.csv` / `compare_volumes.py` loop —
adding it there was noted as a planned-but-postponed step. Unlike AdTree
and TreeQSM, its purpose is not "get the closest total volume" but to
guarantee three things ANSYS structural analysis actually needs:

1. **Geometric connectivity** — no floating/disconnected cylinders (the
   same class of problem found in the TreeQSM branch-island issue, §4.9
   of the historical notes below).
2. **Simplifiability** — a reducible cylinder count.
3. **No physically absurd branch-size jumps** — no child cylinder more
   than **2× its parent's radius**. This 2× rule is not arbitrary: the
   TreeGraph paper itself derives it empirically from manual branch
   measurements (daughter branches were occasionally thicker than their
   parent, but never more than 2×) and applies it as an upper-bound
   radius correction.

Current focus: diagnosing and fixing **260 radius-ratio violations** on
specimen `IND01_054` (config `MIN_PTS=5, CLUSTER_SIZE=0.10, k=1`), coming
from a block-constant radius bug in `estimate_radius.py`'s `p2a_2()`
function and a one-sided (upper-bound-only) cap in `taper.py`'s
`radius_correct()`. A three-layer repair (log-linear interpolation between
reliable neighbours, a global top-down 2× cap, and a "dominant child"
anchor rule for large furcation nodes) has been dry-run tested and is
being turned into production code (`treegraph_convert.py`) on branch
`treegraph-radius-fix-and-export`. Once stable, it should get its own
entry in this guide's parameter table and its own row(s) in
`volume_results.csv`.

---

## 8. File-dependency quick reference

| Script | Reads | Writes |
|---|---|---|
| `adtree_reconstruct_compare.py` | `<TREE_NAME> - Cloud_skeleton.ply`, AdQSM `taper.txt` / `BranchStructure.txt` / `TreesParams.txt` | `calib_*.npz`, rows in `volume_results.csv` |
| `export_geom_ansys.py` | `calib_*.npz` | `geom_*.txt` (for ANSYS `buk.mac`) |
| `runsken.m` (MATLAB) | `<tree_id>.txt` point cloud | `volumes_*.csv` + sidecar `.txt` files, optional `geom_*.txt` |
| `import_matlab_results.py` | `volumes_*.csv` + sidecar files from `runsken.m` | Rows in `volume_results.csv` |
| `qsm_volume_mean.py` | 20 published TreeQSM realization files (`qsm_tanago/`) | Optional per-realization CSV, one row in `volume_results.csv` |
| `reference_volume.py` | `IND.h.trees.txt` | One `"10cm"` row in `volume_results.csv` |
| `compare_volumes.py` | `volume_results.csv` | Console report only |
| `plot_volumes.py` | `volume_results.csv` | PNGs in `plots/` |
| `check_adqsm_thin_branch_length.py` | `BranchStructure.txt`, `TreesParams.txt` | One diagnostic PNG in `plots/` (standalone; does not touch `volume_results.csv`) |
| `tree_geom_utils.py` | — (shared library, no top-level code) | — (imported by the other Python scripts) |

---

## 9. Open items / things to double check next session

- TreeGraph is not yet feeding `volume_results.csv` — decide the method
  label and whether/when to add it once the radius-fix branch lands.
- Crown reconstruction quality above ~16 m remains an acknowledged open
  problem across methods (point density drops, branches thin out below
  what any method can reliably resolve).
- Keep `TAPER_H_LOWER = 1.3 m` / `TAPER_H_UPPER = 10.0 m` and
  `THIN_BRANCH_CUT_CM = 10.0` consistent across *every* script — a
  mismatch here silently makes taper/branch_filter comparisons invalid
  without any error being raised.

---

*Last updated: reflects the AdTree/AdQSM/TreeQSM/Python pipeline as of the
scripts reviewed on 2026-09-01. Update this file whenever a script's
parameters, file formats, or the TreeGraph integration status change.*
