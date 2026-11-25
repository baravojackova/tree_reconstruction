# Handoff — TreeQSM/AdTree Comparison Pipeline (Session 2026-09)

> For starting a fresh chat with continuity. Paste this in, or attach it,
> at the start of the new conversation.

## Project context

Same project as `handoff_master.md` (AdTree/AdQSM/TreeQSM volume comparison
for Babice dataset, test case IND01_054, eventual ANSYS export goal). This
handoff covers a long follow-up session that touched `runsken.m` (MATLAB,
`scripts/matlab/`) heavily, plus `plot_volumes.py`, `plot_box.py` (new
script), `compare_volumes.py`, `import_matlab_results.py`,
`tree_geom_utils.py`, `qsm_volume_mean.py`, `reference_volume.py`
(all Python, `scripts/python/`).

## Shared CSV schema change (affects 4 files)

`volume_results.csv` gained 8 new columns, appended after `n_cylinders`:
`mode, pd1_m, pd2min_m, pd2max_m, mincylrad_m, simp_maxorder,
simp_smallradii, simp_replaceiterations`. These record TreeQSM's actual
reconstruction parameters (manual vs. auto, patch diameters, simplification
settings) — previously NOT captured anywhere except a hand-typed run tag.

**All four scripts that write this CSV must stay schema-identical**:
`tree_geom_utils.py`, `import_matlab_results.py`, `qsm_volume_mean.py`,
`reference_volume.py` — each has its own local `upsert_result()` copy; all
four were updated this session to the same 20-column header. AdTree/AdQSM/
reference/TreeQSM-published rows leave the 8 new columns blank (they have
no reconstruction-run metadata); only `import_matlab_results.py`'s
"TreeQSM mine (...)" rows populate them, sourced from a new
`params_<tree>_<run_tag>.csv` sidecar file that `runsken.m` now writes.

**Naming**: the volume field is now `trunk_m3`/`"Trunk volume [m^3]"`
everywhere (renamed from `stem_m3`/"Stem volume" — matches the pre-existing
`trunk_len`/`trunk_vol` internal naming). `stem` still appears as external
vocabulary inside `reference_volume.py` (source file's own `Fraction`
column) and in a few internal diagnostic variable names — intentionally
untouched, out of scope.

## runsken.m changes, in order

1. **Figure numbering fix**: point cloud plot moved to an explicit
   non-colliding figure number (avoids colliding with `simplify_qsm.m`'s
   hard-coded figure(1)/figure(2), which are library code and were
   deliberately left untouched).
2. **`run_tag` split into `recon_tag` + `run_tag`**:
   `recon_tag` = mode + PD values only (drives `res_name`/`res_new_name`,
   i.e. the expensive reconstruction .mat files).
   `run_tag` = `recon_tag` + simp settings (drives everything else: `vol_file`,
   `export_prefix`, dbh/height/params files, `simplified_file`).
   Both auto-computed from current settings via a shared local function
   `compute_run_tag()` (used both in the "DERIVED NAMES" section and in a
   section-16 sanity check — see below). **This split was necessary**:
   changing only `simp_ReplaceIterations` between variants must NOT trigger
   a full expensive reconstruction re-run, only `recon_tag` governs that.
3. **Skip-if-exists reconstruction guard**: new switch
   `overwrite_existing_reconstruction` (default `false`). If
   `res_file`/`res_new_file` already exist under the current `recon_tag`
   (checked via `isfile(fullfile('results', res_file))` — note the
   `results/` subfolder prefix is REQUIRED, `make_models`/
   `make_models_parallel` save there, not to a bare filename — a bare
   `isfile(res_file)` would incorrectly report false), the expensive
   `make_models`/`make_models_parallel` calls are skipped and the existing
   file is loaded instead.
4. **Section 16 sanity check**: at the top of section 16 (right after
   `clear QSM_simple_clean`, itself a fix for stale-variable leakage across
   repeated in-session simplification passes — see "Workflow: trying
   multiple simplification variants" below), `run_tag` is recomputed from
   current settings and compared against the live `run_tag` — errors out
   loudly if they don't match (catches "changed simp_*/manual_patchdiam
   but forgot to re-run section 2 (DERIVED NAMES)").
5. **Gated ANSYS export**: `export_to_ansys` (default `false`) gates
   section 20 entirely; `ansys_export_name` (default blank) optionally
   overrides the auto-generated tag used for `geom_*.txt` filenames.
6. **Simplified model auto-saved to disk**: every simplification pass now
   unconditionally saves `simplified_<tree_id>_<run_tag>.mat` (previously
   only lived in the MATLAB workspace, lost when trying a second variant).
   `EXPORT_FROM_SAVED_RUN_TAG` (blank by default) lets section 20 export
   from ANY previously-saved variant instead of whatever's currently live
   in memory.
7. **`clean_start` archive-or-delete choice**: new switch
   `clean_delete_permanently` (default `false` = safe/archive mode — moves
   matched files into a single timestamped `archive/<timestamp>/` folder
   instead of deleting; `true` = old permanent-delete behavior, useful for
   genuine trial-run cleanup). Typed confirmation word changes to match
   ("ARCHIVOVAT" vs "SMAZAT").
8. **`params_*.csv` export bug (found + fixed)**: section 19 originally read
   `inputs.PatchDiam1/2Min/2Max` directly — these can be MULTI-ELEMENT
   ARRAYS whenever `nPD1`/`nPD2Min`/`nPD2Max` > 1 (a real, deterministic
   `define_input.m` behavior, not an edge case), corrupting the fprintf
   output (scrambled scientific notation, stray extra line) for auto-mode
   runs specifically. **Fixed** by reading from `OptInputs` instead (the
   `select_optimum()`-picked WINNING scalar values — confirmed via full
   source trace through `select_optimum.m`/`treeqsm.m` that
   `OptInputs.PatchDiam1/2Min/2Max/MinCylRad` are always genuine scalars,
   for both auto and manual mode, and are the semantically-correct
   "actually used" values, not just the first candidate tried).
9. **PROPOSED BUT NOT YET CONFIRMED DONE**: auto-saving all 5 runsken.m
   figures (point cloud, optimal model, simplified model, islands,
   simplified-no-islands) as PNG via a `save_figures_png` switch. A full
   prompt was drafted (naming convention: `pointcloud_<tree_id>.png`,
   `optimal_<tree_id>_<recon_tag>.png`, `simplified_<tree_id>_<run_tag>.png`,
   `islands_<tree_id>_<run_tag>.png`,
   `simplified_no_islands_<tree_id>_<run_tag>.png`) but **Bara has not yet
   confirmed running it** — check this first in the new chat before
   assuming it's done.

### Workflow: trying multiple simplification variants without re-running reconstruction

Because of point 3 above (skip-if-exists), Bara's normal workflow for
trying several `simp_*` variants on the SAME reconstruction (same
mode+PD) is:
1. Edit `simp_MaxOrder`/`simp_SmallRadii`/`simp_ReplaceIterations` in
   section "1c) USER SETTINGS".
2. Re-run ONLY section "2) DERIVED NAMES" (recomputes `run_tag`).
3. Skip sections 3–15 (reconstruction) entirely.
4. Continue from section 16 onward.
This works within one continuous MATLAB session (relies on `QSM_opt`/`res`
staying in the workspace). **Never run section "1b) CLEAN START"** in the
middle of this — it wipes/archives everything.

## import_matlab_results.py changes

1. **`IMPORT_GROUPS`** (was `IMPORT_GROUP`, singular): now a list, imports
   several named groups per file in one run (default:
   `["Optimal", "Simplified (no islands)", "Filtered <10cm"]`). Each
   (file, group) pair becomes its own CSV row.
2. **`if __name__ == "__main__":` guard added** — this file previously
   executed its RUN section unconditionally even on a bare `import`, which
   caused ONE confirmed accidental real write during unrelated testing
   (data was fine, just unintended timing) before the guard was added.
   **This was also, separately, the root cause of a real data-loss
   incident** (see below) via a DIFFERENT mechanism (schema mismatch, not
   the import guard) — don't conflate the two.
3. **`read_params_file()`** added — reads the new `params_<tree>_<run>.csv`
   sidecar once per file (same "read once, reuse across every group from
   that file" pattern as the existing dbh/height/taper reads).

## plot_volumes.py changes

- Earlier small fixes: shortened method labels (regperorder→reg),
  annotation top-margin/headroom fix (scales with `ANNOTATION_FONTSIZE`),
  bar width/spacing parameters, `SELECTED_METHODS` (filter which methods
  plot_tree_overview shows — `None` = all), `PNG_FILENAME_SUFFIX`.
- `FAMILY_GRADIENTS` (dict of family→hex gradient, including a "Reference"
  pink family) and `classify_family()` extracted from
  `build_method_color_map()` as module-level, importable pieces (used by
  `plot_box.py`).
- **`shorten_method_label()` redesigned for TreeQSM**, now accepts optional
  `pd1, pd2min, pd2max, simp_maxorder, simp_smallradii,
  simp_replaceiterations, mode` kwargs (default `None`, AdTree/AdQSM
  unaffected). New compact format when all params are known:
  `"TQ_%s_%s_m%d_s%d_r%d_%s" % (mode_short, pd_token, maxorder, smallradii*1000,
  replaceiter, stage_short)` → e.g. `"TQ_man_p8-2-7_m8_s5_r0_Filt"`,
  `"TQ_aut_p7-2-7_m8_s5_r1_Filt"`. Mode is ALWAYS shown explicitly (both
  "man" and "aut" — an earlier implicit-from-PD-presence design was
  deliberately rejected in favor of always-explicit). Falls back to the
  old `"TQ_{run}_{stage_short}"` form when params aren't available (old
  pre-`params_*.csv` rows). New helpers: `treeqsm_stage_short()`,
  `treeqsm_pd_token()`, `treeqsm_mode_short()`. All 8 call sites across the
  file updated for full consistency (not just `plot_tree_overview` —
  `plot_total_volume_by_tree`, `plot_error_boxplot`,
  `plot_error_metrics_bar`, and the RUN diagnostic printout too).

## plot_box.py (new script this session)

Box plots comparing GROUPS of related method variants (not individual
methods) per metric field (same 9 fields as `plot_tree_overview`). Built
up over many iterations:

- `GROUP_RULES`: list of `(regex, group_label)` tuples, first-match-wins,
  matched against `shorten_method_label(method)` — user-editable, drives
  AdTree/AdQSM grouping. Currently only 4 of the originally-sketched 7
  rules are active (AdTree raw + 3 seg-calib groups) — AdQSM/TreeQSM rules
  were superseded by the newer, separate TreeQSM mechanism below.
- Reference row (`REFERENCE_METHOD`) and TreeQSM-published-mean row(s)
  (`TREEQSM_REF_METHODS` — two mutually-exclusive strings depending on
  `BRANCH_FILTER`) are pulled OUT of grouping entirely and drawn as
  horizontal reference lines (`SHOW_REFERENCE_LINE`,
  `SHOW_TREEQSM_REF_LINE`), each in its own distinct color, with a legend.
- Boxes use `whis=(0,100)` (whiskers = true min/max, no IQR-outlier dots —
  deliberate, since n per group is small). Individual data points are
  jittered and overlaid on top (`JITTER_POINT_SIZE`, `JITTER_POINT_ALPHA`,
  `JITTER_DARKEN_FACTOR`, `JITTER_RANGE`). Per-point on-plot labels show the
  varying parameter (`POINT_LABEL_PATTERNS` for AdTree, regex-based;
  structured-column-based for TreeQSM — see below), with vertical stagger
  (`LABEL_VERTICAL_STAGGER`) to reduce overlap when values cluster. Colors:
  box edge/whisker/median/jitter/label all use `darken_color()` with
  independently-tunable factors — **note the historical gotcha**: this
  function's `factor` parameter was ONCE inverted relative to its own doc
  comments (factor=1 meant "unchanged", not "black") — this was fixed; if
  you ever see a darken factor behaving backwards again, check this first.
- **`overwrite_existing_reconstruction`-style structured TreeQSM grouping**
  (separate from `GROUP_RULES`): `TREEQSM_STAGE_FILTER` (one of "Optimal",
  "Simplified", "Simplified (no islands)", "Filtered <10cm", or `None`),
  `TREEQSM_VARY_BY` (list of short param names that become the per-point
  spread inside a box — everything else becomes the box-defining group
  key), `TREEQSM_PARAM_FILTERS` (optional dict restricting to specific
  parameter values, e.g. `{"pd1": 0.08}`, matched with `1e-4` float
  tolerance), `TREEQSM_PARAM_SHORT_NAMES` (short name ↔ CSV column
  mapping). Routing: any row with `method.startswith("TreeQSM mine (")`
  goes through this new mechanism instead of `GROUP_RULES`; both
  mechanisms' resulting groups are merged into ONE combined chart (AdTree/
  AdQSM groups first, TreeQSM groups appended, ordered by first CSV
  appearance) — confirmed working via direct testing, not just "no
  crash". Old pre-`params_*.csv` TreeQSM rows (`mode == ""`) are excluded
  from this mechanism with a batched warning (not a regression — they
  never matched any `GROUP_RULES` pattern either).
- Box/point label format for TreeQSM matches `plot_volumes.py`'s new style:
  `"TQ_" + "_".join(fixed tokens)` for box labels (e.g.
  `"TQ_man_p7-2-7_m8_s5"`, no `r` token since that's the varying dimension
  shown per-point), plain short tokens (e.g. `"r0"`) for point labels.

## Data-loss incident (fully resolved, but worth knowing about)

Mid-session, `volume_results.csv` was found wiped to just its header row.
Root cause: `qsm_volume_mean.py`/`reference_volume.py` still had their OLD
12-column `upsert_result()` (from before this session's 8-new-column
schema change) — `csv.DictWriter`'s default `extrasaction='raise'` crashed
on the FIRST `writerows()` call against a file that now had 20-column rows,
but `csv.DictWriter` opens in `"w"` mode (truncates immediately) — so the
crash left the file with only its header, zero data rows. **Both scripts'
`upsert_result()` are now fixed** (20-column schema, matching the other two
files). Recovered via: restore from most recent good backup
(`volume_results.csv.pre_diag_backup`, predating the crash) → discovered
that backup was itself slightly stale (missing 18 rows from a later
same-day import) → re-ran `import_matlab_results.py` against the restored
file to regenerate those 18 rows from still-present source files → full
recovery confirmed via row-count/diff check. **Lesson reinforced**: always
verify a restored backup against what SHOULD be present (check timestamps/
row counts), don't assume the most recent backup is complete.

## Two further bugs found, one fixed, one superseded by a rebuild decision

1. **Bug 1 (fixed)**: `params_*.csv` corruption for auto-mode runs — see
   "runsken.m changes" point 8 above.
2. **Bug 2 (superseded)**: 3 stale `"TreeQSM mine (v1aut, ...)"` rows in
   `volume_results.csv` had mismatched mode/PD values (traced to a
   pre-session artifact, confirmed NOT a live/recurring bug in current
   `read_params_file()` logic — exact origin never pinned down, but ruled
   out as a current-code risk). Rather than patching those 3 rows, **Bara
   opted for a full rebuild** of `volume_results.csv` instead (see below).

## PENDING: full `volume_results.csv` rebuild — check status first in new chat

Bara decided to fully rebuild `volume_results.csv` from scratch rather than
patch individual bad rows, for certainty. A read-only reconnaissance
confirmed **92 of 98 rows are fully reproducible** by re-running, in any
order: `adtree_reconstruct_compare.py`, `import_matlab_results.py`,
`qsm_volume_mean.py`, `reference_volume.py`. **6 rows are NOT
reproducible** (their source `volumes_IND01_054_v1aut.csv`/`v1man.csv` no
longer exist on disk) — 3 `v1aut` rows (already known-bad, losing them is
a feature not a cost) and 3 `v1man` rows (confirmed CORRECT, legitimate
data — Bara explicitly signed off on losing these too, in favor of a clean
rebuild). **Bara said she would run this herself** (not via a Claude Code
prompt) — **status unconfirmed as of this handoff**. Check in the new chat
whether `volume_results.csv` now has ~92 rows (rebuilt) or still has the
old 98 (not yet done) before doing further TreeQSM-related analysis, since
downstream charts/groupings will look different either way (the 3
`aut_mo8_sr005_*` rows currently show blank structured params until BOTH
the rebuild happens AND those 3 MATLAB runs are re-executed with the
now-fixed Bug 1 `params_*.csv` export).

## Working style reminders (carried over, still true)

- Bara is experienced in ANSYS, basic in Python — every script/step needs
  plain-language explanation of what it does, not just what changed.
- Strict staged workflow: read-only diagnostic → explicit confirmation →
  implementation with diff summary → visual verification reminder.
- Confirm design/methodology BEFORE drafting a Claude Code implementation
  prompt, especially for anything touching shared file schemas or
  multi-file consistency (the `upsert_result()` header, `GROUP_RULES`-style
  mechanisms, etc.).
- Never silently touch `simplify_qsm.m`/`plot_cylinder_model.m` (TreeQSM
  library code) — work around it in `runsken.m` instead, as already done
  for the figure-numbering fix.
- Backups (`*.pristine_backup`, incrementing numbers) before any edit to
  the live `volume_results.csv`.
- Conversation in Czech; code/comments/changelogs/handoffs in English.
