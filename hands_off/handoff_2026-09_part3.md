# Handoff — 2026-09, part 3

## Purpose & context

Continuation of the Babice tree reconstruction/volume-comparison project (TreeQSM, AdQSM,
AdTree, TreeGraph → ANSYS structural export). This session was dominated by a deep-dive
into **AdQSM segmentation (Height_Segmentation) instability** on reference tree `B21_S01`,
plus several infrastructure fixes surfaced along the way. Ends with a fresh TreeQSM finding
that reframes the next reconstruction pass.

---

## Immediate next step (stated goal for the next chat)

Bára will run **further TreeQSM reconstructions** for beech production trees, tuning
**PD2Min (PatchDiam2Min)** specifically (see finding below), and will **evaluate results
using `pmdist_mean` / `pmdist_trunk_mean` / `pmdist_branch_mean`** — the point-to-cylinder
distance columns added to `params_*.csv` this session (see Infrastructure below).

---

## Newest finding (this session, not yet folded into a written decision doc)

- **TreeQSM `MinCylRad` is already at the right value and should NOT be tuned further.**
  It's currently set at the point-cloud resolution boundary (~2.5 cm, matching the
  subsampling). Empirically: increasing it visibly degrades the reconstruction; decreasing
  it has no effect at all. Conclusion: **leave `MinCylRad` as-is** for TreeQSM going
  forward — this closes out the "two MinCylRad candidates" exploration plan discussed
  earlier this session (default 0.0025 vs. a `1.5×Res`-derived alternative); that plan is
  now superseded for TreeQSM by this direct empirical result. (Note: this finding is
  reported for the tree(s) just tested this session — treat as a working hypothesis to
  re-confirm if it's applied to a materially different tree/point-cloud density.)
- **`PD2Min` (PatchDiam2Min) is the parameter that meaningfully moves the reconstruction**
  for TreeQSM — this is now the priority parameter for the manual PatchDiam search on
  beech production trees, more so than `MinCylRad`.

---

## Open item carried over, still unresolved

- **B21_S01's manual PatchDiam values were never finalized.** An AUTO-mode probe was run
  and produced a console "best three cases" table + `Optimal model: 58`, but the console
  table was initially misread (columns didn't line up with the `PatchDiam1/2Min/2Max`
  auto-generated candidate lists shown just above it) and self-corrected mid-session. The
  reliable source — `params_B21_S01_aut_mo8_sr005_ri0.csv`'s `pd1_m`/`pd2min_m`/`pd2max_m`
  columns (built from `OptInputs`, always scalar and correct for the winning model) — was
  requested but **never actually supplied/confirmed**, because the conversation pivoted
  into the AdQSM investigation. **Before starting new manual TreeQSM runs, re-pull those
  three values from that CSV (or re-run the AUTO probe) rather than reusing the
  console-table numbers.**

---

## This session's major body of work: AdQSM segmentation (Height_Segmentation) sensitivity

Full findings are written up in **`AdQSM_Segmentation_Sensitivity_B21_S01.docx`** (created
this session, in the project). Summary:

- AdQSM's "segmentation variant" folders (`040`–`100`) encode **Height_Segmentation (HS)**
  × 100 — confirmed against AdQSM's own `Introduction.pdf`/GitHub README, whose worked
  example (HS = 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) matches variant numbers 50–100 exactly.
- Dense 61-variant sweep on `B21_S01` (reference DBH 56.0 cm): **44/61 cluster near the
  reference, 17/61 are catastrophic outliers** (up to ~40× normal volume). No smooth
  relationship to HS value — adjacent settings can span the full range (e.g. HS 0.73/0.74/
  0.75 = 10%-off / best-match / worst-outlier-in-dataset).
- Even the 44 "close" ones are **systematically biased ~10% low** vs. the measured
  reference (mean −9.9%), not randomly scattered — only 9/44 fall inside a strict ±10%
  band.
- Structural relevance: a 10% DBH underestimate → ~19% underestimate in cross-section
  (D²) → **~34% underestimate in bending stiffness (D⁴)** → ~1.5× overestimated deflection
  in ANSYS if uncorrected. This is bigger than the volume effect alone.
- Root cause: AdQSM's own paper (Fan et al. 2020, *Remote Sensing* 12(18):3089) names
  initial-trunk-cylinder point selection as its **largest source of error**, and recommends
  repeated modeling — same mitigation as the HS-averaging guidance.
- Mechanism, confirmed via `taper.txt` comparison (worst outlier variant 040 vs. best
  match 091): the two curves are the **same shape, scaled by a near-constant ~3.105×**
  across 0–24.6 m (95% of trunk height) — points to a single anchoring-value error (likely
  the initial cylinder radius), not shape distortion or accumulating error. (A crown-tip
  divergence initially flagged as a second anomaly was later shown to be a **ratio
  artifact** of both curves' diameters approaching zero there — not a real defect.)
- **Existing safety net doesn't catch this**: `tree_geom_utils.py`'s
  `_reject_taper_spikes()` only drops an isolated single-row spike >2× both neighbors; it
  structurally excludes first/last rows and never fires on undersized values or a
  uniformly-scaled whole section. Connects to a pre-existing, still-open
  `TODO_investigations.md` item 8 (same class of problem previously confirmed on
  `IND07_083`, an even longer ~22.6 m flat/wrong section).
- **Outlier-correction comparison**: with 7 candidate values (HS 0.4–1.0) and up to 3/7
  being outliers, classic Tukey-fence removal (Q3+1.5×IQR) **removes 0 of 7** (outliers are
  spread out enough to inflate the IQR past their own values) — worse than doing nothing.
  Plain **median remains the most robust simple choice**, though only one bad value away
  (4/7) from failing the same way.
- **Decision going forward**: per tree, run HS 0.4–1.0 (7 settings); select the variant
  closest to that tree's own measured DBH as the primary reconstruction; **separately**
  build an independent median variant (see below) as a comparison point — not by copying
  whichever real variant happens to be closest to the median (that would make the
  comparison meaningless).
- **`adqsm_build_median_variant.py`** (new script, built & spot-validated this session):
  builds a synthetic variant folder (`999`) from N source HS variants — per-cylinder
  median for `taper.txt` and `BranchStructure.txt` (row-correspondence across variants
  verified byte-identical on all 21,263 `B21_S01` rows for order/parent/length/height/
  angle/azimuth/zenith, so per-cylinder — not just per-order — median is valid), and
  independent per-field median for `TreesParams.txt` (deliberately NOT copying a single
  "closest" row, so `TreeVolume` may not exactly equal `Trunk+Branch` — expected and
  documented). **Still pending**: loading `"999"` into `adtree_reconstruct_compare.py`'s
  `ADQSM_VARIANTS` and running the full pipeline end-to-end (Bára's step, not yet done as
  of this handoff).
- Unexplained, flagged but not investigated: **exact byte-identical duplicate outputs**
  for HS pairs 0.52/0.53 and 0.58/0.59.
- The ~9.6 m near-flat region at the trunk base (present in both good and bad
  reconstructions) is a known/expected artifact of AdQSM's initial-cylinder averaging —
  deliberately deprioritized, not today's main problem.

---

## Infrastructure fixes made this session (all validated by Bára unless noted)

- **`pmdistance` export**: `runsken.m` section 19 now writes `pmdist_mean`,
  `pmdist_trunk_mean`, `pmdist_branch_mean` (from `QSM.pmdistance`, already computed
  per-model inside `treeqsm.m`, independent of `select_optimum`) into `params_*.csv`, and
  all four `upsert_result()` copies (`tree_geom_utils.py`, `import_matlab_results.py`,
  `qsm_volume_mean.py`, `reference_volume.py`) carry the 3 new columns in sync. Confirmed
  safe for AdTree rows (blank via existing optional-field fallback) and for
  `compare_volumes.py`/`qsm_volume_mean.py`/`reference_volume.py` reads (column-name-based
  or not reading `volume_results.csv` directly). **This is the metric the next session's
  pmdist evaluation will use.**
- **`.py.bak` cleanup**: 5 stale files archived (not deleted) to
  `scripts/python/bak_archive_2026-09-06/`, after confirming each had a committed,
  tracked counterpart.
- **`WRITE_THIN_BRANCH_FILTERED_ROW`** (`adtree_reconstruct_compare.py`) now also skips
  the underlying >=10cm-filter *computation* (not just the CSV write) for both AdTree's
  own thin-branch stats and the separate AdQSM-side `report_adqsm_thin_branch()` console
  report — confirmed no other code path reads the now-`None` values.
- **`MinCylRad` run_tag/label/grouping fix**: `compute_run_tag()` now encodes
  `man_MinCylRad` (byte-identical tag at the historical default 0.0025; new `_mcr` token
  otherwise); `plot_volumes.py` label and `plot_box.py`'s **actual box-grouping key**
  (previously silently pooling different `MinCylRad` values into one box — a real
  data-correctness bug, now fixed) both updated; `parameter_sensitivity.py`'s
  `TREEQSM_PARAMS_FOR_CORR` now includes `mincylrad_m`.
- **PD1/PD2Min/PD2Max precision fix**: 2-decimal (cm) rounding caused real collisions
  (e.g. `0.005` and `0.01` both → `"01"`, same run_tag → silent reconstruction reuse).
  Fixed with a conditional cm/mm-width token (`pd_token()` in `runsken.m`,
  `_pd_field_token()` mirrored in `plot_volumes.py`/`plot_box.py`) — clean cm values render
  unchanged, non-clean values get a wider, unique token. Note: this **does** visibly widen
  existing AUTO-mode chart labels (AUTO's PD values are almost never clean cm) — accepted
  tradeoff, confirmed not a grouping/data bug (`plot_box.py`'s grouping key already used
  raw floats).
- **New shared `scripts/python/plot_style.py`**: centralizes color palette
  (`FAMILY_GRADIENTS`, now also a `family_shades(family, n)` helper for gradient-sampling
  N distinct shades when a chart needs more lines than the base palette has stops — used
  by `taper_curve_compare.py` for its up-to-~21-variant chart), line widths, marker sizes,
  and font sizes. `plot_volumes.py`, `plot_box.py`, `parameter_sensitivity.py` migrated as
  a near-pure refactor (one approved deliberate change: axis-label font unified 8→9pt,
  sourced from `plot_volumes.py`'s side only). `adqsm_variant_sensitivity.py` given its own
  heavier local style tier (`linewidth=3`, `markersize=5`, legend `fontsize=9`) plus
  `FAMILY_GRADIENTS` coloring, per Bára's explicit numbers.
  **Standing rule now in memory**: all plotting scripts should keep style settings as
  labeled, configurable constants near the top, preferably sourced from this shared module.
- **`taper_curve_compare.py` style update**: per-variant curves now colored via
  `family_shades("AdQSM", n)` (clamped away from near-white shades) instead of the
  matplotlib default cycle; font sizes unified to `plot_style` constants; per-variant
  `markersize`/`linewidth` deliberately left unchanged (readability at up to 21 lines).
  **Validation from Bára still pending as of this handoff.**
- **New `runsken.m` section "1B) ARCHIVE SELECTED TAGS (OPTIONAL)"**: same safety
  architecture as `CLEAN START` (master switch, dry-run default `true`, typed
  ARCHIVOVAT/SMAZAT confirmation, archive-not-delete default, reuses `clean_patterns`/
  search dirs/`Cloud`-name guard), but filters to files whose name contains one or more
  Bára-specified tag substrings — for archiving specific unwanted MATLAB variants without
  touching everything. Note: the existing `"1B) USER SETTINGS"` cell was renamed to
  `"1C)"` to make room (comment-text-only rename, no logic change).
  **Validation from Bára still pending as of this handoff.**
- **Usage note, not a bug**: `plot_volumes.py`'s `METHOD_FILTERS` list requires every
  entry to have its own `"family"` key (`resolve_method_filters()` reads it
  unconditionally) — a missing key raises `KeyError: 'family'`. Bára hit this from her own
  edit and fixed it; documenting so it's not re-diagnosed from scratch next time.

---

## Reminders carried forward

- Manual backup of `plots/` before any Claude Code session touching the filesystem —
  used consistently this session, keep it up.
- Anti-deletion safety line at the start of every filesystem-touching Claude Code prompt —
  used consistently this session.
- Staged workflow (diagnostic → confirm → implement → validate with real printed numbers)
  — used throughout; several diagnostics this session caught real design flaws before
  implementation (e.g. the `treeqsm.m` cross-product-vs-paired-candidates issue, the
  `BranchStructure.txt` row-correspondence check, the Python PD-token collision check).
