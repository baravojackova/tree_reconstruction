# Handoff — AdTree/AdQSM/TreeQSM comparison pipeline (update to handoff_2026-09.md)

## Immediate next step for the new chat

Bára is about to configure `runsken.m` (MATLAB, TreeQSM) to try **several
manual `PatchDiam` (pd1/pd2min/pd2max) combinations** — the goal is to find
values appropriate for the **beech (buk) production trees**, since the
IND-group's published manual values (`pd1=0.08, pd2min=0.02, pd2max=0.07`)
were confirmed this session to be **specific to that forest/species/scan
group**, not universally transferable (see "Key learnings" below). Walk
through `runsken.m`'s PatchDiam section step by step (Bára is
comfortable with MATLAB syntax but wants each parameter's role explained
before changing it) — do **not** assume the IND values are a safe
starting point without discussing scope (how many combinations, whether
to test on a beech tree yet or continue using the 3 IND reference trees
for now) first.

## Purpose & context

Bára is building a research pipeline for 3D tree reconstruction from
terrestrial laser scanning (TLS) point clouds (Babice dataset),
cross-comparing wood-volume estimates across independent reconstruction
methods (AdTree, AdQSM, TreeQSM, TreeGraph), calibrating AdTree's
cylindrical models against AdQSM and a destructive field reference, and
exporting calibrated models to ANSYS for structural analysis.

**Reference trees** (destructive volume data available): `IND01_054`
(~24.9 m), `IND03_088`, `IND07_083` — used for every calibration/
sensitivity/reliability decision so far.

**Production trees**: 12 beech (buk) trees — **no destructive volume
reference**, but **do have their own structural reference measurements**
(first natural frequency, tilt, trunk deflection). The plan (agreed this
session) is to use the IND-tree-based screening to narrow parameter
grids, then run a smaller, targeted validation pass directly on the
beech trees against their own structural reference before final ANSYS
candidate selection.

## Current state — this session's work (large session, 4 new/updated documents)

### Phase B (AdTree `calmethod`) — closed out, deeply validated
`calmethod = "regression-perorder"` was already the working assumption;
this session built the full evidentiary/documentation trail:
- `calmethod_decision_summary.py` (new) — Bias/MAE/RMSE/tree-spread for
  `min5mm` vs. `regression-perorder`, pooled across all radius_threshold/
  adqsm_variant combos, for total/trunk/branch, both branch_filter modes.
- **`AdTree_Calmethod_Decision.docx`** — full write-up. Regression wins
  clearly on branch (lower bias+RMSE+tree-spread); trunk is byte-identical
  between the two methods (trunk calibration doesn't depend on
  `calmethod` at all — see "Key learnings").
- **`Understanding_Bias_and_RMSE.docx`** — plain-language explainer (dart-
  target analogy) requested separately, ties back to the calmethod numbers.
- **IND07_083 root-cause investigation** (triggered by "why is trunk
  calibrated *worse* than raw?"): AdTree's trunk radius comes **entirely**
  from AdQSM's own `taper.txt` curve (none of AdTree's own measurement is
  used — see `apply_order_calibration_factors()`/
  `apply_radius_regression_per_order()` in `tree_geom_utils.py`).
  IND07_083's `taper.txt` is implausibly **flat** (~0.708 m diameter from
  0–22.6 m) — a real AdQSM data-quality issue, confirmed via manual
  inspection of the raw file (not a pipeline bug). **New script
  `taper_curve_compare.py`** overlays every AdQSM variant's taper curve
  for a tree — reusable for beech trees later. Discovered along the way:
  **`TrunkVolume ratio ≈ (DBH ratio)²` almost exactly**, for every
  variant pair on all 3 trees (error < 0.15 percentage points) — AdQSM
  variants scale near-uniformly.
- Amended into **`AdTree_Calmethod_Decision.docx`** (Section 5/8 caveat)
  and written up in full in **`AdTree_Smoothing_Resampling_Summary.docx`
  Section 4** (new): explains the *mechanism* — `convert()` collapses each
  cylinder to `r = 0.5*(r_start+r_end)` and treats it as a constant-radius
  cylinder instead of a true frustum, which **always underestimates**
  volume (0% bias at no taper, up to −25% at a segment tapering to a
  point) — worse for `SEG_LEN_MIN`-floored thin twigs. Plausible mechanism
  behind the previously-observed non-monotonic `seg_min_mm` effect.

### AdTree parameter sensitivity ("Krok 1") — fully done
Extended `parameter_sensitivity.py` from single-tree to multi-tree
(`REFERENCE_TREES`), for both `branch_filter` modes
(`BRANCH_FILTERS_TO_RUN`), with **5 chart types**: facet grid, ranked
deviation (pooled, tree-marker-coded), correlation heatmap (**two modes**
— `raw` vs. `pct_dev`, the latter avoiding a real, measured cross-tree
scale confound), tornado, and ANOVA (`eta²`, main-effects only).
Documented in **`AdTree_Parameter_Sensitivity_Summary.docx`** (+ later
extended with the taper-variant-comparison section 6.1/6.2 above). Key
findings:
- **`radius_threshold_mm`** flips from irrelevant (`"10cm"` mode — masked
  by the bigger ≥10cm reference cut) to **dominant** (`"none"` mode — the
  actual geometry that goes into ANSYS) for both `branch_m3` and
  `n_cylinders`.
- **`seg_min_mm`** has a real but **non-monotonic ("sawtooth")** effect on
  `branch_m3` (Spearman ≈ 0, ANOVA `eta²` up to 0.24) — Spearman missed
  it because it isn't a smooth trend; confirmed via a dedicated facet grid.
- **`seg_max_mm`**: negligible everywhere — safe to fix at one value.
- **`adqsm_variant`**: dominates `trunk_m3`/`total_m3` (mechanistic —
  trunk = that variant's taper curve) but the "winning" variant is
  **tree-specific, not transferable** (confirmed both via Spearman-vs-ANOVA
  divergence and via `taper_curve_compare.py`'s direct visual/numeric check).
- **DBH-based variant selection method validated**: comparing each AdQSM
  variant's own `dbh_m` against the real measured DBH reliably picks a
  good-matching variant (re-validated later this session on all 3 trees,
  errors 0.7–4.0 cm on 35–84 cm-diameter trees).

### TreeQSM Part A — `simp_smallradii`/`simp_replaceiterations` sensitivity
Restricted to the `"Simplified (no islands)"` stage (the one that goes
to ANSYS) — this **is** a genuine 6-combination sweep (unlike
`pd1/pd2min/pd2max`, which are single-valued per tree/mode). Run
separately for `manual`-only, `auto`-only, and pooled (to rule out a mode
confound — turned out **not** to matter here, direction agreed in every
case).
- **Bug found and fixed mid-session**: filtering rows to the
  `"Simplified"` stage *before* resolving the reference row silently
  stripped the reference too (it's never itself a `"Simplified"` TreeQSM
  row) — tornado/ANOVA/`pct_dev`-heatmap silently fell back to a
  "tree's own mean" baseline for 4 of 5 metrics. **Fixed via a
  `reference_rows=` parameter** threaded through `_pct_dev_pooled()` and
  its callers, resolving the reference from the *unfiltered* dataset
  while the analyzed population stays stage-filtered. **This is now a
  standing rule for this whole project — see "Key learnings."**
- Corrected result: `simp_replaceiterations` controls `n_cylinders`
  almost entirely (`eta² = 0.86–0.98`), **negligible** effect on
  total/trunk/branch/DBH (< 0.2%) — a near-free cost lever.

### TreeQSM Part B — reconstruction-method reliability (AdTree vs. TreeQSM)
New script **`reconstruction_method_decision_summary.py`** — compares 5
groups (AdTree raw, AdTree calibrated, TreeQSM manual, TreeQSM auto,
TreeQSM de Tanago — a published external reference, confirmed complete
for all 3 trees/both stages it applies to) across **3 blocks** (different
TreeQSM stage + matching reference: `Filtered <10cm` vs. Destructive,
`Optimal` vs. AdQSM, `Simplified (no islands)` vs. AdQSM), later extended
with:
- **`AdQSM (other variants)`** group (`"none"`-mode blocks only) — how
  much AdQSM's own variants disagree with each other, as a noise-floor
  reference.
- **"Best-of-tree" selection** (per tree, the single AdTree/AdQSM
  combination matching all of total/trunk/branch best *simultaneously*)
  — Bára correctly flagged this is an **oracle/upper-bound** (uses
  knowledge of the answer), not achievable on a tree with no known
  reference.
- **"Realistic" AdTree pick** added alongside it: `radius_threshold_mm`
  fixed at 5 mm (matching TreeQSM's own finest baseline) +
  `adqsm_variant` chosen via the DBH-matching method — achievable
  without knowing the answer. Consistently **~2.3–2.4× worse** than the
  oracle for `total_m3`.
- Final chart styling: per-tree bars coloured by family×tree-shade
  (matching `plot_volumes.py`'s `FAMILY_GRADIENTS` convention), winning
  configs shown in a legend-style text box (not per-bar rotated labels).

**Headline result** (documented in **`TreeQSM_Reconstruction_Decision_Summary.docx`**,
new, covers both Part A and Part B):
- **`TreeQSM manual` wins on `total_m3`/`trunk_m3` in every block**,
  including against AdTree's *oracle* best-of-tree pick (not just the
  realistic one) — driven by AdTree's trunk-calibration weakness above.
- **`AdTree calibrated` wins specifically on `branch_m3`** (both oracle
  and realistic configurations) — the per-order regression calibration is
  genuinely good; only trunk is compromised.
- Practical implication flagged, **not built**: a **hybrid** export
  (TreeQSM-quality trunk + AdTree-calibrated branch) could beat either
  method alone on total volume.
- **`TreeQSM manual`'s published `pd1/pd2min/pd2max` values are
  IND-group-specific** (confirmed with Bára — not tuned to these 3
  trees, so no circularity in the reference-tree comparison, but **not
  transferable to beech** either) — this is exactly what tomorrow's
  `runsken.m` session needs to address.

### TODO_investigations.md
A prompt was just sent (end of previous message, **not yet confirmed
run** by Bára) to append 4 new items: the IND07_083-only taper.txt
check (open on the other 2 trees), the best-of-vs-realistic reporting
rule, the hybrid-export idea, and the beech-manual-parameters gap.
**Confirm this ran before assuming those items exist in the file.**

## Key learnings & principles (new this session, additive to prior lists)

- **Reference-row-stripped-by-filter bug pattern (now a standing rule)**:
  whenever rows are filtered down to a specific stage/group/family
  *before* looking up a reference row, the reference gets silently
  stripped too if it doesn't belong to that stage/group. **Always
  resolve the reference from the full, unfiltered tree/branch_filter
  population — never from an already-restricted subset.** Caught once
  this session (TreeQSM Part A); costly enough in wasted/wrong charts
  that it's worth checking explicitly any time a new stage/group filter
  is introduced.
- **Oracle ("best-of") vs. realistic — always report both.** A
  retrospective "best-matching configuration per tree" is an upper bound,
  not a deployable choice, on a tree with no known reference (like the
  beech trees). This project now computes both wherever "best-of"-style
  selection is used.
- **AdTree trunk calibration depends 100% on AdQSM's `taper.txt`** — no
  AdTree measurement survives calibration for the trunk. A bad
  `taper.txt` (confirmed: IND07_083; **not yet checked**: IND01_054,
  IND03_088) makes calibrated trunk *worse* than raw, not better.
- **`DBH ratio² ≈ TrunkVolume ratio`**, essentially exactly, across AdQSM
  variant pairs on all 3 trees — variants scale the whole tree by a
  near-uniform factor, which is *why* DBH-matching is such a reliable,
  cheap way to pick a good variant (Section 3.6 of the new doc).
- **DBH-based AdQSM variant selection is legitimate** (uses an
  independent, physically-measured ground truth, not the volume
  reference) **but is tree-specific** — must be re-run per tree
  (including each beech tree), never reused across trees.
- **The average-radius-as-cylinder volume approximation always
  underestimates**, worse for segments that taper more within their own
  length (0% at no taper, up to −25% tapering to a point) — a plausible
  mechanism behind `seg_min_mm`'s non-monotonic branch-volume effect,
  and the documented rationale for the still-deferred "frustum-accurate
  cylinder volume" TODO item.
- **`radius_threshold_mm`'s importance is mode-dependent** — irrelevant
  under the `"10cm"` reference cut, dominant in `"none"` mode (the mode
  that matters for the beech ANSYS geometry). Don't judge a parameter's
  importance from `"10cm"`-mode evidence alone.
- **Raw cross-tree pooling for correlation/ANOVA is a real, bidirectional
  confound** (can overstate OR understate a true effect depending on the
  metric) — `pct_dev` (normalize to each tree's own reference first) is
  the trustworthy default; `raw` is kept only for side-by-side comparison.
- **Tornado (extremes only) and ANOVA (all groups' variance) can
  legitimately disagree** — not a bug, they answer different questions;
  keep both rather than picking one.
- **TreeQSM `manual` parameter values are forest/species/scan-group-
  specific, not universal** — confirmed this session. Do not assume any
  published or previously-used `pd1/pd2min/pd2max` set transfers to a
  new group of trees without its own validation.

## Approach & patterns (reaffirmed, unchanged)

- Strict staged workflow: read-only diagnostic → explicit confirmation →
  implementation with diff summary → visual verification. No step
  skipped, especially for filesystem-touching changes.
- **Every Claude Code prompt touching the filesystem starts with an
  explicit anti-deletion safety line** — strengthened this session after
  a *third* near-incident risk was identified (the line now explicitly
  forbids "clean test environment" deletions, not just unscoped
  user-facing cleanup); Bára also now keeps manual backups of `plots/`
  before sessions, since prompt-text alone has failed before.
- English for code/comments/changelogs/prompts-to-Claude-Code; Czech for
  conversation.
- Additive changes only; new methods/groups/parameters added alongside
  existing ones until a decision is confirmed, never silently replacing.
- Sanity checks demand **real numbers from real runs**, not "it ran
  without crashing" — this project's Claude Code sessions consistently
  self-verify with manual hand-computations, MD5/byte-identity checks on
  unrelated output, and explicit before/after file counts.
- Handoff documents prepared at natural breakpoints or when Bára feels
  she's "getting lost in the amount of information" — this document is
  exactly that, for a genuinely large session.

## Tools & resources (updates only)

- **New scripts this session**: `calmethod_decision_summary.py`,
  `taper_curve_compare.py`, `reconstruction_method_decision_summary.py`.
- **Substantially extended**: `parameter_sensitivity.py` (multi-tree,
  tornado, ANOVA, raw/pct_dev heatmap modes, seg_min-view facet grid,
  dbh_m ranked deviation, TreeQSM Part A simp_* analysis with the
  reference-row fix).
- **Modified**: `plot_volumes.py` (`PLOT_DPI`/`TREES_TO_PLOT`/
  `BRANCH_FILTERS_TO_PLOT` params, `error_boxplot` label-rotation fix).
- **New documents**: `AdTree_Calmethod_Decision.docx` (updated with
  taper.txt caveat), `Understanding_Bias_and_RMSE.docx`,
  `AdTree_Parameter_Sensitivity_Summary.docx` (updated with taper-variant
  section), `AdTree_Smoothing_Resampling_Summary.docx` (updated with new
  Section 4 on frustum volume bias), `TreeQSM_Reconstruction_Decision_Summary.docx`.
- **New/updated data files**: `calmethod_decision_summary.csv`,
  `reconstruction_method_decision_summary.csv`,
  `reconstruction_method_best_picks.csv`.
- **`TODO_investigations.md`**: 4 new items queued (confirm the append
  ran — see "Current state" above).

## On the horizon (priority order, unchanged phases + this session's additions)

- **Immediate**: `runsken.m` manual PatchDiam configuration (see top of
  this document).
- Check `taper_curve_compare.py` output for IND01_054/IND03_088 (is the
  flat-taper.txt issue IND07_083-specific, or universal?).
- Determine beech-appropriate TreeQSM manual parameters (literature
  search or small testing pass) — blocks defaulting to `TreeQSM manual`
  for beech production reconstructions.
- Prototype the hybrid TreeQSM-trunk + AdTree-branch export (not
  started).
- **Phase C** (deferred from earlier in the project): AdTree vs. TreeQSM
  final-candidate matching for ANSYS — comparing *similarly-configured*
  models against each other, not just against references. Not started.
- Eventually: apply the narrowed AdTree parameter grid + validated
  TreeQSM manual params to reconstruct the 12 beech trees; run a
  targeted structural-sensitivity pass using the beech trees' own
  frequency/tilt/deflection reference; select final ANSYS candidates.
