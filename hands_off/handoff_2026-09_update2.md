# Handoff — AdTree/AdQSM/TreeQSM comparison pipeline (update to handoff_2026-09.md)

## Immediate next step for the new chat

Bára just asked Claude Code to fix two things found while testing Phase A
(see "Phase A" below) — a documentation clarification in `plot_volumes.py`'s
new `METHOD_FILTERS` mechanism, and extending the same method-filtering to
`plot_error_boxplot()`/`plot_error_metrics_bar()`, which currently show
EVERY method unfiltered. That fix was just handed to Claude Code; **confirm
it ran successfully before starting Phase B.** After that, proceed through
Phases B → C → D in order (see "On the horizon" below) — do not skip ahead
to C or D before B is resolved, and do not assume a Spearman-correlation
methodology for Phase C without discussing it first (explicitly undecided,
see below).

## Purpose & context

Bára is building a research pipeline for 3D tree reconstruction from
terrestrial laser scanning (TLS) point clouds, cross-comparing wood-volume
estimates across independent reconstruction methods (TreeQSM/AdQSM,
AdTree, TreeGraph), calibrating AdTree's cylindrical models against AdQSM
and a destructive field reference, and exporting calibrated models to
ANSYS for structural analysis (static tilt/displacement, first natural
frequency — this last part is not yet started, see below).

**Reference trees (destructive volume data available):** `IND01_054`
(~24.9m), `IND03_088`, `IND07_083` — all three now have COMPLETE data in
`volume_results.csv` for `branch_filter="10cm"` (AdTree calibrated + raw,
AdQSM direct, TreeQSM, and the real `"Reference (destructive)"` row).

**Production trees (no destructive reference):** 12 trees still to be
reconstructed and exported to ANSYS, once a parameter configuration is
chosen using the 3 reference trees. B21_S01 and Ave_B, mentioned in
earlier sessions as secondary test cases, turned out NOT to be the trees
Bára actually had destructive data ready for this round — `IND03_088`/
`IND07_083` were used instead. Confirm with Bára whether B21_S01/Ave_B are
still planned as additional reference trees, or are among the 12
production trees, if this becomes ambiguous later.

## Current state

### CSV schema (`volume_results.csv`, 27 columns)
```
tree, method, total_m3, trunk_m3, branch_m3, std_m3, dbh_m, height_m,
taper_cm_per_m, trunk_len_m, branch_len_m, branch_filter, n_cylinders,
mode, pd1_m, pd2min_m, pd2max_m, mincylrad_m, simp_maxorder,
simp_smallradii, simp_replaceiterations, adqsm_variant,
radius_threshold_mm, seg_min_mm, seg_max_mm, seg_k_pct, calmethod
```
`calmethod` (values: `"regression-perorder"` or `"min5mm"`, blank for
non-AdTree-calibrated rows) was JUST added this session, to ALL FOUR
`upsert_result()` copies (`tree_geom_utils.py` canonical + the three
sibling scripts `import_matlab_results.py`/`qsm_volume_mean.py`/
`reference_volume.py`) AND the read side (`load_results()` in
`compare_volumes.py`) in a single pass, with a mandatory backup taken
first — verified via full diff against the backup: 1958 rows before and
after, zero non-`calmethod` value changes. **This was necessary because
without a structured `calmethod` column, turning on
`COMPUTE_CALREF_MIN5MM` would have caused `select_rows()`/
`assign_adtree_groups()` to silently merge `regression-perorder` and
`min5mm` rows into the same group** (same bug class as the earlier
"AdTree raw wrongly carries adqsm_variant" issue). `COMPUTE_CALREF_MIN5MM`
is currently still `False` — no `min5mm` rows exist yet for any of the 3
reference trees. Flipping it on is Phase B's decision, not yet made.

**This CSV schema has a known fragility, treat with care going forward:**
any future new column must be added to all four write copies AND
`load_results()` in the same pass, with a backup taken first. This has
caused near-incidents twice this session already (both caught before real
data loss, but only through careful staging).

### Scripts/tools status
- `parameter_sensitivity.py`: three chart types (facet grid, ranked
  deviation dot plot, Spearman correlation heatmap) for AdTree/TreeQSM
  parameter sensitivity vs. a resolved reference (destructive for
  `branch_filter="10cm"`, dynamically-resolved smallest-numbered AdQSM
  variant for `"none"`). Currently single-tree only (`SELECT_TREE`) — NOT
  yet extended to compare across the 3 reference trees (that's Phase C).
  Method labels now use `shorten_method_label()` (matches `plot_volumes.py`
  style). `FACET_Y_MODE` toggle added: `"absolute"` (raw metric value) or
  `"relative"` (% deviation from reference, same formula as the ranked
  plot).
- `tree_coverage_check.py`: deliberately SIMPLIFIED from an earlier
  "coverage gap" cartesian-product design (which produced confusing
  false-positive "missing combinations" for TreeQSM, since `mode` and
  `(pd1, pd2min, pd2max)` are coupled, not independent axes) down to a
  plain "list the parameter values used for this tree, copy-paste ready"
  report — plus a CSV export (`coverage_<tree>.csv`, one row per
  parameter, values joined with `"; "`). Use this before reconstructing
  each new tree to keep parameter grids comparable.
- `plot_volumes.py` — **Phase A, in progress:**
  - DONE: `METHOD_FILTERS` (structured filter dicts, e.g. `{"family":
    "adtree", "calmethod": "regression-perorder", "adqsm_variant": [...],
    ...}`) added alongside the existing literal-string `SELECTED_METHODS`
    (the two UNION together) for `plot_tree_overview()`. Verified: the
    example filter resolved to the exact same 8 methods as the old
    hand-typed list, chart visually identical.
  - PENDING (prompt just given, not yet confirmed run): (1) clarify in the
    `METHOD_FILTERS` comment that `"family"` only accepts `"adtree"`/
    `"treeqsm"` — `Reference (destructive)`/AdQSM-direct/TreeQSM-mean rows
    must go in `SELECTED_METHODS` instead (they union automatically); (2)
    `plot_error_boxplot()` and `plot_error_metrics_bar()` currently show
    EVERY method across ALL trees with NO filtering at all — extend the
    same `SELECTED_METHODS`/`METHOD_FILTERS` restriction to both, via a
    new tree-agnostic `resolve_selected_methods_multi_tree()` wrapper
    (these two charts compare across all trees at once, unlike
    `plot_tree_overview()` which renders one tree at a time).
  - NOTE: `plot_error_boxplot()` is ALREADY a cross-tree comparison chart
    (one box per method, one point per tree) — this may turn out to
    overlap significantly with Phase C's goals once that starts; worth
    checking whether Phase C extends this existing chart or builds
    something new in `parameter_sensitivity.py` instead of duplicating.

### Documentation
- `CHANGELOG_adtree.md`: new entry documenting the `SEG_LEN_MIN` floor
  finding — the diameter below which the fixed-length floor applies is
  `2 * SEG_LEN_MIN / SEG_LEN_K`, which reaches much further up the tree
  than the parameter's name suggests (e.g. `SEG_LEN_MIN=200mm` at the
  project's typical `SEG_LEN_K=0.5` floors almost the ENTIRE tree below
  80cm diameter, not just "thin twigs"). Two candidate mechanisms for why
  this shifted `total_m3` so much in `branch_filter="10cm"` data were
  discussed — leading hypothesis: it refits different per-order
  `(a, b)` regression coefficients (`[calmethod=regression-perorder]`),
  rescaling every cylinder of that order, not just ones near the 10cm
  boundary — but this is **UNVERIFIED**, see `TODO_investigations.md`
  item 1.
- `TODO_investigations.md` (new file, consolidating previously-scattered
  open items):
  1. Regression coefficients vs. `SEG_LEN_MIN` — unverified hypothesis,
     next step defined (compare fitted `(a,b)` per order across two
     `SEG_LEN_MIN` values).
  2. Frustum-accurate cylinder volume — deferred, needs coordinated ANSYS
     export format + TreeQSM equivalence check first.
  3. TreeGraph B21_S01 native crash (`0xC06D007F`) — unresolved.
  4. TreeGraph crown reconstruction quality (2x radius safeguard
     violations above ~16m on IND01_054) — unresolved.
  5. `COMPUTE_CALREF_MIN5MM` — status uncertain before this session
     (handoff said "prepared, not implemented"); now clarified: it's
     `False`, zero `min5mm` rows exist for the 3 reference trees, and
     `calmethod` column exists to keep it safe when it IS turned on —
     that decision is Phase B.
  6. Chart group ordering / reference-family colour refinement for
     `"none"` mode — still pending, no new detail.
  7. Interactive (Plotly/HTML) `parameter_sensitivity.py` charts —
     deliberately deferred (static facet grid judged sufficient for now).
- `TreeQSM_vs_AdQSM_Radius_Summary.docx` (new): documents a completed
  investigation — TreeQSM's own `MinCylRad` parameter (a hard floor on
  cylinder radius, e.g. `0.0025m` for one tested run) was identified as
  the likely dominant driver of TreeQSM reporting ~110% more total volume
  than AdQSM in `branch_filter="none"` mode, via a quantile-matched
  per-branch-order radius comparison (ratio grows from 1.8x at order 1 to
  8x at order 8 — consistent with a floor effect, not uniform bias).
  Caveats noted: based on ONE TreeQSM run and ONE AdQSM variant; compares
  TreeQSM against AdQSM, not against destructive truth (no reference
  exists for `"none"` mode). A concrete testable next step (re-run with
  lower `MinCylRad`) was proposed but NOT started.
- `AdTree_Smoothing_Resampling_Summary.docx` and
  `AdTree_Calibration_Summary.docx`: pre-existing background docs, still
  current, no changes this session.

## On the horizon — the 4-phase roadmap (current priority order)

Agreed this session, when Bára said she was "getting lost in the amount
of information." Do NOT reorder without her explicit input — the order
matters (each phase's output feeds the next).

- **Phase A — simplify `plot_volumes.py` method selection.** Mostly done
  (see "Scripts/tools status" above) — one follow-up fix pending
  confirmation. Finish this first.
- **Phase B — decide AdTree `calmethod`: `"regression-perorder"` vs.
  `"min5mm"`, going forward.** NOT STARTED. This decision should be made
  using the 3 reference trees' data (now that `calmethod` is a safe,
  structured column). Turning on `COMPUTE_CALREF_MIN5MM` to generate the
  `min5mm` comparison data is part of this phase, not before.
- **Phase C — extend `parameter_sensitivity.py` across the 3 reference
  trees; general cross-tree comparison; sensitivity analysis (optionally
  per-tree or combined).** NOT STARTED. Two things Bára explicitly wants
  input on, not yet resolved:
  - Whether/how to use Spearman correlation here — **UNDECIDED, discuss
    with Bára before implementing**, do not assume the same per-tree
    heatmap approach from `parameter_sensitivity.py` is automatically the
    right multi-tree design.
  - Whether an x/y scatter plot with a fitted regression line is wise —
    Bára was advised (this session) that it's reasonable ONLY for
    genuinely continuous, roughly-monotonic parameters (e.g.
    `radius_threshold_mm`, `seg_min_mm`), NOT for categorical-ish ones
    (`adqsm_variant`, `mode`) where the facet grid/ranked plot remain
    primary — re-confirm this framing still holds before building it.
- **Phase D — recommend an optimal AdTree parameter combination for the
  12 production trees**, trading off `n_cylinders` (target <=8000,
  ideally <=4000 — Bára noted `radius_threshold_mm=20` is "already quite
  a lot" on the coarse side) against retained detail/mass (deviation from
  reference). NOT STARTED, depends on B and C. Suggested approach discussed
  but not built: a scatter of `n_cylinders` (x) vs. `%` deviation from
  reference (y), with vertical reference lines at 8000/4000, to visually
  identify Pareto-efficient combinations rather than picking one number.

**Deferred, not on the critical path:** Bára raised the idea of comparing
ANSYS structural OUTPUTS (max tilt/displacement under a defined load,
first natural frequency) across reconstruction methods/parameters, since
`branch_filter="10cm"`'s destructive-reference comparison says nothing
about how much the (unreferenced) thin crown branches affect structural
response. This was deliberately set aside — "ANSYS zatím odložíme, pojďme
na krok 1-3" — to revisit only after the volume-based parameter decision
(A-D) is closed out. When revisited: the plan discussed was to run ANSYS
on only a SHORT LIST of candidate configs (2-4) surviving the volume-based
selection, not the full parameter sweep, and compare just those two output
numbers across methods. How Bára currently extracts those two numbers
from ANSYS (macro output file vs. manual GUI read) was never established
— ask before assuming a format.

## Key learnings & principles (new this session, additive to handoff_2026-09.md's list)

- **`SEG_LEN_MIN`'s floor reaches much higher up the tree than its name
  suggests:** the diameter below which the fixed-length floor applies is
  `2 * SEG_LEN_MIN / SEG_LEN_K`, not just "thin twigs" — always compute
  this explicitly before assuming a `SEG_LEN_MIN` change only affects
  small branches.
- **The 4-copy `upsert_result()` CSV schema is fragile by construction:**
  any new column needs all four write copies + `load_results()` updated
  in the same pass, backup first. Caused two near-incidents this session.
- **TreeQSM's `MinCylRad` creates a hard radius floor**, causing
  disproportionate volume inflation specifically for thin/high-order
  branches relative to AdQSM — not a uniform bias, worse for thinner
  branches.
- **Destructive reference only exists for `branch_filter="10cm"`**; the
  `"none"`-mode reference is dynamically resolved
  (`resolve_reference_method_none()`) to the smallest-numbered AdQSM
  variant PRESENT for that tree — this number is NOT stable across
  re-runs/trees (was "05", later "04" for the same tree after more data
  existed) — never hardcode an expected variant number.
- **Spearman correlation understates non-monotonic/categorical parameter
  effects** (e.g. `adqsm_variant` identity isn't a smooth numeric axis) —
  treat the facet grid and ranked-deviation plot as the primary sensitivity
  tools, correlation heatmap as a numeric-screening supplement only.
- **"AdTree raw" rows must never carry `adqsm_variant`** (already fixed,
  `None`/blank) — raw AdTree geometry doesn't depend on which AdQSM
  variant was active during that reconstruction loop iteration.
- **A geometric/volume match to a destructive reference does not by
  itself guarantee similar downstream ANSYS structural behaviour** — full
  crown mass/stiffness distribution (including branches the `"10cm"`
  filter excludes entirely) can matter for tilt/displacement and natural
  frequency even without a volume ground truth for that part of the tree.

## Approach & patterns (reaffirmed, unchanged from handoff_2026-09.md)

- Strict staged workflow: read-only diagnostic → explicit confirmation →
  implementation with diff summary → visual verification. No step skipped,
  especially for anything touching the shared CSV schema.
- English for code/comments/changelogs/prompts-to-Claude-Code; Czech for
  conversation.
- Prefer smaller, safer diffs; new methods added alongside existing ones
  until a decision is confirmed.
- Never guess at the contents of a file not yet shared this session — ask
  for it or request an upload rather than assuming its structure.
- Handoff documents prepared at the end of a session (or when a person
  feels they're "getting lost in the amount of information") to keep
  continuity across new chats — this document is exactly that.

## Tools & resources (updates only — see handoff_2026-09.md for the full original list)

- New scripts this session: `parameter_sensitivity.py`,
  `tree_coverage_check.py`.
- New/updated docs: `CHANGELOG_adtree.md` (new entry),
  `TODO_investigations.md` (new file),
  `TreeQSM_vs_AdQSM_Radius_Summary.docx` (new).
- `volume_results.csv`: 27 columns (added `calmethod`), 3 reference trees
  (`IND01_054`, `IND03_088`, `IND07_083`) with complete `branch_filter=
  "10cm"` data including a real destructive reference row each.
