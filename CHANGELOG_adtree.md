# Changelog - AdTree/AdQSM calibration investigation (2026-09-02)

Scope: `scripts/python/adtree_reconstruct_compare.py` and
`scripts/python/tree_geom_utils.py` (the AdTree -> AdQSM radius
calibration step of the pipeline). Unlike the existing TreeGraph
changelog (`SHRNUTI_zmen.md`, in Czech), this file is in English.

---

## Step 1 - diagnosed the self-referencing calibration bug

**Problem found:** `calibrate_cylinder_radii()` computed each branch
order's AdTree median radius from the SAME, already-pruned cylinder set
it was about to calibrate. Since pruning (`RADIUS_THRESHOLDS`) removes
the thinnest cylinders of each order first, a higher threshold
mechanically raised that order's AdTree median, shrinking the
calibration factor and over-rescaling even the thick, never-pruned
cylinders of that order.

**Consequence:** "AdTree calibrated" volumes - especially the
`>=10cm`-filtered branch volume used for the field-reference comparison
- ended up depending on `RADIUS_THRESHOLDS`, a parameter that is only
supposed to control which branches survive in the exported model, not
how they get calibrated.

---

## Step 2 - split the calibration into fixed-reference variants

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`

- Split `calibrate_cylinder_radii()` into
  `compute_order_calibration_factors()` (factor computation) and
  `apply_order_calibration_factors()` (factor application), so a factor
  dict computed from one cylinder set can be applied to a different one.
  `calibrate_cylinder_radii()` itself kept its exact old
  signature/behaviour, now implemented as a thin wrapper around the two
  new functions.
- Added two fixed reference cylinder sets, computed once per AdQSM
  variant (not once per threshold): `calref=unpruned` (from
  `convert(..., threshold=0.0, ...)`, the fully unpruned skeleton) and
  `calref=min5mm` (from a fixed 5 mm reference threshold, later
  generalized - see Step 3).
- Added corresponding `[calref=unpruned]` / `[calref=min5mm]` rows to
  `volume_results.csv`, purely additive alongside the existing
  self-referencing rows.

---

## Step 3 - tested calibration sensitivity to the reference threshold

**Files:** `adtree_reconstruct_compare.py`

- `calref=unpruned` turned out to massively overinflate branch volume -
  traced to noise contamination from AdTree's raw, unpruned
  micro-segments (lots of spurious very-thin cylinders inflating the
  order's radius distribution).
- Generalized the single fixed `min5mm` reference into a list
  (`CALIBRATION_REF_THRESHOLDS_MM = [0.002, 0.003, 0.004, 0.005]`),
  producing one `[calref=min2mm]` / `[calref=min3mm]` / `[calref=min4mm]`
  / `[calref=min5mm]` row family each.
- Result: the fixed-reference approach IS reproducible across
  `RADIUS_THRESHOLDS` for a given reference choice, but is still highly
  sensitive to WHICH reference threshold is picked - branch volume
  (`branch_m3`) ranged roughly **0.19-0.65 m3** depending on reference
  choice alone, for the same tree/AdQSM variant.
- Also added three finer pruning thresholds to `RADIUS_THRESHOLDS`
  itself (`0.002, 0.003, 0.004`, in addition to the existing `0.005,
  0.010, 0.020`) to widen the comparison sweep.

---

## Step 4 - new calibration method: global log-log regression

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`

Started a fourth calibration candidate, intended to eventually replace
the discrete per-order median-ratio approach entirely, once validated:
a single global power-law regression `AdQSM_radius = a *
AdTree_radius^b`, fitted with a robust Theil-Sen regression on
log-log, quantile-matched AdTree/AdQSM radius pairs (per branch order,
matching the two radius DISTRIBUTIONS rather than a single median
point).

- `parse_adqsm_branch_file_raw()`: new function exposing AdQSM's raw
  per-order radius lists (`parse_adqsm_branch_file()` now just reduces
  this to medians - unchanged behaviour/signature for existing callers).
- `build_quantile_matched_pairs()`: builds quantile-matched
  (AdTree, AdQSM) radius pairs per order, reusing the same "unpruned"
  reference cylinder set already computed for `calref=unpruned`.
- `fit_radius_regression()`: robust Theil-Sen fit of `(a, b)`, with an
  R^2 diagnostic printed (OLS, for fit-quality reporting only).
- `plot_radius_regression()`: saves a log-log scatter + fitted curve
  diagnostic PNG to `plots/adtree_adqsm_radius_regression_<tree>_
  AdQSM<variant>.png`.
- `apply_radius_regression()`: applies the fitted power law to
  calibrate non-trunk cylinders (trunk still comes from the taper
  curve, unchanged).
- New `[calmethod=regression]` rows added to `volume_results.csv`
  (`"none"` and `"(>=10cm only)"` variants), purely additive alongside
  every existing row.

---

## Step 5 - diagnosed order-dependent bias; added per-order (grouped) regression

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`

**Diagnostic findings that motivated this step** (from `IND01_054`,
AdQSM variant `08`):

- The global regression (`[calmethod=regression]`, Step 4) produced a
  `>=10cm`-filtered branch volume of **0.037 m3** - worse than
  `calref=min5mm`'s 0.196 m3, both far from the destructive field
  reference's **0.461 m3**. A single global `(a, b)` was suspected of
  averaging away real differences between branch orders.
- A read-only per-order diagnostic table (order, pair count, median
  AdTree/AdQSM radius, per-order ratio, AdTree radius range) confirmed
  it: order 1's own ratio was **~1.6**, clearly separate from every
  other order's **~2.1-2.35**. For this tree/variant, order 1 had 30
  quantile-matched pairs - above the `MIN_PAIRS_PER_ORDER = 15`
  threshold introduced in this step (see below).
- A closer look printed every one of order 1's 30 individual matched
  pairs (sorted by AdTree radius). The ratio trend across those 30
  points was smooth and monotonic, with no visible jump - so there is
  **no confirmed evidence of a hard mixed-population split** within
  order 1 itself. The check was originally motivated by a visual hint
  from a reconstruction render showing what looked like a codominant
  leader classified as branch order 1; that hint was not corroborated
  by the pair data, but 30 points is still a thin sample for a
  two-parameter power-law fit, so the underlying caution (order 1
  deserves its own fit, and any future merge should be loud, not
  silent) stands regardless.

**Change made:** a fourth calibration candidate,
`[calmethod=regression-perorder]`, extending the global regression with
per-order grouping instead of abandoning the regression approach:

- `MIN_PAIRS_PER_ORDER = 15`: new parameter controlling how many
  quantile-matched pairs an order needs before it gets its own fit.
- `group_orders_for_fitting()`: greedy upward merge of sparse orders
  (walking ascending/thickest-first, closing a group once it reaches
  `MIN_PAIRS_PER_ORDER`, folding any sparse leftover tail into the
  previous group). Prints the raw per-order pair counts and the final
  groups, and fires a loud, clearly-flagged console warning if branch
  order 1 ever ends up merged with another order instead of getting its
  own solo group.
- Per-group fitting reuses the existing, unmodified
  `fit_radius_regression()` on each group's own pooled pairs -
  `order_to_ab = {order: (a, b)}` maps every order to its group's fit.
- `apply_radius_regression_per_order()`: same contract as
  `apply_radius_regression()`, but looks up each cylinder's own order in
  `order_to_ab` instead of using one global `(a, b)`; leaves a
  cylinder's radius unscaled (with a one-time warning) if its order has
  no fit at all.
- `plot_radius_regression_per_order()`: same quantile-matched scatter as
  the existing regression plot, but with one fitted line per group
  (each spanning only that group's own AdTree radius range) instead of
  a single global line; adds an on-plot annotation if order 1 was
  merged. Saved separately as
  `plots/adtree_adqsm_radius_regression_perorder_<tree>_AdQSM<variant>.png`
  - the existing single-line plot/file is untouched.
- New `[calmethod=regression-perorder]` rows added to
  `volume_results.csv` (`"none"` and `"(>=10cm only)"` variants), purely
  additive alongside every existing row (self-referencing,
  `calref=unpruned`/`minXmm`, and the global `[calmethod=regression]`).

---

## Step 6 - per-order power-law tail-underestimation; added PCHIP interpolation

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`

**Diagnosed problem:** the per-order power-law fit from Step 5
(`[calmethod=regression-perorder]`) systematically UNDERESTIMATES at
both tails for order 1's own group. Order 1 has only 30
quantile-matched pairs, spanning roughly 3 orders of magnitude, with
visible ties/plateaus in AdQSM's own reported values - a 2-parameter
power-law is a poor fit for that shape. Concretely (order 1, sorted by
AdTree radius): the smallest matched pair (idx 0) is underestimated by
**-51%**, and the largest (idx 29) by **-42%** - the large end is
exactly the size range that matters most for the `>=10cm`-diameter
comparison against the field reference.

**Change made:** a fifth calibration candidate,
`[calmethod=interpolation-perorder]`, replacing the fitted power-law
curve with a **monotonic PCHIP interpolator** (fit in log-log space,
same convention as the power-law fit) - applied UNIFORMLY to every
order-group produced by the EXISTING `group_orders_for_fitting()`
grouping (reused unchanged, not recomputed; `MIN_PAIRS_PER_ORDER` also
unchanged), not just order 1, to keep the mechanism simple and
consistent:

- `fit_radius_interpolation()`: builds a `scipy.interpolate.
  PchipInterpolator` on log(AdTree radius) -> log(AdQSM radius) for one
  group's pooled pairs; de-duplicates tied x-values by averaging (PCHIP
  requires strictly increasing x); returns the interpolator plus the
  log-x range it was actually fitted over.
- `apply_radius_interpolation_per_order()`: same contract as
  `apply_radius_regression_per_order()`, but evaluates each cylinder's
  own order's interpolator instead of a power law - CLAMPS
  log(radius) to the fitted range first (never extrapolates), with a
  one-time warning per order if clamping was needed.
- `plot_radius_interpolation_per_order()`: same quantile-matched
  scatter as the regression-perorder plot, with the interpolated curve
  per group instead of the fitted line. Saved separately as
  `plots/adtree_adqsm_radius_interpolation_perorder_<tree>_AdQSM<variant>.png`
  - every existing plot/file from Steps 4-5 is untouched.
- New `[calmethod=interpolation-perorder]` rows added to
  `volume_results.csv` (`"none"` and `"(>=10cm only)"` variants), purely
  additive alongside every existing row.

---

## Step 7 - Decision and cleanup

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`,
`volume_results.csv`

**Decision:** the calibration-method investigation (Steps 1-6) is
concluded. Per-order regression (`[calmethod=regression-perorder]`) is
adopted as the **PRIMARY** calibration method going forward - its
cylinders now become the final `cyl` used for the exported `.npz`/
`geom_*.txt` geometry and the threshold's "processed, CALIBRATED"
report, replacing the old, buggy self-referencing
`calibrate_cylinder_radii()` call site. The fixed 5mm-reference
median-ratio calibration (`[calref=min5mm]`) is kept only as a
**SECONDARY/backup** reference point, computed and written alongside
the primary method for comparison. A full write-up with diagnostic
figures for this investigation exists as a separate presentation
document (not tracked in this repo).

**Removed** (superseded candidates from Steps 1-6, reasons documented
in the corresponding Step entries above):
- The old self-referencing `calibrate_cylinder_radii()` call site and
  its two untagged upsert rows (the original diagnosed bug). The
  function itself was also removed from `tree_geom_utils.py`, since it
  had no other callers left and leaving a known-buggy function around
  as dead code seemed more likely to cause confusion (or accidental
  reuse) than a clean removal - `compute_order_calibration_factors()`/
  `apply_order_calibration_factors()`, the two functions it wrapped,
  are both still used directly (by `calref=min5mm`) and remain.
- `calref=unpruned` and its rows (massively overinflated volume - noise
  contamination from unpruned micro-segments, Step 3).
- `calref=min2mm`/`min3mm`/`min4mm` and their rows -
  `CALIBRATION_REF_THRESHOLDS_MM` collapsed back to `[0.005]` (the
  list-based mechanism itself was left in place, not reverted to
  single-value code, since it already handles one entry with no extra
  complexity).
- The global (non-grouped) regression, `[calmethod=regression]`, and
  its rows - `plot_radius_regression()` and `apply_radius_regression()`
  were also removed from `tree_geom_utils.py` (grep-confirmed no other
  callers); `build_quantile_matched_pairs()` and
  `fit_radius_regression()` themselves are still reused by the primary
  per-order regression method and remain.
- Per-order PCHIP interpolation, `[calmethod=interpolation-perorder]`,
  and its rows - `fit_radius_interpolation()`,
  `apply_radius_interpolation_per_order()`, and
  `plot_radius_interpolation_per_order()` all removed from
  `tree_geom_utils.py` (no other callers).
- The two DIAGNOSTIC-ONLY console printouts added earlier to inform this
  decision (the per-order table and the order-1 individual-pairs dump) -
  removed for a cleaner console output now that the investigation is
  concluded and its findings are documented above; the same technique
  can be reapplied by hand if a future tree needs the same kind of
  diagnosis.
- `volume_results.csv`: 252 superseded rows removed (98 kept), across
  every tree in the file - see the row counts and surviving method list
  in this cleanup's own verification pass.

**Kept unchanged:** `compute_order_calibration_factors()`,
`apply_order_calibration_factors()`, `build_quantile_matched_pairs()`,
`fit_radius_regression()`, `group_orders_for_fitting()`,
`apply_radius_regression_per_order()`, `plot_radius_regression_per_order()`,
`MIN_PAIRS_PER_ORDER`, and every `[calref=min5mm]` /
`[calmethod=regression-perorder]` row already in `volume_results.csv`.

---

## Step 8 - AdQSM now has a ">=10cm only" comparison row

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`

`report_adqsm_thin_branch()` was already computing a `>=cut_cm`-kept
stem/branch volume (and now also length) from `BranchStructure.txt`,
but only printing it, not returning it - so AdQSM had no row at all in
the `branch_filter="10cm"` comparison against the destructive field
reference, unlike every AdTree-derived method. The function now returns
its results (or `None` in the no-length-column proxy-only fallback,
where there is genuinely no volume to report); `adtree_reconstruct_compare.py`
upserts one new `"AdQSM (BranchStructure, cyl. approx.)"` row per AdQSM
variant from it, independent of `RADIUS_THRESHOLDS`.

This is a **compromise value, not an official AdQSM total**: it's the
same single-diameter cylinder approximation `report_adqsm_thin_branch()`
already warned about (over-estimates volume for tapering branches - see
that function's own docstring), and `plot_volumes.py`'s existing
"AdQSM values in this >=10cm subset are approximate" warning (unchanged)
already covers it. `compare_volumes.py` and `plot_volumes.py` needed NO
changes - both already fold in any `branch_filter="10cm"` row generically.

---

## Step 9 - reverted the AdQSM BranchStructure approximation row

**Files:** `adtree_reconstruct_compare.py`, `volume_results.csv`

**Reverted** the `"AdQSM (BranchStructure, cyl. approx.)"` comparison row
added in Step 8: the underlying cylinder approximation was shown to
badly OVERESTIMATE, by a logically impossible amount - for all three
AdQSM variants tested (05, 06, 08), the `>=10cm`-filtered approximation
(`stem_vol_kept + branch_vol_kept`) exceeded even AdQSM's own OFFICIAL,
UNFILTERED whole-tree `BranchVolume` from `TreesParams.txt`. A filtered
subset can never legitimately exceed its own unfiltered total, so this
approximation cannot be trusted as a comparison value.

- `adtree_reconstruct_compare.py`: removed the `upsert_result()` call
  that wrote this row. `report_adqsm_thin_branch()` itself (and its
  Step-8 return value) is untouched - it's still called and still
  prints to console for manual reference, just no longer wired into
  `volume_results.csv`.
- `volume_results.csv`: removed the 3 existing
  `"AdQSM (BranchStructure, cyl. approx.) (AdQSM {variant})"` rows left
  over from Step 8's test runs (one per AdQSM variant: 05, 06, 08).
- `plot_volumes.py` was NOT touched: its unconditional "NOTE: AdQSM
  values in this >=10cm subset are approximate..." warning in
  `plot_tree_overview()` (`branch_filter=="10cm"` mode) already fires
  regardless of whether an AdQSM row is present, and correctly explains
  AdQSM's ABSENCE from this view again now that the row is gone.

---

## Step 10 - Open TODO: frustum-accurate cylinder volume

**Files:** none (documentation only - no code changed)

**Finding:** `convert()` (`tree_geom_utils.py`, ~line 1150) collapses
each cylinder's two endpoint radii into a single average radius
(`r = 0.5*(rad[anchor]+rad[n])`). `volume_stats()`/
`report_thin_branch_volume()` then treat every cylinder as
CONSTANT-radius (`V = pi*r^2*h`), not as a true frustum of a cone
(`V = pi*h/3*(r1^2 + r1*r2 + r2^2)`).

**Quantified impact for THIS specific simplification: small.** Even a
2x radius taper within one (adaptively short) segment underestimates
volume by only ~3.7% - much smaller than the AdQSM BranchStructure.txt
single-diameter-per-whole-branch issue (~100%+ overestimate, see Step 8/
Step 9 above), because AdTree segments are short
(`SEG_LEN_MIN`/`SEG_LEN_MAX`-bounded) relative to how much taper occurs
within any one of them.

**Why deferred (TODO, not implemented now):**
1. Preserving both endpoint radii ripples through many functions that
   currently assume one radius per cylinder (the calibration functions,
   the `>=10cm` diameter filter, the `geom_*.txt` export) - a wider
   change than first estimated.
2. Fixing only the reported volume numbers without also changing the
   ANSYS-facing export (`geom_*.txt` / `bk1.mac`) would create a
   mismatch between what's reported here and what's actually simulated
   in ANSYS - the ANSYS macro itself would need updating to use tapered
   (not constant-radius) beam/pipe elements, which is a
   structural-modelling decision outside this Python pipeline's scope.
3. For a fair three-way comparison (AdTree/AdQSM/TreeQSM), an
   equivalent frustum-accurate volume would need to be confirmed
   available from TreeQSM's own cylinder model/MATLAB export too -
   otherwise "fixing" only AdTree's number would introduce a NEW,
   asymmetric bias across methods, the mirror of the problem being
   solved.

**Status: deferred.** The current constant-radius approximation is
accepted as a small, bounded simplification (~3-4% max for realistic
taper) for present comparison purposes.

---

## Step 11 - Open TODO: SEG_LEN_MIN reaches much further up the tree than its name suggests

**Files:** none (documentation only - no code changed)

**Finding:** the resampling target length is `target_length = clamp(SEG_LEN_K
* r, SEG_LEN_MIN, SEG_LEN_MAX)` - a per-cylinder radius-proportional target,
clamped at both ends. The `SEG_LEN_MIN` floor takes over below
`r_floor = SEG_LEN_MIN / SEG_LEN_K`, i.e. a DIAMETER floor of
`2 * SEG_LEN_MIN / SEG_LEN_K`. At this project's typical `SEG_LEN_K = 0.5`,
that floor reaches much higher up the tree than the parameter's own name/
comment ("shortest segment for thin twigs") suggests:
`SEG_LEN_MIN = 50mm` -> floor applies below 20cm diameter;
`SEG_LEN_MIN = 200mm` -> floor applies below 80cm diameter - i.e. almost
the ENTIRE tree, since DBH here is only ~35-40cm. Cylinders inside the
floor zone are close to `SEG_LEN_MIN` long (not exactly - a branch point or
tip always ends a cylinder early, so most floor-zone cylinders are near
`SEG_LEN_MIN` but some are shorter, never longer).

**Practical consequence discovered this session:** varying `SEG_LEN_MIN`
alone (holding `radius_threshold_mm`/`seg_k_pct` fixed) caused a large,
systematic shift in `total_m3` for the `>=10cm`-filtered AdTree-calibrated
rows - confirmed via `parameter_sensitivity.py`'s facet grid (one panel per
`adqsm_variant`, x=`radius_threshold_mm`, coloured by `seg_min_mm`): the
`seg_min_mm` colour groups were clearly separated by roughly 10-15
percentage points, while `radius_threshold_mm` itself had almost no effect
within any one panel.

**Two candidate mechanisms discussed, NEITHER confirmed yet** (tracked as
item 1 in `TODO_investigations.md`):
1. Per-cylinder radius averaging over a coarser `SEG_LEN_MIN` could push
   some cylinders straddling the 10cm diameter cut-off across the boundary
   either way - not a one-directional bias by itself, so unlikely to fully
   explain a broad systematic shift on its own.
2. Changing `SEG_LEN_MIN` changes the raw AdTree radius distribution per
   branch order, which feeds `build_quantile_matched_pairs()`/
   `fit_radius_regression()` for `[calmethod=regression-perorder]` - a
   different `SEG_LEN_MIN` could therefore refit different `(a, b)`
   coefficients per order, rescaling EVERY cylinder of that order, not just
   ones near the 10cm boundary.

Mechanism 2 is the leading hypothesis, since the observed effect is a
broad, systematic shift rather than boundary noise - but this is
**UNVERIFIED**.

**Status: open, unverified.** See `TODO_investigations.md` (item 1) for the
concrete next step (compare fitted `(a, b)` coefficients across two
`SEG_LEN_MIN` values directly).

---

## Step 12 - Confirmed Step 11's mechanism 2; added optional length-weighted radius statistics

**Files:** `tree_geom_utils.py`, `adtree_reconstruct_compare.py`

**Observation:** across a 9-combination `SEG_LEN_MIN`/`SEG_LEN_K` sweep, raw
AdTree volume varied by only **0.34%** between the extremes, while
`[calmethod=regression-perorder]` calibrated volume varied by **33.6%** over
the same sweep - confirming Step 11's open item: `SEG_LEN_MIN`/`SEG_LEN_K` (a
purely representational resampling-density choice) was changing the
calibrated answer even though the underlying tree did not change.

**Cause, confirming Step 11's mechanism 2:** `build_quantile_matched_pairs()`
and `compute_order_calibration_factors()` both compute their per-order radius
statistic (a quantile-matched value, a median) treating each AdTree cylinder
as one sample regardless of its length. The adaptive resampling target
length is proportional to radius in the linear region (`L = SEG_LEN_K * r`,
so cylinders per metre of branch scales as `1/r`) and clamped to a constant
floor (`L = SEG_LEN_MIN`) below a break radius - so changing
`SEG_LEN_MIN`/`SEG_LEN_K` changes how densely different radius ranges get
sampled into cylinders, reweighting the per-order radius distribution the
statistic sees, without changing the branch it describes.

**Evidence:** between two `SEG_LEN` cases whose per-order total cylinder
length and volume agreed to within 1% (confirming the geometry itself was
unchanged), order 1's median AdTree radius shifted by **39%** unweighted vs.
**0.15%** length-weighted - isolating the sampling-density artifact from any
real geometric difference.

**Fix:** both functions gained an optional `weighting="none"|"length"`
keyword argument (default `"none"`, BYTE-IDENTICAL to before this change),
plus a new `weighted_quantile()` helper (plain cumulative-weight
interpolation, no new dependency - ported from the diagnostic script that
first demonstrated the effect, `scratch_calib_sensitivity.py`).
`adtree_reconstruct_compare.py`'s `RADIUS_STAT_WEIGHTING_LIST = ["none",
"length"]` computes BOTH variants in every run (looped around the fit
computation/application, NOT around `convert()` - geometry does not depend
on weighting), each writing to its own row/filename via a `"-wlen"` suffix
tag so neither overwrites the other.

**Result:** over the same two extreme `SEG_LEN` cases,
`[calmethod=regression-perorder]` spread fell from **29.5% to 0.61%** total
(branch volume alone: **50.4% to 0.95%**); `[calref=min5mm]` fell from
**15.4% to 0.17%**. The raw AdTree control (unaffected by calibration) sits
at **0.10%** - length weighting brings both calibration methods to the same
order of magnitude as that control, down from 30-150x larger.

**Limitation, stated plainly:** weighting does NOT improve every order.
Comparing fitted per-order regression coefficients between the same two
cases, the ratio moves toward 1 for orders 1-5 under weighting, but AWAY from
1 for orders 6-7 (order 7: 1.053 unweighted -> 1.105 weighted). Those two
orders carry about 3.2% of branch volume, so the aggregate result above is
unaffected, but the mechanism is not uniform across orders.

**Incidental finding:** branch orders 8 and above receive NO calibration at
all, from either method - this AdQSM variant's `BranchStructure.txt` has no
rows for those orders, so `build_quantile_matched_pairs()` has nothing to
match them against; those cylinders keep their raw AdTree radii via the
existing "leave unscaled" fallback in `apply_radius_regression_per_order()`.
They are **0.35%** of total volume but **22.97%** of total branch length.
Whether AdQSM and AdTree define branch order equivalently is an **OPEN
QUESTION** - not investigated here, recorded for future follow-up.

---

## Open items

- The remaining gap between the primary calibration's `>=10cm`-filtered
  branch volume and the destructive field reference is attributed to
  reconstruction completeness (AdTree/AdQSM not resolving every branch
  the field measurement captured), not a calibration-method problem -
  tracked as a separate open item, not part of this investigation.
- See `TODO_investigations.md` for this and every other currently open
  investigation across the whole project (not just this AdTree/AdQSM
  calibration scope) - added as a consolidated cross-project list so open
  items stop being scattered across handoffs and chat history.
