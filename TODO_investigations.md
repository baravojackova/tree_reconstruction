# Open investigations

Consolidated list of every currently open/pending investigation across the
project, so they stop being scattered across handoffs and chat history.
Each item: a short title, current status, and next concrete step. This
list reflects what's known as of this session (2026-09-04) - it is not
necessarily exhaustive, and nothing here is marked resolved; resolving any
of these is separate future work. When an item IS resolved, move it out of
this file (into the relevant changelog) rather than just deleting it.

---

## 1. Regression coefficients vs. SEG_LEN_MIN

**Status:** hypothesis proposed, not yet checked. See
`CHANGELOG_adtree.md`, Step 11, for the full finding this comes from.

**Next step:** run `adtree_reconstruct_compare.py` for at least two
`SEG_LEN_MIN` values (e.g. 50mm and 200mm), with the fitted `(a, b)` per
branch order logged/printed by `fit_radius_regression()`/
`apply_radius_regression_per_order()` if not already printed, then compare
the two coefficient sets directly. If they differ substantially, this
confirms mechanism (b) from Step 11 (a different `SEG_LEN_MIN` refits
different per-order coefficients, rescaling every cylinder of that order);
if the two sets are nearly identical, look harder at mechanism (a)
(cylinders straddling the 10cm boundary) or a third explanation not yet
considered.

---

## 2. Frustum-accurate cylinder volume

**Status:** deferred, untouched. Fully documented in `CHANGELOG_adtree.md`,
Step 10, including the quantified impact (~3-4% max) and why it was
deferred rather than fixed immediately.

**Next step:** requires coordinated changes to the ANSYS export format
(`write_geom()` / the `bk1.mac` macro) and confirmation of equivalent
accuracy from TreeQSM's own cylinder model before changing AdTree's -
otherwise "fixing" only AdTree's number introduces a new, asymmetric bias
across methods.

---

## 3. TreeGraph - B21_S01 native crash

**Status:** unresolved. `identify_stem()` crashes with `0xC06D007F`
(Windows MKL/BLAS/OpenMP conflict class). A float32 precision bug in the
`.xyz -> .ply` conversion was already fixed (recentering applied), but the
crash persists.

**Next step:** not yet identified - likely needs isolating which specific
MKL/BLAS call triggers it (env var toggles, or a minimal repro script).

---

## 4. TreeGraph - crown reconstruction quality

**Status:** unresolved. 3.15% of cylinders on IND01_054 above ~16m violate
the 2x parent/child radius safeguard.

**Next step:** cause (a data limit vs. an algorithmic issue in
`generate_path()` or `build_skeleton.py`) not yet investigated.

---

## 5. COMPUTE_CALREF_MIN5MM toggle

**Status:** UNCERTAIN - re-check before acting. A handoff prompt was
prepared earlier, but implementation status was never confirmed in this
session.

**Next step:** check `adtree_reconstruct_compare.py`'s current
`COMPUTE_CALREF_MIN5MM` value and whether `[calref=min5mm]` rows exist in
the current `volume_results.csv` before assuming either way.

---

## 6. Chart group ordering / reference-family colour refinement for "none" mode

**Status:** handoff prompt prepared, not yet implemented. No further detail
available from prior sessions - revisit the original handoff note if this
becomes a priority.

---

## 7. Interactive (Plotly/HTML) version of parameter_sensitivity.py charts

**Status:** discussed, deliberately NOT pursued yet. Bara decided the
existing static facet grid already answers "which combination is min/max
for a given metric" well enough via axis position, and an interactive
version is a bigger investment (new dependency, `.html` output instead of
`.png`) only worth it if this kind of exploration becomes a repeated need.

**Next step:** none planned - revisit only if the static charts stop being
enough.

---

## 8. IND07_083/04 taper.txt - implausibly flat trunk diameter section

**Status:** confirmed AdQSM export data-quality issue for this specific
tree/variant, not a pipeline bug. `data/IND07_083/04/taper.txt` reports an
almost exactly constant diameter (~0.708 m) from height 0 m to 22.6 m -
real trunk taper should decrease measurably over 22+ metres, so this is
not physically plausible - then becomes erratic (0.49/0.71/0.55/0.62/0.63
m) from 23.6-27.6 m, beneath the already-handled 28.6 m spike (6.3952 m
diameter, correctly dropped by `_reject_taper_spikes()`). This inflates
that tree's calibrated trunk volume (~16.47 m3) well past AdQSM's own
official TrunkVolume (10.535 m3, `TreesParams.txt`).

Not a pipeline bug: trunk (order 0) radius calibration reads ONLY AdQSM's
taper.txt (never AdTree's own measured radius) and is calmethod-
independent - confirmed by `calmethod=min5mm` and
`calmethod=regression-perorder` producing byte-identical trunk numbers for
this tree. It also does NOT affect the `calmethod=regression-perorder`
recommendation (`AdTree_Calmethod_Decision.docx`) - the effect is
identical for both candidate methods, so it favours neither one.

**Next step:** check the other reference/production trees' `taper.txt`
files for similarly flat or erratic multi-row sections before trusting
their calibrated trunk volumes. Also consider whether
`_reject_taper_spikes()` (`tree_geom_utils.py`) should be extended to also
flag suspiciously flat or erratic multi-row SPANS, not just isolated
single-row spikes - not implemented yet, just recorded as a future option
here.

---

## 9. AdTree calibrated trunk under-performance: confirmed on IND07_083 only, not yet checked on other trees

**Status:** partially explained, not yet generalized. AdTree's trunk
calibration replaces trunk radius entirely with AdQSM's own taper.txt
curve (see `apply_order_calibration_factors()`/
`apply_radius_regression_per_order()` in `tree_geom_utils.py`) - it uses
none of AdTree's own measurement. For IND07_083, that taper.txt curve was
shown to be implausibly flat (~0.708m diameter, 0-22.6m height) - a
data-quality issue in AdQSM's own export, not a pipeline bug (see
`AdTree_Smoothing_Resampling_Summary.docx` Section 4 and the
`taper_curve_compare.py` diagnostic script).

Open question: has this same taper.txt flatness/implausibility been
checked on IND01_054 and IND03_088? Not yet done. Until checked, it's
unknown whether AdTree calibrated's trunk/total RMSE disadvantage vs.
TreeQSM manual (see `TreeQSM_Reconstruction_Decision_Summary.docx`) is a
universal AdTree limitation or specific to trees with a flawed taper.txt.

**Next step:** run `taper_curve_compare.py` for all 3 reference trees,
visually/numerically check for flat or erratic multi-row spans the same
way IND07_083 was checked.

---

## 10. Best-of-tree vs. realistic AdTree picks - always report both

**Status:** established convention, apply going forward.
`reconstruction_method_decision_summary.py`'s "best-of-tree" AdTree
calibrated/AdQSM-other-variants groups select the single configuration
matching the reference best PER TREE - this uses knowledge of the answer
and is an oracle/upper-bound, not a configuration achievable on a tree
with no known reference (like the 12 production beech trees).

The "realistic" variant (radius_threshold_mm fixed at 5mm, adqsm_variant
chosen via DBH-matching only) is the fair, achievable comparison point -
it is consistently ~2.3-2.4x worse than the oracle for total_m3 across
all 3 comparison blocks.

**Next step:** rule going forward - any time a "best-of"/oracle-style
number is quoted from this pipeline, quote the realistic number alongside
it. Do not present best-of-tree numbers alone as if they were achievable
in production.

---

## 11. Hybrid TreeQSM-trunk + AdTree-branch export - not yet built

**Status:** open idea, not implemented. TreeQSM manual wins clearly on
trunk/total volume (driven by the taper.txt-dependency issue in item 9
above); AdTree calibrated wins clearly on branch volume specifically (its
per-order regression calibration against AdQSM's branch data is robust
even in the "realistic", non-cherry-picked configuration).

**Next step:** not started. Combine TreeQSM's trunk geometry with
AdTree's calibrated branch geometry into one hybrid export - could in
principle beat either method alone on total volume. Would need a new
export path merging two different reconstruction pipelines' outputs - not
attempted, no prototype exists.

---

## 12. Beech-appropriate TreeQSM manual parameters not yet determined

**Status:** open, not started. TreeQSM manual's pd1=0.08/pd2min=0.02/
pd2max=0.07 are published values specific to the IND reference-tree
group's forest/species/scan conditions (confirmed with Bára - not tuned
to these 3 trees, so no circularity in the reference-tree comparison) -
but NOT directly transferable to the beech production trees (different
forest/group).

**Next step:** find beech-appropriate published values, or run a small
manual-parameter testing pass specifically on the beech trees (using
their own structural reference - frequency/tilt/deflection - once that
phase begins), before defaulting to TreeQSM manual for production
reconstructions.

---

## 13. Branch-order correspondence between AdQSM and AdTree

**Status:** unresolved, not blocking. A fixed order shift is ruled out (see
`CHANGELOG_adtree.md`, Step 13), but the two tools' per-order diameter
ratios differ twofold between trees and the residual monotone per-order
trend after length weighting is unexplained. Not blocking: with length
weighting the calibrated volume spread across `SEG_LEN` settings is 0.61%
against a 0.10% raw control, orders 8+ are already uncalibrated by
construction, and the trunk is replaced wholesale by the AdQSM taper curve
so per-order regression does not affect it.

**Next step:** revisit when all 12 beech trees are reconstructed, and
specifically check the calibrated volumes of the three IND trees against
their destructive reference - that is the only available control that
would reveal a systematic bias, as opposed to the instability that length
weighting fixed.

---

## 14. bk_kp.mac - Poisson ratio block (Zlamal et al. material model)

**Status:** found during the PC/Mac merge review (2026-11-25), not yet
fixed. The new orthotropic material block in `scripts/ansys/bk_kp.mac`
(Zlamal et al., added after commit `1305fad`) has two issues in the
"conditions for proper setting of poisson numbers" section (around
lines 266-300):

1. **Syntax error (double `=`):** `pr_xy = = pr_rt*e_r/e_t` in the
   `*elseif,e_t,gt,e_r` branch. For beech (E_R = 1245 MPa > E_T = 474 MPa)
   this branch is never executed, so it does not crash today - but APDL
   would fail on it if E_T > E_R ever occurs (other species, other data).
2. **Computed values are overwritten:** the `*if` block computes
   `pr_xz`, `pr_yz`, `pr_xy` from the E-ratios, but a few lines later
   `pr_xy=pr_lr`, `pr_yz=pr_rt`, `pr_xz=pr_lt` replace them with the raw
   input values. The `*if` result is therefore never used. Also check the
   direction mapping in that block (comments say "z gt x" next to
   `e_l gt e_r`, while e_x=e_l, e_y=e_r, e_z=e_t).

**Next step:** decide which Poisson ratios (major/minor, nu_ij vs nu_ji)
ANSYS `mp,prxy/pryz/prxz` should receive for the L=X, R=Y, T=Z mapping,
then either fix the `*if` block and delete the overwrite, or delete the
`*if` block. Rerun one tree (modal) to confirm frequencies before and
after.

---
