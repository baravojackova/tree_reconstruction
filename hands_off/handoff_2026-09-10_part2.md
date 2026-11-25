# Handoff — 2026-09-10, part 2

Supplement to `handoff_2026-09-10.md`. That file describes the evaluation phase and
the original run plan. **This one supersedes its section 5 (the run plan) and adds a
newly found data defect.**

---

## 1. Final model selection — decided

**One AdTree model and one TreeQSM model per tree.**

```
AdTree   r5mm, regression-perorder, -wlen, selected AdQSM variant, seg50-500-k80
TreeQSM  pd2min = 0.002, sr010, ri1, Simplified (no islands)
```

File name patterns:

```
geom_<TREE>_r5mm_adqsm<VAR>_seg50-500-k80-wlen.txt
scripts/matlab/archive/<date>_<TREE>/simplified_<tree>_man_pd<PD1>-002-<PD2MAX>_mo9_sr010_ri1.mat
```

Per-tree parameters, ordered by field DBH ascending:

| tree | field DBH | AdQSM variant | AdTree r5mm cyl | TreeQSM tag | TreeQSM cyl | total | time |
|---|---|---|---|---|---|---|---|
| S12 | 30 cm | 100 | 4 173 | `man_pd04-002-05` | 4 710 | 8 883 | 16.9 min |
| S16 | 33 | 090 | 4 072 | `man_pd045-002-05` | 2 611 | 6 683 | 12.7 |
| S09 | 36 | 060 | 10 332 | `man_pd05-002-08` | 11 759 | 22 091 | 42.1 |
| S10 | 37 | 040 | 27 624 | `man_pd05-002-09` | 9 765 | 37 389 | 71.3 |
| S08 | 39 | 060 | 9 634 | `man_pd05-002-09` | 8 413 | 18 047 | 34.4 |
| S14 | 41 | 080 | 15 963 | `man_pd05-002-09` | 14 468 | 30 431 | 58.0 |
| S04 | 43 | 090 | 13 956 | `man_pd05-002-10` | 10 212 | 24 168 | 46.1 |
| S11 | 46 | 070 | 16 389 | `man_pd06-002-11` | 12 688 | 29 077 | 55.4 |
| S13 | 49 | 070 | 12 233 | `man_pd06-002-10` | 9 711 | 21 944 | 41.8 |
| S21 | 49 | 070 | 29 338 | `man_pd066-002-099` | 15 435 | 44 773 | 85.4 |
| S01 | 56 | 080 | 45 461 | `man_pd07-002-10` | 21 864 | 67 325 | 128.4 |
| S07 | 64 | 080 | 35 387 | `man_pd08-002-14` | 19 220 | 54 607 | 104.1 |

**Nine trees fit in 7 hours** (S12 through S13, cumulative 6.31 h). The three largest —
S21, S01, S07 — are deferred, exactly as before, and for the same reason: the selection
is by computational cost, not by result, and that must be stated in the text. All
twelve would be 11.6 h.

**The TreeQSM tag is not unique across trees** — `pd05-002-09` belongs to S08, S10 and
S14 alike. The tree is identified by the archive folder, never by the filename. Put the
tree name in every exported filename.

### Why r5mm and not r10mm

AdTree's radius threshold **discards** cylinders thinner than the threshold; TreeQSM's
simplification **merges** cylinders along a branch and discards only ~13 % of length.
The two are not comparable operations. Measured effect on AdTree branch length:

```
r5mm  100 %      r10mm  48 %      r20mm  23 %
```

r10 mm would enter ANSYS missing half of AdTree's branch skeleton, purely because of a
threshold we chose. Ratio of TreeQSM to AdTree branch length:

```
TreeQSM vs AdTree r10mm   min 2.61   median 4.55   max 6.99
TreeQSM vs AdTree r5mm    min 1.43   median 2.49   max 3.92
```

Cost of r5mm over r10mm: about +1.8 h across eight trees. **Even at r5mm the models are
not equivalent** — a factor 2.5 in branch length remains, and branch *volume* differs
by 1.5× to 17.7× (median 3.4×, extreme on S12). That residual difference is handled by
the density-decomposition experiment, not by further matching.

### Why pd2min = 0.002 — and on what basis

**No geometric criterion separates the pd2min levels.** Deviation of the trunk taper
curve from the field measurement, median across 12 trees:

```
0.002  1.43 %     0.005  1.33 %     0.010  1.34 %     0.020  1.67 %     AUTO  1.59 %
```

The spread between best and worst is **0.34 percentage points**, against 5.05 points
between TreeQSM and the AdQSM taper curve — roughly fifteen times smaller than the
difference between methods. And the trees disagree completely: best level per tree is
0.002 on 2, 0.005 on 4, 0.010 on 2, 0.020 on 2, AUTO on 2. The median favouring 0.005
is not a consensus.

Cylinder count differs by 1.5 %, branch length by 0.1 %, DBH by 0.2 points — all noise.
The only quantity that moves measurably is **branch volume, median +8 % at 0.005**, and
that must not be used to select: the beeches have no volume reference, and the one
crude estimate available (crown mass from an ellipsoid times density) cannot rank
methods that differ by less than a factor of two.

**Decision: pd2min = 0.002**, chosen as the finest available reconstruction when no
criterion distinguishes the alternatives. **Write it up that way** — state that the
field taper does not discriminate between levels (0.34 pp, trees inconsistent), and
that the finest setting was retained in the absence of evidence to depart from it. A
documented convention is defensible; a silent choice presented as a result is not.

### One measurement that confirmed an expectation

`pd2min` was expected not to affect the trunk, since it is the *minimum* patch diameter
while the trunk is covered by patches near `PatchDiam2Max`. Measured across the four
levels:

```
DBH at 1.3 m      range: median  2.2 %,  max  7.3 %
trunk volume      range: median 13.4 %,  max 27.7 %
trunk taper error range: 0.34 pp total
```

So the parameter barely moves the diameter and does not move the accuracy of the trunk
at all. What it does move is trunk *volume*, and that happens higher up, where
segmentation decides how much of a fork belongs to the trunk rather than to a branch.
The two statements are consistent: `pd2min` shifts the trunk/branch boundary, not the
trunk's shape where the field measured it.

---

## 2. Data defect found — `branch_len_m` and `trunk_len_m` on TreeQSM rows

**`volume_results.csv` carries the pre-simplification length on the `Simplified` rows.**

```
branch_len_m   identical across both stages : 602 of 602 pairs (100 %)
trunk_len_m    identical across both stages : 602 of 602 pairs (100 %)
n_cylinders    identical across both stages :   0 of 602
branch_m3      identical across both stages :   0 of 602
```

Volumes and cylinder counts differ correctly between stages. **Both length columns do
not** — the `Simplified` row carries the `Optimal` model's length. On `B21_S07`
pd2min 0.002 sr010 ri1:

```
              n_cylinders   trunk_len_m   branch_len_m   trunk_m3   branch_m3
Optimal            35 495     27.296671    4373.008301   4.202176   10.050113
Simplified         19 220     27.296671    4373.008301   3.908053   10.725138
```

The true Simplified branch length, computed from the archived cylinders, is **3752.90 m**
against the 4373.01 m in the CSV — the ratio is a stable **0.84 to 0.90** across all 24
models checked, median 0.87.

### Fix — no MATLAB re-run required

The Simplified cylinder data is already in the archive:
`simplified_*.mat` → `QSM_final.cylinder.{length, BranchOrder}`. Verified by the fact
that `n_cylinders` and `branch_m3` computed from those files match the CSV exactly.

The repair is Python: walk the archives, sum cylinder lengths by branch order, and
write the two columns back through `upsert_result()`. Two things to settle first:

1. **Are the `Optimal` rows correct?** They may be, or the value may come from a
   `treedata` field belonging to neither stage. Check whether an optimal-stage `.mat`
   exists in the archive to compare against. If it does not, the `Optimal` length is
   unverified and must be labelled as such rather than assumed good.
2. **Fix the source too.** `runsken.m` exports the wrong `treedata` field. Repairing the
   CSV without repairing the export means the defect returns on the next import. This
   is a MATLAB change and is lower priority than the CSV fix, but it should not be
   forgotten.

### What this affects

- Any statement made about TreeQSM branch or trunk **length**. The method comparison in
  `Vyhodnoceni_metod_12_buku.docx` uses lengths only descriptively, but the numbers
  there are the pre-simplification ones.
- The `eta²` values for `branch_len_m` and `trunk_len_m` in `method_eval_eta2.csv`.
- **Not** affected: volumes, DBH, cylinder counts, the taper comparison, the base/DBH
  comparison, and every accuracy conclusion. Those use volume and radius, not length.

---

## 3. Correction to the previous handoff

`handoff_2026-09-10.md` lists "missing `sr010/ri1` exports at `pd2min = 0.020` and AUTO"
as an open item. **That was wrong — remove it.** Those combinations were never intended
to exist. From the reconstruction phase: `pd2min` 0.002, 0.005 and 0.010 each carry the
full 4 × 4 simp grid; `pd2min = 0.020` is a single probe at `sr005/ri0` and AUTO is a
single model per tree at the same setting. Confirmed by counting rows.

Side effect worth knowing: in the per-`pd2min` error figure, the AUTO and 0.020 points
use a different simplification setting from the other three levels, so they are not
strictly comparable to them. It does not change the conclusion, since all five sit
within 0.34 points.

---

## 4. Note on the taper error definition

Sanity check 2 of the per-`pd2min` extension came out with a deviation of up to 0.51
percentage points against the existing summary. **This is not an error.** The two
numbers are different quantities:

- `trunk_taper_error_summary.csv` takes the median TreeQSM diameter across the five
  models at each height, builds one curve, and computes its deviation — the **error of
  the median model**.
- The new per-model breakdown computes each model's deviation and takes the median of
  those — the **median of the errors**.

Median and absolute value do not commute. Neither is wrong, but the existing column
should be labelled "error of the median model" in the table caption, since the two
differ by up to half a point.

All originally verified numbers are unchanged: 1.4368 % vs 6.4930 %, 12 of 12, and the
per-height deviations on `B21_S04`. Sanity check 1 passed byte-identically on all four
original outputs.

---

## 5. Run plan

Nine trees, two models each, in field-DBH order: S12, S16, S09, S10, S08, S14, S04,
S11, S13. **6.31 hours.** Then, if time allows, S21, S01, S07 (a further 5.3 h).

**Start with S16 (12.7 min) and compare the actual time against the estimate before
committing the rest.** The whole cost model rests on one measured run — `B21_S07`
AdTree r5mm, 35 387 cylinders, 67.5 minutes — and assumes linear scaling at
**0.1144 s per cylinder**. That holds for geometry and constraint equations, which are
87 % of the time; it may not hold for the solver.

Cost relations measured on that run, all exactly linear:

```
2.00 keypoints per cylinder      6.00 constraint equations per cylinder
3.08 elements per cylinder       effective element size 1.88 cm (e_sizetree = 0.1 m)
```

### Control pairs

**Simplification, `ri0` vs `ri3` on B21_S16 at pd2min 0.002** — 4 925 and 861 cylinders,
9.4 + 1.6 = 11 minutes. Volume differs by 3 %, cylinder count fivefold. Measures
topology at near-constant mass.

**Crown mass, density scaling** — method in
`Rozklad_frekvence_hustota_vetvi.docx`. `phi = V_branch,AdTree / V_branch,TreeQSM`,
`rho*_br = rho_br * phi`. **Verify first whether the macro can restart from a built
model**; if it can, only the material constant changes and 87 % of the time is avoided.
Branch/trunk separation is column index 10 of `geom_*.txt` (branch order, 0 = trunk).

Note that `phi` for this pairing is extreme on some trees — branch volume ratios run
from 1.5× (S07, S10) to 17.7× (S12). On S12, `phi ≈ 0.06`, i.e. branch density scaled
to 6 % of nominal. That is itself informative and should be reported, not smoothed.

---

## 6. Macro reminder

`bk_kp.mac`: `n_rad` on line 26 and the geometry filename on line 29 **must be changed
together**. A mismatch does not fail immediately — the array is dimensioned to `n_rad`,
filled only as far as the file goes, and the loop then computes `sqrt(0+0+0)` and
divides by it partway through. This was hit once already this session.

The keypoint collision fix is confirmed working at 70 811 keypoints, seven times the old
10 000 limit, geometry checked visually.
