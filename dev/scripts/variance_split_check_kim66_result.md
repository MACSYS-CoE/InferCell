# KIM-66: the variance-split estimator, checked on known splits

**Run:** `dev/scripts/variance_split_check_kim66.sh` at `a742224` (pushed),
job 19560131, into `/fred/oz022/tkimpson/infercell-m4-kim66-2026-10-30/`
(`SHA256SUMS`). R = 2,000 replicates of each configuration.

**Design** (the script's header, declared before the run):
- three shapes of μ: balanced; held-dominant; and balanced moved by 10⁵;
- three masks: none, independent (π = 0.7), and θ-dependent
  (π = logistic(1 + x + y));
- both conditioned weight sets under a mask;
- two sizes: (M, K) = (100, 20) and (30, 10);
- the truth under each law computed on a grid.

**Pass:** each of the five components is within ±10% of the truth, the 10%
taken of max(|truth|, 0.05 × total).
- **That floor was added in the script, not in the approved plan.** All 150 rows
  also pass at 10% of |truth|.
- **The archive's `SHA256SUMS`** has SHA256 `7e3fa6fc104e12a7…`.

## Result: all 150 component checks pass

- **At the planned size (100 × 20):** the worst bias is 5.6% of the truth, 56% of
  the tolerance. It is in Var_I[E(μ | θ_I)] under the θ-dependent mask with the
  held parameters dominant (per_theta_I: 0.0643 against 0.0681). It is real, at
  4.3 MC SE.
  - It is consistent with the O(1/K) ratio bias the estimator's docstring names.
  - It is also consistent with the shortfall from rows with fewer than two viable
    columns.
  - From (30, 10) to (100, 20) it shrinks from 7.0% to 5.6%, less than halving.
    So the cause is not pinned down.
- **At the small size (30 × 10):** the worst is 7.0% of the truth.
- **Without a mask:** the worst is 2.4% of the truth (24% of its tolerance),
  which is MC noise, as expected of an unbiased estimator.
- **The offset shape** passes like the balanced one, so the centring holds at
  10⁵.
- **Some component came out negative** in up to 19% of replicates at the small
  size; it is reported as zero and flagged. The TSV records only whether any
  component was flagged, not which one.

**Consequence:** the H-first order passes under the θ-dependent mask, so the
fallback KIM-66 declared (I-first only) is not needed.

**What it does not show:** precision.
- At 100 × 20, one run's components vary by about 25–40%, and up to 64% for a
  small Var_I term (from `mc_se × √2000`). That is mostly the √(2/(K−1)) noise
  of 20 θ_H draws.
- The held share, a ratio, largely cancels it.
- The smallest component tested is 5–8% of the total. KIM-35 puts θ_I's share
  near 3%.
