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

## Result: all 150 component checks pass

- **At the planned size (100 × 20):** the worst bias is 5.6% of the truth, 56% of
  the tolerance. It is in Var_I[E(μ | θ_I)] under the θ-dependent mask with the
  held parameters dominant (per_theta_I: 0.0643 against 0.0681). It is real, at
  4.3 MC SE: the O(1/K) ratio bias the estimator's docstring names.
- **At the small size (30 × 10):** the worst is 7.0% of the truth.
- **Without a mask:** the worst is 2.4% of the truth (24% of its tolerance),
  which is MC noise, as expected of an unbiased estimator.
- **The offset shape** passes like the balanced one, so the centring holds at
  10⁵.
- **Negative components** (reported as zero, flagged) occur in up to 19% of
  replicates at the small size. They are mainly the small Var_I term when the
  held parameters dominate.

**Consequence:** the H-first order passes under the θ-dependent mask, so note
§5.2's fallback (I-first only) is not needed.
