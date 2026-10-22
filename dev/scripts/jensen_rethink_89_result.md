# Issue #89, after B3: what makes the GTP gap, and whether an offset fixes it

**Rules:** fixed before the runs in `spec/phases/approximate-method.md` §12
(2026-10-22). §8 KG is unchanged and gives the verdicts of record. The KG null and
its quantiles are information only (§12 2026-10-23).
**Script:** `dev/scripts/jensen_rethink_89.jl`, with `.slurm`, `.sh` (stages 1 and
2A) and `_transfer.sh` (2A on independent seeds).
**Runs:** every job was submitted from a detached worktree at a pushed commit.

| Jobs | Commit | What | CPU-h |
|---|---|---|---|
| 19012789–19012801 | `bf54df3` | Smoke: 4 cells, 120 s, every stage (the score failed on a `Dict` key type, fixed in `29781f6`) | 1.11 |
| 19013270–19013346 | `29781f6` | Stage 1 arms; 2A ensembles at seeds 50001…; 31 mean fields; score | 21.42 |
| 19014973–19014975 | `bc71511` | Seed-range check: f10_ref, m0_low and m0_high at seeds 50001… | 1.21 |
| 19015472–19015544 | `b2f6b4c` | 2A on independent seeds, and twenty posterior-draw mean fields and ensembles. 76 tasks failed to load a precompile cache built for another CPU; 5 ran past the dependent stages | 18.50 |
| 19017428, 19017558–19017599 | `b2f6b4c` (fresh worktree) | Precompile on a milan node; the missing pieces rerun on milan only; then transfer (19017597), calibrate (19017598) and score (19017599) | 5.75 |
| **Total** | | | **48.0** of the 100 allowed |

**Raw outputs:** `/fred/oz022/tkimpson/infercell-89-rethink-archive-2026-10-22/`.
It holds `scores.md`, `transfer.md` and `calibrate.md`, every ensemble and mean
field, and the logs.
**Regenerate:**
- `dev/scripts/jensen_rethink_89.sh <outdir>`;
- then `dev/scripts/jensen_rethink_89_transfer.sh <outdir>`;
- both from a worktree at the commits above. Run on milan nodes only (below).

## Verdicts

1. **The cause is cell-to-cell variation in the drains.** Give every cell the mean
   field's drains, and the 200-cell mean *is* the mean field: KG totals of 70.8,
   12.3 and 25.9, against the mean field's 70.7, 12.3 and 25.9. Enzyme counts,
   the per-cell rebuild feedback and the clip are ruled out by the pre-registered
   rules. Shot noise is not ruled out: it accounts for up to 0.44 of the offset at
   m0_low.
2. **Under KG as agreed, no description passes 2A's rule** (≥ 7 of 8 targets, ≥ 9 of
   10). A fixed offset δ from the reference point passes KG at 5 of 8 F10 posterior
   draws, 4 of 10 targets around m0_low and 6 of 10 around m0_high. The
   uncorrected mean field passes at 0 of 28.
3. **KG as agreed fails a description known to be right about half the time.** Take
   a 600-cell mean at the same θ, scored as a mean field against 200 other cells: it
   passes KG at 46%, 36% and 56% of 200 random partitions (f10_ref, m0_low,
   m0_high).
   - So 2A's rule could not be met by a perfect description. Its expected pass count
     is about 3.7 of 8 at F10, and 3.6 and 5.6 of 10 at m0_low and m0_high.
   - δ's pass rates (63%, 40%, 60%) sit at the null's.
4. **Placed in that null (information only), δ behaves like a right description.**
   Both of KG's statistics fall within the null's 95% at 26 of 28 independent
   targets, against 5 of 28 for the uncorrected mean field.
   - **What a pass in that null means:** the residual gap is below what a 200-cell
     reference can detect. The null's power grows with the reference's size, so a
     larger reference would sharpen the test and could fail δ.
   - **The test does have power at this size:** it fails the uncorrected mean field
     at 23 of 28 targets.

**Recommendation.**
- **Nothing is certified under the gate as agreed.** The gate must be recalibrated
  before any description can pass it. The proposal is in the sub-spec's §12, for
  your approval; it is not applied here.
- **If you adopt a calibrated gate,** the fixed offset is the second description on
  M0, as a declared reduction.
  - It is measured once per dataset, from the same ensemble that gives S.
  - It cost under 1 CPU-h of ensemble here (800 M0 cells, about 0.5 s each after
    compilation).
- **Before B4,** the next runs are listed at the end.

## Stage 0: the impulse against the pool (B3's archive)

- **Pool against impulse.** M0's GTP pool is 14k–22k particles per cell (mean over
  saves and cells). One GAPD completion debits about 670, which is 0.03–0.05 of
  the pool.
- **Yet the spread is large.** GTP's per-cell CV grows from 0 at t = 0 to about 0.5
  by 600 s, and 3PG's exceeds 1.
- **What changed in the plan:** the variation is slow, not shot noise, so Stage 1's
  Poisson-completion arm was replaced by the four arms below (§12).

## Stage 1: which variation makes the gap (M0, 200 cells per arm)

**Check:**
- The full arm reproduces B3's archived cells exactly: 200 of 200 cells per truth
  matched by seed, worst relative difference 0.
- The clip is not a cause: 2 drains clipped out of 360,000, both GTP_translat in
  one cell at m0_high.

**Share** is the arm's GTP offset (the mean over saves of log x̄_arm/x̄_full) over the
mean field's.
- Share ≈ 1: the arm removed the variation that causes the gap.
- Share ≈ 0: what it removed is irrelevant.

| Arm: what it removes | m0_v6 | m0_low | m0_high | Rule (§12) |
|---|---|---|---|---|
| `drains_mf`: all drain variation | **1.00** | **1.00** | **1.00** | ≥ 0.7 at all three: **named cause** |
| `drains_own`: drain shot noise only | 0.08 | 0.44 | 0.13 | Neither named nor ruled out (m0_low is 0.44) |
| `enzymes_mf`: catalytic and volume counts | 0.14 | −0.12 | −0.29 | < 0.3 at all three: ruled out for GTP |
| `popmean`: per-cell rebuild feedback, read live | −0.01 | 0.18 | −0.05 | < 0.3 at all three: ruled out |

The mean field's GTP offsets are −0.172, −0.107 and −0.107.

**`drains_mf` reproduces the mean field to three figures:**

| Truth | Total | Noise |
|---|---|---|
| m0_v6 | 70.8, against the mean field's 70.7 | 25.0 |
| m0_low | 12.3 (12.3) | 10.3 |
| m0_high | 25.9 (25.9) | 21.3 |

With the drains fixed, nothing else in M0 produces cell-to-cell variation in the
metabolism that the mean misses.

**The mechanism**, from the full arm's saved gene expression (200 cells, a
login-node read of the archive):
- A cell's GTP at 600 s falls with the number of proteins it made. The correlation
  with ptsG protein gain is −0.78 at m0_v6, −0.18 at m0_low and −0.48 at m0_high.
- Protein output varies about 2-fold across cells (gain CV 0.3–0.8), so GTP spend
  does too.
- *Interpretation, not measured:* GTP's only regeneration is PGK3 and PYK3, which
  are saturable. If GTP is convex in the spend, the population mean sits above GTP
  at the mean spend. That is the Jensen gap. What is measured is that the drains
  cause the gap (`drains_mf`), and that per-cell GTP falls with protein output.
- The rest of the imbalance follows from the same spend. GDP is GTP's partner, and
  ATP, AMP and GMP couple to it through GK1 and ADK1.

**What the arms do not separate.**
- `enzymes_mf` moves pools other than GTP at m0_v6: 3PG, 2PG, Pi and lactate shift
  by 2–3× their noise, while GTP's offset is small.
- Fixing enzymes at the mean also breaks their covariance with the drains, since
  both follow the same mRNA path. The decomposition is not additive.
- Any second-order correction must carry the drain-rate variance *and* its
  covariance with the enzyme counts.

## Stage 2A: an offset at a posterior's width

**What δ is.** δ = log x̄_ref − log x̄_MF(ref), per pool and save, from 800 full cells
at the reference point (seeds 40001…). It is added to the mean field at each target.
Each target is scored by KG against its own 200 full cells.

**The first targets shared seeds, which made them a single test.**
- At `29781f6` every target ran at seeds 50001….
- f10_ref itself, at those seeds, gives a GTP offset of +0.108, against +0.147 from
  its 800 cells at 40001…. The four 200-cell groups there range from +0.119 to
  +0.169.
- So the eight F10 targets' common offset (+0.08 to +0.13) was the noise of one
  shared seed range, not a θ-dependence of δ.
- That run's `scores.md` (δ_ref passes at 2 of 8, 9 of 10 and 4 of 10) and
  `calibrate.md`'s "2A's targets in the null" are **superseded**. The
  leave-one-out δ̄ there also leaks the target's own noise through the shared
  seeds.
- **The rerun:** each target got seeds of its own (`transfer.md`).
- **Order of events:** the shared-seed verdicts above were read first. Diagnosing
  them exposed the fault, and the seed check (19014973–5) confirmed it. The rerun's
  design, δ̄ included, was then fixed before its verdicts were read.

**Results on independent seeds** (`transfer.md`):

| Source of δ | Targets | KG passes: uncorrected | KG passes: δ_ref | KG passes: δ̄ (posterior spread) | Within the null's 95%: uncorrected / δ_ref / δ̄ |
|---|---|---|---|---|---|
| f10_ref (F10's posterior mean, V6's dataset) | 8 posterior draws at Mahalanobis quantiles | 0 | 5 | 4 | 3 / 7 / 7 |
| m0_low | ±1 posterior SD in one target at a time | 0 | 4 | – | 2 / 9 / – |
| m0_high | the same | 0 | 6 | – | 0 / 10 / – |
| **Pre-registered rule** | | | **≥ 7 of 8, ≥ 9 of 10: not met** | | |

- **δ̄** is spread over the posterior: 40 cells at each of twenty random F10 draws,
  as log mean_k[x̄_k / x_MF(θ_k)]. It does no better than δ_ref, consistent with δ
  barely varying across the posterior.
- **The worst target** is f10_d4: total 55 against noise 22.8, with FDP, external
  lactate and ADP above 3× their noise. Its total sits at the null's 0.96 quantile.
- **The pools that fail under δ** are mostly ones whose split-half noise is tiny:
  FDP at about 1e-4 nats, external lactate and ADP. `transfer.md` does not report
  the 0.1-nat floor variant; `scores.md`'s does, but only on the superseded
  shared-seed targets.

## KG's null (information only)

**How it is built.** At each 800-cell point, 200 random partitions score the
600-cell mean, in place of the mean field, against the other 200 cells, exactly as
KG scores a mean field (`calibrate.md`, job 19017598, `Xoshiro(89)`).

| Point | KG passes | total/noise: median (90%) | max pool shift/noise: median (90%) |
|---|---|---|---|
| f10_ref | 93 of 200 (46%) | 0.65 (1.55) | 2.5 (14.9) |
| m0_low | 72 of 200 (36%) | 0.64 (1.37) | 3.5 (12.1) |
| m0_high | 112 of 200 (56%) | 0.69 (1.37) | 2.0 (7.0) |

- **The total clause is about right.** A right description's total/noise has a
  median of about 0.65.
- **The per-pool clause is what fails.** The largest pool's shift/noise has a
  median of 2–3.5 and a 90% point of 7–15.
  - Each pool's noise is one split-half estimate from about ten correlated saves.
    It has few degrees of freedom, so "≤ 3× noise in every one of 17 pools" fails
    a right answer often.
  - FDP, external lactate and ADP, whose noise is ~1e-4 to 1e-3 nats, trip it
    most.
- **This null is slightly generous to the mean field.** Its M is a 600-cell mean,
  whose own variance is a third of the target's. A deterministic mean field with
  δ from 800 cells carries a quarter.

## Not done, and why

- **Stage 2B (LNA or second order).** Stage 0's impulse/pool ratio (0.03–0.05)
  allowed it under the rule, but it was not built. Stage 1 shows what it must carry:
  - the drain-rate variance;
  - its covariance with the enzyme counts;
  - its propagation through pools with a per-cell CV above 1 (3PG, 2PG, Pi), where
    a Gaussian closure is doubtful.

  The full arm's saved gene-expression rows are its validation data.
- **Stage 3 (the reduced model).** Its rule was "only for whatever passed on M0".
  Nothing passed KG's rule, so it was not run. It is affordable in about 50 CPU-h
  once a gate is agreed.
- **The reference θ is F10's posterior mean.** That is known here from the exact
  sampler, not from an approximate fit. The realistic two-pass pipeline (fit with
  the uncorrected mean field, measure δ there, refit) is untested. B7 (AV6) is
  where it would be tested.

## Next runs this points to

1. **Agree a calibrated gate.** §12 has the proposal. Then re-read this file's
   tables under it; that needs no new compute.
2. **The reduced model** at A3's truths, with the same independent-seed design: δ
   from 400–800 cells, 4–6 targets at the information bound's width. About 50 CPU-h.
3. **The two-pass pipeline on M0:** δ at the uncorrected mean field's mode, not at
   F10's mean. It needs B5's likelihood, or a cheap optimiser over the mean field
   with a fixed S.

## A cluster note for anyone rerunning

- **Symptom:** a fresh worktree's first jobs, spread over node types (dave milan,
  gina, turin), wrote a precompile cache that milan nodes then refused: "Unable to
  find compatible target in cached code image", at `using InferCell`.
- **Fix:** precompile once in a single job on a milan node (job 19017428), and keep
  the array on `--partition=milan --exclude='dave[148-149]'`.
