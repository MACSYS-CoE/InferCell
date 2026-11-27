# KIM-178 result: ten replicate truths on the reduced model, by MAP and Laplace

The design was frozen at tag `prereg/kim-178` (2050f5c) before any replicate existed. Every
job of record ran at commit `f66e0ec`. The archive is
`/fred/oz022/tkimpson/infercell-kim178-archive-2026-11-28`, 80 files with a `SHA256SUMS`.
`sha256sum -c` of that file passes in the archive, and also in the run directory
`/fred/oz022/tkimpson/infercell-kim178-2026-11-27` it was copied from.

To regenerate the numbers from the archive:

    julia --project dev/scripts/corea_replicates/corea_replicates_reads.jl <archive>     # §1–§3
    julia --project dev/scripts/corea_replicates/corea_replicates_diagnose.jl <archive>  # §4

The per-coordinate z for all 21 coordinates of each replicate is in the archive's
`seed_<seed>/fit.md`, and KIM-92's is in `kim92/kim92.md` (prereg §6). The times below are
Slurm's (`sacct`), and the per-stage seconds are printed by the reads script.

| job | what | time |
|---|---|---|
| 21715763–65 | first smoke, at 83a7e40: 21715764 failed on a Symbol name, fixed in f66e0ec; 21715765 never ran | — |
| 21719135–37 | smoke at f66e0ec, seed 17_200_199 at tiny sizes (plumbing only) | 33 min |
| 21727400 | prereg §5: fitted MAP and Laplace on KIM-92's dataset | 2.06 h |
| 21770002_[1-10] | the ten replicates, task i is seed 17_200_001 + i | 7.8–12.2 h each, 96.1 CPU-h |

## 1. KIM-92's dataset first (prereg §5, job 21727400)

- **The MAP check holds.** The fitted log density is −123852.9728 at the MAP against
  −123856.6815 at θ_ref, a gain of 3.71. −H is positive definite, with eigenvalues 128.4 to 5.857e4.
- **Laplace matches NUTS here.** Laplace SD / NUTS SD (1,600 pooled draws) runs 0.960 to 1.025
  across the 21 coordinates at the fitted MAP, and 0.960 to 1.083 at θ_ref.
- **The mean z agrees.** The promoters' mean z is −0.797 by MAP and Laplace, against −0.808
  by NUTS, with 14 of 17 promoters below the truth. At θ_ref, on θ_ref's own Laplace SD, it is −0.855.

So on KIM-92's dataset, MAP and Laplace stand in for NUTS. This is the only dataset with
NUTS draws. Several replicates are less well identified: the smallest eigenvalue of −H is
10.2 on 003, 20.7 on 007, 26.2 on 010 and 34.2 on 005, against 128.4 here. Laplace is
not checked there.

## 2. The ten replicates (array 21770002)

All ten fits succeed under prereg §3, so no spares were used.

"Mean z at θ_ref" divides θ_ref's error by the fitted MAP's Laplace SD. Prereg §6 asked for
θ_ref's own SD, but no Hessian is computed at θ_ref on a replicate. z and z at θ_ref share
the truth and the denominator, so their agreement is partly built in.

| seed | promoters' mean z | below the truth | correlation-aware SE | mean z / SE | mean z at θ_ref |
|---|---|---|---|---|---|
| 17200002 | −2.156 | 16 of 17 | 0.246 | −8.8 | −2.913 |
| 17200003 | +0.028 | 8 of 17 | 0.256 | +0.1 | +0.031 |
| 17200004 | −0.321 | 13 of 17 | 0.268 | −1.2 | +0.009 |
| 17200005 | −2.721 | 17 of 17 | 0.274 | −9.9 | −3.021 |
| 17200006 | +0.434 | 7 of 17 | 0.223 | +1.9 | +1.274 |
| 17200007 | +0.064 | 7 of 17 | 0.255 | +0.3 | +0.069 |
| 17200008 | −0.325 | 12 of 17 | 0.249 | −1.3 | −0.258 |
| 17200009 | −0.361 | 11 of 17 | 0.256 | −1.4 | −0.369 |
| 17200010 | −1.547 | 17 of 17 | 0.263 | −5.9 | −1.564 |
| 17200011 | −2.945 | 17 of 17 | 0.271 | −10.9 | −3.840 |

## 3. Gate 2 (prereg §4): passes as frozen

- (a) The mean z is below zero in **7 of 10**. The gate fails at 9 or more, so (a) does not fire.
- (b) The pooled mean is **−0.985**, with SD 1.245 and SE 0.394, giving **t = −2.50** on 9 df.
  The gate fails at |t| > 3, so (b) does not fire.
- **Verdict as pre-registered: the transcript factor passes gate 2.**

**The pass does not mean the factor is calibrated.** Under calibration, each replicate's
mean z would scatter about zero with roughly its own correlation-aware SE, about 0.25. Instead:

- The ten means have an SD of 1.245, which is 4.9 times the mean SE.
- Four replicates (002, 005, 010, 011) are individually 5.9 to 10.9 SE below zero, and all
  four are negative. Three of them have all 17 promoters below the truth.
- Σ(mean z / SE)² = 337.3 on 10 df.

These SEs come from the Laplace approximation on each replicate. If Laplace understates the
SD on a replicate, its mean z / SE is overstated too. It would have to be several-fold too
narrow on 005 and 011 to account for them. On KIM-92's dataset it is within 4% of NUTS.

Test (b) takes its SE from the spread between replicates. A bias that varies with the truth
widens that spread and lets the test pass. So (b) asks whether the average bias is zero, not
whether the factor is calibrated. The gate's definition did not anticipate this. It is
recorded in the Decision log entry of 2026-11-28 as a post-hoc observation, not a change to
the gate.

The fitted-likelihood MAP is not the cause. The offset is already present at θ_ref, the MAP
with δ = S = 0, and it is larger there on the worst replicates (−3.84 against −2.95 on 011).

## 4. What the offset truths share (post hoc, after gate 2 was read; n = 10)

These come from `corea_replicates_diagnose.jl` on the archive. "Truth X" is the truth's
coordinate in prior SDs from the prior median. r is the Pearson correlation, across the ten,
of a quantity with the promoters' mean z.

At n = 10 each r is wide, and the two worst replicates (005 and 011) carry much of it. The
script prints a Fisher 95% interval, the rank correlation and r without those two.

- **The offset is mostly a shift shared by the promoters.** On the four bad replicates it is
  −0.046 to −0.086 in log units, so the promoters are 4.5 to 8.2% low. Their Laplace SD is
  about 0.03, so a small absolute bias is many SDs. Within a replicate, the promoters' z
  have an SD of 0.78 to 1.84, against a mean of −1.5 to −2.9.
- **It is associated with the size of δ.** δ is the reference offset of the stochastic
  ensemble from the mean field at θ_ref. For max |δ|, r = −0.84 (95% interval −0.96 to
  −0.46; rank r −0.68). It stays at −0.78 without 005 and 011. The two worst replicates have
  max |δ| of 0.90 and 1.04, against 0.09 to 0.28 for the other eight. For tr S, r = −0.71,
  but this falls to −0.04 without those two, so it rests on them.
- **It is associated with slow mRNA degradation in the truth.** For truth krnadeg, r = +0.78
  (interval +0.29 to +0.94; rank r +0.90; +0.54 without 005 and 011). The worst two have
  truths of −1.35 (005) and −1.02 (011) prior SDs. 010 sits at −0.09 with a −1.5 offset,
  so krnadeg does not account for every bad replicate.
- **A bulk-side coordinate is off as well.** z for kcat FBA has r = −0.72 (interval −0.93
  to −0.16). It is +4.45 on 002 and +2.85 on 011. z for σ_b has r = +0.53, with an interval
  that crosses zero. So the miscalibration is not confined to the promoters.
- **Weak or absent associations.**
  - The promoters' own truth: r = −0.61, with an interval that crosses zero, and all ten
    means lie within ±0.41 prior SD of the median.
  - The MAP's gain in log density: r = +0.05.
  - σ_b's truth: r = −0.22.
- **The promoters–krnadeg trade-off is modest.** The posterior correlation of the promoters'
  mean with log_krnadeg is +0.12 to +0.45, and it is smallest on the two worst (005 and 011).
  The krnadeg error is at most 0.014 in log units.

| seed | mean z | truth krnadeg | max abs δ | tr S | promoter error (log) | z FBA |
|---|---|---|---|---|---|---|
| 17200002 | −2.156 | −0.59 | 0.276 | 0.415 | −0.0529 | +4.45 |
| 17200003 | +0.028 | +0.50 | 0.123 | 0.355 | +0.0040 | −1.74 |
| 17200004 | −0.321 | +1.94 | 0.122 | 0.518 | −0.0115 | +0.87 |
| 17200005 | −2.721 | −1.35 | 0.898 | 0.644 | −0.0734 | +1.33 |
| 17200006 | +0.434 | +0.57 | 0.115 | 0.354 | +0.0191 | −0.30 |
| 17200007 | +0.064 | +0.44 | 0.118 | 0.312 | +0.0029 | −0.10 |
| 17200008 | −0.325 | −0.04 | 0.104 | 0.272 | −0.0078 | −0.14 |
| 17200009 | −0.361 | −0.06 | 0.0898 | 0.269 | −0.0129 | −2.40 |
| 17200010 | −1.547 | −0.09 | 0.132 | 0.259 | −0.0461 | −0.14 |
| 17200011 | −2.945 | −1.02 | 1.04 | 0.785 | −0.0856 | +2.85 |

**A hypothesis, untested.**
- **Where δ and S act.** In the mixed route, δ and S enter only the bulk factor
  (`bulk_loglik`, `src/mixed_likelihood.jl`). The transcript factor builds its rates from
  route (ii)'s mean-field states and never reads them.
- **Where the offset lives.** The offset is already present at θ_ref, where δ = S = 0, and
  KIM-116 placed it in the transcript factor. So δ and S do not cause it.
- **What a large δ likely signals.** It marks a truth where the mean field sits far from the
  stochastic mean, as it does with slow mRNA turnover. The transcript factor's mean-field
  rates are uncorrected there, and the error would surface as a shared shift in the promoter
  strengths.
- **How to test it.** Recomputing δ and S and refitting would change only the bulk factor,
  so it could come back null. The direct test is a transcript-only fit (KIM-116's cut,
  stage 1) on 005 and 011, set beside the transcript rates from the reference ensemble.
  It is not run here.

## 5. For KIM-187

Gate 2 passes as frozen. On these ten truths, though, the transcript factor's bias depends on
the truth. It is near zero on about half, and 6 to 11 posterior SEs on four, as Laplace
measures them. The posterior's stated uncertainty does not cover it. The bias is associated
(n = 10) with a large gap between the stochastic ensemble and the mean field. The full cell
inherits the factor, so the route decision should take the bias as real and truth-dependent,
not as one draw's noise.
