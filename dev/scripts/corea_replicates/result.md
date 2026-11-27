# KIM-178 result: ten replicate truths on the reduced model, by MAP and Laplace

The design was frozen at tag `prereg/kim-178` (2050f5c) before any replicate existed. Every
job ran at commit `f66e0ec`. The archive is `/fred/oz022/tkimpson/infercell-kim178-archive-2026-11-28`
(80 files). `sha256sum -c SHA256SUMS` passes on both the archive and the run directory
`/fred/oz022/tkimpson/infercell-kim178-2026-11-27`.

To regenerate the numbers from the archive:

    julia --project dev/scripts/corea_replicates/corea_replicates_reads.jl <archive>     # §1–§3
    julia --project dev/scripts/corea_replicates/corea_replicates_diagnose.jl <archive>  # §4

| job | what | time |
|---|---|---|
| 21719135–37 | smoke, seed 17_200_199 at tiny sizes (plumbing only) | 29 min |
| 21727400 | prereg §5: fitted MAP and Laplace on KIM-92's dataset | 2.06 h |
| 21770002_[1-10] | the ten replicates, task i is seed 17_200_001 + i | 7.8–12.2 h each, 96.1 CPU-h |

## 1. KIM-92's dataset first (prereg §5, job 21727400)

- **The MAP check holds.** The fitted log density is −123852.9728 at the MAP against
  −123856.6815 at θ_ref, a gain of 3.71. −H is positive definite, with eigenvalues 128.4 to 5.857e4.
- **Laplace matches NUTS.** Laplace SD / NUTS SD (1,600 pooled draws) runs 0.960 to 1.025
  across the 21 coordinates at the fitted MAP, and 0.960 to 1.083 at θ_ref.
- **The mean z agrees.** The promoters' mean z is −0.797 by MAP and Laplace, against −0.808
  by NUTS, with 14 of 17 promoters below the truth.

So MAP and Laplace stand in for NUTS on this model. Whatever the ten show belongs to the
posterior, not to the approximation.

## 2. The ten replicates (array 21770002)

All ten fits succeed under prereg §3, so no spares were used.

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

- The ten means have an SD of 1.245, about five times that.
- Four replicates (002, 005, 010, 011) are individually 5.9 to 10.9 SE below zero, and all
  four are negative. Three of them have all 17 promoters below the truth.
- Σ(mean z / SE)² = 337 on 10 df. Calibration cannot produce that.

Test (b) takes its SE from the spread between replicates. A bias that varies with the truth
widens that spread and lets the test pass. So (b) asks whether the average bias is zero, not
whether the factor is calibrated. The gate's definition did not anticipate this. It is
recorded in the Decision log entry of 2026-11-28 as a post-hoc observation, not a change to
the gate.

The mean z at θ_ref (r = +0.99 with the MAP's) shows that the fitted-likelihood MAP is not
the cause. The bias is already present at the δ = S = 0 MAP.

## 4. What the offset truths share (exploratory, after gate 2 was read; n = 10)

These come from `corea_replicates_diagnose.jl` on the archive. "Truth X" is the truth's
coordinate in prior SDs from the prior median. r is the correlation, across the ten, of a
quantity with the promoters' mean z.

- **The offset is common to all the promoters.** In log units it is −0.05 to −0.09 on the
  four bad replicates, which is 5 to 9% low. The promoters' Laplace SD is about 0.03, so a
  small absolute bias is many SDs.
- **It follows the size of δ.** The reference offset δ (stochastic ensemble against mean
  field at θ_ref) gives r = −0.84 for max |δ| and −0.71 for tr S. The two worst replicates
  have max |δ| of 0.90 and 1.04, against 0.09 to 0.28 for the other eight.
- **It follows slow mRNA degradation in the truth.** Truth krnadeg gives r = +0.78. The
  worst replicates have truths of −1.35 (005) and −1.02 (011) prior SDs.
- **It follows the promoters' own truth only weakly.** Their mean in prior SDs gives
  r = −0.61, but all ten lie within ±0.41 of the prior median. It does not follow the
  MAP's gain in log density (r = +0.05) or σ_b's truth (r = −0.22).
- **The promoters–krnadeg trade-off is modest.** The posterior correlation of the promoters'
  mean with log_krnadeg is +0.12 to +0.45, and it is smallest on the two worst (005 and 011).
  The krnadeg error is at most 0.014 in log units.

| seed | mean z | truth krnadeg | max abs δ | tr S | promoter error (log) |
|---|---|---|---|---|---|
| 17200002 | −2.156 | −0.59 | 0.276 | 0.415 | −0.0529 |
| 17200003 | +0.028 | +0.50 | 0.123 | 0.355 | +0.0040 |
| 17200004 | −0.321 | +1.94 | 0.122 | 0.518 | −0.0115 |
| 17200005 | −2.721 | −1.35 | 0.898 | 0.644 | −0.0734 |
| 17200006 | +0.434 | +0.57 | 0.115 | 0.354 | +0.0191 |
| 17200007 | +0.064 | +0.44 | 0.118 | 0.312 | +0.0029 |
| 17200008 | −0.325 | −0.04 | 0.104 | 0.272 | −0.0078 |
| 17200009 | −0.361 | −0.06 | 0.0898 | 0.269 | −0.0129 |
| 17200010 | −1.547 | −0.09 | 0.132 | 0.259 | −0.0461 |
| 17200011 | −2.945 | −1.02 | 1.04 | 0.785 | −0.0856 |

**The hypothesis, untested.** The mixed route fixes δ and S at θ_ref. Where the truth puts the
cell in a strongly stochastic regime (slow mRNA turnover, so a large gap between the ensemble
and the mean field), δ is large. Any error in δ then passes straight into a common shift of
the promoter strengths, which is the transcript factor. Replicate 010 has an ordinary δ
(0.132) and a −1.5 offset, so δ's size is not the whole story.

The direct test is to recompute δ and S at the fitted MAP and refit, on 005 and 011 (about
2 × 8 h). It is not run here.

## 5. For KIM-187

Gate 2 passes as frozen. On these ten truths, though, the transcript factor's bias depends on
the truth: zero on about half, and 6 to 11 posterior SEs on four. The posterior's stated
uncertainty does not cover it. The full cell inherits the factor, so the route decision
should take the bias as real and truth-dependent, not as one draw's noise.
