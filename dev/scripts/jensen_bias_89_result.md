# Issue #89, B3 close-out: the parameter bias the Jensen gap implies

**Rule:** sub-spec §12 2026-10-23 ("the gap is a declared reduction"). The bias is
reported with no pass line. AV6 (B7) is the test of whether the gap moves the
posterior; this is a linearised estimate.
**Script:** `dev/scripts/jensen_bias_89.jl`, with `.slurm`. It reads the rethink's
archive (`/fred/oz022/tkimpson/infercell-89-rethink-archive-2026-10-22/`, checked
unchanged against its `SHA256SUMS` after the runs) and simulates nothing.

| Job | Commit | What | CPU-h |
|---|---|---|---|
| 19038595 | `774ba2d` | First run. The m0_v6 section lacked SEs and noise; superseded, kept in `superseded_19038595/`. The rest is identical to the run of record | 0.03 |
| 19038794 | `c9184c2` | **Of record** | 0.02 |

**Outputs:** `/fred/oz022/tkimpson/infercell-89-bias-2026-10-23/` (`bias.md`,
`bias.jls`, the logs, `SHA256SUMS`).
**Regenerate:** `sbatch --output=<outdir>/%x_%j.log dev/scripts/jensen_bias_89.slurm
<outdir>` from a worktree at `c9184c2`. It takes about 1 minute.

## What was computed

bias ≈ (JᵀWJ + P)⁻¹ JᵀW · gap, in log θ, divided by F10's posterior SD per
parameter.
- **gap:** at each of 2A's 28 targets, log x̄ of its 200 cells on independent seeds,
  minus (m + δ_ref) or minus m (uncorrected). The KG totals this reproduces match
  `transfer.md` (job 19018422) to the printed figure.
- **J:** from the archived mean fields. The prototype rounds counts inside the
  driver, so ForwardDiff (AV3) does not reach it.
  - At m0_low and m0_high: central differences over the ±1 posterior SD targets.
    The curvature ratio ‖m₊ + m₋ − 2m₀‖ / ‖m₊ − m₋‖ is at most 0.12.
  - At f10_ref: least squares over the 28 F10 draws' mean fields, R² = 0.989.
  - So m(θ) is close to linear over the posterior's width.
- **W = (σ_b² I + S)⁻¹:** S = Cov(log x̄₅₀) at the source, by the delta method from
  its 800 cells. σ_b is KG's 0.10 and, as a check, F10's posterior median, 0.072.
  The two give the same picture; the numbers below are at 0.10.
- **P:** the prior precisions in log θ (all LogNormal).
- **SE:** the gap's sampling error, from the target's 200 cells and δ's 800,
  propagated linearly.

**The formula leaves out the transcripts.** (JᵀWJ + P)⁻¹ is the bulk-plus-prior
posterior. For S_0607, S_0779 and krnadeg its SD is 0.8–5.5× F10's, so it
overstates the shift in those parameters. Beside it, at the 8 F10 targets, is the
linear response with F10's posterior covariance, Σ_F10 JᵀW · gap, which carries the
transcripts' information (information only).

## Results (σ_b = 0.10)

**Uncorrected, the mean field biases krnadeg and S_0607 by several posterior SDs.**

| Source | Read | S_0607 | S_0779 | krnadeg | kcatF_ENO | kcatF_FBA |
|---|---|---|---|---|---|---|
| f10_ref (8 targets) | spec formula, RMS | 7.37 | 0.48 | 5.98 | 0.17 | 0.31 |
| f10_ref (8 targets) | Σ_F10, range | +0.87 to +1.63 | −0.41 to +0.48 | **+1.91 to +2.40** | −0.43 to −0.03 | −0.08 to +0.24 |
| m0_low (10) | spec formula, RMS | 2.95 | 0.59 | 6.91 | 0.32 | 0.23 |
| m0_high (10) | spec formula, RMS | 0.86 | 1.62 | 6.23 | 0.39 | 0.72 |

- |bias| > 2 SE at 99 of 140 target × parameter cells.
- With the transcripts' information (Σ_F10), krnadeg still moves by about 2
  posterior SDs at every F10 target, and S_0607 by about 1.

**With δ_ref, the residual bias is within about 1 SD everywhere, and mostly within
its noise.**

| Source | Read | S_0607 | S_0779 | krnadeg | kcatF_ENO | kcatF_FBA |
|---|---|---|---|---|---|---|
| f10_ref (8) | spec formula, RMS | 1.02 | 0.53 | 0.49 | 0.17 | 0.09 |
| f10_ref (8) | Σ_F10, max \|bias\| | 0.40 | 0.45 | 0.33 | 0.25 | 0.22 |
| m0_low (10) | spec formula, RMS | 0.45 | 0.11 | 1.13 | 0.31 | 0.03 |
| m0_high (10) | spec formula, RMS | 0.79 | 0.17 | 0.74 | 0.21 | 0.03 |

- |bias| > 2 SE at 10 of 140 cells. About 7 would be expected by chance.
- **Resolution:** the SEs are 0.4–0.9 SD on S_0607 and 0.5–1.2 SD on krnadeg, from
  the 200-cell targets. So δ_ref's residual is not resolved below about 1 SD there.
  kcatF_FBA is resolved (SE ≤ 0.06) and shows at most 0.17 SD.
- **One resolved residual on ENO:** at low_kcatF_R_ENO_−, kcatF_ENO's bias is
  −0.80 ± 0.24 under δ_ref, the same as uncorrected (−0.82 ± 0.21). δ does not
  remove that component there. It is the only target where δ leaves a bias above
  3 SE.
- **At σ_b = 0.072** (`bias.md`):
  - δ_ref's RMS values move by at most 0.12 SD, and its count above 2 SE is 9 of
    140.
  - The uncorrected RMS values grow by up to 0.84 SD (krnadeg at m0_low, 7.75),
    and the count above 2 SE stays at 99.

## δ at m0_v6, re-read on independent seeds (information only)

This corrects what #97 added to the rethink. That read used the shared-seed
targets, and on independent seeds it does not hold up.
- **m0_v6 is a typical point of F10's posterior:** Mahalanobis 2.48, against the
  draws' median of 2.06 and 95% point of 3.29. It is not "1–3 posterior SDs away".
- **δ_v6 is not distinguishable from δ_ref:** ‖δ_v6 − δ_ref‖ = 0.642 over all pools
  and saves, against 0.744 expected from the two δs' sampling noise alone.
- **KG totals on independent seeds:** δ_v6 does better than uncorrected at 5 of 8
  draws and worse than δ_ref at 7 of 8. That is what its 200 cells predict: a
  right δ from 200 cells should score about the noise, and one from 800 about 0.6
  of it.
  - δ_v6's totals ÷ noise are 0.5–1.1 at six draws.
  - The exceptions are f10_d4 (3.7) and f10_d6 (1.9).
- **Its implied biases** are mostly within 2 SE. S_0607 at f10_d7 and f10_d8 is
  −2.1 ± 0.9 and −2.9 ± 0.9, and kcatF_FBA at f10_d2 and f10_d5 is about −0.2 ±
  0.07.
- **So:** whether δ depends on θ across the posterior is neither established nor
  ruled out at these ensemble sizes. The shared-seed table in `scores.md` stays as
  history.

## What this says for B4–B7

- **δ is needed.** Without it the bulk factor would move krnadeg by about 2
  posterior SDs, even with the transcripts' information.
- **With δ held at one θ, the linearised residual is within the noise of this
  test.** That noise is about 1 SD on S_0607 and krnadeg. AV6 measures the real
  posterior shift against the exact sampler, at V7's 50 datasets.
- **What AV6 should look at:** krnadeg and S_0607, where this test is weakest, and
  kcatF_ENO at low-ENO truths.
- **The linearisation:** it is good over the posterior's width (R² 0.989,
  curvature ≤ 0.12), so the gap is unlikely to act nonlinearly near the posterior.
