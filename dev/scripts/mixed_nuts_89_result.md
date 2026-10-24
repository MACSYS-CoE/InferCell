# Issue #89 B5: the mixed likelihood and NUTS on m0_v6

**Built:** `src/mixed_likelihood.jl`, the approximate posterior over
(log θ, log σ_b) of §3. It has three factors:
- B1's exact transcript factor, at the transcription rates that route (ii)
  rebuilds;
- a Gaussian bulk factor N(log D_met; m(θ) + δ, σ_b² I + S);
- the priors.

δ and S come from an 800-cell reference ensemble at θ_ref, the bounded MAP, by the
delta method (§12 2026-10-23, "B5 as planned"). The density is −∞ outside the
prior ± 6 SD box and when route (ii) needs more than 10⁴ steps (§12 2026-10-23,
"B5's support box and step cap").

**Suite:** 28791 of 28791 at `469ca26` (job 19084578). It includes B5's
verification clause: the clock identity, the bulk reduction, Haldane on duals, the
delta-method S, and AV3 on the whole density to 1e-6 relative. The clock identity
holds the last column to 1e-12 and the earlier ones to 1e-6 against separate solves
(§12 2026-10-24, "B5 as run").

**The run:** every stage ran. The four chains' draws have no divergent or
non-finite steps, and R̂ ≤ 1.012. It cost **41.5 CPU-h**, inside KC's bound of
50.75 (a tenth of F10's 507.5, job 18751635) on this dataset at this run length.
The margin is 18% and depends on tree depth; B7's refits read KC in general. B5 has no pass line on the posterior. AV6 (B7) judges
it. The comparison with F10 below is information only.

## Runs

- **Commit:** every stage ran at `469ca26`, the PR branch's tip when it was
  launched.
- **Archive:** `/fred/oz022/tkimpson/infercell-89-b5-archive-2026-10-23/`, with
  `SHA256SUMS`.
- **Rerun:** `NADAPT=300 NDRAW=300 dev/scripts/mixed_nuts_89.sh <outdir>` reruns
  every stage.
- **Header:** the merge's `summary.md` is headed "B5 smoke", the script's label.
  `SMOKE = false` in that file confirms this is the full run. Its closing line,
  that "a smoke's chains are short, so this is a lower bound", is also the
  script's smoke wording and does not apply here.
- **Reproducibility:** the chains are seeded (`Xoshiro(1000 + c)`), but a rerun
  on other nodes will agree statistically, not to the digit.

| What | Job | Elapsed |
|---|---|---|
| Precompile | 19091859 | 0:02:22 |
| `theta_ref` (bounded MAP at δ = S = 0) | 19091860 | 2:02:13 |
| Reference ensemble, 16 chunks of 50 cells | 19091861–19091876 | 0:04:19–0:06:36 each |
| NUTS chains 1–4, 300 adaptation + 300 draws each | 19091877–19091880 | 9:24:33, 9:13:34, 10:27:54, 8:53:39 |
| Merge (`summary.md`) | 19091881 | 0:01:31 |

## The data and the reference θ

- **Data:** m0_v6, V6's replicate 16085: 50 cells, a 600 s horizon, 10 saves,
  2 genes and a 17-pool panel.
- **θ_ref:** the bounded MAP converged at log density −591.04 in 6793 s. The
  prior-median fallback was not needed.

| coordinate | θ_ref (log) | truth (log) |
|---|---|---|
| log_S_JCVISYN3A_0607 | 1.6772 | 1.4892 |
| log_S_JCVISYN3A_0779 | 0.7434 | 1.0121 |
| log_krnadeg | 2.3568 | 2.2433 |
| log_kcatF_R_ENO | 4.2467 | 4.1818 |
| log_kcatF_R_FBA | 4.3172 | 4.2066 |
| log_sigma_b | −1.9053 | −2.6220 |

**θ_ref's σ_b (0.149) is about twice the truth's (0.073).** That is consistent
with the MAP being taken at δ = 0 and S = 0: σ_b² rises from 0.0053 to 0.0221,
about what S's mean diagonal (≈ 0.010) and δ² (≈ 0.003) would add. This is an
explanation, not a test. The chains bring it back to −2.621.

## δ and S

δ is small on M0 at θ_ref:
- every pool's mean over the saves lies in [−0.084, +0.114] nats;
- the largest single entry is M_gtp_c's +0.235;
- the full per-pool table is in `summary.md`.

S, with C = 50:
- trace 1.718;
- the five largest eigenvalues are 0.411, 0.215, 0.165, 0.154 and 0.134;
- the smallest is 7.4e-11.

So S is full rank in exact arithmetic, though numerically nearly singular in some
directions. The σ_b² I term
(0.022 at θ_ref) keeps the covariance well conditioned.

## The chains

| coordinate | mean | sd | ESS (bulk) | R̂ |
|---|---|---|---|---|
| log_S_JCVISYN3A_0607 | 1.5161 | 0.1057 | 1291 | 1.000 |
| log_S_JCVISYN3A_0779 | 0.6215 | 0.1835 | 1642 | 0.999 |
| log_krnadeg | 2.2636 | 0.0595 | 1319 | 1.000 |
| log_kcatF_R_ENO | 4.2585 | 0.0430 | 1812 | 1.012 |
| log_kcatF_R_FBA | 4.3977 | 0.1877 | 1905 | 0.998 |
| log_sigma_b | −2.6212 | 0.0603 | 1706 | 1.000 |

- **Divergences:** 0 divergent or non-finite steps in the 1200 kept draws, so no
  −∞ (box, step cap, support or overdraw) was hit there. The 1200 adaptation
  iterations are discarded by Turing and not counted.
- **Starts:** all four chains start at θ_ref, so R̂ here checks agreement from a
  common start, not from overdispersed ones. B7 should jitter the starts.
- **Trees:** depth median 3 and max 4, with a median of 7 leapfrog steps per draw.
- **ESS:** the bulk ESS is above the 1200 draws in all six coordinates. NUTS draws
  can be anticorrelated, which would do this; the autocorrelation was not
  measured, and the lowest values (1291, 1319) are within the estimator's noise of
  1200.

## Cost (KC)

- **Per gradient:** 6.13 s and 5.22 GiB allocated (median of 3, at the prior
  median). One density costs 2.52 s.
- **By stage:** the chains took 37.9 CPU-h for 7900 leapfrog steps in their draws.
  The MAP took 1.9 and the ensemble 1.0.
- **Total by `sacct`:** 41.5 CPU-h. That includes precompile, the merge and each
  job's start-up, which the script's own count of 40.74 leaves out.
- **Against the projection:** the projection was about 59 CPU-h, from smoke r3's
  15 leapfrog steps per iteration. At 6.13 s per gradient the chains' 37.9 CPU-h
  is about 22,200 gradients. The draws' 7900 leave about 14,300 for the 1200
  adaptation iterations, about 12 each. So the draws' shallower trees (7 steps)
  explain most of the saving; adaptation was below the projection's 15 but its steps are
  not recorded directly.

The gradient is still the cost. A login-node profile (not a Slurm job) put the dual
solve of route (ii) at about 5.1 s of a 5.6 s gradient. Smoke r3 (3 cells, 120 s)
cost nearly the same per gradient, 5.35 s against 6.13 s. The handoff recommends
profiling it first in B7.

## Beside F10's exact posterior (information only)

F10's draws are its reference chains 1–4, sweeps 101–1,650
(`dev/scripts/f10_16a10_draws.tsv`). Both are on the log scale.

| coordinate | mixed: mean ± sd | F10: mean ± sd | difference / F10 sd |
|---|---|---|---|
| log_S_JCVISYN3A_0607 | 1.5161 ± 0.1057 | 1.5454 ± 0.1159 | −0.25 |
| log_S_JCVISYN3A_0779 | 0.6215 ± 0.1835 | 0.6999 ± 0.1995 | −0.39 |
| log_krnadeg | 2.2636 ± 0.0595 | 2.2246 ± 0.0763 | +0.51 |
| log_kcatF_R_ENO | 4.2585 ± 0.0430 | 4.2423 ± 0.0419 | +0.39 |
| log_kcatF_R_FBA | 4.3977 ± 0.1877 | 4.4077 ± 0.1844 | −0.05 |
| log_sigma_b | −2.6212 ± 0.0603 | −2.6250 ± 0.0601 | +0.06 |

Every mean lies within 0.51 of F10's posterior SD of F10's own mean. The SDs agree
to within 22%: krnadeg's is the narrowest, at 0.78 of F10's. This is one dataset
and has no pass line. AV6 in B7 is the test.
