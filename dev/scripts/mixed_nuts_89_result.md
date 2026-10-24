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
delta-method S, and AV3 on the whole density to 1e-6 relative.

**The run:** every stage ran. The four chains have no divergent or non-finite steps
and R̂ ≤ 1.012. It cost **41.5 CPU-h**, inside KC's bound of 50.75 (a tenth of
F10's 507.5, job 18751635). B5 has no pass line on the posterior. AV6 (B7) judges
it. The comparison with F10 below is information only.

## Runs

- **Commit:** every stage ran at `469ca26`, the PR branch's tip when it was
  launched.
- **Archive:** `/fred/oz022/tkimpson/infercell-89-b5-archive-2026-10-23/`, with
  `SHA256SUMS`.
- **Rerun:** `NADAPT=300 NDRAW=300 dev/scripts/mixed_nuts_89.sh <outdir>` reruns
  every stage.
- **Header:** the merge's `summary.md` is headed "B5 smoke", the script's label.
  `SMOKE = false` in that file confirms this is the full run.

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

**θ_ref's σ_b (0.149) is about twice the truth's (0.073).** That is expected. The
MAP is taken at δ = 0 and S = 0, so σ_b absorbs the Jensen gap and the cell-to-cell
spread that δ and S later carry. The chains bring it back to −2.621.

## δ and S

δ is small on M0 at θ_ref:
- every pool's mean over the saves lies in [−0.084, +0.114] nats;
- the largest single entry is M_gtp_c's +0.235;
- the full per-pool table is in `summary.md`.

S, with C = 50:
- trace 1.718;
- the five largest eigenvalues are 0.411, 0.215, 0.165, 0.154 and 0.134;
- the smallest is 7.4e-11.

So S is full rank, though nearly singular in some directions. The σ_b² I term
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

- **Divergences:** 0 divergent or non-finite steps in 1200 draws, so the box and
  the step cap never bound.
- **Trees:** depth median 3 and max 4, with a median of 7 leapfrog steps per draw.
- **ESS:** the bulk ESS is above the 1200 draws in five coordinates. NUTS draws can
  be anticorrelated, and these are.

## Cost (KC)

- **Per gradient:** 6.13 s and 5.22 GiB allocated (median of 3, at the prior
  median). One density costs 2.52 s.
- **By stage:** the chains took 37.9 CPU-h for 7900 leapfrog steps in their draws.
  The MAP took 1.9 and the ensemble 1.0.
- **Total by `sacct`:** 41.5 CPU-h. That includes precompile, the merge and each
  job's start-up, which the script's own count of 40.74 leaves out.
- **Against the projection:** I projected about 59 CPU-h from smoke r3's 15
  leapfrog steps per draw. These chains took 7, so they cost less.

The gradient is still the cost, and the dual solve of route (ii) is almost all of
it. Smoke r3 (3 cells, 120 s) cost the same per gradient. Profiling it is B7's
first job.

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
