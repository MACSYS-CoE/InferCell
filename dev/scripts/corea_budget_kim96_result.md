# KIM-96: the reduced model's uncertainty budget once every parameter is inferred

**What this is.** KIM-42 split Q's parameter variance on the reduced model between θ_I, from
M0's posterior, and θ_H, the 15 promoter strengths held at their priors. It found 91% in θ_H.
KIM-96 was to repeat that split with every parameter inferred. **There is nothing left to
split:**
- KIM-92's joint posterior reaches all 21 coordinates.
- So the planned split by contraction (informed against prior-dominated) has an empty second
  group, and the held share is 0 by construction.

Instead, this note sets Q's parameter variance side by side in three states (Decision log
2026-11-19, Tom's call: no new compute).

**No new simulation.** The numbers come from three archives. One short Slurm job reads them,
and each input is checked against its archive's `SHA256SUMS` first.

**Run:**
- `dev/scripts/corea_budget_kim96.sh` at c507315 ran job **21188187**, which ended COMPLETED 0:0 after 14 s on dave24 (`sacct.tsv`).
- It wrote to `/fred/oz022/tkimpson/infercell-kim96-archive-2026-11-19/`:
  - `budget.md`, `budget_rows.tsv` and `inputs.tsv`;
  - the script that ran, and `COMMIT`;
  - `jobs.tsv`, `sacct.tsv` and `logs/`.
- `SHA256SUMS` covers the outputs, the script and `COMMIT`, written by the job. `jobs.tsv` and `sacct.tsv` were appended after it finished. All 7 lines pass `sha256sum -c`.
- The tables below are `budget.md`'s.

**Inputs** (`inputs.tsv` lists their SHA-256):
- **KIM-92** `pilot_coordinates.tsv`, from `infercell-kim92-archive-2026-11-17`, merge job 21052425.
- **KIM-95** `predictive_{prior,posterior}_{mu,s2,n}.tsv`, from `infercell-kim95-archive-2026-11-19`.
  - The runs were at baf5331: posterior 21154498, prior 21154499, collect 21154500.
- **KIM-42** `budget_rows.tsv`, from `infercell-m5-kim42-2026-11-05`, at f23fd06.
  - The jobs were 20041038–40; width 1 is task 1, raw id 20041246.

**Q is the same quantity in all three:** `:lactate_out`, the population mean over cells at
6,300 s, on the hybrid. KIM-42 passes `quantity = :lactate_out` and the 6,300 s horizon
explicitly (`reduced_budget_kim42.jl`), and KIM-95 passes the horizon
(`corea_predict_kim95.jl`). Everything else is `predict`'s default in both: Q is
`:lactate_out`, `forward` is `:hybrid` and `estimand` is `:population_mean`.

## There is no prior-dominated set

The threshold is frozen in `spec/prereg/kim-92.md` §3, which says "the data reach" a
parameter at contraction ≥ 0.5.
- On KIM-92's 1,600 pooled draws, **21 of 21** coordinates reach it.
- The lowest contraction is 0.919, on kcatF_R_ENO, which has the tightest prior (SD 0.157).
- Every promoter's contraction is 0.993 or more; the lowest is S_0651 at 0.99382. `budget.md` has the full table, at three decimals.

No threshold between 0.5 and 0.919 changes this. A threshold above 0.919 would be chosen after
the contractions were read (2026-11-17), so none is used.

## Q's parameter variance in three states

| state | what is uncertain | mean Q ± SE | Var_θ ± SE | √Var_θ | √Var_θ / mean | MC share |
|---|---|---|---|---|---|---|
| prior | all 20 drawn parameters at their priors (KIM-95 prior arm, 100 draws × 200 cells) | 1.1388e+07 ± 3.4e+05 | 1.153e+13 ± 1.37e+12 | 3.39e+06 | 0.298 | 0.0008 |
| KIM-42 | M0's posterior on 5 (V7 replicate 1); 15 promoters at their priors (width 1) | 9.5765e+06 ± 3.3e+05 | 2.169e+12 ± 5.93e+11 | 1.47e+06 | 0.154 | — |
| posterior | KIM-92's joint posterior on all 21 (KIM-95 posterior arm, 100 draws × 973 cells) | 1.5013e+07 ± 1.3e+04 | 1.291e+10 ± 2.25e+09 | 1.14e+05 | 0.0076 | 0.2108 |

The variance ratios use the delta method, on runs that are independent. With KIM-42's
relative SE at 27%, a ratio's sampling distribution is skewed, so `budget.md`'s linear ± SE
should not be read as a symmetric interval. The last column is the log-scale interval,
r·exp(±2·√((se_a/a)² + (se_b/b)²)), derived from `budget_rows.tsv`:

| ratio of Var_θ | estimate ± SE (`budget.md`) | ±2 SE on the log scale (derived) |
|---|---|---|
| prior / posterior | 893 ± 189 | 585–1,363 |
| KIM-42 / posterior | 168 ± 54 | 88–321 |
| prior / KIM-42 | 5.3 ± 1.6 | 2.9–9.6 |

**These ratios hold at the default prior width only** (geometric SD 2, KIM-42's width 1, which
is also KIM-95's prior arm). KIM-42's own note says its budget "depends strongly on the prior
width". Its Var(μ) is 5.65e11, 2.17e12 and 1.47e13 at widths ½, 1 and 2 (`budget_rows.tsv`,
`total`), so against the same posterior:
- KIM-42 / posterior runs about 44, 168 and 1,140;
- prior / KIM-42 runs about 20, 5.3 and 0.78, if the prior state stays at width 1.

That spread, set by the choice of width, is more than an order of magnitude, far larger than
the Monte Carlo SEs above. The prior state's 1.15e13 also rests on the chosen geometric SD of
2. The posterior's 1.29e10 does not, because the data dominate it.

In KIM-42, 0.908 ± 0.041 of the 2.17e12 sat in the 15 held promoters (I first; 0.898 ± 0.042
with H first).

**The sentence F6 can carry** (Decision log 2026-11-19 has the template):

*Once every parameter is inferred, none stays prior-dominated (contraction ≥ 0.919 on all
21). At the default prior width (geometric SD 2), Q's parameter variance falls in three steps:*
- *1.15e13 with every parameter at its prior;*
- *2.17e12 when M0's posterior constrains five and the 15 promoters stay at their priors, with
  91% of it in those 15;*
- *1.29e10 under the joint posterior.*

*The joint posterior is about 170 times below KIM-42's state (88–321) and about 900 times below
the prior's (585–1,363). As a standard deviation, Q's spread falls from 30% of Q, through 15%,
to 0.8%. KIM-42's step depends on its prior width: at widths ½ and 2 the factor is about 44 and
1,140.*

## How to read it

- **The 91% does not "fall to X%".** Once all 15 promoters are inferred, the question of how
  much variance the held parameters carry has no subject. What the data remove is most of the
  variance itself. The prior → KIM-42 step takes out a factor of about 5, and the step from
  KIM-42 to the joint posterior takes out about 170 more.
- **That second step is not only the promoters.** KIM-42's five M0-inferred parameters alone
  carried 2.0e11 ± 5.6e10 of Var(μ) at width 1. That is `i_between`, the variance of the mean
  over θ_H; with H first, `h_within` is 2.2e11. This is about 15 times the joint posterior's
  whole Var_θ. Inferring the 15 promoters while keeping M0's θ_I would leave about 2e11, a
  factor of about 11 below 2.17e12. The remaining factor of about 15 comes from KIM-92's data
  also tightening those five: krnadeg, for example, has posterior SD 0.0067 on the log scale
  and contraction 1.000 (`pilot_coordinates.tsv`). The joint posterior does both at once, so
  the two are not separated here.
- **The three states have different θ and different means.**
  - KIM-42's θ_I comes from V7 replicate 1, an M0 fit carried to the reduced model, and its
    mean Q is 9.6e6.
  - The posterior's mean is 1.50e7, against a truth of 1.498e7 ± 1.8e4 (KIM-95).
  - So the ratios compare variance around different centres. They are not a decomposition of
    one distribution.
- **The posterior's Var_θ is resolved, but only just.**
  - The cell noise is 21% of var(μ̂) across the 100 draws (the MC share), so Var_θ is 79% of
    it.
  - The jackknife SE is 17% of the estimate.
  - KIM-95's sizing pass, at 100 cells per draw, could not resolve it (job 21126961). The
    full run's 973 cells can.
- **V.** KIM-42 conditions on its viability rule, and KIM-95 applies none. At width 1, V kept
  every pair in KIM-42 (acceptance 1.0000), so its `per_theta_I` and `joint` rows are equal,
  and this source of difference vanishes.
- **What the KIM-42 / posterior ratio's SE omits.** It is dominated by KIM-42's jackknife
  (27%) and treats the runs as independent, which they are: different θ draws and cell
  seeds. It does not cover KIM-42's dependence on replicate 1's θ_I, which its own note
  flags.
- **A well-specified, synthetic case.** KIM-92's data were simulated by this same model at
  θ*, so there is no model-form error. The 0.8% is precision under a correct model.
  - It is not a calibrated uncertainty on every output: KIM-95 found 10 of the 17 panel pools
    offset by 1–1.4 posterior SD on the forecast window, and those offsets are θ-driven.
  - Q itself is covered: the truth sits 0.26 posterior SD below the median (KIM-95).
- **Parameter variance only.** A new cell's spread around Q adds the cells' own variance.
  On the posterior arm that is mean s² = 3.35e12 per cell (derived from `predictive_posterior_s2.tsv`),
  so Var_θ/σ² ≈ 0.0038. For one cell, the inferred θ's share of the variance is about 0.4%.

## Reproducing

```
dev/scripts/corea_budget_kim96.sh <new outdir>   # Slurm, from a clean worktree at a pushed commit
julia --project dev/scripts/corea_budget_kim96.jl <outdir>   # or directly: about 4 s, reads only
```

The direct run takes `KIM92_DIR`, `KIM95_DIR` and `KIM42_DIR` to point it at copies of the
archives. The archive's `sacct.tsv` was exported after the job with
`sacct -j 21188187 -o JobID,JobName,State,ExitCode,Elapsed,NodeList,TotalCPU,MaxRSS,Submit,Start,End -P`.
The archive holds the inputs' hashes (`inputs.tsv`) but not copies of the inputs, so off the
cluster the result can be checked against those hashes but not regenerated.

**Check:** the script was rerun on the login node, on copies of the three archives' inputs
and their `SHA256SUMS` (KIM92_DIR, KIM95_DIR and KIM42_DIR point it at the copies). It exited
0 and wrote a `budget.md` byte-identical to the job's (`cmp`).
