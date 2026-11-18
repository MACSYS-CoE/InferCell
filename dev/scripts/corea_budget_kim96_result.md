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
6,300 s, on the hybrid. In KIM-42 it is passed explicitly (`reduced_budget_kim42.jl`). In
KIM-95 it is `predict`'s default (`corea_predict_kim95.jl`).

## There is no prior-dominated set

The threshold is frozen in `spec/prereg/kim-92.md` §3, which says "the data reach" a
parameter at contraction ≥ 0.5.
- On KIM-92's 1,600 pooled draws, **21 of 21** coordinates reach it.
- The lowest contraction is 0.919, on kcatF_R_ENO, which has the tightest prior (SD 0.157).
- Every promoter's contraction is 0.994 or more (`budget.md` has the full table).

No threshold between 0.5 and 0.919 changes this. A threshold above 0.919 would be chosen after
the contractions were read (2026-11-17), so none is used.

## Q's parameter variance in three states

| state | what is uncertain | mean Q ± SE | Var_θ ± SE | √Var_θ | √Var_θ / mean | MC share |
|---|---|---|---|---|---|---|
| prior | all 20 drawn parameters at their priors (KIM-95 prior arm, 100 draws × 200 cells) | 1.1388e+07 ± 3.4e+05 | 1.153e+13 ± 1.37e+12 | 3.39e+06 | 0.298 | 0.0008 |
| KIM-42 | M0's posterior on 5 (V7 replicate 1); 15 promoters at their priors (width 1) | 9.5765e+06 ± 3.3e+05 | 2.169e+12 ± 5.93e+11 | 1.47e+06 | 0.154 | — |
| posterior | KIM-92's joint posterior on all 21 (KIM-95 posterior arm, 100 draws × 973 cells) | 1.5013e+07 ± 1.3e+04 | 1.291e+10 ± 2.25e+09 | 1.14e+05 | 0.0076 | 0.2108 |

The variance ratios use the delta method, on runs that are independent:

| ratio of Var_θ | estimate ± SE |
|---|---|
| prior / posterior | 893 ± 189 |
| KIM-42 / posterior | 168 ± 54 |
| prior / KIM-42 | 5.3 ± 1.6 |

In KIM-42, 0.908 ± 0.041 of the 2.17e12 sat in the 15 held promoters (I first; 0.898 ± 0.042
with H first).

**The sentence F6 can carry** (Decision log 2026-11-19 has the template):

*Once every parameter is inferred, none stays prior-dominated (contraction ≥ 0.919 on all
21). Q's parameter variance falls from 1.15e13 with every parameter at its prior, through
2.17e12 when M0's posterior fixes five and the 15 promoters stay at their priors (91% of it in
those 15), to 1.29e10 under the joint posterior. That is a factor of about 170 ± 50 below
KIM-42's, and about 900 ± 190 below the prior's. As a standard deviation, it falls from 30% of
Q, through 15%, to 0.8%.*

## How to read it

- **The 91% does not "fall to X%".** Once all 15 promoters are inferred, the question of how
  much variance the held parameters carry has no subject. What the data remove is most of the
  variance itself: the prior → KIM-42 step takes out a factor of about 5, and the step that
  infers the promoters takes out about 170 more.
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
- **Parameter variance only.** A new cell's spread around Q adds the cells' own variance.
  On the posterior arm that is mean s² = 3.35e12 per cell (derived from `predictive_posterior_s2.tsv`),
  so Var_θ/σ² ≈ 0.0038. For one cell, the inferred θ's share of the variance is about 0.4%.

## Reproducing

```
dev/scripts/corea_budget_kim96.sh <new outdir>   # from a clean worktree at a pushed commit
```

**Check:** the script was rerun on the login node, on copies of the three archives' inputs
and their `SHA256SUMS` (KIM92_DIR, KIM95_DIR and KIM42_DIR point it at the copies). It exited
0 and wrote a `budget.md` byte-identical to the job's (`cmp`).
