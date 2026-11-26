# KIM-184: a gradient of the full cell's first-order mean field costs about 1.3 s

**This is a lower bound on cost, and it says nothing about accuracy.** The system is
the full syn3A CME table's 2,002 first-order reactions alone. Those make it linear,
with constant rate constants. The full cell's real mean field also has:

- the 363 second-order reactions;
- the 175 non-linear metabolic ODE fluxes;
- the rate constants rebuilt from the pools every minute.

Each of those makes the Jacobian depend on the state. Each is likely to add solver
steps, and none is in this timing. Per-step costs are given below so that KIM-187
and KIM-194 can rescale.

## Runs

All runs use `dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh`,
from this worktree at a pushed commit. They ran on milan (AMD EPYC 7543 32-Core),
Julia 1.10.5, `OPENBLAS_NUM_THREADS=1`.

| what | jobs | commit | outdir (under `/fred/oz022/tkimpson/infercell-kim184-2026-11-26/`) |
|---|---|---|---|
| **the run of record**: check, time (serial), threaded (4 threads, `--gcthreads=4,1`) | 21666251, 21666252, 21666253 (precompile 21666250) | 74755da | `run1/` (`check.md`, `time.md`, `threaded.md`, `.jls`, `logs/`) |
| smoke 2, 300 s window, 2 repeats (plumbing only) | 21664745–48 | 74755da | `smoke2/` |
| smoke 1, which found the two configuration faults below | 21663390–93 | 502aea4 | `smoke1/` |

| job | node |
|---|---|
| check 21666251 | dave5 |
| time 21666252 | dave28 |
| threaded 21666253 | dave64 |

**Compute used:** 2.6 CPU-h for every job of the three rounds. Most of it is smoke 1's
two jobs, which ran on the faulty configuration until they hit their 1 h limit.

**Inputs:** `/fred/oz022/tkimpson/syn3a-ensemble` at commit
`2075809cbb825140cb31a7bec8bc3529770c8e13`, read only.

| table | sha256 |
|---|---|
| `data/model/cme/reactions.csv` | `63fbeca180fe8bdebec646df0ddbffb0149958fc2ff9fc3e6ec8803cd24c6aaf` |
| `data/model/cme/species.csv` | `c1acb80878e8a500a86e0e1d9e195a491b46a01c337788c0a6554df19cbe6ec9` |

The ODE tables (`data/model/ode/`) are not used. Every ODE species is also a CME
species. The ODE's 175 fluxes are Michaelis–Menten forms, not first-order reactions,
so they fall outside "first-order reactions only".

**Environment:** `dev/scripts/fullcell_timing_meanfield/Project.toml`. Its Manifest is
the root Manifest's, resolved offline. It uses OrdinaryDiffEq 6.107.0, LinearSolve
3.73.0, ForwardDiff 0.10.39 and SciMLBase 2.153.1, the versions InferCell itself runs.

## The system (checklist item 1)

| quantity | value |
|---|---|
| CME species in the table | 2,661 |
| **states**: the species a first-order reaction reads or changes | **2,112** |
| **first-order reactions** (LM type 1, propensity k·x) | **2,002** |
| second-order reactions (type 2), left out | 363 |
| **Jacobian non-zeros** from the table: distinct (changed species, reactant) pairs | **20,044** |
| … with a structural diagonal added, as stored | 21,095 (0.47% dense) |
| reactions with non-zero flux over the window | 1,511 |

**The first-order count is 2,002, as the issue and `inventory.md` say.** The issue's
"about 2,800 states" is not what the table gives. The table has 2,661 species. Of
those, 549 take part in no first-order reaction and are left out, since they are
constant here. The 2,800 may have been 2,661 CME + 146 ODE species = 2,807, but the
146 ODE species are all CME species already.

The rate constants are the table's minute-1 values (`rate`). Under these, the mean
field is dx/dt = S·(k ∘ x[r]). Its Jacobian is S·diag(k)·E, which is constant in x and
depends on θ only through k. It is written analytically into the sparse prototype:
- no sparse AD;
- no colouring;
- `tgrad` is zero, because the system is autonomous.

**Check: the RHS at the initial state.** The RHS at procid 1's initial counts was
checked against the table's propensities, Σ_j ν_ij k_j x_j. The sum was formed
independently: in BigFloat at 256 bits, by species name, straight from the CSV rows.
- 677 states are non-zero in the reference. Their largest relative error is **2.6e-15**
  (bound 1e-10).
- 1,435 are exactly zero in the reference, and every one is exactly zero in the RHS.
- J·x₀ matches the RHS to 1.5e-16 relative.

**PASS** (job 21666251).

## The setup

The scalar and the solver follow route (ii) on the reduced model (`src/mean_field_route.jl`,
KIM-99), where they carry over.

| | |
|---|---|
| solver | `Rodas5P`, as `mean_field_solver` uses |
| tolerances | abstol 1e-10, reltol 1e-8, maxiters 10⁵: `solve_mean_field`'s pinned values |
| linear solver, primal | LinearSolve's `KLUFactorization` (symbolic factorisation reused) |
| linear solver, dual | `DualSparseKLU` (below): KLU on the value part, with all partials as one n×N right-hand side |
| sparsity | analytic Jacobian into a `SparseMatrixCSC` prototype, 21,095 non-zeros, `FullSpecialize` |
| fit window | 0 to `COREA_FIT_S` = 3,150 s, saved at 60, 120, …, 3,120 s (52 saves), as in KIM-92 and KIM-95 |
| parameters | θ_b = log(k_b / k_b⁰), the B reactions with the largest flux integrated over the window at θ = 0, ties broken by reaction index. The blocks of 10, 20 and 40 are nested. |
| evaluated at | θ*_b = 0.05·(−1)^b |
| scalar | Gaussian misfit over all 2,112 states at the 52 saves, against the θ = 0 solution, with variance max(d, 1) |
| gradient | forward mode, chunk 10 (`ForwardDiff.pickchunksize` at B = 10, 20 and 40), so 1, 2 and 4 passes. One pass is one dual solve. |
| threaded | the passes run one per thread (`Threads.@threads :static`), 4 threads |

**How the block was chosen.** The block is the 40 reactions with the largest
integrated flux. The integral is 3-point Gauss–Legendre on each of the solver's steps,
through its dense output.
- 39 of the 40 are translation reactions (k ≈ 0.06–0.12 s⁻¹). Each reads an mRNA
  without consuming it, and makes the gene's product plus the cost counters.
- Rank 20 is reaction 25, the translocation `ptsg_cyto` → `ptsg`.
- The list is in `run1/check.md`.

This is a plausible block for cost: each parameter moves the states it feeds over the
whole window. Which block it is should not move the cost much, because every pass
solves the same system. The scientific block is KIM-187's to pick, by doubling time.

## The dual gradient against central differences (checklist item 2)

The setup is as AV3 (`test/test_mixed_likelihood.jl`): B = 20, tight tolerances
(abstol 1e-14, reltol 1e-12) and h = 1e-4 in log k. The three parameters are ranks
1, 10 and 20.

| θ | reaction | dual | central FD | relative |
|---|---|---|---|---|
| θ_1 | 1164 | −2.607077870429e+04 | −2.607078230273e+04 | 1.4e-07 |
| θ_10 | 699 | −9.774521756389e+03 | −9.774522260676e+03 | 5.2e-08 |
| θ_20 | 25 | −1.601730910017e+01 | −1.601731171263e+01 | 1.6e-07 |

**PASS** at 1e-6 relative (job 21666251).

Two more equality checks hold at the timing tolerances:
- The pass-wise gradient equals `ForwardDiff.gradient` with `Chunk{10}` bitwise
  (max |Δ| = 0).
- The threaded gradient equals the serial one bitwise at B = 10, 20 and 40.

## Timings (checklist items 3–5)

Each figure is the median of 5 repeats after a warm-up. Each repeat is listed in
`run1/time.md` and `run1/threaded.md`.

| B | passes | primal solve | gradient, serial (21666252) | gradient, 4 threads (21666253) | threaded speed-up |
|---|---|---|---|---|---|
| 10 | 1 | 0.110 s | 0.649 s | 0.581 s | 1.12× (one pass: no parallel work, a different node) |
| **20** | 2 | **0.109 s** | **1.277 s** | **0.620 s** | **2.06×** |
| 40 | 4 | 0.109 s | 2.551 s | 0.760 s | 3.36× |

**The split, at B = 20, serial.**
- The primal solve takes 0.109 s, over 174 steps.
- The gradient is 2 dual passes of 0.636 s each, over 144 steps each.
- So one pass of 10 partials costs 5.8× a primal solve.
- Per step, that is 4.4 ms against 0.63 ms, or 7×.
- The dual solve takes fewer steps than the primal (144 against 174, none rejected) at
  the same tolerances. The partials enter the step-size control.

**Block size (item 5).** Serial, B = 10 : 20 : 40 costs 0.649 : 1.277 : 2.551 s, a ratio
of **1 : 1.97 : 3.93**. That is linear in the number of passes: each pass costs the same
0.64 s whatever B is. Threaded with one pass per thread, B = 40 costs 1.3× B = 10. The
block's size is then mostly a matter of cores, up to one thread per 10 parameters.

**Days per chain**, at 8,600 gradients (the reduced model's chain length, from the
issue):

| B | serial | 4 threads |
|---|---|---|
| 10 | 1.55 h (0.065 d) | 1.39 h (0.058 d) |
| **20** | **3.05 h (0.127 d)** | **1.48 h (0.062 d)** |
| 40 | 6.09 h (0.254 d) | 1.82 h (0.076 d) |

That counts only the forward model's gradient. A likelihood on top, such as
transcript factors, adds to it.

## Two configuration faults, fixed before any timing was read

1. **LinearSolve's ForwardDiff extension densifies a sparse dual W.**
   - Under a dual state, Rodas5P hands a `SparseMatrixCSC{Dual}` W to
     `LinearSolveForwardDiffExt`'s `DualLinearCache`. That code solves the value part,
     then the partials against the same factorisation. The algebra is right.
   - But `partials_to_list` stores each partial of W as a dense n×n matrix. It fills
     that matrix by sparse `getindex` over all n² entries, and then multiplies by it
     densely.
   - On this system, the check stage measured it on one pass of 10 partials over 30 s:
     **297.1 s against 0.353 s, over the same 82 steps.** The two agree bitwise, both
     in partials and in values (job 21666251).
   - In smoke 1 (job 21663392, 300 s window), the library path took 437 s per pass
     against a 0.08 s primal.
   - **Extrapolation, not measured:** at the run's 144 steps per pass and 3.6 s per step,
     the library path would take about 520 s per pass. That is about 17 min per
     20-parameter gradient, or about 100 days per chain.
   - `DualSparseKLU` (in the script) does the same algebra on the sparse pattern:
     - KLU on the value part, refactored only when the values change;
     - R = ∂b − ∂W·x₀, formed from the non-zeros;
     - one multi-right-hand-side solve.
   - It is wired in as a `LinearSolveFunction`, kept off `DualLinearCache` as LinearSolve
     does for `SparspakFactorization`, with `needs_concrete_A = true`.
   - This is the default path's cost on any sparse full-cell solve under duals, so it
     belongs in KIM-194's options.
2. **The first ranking used a 10 s trapezoid.** It put 22 tRNA-charging drains first,
   with an "integrated flux" of 150 each. Each is one molecule that drains in about
   0.03 s at k = 30 s⁻¹, so its true flux is 1. Against that block the misfit hardly
   moved (|g| ≈ 3e-27). The ranking now integrates on the solver's steps.

## What this timing does not cover

- **Second-order reactions (363), metabolism (175 non-linear fluxes), and the per-minute
  rate rebuilds.** With these in, the Jacobian depends on the state. That means:
  - the analytic constant-pattern Jacobian would have to come from AD or from symbolic
    code;
  - the system's stiffness likely rises (metabolism runs on faster time scales than
    gene expression);
  - the step count, here 174 primal and 144 dual over 3,150 s, could rise by orders of
    magnitude.

  For rescaling, use 0.63 ms per primal step and 4.4 ms per dual step of 10 partials.
  The second-order terms would also make the mean field a closure, not exact.
- **Accuracy.** Nothing here says the first-order mean field approximates the full
  cell. It tests cost only.
- **Discontinuities.** A model with minute boundaries needs `tstops` or callbacks. Each
  restarts the stiff solver, and there are 52 of them in the window.

## Reproduce

```bash
# from a clean worktree at a pushed commit
dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh <outdir>        # all stages
SMOKE=1 dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh <outdir> # plumbing pass
```
