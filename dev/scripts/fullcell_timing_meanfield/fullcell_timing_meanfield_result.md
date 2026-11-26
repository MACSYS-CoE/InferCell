# KIM-184: a gradient of the full cell's first-order mean field: 1.3 s with a sparse dual solve, about 18 min with LinearSolve 3.73.0's

**Two numbers, and which one holds depends on the linear solve.** A 20-parameter
dual gradient over the fit window costs:

- **1.28 s serial and 0.62 s on 4 threads** with `DualSparseKLU`, a custom sparse dual
  linear solve that lives in this script (jobs 21666252 and 21666253);
- **about 1,108 s (18.5 min)** with the stock path in **LinearSolve 3.73.0**, its
  ForwardDiff extension's `DualLinearCache`, behind `Rodas5P(linsolve = KLUFactorization())`.
  This is 2 × 554.1 s, from one measured pass over the whole window (job 21667143),
  860× the sparse solve's pass. The second pass was not run on the stock path.

On the one pass compared, the two solutions are bitwise equal. That pass is partials
1–10 at B = 20 and θ*, at the 52 saves; a 30 s pass in the check stage agrees too.
Equality is expected here: each ∂W/∂θ_b is non-zero in one column only, so every row
of ∂W·x₀ is a single product. A state-dependent Jacobian will need a tolerance instead.

**The 860× is specific to LinearSolve 3.73.0**, the version InferCell's Manifest
resolves. LinearSolve 3.87.0, which is in the local depot, adds `SparseMatrixCSC`
methods for the partials and separate A and b cache flags. **3.87 is untimed here.**
- KIM-194 should time this script's `library` stage under 3.87 before choosing
  between a version bump and `DualSparseKLU`.
- That run should keep `KLUFactorization` explicit. 3.87's default sparse choice is
  `PureKLUFactorization`, which 3.87 sends down a different (direct-dual) path.
- Moving to 3.87 needs a re-resolve. SciMLLogging must be at least 1.10.1, and the
  Manifest has 1.9.1.

`DualSparseKLU` relies on two method definitions for its own `LinearSolveFunction`
type:
- `SciMLBase.init`, which keeps the dual problem off `DualLinearCache`, as LinearSolve
  does for `SparspakFactorization`;
- `LinearSolve.needs_concrete_A`.

It also reaches into three LinearSolve internals:
- `Base.get_extension(LinearSolve, :LinearSolveSparseArraysExt).KLU`, the vendored KLU;
- `LinearSolveForwardDiffExt.DualAbstractLinearProblem`;
- `LinearSolve.__init`.

Anything that adopts the 1.3 s figure (KIM-187, KIM-194) inherits that code, or
replaces it with a LinearSolve version whose dual path is measured to be sparse.

**This is a lower bound on cost, and it says nothing about accuracy.** The system is
the full syn3A CME table's 2,002 first-order reactions alone. Those make it linear,
with constant rate constants. The full cell's real mean field also has:

- the 363 second-order reactions;
- the 175 metabolic ODE fluxes, with non-linear or zero-order rate laws;
- the rate constants rebuilt from the pools every minute.

The second-order reactions and most metabolic rate laws make the Jacobian depend on
the state, and the rebuilds make the rates jump each minute. None is in this timing.
The per-step costs below are floors for this system, not constants to rescale by.

## Runs

All runs use `dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh`,
from this worktree at a pushed commit. They ran on milan (AMD EPYC 7543 32-Core),
Julia 1.10.5, `OPENBLAS_NUM_THREADS=1`. The outdir is
`/fred/oz022/tkimpson/infercell-kim184-2026-11-26/`. It holds:
- `sacct.tsv`: all 14 jobs, with submit, start and end times, node, CPUs and elapsed time;
- `SHA256SUMS`: every file in the outdir (35 files; the sum file's own sha256 is
  `1cf773f0…`).

| what | jobs | commit | outdir |
|---|---|---|---|
| **the run of record**: check, time (serial), threaded (4 threads, `--gcthreads=4,1`) | 21666251, 21666252, 21666253 (precompile 21666250) | 74755da | `run1/` |
| the stock dual path over the whole window (`library` stage) | 21667143 (precompile 21667142) | 6a20e2c | `library/` |
| smoke 2, 300 s window, 2 repeats (plumbing only) | 21664745–48 | 74755da | `smoke2/` |
| smoke 1, which found the two configuration faults below | 21663390–93, submitted at 502aea4 | see note | `smoke1/` |

**Smoke 1's commit.** Its check and time jobs ran the faulty code at 502aea4. Its
threaded job (21663393) started after 74755da was committed, so it ran the fixed
code: 0.5 s passes, and |g| ≈ 61 against 3e-27 in the time job. The tree moved under
a queued job. Smoke 1 is plumbing only, and none of its numbers is used except the
437 s that found fault 1.

| job | node |
|---|---|
| check 21666251 | dave5 |
| time 21666252 | dave28 |
| threaded 21666253 | dave64 |
| library 21667143 | dave1 |

**The stage outputs are committed** in `results/`: `check.md`, `time.md` and
`threaded.md` from `run1/`, and `library.md` from `library/`. They are byte-identical
to the outdir copies.

**The committed driver is not quite the code that ran.** It differs from 74755da and
6a20e2c in two ways:
- guards: every solve's retcode is checked, and the logged bitwise comparisons and
  the linear-cache path are asserted, where the runs only logged them;
- comments and output strings.

No number changes. In the runs of record, every logged comparison was exactly 0, and
the linear-cache path was `LinearSolveFunction`.

**Compute used:** about 2.9 CPU-h over every job, the library stage's 0.26 CPU-h
included. Most of it is smoke 1's check and time jobs, which ran on the faulty
configuration until they hit their 1 h limit.

**Inputs:** the syn3a-ensemble repo, `/fred/oz022/tkimpson/syn3a-ensemble`, at commit
`2075809cbb825140cb31a7bec8bc3529770c8e13`, read only.

| table | sha256 |
|---|---|
| `data/model/cme/reactions.csv` | `63fbeca180fe8bdebec646df0ddbffb0149958fc2ff9fc3e6ec8803cd24c6aaf` |
| `data/model/cme/species.csv` | `c1acb80878e8a500a86e0e1d9e195a491b46a01c337788c0a6554df19cbe6ec9` |

The ODE tables (`data/model/ode/`) are not used. Every ODE species is also a CME
species. The ODE's 175 fluxes use 168 distinct rate laws (the `rate_form` column):
bespoke reversible Michaelis–Menten forms, `zeroOrderOnOff`, the `FAFBA` constant
uptakes and `RateForm0–5`, per `inventory.md`. None is an LM first-order reaction, so
they fall outside "first-order reactions only".

**Environment:** `dev/scripts/fullcell_timing_meanfield/Project.toml`. Its Manifest is
the root Manifest's, resolved offline. `[compat]` pins OrdinaryDiffEq 6.107.0,
LinearSolve 3.73.0, ForwardDiff 0.10.39 and SciMLBase 2.153.1, the versions InferCell
itself runs. Offline `Pkg.instantiate()` still works with the pins; it warns that the
project hash has changed.

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

The rate constants are the table's `rate` column: the first minute's values, in s⁻¹.
In the reference, 1,434 of the 2,365 rates change at each minute boundary
(`inventory.md`); here all are held at their first-minute values. Under these, the
mean field is dx/dt = S·(k ∘ x[r]). Its Jacobian is S·diag(k)·E, which is constant in
x and depends on θ only through k. It is written analytically into the sparse
prototype:
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

**PASS** (job 21666251; `results/check.md`).

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
- The full list is in `results/check.md`.

This is a plausible block for cost: each parameter moves the states it feeds over the
whole window. The scientific block is KIM-187's to pick, by doubling time.

**Does the block move the cost?** Here it did not. Every pass, in every repeat, at
B = 10, 20 and 40, serial and threaded, took 144 steps. The serial per-pass medians
agree within 4%. That is an observation, not a guarantee: the step count depends on
the partials through the error norm (below), so another block could take a different
number of steps.

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

The two checks of `DualSparseKLU` that are independent of it are:
- these Float64 central differences;
- the bitwise comparison with LinearSolve's `DualLinearCache` (above, and fault 1
  below).

Two more equalities hold at the timing tolerances, but both compare `DualSparseKLU`
with itself:
- The pass-wise gradient equals `ForwardDiff.gradient` with `Chunk{10}` bitwise
  (max |Δ| = 0). Both go through `DualSparseKLU`, so this checks only how the passes
  are seeded and assembled.
- The threaded gradient equals the serial one bitwise at B = 10, 20 and 40.

The gradient at the timing tolerances is not checked against differences. That is
fine for a cost measurement, but KIM-187 should know it.

## Timings (checklist items 3–5)

Each figure is the median of 5 repeats after a warm-up. Every repeat is in
`results/time.md` and `results/threaded.md`.

| B | passes | primal solve, serial | gradient, serial (21666252) | gradient, 4 threads (21666253) | speed-up, same node | speed-up, across jobs |
|---|---|---|---|---|---|---|
| 10 | 1 | 0.110 s | 0.649 s | 0.581 s | 1 (the baseline) | 1.12× |
| **20** | 2 | **0.109 s** | **1.277 s** | **0.620 s** | **1.87×** | 2.06× |
| 40 | 4 | 0.109 s | 2.551 s | 0.760 s | 3.06× | 3.36× |

The threaded job's own split:

| B | primal solve | per pass, median | per pass, mean | gradient (wall) |
|---|---|---|---|---|
| 10 | 0.107 s | 0.581 s | 0.591 s | 0.581 s |
| 20 | 0.101 s | 0.603 s | 0.607 s | 0.620 s |
| 40 | 0.102 s | 0.680 s | 0.701 s | 0.760 s |

**Read the same-node column.** The serial and threaded jobs ran on different nodes
(dave28 and dave64), and with different GC settings (`--gcthreads=4,1` against the
default). B = 10 has one pass and so no parallel work, yet it is 12% faster in the
threaded job. That gap is most likely the node and GC difference; it was not
isolated. The across-jobs speed-up at B = 20 (2.06×) even exceeds the 2× that two
passes allow. The same-node column takes the threaded job's own one-pass time
(B = 10, 0.581 s) as the serial cost of a pass. Its speed-up is the number of passes
× 0.581 s, divided by the wall time.

**Threads contend.** With 4 passes on 4 threads (B = 40), the mean pass is 21% slower
than a lone pass (0.701 against 0.581 s), and the wall time is 31% longer. Scaling past
4 threads, or several chains on one node, is unmeasured.

**The split, at B = 20, serial.**
- The primal solve takes 0.109 s, over 174 steps.
- The gradient is 2 dual passes of 0.636 s each, over 144 steps each.
- So one pass of 10 partials costs 5.8× a primal solve. That ratio is **not at matched
  accuracy**, because the dual solve takes fewer steps (see below).
- **Per step, the like-for-like figure:** 4.4 ms against 0.63 ms, or **7×**.

**Why the dual solve takes fewer steps.** The dual solve took 144 steps against the
primal's 174, none rejected, at the same tolerances. The reviewer of PR #207
identified the mechanism (the reviewer's arithmetic, not re-derived here):
- DiffEqBase's error norm for dual states is sqrt(Σ sse / (n(1+N))). Its denominator
  counts the N = 10 partials as well as the values.
- Small errors in the partials therefore dilute the RMS by about √11.
- That predicts 137–142 steps, against the 144 observed. Smoke 2 shows the same ratio
  (129/107 = 1.21).

At matched accuracy (174 steps) a pass would take about 0.77 s, and a B = 20 serial
gradient about **1.54 s**. That is the reviewer's estimate; **1.277 s is the measured
figure**. The step count also depends on the chunk size and the block.

**Block size (item 5).** Serial, B = 10 : 20 : 40 costs 0.649 : 1.277 : 2.551 s, a ratio
of **1 : 1.97 : 3.93**. That is linear in the number of passes: each pass costs about
0.64 s whatever B is. Threaded with one pass per thread, B = 40 costs 1.3× B = 10, at
up to 4 threads.

**Days per chain.** These are for the first-order linear surrogate, which has a
constant Jacobian. They are a lower bound on both the cost per step and the step
count. The assumptions:
- 8,600 gradients per chain, the reduced model's chain length, from the issue;
- a constant cost per gradient, as measured at θ*;
- no likelihood, no warm-up and no primal-only evaluations;
- at B = 20 only 2 of the 4 threads have work.

| B | serial | 4 threads |
|---|---|---|
| 10 | 1.55 h (0.065 d) | 1.39 h (0.058 d) |
| **20** | **3.05 h (0.127 d)** | **1.48 h (0.062 d)** |
| 40 | 6.09 h (0.254 d) | 1.82 h (0.076 d) |

At the matched-accuracy estimate of 1.54 s, the B = 20 serial figure would be about
3.7 h.

## Two configuration faults, fixed before the run of record

Smoke 1 found both. Each was fixed in 74755da, before run1 was submitted.

1. **LinearSolve 3.73.0's ForwardDiff extension densifies a sparse dual W.**
   - Under a dual state, Rodas5P hands a `SparseMatrixCSC{Dual}` W to
     `DualLinearCache`. That code solves the value part, then the partials against
     the same factorisation. The algebra is right.
   - But 3.73.0 has only the generic dense `partials_to_list` and
     `update_partials_list!`. Each partial of W becomes a dense n×n matrix, filled by
     sparse `getindex` over all n² entries, and is multiplied densely.
   - The recurring cost is `update_partials_list!`. One `rhs_cache_valid` flag is
     cleared by both `setA!` and `setb!`, and `setb!` runs at every stage solve. So the
     dense rebuild happens at all 8 stage solves of each Rodas5P step, not once per W.
     This diagnosis is from reading the source, not from a profile.
   - **The check stage measured it on one pass of 10 partials over 30 s:** 297.1 s
     against 0.353 s, over the same 82 steps. The two agree bitwise (job 21666251).
   - **Over the whole window** (the `library` stage, job 21667143, one run): one pass took
     **554.1 s against 0.645 s** (median of 5), over the same 144 steps. That is
     3.85 s per step against 4.5 ms. The pass is again bitwise equal.
   - That is about 18.5 min per 20-parameter gradient (2 passes), or about 110 days
     per chain at 8,600 gradients, serial.
   - `DualSparseKLU` (in the script) does the same algebra on the sparse pattern:
     - KLU on the value part, refactored only when the values change;
     - R = ∂b − ∂W·x₀, formed from the non-zeros;
     - one multi-right-hand-side solve.
   - It is wired in as a `LinearSolveFunction`, kept off `DualLinearCache` as LinearSolve
     does for `SparspakFactorization`, with `needs_concrete_A = true`.
   - This is LinearSolve 3.73.0's cost on a sparse dual solve. 3.87.0 may not have it
     (see the top of this note), so the choice between a version bump and
     `DualSparseKLU` belongs in KIM-194's options, after 3.87 is timed.
2. **The first ranking used a 10 s trapezoid.** It put 22 tRNA-charging drains first,
   with an "integrated flux" of 150 each. Each is one molecule that drains in about
   0.03 s at k = 30 s⁻¹, so its true flux is 1. Against that block the misfit hardly
   moved (|g| ≈ 3e-27). The ranking now integrates on the solver's steps.

## What this timing does not cover

- **Second-order reactions (363), metabolism (175 fluxes with non-linear or zero-order
  rate laws), and the per-minute rate rebuilds.** With these in, the Jacobian depends
  on the state. That means:
  - the analytic constant-pattern Jacobian would have to come from AD or from
    symbolic code;
  - the cost per step rises. The Jacobian becomes state-dependent, second-order terms
    add fill-in, and fewer states sit at zero. Here 501 states are zero at the start
    and at every save;
  - the step count is unmeasured for the full system. Two hints that it rises: the
    reference rebuilds its ODE every 1 s, and the steps here are front-loaded (129
    primal steps in the first 300 s, smoke 2, against 174 over all 3,150 s). How far
    it rises is not known.

  The per-step costs here are 0.63 ms per primal step and 4.4 ms per dual step of 10
  partials. They are floors for this system (n = 2,112, nnz = 21,095, a free
  Jacobian). The full cell raises both the cost per step and the step count. The
  second-order terms would also make the mean field a closure, not exact.
- **Accuracy.** Nothing here says the first-order mean field approximates the full
  cell. It tests cost only.
- **Discontinuities.** Each discontinuity needs `tstops` or a callback, and each
  restarts the stiff solver. How many there are depends on which ones a mean field
  keeps:
  - **52 in the window** if only the CME's minute rebuilds count;
  - **about 3,150** if the reference's 1 s ODE hook is kept.

## Reproduce

```bash
# from a clean worktree at a pushed commit
dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh <outdir>          # check, time, threaded
dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh <outdir> library  # the LinearSolve 3.73.0 dual path
SMOKE=1 dev/scripts/fullcell_timing_meanfield/fullcell_timing_meanfield.sh <outdir>  # plumbing pass
```
