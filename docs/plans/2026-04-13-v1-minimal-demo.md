# v1 Minimal Demo Plan

> Date: 2026-04-13

## Goal

Demonstrate the core InferCell value proposition end-to-end: **hybrid Bayesian inference across heterogeneous blocks (NUTS + ABC-SMC) with uncertainty propagation through a boundary protocol**, all via a single composable API in pure Julia.

## Current state

Phase 1 (TX/TL ODE + NUTS) and Phase 2 (Stochastic GE + ABC-SMC) are complete and merged to `main`. The unified `infer()` API dispatches correctly based on model formalism. Two independent inference blocks work; they are not yet connected.

| What works | What's missing |
|-----------|---------------|
| ODE model (TX/TL) with NUTS parameter recovery | Multi-module ODE composition |
| SSA model (Stochastic GE) with ABC-SMC | Level 1 boundary protocol |
| Unified `infer()` dispatch | Shared parameters across modules |
| Posterior predictive for both backends | End-to-end pipeline figure |
| Structural identifiability check | |

## Implementation plan

### Step 1 — Light metabolism module + shared parameters

**What:** A minimal ODE metabolism module (3-5 reactions) that couples to TX/TL through shared resource pools (ATP, nucleotides, amino acids).

**Biology (kept deliberately simple):**
- Energy source → ATP (production)
- ATP + precursors → nucleotides (feeds transcription)
- ATP + precursors → amino acids (feeds translation)
- ATP consumption by TX/TL (coupling point)

**Parameters (~5 free):** metabolic rate constants, pool sizes. Shared parameters with TX/TL (e.g. energy cost per transcription event) exercise the composition machinery.

**Engineering:**
- New file: `src/models/light_metabolism.jl`
- Solve open design question #1: shared parameters across modules. Current `param_idxs` is `UnitRange{Int}` (contiguous, non-overlapping). Options:
  - **(a)** Change to `Vector{Int}` — modules can reference the same slot in the flat parameter vector
  - **(b)** Aliasing at orchestrator level — one slot, multiple modules get the same index via `input_map`-like mechanism
  - Option (a) is simpler and sufficient for v1
- Update `build_problem` to handle shared parameter indices
- Tests: unit tests for the new model, composition tests for `build_problem([txl, metabolism])`

**Done when:** `build_problem([TranscriptionTranslation(), LightMetabolism()])` returns a single `ODEProblem` with ~4 states and ~8-10 free parameters. Forward simulation produces biologically reasonable trajectories (ATP consumed, mRNA/protein produced).

### Step 2 — Joint ODE inference

**What:** NUTS inference on the composed TX/TL + metabolism system.

**Engineering:**
- The existing `_infer_nuts` path should work with minimal changes if step 1 is done correctly — the Turing model already reads `InferParameter` metadata generically
- Validate that `ForwardDiffSensitivity` scales acceptably at ~10 parameters (benchmark wall-clock)
- Synthetic twin experiment: generate data from known parameters, recover posteriors

**Done when:** Joint NUTS recovers all ~8-10 parameters at 90% CI from synthetic data. Wall-clock is under 10 minutes on a single HPC node.

### Step 3 — Level 1 boundary protocol (sequential conditioning)

**What:** The architectural centrepiece. Infer the ODE block → take posterior samples → condition the SSA block on those posteriors → ABC-SMC for remaining stochastic parameters.

**Design:**
```
ODE block (TX/TL + metabolism)
    ↓ NUTS
posterior samples for shared rates (k_tx, k_tl, γ_mRNA, γ_protein)
    ↓ sequential conditioning
SSA block (Stochastic GE)
    ↓ ABC-SMC with informed priors from ODE posterior
posterior samples for all parameters
```

Information flows one direction only. This produces a **cut posterior** — the SSA block does not feed back into the ODE block. Valid when blocks are conditionally independent given shared boundary variables.

**Engineering:**
- New function: `boundary_condition(ode_chain, ssa_model; method=:sequential)` or extend `infer()` to accept a prior source
- The key mechanism: replace ABC-SMC priors for shared parameters with KDE or empirical distribution fitted to ODE posterior samples
- Alternatively: run ABC-SMC once per ODE posterior draw (embarrassingly parallel, conceptually simpler, more expensive)
- New file or extend: `src/boundary.jl`
- Tests: verify that SSA posteriors tighten when conditioned on ODE posteriors vs uninformed priors

**Done when:** The full pipeline runs: ODE inference → boundary conditioning → SSA inference. Posteriors for shared parameters are consistent across blocks. Conditioning demonstrably improves SSA parameter recovery compared to uninformed priors.

### Step 4 — End-to-end demo + figure

**What:** One script that runs the full v1 pipeline and produces a publication-quality multi-panel figure.

**Figure layout (single figure, 4 panels):**
- **(A)** Architecture schematic — ODE block → boundary → SSA block (hand-drawn or TikZ, not generated)
- **(B)** Joint ODE posteriors — TX/TL + metabolism parameters recovered, ground truth overlaid
- **(C)** Boundary protocol effect — SSA posteriors with vs without ODE conditioning, showing uncertainty propagation
- **(D)** Posterior predictive — forward simulations from the full pipeline posterior, data overlaid

**Engineering:**
- New file: `examples/v1_demo.jl`
- Slurm submission script for HPC
- Figure output to `examples/figures/`

**Done when:** `v1_demo.jl` runs end-to-end on one HPC node and produces a single figure that tells the complete InferCell story.

## Tier 2 — Strong but deferrable

These strengthen the narrative but aren't required for the minimal demo.

### Bursting gene expression (two-state promoter)

Two-state promoter model: gene switches between ON (transcribing) and OFF states stochastically. mRNA is produced in bursts during ON periods. This has no ODE equivalent — the dynamics are inherently discrete — making it the strongest argument for why likelihood-free inference is necessary alongside gradient-based methods.

**Parameters:** k_on, k_off (switching rates), k_tx_burst (transcription rate when ON), γ_mRNA, k_tl, γ_protein.

Would replace or supplement StochasticGE in the SSA block. Enables the model selection figure (constitutive vs bursty, Bayes factors).

### Model selection via Bayes factors

Compare constitutive transcription (current StochasticGE) vs bursty transcription (two-state promoter) on the same synthetic data. Compute Bayes factors from ABC-SMC evidence estimates or posterior odds.

## AD backend: ForwardDiff → Mooncake migration

The SciML ecosystem is migrating from ForwardDiff.jl to [Mooncake.jl](https://github.com/compintell/Mooncake.jl) as the default AD backend. Mooncake is a reverse-mode AD that should scale better for higher parameter counts (O(1) gradient cost vs ForwardDiff's O(n_params)). This is relevant for step 2 (joint ODE inference at ~10 params) and future scaling.

**Action items (separate PR):**
- Check Turing.jl + SciMLSensitivity.jl Mooncake support status
- Benchmark ForwardDiff vs Mooncake on the joint TX/TL + metabolism system
- If Mooncake is production-ready, switch `sensealg` default
- If not, document the migration path and pin to ForwardDiff for now

## Risks

| Risk | Mitigation |
|------|-----------|
| Shared parameter design blocks composition | Start with option (a) — `Vector{Int}` indices. Simple, sufficient for v1 |
| ForwardDiff too slow at 10 params | Benchmark early in step 2. Fall back to adjoint methods or Mooncake |
| Boundary protocol is conceptually unclear | Start with the simplest version: replace ABC-SMC priors with ODE posterior samples. Iterate |
| ABC-SMC too slow for conditioned inference | Reduce particle count, parallelise across ODE posterior draws |
| Joint figure is hard to produce programmatically | Panel A is hand-drawn. Panels B-D are CairoMakie. Assemble in a layout |

## Non-goals for v1

- Level 2/3 boundary protocols
- Real experimental data
- Spatial dynamics, genome-scale metabolism
- Population-level inference
- Moment closure methods
- Neural posterior estimation (future SBI upgrade)
