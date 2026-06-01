# Phase 1 Implementation Log

**Date:** 2026-04-01
**Branch:** `phase-1`
**Commit:** `6815f29`

## What was done

Implemented the full Phase 1 architecture: a transcription-translation ODE model with end-to-end Bayesian parameter inference via NUTS/Turing.jl.

### Files created

**Package source (`src/`):**
- `InferCell.jl` — module entry point
- `parameters.jl` — `InferParameter` struct and filtering utilities
- `interface.jl` — `AbstractSubModel` abstract type, `SubModelContext`, formalism validation
- `likelihoods.jl` — `ObservedData` and `PosteriorPredictive` data types
- `orchestrator.jl` — `build_problem`: composes sub-models into a single `ODEProblem`
- `models/transcription_translation.jl` — TX/TL sub-model (2 states, 5 free parameters)
- `inference.jl` — generic Turing `@model`, `infer`, `observe`, `posterior_predictive`, `check_identifiability`

**Tests (`test/`):**
- `test_parameters.jl` — InferParameter construction and filtering (14 tests)
- `test_orchestrator.jl` — build_problem, forward solve, formalism validation (16 tests)
- `test_inference.jl` — observe, build_turing_model, identifiability check (15 tests)
- `test_txl.jl` — integration test: 1000 NUTS samples, parameter recovery, posterior predictive (8 tests, gated behind `INFERCELL_INTEGRATION_TESTS=true`)
- `runtests.jl` — test runner

**Demo (`examples/`):**
- `phase1_demo.jl` — end-to-end workflow with CairoMakie diagnostic plots

### Key design decisions

1. **Generic Turing model via metadata.** The `@model` function samples parameters in a loop (`theta[i] ~ priors[i]`) rather than hardcoding names. Chain columns are renamed from `theta[1]` to `:k_tx` etc. after sampling via `MCMCChains.replacenames`.

2. **`u0` is SVector, `p0` is Vector.** Initial conditions use `SVector` (matching the `SA[...]` return type of `dynamics`). Parameters use plain `Vector{Float64}` because `remake(prob, p=...)` must accept ForwardDiff `Dual` numbers during AD — SVector conversion from abstract-typed containers triggers `Float64(::Dual)` errors.

3. **Comprehension `[theta[i] for i in 1:n_ode]` in hot path.** Direct slicing `theta[1:n_ode]` returns a `SubArray{Real}` that `remake` tries to convert to `Vector{Float64}`, failing with Dual numbers. The comprehension creates a new vector with the correct element type.

4. **Julia 1.10.5 via module system.** The cluster has Julia available via `module load julia/1.10.5`. The Manifest.toml was regenerated from the original (which targeted Julia 1.11.6).

### Results

All Phase 1 done criteria pass:

| Criterion | Result |
|-----------|--------|
| Forward simulation produces plausible trajectories | mRNA → 2.0, protein → 39.7 (matches analytic steady states) |
| Structural identifiability (sensitivity matrix rank) | Rank 4/4 — full rank |
| Parameter recovery (NUTS, 1000 samples) | All 5 parameters within 90% CI of true values |
| Posterior predictive checks | Trajectory envelope brackets true solution and data |

Sampling completes in ~2 minutes on a single CPU core. GPU acceleration is not beneficial at this problem scale.

### What this does not yet do

- Multi-module composition (Phase 2: metabolism coupled to TX/TL)
- Stochastic gene expression via SSA/Gillespie (Phase 2: simulation-based inference)
- Boundary protocol for uncertainty propagation across inference blocks
- Real experimental data (synthetic twin experiments only)
