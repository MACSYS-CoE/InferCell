# InferCell Phase 1 Design

> Transcription-translation model with end-to-end Bayesian parameter inference.

For the overall architecture, see [Architecture](2026-03-13-architecture.md).

## Goal

Prove the architecture works end-to-end on the simplest possible biological system. The biology here is deliberately boring; the point is validating that gradients flow from data through inference through the solver through the model and back.

## Biology: Transcription-Translation

A single sub-model with two state variables and four rate parameters:

```
d(mRNA)/dt    = k_tx - γ_mRNA * mRNA
d(protein)/dt = k_tl * mRNA - γ_protein * protein
```

| Parameter | Symbol | Prior | Description |
|-----------|--------|-------|-------------|
| Transcription rate | `k_tx` | `LogNormal(0, 1)` | mRNA production rate |
| Translation rate | `k_tl` | `LogNormal(0, 1)` | protein production per mRNA |
| mRNA degradation | `γ_mRNA` | `LogNormal(0, 1)` | mRNA decay rate |
| Protein degradation | `γ_protein` | `LogNormal(0, 1)` | protein decay rate |
| Observation noise | `σ_obs` | `HalfNormal(1.0)` | Inferred, not fixed |

### Implementation

```julia
struct TranscriptionTranslation <: AbstractSubModel
    params::Vector{InferParameter}
end

formalism(::TranscriptionTranslation) = :ode
states(::TranscriptionTranslation) = [:mRNA, :protein]
inputs(::TranscriptionTranslation) = Symbol[]  # no coupling in phase 1

# Out-of-place dynamics for AD compatibility
function dynamics(u, p, t, ::TranscriptionTranslation)
    mRNA, protein = u
    k_tx, k_tl, γ_mRNA, γ_protein = p
    SA[k_tx - γ_mRNA * mRNA, k_tl * mRNA - γ_protein * protein]
end
```

`parameters(m)` returns all seven `InferParameter` objects:
- Four rate parameters (`role=:rate`) with `LogNormal(0, 1)` priors
- Two initial conditions (`role=:initial_condition`, `fixed=true`): `mRNA₀=0.0`, `protein₀=0.0`
- One observation noise (`role=:observation`): `σ_obs` with `HalfNormal(1.0)` prior

All rate parameters are positive. `LogNormal` priors encode this constraint, and Turing handles unconstrained sampling internally via Bijectors.jl. No explicit transform code needed.

Initial conditions can be made inferrable by flipping `fixed=false`, but are fixed for phase 1.

Phase 1 uses a single `σ_obs` for both mRNA and protein, assuming equal observation noise across species. Per-species noise is a straightforward extension for later phases.

## Formalism: Pure ODE

Phase 1 uses deterministic ODEs only, which gives:
- Differentiability via AD
- Adjoint sensitivity methods via SciMLSensitivity.jl
- Compatibility with NUTS/HMC

The sub-model declares `formalism(m) = :ode`, so the orchestrator returns a plain `ODEProblem{false}` (out-of-place). No noise terms, no jumps.

### Solver and AD choices

| Choice | Value | Rationale |
|--------|-------|-----------|
| Solver | `Tsit5()` | Reliable explicit RK method, well-tested with AD |
| Sensitivity | `ForwardDiffSensitivity()` | Works reliably for small parameter counts (5 params) |
| Dynamics | Out-of-place | Required for broad AD backend compatibility |

These choices are deliberately conservative. Phase 1 validates that the pipeline works; solver/AD optimization is a later concern.

## Data: Synthetic Twin Experiment

1. Choose "true" parameter values for all four rate parameters and `σ_obs`.
2. Solve the forward ODE to generate a ground-truth trajectory.
3. Sample observations at discrete timepoints with additive Gaussian noise (`σ_obs`).
4. Package as `ObservedData` for the inference layer.

### Data format

```julia
struct ObservedData
    times::Vector{Float64}          # observation timepoints
    observations::Matrix{Float64}   # species × timepoints
    species::Vector{Symbol}         # which species are observed
end
```

The `observe()` function (in `inference.jl`) generates `ObservedData` from a solution:

```julia
data = observe(true_sol, timepoints; σ=0.1)
```

### Observation model

The likelihood evaluates per-timepoint:

```julia
for i in 1:length(data.times)
    data.observations[:, i] ~ MvNormal(sol[:, i], σ_obs)
end
```

Where `σ_obs` is inferred (not fixed) with a `HalfNormal(1.0)` prior.

Both mRNA and protein are observed (full observability). Partial observation (e.g. protein only) is a straightforward extension but not needed for phase 1.

## Inference: NUTS via Turing.jl

A single generic Turing `@model` function iterates over the parameter declarations at runtime:

1. Sample each free `InferParameter` from its prior (including `σ_obs`). Turing handles unconstrained sampling internally; the prior encodes the constraint.
2. Collect sampled values into the flat parameter vector for the ODE solver.
3. Solve the ODE with `Tsit5()` and `ForwardDiffSensitivity()`.
4. Evaluate per-timepoint Gaussian likelihood against the observed data.

One reusable `@model` definition that reads `InferParameter` metadata, not macro-generated. Adding parameters to a sub-model automatically extends the inference problem.

NUTS is the sampler. It's the most demanding test that AD actually works through the full pipeline.

### Modular calibration

For phase 1, the model is a single sub-model, so modular vs. joint calibration is the same thing. But the `infer()` function is designed to accept either:

```julia
# These are equivalent in phase 1
chain = infer(txl_model, data)
chain = infer(composed_model, data)
```

This distinction becomes meaningful in phase 2 when multiple sub-models can be calibrated independently or jointly.

## Done Criterion

Phase 1 is complete when:

1. **Forward simulation** produces plausible mRNA and protein trajectories.
2. **Parameter recovery.** NUTS posteriors concentrate around the true parameter values used to generate synthetic data.
3. **Posterior predictive checks.** Trajectories sampled from the posterior reproduce the synthetic data distribution. If posterior predictions match the data, the pipeline works.

## End-to-End Workflow

```julia
using InferCell

# 1. Define sub-model
txl = TranscriptionTranslation()

# 2. Build problem
prob = build_problem([txl]; tspan=(0.0, 100.0))

# 3. Generate synthetic data
true_sol = solve(prob, Tsit5())
data = observe(true_sol, 0.0:5.0:100.0; σ=0.1)

# 4. Infer (σ_obs is inferred alongside rate parameters)
chain = infer(txl, data)

# 5. Posterior predictive check
ppc = posterior_predictive(txl, chain)
```

## Project Structure (Phase 1 Files)

```
src/
├── InferCell.jl                          # main module
├── parameters.jl                          # InferParameter struct
├── interface.jl                           # AbstractSubModel definition
├── orchestrator.jl                        # build_problem
├── inference.jl                           # build_turing_model, infer(), posterior_predictive()
├── likelihoods.jl                         # Gaussian observation model
└── models/
    └── transcription_translation.jl       # TX/TL sub-model

test/
├── runtests.jl
├── test_parameters.jl
├── test_orchestrator.jl
├── test_inference.jl
└── test_txl.jl                            # forward sim + parameter recovery

examples/
└── phase1_demo.jl                         # end-to-end demo script
```

## Dependencies

| Package | Role in Phase 1 |
|---------|-----------------|
| DifferentialEquations.jl | ODE solver (`Tsit5()`) |
| Turing.jl | NUTS sampler (uses Bijectors.jl internally for unconstrained sampling) |
| Distributions.jl | Priors (LogNormal, HalfNormal), likelihood (MvNormal) |
| SciMLSensitivity.jl | `ForwardDiffSensitivity()` for gradients through ODE solve |
| Catalyst.jl | Not required for phase 1, but available |
