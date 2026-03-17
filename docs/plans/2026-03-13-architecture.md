# InferCell Architecture

> Composable, inference-first whole-cell modelling in Julia.

## Vision

InferCell is a minimal whole-cell model designed so that inference works from the start. The model and the architecture are co-developed: biology validates the platform, the platform makes the biology calibratable.

Instead of orchestrating isolated black-box simulators, InferCell compiles hybrid dynamics (ODEs/SDEs + stochastic simulation + discrete events) into a single computational graph in pure Julia. AD, likelihood evaluation, and Bayesian calibration fall out of this design.

## Phased Roadmap

- **Phase 1 Spike:** Validate AD through ODE solve in a single script before building abstractions. See [Spike Plan](2026-03-13-spike-plan.md).
- **Phase 1:** Transcription-translation with ODEs, NUTS inference, synthetic data, parameter recovery + posterior predictive checks. See [Phase 1 Design](2026-03-13-phase1-design.md).
- **Phase 2:** Add coarse metabolism (energy/resource coupling), SDE formalism option, simulation-based inference fallbacks.
- **Phase 3:** Cell growth and division events, hybrid discrete-continuous dynamics, full likelihood-free inference for non-differentiable components.

## Four-Layer Architecture

```
Parameters → SubModels → Orchestrator → Inference
```

### Parameters

Each parameter carries its value, prior distribution, fixed/free flag, and metadata. The inference layer reads these declarations to build the Bayesian model automatically.

```julia
struct InferParameter
    value::Float64        # current/default value
    prior::Distribution   # e.g. LogNormal(0, 1)
    fixed::Bool           # if true, excluded from inference
    name::Symbol          # e.g. :k_transcription
    module_id::Symbol     # which sub-model owns this
    role::Symbol          # :rate, :initial_condition, or :observation
end
```

`InferParameter` is deliberately non-parametric. It is metadata read once during model construction, never in a hot loop. The abstract `prior::Distribution` field lets parameters with different prior types (`LogNormal`, `HalfNormal`, etc.) coexist in a plain `Vector{InferParameter}`. The hot path uses a flat `Vector{Float64}` parameter vector built by the orchestrator.

**No explicit transform field.** When the `@model` function samples `k_tx ~ LogNormal(0, 1)`, Turing internally uses Bijectors.jl to sample in unconstrained space and transform back. The sampled value the model sees is already constrained (positive). No user-facing transform is needed between Turing and the ODE solver. The prior encodes the constraint, Turing handles the rest.

**Initial conditions** are represented as `InferParameter` with `role=:initial_condition` and `fixed=true`. They can be made inferrable by flipping the `fixed` flag. The inference layer treats them identically to rate parameters.

This design means:
- Adding a new sub-model automatically extends the inference problem (its parameters carry their own priors).
- Sensitivity analysis is a flag toggle, not a code change.
- HMC/NUTS works in unconstrained space; Turing handles the mapping and jacobian correction internally.

### SubModels

Each sub-model is a Julia struct implementing a common interface:

```julia
abstract type AbstractSubModel end

states(m::AbstractSubModel)       # → vector of Symbol
parameters(m::AbstractSubModel)   # → vector of InferParameter (rate params, ICs, observation params)
dynamics(u, p, t, m::AbstractSubModel)  # out-of-place rate equations (returns du)
inputs(m::AbstractSubModel)       # → vector of Symbol (states read from other modules)
formalism(m::AbstractSubModel)    # → :ode, :sde, :jump
```

`parameters(m)` returns *all* `InferParameter` objects for a sub-model: rate parameters, initial conditions (`role=:initial_condition`), and observation parameters (`role=:observation`). The orchestrator and inference layer filter by `role` and `fixed` as needed.

**Out-of-place dynamics.** `dynamics` returns `du` rather than mutating it in-place. This is required for reliable AD. ForwardDiff.jl works with in-place functions but ReverseDiff.jl and Zygote.jl do not. The performance cost is negligible for the small state vectors in phase 1, and SciML supports `ODEProblem{false}` (out-of-place) natively.

**Dynamics wrapping.** Sub-models define 4-argument `dynamics(u, p, t, m)` for dispatch. The orchestrator wraps this into the 3-argument `f(u, p, t)` closure that DifferentialEquations.jl expects.

**Coupling contract:** Sub-models declare inputs, state variables they read but don't own. The orchestrator resolves these at composition time via a `SubModelContext`:

```julia
struct SubModelContext
    state_idxs::UnitRange{Int}    # indices into the flat state vector
    param_idxs::UnitRange{Int}    # indices into the flat parameter vector
    input_map::Dict{Symbol, Int}  # input name → index in flat state vector
end
```

At composition time, the orchestrator builds a `SubModelContext` for each sub-model. The composed dynamics function uses `@views` to slice the flat vectors, passing each sub-model only its own states and parameters plus read-only access to declared inputs. Coupling errors (undeclared inputs, missing providers) are caught at composition time, not runtime.

**Formalism freedom:** Each sub-model declares whether it's `:ode`, `:sde`, or `:jump`. Internally, it can use Catalyst.jl to define reactions, hand-written ODEs, a neural surrogate, or constraint-based methods. The interface doesn't care; only the contract matters.

### Orchestrator

The orchestrator composes sub-models into a single `DEProblem`. Its job is mechanical, no biology lives here. It does not commit to a single problem type; the composed formalism emerges from what the sub-models declare.

1. **Collects states.** Builds a flat state vector with index mapping.
2. **Resolves coupling.** Verifies that every declared input is owned by another sub-model. Errors at composition time.
3. **Determines formalism.** Promotes to the most general type needed:
   - All `:ode` → `ODEProblem`
   - Any `:sde` → `SDEProblem` (continuous sub-models contribute zero noise)
   - Any `:jump` → `JumpProblem` wrapping the continuous base
4. **Builds combined RHS.** Dispatches to each sub-model's dynamics. Also builds noise functions (SDE) or collects jump sets as needed.
5. **Collects parameters.** Flat parameter vector for the inference layer.
6. **Collects events.** Bundles callbacks (e.g. division) into a `CallbackSet`.
7. **Returns a `DEProblem`:**

```julia
function build_problem(models::Vector{<:AbstractSubModel}; tspan=(0.0, 100.0))
    return problem  # ODEProblem, SDEProblem, or JumpProblem
end
```

### Inference

The inference layer builds a Bayesian inference problem from the model's parameter declarations.

1. **Collects free parameters.** Reads priors from all `InferParameter` objects where `fixed=false`.
2. **Generic Turing `@model`.** A single `@model` function iterates over parameter declarations at runtime: samples each free parameter from its prior (Turing handles unconstrained sampling internally via Bijectors.jl), solves the forward model, and evaluates the likelihood. One reusable model definition that reads metadata, not macro-generated.
3. **Observation data format:**
   ```julia
   struct ObservedData
       times::Vector{Float64}          # observation timepoints
       observations::Matrix{Float64}   # species × timepoints
       species::Vector{Symbol}         # which species are observed
   end
   ```
   The `observe()` function (in `inference.jl`) generates `ObservedData` from a solution for synthetic experiments.
4. **Flexible likelihood.** Accepts any callable `log_likelihood(data, solution, params)`. Gaussian for phase 1; Poisson, partial observations, custom noise for later.
5. **Runs inference:**
   ```julia
   function infer(model, data::ObservedData; sampler=NUTS(), n_samples=1000)
       turing_model = build_turing_model(model, data)
       chain = sample(turing_model, sampler, n_samples)
       return chain
   end
   ```
6. **Posterior predictive checks.** Samples from the posterior, runs forward, returns trajectory ensembles.

**Modular calibration:** `infer()` accepts a single sub-model or a composed model. Fit modules independently, then use sub-model posteriors as informed priors for joint calibration.

**Inference strategy adapts to formalism:** Fully continuous models use NUTS. Models with jump processes fall back to simulation-based methods (ABC, synthetic likelihood) or hybrid schemes. The `infer()` signature stays the same.

## Design Principles

- **Information flow as design test.** For every module boundary: "how does inference work across this interface?" If you can't answer, the design needs to change.
- **Interfaces over implementations.** Sub-models are defined by contracts, not internals. Swappable without touching the rest of the system.
- **Single-language computational graph.** Julia removes language boundaries. Solvers, models, and inference share one compiled representation.
- **Synthetic first, real data ready.** Validate on twin experiments where ground truth is known. Interface designed so real data slots in as a new likelihood function.

## Core Dependencies

| Package | Role |
|---------|------|
| DifferentialEquations.jl | ODE/SDE/Jump solvers |
| Turing.jl | Bayesian inference, NUTS sampler (uses Bijectors.jl internally) |
| Distributions.jl | Priors, likelihoods |
| SciMLSensitivity.jl | Adjoint methods for gradients through solvers |
| Catalyst.jl | Optional, for reaction-network sub-models |

## Project Structure

```
InferCell/
├── spike/
│   └── spike_txl_inference.jl # pre-architecture validation (see spike plan)
├── src/
│   ├── InferCell.jl           # main module, exports public API
│   ├── parameters.jl          # InferParameter struct and utilities
│   ├── interface.jl           # AbstractSubModel, required method definitions
│   ├── orchestrator.jl        # build_problem, state/param flattening, formalism promotion
│   ├── inference.jl           # build_turing_model, infer(), observe(), posterior_predictive()
│   ├── likelihoods.jl         # observation models (Gaussian, Poisson, custom)
│   └── models/                # sub-model implementations
├── test/
├── examples/
├── docs/
│   └── plans/
├── Project.toml
└── readme.md
```
