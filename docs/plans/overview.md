# InferCell: Project Overview

> Composable, inference-first whole-cell modelling in Julia.

## What InferCell is

Whole-cell models simulate a living cell from its molecular parts. The field has made real progress wiring sub-models of transcription, translation, metabolism, and replication into integrated simulations. But the software architectures that run these simulations were not designed for inference. Existing platforms treat sub-models as black boxes, orchestrated via message passing or multi-language pipelines. You can't push gradients or likelihoods through boundaries you can't see inside.

InferCell asks a different question: *"Given noisy observations of cellular dynamics, what can we learn about the parameters governing those dynamics, and how certain are we?"*

The Luthey-Schulten Lab's [4D Whole-Cell Model](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM) (MC4D) is a massive forward simulation: 493 genes, 4 coupled solvers, spatial resolution at 10nm voxels, GPU-accelerated, published in Cell (2026). It answers: *"Can we simulate a complete minimal cell cycle in 4D?"* MC4D cannot do Bayesian inference. Its architecture (Python glue orchestrating 4 separate solvers via a shared `sim_properties` dictionary) has no gradient path and no likelihood function. The two projects are complementary, not competing.

Instead of orchestrating isolated black-box simulators, InferCell compiles hybrid dynamics (ODEs/SDEs + stochastic simulation + discrete events) into a single computational graph in pure Julia. This makes AD, likelihood evaluation, and Bayesian calibration possible without extra plumbing.

The core output is a whole-cell model that can do Bayesian inference natively. The specific "lightweight" biological content (which genes, which reactions) is chosen to exercise the architecture. We can then add higher order complexity in at a later date. The framework is organism-agnostic (not committed to JCVI-syn3A), though syn3A should be supported for comparison against MC4D.

## Why this matters

1. **Uncertainty quantification.** No existing whole-cell model quantifies posterior uncertainty over parameters. MC4D uses ~1000 point estimates from literature; InferCell produces full posterior distributions.
2. **Model selection.** Compare competing biological hypotheses (e.g., alternative gene regulation mechanisms) using Bayesian model comparison, not just visual trajectory matching.
3. **Principled calibration.** Learn parameters from data without hand-tuning. As experimental data becomes available, the framework calibrates automatically.

## Architecture

### The inference graph

An InferCell model is a **graph of biological modules**, each with a forward model and an inference backend, connected by typed interfaces that carry uncertainty.

```
Parameters --> SubModels --> Orchestrator --> InferenceGraph
                                                  |
                                         ---------+---------
                                         |        |        |
                                    DiffBlock  SimBlock  BoundaryProtocol
                                    (NUTS)     (SBI)    (uncertainty flow)
```

- **DifferentiableBlock:** Groups ODE/SDE modules. Builds a joint Turing.jl model. Infers with NUTS/HMC or some other fast likelihood based technique.
- **SimulationBlock:** Groups SSA/non-differentiable modules. Runs forward simulations. Infers with SBI (e.g., ABC-SMC, neural posterior estimation).
- **BoundaryProtocol:** Handles uncertainty propagation between blocks (see [Boundary Protocol](#boundary-protocol)).

Each graph node declares:
- **Forward dynamics:** ODE, SDE, SSA, or other formalism
- **Inference mode:** `:differentiable`, `:simulation`, or `:auto`
- **Parameters:** with priors, as `InferParameter` objects
- **Inputs/outputs:** state variables read from or written to other modules

Graph edges represent:
- **Shared parameters:** the same rate constant appearing in multiple modules
- **State coupling:** output of one module feeds into another

The framework partitions the model into **inference blocks** (i.e. groups of modules that share an inference backend) and handles uncertainty propagation across block boundaries.

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

`InferParameter` is deliberately non-parametric. It is metadata read once during model construction, never in a hot loop. The hot path uses a flat `Vector{Float64}` parameter vector built by the orchestrator.

**No explicit transform field.** When the `@model` function samples `k_tx ~ LogNormal(0, 1)`, Turing internally uses Bijectors.jl to sample in unconstrained space and transform back. The prior encodes the constraint; Turing handles the rest.

**Initial conditions** are `InferParameter` with `role=:initial_condition` and `fixed=true`. They can be made inferrable by flipping the `fixed` flag. The inference layer treats them identically to rate parameters.

Adding a new sub-model automatically extends the inference problem. Sensitivity analysis is a flag toggle, not a code change.

### SubModels

Each sub-model is a Julia struct implementing a common interface:

```julia
abstract type AbstractSubModel end

states(m::AbstractSubModel)              # --> vector of Symbol
parameters(m::AbstractSubModel)          # --> vector of InferParameter
dynamics(u, p, t, m::AbstractSubModel)   # out-of-place rate equations (returns du)
inputs(m::AbstractSubModel)              # --> vector of Symbol (states from other modules)
formalism(m::AbstractSubModel)           # --> :ode, :sde, :jump
inference_mode(m::AbstractSubModel)      # --> :differentiable, :simulation, :auto
```

`parameters(m)` returns *all* `InferParameter` objects: rate parameters, initial conditions (`role=:initial_condition`), and observation parameters (`role=:observation`). The orchestrator and inference layer filter by `role` and `fixed` as needed.

`dynamics` returns `du` rather than mutating it in-place. This is required for reliable AD (ReverseDiff.jl and Zygote.jl do not support mutation). Negligible performance cost for small state vectors. Sub-models define the 4-argument `dynamics(u, p, t, m)` for dispatch; the orchestrator wraps this into the 3-argument `f(u, p, t)` closure that DifferentialEquations.jl expects.

When `inference_mode` is `:auto`, the orchestrator selects based on `formalism(m)`:
- `:ode`, `:sde` --> `:differentiable`
- `:jump`, `:ssa` --> `:simulation`

Sub-models declare `inputs()`, the state variables they read but don't own. The orchestrator resolves these at composition time via a `SubModelContext`:

```julia
struct SubModelContext
    state_idxs::UnitRange{Int}    # indices into the flat state vector
    param_idxs::UnitRange{Int}    # indices into the flat parameter vector
    input_map::Dict{Symbol, Int}  # input name --> index in flat state vector
end
```

Coupling errors (undeclared inputs, missing providers) are caught at composition time, not runtime.

### Orchestrator

The orchestrator composes sub-models into a single `DEProblem`. Its job is mechanical --- no biology lives here.

1. **Collects states.** Builds a flat state vector with index mapping.
2. **Resolves coupling.** Verifies that every declared input is owned by another sub-model.
3. **Determines formalism.** Promotes to the most general type needed:
   - All `:ode` --> `ODEProblem`
   - Any `:sde` --> `SDEProblem` (continuous sub-models contribute zero noise)
   - Any `:jump` --> `JumpProblem` wrapping the continuous base
4. **Partitions into inference blocks.** Groups modules by inference mode into differentiable and simulation blocks.
5. **Builds combined RHS.** Dispatches to each sub-model's dynamics using `@views` to slice flat vectors.
6. **Collects parameters.** Flat parameter vector for the inference layer.
7. **Collects events.** Bundles callbacks into a `CallbackSet`.

```julia
function build_problem(models::Vector{<:AbstractSubModel}; tspan=(0.0, 100.0))
    return problem  # ODEProblem, SDEProblem, or JumpProblem
end
```

### Inference layer

The inference layer builds a Bayesian inference problem from the model's parameter declarations.

**Differentiable block (NUTS/HMC):**

1. Collects free parameters from all `InferParameter` objects where `fixed=false`.
2. A single generic Turing `@model` function iterates over parameter declarations at runtime: samples each free parameter from its prior, solves the forward model, evaluates the likelihood. One reusable definition that reads metadata.
3. Runs NUTS via Turing.jl. Posterior predictive checks sample from the posterior and run forward.

Note: we may need to explore if NUTS/HMC is the best way to do this, or if there are other techniques to explore.

**Simulation block (SBI):**

For non-differentiable modules (SSA, jump processes), inference uses simulation-based methods. v1 targets ABC-SMC (simple to implement, enough for v1). More advanced backends (neural posterior estimation) can be swapped in later. StochasticAD.jl may also be useful here. 

**Observation data:**

```julia
struct ObservedData
    times::Vector{Float64}          # observation timepoints
    observations::Matrix{Float64}   # species x timepoints
    species::Vector{Symbol}         # which species are observed
end
```

**API:**

```julia
function infer(model, data::ObservedData; sampler=NUTS(), n_samples=1000)
    turing_model = build_turing_model(model, data)
    chain = sample(turing_model, sampler, n_samples)
    return chain
end
```

`infer()` accepts a single sub-model or a composed model. Fit modules independently, then use sub-model posteriors as informed priors for joint calibration.

### Boundary protocol

This is the part that will require some thought and careful consideration propagating uncertainty across inference blocks. This is an instance of **modular Bayesian inference** --- doing inference in separate modules and stitching the results together. The literature on cut posteriors ([Plummer 2015](https://doi.org/10.1007/s11222-014-9503-z), [Jacob et al. 2017](https://arxiv.org/abs/1708.08719)) and expectation propagation ([Minka 2001](https://dl.acm.org/doi/10.5555/2074022.2074067)) covers the conditions under which modular posteriors are coherent. InferCell implements three levels of sophistication:

**Level 1 --- Sequential conditioning (v1 target):**
Infer differentiable block (NUTS) --> posterior samples --> condition simulation block on those samples --> SBI per sample. Information flows one direction only. This produces a **cut posterior**, not a joint posterior --- the simulation block does not feed back into the differentiable block. Valid when the blocks are conditionally independent given shared boundary variables; biased otherwise.

**Level 2 --- Gibbs-like alternation (v1 demonstration):**
Alternate: fix stochastic module state --> NUTS on differentiable params; fix differentiable params --> SBI on stochastic params. Information flows both ways. Because the ABC-SMC step produces approximate (not exact) conditional samples, this is a pseudo-marginal Gibbs sampler ([Andrieu and Roberts 2009](https://doi.org/10.1214/07-AOS574)). The stationary distribution is biased by the ABC tolerance --- it converges to an approximation of the joint posterior, not the exact joint. The quality of this approximation must be characterised empirically. See [Open Design Questions](#open-design-questions).

**Level 3 --- Particle MCMC / pseudo-marginal (future work):**
Stochastic modules contribute unbiased likelihood estimates via particle filters into a joint MCMC sampler. Asymptotically exact. Good Julia infrastructure exists (AdvancedMH.jl, SequentialMonteCarlo.jl).

The boundary protocol representation (full posterior samples, summary statistics, normalizing flow approximations) is an empirical question. v1 implements Level 1, demonstrates Level 2, defers Level 3.

## v1 scope

### Biological modules

Three modules that exercise different aspects of the inference graph:

| Module | Dynamics | Inference Block | Parameters | Purpose |
|--------|----------|-----------------|------------|---------|
| Transcription-Translation | ODE | Differentiable (NUTS) | k_tx, k_tl, gamma_mRNA, gamma_protein, sigma_obs | Baseline. Spike-validated. |
| Light Metabolism | ODE | Differentiable (NUTS, joint with TX/TL) | ~5-10 metabolic rate constants | Multi-module composition within a differentiable block. Shared parameters with TX/TL via energy/nucleotide costs. |
| Stochastic Gene Expression | SSA (Gillespie) | Simulation (SBI) | TX/TL rates, burst parameters | Same biology as TX/TL but stochastic. Exercises the inference graph boundary. |

**Light Metabolism:** A minimal metabolic core (5-10 reactions) producing ATP and nucleotides. Simplified glycolysis or energy source --> ATP production, nucleotide synthesis (feeds transcription), amino acid pool (feeds translation), ATP consumption by TX/TL. Not genome-scale; the point is demonstrating multi-module composition and parameter sharing.

**Stochastic Gene Expression:** The same TX/TL biology modeled as an SSA process. This creates a direct comparison: same biology modeled two ways, inferred two ways. If the posteriors agree, the framework works. Also shows that SBI can recover parameters from stochastic trajectories where NUTS can't be applied.

### What we take from MC4D

Reuse where appropriate:
- **Biological coupling graph:** what talks to what (metabolism feeds nucleotides to TX, TX produces mRNAs consumed by TL)
- **Kinetic parameters and rate laws:** from `kinetic_params.xlsx` and `GIP_rates.py` (curated from literature)
- **Organism data:** genome (`syn3A.gb`), metabolic model SBML (`Syn3A_updated.xml`) as optional inputs
- **Validation targets:** growth rate, mRNA/protein counts, metabolic flux distributions

Do NOT inherit:
- Hook interrupt architecture, `sim_properties` dictionary pattern, separate solver per formalism, Python/Cython/CUDA stack

### Key figures

1. Architecture diagram: the inference graph with differentiable and simulation blocks, boundary protocol.
2. Parameter recovery: posterior distributions for all ~15-30 parameters, ground truth overlaid.
3. ODE vs SSA comparison: same biology modeled two ways, consistent posteriors across formalisms and inference backends.
4. Uncertainty propagation: posterior uncertainty flowing through the boundary protocol between inference blocks.
5. Model selection: Bayes factors comparing competing hypotheses (e.g., constitutive vs bursty transcription).

### Deferred scope

Explicitly out of scope for v1:

- Spatial dynamics (RDME)
- Chromosome dynamics (Brownian Dynamics)
- Cell division and growth
- DNA replication
- Genome-scale metabolism
- Real experimental data (synthetic twin experiments first)
- Population-level inference (single-cell for v1)
- Moment closure (available as optional tool, not the core approach)

## Open design questions

These need answers before or during implementation:

1. **Shared parameters vs `SubModelContext.param_idxs`.** `param_idxs` is a `UnitRange{Int}`, assuming each sub-model's parameters occupy a contiguous, non-overlapping slice. But v1 requires shared parameters across modules (e.g., the same rate constant in TX/TL and metabolism). Resolution options: aliasing at orchestrator level (one slot in flat vector, multiple modules get the same index), or change to `Vector{Int}`.

2. **`ObservedData` is too rigid for SBI.** The struct assumes `Matrix{Float64}` at fixed timepoints. ABC-SMC typically operates on summary statistics with a distance metric, not raw time-series observations. The simulation block needs a different or more general observation/comparison interface.

3. **No summary statistics abstraction.** The current `log_likelihood` callable assumes evaluable likelihoods. SBI cannot provide this. The simulation block needs its own interface: summary statistics extraction, distance metric, acceptance threshold.

4. **`formalism()` returns bare Symbols.** `:ode`, `:sde`, `:jump`, `:ssa` are unchecked. A typo like `:ODE` fails silently. Consider an enum or type hierarchy for type safety.

5. **Level 2 boundary protocol convergence.** Gibbs-like alternation between NUTS and ABC-SMC is heuristic --- the ABC-SMC step does not produce exact conditional samples, so standard Gibbs convergence guarantees do not apply. This should be presented as experimental, with Level 3 (particle MCMC) as the theoretically complete solution.

6. **Parameter identifiability.** With 15-30 parameters and potentially correlated modules, structural and practical identifiability is a concern. Synthetic twin experiments provide empirical checks, but identifiability should be verified analytically (e.g., sensitivity matrix rank) before running any sampler.

7. **SBI risk is higher than rated.** The stochastic gene expression module is what makes the inference graph real --- it's the whole point of the mixed-backend architecture. If Julia SBI tooling proves insufficient, the central demonstration is compromised. The PythonCall fallback introduces the Python dependency that the single-language design principle exists to avoid. Mitigations: keep the stochastic module's parameter count low (~5-8), start ABC-SMC implementation early to surface issues.

8. **ForwardDiffSensitivity scaling.** Phase 1 uses `ForwardDiffSensitivity()` which costs O(n_params) per gradient evaluation. At 5 parameters this is fine; at 15-30 parameters across coupled ODE modules, gradient cost scales linearly and NUTS wall-clock may become impractical. Benchmark at increasing parameter counts before committing to v1 scope. If scaling is poor, switch to adjoint sensitivity methods (`InterpolatingAdjoint` or `BacksolveAdjoint` from SciMLSensitivity.jl).

9. **Model misspecification.** All v1 validation uses synthetic twin experiments where the data-generating model matches the inference model. This tests self-consistency, not robustness. Real biological models are always wrong to some degree, and Bayesian inference on a misspecified model can produce posteriors that concentrate on incorrect parameter values ([Kleijn and van der Vaart 2012](https://doi.org/10.1214/12-EJS675)). v1 should include at least one experiment where the data-generating process differs from the inference model (e.g., data from Michaelis-Menten kinetics, inference with mass-action kinetics) to characterise how the framework behaves under misspecification.

10. **Observation model specification.** The `ObservedData` struct defines the format but not the content. Identifiability, posterior geometry, and computational cost all depend on what is observed, at what temporal resolution, and with what noise structure. Observing mRNA + protein at high resolution is a fundamentally different inference problem from observing protein-only at low resolution. The observation model for each module should be specified explicitly before running inference.

## Design principles

- **Information flow as design test.** For every module boundary: "how does inference work across this interface?" If you can't answer, the design needs to change.
- **Interfaces over implementations.** Sub-models are defined by contracts, not internals. Swappable without touching the rest of the system.
- **Single-language computational graph.** Julia removes language boundaries. Solvers, models, and inference share one compiled representation.
- **Synthetic first, real data ready.** Validate on twin experiments where ground truth is known. Interface designed so real data slots in as a new likelihood function.

## Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| SBI tooling in Julia is immature | **High** | Start ABC-SMC early. Keep stochastic module params low (~5-8). Wrap Python `sbi` via PythonCall only as last resort. |
| Boundary protocol design is unclear | High | Start with sequential conditioning (Level 1). Validate empirically. Iterate. |
| Joint inference across blocks may not converge | High | Synthetic twin experiments provide known ground truth for diagnosing failures. |
| Light metabolism scope is vague | Low | Define 5-10 reactions producing ATP and nucleotides. Couple to TX/TL. Keep minimal. |
| "Why not SBI on everything?" | Medium | Demonstrate empirically: NUTS on differentiable modules is orders of magnitude more sample-efficient. Hybrid is strictly better. |
| Parameter dimensionality too high for SBI | Medium | Keep stochastic module parameters low (~5-8). ABC-SMC works well in this regime. |
| ForwardDiffSensitivity scaling | Medium | Benchmark gradient cost at 5, 10, 20, 30 params. Switch to adjoint methods if wall-clock becomes impractical. |

## Dependencies

| Package | Role |
|---------|------|
| OrdinaryDiffEq.jl | ODE solvers (Tsit5) |
| SciMLBase.jl | Problem/solution types |
| SciMLSensitivity.jl | Adjoint methods for gradients through solvers |
| Turing.jl | Bayesian inference, NUTS sampler (uses Bijectors.jl internally) |
| Distributions.jl | Priors, likelihoods |
| StaticArrays.jl | Fast small-array allocations for out-of-place dynamics |
| CairoMakie.jl | Plotting |
| Catalyst.jl | *Planned.* Optional, for reaction-network sub-models |

## Project structure

```
InferCell/
├── src/
│   ├── InferCell.jl           # main module, exports public API
│   ├── parameters.jl          # InferParameter struct and utilities
│   ├── interface.jl           # AbstractSubModel, required method definitions
│   ├── orchestrator.jl        # build_problem, state/param flattening, formalism promotion
│   ├── inference.jl           # build_turing_model, infer(), observe(), posterior_predictive()
│   ├── likelihoods.jl         # observation models (Gaussian, Poisson, custom)
│   └── models/                # sub-model implementations
├── test/
│   └── spike-test/            # pre-architecture spike validation (completed)
├── examples/
├── docs/
│   └── plans/
├── Project.toml
└── README.md
```

## Related documents

- [Phase 1 Design](2026-03-13-phase1-design.md) --- detailed TX/TL implementation spec, solver choices, done criteria, end-to-end workflow
- [Spike Plan + Results](../../test/spike-test/) --- AD-through-ODE validation (completed, all success criteria passed)
