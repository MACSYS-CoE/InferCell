```@meta
CurrentModule = InferCell
```

# ODE inference (NUTS)

!!! note "Stub page"
    This tutorial is a placeholder for v0.0.1. The narrative will fill out as the framework stabilises; for now, the runnable example below mirrors the [Getting started](../getting-started.md) ODE path.

InferCell dispatches sub-models with `inference_mode = :differentiable` to the No-U-Turn Sampler (NUTS) via [Turing.jl](https://turing.ml). Gradients flow through DifferentialEquations.jl via `ForwardDiffSensitivity`.

## Minimal example

```julia
using InferCell

# 1. Build a sub-model (priors and defaults live in the constructor)
model = TranscriptionTranslation()

# 2. Compose into an ODEProblem
prob = build_problem(model; tspan=(0.0, 50.0))

# 3. Generate synthetic data
sol  = solve(prob, Tsit5(); saveat=0:2.5:50)
data = observe(sol, 0:2.5:50, model; sigma=0.3)

# 4. Run NUTS
chain = infer(model, data; n_samples=1000)
```

## What to look at next

- [`build_problem`](@ref) — composing multiple sub-models.
- [`build_turing_model`](@ref) — when you need to override priors directly.
- [`check_identifiability`](@ref) — Jacobian-rank diagnostic before sampling.
