# Sub-model interface

Source: [`src/interface.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/interface.jl).

Every block in InferCell — TX/TL, metabolism, stochastic gene expression — is a subtype of `AbstractSubModel` that implements a small protocol. The orchestrator only ever talks to this interface; sub-model internals stay private.

## `AbstractSubModel`

Abstract supertype. Concrete sub-types implement the protocol functions below.

## Required protocol functions

```julia
states(m)        # Vector{Symbol} — ordered state-variable names
parameters(m)    # Vector{InferParameter} — every parameter, all roles
dynamics(u, p, t, m[, u_inputs])   # out-of-place ODE RHS (required if formalism == :ode)
reactions(m)     # Vector{ConstantRateJump} (required if formalism == :jump)
```

## Optional protocol functions (with defaults)

```julia
inputs(m)           # Symbol[]              — coupling inputs from other sub-models
formalism(m)        # :ode | :sde | :jump   — defaults to :ode
inference_mode(m)   # :differentiable | :simulation — defaults to :differentiable
```

`formalism` determines which problem type the orchestrator builds; `inference_mode` determines which sampler `infer()` dispatches to.

## `SubModelContext`

```julia
SubModelContext(state_idxs, param_idxs, input_map)
```

Per-sub-model bookkeeping built by the orchestrator: where the sub-model's states live in the global state vector, which global parameter indices it reads, and how its declared input symbols map onto global state indices. End users rarely touch this directly; it's the wiring the orchestrator threads through `dynamics()`.
