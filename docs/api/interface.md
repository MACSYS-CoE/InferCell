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
coupling(m)         # CouplingEdge[]        — typed coupling declarations
contributed_states(m)                        # Symbol[] — registry states of *other* modules this one adds derivative terms to
contributions(u, p, t, m[, u_inputs])        # static vector, one signed term per contributed state
module_id(m)        # Symbol                — identity, defaults to the type name
reduction_notes(m)  # String[]              — simplifications this module introduces
formalism(m)        # :ode | :sde | :jump   — defaults to :ode
inference_mode(m)   # :differentiable | :simulation — defaults to :differentiable
```

`formalism` determines which problem type the orchestrator builds; `inference_mode` determines which sampler `infer()` dispatches to.

`inputs` says *which* state a sub-model reads from another; `coupling` says *how* it crosses, as one of seven typed edge kinds. Both default to empty, so a sub-model written before the Core A′ contract needs no change. See [Core A′ interface contract](corea-interface.md).

`contributed_states` and `contributions` are the outbound counterpart of `inputs` and `dynamics`: where `inputs` lets a module *read* a state another module integrates, `contributions` lets it *write* a derivative term into one. The orchestrator adds the term to the owner's derivative at the owner's global index; a module never writes outside its own slice directly. The species must be Core A′ registry names and may not be chemostatted (checked by the resolver), and must be owned by some module in the composition (checked by the orchestrator at build time). For an ODE module, every mass or currency edge (either direction) on a dynamic registry species it does not own must have a matching contribution, and every contribution a matching edge — the resolver checks both directions, as it does for `inputs`, and a module listing contributions with no edges is checked rather than skipped. A module that consumes a foreign pool contributes a negative term; one that produces into it contributes a positive one. A `:jump` module may not list contributions; its writes to a peer's state go through its reactions.

## `SubModelContext`

```julia
SubModelContext(state_idxs, param_idxs, input_map[, contrib_idxs])
```

Per-sub-model bookkeeping built by the orchestrator: where the sub-model's states live in the global state vector, which global parameter indices it reads, how its declared input symbols map onto global state indices, and the global index of each state in `contributed_states(m)`, in that order (`contrib_idxs`, empty by default). End users rarely touch this directly; it's the wiring the orchestrator threads through `dynamics()`.
