# Parameters

Source: [`src/parameters.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/parameters.jl).

`InferParameter` is the unit of parameter declaration shared by every sub-model. The role-filtering helpers below let the inference layer separate kinetic rates, initial conditions, and observation parameters without case-by-case logic.

## `InferParameter`

```julia
InferParameter(value, prior, fixed, name, module_id, role)
```

A single parameter declared by a sub-model: nominal value, prior, identity, and role.

- `value::Float64` — nominal point value used for forward simulation and as the starting point for inference.
- `prior::Distribution` — Distributions.jl prior, sampled by NUTS / ABC-SMC.
- `fixed::Bool` — when `true`, the parameter is held at `value` and not inferred.
- `name::Symbol` — global identity. Parameters with the same name across sub-models are treated as **shared**; the orchestrator validates that their `value`, `prior`, and `fixed` agree.
- `module_id::Symbol` — owning sub-model (`:txl`, `:metab`, `:bge`, …).
- `role::Symbol` — `:rate`, `:initial_condition`, or `:observation`.

## Role-filtering helpers

```julia
free_params(params)         # !fixed
rate_params(params)         # role == :rate
ic_params(params)           # role == :initial_condition
obs_params(params)          # role == :observation
model_free_params(params)   # free and role != :observation  (dynamics RHS params)
obs_free_params(params)     # free and role == :observation  (likelihood params)
ode_free_params             # backward-compat alias for model_free_params
unique_params(params)       # dedupe by name, first-occurrence wins
```

Each takes a `Vector{InferParameter}` and returns the filtered subset; chain them as needed.
