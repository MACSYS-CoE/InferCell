# Parameters

Source: [`src/parameters.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/parameters.jl).

`InferParameter` is the unit of parameter declaration shared by every sub-model. The role-filtering helpers below let the inference layer separate kinetic rates, initial conditions, and observation parameters without case-by-case logic.

## `InferParameter`

```julia
InferParameter(value, prior, fixed, name, module_id, role)
InferParameter(value, prior, fixed, name, module_id, role, provenance)
```

A single parameter declared by a sub-model: nominal value, prior, identity, and role.

- `value::Float64` — nominal point value used for forward simulation and as the starting point for inference.
- `prior::Distribution` — Distributions.jl prior, sampled by NUTS / ABC-SMC.
- `fixed::Bool` — when `true`, the parameter is held at `value` and not inferred.
- `name::Symbol` — global identity. Parameters with the same name across sub-models are treated as **shared**; the orchestrator validates that their `value`, `prior`, and `fixed` agree.
- `module_id::Symbol` — owning sub-model (`:txl`, `:metab`, `:bge`, …).
- `role::Symbol` — `:rate`, `:initial_condition`, or `:observation`.
- `provenance::Union{ParameterSource, Nothing}` — optional record of the file the value was imported from. Defaults to `nothing`, so the six-argument form above remains valid.

## `ParameterSource`

```julia
ParameterSource(file; table, identifier, informedness, alternatives)
```

Where an imported value came from. Carried per value rather than resolved once by hand, because the same cross-file ambiguity has already produced two errors of record in the source model — see [Core A′ interface contract](corea-interface.md) for the loader that produces these.

- `file::String` — the source file, as provenance records it.
- `table`, `identifier` — the table within the file and the identifier the value appeared under.
- `informedness::Symbol` — `:balanced`, `:prior_default`, `:asserted`, or `:not_imported`.
- `alternatives::Vector{Pair{String,Float64}}` — the same identifier's value in every *other* file it appears in, so a value chosen in the presence of an alternative says so.

## Provenance helpers

```julia
provenance_of(p)              # the ParameterSource, or nothing
source_file(p)                # the file name, or nothing
informedness(p)               # :balanced / :prior_default / :asserted / :not_imported
asserted_prior_params(params) # priors this project asserted, not inherited
uninformed_params(params)     # sitting at the prior median at prior width
provenance_conflicts(params)  # names arriving from more than one source file
```

## Role-filtering helpers

```julia
free_params(params)         # !fixed
rate_params(params)         # role == :rate
ic_params(params)           # role == :initial_condition
obs_params(params)          # role == :observation
model_free_params(params)   # free and role != :observation  (dynamics RHS params)
obs_free_params(params)     # free and role == :observation  (likelihood params)
ode_free_params             # backward-compat alias for model_free_params
unique_params(params)       # dedupe by name; a provenance-carrying copy wins over a bare one
```

Each takes a `Vector{InferParameter}` and returns the filtered subset; chain them as needed.
