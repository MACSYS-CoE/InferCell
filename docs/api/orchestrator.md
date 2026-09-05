# Orchestrator

Source: [`src/orchestrator.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/orchestrator.jl).

The orchestrator composes one or more sub-models into a single problem object. Shared parameters are deduplicated by name; cross-block coupling is resolved through each sub-model's declared inputs.

## `build_problem`

```julia
build_problem(models; tspan=(0.0, 100.0), kwargs...)
build_problem(model;  tspan=(0.0, 100.0), kwargs...)
```

Returns:

- an `ODEProblem` when every sub-model has `formalism = :ode`
- a `JumpProblem` when every sub-model has `formalism = :jump`
- a [`HandshakeDriver`](handshake.md) for a mixed ODE/jump composition
- an error for `:sde`, which has no build path in any composition

Keyword arguments other than `tspan` are the handshake driver's declared policy (`interval`, `rounding`, `radius_nm`, `ode_solver`, `abstol`, `reltol`) and are refused on a homogeneous composition rather than silently ignored.

### Shared-parameter handling

Parameters with the same `name` across sub-models are treated as shared and stored once in the global parameter vector. Definitions must agree: same `value`, same `prior` type, same `fixed` flag. The orchestrator raises an error on inconsistency, naming the two offending modules.

### Coupling

Each sub-model declares its inputs via [`inputs`](interface.md). The orchestrator resolves these against the global state ownership table — every declared input must be a `state` of some sub-model in the composition — and threads the resolved values into `dynamics(u, p, t, m, u_inputs)`. Models with no inputs use the simpler `dynamics(u, p, t, m)` overload.

### Example

```julia
using InferCell

# Composed ODE: TX/TL feeds enzyme to metabolism; metabolism shares k_tx, k_tl with TX/TL
prob = build_problem([TranscriptionTranslation(), LightMetabolism()];
                     tspan=(0.0, 100.0))
sol  = solve(prob, Tsit5())
```
