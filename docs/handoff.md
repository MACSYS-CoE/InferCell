# Handoff

## Current state

Branch `step1-metabolism-shared-params` implements Step 1 of the [v1 minimal demo plan](plans/2026-04-13-v1-minimal-demo.md).

### What's working

- **TX/TL module** — ODE model with NUTS inference, posterior predictive, identifiability
- **Stochastic gene expression** — SSA model with ABC-SMC, automatic dispatch via `infer()`
- **LightMetabolism module** — minimal ODE metabolism (ATP, NTP, AA pools) coupled one-way to TX/TL
- **Shared parameter composition** — `build_problem([TranscriptionTranslation(), LightMetabolism()])` produces a single ODEProblem with 5 states and 7 deduplicated free parameters
- **121 tests pass** including unit and composition tests for the new module

### Key infrastructure added in Step 1

- `param_idxs` changed from `UnitRange{Int}` to `Vector{Int}` — enables shared parameter slots across modules
- Name-based parameter deduplication in orchestrator (`_build_contexts`, `_build_p0`) and inference (`build_turing_model`, `_infer_abc`)
- 5-arg `dynamics(u, p, t, m, u_inputs)` with fallback — enables cross-module state coupling (metabolism reads `mRNA` from TX/TL)
- `_validate_shared_params` — errors on inconsistent shared parameter definitions across modules

### Files changed/added

- `src/interface.jl` — `SubModelContext.param_idxs` type change, dynamics fallback
- `src/orchestrator.jl` — shared param dedup in `_build_contexts`/`_build_p0`, input passing in `_build_rhs`, validation
- `src/parameters.jl` — `unique_params` helper
- `src/inference.jl` — dedup in `build_turing_model` and `_infer_abc`
- `src/models/light_metabolism.jl` — new model
- `test/test_light_metabolism.jl`, `test/test_composition.jl` — new tests

## Next step

**Step 2 — Joint ODE inference.** Run NUTS on the composed TX/TL + metabolism system. The existing `_infer_nuts` path should work with minimal changes since the Turing model already reads `InferParameter` metadata generically. Main tasks:

- Synthetic twin experiment: generate data from known parameters, recover posteriors
- Validate ForwardDiffSensitivity scales at 7 ODE + 1 obs params
- Wall-clock benchmark on one HPC node

See `docs/plans/2026-04-13-v1-minimal-demo.md` for the full 4-step plan.
