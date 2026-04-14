# Handoff

## Current state

Branch `main` implements Steps 1-2 of the [v1 minimal demo plan](plans/2026-04-13-v1-minimal-demo.md).

### What's working

- **TX/TL module** — ODE model with NUTS inference, posterior predictive, identifiability
- **Stochastic gene expression** — SSA model with ABC-SMC, automatic dispatch via `infer()`
- **LightMetabolism module** — minimal ODE metabolism (ATP, NTP, AA pools) coupled one-way to TX/TL
- **Shared parameter composition** — `build_problem([TranscriptionTranslation(), LightMetabolism()])` produces a single ODEProblem with 5 states and 7 deduplicated free parameters
- **Joint ODE inference** — `infer([txl, metab], data)` runs NUTS on the composed 5-state system, recovering all 8 parameters (7 ODE + 1 obs noise)
- **Multi-model utilities** — `observe`, `posterior_predictive`, `check_identifiability` all accept model vectors with proper parameter deduplication
- **121+ tests pass** including unit, composition, and joint inference integration tests

### Key infrastructure added in Step 2

- Multi-model `observe(sol, times, models::Vector)` — concatenates species from all models
- Multi-model `posterior_predictive(models::Vector, chain)` — uses deduplicated params via `unique_params`
- Multi-model `check_identifiability(models::Vector, prob, times)` — validates structural identifiability of composed system
- Refactored `_run_posterior_predictive` to accept model vectors; single-model callers delegate to vector version
- Joint inference integration test — synthetic twin on TX/TL + metabolism, 8-parameter recovery at 90% CI
- Step 2 demo script with wall-clock benchmark and diagnostic figures

### Files changed/added

- `src/inference.jl` — multi-model dispatches for `observe`, `posterior_predictive`, `check_identifiability`, `_run_posterior_predictive` refactor
- `test/test_inference.jl` — unit tests for composed-model `observe` and `check_identifiability`
- `test/test_joint_inference.jl` — new integration test for joint NUTS recovery
- `test/runtests.jl` — registered joint inference integration test
- `examples/step2_joint_demo.jl` — end-to-end demo script
- `examples/run_step2_demo.slurm` — Slurm submission script

## Next step

**Step 3 — Level 1 boundary protocol (sequential conditioning).** Infer the ODE block via NUTS, take posterior samples, condition the SSA block on those posteriors, run ABC-SMC for remaining stochastic parameters. This produces a cut posterior where information flows one direction only.

Main tasks:

- New function: `boundary_condition(ode_chain, ssa_model)` or extend `infer()` to accept a prior source
- Replace ABC-SMC priors for shared parameters with KDE/empirical distribution from ODE posterior
- Test: SSA posteriors tighten when conditioned on ODE posteriors vs uninformed priors

See `docs/plans/2026-04-13-v1-minimal-demo.md` for the full 4-step plan.
