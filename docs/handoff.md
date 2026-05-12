# Handoff

## Current state

Branch `step3-boundary-protocol` implements Steps 1-3 of the [v0.0.1 plan](plans/2026-05-12-v0.0.1-scoping.md).

### What's working

- **TX/TL module** — ODE model with NUTS inference, posterior predictive, identifiability
- **Stochastic gene expression** — SSA model with ABC-SMC, automatic dispatch via `infer()`
- **LightMetabolism module** — minimal ODE metabolism (ATP, NTP, AA pools) coupled one-way to TX/TL
- **Shared parameter composition** — `build_problem([TranscriptionTranslation(), LightMetabolism()])` produces a single ODEProblem with 5 states and 7 deduplicated free parameters
- **Joint ODE inference** — `infer([txl, metab], data)` runs NUTS on the composed 5-state system, recovering all 8 parameters (7 ODE + 1 obs noise)
- **Multi-model utilities** — `observe`, `posterior_predictive`, `check_identifiability` all accept model vectors with proper parameter deduplication
- **Level 1 boundary protocol** — `sequential_infer(ode_model, ode_data, ssa_model, ssa_data)` runs the full pipeline: ODE NUTS → KDE boundary conditioning → SSA ABC-SMC with informed priors
- **KDEPrior** — custom `ContinuousUnivariateDistribution` wrapping KDE-fitted posteriors, compatible with ABC-SMC's `rand`/`pdf` interface
- **boundary_condition** — extracts ODE posterior samples, fits KDE (or Normal) priors for shared parameters, passes through original priors for unmatched parameters
- **121+ unit tests pass** plus integration tests for joint inference and sequential inference

### Key infrastructure added in Step 3

- `KDEPrior <: ContinuousUnivariateDistribution` — wraps KernelDensity.jl KDE with `rand()`, `pdf()`, `logpdf()` methods for use as ABC-SMC priors
- `boundary_condition(chain, ssa_models; method=:kde)` — builds conditioned priors from ODE posterior, supports `:kde` and `:normal` methods, handles partial conditioning
- `sequential_infer(ode_models, ode_data, ssa_models, ssa_data)` — end-to-end two-stage pipeline returning `(ode_chain, ssa_posterior, boundary)` named tuple
- Modified `_infer_abc` to accept optional `priors` and `param_names` kwargs for external prior injection
- Sequential inference integration test — synthetic twin verifying parameter recovery and posterior tightening vs uninformed baseline
- Step 3 demo script with conditioned vs unconditioned comparison and diagnostic figures

### Files changed/added

- `src/boundary.jl` — new: KDEPrior, boundary_condition, sequential_infer
- `src/inference.jl` — modified: `_infer_abc` accepts optional `priors`/`param_names` kwargs
- `src/InferCell.jl` — added KernelDensity import, include boundary.jl, new exports
- `Project.toml` — added KernelDensity dependency
- `test/test_boundary.jl` — new: unit tests for KDEPrior, boundary_condition, partial conditioning
- `test/test_sequential_inference.jl` — new: integration test for sequential inference + tightening comparison
- `test/runtests.jl` — registered new test files
- `examples/step3_boundary_demo.jl` — new: end-to-end demo script
- `examples/run_step3_demo.slurm` — new: Slurm submission script

## Next step

**Step 4 — Level 2 boundary protocol (message-passing / bidirectional).** Close the loop so SSA posteriors can feed back to refine ODE inference, enabling iterative convergence of the cut posterior.

Main tasks:

- Implement iterative message-passing loop between ODE and SSA blocks
- Define convergence criterion (e.g., KL divergence between successive posteriors)
- Test: demonstrate convergence and improved parameter recovery over single-pass sequential conditioning

See `docs/plans/2026-05-12-v0.0.1-scoping.md` for the full plan (now 6 steps, with Step 4 next).
