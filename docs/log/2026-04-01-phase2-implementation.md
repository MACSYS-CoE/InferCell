# Phase 2 Implementation Log: Stochastic Gene Expression + ABC-SMC

**Date:** 2026-04-01
**Branch:** `phase-2-stochastic`
**Status:** Implementation complete, pending full validation

## What was built

### New capabilities
1. **StochasticGeneExpression sub-model** — Same TX/TL biology as Phase 1, implemented as Gillespie SSA via JumpProcesses.jl (4 reactions: mRNA production, mRNA degradation, translation, protein degradation)
2. **ABC-SMC inference engine** — Hand-rolled sequential Monte Carlo for simulation-based inference (~100 lines). Supports adaptive tolerance, Gaussian perturbation kernels, importance weighting.
3. **Summary statistics framework** — Computes mean, variance, Fano factor from replicate SSA trajectories. Euclidean distance metric.
4. **Inference mode dispatch** — `infer()` now routes to NUTS or ABC-SMC based on `inference_mode(model)` (`:differentiable` or `:simulation`).
5. **GitHub Actions CI** — Unit tests + integration tests (gated) on Julia 1.10.

### Architecture changes
- `build_problem` dispatches on collective formalism (`:ode` → `ODEProblem`, `:jump` → `JumpProblem`). Mixed formalism errors explicitly (deferred to v0.0.2).
- `reactions(m::AbstractSubModel)` added to interface for jump models.
- `ode_free_params` renamed to `model_free_params` (alias preserved for backward compat).
- `ABCPosterior` type for simulation-based inference results (weighted particles, not MCMC chains).
- `posterior_predictive` dispatches on both `Chains` and `ABCPosterior`.

### Files added
- `src/models/stochastic_gene_expression.jl` — SSA sub-model
- `src/summary_statistics.jl` — Summary stats computation + distance
- `src/abc_smc.jl` — ABC-SMC algorithm + ABCPosterior type
- `test/test_stochastic_ge.jl` — Unit tests for SSA model
- `test/test_summary_stats.jl` — Summary stats tests
- `test/test_abc_smc.jl` — ABC-SMC algorithm tests (toy Gaussian problem)
- `test/test_stochastic_inference.jl` — Integration test: ABC-SMC parameter recovery
- `examples/phase2_demo.jl` — ODE/NUTS vs SSA/ABC-SMC comparison demo
- `.github/workflows/CI.yml` — GitHub Actions CI

### Files modified
- `src/InferCell.jl` — New imports, includes, exports
- `src/interface.jl` — Added `reactions()` default
- `src/parameters.jl` — `model_free_params` + alias
- `src/orchestrator.jl` — Formalism dispatch, `_build_jump_problem`, `_build_u0_integer`
- `src/inference.jl` — Inference mode dispatch, `_infer_abc`, ABC `posterior_predictive`
- `Project.toml` — JumpProcesses dep, [compat] section
- `README.md` — Quickstart section
- `test/runtests.jl` — Include new test files

### Key design decisions
1. **JumpProcesses.jl over hand-rolled Gillespie** — SciML integration, multiple algorithms, battle-tested. AD irrelevant for SBI path.
2. **Hand-rolled ABC-SMC over library** — Full control, understanding of every line. ApproxBayes.jl/KissABC.jl are limited. No Python deps.
3. **Same biology in both formalisms** — Controlled comparison: same ground truth, two inference backends. Most compelling for validation.
4. **ABCPosterior as separate type** — Honest about what ABC-SMC produces (weighted particles) vs MCMC chains.

## What this proves

- The InferCell architecture supports both differentiable (NUTS) and simulation-based (ABC-SMC) inference through the same `infer()` API.
- The sub-model interface (`AbstractSubModel`) is flexible enough to express both ODE and jump process dynamics.
- The orchestrator dispatches correctly on formalism without breaking existing ODE functionality.

## What this doesn't prove (deferred)

- Multi-module composition (ODE + SSA modules in the same inference call)
- Boundary protocol (uncertainty propagation across blocks)
- Scaling to >4 parameters
- Real experimental data
