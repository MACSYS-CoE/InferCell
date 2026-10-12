```@meta
CurrentModule = InferCell
```

# Boundary protocol

!!! note "Toy models"
    This page summarises the boundary-protocol API. For the calibration behaviour of the two levels, see `dev/notes/2026-05-22-calibration-curve.md` in the repo; for the inference design they sit inside, `dev/notes/modular-bayesian-inference-heterogeneous-modules.md`.

NUTS-ODE and ABC-SMC-SSA cannot be jointly composed in a single sampler. InferCell instead **propagates posteriors across the heterogeneous block boundary** as an explicit step.

## Levels

- **Level 1 — sequential conditioning** ([`sequential_infer`](@ref)). One forward pass: NUTS on the ODE block → fit KDE priors on the shared parameters → ABC-SMC on the SSA block with the conditioned priors. Cheap, robust, sufficient when shared parameters are well-identified by the ODE data.
- **Level 2 — iterative message-passing** ([`iterative_infer`](@ref)). Iterates ODE ↔ SSA until the ODE posterior stops changing (summed KL divergence across shared parameters drops below `kl_tol`, or `max_iters` is reached). Tightens shared-parameter posteriors but can slightly over-propagate confidence.

## Minimal example

```julia
using InferCell

# Composed ODE block (e.g. TX/TL + metabolism); composed SSA block (bursty regulator)
ode_models = [TranscriptionTranslation(), LightMetabolism()]
ssa_models = [BurstyGeneExpression()]

# Generate ode_data and ssa_data via the usual `observe()` workflow

result = sequential_infer(ode_models, ode_data, ssa_models, ssa_data;
                          ode_kwargs=(; n_samples=1000),
                          ssa_kwargs=(; n_samples=300, n_populations=5),
                          verbose=true)

# result.ode_chain — MCMCChains.Chains
# result.ssa_posterior — ABCPosterior
# result.boundary — (priors, param_names, conditioned)
```

For the iterative version, swap `sequential_infer` for `iterative_infer` and add `max_iters=5, kl_tol=0.05`.

## Calibration

Both levels are proofs of concept of passing a posterior between the two toy
backends, not calibrated methods. A single-seed run of `iterative_infer`
over-tightened at 95% (two of four ODE parameters missed truth at 95% and were
covered at 99%), and no multi-seed calibration has been run. The reduced syn3A
model has no shared parameter to pass, and uses the
[conditional sampler](../api/conditional-sampler.md) instead.
