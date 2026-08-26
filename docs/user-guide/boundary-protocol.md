# Boundary protocol

!!! note "Stub page"
    This page summarises the boundary-protocol API. For the calibration behaviour of the two levels, see `dev/notes/2026-05-22-calibration-curve.md` in the repo; for the inference design they sit inside, `dev/notes/modular-bayesian-inference-heterogeneous-modules.md`.

NUTS-ODE and ABC-SMC-SSA cannot be jointly composed in a single sampler. InferCell instead **propagates posteriors across the heterogeneous block boundary** as an explicit step.

## Levels

- **Level 1 — sequential conditioning** ([`sequential_infer`](../api/boundary.md)). One forward pass: NUTS on the ODE block → fit KDE priors on the shared parameters → ABC-SMC on the SSA block with the conditioned priors. Cheap, robust, sufficient when shared parameters are well-identified by the ODE data.
- **Level 2 — iterative message-passing** ([`iterative_infer`](../api/boundary.md)). Iterates ODE ↔ SSA until the ODE posterior stops changing (summed KL divergence across shared parameters drops below `kl_tol`, or `max_iters` is reached). Tightens shared-parameter posteriors but can slightly over-propagate confidence.

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

A single-seed run of `iterative_infer` was found to over-tighten by a small amount at 95% CI (2/4 ODE parameters miss truth at 95% but cover at 99%). A multi-seed calibration sweep is tracked in [Issue #11](https://github.com/MACSYS-CoE/InferCell/issues/11) and should land before the v0.0.1 publication figure.
