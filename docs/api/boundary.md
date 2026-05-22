# Boundary

Source: [`src/boundary.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/boundary.jl).

The boundary module implements posterior propagation across heterogeneous block boundaries (ODE ↔ SSA). See [User guide → Boundary protocol](../user-guide/boundary-protocol.md) for the design narrative.

## `sequential_infer`

```julia
sequential_infer(ode_models, ode_data, ssa_models, ssa_data;
                 ode_kwargs=(;), ssa_kwargs=(;),
                 boundary_method=:kde, verbose=false)
```

Level-1 boundary protocol (cut posterior, single-pass). Three stages:

1. NUTS on the ODE block.
2. Fit `boundary_method` priors (`:kde` or `:normal`) from the ODE posterior for shared parameters.
3. ABC-SMC on the SSA block with conditioned priors.

Returns `(ode_chain, ssa_posterior, boundary)`.

## `iterative_infer`

```julia
iterative_infer(ode_models, ode_data, ssa_models, ssa_data;
                ode_kwargs=(;), ssa_kwargs=(;),
                boundary_method=:kde,
                max_iters=5, kl_tol=0.05,
                n_backward_samples=1000,
                keep_history=false, verbose=false, rng=...)
```

Level-2 boundary protocol. Iterates ODE ↔ SSA message-passing until the ODE posterior stops changing (summed KL across shared parameters drops below `kl_tol`) or `max_iters` is reached.

Returns `(ode_chain, ssa_posterior, n_iters, kl_trace, ode_history, ssa_history, boundary)`. Histories are empty when `keep_history=false`.

## `boundary_condition`

```julia
boundary_condition(chain, ssa_models; method=:kde)                              # forward (ODE → SSA)
boundary_condition(posterior::ABCPosterior, ode_models; method=:kde,
                   n_samples=1000, rng=...)                                     # backward (SSA → ODE)
```

Build conditioned priors for one block from the other's posterior. Returns `(priors, param_names, conditioned)` — a parallel arrays NamedTuple where `conditioned[i]` indicates whether parameter `param_names[i]` was actually found in the upstream posterior (un-found parameters keep their model-declared prior).

## `KDEPrior`

```julia
KDEPrior(samples::Vector{Float64})
```

Univariate `ContinuousUnivariateDistribution` fitted to posterior samples via KernelDensity.jl. Satisfies the `rand(rng, d)` and `pdf(d, x)` interface needed by ABC-SMC, plus the `minimum` / `maximum` / `insupport` interface needed by Turing's `Bijectors.jl` for unconstrained sampling.

## KL diagnostics

```julia
kl_divergence(samples_p, samples_q; n_grid=200, q_floor=1e-12)
chain_kl(chain_new, chain_old, param_names)
```

`kl_divergence` estimates KL(p ‖ q) by KDE-smoothing both sample sets on a shared grid and integrating numerically; `q_floor` keeps the result finite when supports don't overlap. `chain_kl` sums per-parameter KLs across `param_names`, skipping any parameter not present in both chains — used as the convergence criterion in [`iterative_infer`](#iterative_infer).
