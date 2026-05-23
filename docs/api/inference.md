# Inference

Source: [`src/inference.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/inference.jl), [`src/abc_smc.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/abc_smc.jl), [`src/summary_statistics.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/summary_statistics.jl), [`src/likelihoods.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/likelihoods.jl).

## `infer`

```julia
infer(models, data; sampler=NUTS(), n_samples=1000, kwargs...)
infer(model,  data; kwargs...)
```

Dispatch entry point. Selects NUTS for sub-models with `inference_mode = :differentiable` and ABC-SMC for `:simulation`. Returns an `MCMCChains.Chains` (NUTS) or [`ABCPosterior`](#abcposterior) (ABC-SMC).

`priors_override::Dict{Symbol, Distribution}` is forwarded to [`build_turing_model`](#build_turing_model) on the NUTS path — the hook the boundary protocol uses to feed conditioned priors back into NUTS.

## `build_turing_model`

```julia
build_turing_model(models, data, prob; solver=Tsit5(), sensealg=ForwardDiffSensitivity(),
                   priors_override=nothing)
```

Construct the underlying Turing.jl `@model` used by the NUTS path. Useful when you want to call `sample()` directly, mix in custom samplers, or interrogate the prior structure.

Returns `(turing_model, param_names::Vector{Symbol})`. The order in `param_names` is: model free params first, then observation free params.

## Synthesising data

```julia
observe(sol, times, model_or_models; sigma=0.1, rng=...)         # ODE: noisy point sample
observe(trajectories, times, model_or_models)                    # SSA: mean of replicates
```

Builds an [`ObservedData`](#observeddata).

## Posterior predictive

```julia
posterior_predictive(model_or_models, chain_or_posterior;
                     n_samples=100, tspan, saveat=nothing)
```

Draws forward simulations from the posterior. Accepts either an MCMCChains chain (NUTS) or an [`ABCPosterior`](#abcposterior) (resampled by particle weights). Returns a [`PosteriorPredictive`](#posteriorpredictive). `tspan` is required.

## Identifiability

```julia
check_identifiability(model_or_models, prob, times; solver=Tsit5())
```

Local structural-identifiability diagnostic. Builds the forward-map Jacobian via `ForwardDiff` and returns `(rank, n_params, n_obs, full_rank, jacobian)`. `full_rank == true` means every parameter is locally identifiable from the chosen observation pattern.

## ABC-SMC

### `ABCPosterior`

Weighted-particle posterior returned by ABC. Fields:

- `particles::Matrix{Float64}` — shape `(n_params, n_particles)`
- `weights::Vector{Float64}` — sums to 1
- `param_names::Vector{Symbol}` — orders the rows of `particles`
- `tolerance::Float64`, `n_populations::Int` — final SMC state

### `abc_smc`

```julia
abc_smc(simulate, observed_stats, priors, param_names;
        n_particles=1000, n_populations=10, alpha=0.5, verbose=false, rng=...)
```

ABC with adaptive Sequential Monte Carlo. `simulate(theta)` must return a summary-statistic vector comparable to `observed_stats`. Tolerance schedule is adaptive: the new tolerance is the `alpha`-quantile of the previous population's distances.

## Summary statistics

```julia
compute_summary_stats(trajectories, species; times=nothing)
summary_distance(s1, s2)
```

`compute_summary_stats` reduces an SSA ensemble to a single statistic vector: with `times=nothing` returns `[means; vars; fanos]` of the final state; otherwise returns concatenated per-time-point means. `summary_distance` is Euclidean.

## Data containers

### `ObservedData`

```julia
ObservedData(times, observations, species)
```

`times::Vector{Float64}`, `observations::Matrix{Float64}` (`species × T`), `species::Vector{Symbol}`. Produced by [`observe`](#synthesising-data), consumed by [`infer`](#infer).

### `PosteriorPredictive`

```julia
PosteriorPredictive(solutions, times, species)
```

Container for forward simulations drawn from a posterior. Built by [`posterior_predictive`](#posterior-predictive).
