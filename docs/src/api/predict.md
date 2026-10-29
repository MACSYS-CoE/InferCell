```@meta
CurrentModule = InferCell
```

# Prediction

Source: [`src/posterior_draws.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/posterior_draws.jl), [`src/predict_completion.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_completion.jl), [`src/predict_forward.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_forward.jl), [`src/predict_viability.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_viability.jl).

## Reading a posterior

[`posterior_draws`](@ref) reads any posterior InferCell writes into one
[`PosteriorDraws`](@ref): the exact route's V7 replicate files, the mixed route's
B7 chain files, the `Chains` that `infer(models, ds; route)` returns, and the toy
models' `Chains` and `ABCPosterior`. Every column is mapped to its parameter by
name, and every value is held on the natural scale, ready to write into a model.
The bulk noise σ_b is kept in its own field: it is observation noise, not a cell
parameter, and a forward simulation never reads it.

## Completing the held parameters

On the reduced model, M0's posterior covers five parameters, and every other free
parameter is held: no data reaches it. [`held_parameters`](@ref) lists them for a
composition, skipping the slots the driver overwrites (the `enz_R_*` catalytic
slots), the reverse constants, and point priors, and says why each was skipped.
[`draw_held`](@ref) draws them from their priors, which `width` widens or narrows
about their medians for the sensitivity runs. Reverse constants are derived from
the whole set that runs, posterior and held values together, by
[`haldane_values`](@ref), so every reaction keeps its equilibrium constant.

## Forward descriptions and δ

[`forward_cells`](@ref) runs a composition at given parameter writes, read for
the predicted quantity Q:
- `forward = :hybrid` runs the reference composition, one stochastic cell per
  seed, and each seed fixes its cell;
- `forward = :mean_field` runs route (ii), one ODE, with no cells and no
  single-cell distribution.

Q is read from the run's metered lactate export, so every composition is
[`metered`](@ref) first: `:lactate_out` is the lactate exported in particles and
`:exported` is the fraction of the lactate formed that left the cell. A function
of the run is accepted too. Each run is kept as a [`CellRun`](@ref): its Q, its
clips, the smallest particle count of each state, and its failure, if it failed,
with the reason.

The mean field has no δ for Q of its own. The mixed route's δ covers only its 17
panel pools. [`delta_correction`](@ref) computes one from a hybrid ensemble at a
reference point, once, and keeps it as a [`DeltaCorrection`](@ref) to reuse.

## Viability

A viability rule V is data, a [`ViabilityRule`](@ref), and nothing in the code
fixes its thresholds. V is frozen in M4, before any prediction is read.
[`cell_status`](@ref) puts each run in one of five categories, which are not
interchangeable: failed, starved, below a floor, clipped, or clean. V is a
property of θ: [`assess_viability`](@ref) accepts a parameter set when enough of
its cells are viable, and otherwise keeps its count of cells per category as the
reason. Over the crossed matrix of posterior draws × held draws,
[`viability_weights`](@ref) gives both conditioned variants of the methods note
as weights: (b) keeps θ_I's marginal as the posterior's, and (a) lets the filter
act on θ_I. Each comes with its acceptance rate and effective sample size.

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["posterior_draws.jl", "predict_completion.jl", "predict_forward.jl", "predict_viability.jl"]
```
