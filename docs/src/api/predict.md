```@meta
CurrentModule = InferCell
```

# Prediction

Source: [`src/posterior_draws.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/posterior_draws.jl), [`src/predict_completion.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_completion.jl).

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

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["posterior_draws.jl", "predict_completion.jl"]
```
