```@meta
CurrentModule = InferCell
```

# Prediction

Source: [`src/posterior_draws.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/posterior_draws.jl).

## Reading a posterior

[`posterior_draws`](@ref) reads any posterior InferCell writes into one
[`PosteriorDraws`](@ref): the exact route's V7 replicate files, the mixed route's
B7 chain files, the `Chains` that `infer(models, ds; route)` returns, and the toy
models' `Chains` and `ABCPosterior`. Every column is mapped to its parameter by
name, and every value is held on the natural scale, ready to write into a model.
The bulk noise σ_b is kept in its own field: it is observation noise, not a cell
parameter, and a forward simulation never reads it.

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["posterior_draws.jl"]
```
