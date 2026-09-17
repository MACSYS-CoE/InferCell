```@meta
CurrentModule = InferCell
```

# Models

Built-in sub-models. Each is a concrete subtype of [`AbstractSubModel`](@ref). Source files under [`src/models/`](https://github.com/MACSYS-CoE/InferCell/tree/main/src/models).

## Deterministic blocks

## Stochastic blocks

`BurstyGeneExpression` is load-bearing for the boundary protocol because it has no ODE equivalent.

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["models/transcription_translation.jl", "models/stochastic_gene_expression.jl", "models/bursty_gene_expression.jl", "models/light_metabolism.jl", "models/tier_b_metabolism.jl"]
```
