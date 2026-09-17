```@meta
CurrentModule = InferCell
```

# Parameters

Source: [`src/parameters.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/parameters.jl).

`InferParameter` is the unit of parameter declaration shared by every sub-model. The role-filtering helpers below let the inference layer separate kinetic rates, initial conditions, and observation parameters without case-by-case logic.

## The declaration

## Provenance

`ParameterSource` is carried per value rather than resolved once by hand, because the same cross-file ambiguity has already produced two errors of record in the source model — see [Core A′ interface contract](corea-interface.md) for the loader that produces these.

## Role filtering

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["parameters.jl"]
```
