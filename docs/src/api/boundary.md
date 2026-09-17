```@meta
CurrentModule = InferCell
```

# Boundary

Source: [`src/boundary.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/boundary.jl).

The boundary module implements posterior propagation across heterogeneous block boundaries (ODE ↔ SSA). See [User guide → Boundary protocol](../user-guide/boundary-protocol.md) for the design narrative.

## Level 1 — cut posterior, single pass

## Level 2 — iterative message passing

[`iterative_infer`](@ref) iterates ODE ↔ SSA until the ODE posterior stops changing (summed KL across shared parameters drops below `kl_tol`) or `max_iters` is reached.

## Conditioning

## KL diagnostics

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["boundary.jl"]
```
