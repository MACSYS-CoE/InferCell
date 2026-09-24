```@meta
CurrentModule = InferCell
```

# Inference

Source: [`src/inference.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/inference.jl), [`src/abc_smc.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/abc_smc.jl), [`src/summary_statistics.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/summary_statistics.jl), [`src/likelihoods.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/likelihoods.jl).

## Entry point

[`infer`](@ref) selects NUTS for sub-models with `inference_mode = :differentiable` and ABC-SMC for `:simulation`. `priors_override` is forwarded to [`build_turing_model`](@ref) on the NUTS path — the hook the boundary protocol uses to feed conditioned priors back into NUTS. The mode is read from the whole composition: every module must declare the same one, and a hybrid ODE/jump composition is refused by name, because it builds a handshake driver that neither backend can run.

## Synthesising data

## Posterior predictive

## Identifiability

## ABC-SMC

## Summary statistics

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["inference.jl", "abc_smc.jl", "summary_statistics.jl", "likelihoods.jl"]
```
