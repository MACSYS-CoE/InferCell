# SSA inference (ABC-SMC)

!!! note "Stub page"
    This tutorial is a placeholder for v0.0.1. The narrative will fill out as the framework stabilises; for now, the runnable example below mirrors the [Getting started](../getting-started.md) SSA path.

Sub-models with `inference_mode = :simulation` are routed to an Approximate Bayesian Computation Sequential Monte Carlo (ABC-SMC) sampler. Summary statistics are computed across replicate trajectories rather than from a single solve.

## Minimal example

```julia
using InferCell

model  = StochasticGeneExpression()
prob   = build_problem(model; tspan=(0.0, 50.0))
trajs  = [solve(prob, SSAStepper(); saveat=0:2.5:50) for _ in 1:200]
data   = observe(trajs, 0:2.5:50, model)

result = infer(model, data;
               n_samples=300,        # ABC particle count
               n_populations=5,      # SMC tolerance schedule
               tspan=(0.0, 50.0))
```

`result` is an [`ABCPosterior`](../api/inference.md) — a weighted-particle posterior.

## What to look at next

- [`abc_smc`](../api/inference.md) — the underlying sampler if you want to bypass `infer()`.
- [`compute_summary_stats`](../api/inference.md) — the default summary used for distance computation.
- [Boundary protocol](boundary-protocol.md) — how to condition this block on an ODE posterior.
