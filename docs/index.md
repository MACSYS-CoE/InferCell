# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

InferCell is a minimal whole-cell modelling framework designed so that **Bayesian inference works from the start**. Sub-models — transcription/translation, metabolism, stochastic gene expression — are composed into a single computational graph in pure Julia. AD, likelihood evaluation, and posterior calibration fall out of this design rather than being bolted on.

## What's here

- **One `infer()` call** that dispatches to NUTS for differentiable blocks and a simulation-based backend for stochastic ones (ABC-SMC today, moving to a learned differentiable likelihood). The user writes the same code regardless of formalism.
- **A boundary protocol** that propagates posteriors across heterogeneous block boundaries with explicit uncertainty. Sequential conditioning (Level 1) and iterative message-passing (Level 2) are both supported.
- **Closed-loop biology**: a syn3A-flavoured toy cell with ODE metabolism, ODE bulk gene expression, and an SSA bursty regulator coupled by shared ATP/NTP/AA pools and a metabolic-enzyme feedback.

## Where to go next

- New to the project? Start with [Getting started](getting-started.md) for installation and a runnable example.
- Curious about the framing and why inference, not just simulation? The design log and research notes live under `dev/` in the repository.
- Building a model? Skim the [User guide](user-guide/ode-inference.md) tutorials, then jump into the [API reference](api/index.md).

## Status

Pre-release. The framework runs end-to-end on synthetic data (Steps 1–5 of the v0.0.1 plan complete on `main`). Real-data calibration against syn3A datasets is v0.0.2.
