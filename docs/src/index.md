# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

Whole-cell models mix formalisms. Metabolism is a system of ordinary differential
equations, gene expression is a stochastic jump process, and the two exchange
state on a clock. InferCell asks whether parameter uncertainty can be propagated
*correctly* across that seam: when the model is fitted to data, do the posteriors
mean what they claim to mean?

It answers on a **reduced JCVI-syn3A model**, cut from the published syn3A
whole-cell model of Thornburg et al. (2022). The model is built from modules that
touch only through a typed coupling boundary, and fitted by a sampler that treats
the stochastic path as a latent variable.

## What's here

- **[The reduced syn3A model](user-guide/reduced-syn3a.md).** Twenty-two metabolic
  reactions in four ODE modules and seventeen genes in three jump modules, which
  swap state every simulated second and rate constants every sixty.
- **A composition framework.** Each module declares the states it owns and its
  coupling edges against one [interface contract](api/corea-interface.md), and
  the resolver checks the composition before anything runs. Parameters carry
  their upstream provenance, and every departure from the published model is
  enumerated.
- **[The hybrid driver](api/handshake.md).** It steps the ODE and jump blocks in
  lock-step, reports conservation, clipping and growth from the outside, and
  replays a recorded jump path exactly.
- **[A conditional sampler](api/conditional-sampler.md).** A three-block
  Gibbs chain over the ODE parameters, the stochastic rate constants and the jump
  path. The path is updated by conditional SMC with exact transcript bridges.
- **[Synthetic data](api/corea-observables.md).** Observables, noise models and a
  dataset generator that records its truth, seeds and reductions.
- **Toy models.** The [ODE](user-guide/ode-inference.md) and
  [stochastic](user-guide/ssa-inference.md) models the package started from, with
  a single `infer()` call.

## Status

Pre-release research code.

- **Done and tested:** the model, the framework, the hybrid driver and each
  sampler block.
- **Running:** the long reference chain and the simulation-based calibration of
  the sampler on a two-gene reference system. Calibration has not yet passed, and
  this release makes no claim that the posteriors are calibrated.
- **Not started:** recovery and coverage on the full reduced model.
- **Not in scope:** real data. The reduced model does not reproduce syn3A's
  doubling time, metabolite levels or proteome, and is not to be compared with
  them.

## Where to go next

- [Getting started](getting-started.md) installs the package and runs the model.
- [The reduced syn3A model](user-guide/reduced-syn3a.md) says what the model
  contains and what it leaves out.
- The [API reference](api/index.md) has every exported symbol's docstring.
- The research specification, `spec/spec.md` in the repository, carries the
  claim, the validity checks and the decisions of record.
