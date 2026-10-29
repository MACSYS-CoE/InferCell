# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

Whole-cell models mix formalisms. Metabolism is a system of ordinary differential
equations, gene expression is a stochastic jump process, and the two exchange
state on a clock. InferCell asks whether parameter uncertainty can be propagated
*correctly* across that seam: when the model is fitted to data, do the posteriors
mean what they claim to mean?

InferCell is a compositional inference framework, demonstrated on a **reduced
JCVI-syn3A model** cut from the published syn3A whole-cell model of Thornburg et
al. (2022). The model is built from modules that touch only through a typed
coupling boundary. It is fitted two ways: exactly, by a sampler that treats the
stochastic path as a latent variable, and approximately, with gene expression
replaced by its mean field.

## What's here

- **[The reduced syn3A model](user-guide/reduced-syn3a.md).** Twenty-two
  reactions in four ODE modules (the 21 published metabolic reactions plus a
  lumped tRNA charging step) and seventeen genes in three jump modules, which
  swap state every simulated second and rate constants every sixty.
- **A composition framework.** Each module declares the states it owns and its
  coupling edges against one [interface contract](api/corea-interface.md), and
  the resolver checks the composition before anything runs. Parameters carry
  their upstream provenance, and every departure from the published model is
  enumerated.
- **[The hybrid driver](api/handshake.md).** It steps the ODE and jump blocks in
  lock-step, reports conservation, clipping and growth from the outside, and
  replays a recorded jump path exactly.
- **The exact route, [a conditional sampler](api/conditional-sampler.md).** A
  three-block Gibbs chain over the ODE parameters, the stochastic rate constants
  and the jump path. The path is updated by conditional SMC with exact
  transcript bridges.
- **The mixed route, [a mean-field description](api/corea-gene-expression.md).**
  Gene expression swapped for its mean field and the whole model run as one ODE,
  with exact transcript likelihoods, an ensemble-corrected bulk term and NUTS.
  `infer(models, ds; route = :exact)` or `route = :mixed` runs either route on an
  M0 dataset ([Inference](api/inference.md)).
- **[Synthetic data](api/corea-observables.md).** Observables, noise models and a
  dataset generator that records its truth, seeds and reductions.
- **Toy models.** The [ODE](user-guide/ode-inference.md) and
  [stochastic](user-guide/ssa-inference.md) models the package started from, with
  a single `infer()` call.

## Status

Pre-release research code.

- **Done and tested:** the model, the framework, the hybrid driver and both
  inference routes.
- **Calibrated on M0,** the two-gene reference system, with six parameters.
  - **The exact route** passes simulation-based calibration (SBC) on all 6
    parameters over 50 prior-drawn datasets (V7).
  - **The mixed route:** On M0 the mixed route passes SBC on all 6 parameters
    over V7's 50 datasets (job 19292232), at a median 6.57 CPU-h per posterior
    against 507.5–678 for the exact sampler. Its largest mean bias against the
    exact posteriors (V7) is 0.095 ± 0.050 SD (krnadeg). On the single dataset
    m0_v6, S_0779, krnadeg and ENO sit 0.42–0.56 SD from the exact reference
    (V6), about 0.3 SD above the matched noise floor, and the cause is open.
- **In progress, for v0.3.0:** `predict`, one entry point for forward and
  intervention predictions. When it lands, v0.3.0 will claim (spec §1):

  > InferCell supports exact stochastic and mean-field descriptions across a
  > coupled ODE–jump boundary. On M0, both inference routes pass
  > simulation-based calibration, with small declared discrepancies.
  > Posterior-updated parameters can be combined modularly with
  > viability-conditioned uncertainty in the remaining parameters to produce
  > forward and intervention predictions on the reduced model.

  The first two sentences hold now; the third waits on `predict`.
- **Not claimed:** an all-parameter posterior, calibrated inference on the full
  reduced model, or a biological finding. Predictions on the reduced model are
  unscored demonstrations of the method; only M0 is scored.
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
