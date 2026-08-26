# Positioning

> Date: 2026-05-12

## What InferCell is

A composable, inference-first whole-cell modelling framework in Julia. The goal is a working hybrid continuous-discrete model of a minimal cell where Bayesian parameter inference works end-to-end across heterogeneous blocks (ODE, SDE, SSA, discrete events) in one compiled computational graph.

The framing is sequenced: every shipped version through v0.0.x is a methods contribution proven on a deliberately minimal cell, walking toward the longer-term vision of a whole-cell model that does inference — concretely, calibrated UQ over a well-mixed reduction of the Luthey-Schulten / Thornburg minimal cell. See the research notes under `dev/notes/` in the repository for the full vision and current scope.

## The design discipline

For every module boundary, ask: **how does inference cross this interface?** If you can't answer that, the design needs to change.

This single question constrains a great deal. It rules out cross-language orchestration (gradients don't cross language boundaries). It rules out black-box sub-models (boundary posteriors need explicit structure). It rules out point-estimate handoffs between blocks (uncertainty has to propagate). It demands that every block declare its parameters, observables, and state in a form the inference machinery can read.

There is a second, sharper form of the same question: **how does inference cross this interface when the module on the other side might be the wrong model?** A boundary that hands over a point estimate fails the first version. A boundary that hands over a posterior conditioned on one model structure, with no record of which structure, passes the first and fails the second.

Most existing WCM frameworks (Karr 2012, CovertLab, Vivarium, Lattice Microbes) are simulation-first: forward dynamics work well, inference is bolted on later as a bespoke project per model. InferCell inverts the priority. The cost is that some simulation patterns become harder; the payoff is that calibration, sensitivity analysis, model reduction, and model comparison are first-class rather than per-model engineering efforts.

## What falls out of the discipline

Four named capabilities are downstream consequences of inference-first design — not separate things we have to build.

### Parameter inference

Bayesian posteriors over kinetic parameters, pool sizes, and initial conditions from time-series data. The `infer()` API dispatches to NUTS (gradient-based) for differentiable models and a simulation-based backend for stochastic ones, behind one call — ABC-SMC today, moving to a learned differentiable likelihood (NLE) so the stochastic block can join the graph. Posteriors propagate across the ODE/SSA boundary via explicit boundary protocols: Level 1 (sequential conditioning) and Level 2 (iterative message-passing) are implemented as a proof-of-concept; the EP-over-trajectories production form and a joint-NUTS oracle are the forward direction. Particle MCMC is ruled out.

This is the headline capability and the one v0.0.1 is engineered around.

### Sensitivity analysis

Local and global parameter sensitivities. The AD machinery that powers NUTS gradients also powers forward and adjoint sensitivities for free. Two practical uses:

- **Experimental design.** Which parameters does this observable actually constrain? Answering this *before* committing to a data collection campaign is cheap.
- **Inference diagnostics.** Which inference problems are well-posed before sampling? Sloppy directions in parameter space show up in the sensitivities; you can see them without burning sampler compute.

### Model reduction

Identifiability checks, sloppy-parameter analysis, and posterior-contraction-guided simplification. The framework can tell you which parameters your data can't constrain — and, by extension, where the model has more degrees of freedom than the biology warrants. Reduction by posterior collapse is more principled than reduction by hand-waved time-scale arguments.

The connection to inference is direct: a parameter the posterior doesn't move away from its prior is one the data doesn't see. That's a reduction signal.

### Model comparison

Comparison across competing model structures — constitutive vs bursty transcription, with vs without enzyme feedback, one regulator vs many. Each candidate is an `infer()` call, and the comparison is a calculation on the resulting posteriors and predictive densities.

Two qualifications, because this capability is easier to promise than to deliver.

**The criterion is predictive, not evidential.** Marginal likelihood is acutely sensitive to prior width, and priors on kinetic parameters here are diffuse and semi-arbitrary by construction. Expected log predictive density — cross-validated, or leave-one-stream-out — is the robust criterion and is the same quantity the robustness layer already needs for choosing per-module influence parameters. Bayes factors are not the target.

**Not every backend produces a comparable number.** A NUTS log-density is comparable across structures fitted to the same data; an ABC-SMC posterior is not, because its implied evidence depends on the tolerance and the summary statistics. The design commitment is therefore that each module declares its normalisation semantics, rather than that likelihood values from different backends be assumed differenceable.

Structural uncertainty is not quantified in v0.0.x — parameter inference across heterogeneous blocks is the prerequisite for every model-comparison quantity. What is committed to now is narrower: that structural identity and normalisation semantics be recorded from the start, so the comparison layer can be added later without rebuilding the module interfaces. The full argument is in `dev/notes/model-uncertainty-and-selection.md`.

## What InferCell is not

- **Not a simulation platform.** Forward simulation is necessary but not sufficient. If forward simulation is the only thing you need, Lattice Microbes, CobraToolbox, or Catalyst.jl directly are better choices.
- **Not a genome-scale model.** v0.0.1 has ~5 named genes and ~5 metabolic reactions. The architectural claim is that the inference machinery scales with the model; the model itself stays deliberately minimal until each addition is justified by data.
- **Not yet validated against real data.** Every shipped version so far is synthetic. The route to real syn3A is the three-step ladder scoped in `dev/notes/reduced-syn3a-scoping.md` (Step 1: the Core A' reduction; Step 2: the published well-stirred model; Step 3: the 4DWCM). Until then, the project is a methods contribution, not a biological one.
- **Not a replacement for domain expertise.** The framework makes inference easy; deciding *which inference question to ask* is still the user's job.

## Who this is for

- **Computational cell biologists** who want to do Bayesian parameter estimation on their models without writing per-model inference code.
- **Methods researchers** working on hybrid inference, boundary protocols, or simulation-based inference, who want a non-trivial test case grounded in real biology.
- **WCM developers** who want a modular, AD-friendly substrate for composing sub-models from across the literature.

## Related work

See README for primary WCM citations. In short:

- **Karr et al. 2012** — established the sub-model composition paradigm InferCell follows. Not inference-first.
- **Macklin et al. 2020 (CovertLab)** — scale and complexity benchmark for E. coli WCMs. Not inference-first.
- **Thornburg et al. 2022, Luthey-Schulten et al. 2026 (4DWCM)** — spatial state-of-the-art for syn3A. Not inference-first.
- **Agmon et al. 2022 (Vivarium)** — composability via orchestration across languages rather than a shared compiled graph. InferCell departs in choosing one Julia compilation unit over message-passing across processes.

The repeated "not inference-first" is the gap in the field this project exists to address.
