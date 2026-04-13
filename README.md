# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

## Quick start

```julia
using InferCell

# ODE model → automatically uses NUTS (gradient-based)
model = TranscriptionTranslation()
prob  = build_problem(model; tspan=(0.0, 50.0))
sol   = solve(prob, Tsit5(); saveat=0:2.5:50)
data  = observe(sol, 0:2.5:50, model; sigma=0.3)
chain = infer(model, data; n_samples=1000)

# SSA model → automatically uses ABC-SMC (likelihood-free)
model  = StochasticGeneExpression()
prob   = build_problem(model; tspan=(0.0, 50.0))
trajs  = [solve(prob, SSAStepper(); saveat=0:2.5:50) for _ in 1:200]
data   = observe(trajs, 0:2.5:50, model)
result = infer(model, data; n_samples=300, n_populations=5,
               tspan=(0.0, 50.0))
```

Same four steps, same `infer()` call. The framework dispatches to the right backend based on the model's formalism.

Run tests:
```bash
julia --project -e 'using Pkg; Pkg.test()'

# Include integration tests (~2 min, NUTS sampling):
INFERCELL_INTEGRATION_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
```

**Status:** Phase 2 complete. The `infer()` API dispatches to NUTS for differentiable models (ODE) or ABC-SMC for simulation-based models (SSA), selected automatically by model formalism.

## Roadmap

### Done

- **TX/TL module** -- ODE model with NUTS inference, posterior predictive checks, identifiability analysis
- **Stochastic gene expression** -- SSA model with ABC-SMC inference, automatic backend dispatch via `infer()`

### Next

- **Bursting gene expression** -- two-state promoter model (ON/OFF switching) with no ODE equivalent, demonstrating why likelihood-free inference is necessary
- **Light metabolism** -- ODE module (~5-10 reactions) composed jointly with TX/TL, exercising multi-module parameter sharing

### v1 target

Three biological modules connected through an inference graph:
- **DifferentiableBlock** (ODE modules) inferred via NUTS/HMC
- **SimulationBlock** (SSA modules) inferred via ABC-SMC
- **Boundary protocol** (Level 1): sequential conditioning for uncertainty propagation between blocks
- Model selection via Bayes factors (e.g. constitutive vs bursty transcription)

See [`docs/plans/overview.md`](docs/plans/overview.md) for full architecture and design decisions.

## Motivation

Whole-cell models simulate a living cell from its molecular parts. The field has made real progress wiring sub-models of transcription, translation, metabolism, and replication into integrated simulations. But the software architectures that run these simulations were not designed for inference. Existing platforms treat sub-models as black boxes, orchestrated via message passing or multi-language pipelines. Good for running forward simulations. Opaque to parameter estimation, automatic differentiation, and uncertainty quantification. You can't push gradients or likelihoods through boundaries you can't see inside.

The problem is architectural. Simulation-first frameworks treat inference as an afterthought, and hard language boundaries make Bayesian calibration a bespoke, fragile process for every new model configuration.

## What we're building

InferCell is a minimal whole-cell model designed so that inference works from the start. 

Instead of orchestrating isolated black-box simulators, InferCell compiles hybrid dynamics (ODEs/SDEs + stochastic simulation + discrete events) into a single computational graph in pure Julia. AD, likelihood evaluation, and Bayesian calibration fall out of this design. The `infer()` API selects the appropriate backend automatically -- NUTS for differentiable models, ABC-SMC for simulation-based models -- so the user writes the same code regardless of the underlying formalism. The SciML and Turing.jl ecosystems provide the numerical and probabilistic foundations; InferCell composes them into one inference-ready whole-cell model.

## Key design ideas

- **Information flow as design test.** For every module boundary, ask: how does inference cross this interface? If you can't answer that, the design needs to change.
- **Interfaces over implementations.** Sub-models are defined by contracts (state variables, time-stepping, event handling), not internals. A metabolism module can be an ODE solver, a constraint-based method, or a neural surrogate, swappable without touching the rest of the system.
- **Single-language computational graph.** Julia removes the language boundaries that make existing WCM stacks opaque. Solvers, models, and inference algorithms share one compiled representation with native interop to DifferentialEquations.jl, Catalyst.jl, and Turing.jl.

## What success looks like

A working hybrid continuous-discrete model of a minimal biological subsystem that runs forward simulations and does Bayesian parameter inference end-to-end in pure Julia. Gradients flow where components are smooth, likelihoods where available, simulation-based methods everywhere else. Validated by calibration against synthetic or experimental data.

## References

### Whole-cell modelling

- Karr et al. (2012), [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://doi.org/10.1016/j.cell.2012.05.044). The first complete whole-cell model (M. genitalium). Established the sub-model composition paradigm that subsequent WCMs follow.
- Macklin et al. (2020), [Simultaneous cross-evaluation of heterogeneous E. coli datasets via mechanistic simulation](https://doi.org/10.1126/science.aav3751). E. coli whole-cell model. Demonstrates the scale and complexity of current WCM efforts.
- Goldberg et al. (2018), [Emerging Whole-Cell Modeling Principles and Methods](https://doi.org/10.1016/j.copbio.2017.12.013). Review of WCM methodology and open challenges.

### Simulation platforms

- Agmon et al. (2022), [Vivarium: an interface and engine for integrative multiscale modeling in computational biology](https://doi.org/10.1093/bioinformatics/btac049) ([GitHub](https://github.com/vivarium-collective/vivarium-core)). Composable multi-scale simulation framework. Demonstrates the orchestration approach InferCell departs from.
- Thornburg et al. (2022), [Fundamental behaviors emerge from simulations of a living minimal cell](https://doi.org/10.1016/j.cell.2021.12.025). Lattice Microbes applied to JCVI-syn3A. Spatial stochastic simulation at whole-cell scale.
- Luthey-Schulten et al. (2026), [Bringing the genetically minimal cell to life on a computer in 4D](https://www.cell.com/cell/fulltext/S0092-8674(26)00174-1). Full cell-cycle simulation of JCVI-syn3A at nanoscale resolution integrating metabolism, genetic information processing, and morphological dynamics. Current state of the art in simulation-first whole-cell modelling.

