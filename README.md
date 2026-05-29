# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

## Quick start (zero to running)

New to Julia? This walks you from nothing installed to running your first
inference. It takes about 10–15 minutes, most of which is Julia downloading
and precompiling dependencies the first time.

### 1. Install Julia

The recommended way is [`juliaup`](https://github.com/JuliaLang/juliaup), the
official version manager — it installs Julia and keeps it up to date.

```bash
# macOS / Linux
curl -fsSL https://install.julialang.org | sh

# Windows (PowerShell)
winget install julia -s msstore
```

Restart your shell, then confirm it works (InferCell needs Julia ≥ 1.10):

```bash
julia --version
```

### 2. Get the code

```bash
git clone https://github.com/MACSYS-CoE/InferCell.git
cd InferCell
```

### 3. Install the dependencies

`Manifest.toml` is committed, so this reproduces the *exact* package versions
CI uses — no version guessing. Run it once after cloning:

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

The first run downloads and precompiles everything (Turing, the SciML stack,
etc.), so expect a few minutes. Subsequent starts are fast.

> **What does `--project` mean?** It tells Julia to use this repository's
> environment (its `Project.toml`/`Manifest.toml`) instead of your global one,
> so InferCell's dependencies stay isolated. Run every command below from the
> `InferCell` directory so `--project` picks it up.

### 4. Run your first inference

Save this as `first_run.jl` in the repo root:

```julia
using InferCell

# 1. Pick a model        2. Build the simulation problem
model = TranscriptionTranslation()
prob  = build_problem(model; tspan=(0.0, 50.0))

# 3. Simulate, then fake some noisy "experimental" data from the solution
sol  = solve(prob, Tsit5(); saveat=0:2.5:50)
data = observe(sol, 0:2.5:50, model; sigma=0.3)

# 4. Infer the parameters back from the data (NUTS, gradient-based)
chain = infer(model, data; n_samples=1000)

# Look at the posterior summary
display(chain)
```

Run it:

```bash
julia --project first_run.jl
```

You'll see a posterior summary table for the model's parameters — that's a
full simulate → observe → infer loop. Prefer working interactively? Start a
REPL with `julia --project`, then paste the lines in one at a time (results
stay in memory between lines, which is handy for poking at `chain`).

### 5. Verify your setup with the test suite

```bash
# Unit tests only (~1 min)
julia --project -e 'using Pkg; Pkg.test()'

# Add the slower Bayesian integration tests (~2 min, runs NUTS sampling)
INFERCELL_INTEGRATION_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
```

If these pass, you're fully set up. 🎉

## The InferCell workflow at a glance

The same four steps — pick a model, build a problem, generate data, call
`infer()` — work for both differentiable (ODE) and simulation-based (SSA)
models. `infer()` dispatches to the right backend based on the model's
declared formalism:

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

Same four steps, same `infer()` call — no separate "ABC API" or "NUTS API"
exposed to the user.

## Where to go next

- [`docs/getting-started.md`](docs/getting-started.md) — the same setup with
  links into the detailed user guide.
- [`docs/plans/`](docs/plans/) — current status, roadmap, and design decisions.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — local workflow, CI checks, and branch
  conventions if you want to contribute.

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

## Capabilities

Inference is the design discipline; the following fall out as named framework outputs, not separate features that need building. See [`docs/positioning.md`](docs/positioning.md) for the longer argument.

- **Parameter inference.** Bayesian posteriors over rates, pool sizes, and initial conditions from time-series data. `infer()` dispatches to NUTS for differentiable models and ABC-SMC for stochastic ones behind one call. Posteriors propagate across heterogeneous block boundaries via explicit boundary protocols.
- **Sensitivity analysis.** Local and global parameter sensitivities. The AD machinery powering NUTS gradients powers forward and adjoint sensitivities for free — useful for asking which parameters matter for which observables before committing to data collection.
- **Model reduction.** Identifiability checks and sloppy-parameter analysis fall out of the same gradient and posterior machinery. The framework can tell you which parameters your data can't constrain, and where the model has more degrees of freedom than the biology warrants.
- **Model comparison.** Bayes factors and posterior odds across competing model structures (constitutive vs bursty transcription, with vs without enzyme feedback, etc.). Each candidate is one `infer()` call; the comparison is downstream.

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

