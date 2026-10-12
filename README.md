# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://macsys-coe.github.io/InferCell/stable/)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://macsys-coe.github.io/InferCell/dev/)

Whole-cell models mix formalisms. Metabolism is a system of ordinary
differential equations, gene expression is a stochastic jump process, and the
two exchange state on a clock. InferCell asks whether parameter uncertainty can
be propagated *correctly* across that seam: when the model is fitted to data, do
the posteriors mean what they claim to mean?

It answers on a **reduced JCVI-syn3A model**: a small organism cut from the
published syn3A whole-cell model of
[Thornburg et al. (2022)](https://doi.org/10.1016/j.cell.2021.12.025), built
from modules that touch only through a typed coupling boundary, and fitted by a
sampler that treats the stochastic path as a latent variable.

## What is in this release

- **The reduced syn3A model.** Twenty-two metabolic reactions in four ODE
  modules: glycolysis to lactate, sugar transport and lactate export, nucleotide
  recycling, and a lumped tRNA charging step. Seventeen genes in three jump
  modules: transcription, translation and transcript decay. The two blocks swap
  state every simulated second and rate constants every sixty. About 92% of the
  cell is removed, leaving its best-measured region.
- **A composition framework.** Each module declares the states it owns and its
  coupling edges, and the resolver checks the composition before anything runs:
  every state has an owner, every edge a partner, and no moiety is stranded.
  Parameters carry their upstream provenance, and every place the reduction
  departs from the published model is enumerated rather than left in prose.
- **The hybrid driver.** It steps the ODE and jump blocks in lock-step, with
  conservation, clipping and growth reported from the outside, and replays a
  recorded jump path exactly.
- **A conditional sampler.** A three-block Gibbs chain over the ODE parameters,
  the stochastic rate constants and the jump path. The path is updated by
  conditional SMC, with exact transcript bridges between observations. Each block
  is unit-tested against an exact or independent reference.
- **Synthetic data.** Observables, noise models and a dataset generator that
  records its truth, seeds and reductions. Truths are drawn from the prior, and
  a calibration run refuses the nominal parameter values.

## Status

Pre-release research code, and the status is stated plainly so the release
claims no more than it has.

- **Done and tested:** the model, the framework, the hybrid driver and each
  sampler block. The suite has about 28,000 assertions.
- **Running:** the long reference chain and the simulation-based calibration of
  the sampler on a two-gene reference system. **Calibration has not yet passed.**
  This release makes no claim that the posteriors are calibrated.
- **Not started:** recovery and coverage on the full reduced model, which waits
  on a measured cost per posterior.
- **Not in scope:** real data. Two published measurements serve as sanity
  checks only. The reduced model does not reproduce syn3A's doubling time,
  metabolite levels or proteome, and is not to be compared with them. The code
  refuses a doubling time outright.

## Quick start

InferCell targets Julia 1.10. Clone it and instantiate the project:

```bash
git clone https://github.com/MACSYS-CoE/InferCell.git
cd InferCell
julia --project -e 'using Pkg; Pkg.instantiate()'
```

Run the reduced syn3A model forward, then fit a small synthetic dataset:

```julia
using InferCell, Random

Random.seed!(1)
d = build_corea(; tspan = (0.0, 180.0))   # all seven modules, in completeness mode
rec = run_handshake!(d, 180)               # 180 one-second exchanges
rec.t[end], size(rec.jump)

ds = m0_dataset(16085; ncells = 3, horizon = 180)   # truth drawn from the prior
g = gibbs_setup(ds; N = 4, horizon = 180)

rng = Xoshiro(2)
start = Dict(draw_truth(rng, g.models; names = M0_TARGETS, purpose = :calibration).values)
s = gibbs_init(rng, g; cme = [start[n] for n in g.cme],
               u = [log(start[f]) for f in g.forwards], σ = rand(rng, SIGMA_MET_PRIOR))
for _ in 1:5
    gibbs_sweep!(s, g)
end
s.trace[end]
```

The dataset is from the two-gene reference system, where the sampler is
validated. Its five free parameters are the GAPD and ptsG promoter strengths,
the transcript decay constant, and the enolase and aldolase catalytic constants.
The code's `CoreA` and `corea` prefixes name the reduced syn3A model.

The [documentation](https://macsys-coe.github.io/InferCell/stable/) covers the
model, the coupling contract, the driver and the sampler, with every exported
symbol's docstring.

## Running the tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

The suite takes about an hour on four cores. On a Slurm cluster:
`sbatch test/run_tests.slurm`.

## Earlier toy models

The package keeps the toy models it started from. These are a transcription and
translation ODE, stochastic and bursty gene expression, and two small metabolic
models, with a single `infer()` call. That call dispatches to NUTS for an ODE
model and to ABC-SMC for a jump model. They are documented under *Toy models*.
The posterior-passing protocols between those two backends, `sequential_infer`
and `iterative_infer`, are proofs of concept, not calibrated methods.

## Where things are

- [`spec/spec.md`](spec/spec.md) is the research specification: the claim, the
  validity checks, the decisions of record and the task list.
- [`dev/notes/`](dev/notes/) holds the research notes the spec rests on.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) covers the local workflow and branch
  conventions.

## Licence

InferCell's own source, tests and documentation are released under the MIT
licence — see [`LICENSE`](LICENSE).

**That covers this project's code, not the parameter values it imports.** The
extracts under `src/organisms/coreA/data/` are derived from
[`Luthey-Schulten-Lab/Minimal_Cell`](https://github.com/Luthey-Schulten-Lab/Minimal_Cell),
which is public but carries no licence of its own and therefore grants no
redistribution rights. They are reshaped tables of published scientific
parameter values, retained here so the model is reproducible without network
access; `src/organisms/coreA/data/README.md` records the upstream file, commit
and regeneration command for each. If you intend to redistribute them, or to
build on them commercially, take that up with the upstream authors rather than
relying on the licence above.

## References

- Thornburg et al. (2022), [Fundamental behaviors emerge from simulations of a living minimal cell](https://doi.org/10.1016/j.cell.2021.12.025). The syn3A whole-cell model this one is reduced from.
- Luthey-Schulten et al. (2026), [Bringing the genetically minimal cell to life on a computer in 4D](https://www.cell.com/cell/fulltext/S0092-8674(26)00174-1). Full cell-cycle simulation of JCVI-syn3A, the current state of the art in simulation-first whole-cell modelling.
- Karr et al. (2012), [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://doi.org/10.1016/j.cell.2012.05.044). The first complete whole-cell model.
- Macklin et al. (2020), [Simultaneous cross-evaluation of heterogeneous E. coli datasets via mechanistic simulation](https://doi.org/10.1126/science.aav3751). The E. coli whole-cell model.
- Agmon et al. (2022), [Vivarium: an interface and engine for integrative multiscale modeling in computational biology](https://doi.org/10.1093/bioinformatics/btac049). The orchestration approach InferCell departs from.
- Andrieu, Doucet and Holenstein (2010), [Particle Markov chain Monte Carlo methods](https://doi.org/10.1111/j.1467-9868.2009.00736.x). Particle Gibbs, which the path update is.
- Hobolth and Stone (2009), [Simulation from endpoint-conditioned, continuous-time Markov chains on a finite state space](https://doi.org/10.1214/09-AOAS247). The transcript bridges.
- Talts et al. (2018), [Validating Bayesian inference algorithms with simulation-based calibration](https://arxiv.org/abs/1804.06788). The calibration check.
