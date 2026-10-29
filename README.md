# InferCell.jl

> Composable, inference-first whole-cell modelling in Julia.

[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://macsys-coe.github.io/InferCell/stable/)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://macsys-coe.github.io/InferCell/dev/)

Whole-cell models mix formalisms. Metabolism is a system of ordinary
differential equations, gene expression is a stochastic jump process, and the
two exchange state on a clock. InferCell asks whether parameter uncertainty can
be propagated *correctly* across that seam: when the model is fitted to data, do
the posteriors mean what they claim to mean?

InferCell is a compositional inference framework, demonstrated on a **reduced
JCVI-syn3A model**. That model is a small organism cut from the published syn3A
whole-cell model of
[Thornburg et al. (2022)](https://doi.org/10.1016/j.cell.2021.12.025), and
built from modules that touch only through a typed coupling boundary. It is
fitted two ways: exactly, by a sampler that treats the stochastic path as a
latent variable, and approximately, with gene expression replaced by its mean
field.

## What is here

- **The reduced syn3A model.** Twenty-two reactions in four ODE modules:
  glycolysis to lactate, sugar transport and lactate export, nucleotide
  recycling, and a lumped tRNA charging step. Twenty-one of them are published
  reactions; the charging step stands in for a hundred. Seventeen genes in three
  jump modules: transcription, translation and transcript decay. The two blocks
  swap state every simulated second and rate constants every sixty. About 92% of
  the cell is removed, leaving its best-measured region.
- **A composition framework.** Each module declares the states it owns and its
  coupling edges, and the resolver checks the composition before anything runs:
  every state has an owner, every edge a partner, and no moiety is stranded.
  Parameters carry their upstream provenance, and every place the reduction
  departs from the published model is enumerated rather than left in prose.
- **The hybrid driver.** It steps the ODE and jump blocks in lock-step, with
  conservation, clipping and growth reported from the outside, and replays a
  recorded jump path exactly.
- **The exact route: a conditional sampler.** A three-block Gibbs chain over the
  ODE parameters, the stochastic rate constants and the jump path. The path is
  updated by conditional SMC, with exact transcript bridges between
  observations. Each block is unit-tested against an exact or independent
  reference.
- **The mixed route: a mean-field description.** `MeanField` swaps the jump
  modules for their mean equations, and route (ii), `build_mean_field_problem`,
  runs the whole model as one ODE that executes the coupling edges. Transcripts
  are scored by their exact likelihood at the mean-field rate constants, and
  bulk metabolites by a Gaussian whose offset and covariance come from a hybrid
  ensemble at a reference point. The posterior is sampled by NUTS.
  `infer(models, ds; route = :exact)` or `route = :mixed` runs either route on
  an M0 dataset.
- **Synthetic data.** Observables, noise models and a dataset generator that
  records its truth, seeds and reductions. Truths are drawn from the prior, and
  a calibration run refuses the nominal parameter values.

## Status

Pre-release research code. The status is stated plainly so that no release
claims more than it has.

- **Done and tested:** the model, the framework, the hybrid driver and both
  inference routes. The suite has about 29,000 assertions.
- **Calibrated on M0,** the two-gene reference system, with six parameters:
  - **The exact route** passes simulation-based calibration (SBC) on all 6
    parameters over 50 prior-drawn datasets (V7, merge job 18689602;
    `dev/scripts/v7_16a9b_result.md`).
  - **The mixed route,** in the wording of `spec/phases/approximate-method.md`
    §12 (2026-10-29): On M0 the mixed route passes SBC on all 6 parameters over
    V7's 50 datasets (job 19292232), at a median 6.57 CPU-h per posterior
    against 507.5–678 for the exact sampler. Its largest mean bias against the
    exact posteriors (V7) is 0.095 ± 0.050 SD (krnadeg). On the single dataset
    m0_v6, S_0779, krnadeg and ENO sit 0.42–0.56 SD from the exact reference
    (V6), about 0.3 SD above the matched noise floor, and the cause is open.
- **In progress, for v0.3.0:** `predict`, one entry point for forward and
  intervention predictions. When it lands, v0.3.0 will make the claim in
  [`spec/spec.md`](spec/spec.md) §1:

  > InferCell supports exact stochastic and mean-field descriptions across a
  > coupled ODE–jump boundary. On M0, both inference routes pass
  > simulation-based calibration, with small declared discrepancies.
  > Posterior-updated parameters can be combined modularly with
  > viability-conditioned uncertainty in the remaining parameters to produce
  > forward and intervention predictions on the reduced model.

  The first two sentences hold on `main` now; the third waits on `predict`.
- **Not claimed:** an all-parameter posterior, calibrated inference on the full
  reduced model, or a biological finding. The data are synthetic and the
  priors asserted. Predictions on the reduced model are unscored
  demonstrations of the method; only M0 is scored.
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

The dataset is from M0, the two-gene reference system on which both routes are
calibrated. Five of its six parameters are the GAPD and ptsG promoter
strengths, the transcript decay constant, and the enolase and aldolase
catalytic constants; the sixth is the bulk-assay noise σ_b. The code's `CoreA`
and `corea` prefixes name the reduced syn3A model.

The [documentation](https://macsys-coe.github.io/InferCell/stable/) covers the
model, the coupling contract, the driver and both routes, with every exported
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
- Thornburg et al. (2026), [Bringing the genetically minimal cell to life on a computer in 4D](https://doi.org/10.1016/j.cell.2026.02.009), Cell 189(9). Full cell-cycle simulation of JCVI-syn3A, the current state of the art in simulation-first whole-cell modelling.
- Karr et al. (2012), [A Whole-Cell Computational Model Predicts Phenotype from Genotype](https://doi.org/10.1016/j.cell.2012.05.044). The first complete whole-cell model.
- Macklin et al. (2020), [Simultaneous cross-evaluation of heterogeneous E. coli datasets via mechanistic simulation](https://doi.org/10.1126/science.aav3751). The E. coli whole-cell model.
- Agmon et al. (2022), [Vivarium: an interface and engine for integrative multiscale modeling in computational biology](https://doi.org/10.1093/bioinformatics/btac049). The orchestration approach InferCell departs from.
- Andrieu, Doucet and Holenstein (2010), [Particle Markov chain Monte Carlo methods](https://doi.org/10.1111/j.1467-9868.2009.00736.x). Particle Gibbs, which the path update is.
- Hobolth and Stone (2009), [Simulation from endpoint-conditioned, continuous-time Markov chains on a finite state space](https://doi.org/10.1214/09-AOAS247). The transcript bridges.
- Talts et al. (2018), [Validating Bayesian inference algorithms with simulation-based calibration](https://arxiv.org/abs/1804.06788). The calibration check.
