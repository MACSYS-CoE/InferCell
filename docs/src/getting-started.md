# Getting started

## Installation

InferCell is a Julia package targeting Julia 1.10. Clone the repository and
instantiate the project:

```bash
git clone https://github.com/MACSYS-CoE/InferCell.git
cd InferCell
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Run the reduced syn3A model

[`build_corea`](@ref) composes all seven modules and checks the composition is
complete: every state has an owner and no moiety is stranded. It returns a
[`HandshakeDriver`](@ref), which [`run_handshake!`](@ref) steps one simulated
second at a time.

```julia
using InferCell, Random

Random.seed!(1)
d = build_corea(; tspan = (0.0, 180.0))
rec = run_handshake!(d, 180)
rec.t[end], size(rec.jump)
```

`rec` holds both blocks' states and parameters at every handshake, with the
cell's geometry, the clipping census and the rate-constant rebuilds. Seeding the
global generator before the build fixes the cell.

## Fit a synthetic dataset

The sampler is validated on M0, a two-gene reference system cut from the same
model. [`m0_dataset`](@ref) draws a truth from the prior and simulates a
population of cells under it. [`gibbs_setup`](@ref) holds what a chain keeps
fixed, and [`gibbs_sweep!`](@ref) is one sweep of the three-block sampler.

```julia
ds = m0_dataset(16085; ncells = 3, horizon = 180)
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

The chain's free parameters are the two promoter strengths and the transcript
decay constant (`g.cme`), the log enolase and aldolase forward constants
(`g.forwards`), and the metabolite noise scale σ. This example is small enough to
run in minutes. A production chain uses 50 cells over 600 s with `N = 20`
particles, and one sweep takes about 13 minutes on five threads. See
[the conditional sampler](api/conditional-sampler.md).

## Running the tests

```bash
julia --project -e 'using Pkg; Pkg.test()'
```

The suite takes about an hour on four cores. On a Slurm cluster, run
`sbatch test/run_tests.slurm`.

## Next steps

- [The reduced syn3A model](user-guide/reduced-syn3a.md): what it contains and
  what it leaves out.
- [The conditional sampler](api/conditional-sampler.md): how the three
  blocks fit together.
- [API reference](api/index.md): every exported symbol.
