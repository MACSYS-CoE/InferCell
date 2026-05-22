# Getting started

## Installation

InferCell is a Julia package targeting Julia ≥ 1.10. Clone the repository and instantiate the project:

```bash
git clone https://github.com/MACSYS-CoE/InferCell.git
cd InferCell
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Quick start

The same four steps work for both differentiable (ODE) and simulation-based (SSA) models — `infer()` dispatches to the right backend based on the sub-model's declared formalism.

### ODE model → NUTS (gradient-based)

```julia
using InferCell

model = TranscriptionTranslation()
prob  = build_problem(model; tspan=(0.0, 50.0))
sol   = solve(prob, Tsit5(); saveat=0:2.5:50)
data  = observe(sol, 0:2.5:50, model; sigma=0.3)
chain = infer(model, data; n_samples=1000)
```

### SSA model → ABC-SMC (likelihood-free)

```julia
using InferCell

model  = StochasticGeneExpression()
prob   = build_problem(model; tspan=(0.0, 50.0))
trajs  = [solve(prob, SSAStepper(); saveat=0:2.5:50) for _ in 1:200]
data   = observe(trajs, 0:2.5:50, model)
result = infer(model, data; n_samples=300, n_populations=5,
               tspan=(0.0, 50.0))
```

Same four steps, same `infer()` call. The framework dispatches to the right backend based on the model's formalism — there's no separate "ABC API" or "NUTS API" exposed to the user.

## Running tests

```bash
julia --project -e 'using Pkg; Pkg.test()'

# Include integration tests (~2 min, NUTS sampling):
INFERCELL_INTEGRATION_TESTS=true julia --project -e 'using Pkg; Pkg.test()'
```

The full test suite — including the iterative-boundary integration test — is typically run on a Slurm cluster via `sbatch test/run_tests.slurm` (~2.5 min) and `sbatch test/run_integ_iterative.slurm` (~30 min).

## Next steps

- [User guide → ODE inference](user-guide/ode-inference.md) walks through the NUTS path in detail.
- [User guide → SSA inference](user-guide/ssa-inference.md) covers the ABC-SMC path.
- [User guide → Boundary protocol](user-guide/boundary-protocol.md) explains how posteriors cross block boundaries.
- [API reference](api/index.md) lists every exported symbol.
