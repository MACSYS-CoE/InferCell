```@meta
CurrentModule = InferCell
```

# Assembly, priors and checks

How the seven modules are put together into the reduced syn3A model
(`corea_models`, `build_corea`), prior draws that keep each reaction's
equilibrium constant, the balance checks on the assembled model, and the
steady-state control analysis. The [overview](corea-modules.md) lists the other
pages.

```@autodocs
Modules = [InferCell]
Pages = [
    "organisms/coreA/assembly.jl",
    "organisms/coreA/prior_draws.jl",
    "organisms/coreA/validation.jl",
    "control_analysis.jl",
]
```
