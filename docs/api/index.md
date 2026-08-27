# API reference

InferCell exports its public API across a dozen source files. The pages in this section are organised by source module; each page lists the public API with the same docstring text returned by `?<symbol>` in the Julia REPL.

| Module | Page | Highlights |
|---|---|---|
| `src/parameters.jl` | [Parameters](parameters.md) | `InferParameter`, role-filtering helpers |
| `src/interface.jl` | [Sub-model interface](interface.md) | `AbstractSubModel`, protocol functions |
| `src/corea/*.jl` | [Core A′ interface contract](corea-interface.md) | `COREA_SPECIES`, the seven `CouplingEdge` kinds, `resolve_coupling`, `load_parameter` |
| `src/orchestrator.jl` | [Orchestrator](orchestrator.md) | `build_problem` |
| `src/inference.jl` | [Inference](inference.md) | `infer`, `observe`, `posterior_predictive`, `check_identifiability` |
| `src/boundary.jl` | [Boundary](boundary.md) | `sequential_infer`, `iterative_infer`, `KDEPrior` |
| `src/models/*.jl` | [Models](models.md) | `TranscriptionTranslation`, `BurstyGeneExpression`, `LightMetabolism`, … |

For the canonical signatures, query the live module:

```julia
using InferCell
?TranscriptionTranslation     # struct docstring
?infer                        # function docstring
```

This site mirrors the docstrings as of the most recent release; the REPL is authoritative.
