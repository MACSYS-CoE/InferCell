# API reference

InferCell exports its public API across twenty source files. The pages in this section are organised by source module; most pages list the public API with the same docstring text returned by `?<symbol>` in the Julia REPL, while the Core A′ page documents its contract as prose.

| Module | Page | Highlights |
|---|---|---|
| `src/parameters.jl` | [Parameters](parameters.md) | `InferParameter`, role-filtering helpers |
| `src/interface.jl` | [Sub-model interface](interface.md) | `AbstractSubModel`, protocol functions |
| `src/edges.jl`, `src/resolver.jl`, `src/loader.jl`, `src/labels.jl`, `src/organisms/coreA/registry.jl` | [Core A′ interface contract](corea-interface.md) | `COREA_SPECIES`, the seven `CouplingEdge` kinds, `resolve_coupling`, `load_parameter` |
| `src/orchestrator.jl` | [Orchestrator](orchestrator.md) | `build_problem` |
| `src/handshake.jl` | [Handshake driver](handshake.md) | `HandshakeDriver`, `run_handshake!`, `clipping_census`, `rebuild_census`, `rate_constant_elasticity`, `driver_declarations`, the count↔concentration conversion |
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
