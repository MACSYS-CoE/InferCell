# API reference

The pages in this section are organised by source module. Each splices the live docstrings in beside the prose that explains why the design is what it is; the Core A′ page is mostly the latter, because a contract is not a signature.

| Module | Page | Highlights |
|---|---|---|
| `src/parameters.jl` | [Parameters](parameters.md) | `InferParameter`, role-filtering helpers |
| `src/interface.jl` | [Sub-model interface](interface.md) | `AbstractSubModel`, protocol functions |
| `src/edges.jl`, `src/resolver.jl`, `src/loader.jl`, `src/labels.jl`, `src/organisms/coreA/registry.jl` | [Core A′ interface contract](corea-interface.md) | `COREA_SPECIES`, the seven `CouplingEdge` kinds, `resolve_coupling`, `load_parameter` |
| `src/orchestrator.jl` | [Orchestrator](orchestrator.md) | `build_problem` |
| `src/handshake.jl` | [Handshake driver](handshake.md) | `HandshakeDriver`, `run_handshake!`, `clipping_census`, `rebuild_census`, `chemostat_census`, `rate_constant_elasticity`, `driver_declarations`, the count↔concentration conversion |
| `src/inference.jl` | [Inference](inference.md) | `infer`, `observe`, `posterior_predictive`, `check_identifiability` |
| `src/boundary.jl` | [Boundary](boundary.md) | `sequential_infer`, `iterative_infer`, `KDEPrior` |
| `src/models/*.jl` | [Models](models.md) | `TranscriptionTranslation`, `BurstyGeneExpression`, `LightMetabolism`, … |

Every signature and field list on these pages is generated from the docstring
attached to the symbol itself, so `?infer` in the REPL and the [Inference](inference.md)
page cannot disagree. The build fails if an exported symbol has no docstring, or if a
cross-reference names a symbol that no longer exists.
