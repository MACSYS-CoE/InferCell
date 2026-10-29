```@meta
CurrentModule = InferCell
```

# Model modules

The seven sub-models that make up the reduced syn3A organism. Each declares its
own coupling edges against the [interface contract](corea-interface.md)
and vendors the parameter extract it imports from, so a module can be read,
tested and composed on its own.

Each module's vendored extracts live under `src/organisms/coreA/data/`, whose
`README.md` records the upstream file, the commit it was taken at and the
command that regenerates it byte for byte.

The docstrings are split over three pages, by block:

| Page | Sources |
|---|---|
| [Metabolism](corea-metabolism.md) | The four ODE modules: central glycolysis, sugar transport and lactate export, nucleotide recycling, and the lumped tRNA charging step |
| [Gene expression](corea-gene-expression.md) | The three jump modules (transcription, transcript decay, translation), their mean field, and route (ii) |
| [Assembly, priors and checks](corea-assembly.md) | `build_corea`, the Haldane-consistent prior draws, the balance checks and the control analysis |
