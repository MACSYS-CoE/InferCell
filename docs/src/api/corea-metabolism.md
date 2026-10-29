```@meta
CurrentModule = InferCell
```

# Metabolism

The reduced syn3A model's four ODE modules: glycolysis from glucose-6-phosphate
to lactate, the phosphotransferase cascade with passive lactate export,
nucleotide recycling, and the lumped tRNA charging step. The
[overview](corea-modules.md) lists the other pages.

```@autodocs
Modules = [InferCell]
Pages = [
    "organisms/coreA/central_glycolysis.jl",
    "organisms/coreA/pts_transport.jl",
    "organisms/coreA/nucleotide_recycling.jl",
    "organisms/coreA/trna_charging.jl",
]
```
