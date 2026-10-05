```@meta
CurrentModule = InferCell
```

# Core A′ observables and synthetic data

What the likelihood may condition on, and the synthetic data it is fitted to
(spec §11 phase 15). This page covers:
- the circularity guard and the truth-drawing rule of spec §4 D8;
- the 60 s observables, fluxes and replicate ensembles;
- the ensemble sensitivity and its noise-floor rank test;
- output F2 on check 8's frozen-expression variant;
- the dataset generator;
- M0, the two-gene reference system of phase 16 (task 16a.8), and its datasets.

```@autodocs
Modules = [InferCell]
Pages = [
    "organisms/coreA/observables.jl",
    "organisms/coreA/control_variant.jl",
    "organisms/coreA/m0.jl",
]
```
