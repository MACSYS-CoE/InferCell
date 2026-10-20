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

**Mean-field gene expression** (track B, task B2; `spec/phases/approximate-method.md`). `MeanField(members)` turns a group of first-order `:jump` modules into one `:ode` module carrying their mean equations, `d E[x]/dt = Σ_r a_r(E[x]) ν_r`. It reads each propensity from the member's own closure and gets each state change by applying the member's own `affect!`, so the two descriptions cannot drift apart. It refuses a propensity that is not affine in the state, since the means are exact only then. `gene_expression_mean_field(; genes)` groups transcription, decay and translation, with the same states, parameters and coupling edges as the jump modules. The ODE path does not yet execute those edges (task B4). With `edges = :none` it builds alone at pinned rate constants.

```@autodocs
Modules = [InferCell]
Pages = [
    "organisms/coreA/central_glycolysis.jl",
    "organisms/coreA/pts_transport.jl",
    "organisms/coreA/nucleotide_recycling.jl",
    "organisms/coreA/trna_charging.jl",
    "organisms/coreA/transcription.jl",
    "organisms/coreA/transcript_decay.jl",
    "organisms/coreA/translation.jl",
    "organisms/coreA/assembly.jl",
    "organisms/coreA/prior_draws.jl",
    "organisms/coreA/validation.jl",
    "control_analysis.jl",
    "mean_field.jl",
]
```
