```@meta
CurrentModule = InferCell
```

# Gene expression

The reduced syn3A model's three jump modules (transcription, transcript decay and
translation), their mean field, and route (ii), which runs the whole model with
gene expression swapped for that mean field. The [overview](corea-modules.md)
lists the other pages.

**Mean-field gene expression** (track B, tasks B2 and B4; `spec/phases/approximate-method.md`). `MeanField(members)` turns a group of first-order `:jump` modules into one `:ode` module carrying their mean equations, `d E[x]/dt = Σ_r a_r(E[x]) ν_r`. It reads each propensity from the member's own closure and gets each state change by applying the member's own `affect!`, so the two descriptions cannot drift apart. It refuses a propensity that is not affine in the state, since the means are exact only then. `gene_expression_mean_field(; genes)` groups transcription, decay and translation, with the same states, parameters and coupling edges as the jump modules. With `edges = :none` it builds alone at pinned rate constants. Route (ii), `build_mean_field_problem`, runs the reduced model as one ODE with the jump modules swapped for one `MeanField`, and executes their coupling edges as part of the right-hand side or as callbacks; the plain `build_problem` does not.

```@autodocs
Modules = [InferCell]
Pages = [
    "organisms/coreA/transcription.jl",
    "organisms/coreA/transcript_decay.jl",
    "organisms/coreA/translation.jl",
    "mean_field.jl",
    "mean_field_route.jl",
]
```
