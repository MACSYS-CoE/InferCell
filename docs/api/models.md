# Models

Built-in sub-models. Each is a concrete subtype of [`AbstractSubModel`](interface.md). Source files under [`src/models/`](https://github.com/MACSYS-CoE/InferCell/tree/main/src/models).

## `TranscriptionTranslation`

```julia
TranscriptionTranslation(; k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1,
                          mRNA0=0.0, protein0=0.0, sigma_obs=0.3)
TranscriptionTranslation(genes::Vector{Symbol}; …, enzyme_slot=auto, overrides=Dict())
```

Mean-field ODE block. The no-arg form is the single-gene Phase-1 model with states `[:mRNA, :protein]`. The `genes`-arg form is the multi-gene v0.0.1 variant: per-gene state names `[:<gene>_mRNA, :<gene>_protein]` and per-gene rate names `[:k_tx_<gene>, …]`. `enzyme_slot` (default `:enzyme` if present in `genes`) marks which gene's protein feeds metabolism. `overrides::Dict{Symbol, NamedTuple}` lets specific genes get bespoke rate values.

Formalism: `:ode`. Inference: NUTS.

## `StochasticGeneExpression`

```julia
StochasticGeneExpression(; k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1,
                          mRNA0=0, protein0=0)
```

Constitutive stochastic transcription + translation as a jump process: four reactions (mRNA production at constant rate `k_tx`, mRNA degradation, translation, protein degradation). States `[:mRNA, :protein]`. No ODE equivalent below the law-of-large-numbers limit.

Formalism: `:jump`. Inference: ABC-SMC.

## `BurstyGeneExpression`

```julia
BurstyGeneExpression(; k_on=0.2, k_off=0.5, k_tx_burst=10.0, k_tl=2.0,
                      gamma_mRNA=0.5, gamma_protein=0.1,
                      promoter0=0, mRNA0=0, protein0=0)
```

Telegraph-promoter stochastic gene expression. A two-state promoter (`:promoter ∈ {0, 1}`) gates mRNA production, producing visible mRNA-count bursting (Fano factor > 1). States `[:promoter, :mRNA, :protein]`; six jump reactions. The v0.0.1 SSA block — load-bearing for the boundary protocol because it has no ODE equivalent.

Formalism: `:jump`. Inference: ABC-SMC.

## `LightMetabolism`

```julia
LightMetabolism(; k_atp=1.0, k_ntp=0.2, k_aa=0.4, k_tx=1.0, k_tl=2.0,
                  ATP_max=10.0, e_tx=0.1, e_tl=0.05, n_tx=0.5, n_tl=0.5,
                  gamma_ntp=0.1, gamma_aa=0.1, K_M_enzyme=5.0,
                  ATP0=5.0, NTP0=2.0, AA0=2.0, sigma_obs=0.3,
                  mRNA_source=:mRNA, enzyme_source=:protein)
```

Lightweight ODE metabolism: ATP production gated by enzyme via Michaelis-Menten (constant `K_M_enzyme`), NTP and AA synthesis, TX/TL-driven sinks. States `[:ATP, :NTP, :AA]`. Couples to a gene-expression block via `mRNA_source` (drives TX consumption) and `enzyme_source` (modulates ATP production). Closes the v0.0.1 autocatalytic loop with [`TranscriptionTranslation`](#transcriptiontranslation).

Formalism: `:ode`. Inference: NUTS.

## `TierBMetabolism`

```julia
TierBMetabolism(; …, k_atp_leak=0.2, n_hill=2.0,
                  mRNA_source=:mRNA, enzyme_source=:protein)
```

"v0.0.1 + ε" metabolism for Tier-B mismatch experiments. Identical interface to [`LightMetabolism`](#lightmetabolism) but with a Hill enzyme response (`n_hill > 1`) and a constant ATP leak `k_atp_leak` — two sources of model error the v0.0.1 fitter must absorb. Not an inference target itself; only used to generate synthetic data.

Formalism: `:ode`. Inference: NUTS (for the mismatch experiment generator).
