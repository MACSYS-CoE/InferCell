# Phase 2 Handoff

## Current state

Phase 2 on branch `phase-2-stochastic` is complete. The demo (`examples/phase2_demo.jl`) runs two independent parameter recovery experiments through the same `infer()` API:

1. **ODE model** (`TranscriptionTranslation`, formalism `:ode`) → dispatches to NUTS. All 4 parameters recovered within 90% CI.
2. **SSA model** (`StochasticGeneExpression`, formalism `:jump`) → dispatches to ABC-SMC. Parameters recovered with wider posteriors (expected for likelihood-free inference).

Figure saved to `examples/figures/phase2_parameter_recovery.png`.

## Known limitation

Both models infer the same 4 parameters (k_tx, k_tl, gamma_mRNA, gamma_protein) for the same TX/TL biology. The SSA model is just a stochastic reformulation of the ODE — there is nothing about it that *requires* likelihood-free inference. This weakens the demo: a reader could reasonably ask "why not just use NUTS for both?"

## Next PR: Bursting gene expression model

Replace `StochasticGeneExpression` in the demo with a model that has **no ODE equivalent**, making the case for ABC-SMC self-evident.

### Proposed model: Transcriptional bursting

A two-state gene model where the promoter switches between OFF and ON states:

```
Gene_OFF  →  Gene_ON     (rate: k_on)
Gene_ON   →  Gene_OFF    (rate: k_off)
Gene_ON   →  Gene_ON + mRNA  (rate: k_tx)
mRNA      →  ∅           (rate: gamma_mRNA)
mRNA      →  mRNA + protein  (rate: k_tl)
protein   →  ∅           (rate: gamma_protein)
```

This produces **transcriptional bursts** — mRNA is produced in discrete pulses when the gene is ON, then silent when OFF. The burst size distribution (geometric, mean = k_tx / k_off) and burst frequency (k_on) are inherently stochastic properties with no smooth ODE limit when gene copy number is 1.

### Parameters to infer

- `k_on`, `k_off` (gene switching rates) — these are the genuinely stochastic parameters
- `k_tx`, `gamma_mRNA` (transcription and degradation)
- Optionally fix `k_tl` and `gamma_protein` to reduce dimensionality

### Implementation steps

1. Create `src/models/bursting_gene.jl` with `formalism = :jump`, `inference_mode = :simulation`
2. Define the 5-6 reactions as `ConstantRateJump`s (3 species: gene_state, mRNA, protein)
3. Update the demo to use `BurstingGene` as the Part 2 model
4. May need to tune ABC-SMC settings (more particles/populations) and summary statistics for the higher-dimensional parameter space
5. Add tests in `test/test_bursting_gene.jl`

### Summary statistics consideration

The current summary statistics (mean trajectories) may not capture burst dynamics well. Consider adding variance-based or distribution-based statistics (e.g. Fano factor, coefficient of variation across replicates at each timepoint) to give ABC-SMC more signal about the switching behaviour.
