```@meta
CurrentModule = InferCell
```

# The reduced syn3A model

A small organism cut from the published JCVI-syn3A whole-cell model of Thornburg
et al. (2022). It keeps what makes a whole-cell model hard to infer on: a genuine
split between deterministic metabolism and stochastic gene expression, coupling in
both directions, low-copy discrete species, and autocatalytic closure. About 92% of
the cell is removed. What is left is the best-measured region of the network.

In the code, the `CoreA` and `corea` prefixes name this model: [`build_corea`](@ref),
[`corea_models`](@ref), `src/organisms/coreA/`.

## What it contains

| | |
|---|---|
| ODE block | 22 reactions in 4 modules: glycolysis through lactate (10), the phosphotransferase cascade plus lactate export (6), nucleotide recycling (5), and a lumped tRNA charging step (1) |
| Stochastic block | 52 reactions in 3 modules: transcription, translation and transcript decay for 17 genes, plus one translocation for ptsG, the only membrane protein |
| State | 32 dynamic species and 5 chemostats, named and ordered once in a registry |
| Copy numbers | mRNA 0–2, protein 266–1,355, metabolites about 200–74,000 particles |
| Coupling | 3 channels exchanged every 1 s, 2 every 60 s, and a volume channel in both directions |
| Growth | Surface-area accounting. The cell reaches about 1.07× its initial volume in one 6,300 s cycle |

The seven modules are [`CentralGlycolysis`](@ref), [`PtsTransport`](@ref),
[`NucleotideRecycling`](@ref) and [`TrnaCharging`](@ref) in the ODE block, and
[`CoreATranscription`](@ref), [`CoreATranslation`](@ref) and
[`CoreATranscriptDecay`](@ref) in the jump block. [`corea_models`](@ref) returns
them as the assembled model composes them.

## How the blocks talk

**Through state, not through parameters.** No rate constant appears in both
blocks. Enzyme copy numbers from translation set metabolic rates. Nucleotide and
charged-tRNA pools from metabolism set transcription and translation rates. This
is a property of the published model, not of the reduction, and it shapes the
sampler: there is no shared parameter whose posterior could be passed across the
seam, so information crosses through the path.

Each module declares the states it owns and its typed coupling edges. The resolver
checks the composition before it runs: every state has an owner, every edge a
partner, and no moiety is stranded. See the
[interface contract](../api/corea-interface.md).

**The handshake.** Every simulated second, counts from the jump block become
rate-law parameters of the ODE block, the ODE block integrates one second, and the
costs the jump block accrued are debited from the metabolite pools. Every sixty
seconds, the jump block's rate constants are rebuilt from those pools. The [handshake driver](../api/handshake.md) runs this, and records it
from the outside so conservation and clipping can be checked.

## Two variants of one model

- **The published model**, `build_corea()`, clamps a drained pool at zero and
  carries the shortfall, and writes pools back as whole particles. It generates
  every synthetic dataset.
- **The smoothed model**, `build_corea(; COREA_SAMPLER...)`, smooths the clamp over
  one particle and carries pools continuously. It is differentiable through the
  handshake, and kept for gradient-based work. The conditional sampler does not
  need it: its parameter block is gradient-free and runs on the published model.

## What it is not

It does not reproduce syn3A's 105-minute doubling time, its metabolite
concentrations or its proteome, and no result is to be compared with them.
[`doubling_time`](@ref) throws rather than return a number; growth is reported as
fractional growth or time to a threshold ([`growth_report`](@ref)).

Several treatments are ours and not the published model's. These include the lumped
charging step, the chemostatted CTP, UTP and amino-acid pools, and exogenous
membrane growth. [`reduction_report`](@ref) enumerates each, so the label can reach
any result that depends on one.

## M0, the two-gene reference system

The sampler is validated on M0, the same model cut to two genes over 600 s: GAPD
(`JCVISYN3A_0607`) and ptsG (`JCVISYN3A_0779`). The other enzyme counts are held at
their proteomics copy numbers ([`HeldEnzymeCounts`](@ref)), so the catalytic
channel runs as it does in the full model with one live count.

M0 has five free parameters ([`M0_TARGETS`](@ref)): the two promoter strengths, the
transcript decay constant, and the enolase and aldolase forward constants. The sixth
is the metabolite noise scale. [`m0_dataset`](@ref) draws a truth from the prior and
simulates a population of cells under it. Metabolites are observed as one bulk
measurement per pool and time, and transcripts exactly per cell.

## Further reading

- `dev/notes/reduced-syn3a-scoping.md` in the repository derives the reduction.
- `spec/spec.md` carries the claim, the validity checks and every decision of
  record.
