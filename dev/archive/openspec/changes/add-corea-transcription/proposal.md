## Why

Core A′'s three ODE modules make energy. This is the first module that spends it,
and the first that is genuinely discrete: transcripts exist at zero to two copies
per cell, so there is no concentration to integrate and the stochastic block is
not an ODE with noise bolted on. It is also where the reduction's central proxy
lives — promoter strength inherited from proteomics — which is the thing the
surrogate work exists to replace and therefore the baseline it will be measured
against.

**Wave:** 1, CME block, per `dev/plans/reduced-syn3a-wave-plan.md`. This is the
`add-corea-transcription` row of that file's wave-1 table.

## Dependencies

- **`establish-corea-interface`** (wave 0). Required and already satisfied. This
  change consumes the registry, the edge kinds and the loader unchanged.

No other change is required. `add-nucleotide-recycling` owns the ATP and GTP
pools this module debits, and `add-central-glycolysis`, `add-pts-transport` and
`add-nucleotide-recycling` between them own the enzymes the seventeen genes
encode — but all of that is declared through typed edges with peers unnamed, and
this module owns no registry state, so it composes and runs alone.

## What Changes

- **Seventeen genes, one constitutive transcription reaction each.** Ten
  glycolytic, four phosphotransferase, three nucleotide-recycling. Each reaction
  is catalytic in its gene and produces one transcript plus the cost counters the
  hook later debits:

  ```
  gene  →  gene + mRNA + n·ATP_trsc + #A·ATP_mRNA + #U·UTP_mRNA
                        + #G·GTP_mRNA + #C·CTP_mRNA
  ```

- **The formalism is a jump process**, and inference on it is simulation-based
  rather than gradient-based. That is not a preference: at fractional copy number
  a continuous approximation has no meaning.

- **This module owns no registry species.** Its states are the seventeen
  transcripts and the five cost counters, none of which the Core A′ registry names
  — the registry is a metabolite and phospho-state vocabulary, and its `:mrna`
  copy-number regime is declared there precisely for this module to use without
  the transcripts themselves being registry entries.

- **NTP costs are deferred counters against the pools `add-nucleotide-recycling`
  owns.** Five of them, because the published model charges transcription twice
  over: `n` ATP hydrolysed per nucleotide polymerised, *and* one NTP monomer
  incorporated per base. The first drains ATP to ADP and phosphate; the four
  monomer counters drain their NTP to pyrophosphate. All five take the published
  clamped-at-zero drain with the deficit carried forward.

- **CTP and UTP are clamped, and that clamp is ours.** Their sources are in the
  nucleotide module that Core A′ cut, so they are chemostatted by this reduction
  and not by the published model — unlike external glucose, which the published
  model clamps too. They are held at the nucleotide balanced file's values, 0.6874
  and 2.7681 mM.

- **Rate constants are settable, not baked in.** Every one of the seventeen is
  recomputed from live NTP pools every 60 s in the published model, so this module
  declares four inbound rate-constant edges at that cadence and exposes a
  recomputation entry point. `add-cme-rebuild-60s` calls it; nothing here does.

- **Promoter strength is the `/180` proteomics proxy, inherited and flagged.**

### Constitutive, not two-state — and why that matters

The published motif has no promoter switching. There is no `k_on`, no `k_off`,
no on-state and off-state to be degenerate with each other; the gene is always
available and the transcription reaction fires at one rate constant. This module
implements that and nothing more.

The consequence is a correction the scoping note already records and this change
carries into the code: **the ISAB bursty-transcription identifiability result does
not transfer.** That result is about a two-state promoter, and the degeneracy it
characterises requires two states to be degenerate between. Importing it here
would be importing a conclusion about a different model.

This is worth stating in the code rather than only in a note, because "17 genes ×
one stochastic transcription reaction" looks exactly like the bursty motif the
project has existing machinery for, and reaching for that machinery is the
natural mistake.

### Promoter strength is a proxy, and conditioning on proteomics would be circular

`getPtnCount` reads each gene's protein copy number from proteomics and the rate
constant scales as that count over 180. So promoter strength is not measured; it
is *inferred from the steady-state protein abundance the model is supposed to
predict*.

Core A′ inherits it, deliberately, as the baseline the surrogate work is measured
against — replacing this proxy is the point of that work, and a baseline you have
not implemented cannot be improved on. But two things follow, and both must be
labelled wherever a result depends on them:

- **It is a proxy, not a parameter.** The seventeen values are a re-encoding of
  proteomics, so recovering them from synthetic data tests the machinery, not the
  biology.
- **Conditioning on proteomics would be circular.** Protein counts are the input
  to these rate constants. An inference that conditions on measured protein
  abundance and reports having recovered promoter strengths has recovered its own
  input. The observable design in wave 3 has to avoid it, and this module records
  the constraint so that it is visible at the point the choice is made.

### The decision of record: a base-to-NTP mapping bug, corrected

`TranscriptRate` charges three of the four bases against the wrong NTP pool. The
local variable names and the base-count dictionary keys are written in different
orders, so the assignment cycles:

| Assigned | Gets the count of | Charged against | Should be |
|---|---|---|---|
| `NMono_A` | A | ATP | A — correct |
| `NMono_U` | **C** | UTP | U |
| `NMono_C` | **G** | CTP | C |
| `NMono_G` | **U** | GTP | G |

It is present on both code paths, including the production 60 s rebuild.

**Core A′ corrects it, and labels the correction as ours.** The effect on the
rate constant itself is about one percent, because the denominator is dominated by
transcript length. The effect on the coupling is not: over the seventeen genes the
base counts are A 7236, C 2078, G 3094, U 5868, so the bug weights the GTP term by
#U rather than #G and makes the transcription rate constant **1.9× more sensitive
to the live GTP pool than it should be**. That is the only ODE→CME channel in the
model, so porting the bug would mean demonstrating bidirectional coupling with a
gain that a typo nearly doubled.

Correcting it is a departure from the published model and is registered as such,
with both mappings recorded so the difference can be measured rather than argued.

### An honest measurement of what channel 4 is worth

Reporting this now rather than at wave 3, because it changes how the result should
be framed and it falls out of the rate law directly.

The transcription rate constant's denominator is
`(1 + K₀/[RNAP])·K_d²/(C₁C₂) + Σ N_x K_d/[NTP_x] + n − 1`. For these seventeen
genes the NTP-dependent term contributes **4.4% to 5.1%** of that denominator; the
rest is transcript length. So the elasticity of every transcription rate constant
with respect to the NTP pools is about 0.045, and with respect to GTP alone about
0.008 to 0.012 once the base mapping is corrected — roughly 0.016 to 0.020 if it
is not.

Channel 4 is therefore **live and bidirectional but weak**: doubling every NTP
pool moves a transcription rate constant by about 3%. That satisfies the design
criterion, which asks whether information flows both ways, and it bounds how much
information can. Wave 3 should size its expectations accordingly rather than
discovering this after a coverage result comes back uninformative.

## Capabilities

### New Capabilities

- `corea-transcription/transcription-reactions`: seventeen genes, one constitutive
  reaction each, the discrete formalism, and the recorded fact that the motif is
  not two-state.
- `corea-transcription/transcription-rate-constants`: the published rate-constant
  formula, the corrected base-to-NTP mapping, the `/180` promoter proxy with its
  circularity warning, and the requirement that every constant be recomputable.
- `corea-transcription/ntp-cost-accounting`: the five cost counters, what each
  drains into, and the clamped drain policy they inherit.
- `corea-transcription/transcription-coupling`: the eleven typed edges this module
  declares, that it owns no registry state, and that it runs standalone.

### Modified Capabilities

None. `corea-interface` is consumed unchanged.

## Impact

**New source**

| Path | Contents |
|---|---|
| `src/organisms/coreA/transcription.jl` | `CoreATranscription <: AbstractSubModel` |
| `src/organisms/coreA/data/transcription_genes.tsv` | per gene: transcript length, base composition, protein count, measured mean mRNA |
| `dev/scripts/extract_transcription_genes.jl` | regenerates it from a source-model checkout |

**Changed source**

- `src/InferCell.jl` — one `include` after the registry.

**New tests**

- `test/test_corea_transcription.jl`, included from `test/runtests.jl`.

**Read-only to this branch**

`src/edges.jl`, `src/resolver.jl`, `src/loader.jl`, `src/labels.jl`,
`src/interface.jl`, `src/orchestrator.jl`, `src/organisms/coreA/registry.jl`.

**The first jump sub-model in Core A′.** `build_problem` already dispatches on
`formalism = :jump`, and the existing stochastic gene-expression models use it, so
this is exercising an existing path rather than adding one. But it is the first
Core A′ module to do so, and the first to declare `inference_mode = :simulation`,
so the composition of a jump module with the three ODE modules is genuinely untried
— and is wave 2's problem, not this change's.

**A gradient obstruction, declared where it arises.** All five deferred counters
take the published clamped-at-zero drain, which is non-differentiable at the
boundary. The wave-0 contract requires that a composition preparing a clamped
counter for a differentiable sub-model reports it. This module's own inference mode
is simulation-based, so nothing is obstructed here; the report will fire when this
module meets the ODE blocks in wave 2, which is the correct place for it to fire.

**Documentation**

`dev/notes/reduced-syn3a-scoping.md` gains the base-mapping bug and its effect,
and the channel-4 elasticity measurement, which bears directly on the note's
argument that channel 4 satisfies design criterion 2.
