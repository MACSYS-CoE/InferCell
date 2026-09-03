## Context

See `proposal.md` — Why. What differs from the three ODE modules already
proposed:

- **This module is a jump process**, so `dynamics` is not implemented at all and
  `reactions` is. `build_problem` already dispatches on `formalism = :jump`; the
  existing stochastic gene-expression models exercise that path.
- **It owns no registry state**, so most of the resolver's machinery — ownership,
  dead ends on states it integrates, inputs consistency — has nothing to act on
  from its side. Its whole boundary is deferred counters, rate-constant edges and
  clamps.
- **Its parameters are computed, not imported.** The seventeen rate constants are
  a formula over per-gene sequence data and four NTP concentrations, not table
  lookups. That is why they can be rebuilt, and why the data extract carries
  sequence composition rather than kinetic constants.
- **The sequence data is not in the balanced tables.** It comes from the genome
  record, the proteomics table and the measured mRNA counts, so the extract has
  three upstream sources and the loader's cross-file machinery has nothing to do
  here.

## Goals / Non-Goals

**Goals**

- Seventeen reactions and their rate constants, simulable alone over a full cycle.
- The corrected base mapping, with the published one still selectable so the
  difference is measurable.
- Five cost counters declared against the pools that pay them.
- Both caveats — the constitutive motif and the proteomics circularity — recorded
  where a later change will trip over them.

**Non-Goals**

- Performing the 60 s rebuild. This module exposes recomputation;
  `add-cme-rebuild-60s` schedules it.
- Debiting the counters. That is the hook's, in `add-hook-1s-coupling`.
- Translation and decay. Same seventeen genes, different changes.
- Choosing observables, or deciding how to avoid the proteomics circularity. This
  module records the constraint; wave 3 acts on it.

## Decisions

### D1 — The reaction, exactly as the published model builds it

`MinCell_CMEODE.py:728` adds, per gene:

```
gene → gene + mRNA + n·ATP_trsc + #A·ATP_mRNA + #U·UTP_mRNA + #G·GTP_mRNA + #C·CTP_mRNA
```

with propensity `k_transcription × [gene]`. The gene is reactant and product, so
its count is invariant, and the propensity is constant — a `ConstantRateJump`.

Genes are carried as a fixed per-gene quantity rather than as seventeen states
pinned at one. Core A′ cuts replication, so nothing can change them, and seventeen
constant states in a jump system is seventeen wasted dimensions in every
trajectory the inverse problem stores. The published model carries them as species
because it has `rep_start.py`; D-record that, so adding replication later knows to
promote them.

States are therefore 17 transcripts + 5 counters = 22, all outside the registry.

### D2 — The rate constant, and the one thing changed in it

```
kcat_mod = min(rnaPolKcat · S_g , 180)                    S_g = ptnCount_g / 180
NMonoSum = K_d · Σ_x  N_x / [NTP_x]
k_g      = kcat_mod / ( (1 + K₀/[RNAP])·K_d²/(C₁C₂) + NMonoSum + n_g − 1 )
```

with `rnaPolKcat = 0.155·187/493·20 = 1.17586` nt/s, `K₀ = 1e-4` mM, `K_d = 0.1`
mM, `[RNAP] = 187` copies = 0.00927 mM, and `C₁`, `C₂` the NTP concentrations of
the transcript's first two bases.

**The one change: `N_x` is the count of the base polymerised from NTP x.** The
published code assigns `NMono_U = baseCount["C"]`, `NMono_C = baseCount["G"]`,
`NMono_G = baseCount["U"]` — the variable names run A, U, C, G while the keys run
A, C, G, U, so three of four are permuted. `MinCell_CMEODE.py:377-383` and
`MinCell_restart.py:547-553`, so the rebuild inherits it.

`base_mapping = :corrected` is the default; `:published` reproduces the permutation
and is recorded as a selection so a trajectory computed under it is identifiable.
Correcting is a departure from the published model, so `:corrected` is what appears
in `reduction_notes` — the unusual case where the default is the labelled one.

The cap never binds: the largest promoter-scaled turnover among these genes is
GAPD's `1.17586 × 7.528 = 8.85` nt/s against a ceiling of 180. The module asserts
that rather than assuming it.

### D3 — NTP concentrations: balanced everywhere, and never the setup constants

Four concentrations enter every rate constant. Core A′ takes all four from
balanced tables:

| | value | gstd | source | in Core A′ |
|---|---|---|---|---|
| ATP | 3.6529 | 1.2825 | central balanced | live, owned by `add-nucleotide-recycling` |
| GTP | 1.6627 | 1.5684 | nucleotide balanced | live, owned by `add-nucleotide-recycling` |
| CTP | 0.6874 | 2.0574 | nucleotide balanced | clamped, ours |
| UTP | 2.7681 | 1.3664 | nucleotide balanced | clamped, ours |

**Not** the 1.04 / 0.68 / 0.34 / 0.68 mM in `MinCell_CMEODE.py:277-280`. Those are
the first-minute rate-constant inputs that the 60 s rebuild replaces with live
pools, and the scoping note already records one error of exactly this shape — ATP's
1.04 mM mistaken for the metabolic pool. The nucleotide file's `Compound` table
repeats 0.34 and 0.68 for CTP and UTP, which makes the wrong choice look
corroborated; it is the same setup constant twice, not two measurements.

ATP and GTP standalone take the registry's values, which are the same balanced
numbers, so before the rebuild exists the module is self-consistent with the pool
owner.

### D4 — Rate constants are state on the struct, rebuilt through one entry point

`recompute_rate_constants!(m; atp, ctp, gtp, utp)` recomputes all seventeen in
place and returns them; `rate_constants(m)` reads them. Nothing inside this module
calls the first — the module declares a 60 s piecewise-constant cadence and
provides the means, and `add-cme-rebuild-60s` does the scheduling.

That split is what the rate-constant edge kind means: "the module whose rate
constants are rebuilt declares it inbound". A test asserts the module does not
rebuild itself, by simulating past the declared cadence and checking the constants
are unchanged.

Mutable state on a sub-model is a departure from the four ODE modules, which are
immutable. It is unavoidable — the constants must change between refreshes without
the sub-model being rebuilt — and it is confined to one vector with one setter.

### D5 — The promoter proxy, inherited and doubly flagged

`S_g = ptnCount_g / 180`, where `ptnCount` is `max(10, round(proteomics))`. The
values reproduce the scoping note's `S = cnt/180` column exactly:

| | PGI | PFK | FBA | TPI | GAPD | PGK | PGM | ENO | PYK | LDH_L |
|---|---|---|---|---|---|---|---|---|---|---|
| copies | 266 | 458 | 775 | 410 | 1355 | 411 | 322 | 998 | 551 | 1100 |
| S | 1.478 | 2.544 | 4.306 | 2.278 | 7.528 | 2.283 | 1.789 | 5.544 | 3.061 | 6.111 |

| | ptsI | Crr | ptsH | ptsG | ADK1 | PPA | GK1 |
|---|---|---|---|---|---|---|---|
| copies | 353 | 314 | 290 | 831 | 213 | 190 | 186 |
| S | 1.961 | 1.744 | 1.611 | 4.617 | 1.183 | 1.056 | 1.033 |

Two `reduction_notes` entries, not one. The first says it is a proxy: promoter
strength is not measured, it is back-calculated from the steady-state protein
abundance the model is meant to predict. The second says conditioning on
proteomics would be circular, and names the seventeen parameters affected — the
trap is not that the proxy is wrong but that an observable chosen later could
close a loop through it without anyone noticing.

The same protein counts feed the metabolic modules' enzyme concentrations, so a
test asserts this module's copy numbers agree with theirs. One number, three
consumers.

### D6 — Five counters, and what each drains into

| Counter | Increment per firing | Debited against | Produces |
|---|---|---|---|
| `ATP_trsc` | `n_g` | `M_atp_c` | `M_adp_c` + `M_pi_c` |
| `ATP_mRNA` | `#A` | `M_atp_c` | `M_ppi_c` |
| `GTP_mRNA` | `#G` | `M_gtp_c` | `M_ppi_c` |
| `CTP_mRNA` | `#C` | `M_ctp_c` | `M_ppi_c` |
| `UTP_mRNA` | `#U` | `M_utp_c` | `M_ppi_c` |

From `in_out.py:176-260`. Two things fall out that are worth stating plainly.

**ATP is debited twice per transcript** — once as polymerisation energy at one per
nucleotide, once as an incorporated monomer. For PGI that is 1284 + 521 = 1805 ATP
per transcript. It looks like double-counting and is not: it is the published
model's accounting, and `reduction_notes` says so, because the alternative is
someone "fixing" it.

**Transcription produces pyrophosphate**, four counters' worth. The scoping note
attributes PPi only to the lumped charging step when it argues that PPA is
required. The argument survives — charging dominates — but the accounting is
incomplete, and `add-nucleotide-recycling`'s phosphate closure will see this flux
in wave 2.

All five take `clip = :clamped_deficit_carried`, the published policy and the
default, so none is a labelled deviation. They are non-differentiable, and the
wave-0 contract fires that report at the join rather than here, because this
module's own `inference_mode` is `:simulation`.

### D7 — Eleven edges, no inputs

| Kind | Species | Dir | Detail |
|---|---|---|---|
| rate constant | `M_atp_c`, `M_ctp_c`, `M_gtp_c`, `M_utp_c` | in | piecewise-constant, 60 s |
| deferred counter | `M_atp_c` | in | counter `ATP_trsc` |
| deferred counter | `M_atp_c` | in | counter `ATP_mRNA` |
| deferred counter | `M_gtp_c` | in | counter `GTP_mRNA` |
| deferred counter | `M_ctp_c` | in | counter `CTP_mRNA` |
| deferred counter | `M_utp_c` | in | counter `UTP_mRNA` |
| clamped | `M_ctp_c` | in | 0.6874 mM, `origin = :ours` |
| clamped | `M_utp_c` | in | 2.7681 mM, `origin = :ours` |

`inputs` is empty. Nothing here reads another module's state at run time: the NTP
pools enter through the rate constants, which are parameters, not through the
propensities. That is also why a chemostatted species can be named by an edge here
without tripping the resolver's rule against chemostatted inputs.

Three kinds coexist on `M_atp_c` inbound (two counters and a rate-constant edge)
and on `M_ctp_c`/`M_utp_c` (counter, rate constant, clamp). The wave-0 contract
requires exactly this: "the published model routes ATP through a currency pool, a
deferred-counter debit and a rate-constant rebuild simultaneously, and the resolver
accepts that composition."

`ClampedEdge` validates a held value against the registry only where the registry
records one, and it records none for CTP or UTP — so the module carries the value
and its provenance itself, which is what the resolver's own error text prescribes.

### D8 — One extract, three upstream sources

`src/organisms/coreA/data/transcription_genes.tsv`, one row per locus:

```
!ID  !Length  !A  !C  !G  !U  !PtnCount  !MeanMRNA  !UpstreamRow
JCVISYN3A_0445  1284  521  125  197  441  266  0.4403  syn3A.gb|proteomics.xlsx|mRNA_counts.csv
```

Three upstream files, no cross-file ambiguity — each column has exactly one source
— so this extract does not go through `load_parameter`'s governing machinery. It is
read directly, and the four base counts are asserted to sum to the length.

The generator parses the genome record for the CDS of each locus, reverse-complements
where the feature is on the complement strand, transcribes, and counts. It fails
loudly on a missing locus rather than defaulting, since a silently absent gene would
show up only as a model with sixteen transcripts.

### D9 — What the numbers say about channel 4

Computed from the formula at the D3 concentrations, corrected mapping, across all
seventeen genes:

| | range |
|---|---|
| rate constant `k_g` | 1.26e-3 to 8.29e-3 s⁻¹ |
| steady-state transcripts, `k_g/k_deg` | 0.44 to 2.87 copies |
| measured mean transcripts | 0.29 to 2.18 copies |
| elasticity of `k_g` to all four NTP pools | 0.044 to 0.051 |
| elasticity to GTP alone, corrected | 0.0079 to 0.0117 |
| elasticity to GTP alone, published mapping | 0.0156 to 0.0196 |

The steady states land within a factor of about 1.5 of the measured counts, which
is the calibration check that the formula and the extract are wired together
correctly — and the test asserts that band rather than exact values.

The elasticities are the finding. `NMonoSum` contributes about 4.5% of the rate
constant's denominator and transcript length contributes the rest, so the only
ODE→CME channel in Core A′ has a gain of a few percent. It is live and it is
bidirectional; it is also weak, and correcting the base mapping halves the GTP half
of it. Recorded here and in the proposal because wave 3 should size its
expectations before it measures, not after.

## Risks / Trade-offs

- **Channel 4 may be too weak to carry the demonstration.** → Measured and
  reported now rather than discovered in wave 3. Nothing in this change can fix it;
  what it can do is make the number visible and reproducible, which the elasticity
  diagnostic does.

- **Correcting a published bug makes Core A′ not-a-port on this one line**, and
  every other value in the project is chosen as what-actually-runs. → The published
  mapping is retained behind a keyword, the correction is in `reduction_notes`, and
  the comparison is one argument away rather than a re-derivation.

- **Mutable rate constants on a sub-model.** → Confined to one vector with one
  setter and one reader; a test asserts nothing else about the sub-model changes
  when they are rebuilt, and that the module never rebuilds itself.

- **The proteomics circularity is a trap that this change cannot close.** Protein
  counts feed promoter strength here and enzyme concentrations in three other
  modules, so almost any protein-level observable touches it. → Recorded as a
  retrievable caveat naming the affected parameters, so wave 3 meets it when
  choosing observables rather than when interpreting a posterior.

- **Transcription's pyrophosphate is not in the scoping note's accounting.** → The
  drain products are recorded per counter, and the handoff flags it for
  `add-nucleotide-recycling`'s phosphate closure, which will see the flux at the
  join.

## Migration Plan

Additive. One `include` in `src/InferCell.jl` after the registry, one in
`test/runtests.jl`. Appends to `src/organisms/coreA/data/README.md` if a sibling
change has created it. No existing sub-model, test or exported name changes.

## Open Questions

- **Whether the elasticity finding changes the wave-3 plan.** A channel with a 4.5%
  gain still satisfies design criterion 2, but it bounds what an inference can
  learn about the parameters that cross it. Deciding what to do about that needs
  translation and the rebuild in place, since translation's charged-tRNA channel
  may be much stronger. It changes nothing here.

- **Whether the published base mapping should be exercised in the suite or only
  available.** Running both doubles the trajectory cost of this module's tests for a
  comparison nothing yet consumes. Leaning toward available-but-untested until wave
  3 asks for the number.
