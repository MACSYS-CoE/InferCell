## Why

Every conserved moiety in Core A′ needs a closed cycle, and without these five
reactions the model is a set of dead ends rather than a sustained pathway. Two of
them make GTP, which is what translation runs on; three of them return spent
adenylate, phosphate and guanylate to the top of the loop.

This is not a hypothetical. The first Core A specification declared costs against
AMP and PPi with nothing consuming either, and the arithmetic is brutal: doubling
the proteome takes 3,484,518 charging events — 553 per second — against a total
adenylate pool of about 79,800 particles. The cell runs out of adenylate after
**144 seconds of a 6,300-second cycle**, or 2.3% of it. PPi, unhydrolysed, reaches
173 mM. Those are the errors this change exists to make impossible.

**Wave:** 1, ODE block, per `dev/plans/reduced-syn3a-wave-plan.md`. This is the
`add-nucleotide-recycling` row of that file's wave-1 table.

## Dependencies

- **`establish-corea-interface`** (wave 0). Required and already satisfied. This
  change consumes the species registry, the seven edge kinds and — more than any
  other wave-1 change — the provenance-carrying loader's governing-file machinery,
  and modifies none of them.

No other change is required. `add-central-glycolysis` and `add-pts-transport` are
parallel siblings, not prerequisites. The coupling to glycolysis is real and
dense — this module shares two genes with it and four metabolites — but it is
declared through typed edges with peers unnamed, so the module resolves and
integrates alone.

## What Changes

- **Five reactions**, as the published reconstruction writes them, water and
  hydrogen ions already absent because Core A′ derives from the `NoH2O` model:

  ```
  PGK3   M_13dpg_c + M_gdp_c  ⇌  M_3pg_c + M_gtp_c
  PYK3   M_gdp_c + M_pep_c    ⇌  M_gtp_c + M_pyr_c
  ADK1   M_amp_c + M_atp_c    ⇌  2 M_adp_c
  PPA    M_ppi_c              ⇌  2 M_pi_c
  GK1    M_atp_c + M_gmp_c    ⇌  M_adp_c + M_gdp_c
  ```

- **Three new genes**: adk `JCVISYN3A_0651`, ppa `JCVISYN3A_0344`, gmk
  `JCVISYN3A_0203`. `PGK3` and `PYK3` are catalysed by `JCVISYN3A_0606` and
  `JCVISYN3A_0221`, the same genes as `PGK` and `PYK`, so they **add reactions but
  no genes**. That shared catalysis restores a little of the enzyme competition
  the reduction otherwise removes, and it means the two modules must derive one
  enzyme concentration the same way.

- **Eight dynamic states owned**: the registry's four `:adenylate` species
  (`M_atp_c`, `M_adp_c`, `M_amp_c`, `M_pi_c`), its three `:guanylate` species
  (`M_gtp_c`, `M_gdp_c`, `M_gmp_c`), and `M_ppi_c`. These are the pools every
  other Core A′ module routes its energy traffic through, so this module is the
  counterpart of every currency and deferred-counter edge in the model.

- **Every import declares its governing file.** All five reactions and all eight
  initial conditions exist in both balanced files. See below — this is the
  substance of the change, not bookkeeping around it.

- **Two conservation checks, both this module's own.** Phosphate closure, counting
  `M_pi_c` plus every phosphorylated species with `M_ppi_c` as two; and adenylate
  and guanylate conservation, run over a full 6,300-second cycle rather than a
  burst, because the failure mode is a slow one-way drain that a short run cannot
  see.

### Why each of the three recycling reactions is required

Recorded with the magnitude that makes it required, because "closes a moiety" is
the kind of justification that gets trimmed:

- **ADK1 is the fatal one.** The lumped charging step turns ATP into AMP, and
  nothing else in Core A′ touches AMP. 553 AMP per second against ~79,800
  adenylate particles exhausts the pool in 144 s. ADK1 returns AMP to ADP, which
  PGK and PYK rephosphorylate.
- **PPA is the same argument in the other product.** One PPi per charging event,
  unhydrolysed, is 173 mM — lactate-scale accumulation against a 0.1 mM starting
  pool.
- **GK1 closes the smaller leak.** Transcription buries GTP as GMP in mRNA and
  decay releases it free; with PGK3 and PYK3 rephosphorylating only GDP, that GMP
  is stranded. Roughly 24,000 GMP over a cycle against a ~39,800-particle
  guanylate pool — about 60%. Not fatal within one cycle, but enough to bias the
  transcription rate constants that the 60 s rebuild computes from live NTP pools,
  which is precisely the coupling Core A′ exists to exercise.

### The decision of record: the cross-file trap is wider than recorded

The scoping note flags `PGK3` and `PYK3` as appearing in both balanced files with
different modes. **All five do**, and three of them disagree by more than the two
already on the record:

| Reaction | nucleotide `Mode` | central `Mode` | ratio | central gstd |
|---|---|---|---|---|
| PGK3 (fwd) | 140.8 | 319.4605 | 2.3× | 2.11e+63 |
| PYK3 (fwd) | 672.84 | 1873.9829 | 2.8× | 2.11e+63 |
| ADK1 (fwd) | 319.156 | 228.0578 | 1.4× | 7.01e+35 |
| **GK1 (fwd)** | **410.227** | **7.5991** | **54×** | 2.11e+63 |
| **PPA (fwd)** | **646.727** | **583611.6071** | **902×** | 1.69e+33 |

The nucleotide file governs for all five: they are nucleotide-module reactions,
and the central file's numbers are not a competing estimate but an artefact of
balancing with no data. The geometric standard deviations say so mechanically —
1e33 to 1e63, against 1.0513 in the nucleotide file — and all seventeen of these
reactions' Michaelis constants sit at exactly 0.1 mM with gstd 10, the prior
median at prior width, in the central file.

**Taking the central file's PPA constant would have run inorganic
pyrophosphatase 902× too fast**, which is the kind of error that hides because it
makes a conservation check pass. So the governing declaration is not defensive
paperwork; on this module it is the difference between a working model and a
plausible-looking broken one.

The same applies to the initial conditions, and cuts both ways:

| Species | governs | value | central file shows |
|---|---|---|---|
| `M_atp_c`, `M_adp_c`, `M_pi_c` | central | 3.6529, 0.2178, 17.8185 | identical |
| `M_ppi_c` | central | 0.1 (gstd 10) | identical, and uninformed in both |
| `M_amp_c` | nucleotide | 0.0832 | 0.1 at gstd 10 |
| `M_gtp_c` | nucleotide | 1.6627 | 0.1 at gstd 10 |
| `M_gdp_c` | nucleotide | 0.2981 | 0.1 at gstd 10 |
| `M_gmp_c` | nucleotide | 0.0117 | 0.1 at gstd 10 |

Four agree across the files and four do not, so this module exercises both halves
of the loader's ambiguity report — the disagreements and the agreements — and all
eight need a governing declaration regardless, because the loader refuses to
choose silently when two files hold one identifier. Reading the guanylate species
from the central file would understate that pool by an order of magnitude:
~39,800 particles rather than the ~4,000 that 0.1 mM implies.

### One thing that is not a trap here

The Michaelis-constant column that split `add-central-glycolysis` — the
`Quantity` table the simulator reads against the balanced `Parameter` table that
carries the priors — **agrees on all seventeen constants in the nucleotide file**.
There is no column choice to make and no deviation to label. That is worth
stating, because it means the two wave-1 ODE modules end up with different
provenance stories through no inconsistency of ours.

## Capabilities

### New Capabilities

- `corea-nucleotide-recycling/recycling-network`: the five reactions, the eight
  states they own, the shared catalysis with glycolysis, and the recorded reason
  each recycling reaction is required.
- `corea-nucleotide-recycling/cross-file-provenance`: every value imported with a
  declared governing file, the disagreements and agreements both enumerable, and
  the informedness of a value distinguishing a balanced estimate from a prior
  default.
- `corea-nucleotide-recycling/moiety-conservation`: adenylate and guanylate
  conservation over a full cycle, phosphate closure, and the demonstration that
  removing ADK1 reproduces the dead end of record.
- `corea-nucleotide-recycling/recycling-coupling`: what this module declares as
  the owner of Core A′'s energy pools, and that it resolves standalone.

### Modified Capabilities

None. `corea-interface` is consumed unchanged.

## Impact

**New source**

| Path | Contents |
|---|---|
| `src/organisms/coreA/nucleotide_recycling.jl` | `NucleotideRecycling <: AbstractSubModel` |
| `src/organisms/coreA/data/nucleotide_recycling.tsv` | the nucleotide-file extract |
| `src/organisms/coreA/data/nucleotide_recycling_central.tsv` | the central file's rival values for the same identifiers, so the ambiguity is real rather than described |
| `dev/scripts/extract_nucleotide_recycling.jl` | regenerates both |

**Changed source**

- `src/InferCell.jl` — one `include` after the registry.

**New tests**

- `test/test_corea_nucleotide_recycling.jl`, included from `test/runtests.jl`.

**Read-only to this branch**

`src/edges.jl`, `src/resolver.jl`, `src/loader.jl`, `src/labels.jl`,
`src/interface.jl`, `src/orchestrator.jl`, `src/organisms/coreA/registry.jl`.

**A coordination point with `add-central-glycolysis`.** `PGK3` shares
`JCVISYN3A_0606` with `PGK`, and `PYK3` shares `JCVISYN3A_0221` with `PYK`. Both
modules therefore carry a nominal concentration for the same enzyme, and if the
two derivations drift apart the assembled model runs one enzyme at two
concentrations. Both use copy number at the registry's initial volume — 411 and
551 copies, 0.020367 and 0.027304 mM — and this change asserts the agreement
rather than assuming it.

**The same orchestrator limitation the sibling changes record.** `_build_rhs`
(`src/orchestrator.jl:203-215`) lets a module write only its own state slice, so
this module's production of `M_3pg_c` and `M_pyr_c` is declared and left
unexecuted; standalone, `M_13dpg_c` and `M_pep_c` arrive as held inputs.

**Documentation**

`dev/notes/reduced-syn3a-scoping.md` records the cross-file inconsistency for
`PGK3` and `PYK3` only, in both its open-questions list and its GTP-regeneration
section. It gains the other three reactions, their ratios, and the observation
that the central file's geometric standard deviations identify the out-of-scope
values mechanically.
