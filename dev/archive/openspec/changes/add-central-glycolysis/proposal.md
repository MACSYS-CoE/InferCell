## Why

Core A′ needs the pathway everything else in it exists to serve. Glycolysis
through lactate is where glucose becomes ATP and GTP, where the enzymes that the
CME block transcribes and translates actually do work, and where the autocatalytic
loop closes. Until it exists there is no flux for the energy interface to carry,
and the other six wave-1 modules have nothing to couple to.

It is also the best-measured region of the reduced model: ten single-gene
reactions whose kinetics carry genuine balanced priors, which makes it the module
where a parameter *should* recover tightly and therefore the natural control on
whether the inference machinery works at all.

**Wave:** 1, ODE block, per `dev/plans/reduced-syn3a-wave-plan.md`. This is the
`add-central-glycolysis` row of that file's wave-1 table.

## Dependencies

- **`establish-corea-interface`** (wave 0). Required and already satisfied: it is
  merged, archived, and its three delta specs are promoted to
  `openspec/specs/corea-interface/`. This change consumes the species registry,
  the seven edge kinds and the provenance-carrying loader; it does not modify any
  of them.

No other change is required. `add-pts-transport` and `add-nucleotide-recycling`
are parallel siblings in wave 1, not prerequisites: this module names the species
it shares with them through the registry and declares typed edges at the
boundary, and resolves standalone with those edges' peers left unnamed. The
execution of that coupling is `add-hook-1s-coupling`'s job in wave 2.

## What Changes

- **Ten reactions, `PGI` → `LDH_L`, one gene each**, ported from the central
  balanced file: `PGI`, `PFK`, `FBA`, `TPI`, `GAPD`, `PGK`, `PGM`, `ENO`, `PYK`,
  `LDH_L`. Every GPR in this set is a single gene, so none of the published
  model's AND/OR effective-enzyme logic is needed.

- **Thirteen dynamic states owned**: the eleven `:glycolytic` registry species
  (`M_g6p_c`, `M_f6p_c`, `M_fdp_c`, `M_dhap_c`, `M_g3p_c`, `M_13dpg_c`,
  `M_3pg_c`, `M_2pg_c`, `M_pep_c`, `M_pyr_c`, `M_lac__L_c`) plus the two
  `:redox` species `M_nad_c` and `M_nadh_c`. No other module in Core A′
  integrates any of these.

- **`M_atp_c`, `M_adp_c` and `M_pi_c` are boundary, not state.** They are
  `add-nucleotide-recycling`'s to integrate. This module declares them as
  `CurrencyEdge`s — ATP and ADP in both directions (`PFK` draws ATP, `PGK` and
  `PYK` supply it), phosphate inbound only (`GAPD` draws it) — and reads them
  through `inputs`.

- **`NOX` is deliberately dropped, and recorded as a decision rather than an
  omission.** `LDH_L` already regenerates NAD⁺, so `NOX` is redundant for redox
  balance; it carries no GPR rule, so cutting it costs no boundary; and its
  `k_cat` was prior-dominated, so it was never an inference target. The scoping
  note makes this call and this change carries it into the code as a labelled
  decision, so that a later reader finds a reason rather than a gap.

- **Kinetics and initial conditions imported through the wave-0 loader** from a
  vendored, derived extract of the central balanced file's `Mode` column: 10
  `kcatF`, 10 `kcatR`, 32 `K_M` and the 13 owned states' concentrations, each
  carrying its geometric standard deviation and its source file.

- **This module's own conservation check: redox balance.** `GAPD` reduces NAD⁺ to
  NADH and `LDH_L` oxidises it back, and no other reaction in Core A′ touches
  either species — so `M_nad_c + M_nadh_c` is exactly conserved by this module
  alone. It is checked here, on this branch, not deferred to wave 3.

### Two decisions of record

**1. The `K_M` values come from the balanced `Parameter` table, not the
`Quantity` table the published simulator reads.** These 32 constants exist twice
in the source file, and 8 of them disagree:

| Identifier | `Quantity` `Value` (what runs) | balanced `Mode` (ported) | gstd |
|---|---|---|---|
| `kmc_R_PGI_M_g6p_c` | 0.28 | 22.9419 | 1.0512 |
| `kmc_R_PGI_M_f6p_c` | 0.15 | 3.3488 | 1.2986 |
| `kmc_R_PFK_M_atp_c` | 0.117 | 0.0255 | 4.8229 |
| `kmc_R_FBA_M_fdp_c` | 0.005 | 0.2447 | 3.573 |
| `kmc_R_FBA_M_dhap_c` | 0.095 | 0.0064 | 4.5985 |
| `kmc_R_FBA_M_g3p_c` | 1.0 | 0.0044 | 8.4265 |
| `kmc_R_GAPD_M_nad_c` | 1.3 | 5.6969 | 1.1808 |
| `kmc_R_PGM_M_2pg_c` | 1.47 | 0.0278 | 6.0168 |

The other 24 agree exactly. Taking the balanced column gives every `K_M` one
provenance and a real inherited prior, which is what the inference problem needs;
the `Quantity` column is a hand-patch layer carrying no uncertainty at all. The
cost is that Core A′ is not, on these eight numbers, running what the published
simulator runs — a departure large enough to change a trajectory (PGI's G6P
constant moves by 82×, FBA's G3P constant by 227×). It is therefore **ours, not
the published model's**, and is registered through `reduction_notes` so that
`reduction_declarations` enumerates it and it reaches any result that depends on
it.

This is a third instance of the column trap the scoping note already records
twice. The note's own summary — "~35 `K_M` from the balanced file (with priors)"
— assumed the balanced column *is* what runs. It is not, and the scoping note
should be corrected to say so.

**2. The parameter table is vendored as a derived extract.** Wave 0 left the
balanced tables unvendored and the loader path-agnostic. That cannot stand for a
module whose tests must run on compute nodes with no network, and the raw
upstream file does not fit the loader anyway: its balanced `Parameter` table keys
each row on `(QuantityType, Reaction id, Compound id)`, while `read_source_table`
assumes a single ID column — and `src/loader.jl` is read-only to a wave-1 branch.
So this change commits a ~65-row extract under `src/organisms/coreA/data/`,
reshaped to the loader's single-ID convention (including the `conc_<species>`
form, so the loader's registry-agreement check actually fires), alongside the
script that regenerates it from a checkout of the source model at commit
`db048ac`.

## Capabilities

### New Capabilities

- `corea-central-glycolysis/reaction-network`: the ten reactions and their
  stoichiometry, the thirteen dynamic states this module owns, the dropped `NOX`
  decision, and the redox conservation invariant those reactions imply.
- `corea-central-glycolysis/glycolytic-kinetics`: the modular rate law each
  reaction takes, the enzyme concentration that scales it, and the import of
  every kinetic constant and initial condition through the wave-0 loader with its
  source file, its uncertainty and its informedness intact.
- `corea-central-glycolysis/glycolysis-coupling`: what this module declares at its
  boundary — typed edges for the registry species it shares with the other Core A′
  modules, untyped inputs for the enzyme counts no edge kind can name, and the
  requirement that the module resolve standalone.

### Modified Capabilities

None. `corea-interface` is consumed unchanged; a module that needs one of its
requirements altered has found an interface bug, which is a conversation on `main`
rather than an edit on this branch.

## Impact

**New source**

| Path | Contents |
|---|---|
| `src/organisms/coreA/central_glycolysis.jl` | `CentralGlycolysis <: AbstractSubModel` — states, parameters, dynamics, coupling, inputs, reduction notes |
| `src/organisms/coreA/data/central_glycolysis.tsv` | the vendored derived extract, ~65 rows |
| `src/organisms/coreA/data/README.md` | upstream file, commit, and what the extract does and does not reshape |
| `dev/scripts/extract_central_glycolysis.jl` | regenerates the extract from a source-model checkout |

**Changed source**

- `src/InferCell.jl` — one `include` after the registry. Include order is a
  reading convention, not a load-time constraint.

**New tests**

- `test/test_corea_central_glycolysis.jl`, included from `test/runtests.jl`.

**Read-only to this branch**

`src/edges.jl`, `src/resolver.jl`, `src/loader.jl`, `src/labels.jl`,
`src/interface.jl`, `src/orchestrator.jl`, `src/organisms/coreA/registry.jl`.

**One interface limitation this change confirms rather than fixes.**
`_build_rhs` in `src/orchestrator.jl` assembles `du` by concatenating each
sub-model's own state slice, so a module can only write derivatives for states it
owns. This module's ATP and ADP production is therefore *declared* — as outbound
currency edges — but not *executed*: standalone, and in any composition built by
the current orchestrator, `M_atp_c` reaches its dynamics as a held input. Closing
that is exactly what `add-hook-1s-coupling` means by "the join", and wave 2 should
know the gap is structural rather than an oversight in this module.

**Documentation**

`dev/notes/reduced-syn3a-scoping.md` gains a correction-of-record row for the
`K_M` column, in the inline table the note already keeps for this purpose.
