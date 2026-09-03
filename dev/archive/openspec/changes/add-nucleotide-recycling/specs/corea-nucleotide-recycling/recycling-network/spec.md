## Purpose

The two reactions that make GTP and the three that return spent adenylate,
phosphate and guanylate to the top of the loop, together with the eight pools
they own — so that every cost Core A′ declares has a state that can pay it and a
reaction that returns it, and so the reason each recycling reaction exists
survives in the code rather than only in a note.

## ADDED Requirements

### Requirement: The five reactions and their stoichiometry

The system SHALL provide a Core A′ sub-model carrying exactly five reactions,
with the stoichiometry of the published reconstruction:

| Reaction | Stoichiometry | Gene (syn3A locus) |
|---|---|---|
| PGK3 | `M_13dpg_c` + `M_gdp_c` ⇌ `M_3pg_c` + `M_gtp_c` | JCVISYN3A_0606 |
| PYK3 | `M_gdp_c` + `M_pep_c` ⇌ `M_gtp_c` + `M_pyr_c` | JCVISYN3A_0221 |
| ADK1 | `M_amp_c` + `M_atp_c` ⇌ 2 `M_adp_c` | JCVISYN3A_0651 |
| PPA | `M_ppi_c` ⇌ 2 `M_pi_c` | JCVISYN3A_0344 |
| GK1 | `M_atp_c` + `M_gmp_c` ⇌ `M_adp_c` + `M_gdp_c` | JCVISYN3A_0203 |

Every species named SHALL be a Core A′ registry species. Neither water nor the
hydrogen ion SHALL appear: Core A′ derives from the `NoH2O` model with
hydrogen-ion accounting removed, which is why pyrophosphatase is carried without
its water substrate.

Each reaction SHALL be reversible with an independent reverse constant, as the
reconstruction writes it. A species appearing with a stoichiometric coefficient
above one SHALL contribute that many terms to the rate law, so that the two
`M_adp_c` of ADK1 and the two `M_pi_c` of PPA enter as repeated factors rather
than as a single term.

#### Scenario: The reaction set is exactly five
- **WHEN** the sub-model's reactions are enumerated
- **THEN** exactly the five reactions above are present, each exactly once

#### Scenario: Water and hydrogen ions are absent
- **WHEN** each reaction's substrates and products are collected
- **THEN** no species denoting water or the hydrogen ion appears
- **AND** every species resolves to a Core A′ registry position

#### Scenario: Repeated stoichiometry enters the rate law twice
- **WHEN** ADK1's rate is evaluated
- **THEN** its product term depends on the square of the ADP saturation ratio
- **AND** the same Michaelis constant is used for both occurrences

### Requirement: Two reactions add no genes

`PGK3` and `PYK3` SHALL be recorded as catalysed by the same genes as `PGK` and
`PYK` respectively. The sub-model SHALL expose which of its reactions share
catalysis with a reaction outside it, and the shared enzyme's concentration SHALL
be derived by the same rule the other module uses, so that one enzyme cannot
appear at two concentrations in a composed model.

#### Scenario: Shared catalysis is declared
- **WHEN** the sub-model is queried for reactions whose gene is also used outside
  it
- **THEN** `PGK3` and `PYK3` are returned with their loci
- **AND** the three recycling reactions are not

#### Scenario: A shared enzyme has one concentration
- **WHEN** a composition contains this sub-model and another declaring the same
  locus
- **THEN** the enzyme concentration each derives for that locus is equal
- **AND** a disagreement is reported rather than silently resolved

### Requirement: Eight dynamic states are owned, and no others

The sub-model SHALL integrate exactly eight of the registry's dynamic states: the
four `:adenylate` species, the three `:guanylate` species, and `M_ppi_c`.

The four glycolytic intermediates its reactions touch — `M_13dpg_c`, `M_3pg_c`,
`M_pep_c` and `M_pyr_c` — SHALL NOT be integrated here. State ordering SHALL
follow the registry's canonical order.

#### Scenario: The owned states are the adenylate, guanylate and pyrophosphate pools
- **WHEN** the sub-model's states are enumerated
- **THEN** they are exactly the registry's `:adenylate` and `:guanylate` groups
  plus `M_ppi_c`
- **AND** `M_pi_c` is among them, because the registry groups phosphate with the
  adenylate species

#### Scenario: Glycolytic intermediates are not owned
- **WHEN** the sub-model's states are queried for `M_13dpg_c`, `M_3pg_c`,
  `M_pep_c` or `M_pyr_c`
- **THEN** none of the four is present

### Requirement: Each recycling reaction records why it is required

`ADK1`, `PPA` and `GK1` SHALL each carry a recorded reason for its presence,
stating the moiety it closes, the rate at which that moiety would otherwise be
lost, and the consequence over a cell cycle. The reasons SHALL be retrievable
from the composed model.

The recorded magnitudes SHALL be at least: for `ADK1`, that around 553 AMP per
second against an adenylate pool of roughly 79,800 particles exhausts it in about
144 seconds of a 6,300-second cycle; for `PPA`, that unhydrolysed pyrophosphate
reaches about 173 mM against a 0.1 mM starting pool; for `GK1`, that roughly
24,000 stranded GMP over a cycle is about 60% of a guanylate pool of some 39,800
particles.

An earlier specification of this model declared costs against AMP and PPi with
nothing consuming either. The record exists so that removing one of these three
requires confronting the number that put it there.

#### Scenario: Each reason is retrievable with its magnitude
- **WHEN** the composed model is queried for this module's recorded rationale
- **THEN** three entries are returned, one per recycling reaction
- **AND** each names the moiety it closes and the magnitude that makes it required

#### Scenario: The GTP-producing reactions are not in that list
- **WHEN** the rationale is enumerated
- **THEN** `PGK3` and `PYK3` do not appear, because they supply a cost rather
  than close a leak
