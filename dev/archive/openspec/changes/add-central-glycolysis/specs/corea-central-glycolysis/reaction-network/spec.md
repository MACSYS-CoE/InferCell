## Purpose

The ten reactions that carry glucose-6-phosphate through to lactate in Core A′,
the thirteen dynamic states they own, and the redox invariant those reactions
imply — so that the pathway is defined in exactly one place and the conserved
moiety it closes is checked on its own branch rather than at the join.

## ADDED Requirements

### Requirement: The ten glycolytic reactions and their stoichiometry

The system SHALL provide a Core A′ sub-model carrying exactly ten reactions,
with the stoichiometry of the published reconstruction:

| Reaction | Stoichiometry | Gene (syn3A locus) |
|---|---|---|
| PGI | `M_g6p_c` ⇌ `M_f6p_c` | JCVISYN3A_0445 |
| PFK | `M_atp_c` + `M_f6p_c` ⇌ `M_adp_c` + `M_fdp_c` | JCVISYN3A_0220 |
| FBA | `M_fdp_c` ⇌ `M_dhap_c` + `M_g3p_c` | JCVISYN3A_0131 |
| TPI | `M_dhap_c` ⇌ `M_g3p_c` | JCVISYN3A_0727 |
| GAPD | `M_g3p_c` + `M_nad_c` + `M_pi_c` ⇌ `M_13dpg_c` + `M_nadh_c` | JCVISYN3A_0607 |
| PGK | `M_13dpg_c` + `M_adp_c` ⇌ `M_3pg_c` + `M_atp_c` | JCVISYN3A_0606 |
| PGM | `M_3pg_c` ⇌ `M_2pg_c` | JCVISYN3A_0729 |
| ENO | `M_2pg_c` ⇌ `M_pep_c` | JCVISYN3A_0213 |
| PYK | `M_adp_c` + `M_pep_c` ⇌ `M_atp_c` + `M_pyr_c` | JCVISYN3A_0221 |
| LDH_L | `M_nadh_c` + `M_pyr_c` ⇌ `M_lac__L_c` + `M_nad_c` | JCVISYN3A_0475 |

Every species named SHALL be a Core A′ registry species, referred to by its
registry name. The sub-model SHALL NOT define its own species names.

Every reaction SHALL be reversible, matching the reconstruction, so that each
carries both a forward and a reverse catalytic constant.

Each reaction's gene-protein-reaction rule is a single gene, so the
gene-to-enzyme-to-reaction map is one-to-one and no AND/OR effective-enzyme logic
is required. The sub-model SHALL expose that map, so that the module supplying
protein counts can be wired to it without re-deriving it.

#### Scenario: The reaction set is exactly ten
- **WHEN** the sub-model's reactions are enumerated
- **THEN** exactly the ten reactions above are present, each exactly once

#### Scenario: Every species is a registry species
- **WHEN** each reaction's substrates and products are collected
- **THEN** every one of them resolves to a Core A′ registry position
- **AND** no species name is defined by this sub-model itself

#### Scenario: Each reaction maps to one gene
- **WHEN** the gene map is queried for a reaction
- **THEN** it returns exactly one syn3A locus, matching the table above
- **AND** the ten loci are distinct

### Requirement: Thirteen dynamic states are owned, and no others

The sub-model SHALL integrate exactly thirteen of the registry's dynamic states:
the eleven species in the registry's `:glycolytic` group and the two in its
`:redox` group.

Every other registry species its reactions touch — `M_atp_c`, `M_adp_c` and
`M_pi_c` — SHALL NOT be integrated by this sub-model. They are owned elsewhere in
Core A′ and reach this sub-model's dynamics as coupling inputs.

The state ordering SHALL follow the registry's canonical order, so that a
position in this sub-model's state vector maps to a registry position by a rule
rather than by a hand-maintained table.

#### Scenario: The owned states are the glycolytic and redox groups
- **WHEN** the sub-model's states are enumerated
- **THEN** they are exactly the registry's eleven `:glycolytic` species and two
  `:redox` species
- **AND** none of them is marked chemostatted in the registry

#### Scenario: The energy currencies are not owned
- **WHEN** the sub-model's states are queried for `M_atp_c`, `M_adp_c` or
  `M_pi_c`
- **THEN** none of the three is present
- **AND** each is instead available to the dynamics as a coupling input

#### Scenario: State order follows the registry
- **WHEN** the sub-model's state names are mapped to registry positions
- **THEN** the positions are strictly increasing

### Requirement: NOX is dropped as a recorded decision

`NOX` SHALL be absent from the reaction set, and its absence SHALL be recorded as
a decision with its reason rather than left as a gap. The recorded reason SHALL
state that `LDH_L` already regenerates NAD⁺, that `NOX` carries no
gene-protein-reaction rule and therefore crosses no boundary, and that its
catalytic constant was prior-dominated and so was never an inference target.

This is a departure from the published reaction list, not from the published
model's behaviour, and SHALL be retrievable from the composed model alongside the
other Core A′ reduction declarations.

#### Scenario: NOX is not in the reaction set
- **WHEN** the reaction set is enumerated
- **THEN** no reaction named `NOX` is present

#### Scenario: The reason for dropping NOX is retrievable
- **WHEN** a composed model containing this sub-model is queried for the
  declarations that are this reduction's rather than the published model's
- **THEN** the dropped `NOX` reaction appears, with its reason

### Requirement: The NAD⁺/NADH pool is conserved by this sub-model alone

`GAPD` reduces NAD⁺ to NADH and `LDH_L` oxidises it back, and no other reaction
in this sub-model, and no coupling edge it declares, touches either species. The
sum of `M_nad_c` and `M_nadh_c` SHALL therefore be invariant under this
sub-model's dynamics.

The system SHALL check this on a nominal trajectory of this sub-model alone —
not deferred to a later validation change — and the check SHALL fail rather than
warn when the sum drifts beyond integrator tolerance.

Because the pool is closed within this one module, the check is meaningful for
the sub-model in isolation, and its passing says nothing about whether the wider
Core A′ boundary is closed.

#### Scenario: Redox balance over a nominal trajectory
- **WHEN** the sub-model is integrated alone from the registry's initial
  conditions at published parameters
- **THEN** `M_nad_c + M_nadh_c` equals its initial value at every reported time,
  to within integrator tolerance

#### Scenario: A broken redox stoichiometry is caught
- **WHEN** the NAD⁺ or NADH stoichiometry of `GAPD` or `LDH_L` is perturbed so
  that the two reactions no longer mirror each other
- **THEN** the redox check fails, naming the conserved pool and the drift

#### Scenario: No coupling edge touches the redox pool
- **WHEN** the sub-model's declared coupling edges are enumerated
- **THEN** none of them names `M_nad_c` or `M_nadh_c`
