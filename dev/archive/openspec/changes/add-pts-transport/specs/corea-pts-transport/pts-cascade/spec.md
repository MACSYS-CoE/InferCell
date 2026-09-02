## Purpose

The five phosphotransferase steps by which glucose enters Core A′, the eight
phospho-carrier states they own, the four conservation invariants those steps
imply, and the membrane-protein flag on ptsG that the growth coupling reads —
so that the one route carbon takes into the model is defined once and checks
itself.

## ADDED Requirements

### Requirement: The five cascade steps and their stoichiometry

The system SHALL provide a Core A′ sub-model carrying exactly five
phosphotransferase reactions, with the stoichiometry of the published transport
reconstruction:

| Reaction | Stoichiometry |
|---|---|
| GLCpts0 | `M_ptsi_c` + `M_pep_c` ⇌ `M_ptsi_P_c` + `M_pyr_c` |
| GLCpts1 | `M_ptsh_c` + `M_ptsi_P_c` ⇌ `M_ptsh_P_c` + `M_ptsi_c` |
| GLCpts2 | `M_crr_c` + `M_ptsh_P_c` ⇌ `M_crr_P_c` + `M_ptsh_c` |
| GLCpts3 | `M_ptsg_c` + `M_crr_P_c` ⇌ `M_ptsg_P_c` + `M_crr_c` |
| GLCpts4 | `M_glc__D_e` + `M_ptsg_P_c` ⇌ `M_g6p_c` + `M_ptsg_c` |

Every species named SHALL be a Core A′ registry species, referred to by its
registry name; the sub-model SHALL NOT define its own species names. Each
reaction SHALL be reversible, carrying an independent forward and reverse
constant.

The cascade SHALL consume phosphoenolpyruvate and not ATP. This is the mechanism
that makes Core A′'s net yield two ATP per glucose rather than three, so a
formulation that debited ATP here would change the model's energy budget.

#### Scenario: The cascade is exactly five steps
- **WHEN** the sub-model's reactions are enumerated
- **THEN** exactly the five reactions above are present, each exactly once

#### Scenario: Every species is a registry species
- **WHEN** each reaction's substrates and products are collected
- **THEN** every one resolves to a Core A′ registry position

#### Scenario: No adenylate species appears in the cascade
- **WHEN** the cascade's stoichiometry is inspected
- **THEN** no reaction names `M_atp_c`, `M_adp_c`, `M_amp_c` or `M_pi_c`
- **AND** `M_pep_c` is consumed by `GLCpts0`

### Requirement: The eight phospho-states are owned, and no others

The sub-model SHALL integrate the registry's eight `:pts` phospho-states. The
registry species its reactions touch but does not own — `M_pep_c`, `M_pyr_c`,
`M_g6p_c` and `M_glc__D_e` — SHALL NOT be integrated here.

State ordering SHALL follow the registry's canonical order.

#### Scenario: The owned phospho-states are the registry's PTS group
- **WHEN** the sub-model's states are filtered to the registry's `:pts` group
- **THEN** all eight are present, each exactly once

#### Scenario: Shared metabolites are not owned
- **WHEN** the sub-model's states are queried for `M_pep_c`, `M_pyr_c`,
  `M_g6p_c` or `M_glc__D_e`
- **THEN** none of the four is present

### Requirement: Four genes, one per carrier protein

The sub-model SHALL record the gene behind each carrier: ptsI
`JCVISYN3A_0233`, Crr `JCVISYN3A_0234`, ptsH `JCVISYN3A_0694`, ptsG
`JCVISYN3A_0779`. Each gene SHALL map to exactly one carrier and therefore to a
pair of phospho-states.

Unlike the glycolytic reactions, the carriers are the reacting species rather
than enzymes scaling a rate law, so the gene map here identifies which pair of
states a translated protein enters, not which rate it catalyses.

#### Scenario: Each gene maps to one phospho-state pair
- **WHEN** the gene map is queried for a locus
- **THEN** it returns exactly two registry species, the unphosphorylated and
  phosphorylated forms of one carrier
- **AND** the four loci are distinct and cover all eight states

### Requirement: ptsG is flagged as Core A′'s only membrane protein

The sub-model SHALL record that `JCVISYN3A_0779` is a membrane protein and that
both `M_ptsg_c` and `M_ptsg_P_c` count toward membrane surface area, together
with the per-protein footprint the published model uses. This flag SHALL be
retrievable from the composed model without reading this module's source.

No other protein in Core A′ carries this flag, so a later change computing
surface area from protein counts SHALL be able to find every contributing state
by querying for it.

#### Scenario: The membrane flag names both phospho-forms
- **WHEN** the composed model is queried for its membrane-protein states
- **THEN** `M_ptsg_c` and `M_ptsg_P_c` are returned, and no other species
- **AND** the per-protein surface-area footprint is available alongside them

#### Scenario: The flag is discoverable from the composed model
- **WHEN** a module that did not declare the flag queries the composition for it
- **THEN** the ptsG states are found without that module knowing which sub-model
  declared them

### Requirement: Each carrier protein is conserved across its two forms

Every cascade step moves a phosphate between adjacent carriers and creates or
destroys no carrier. The sums `M_ptsi_c + M_ptsi_P_c`, `M_ptsh_c + M_ptsh_P_c`,
`M_crr_c + M_crr_P_c` and `M_ptsg_c + M_ptsg_P_c` SHALL therefore each be
invariant under this sub-model's dynamics.

The system SHALL check all four on a nominal trajectory of this sub-model alone
— not deferred to a later validation change — and SHALL fail rather than warn
when any sum drifts beyond integrator tolerance. The four SHALL be checked
independently, so that a failure names which carrier drifted.

#### Scenario: All four carrier pools are conserved
- **WHEN** the sub-model is integrated alone from its initial conditions at
  published parameters
- **THEN** each of the four sums equals its initial value at every reported time,
  to within integrator tolerance

#### Scenario: A failure names the carrier
- **WHEN** one cascade step's stoichiometry is perturbed so that a carrier is
  created rather than transferred
- **THEN** the check fails naming that carrier and the drift
- **AND** the other three carriers' checks still pass, so the failure is localised

#### Scenario: Conservation holds independently of the shared metabolites
- **WHEN** `M_pep_c` and the other non-owned species are held at any values
- **THEN** the four carrier sums remain invariant, because no cascade step
  couples a carrier's total to a metabolite
