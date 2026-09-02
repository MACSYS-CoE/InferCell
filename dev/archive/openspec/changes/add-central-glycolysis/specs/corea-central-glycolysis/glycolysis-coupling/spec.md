## Purpose

What this module declares at its boundary — typed edges for the registry species
it shares with the other Core A′ modules, untyped inputs for the enzyme counts no
edge kind can name — and the requirement that it resolve, integrate and be tested
entirely on its own, before any of its six wave-1 siblings exists.

## ADDED Requirements

### Requirement: The energy currencies are declared as currency edges

`M_atp_c`, `M_adp_c` and `M_pi_c` cross this module's boundary and are integrated
elsewhere in Core A′. Each SHALL be declared as a currency edge, with the
direction its reactions imply:

| Species | Direction | Reactions |
|---|---|---|
| `M_atp_c` | inbound | PFK draws ATP |
| `M_atp_c` | outbound | PGK and PYK supply ATP |
| `M_adp_c` | inbound | PGK and PYK draw ADP |
| `M_adp_c` | outbound | PFK supplies ADP |
| `M_pi_c` | inbound | GAPD draws phosphate |

Each inbound currency edge names a dynamic state this module does not integrate,
so the species SHALL also appear in the module's untyped inputs — that is the only
channel by which it reaches the dynamics.

No deferred-counter, rate-constant or volume edge is declared here: this module
accrues no stochastic cost, rebuilds no rate constants, and drives no volume.

#### Scenario: The five currency edges are declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** exactly the five currency edges above are present, each naming a
  registry species

#### Scenario: Inbound currencies are reachable by the dynamics
- **WHEN** the module's untyped inputs are enumerated
- **THEN** `M_atp_c`, `M_adp_c` and `M_pi_c` are all present
- **AND** composing the module raises no disagreement between its typed edges and
  its untyped inputs

#### Scenario: No edge kind beyond mass and currency is declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** no deferred-counter, rate-constant, volume or clamped edge appears

### Requirement: Shared owned states are declared as mass edges

Four of the thirteen states this module owns are also touched by another Core A′
module, and each such crossing SHALL be declared as a mass edge from this module's
side, with the direction mass flows relative to this module:

| Species | Direction | Counterpart reaction |
|---|---|---|
| `M_g6p_c` | inbound | the PTS cascade's final step supplies it |
| `M_pyr_c` | inbound | the PTS cascade's first step supplies it |
| `M_pep_c` | outbound | the PTS cascade's first step draws it |
| `M_lac__L_c` | outbound | lactate export draws it |

These edges name states this module integrates, so they SHALL NOT require a
matching untyped input — the state is already local.

#### Scenario: The four mass edges are declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** exactly the four mass edges above are present

#### Scenario: A mass edge on an owned state needs no input
- **WHEN** the module is composed
- **THEN** the inbound mass edges on `M_g6p_c` and `M_pyr_c` raise no requirement
  that either be listed as an untyped input

### Requirement: Enzyme counts enter as untyped inputs

Protein counts are not Core A′ registry species, so no edge kind can name them.
The ten enzyme counts this module reads SHALL therefore be declared through the
untyped input channel, which the coupling contract leaves open for exactly this
case, and SHALL pass composition unchecked against the typed declarations.

#### Scenario: Enzyme inputs are outside the typed contract
- **WHEN** the module is composed
- **THEN** its ten enzyme inputs raise no requirement for a matching coupling edge
- **AND** composition succeeds

#### Scenario: No edge names a protein
- **WHEN** the module's coupling declarations are enumerated
- **THEN** every species named resolves to a registry position

### Requirement: The module resolves, integrates and is tested alone

The module SHALL be composable, resolvable and integrable by itself, with none of
its wave-1 siblings present. No coupling edge SHALL name a peer module, so that
each edge resolves against whichever module owns its species and a single-module
composition does not fail for a counterpart that is deliberately absent.

Resolving the module alone SHALL report — not fail on — the states it consumes
with no declared producer and the registry states no module integrates. A
successful resolution of this module alone is evidence that its own declarations
are internally consistent, and SHALL NOT be read as evidence that the Core A′
boundary is closed.

#### Scenario: The module resolves alone
- **WHEN** the module is composed by itself and its coupling is resolved
- **THEN** resolution succeeds
- **AND** no edge fails for a missing counterpart module

#### Scenario: Standalone resolution reports rather than fails
- **WHEN** the module is resolved alone
- **THEN** `M_atp_c`, `M_adp_c` and `M_pi_c` are reported as unowned states
- **AND** the report distinguishes this deliberate partial composition from an
  incomplete one

#### Scenario: The module integrates alone
- **WHEN** the module is integrated by itself from the registry's initial
  conditions at published parameters, with its inbound currencies held at their
  registry values
- **THEN** a trajectory is produced over the requested interval
- **AND** every one of the thirteen owned states remains non-negative

### Requirement: Declared production of a state this module does not own is not executed

This module declares outbound currency edges for the ATP and ADP its reactions
produce, but does not integrate either species, and the composition machinery
gives a module no way to contribute to a derivative it does not own. Those
outbound edges are therefore a declaration of intent that a later change executes.

The module SHALL NOT present its inbound currencies as though the loop were
closed: a trajectory produced by this module alone or by any composition that does
not itself execute the coupling SHALL be labelled as one in which `M_atp_c`,
`M_adp_c` and `M_pi_c` are held rather than integrated.

#### Scenario: Held currencies are labelled as held
- **WHEN** a trajectory is produced from a composition that does not execute the
  currency coupling
- **THEN** the result records that the energy currencies were held at fixed values
- **AND** the record names the three species

#### Scenario: The declaration survives for the later change to consume
- **WHEN** the module's coupling declarations are enumerated
- **THEN** the outbound currency edges on `M_atp_c` and `M_adp_c` are present,
  whether or not any composition executes them
