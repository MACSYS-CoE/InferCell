## Purpose

What the owner of Core A′'s energy pools declares at its boundary: the four
glycolytic intermediates its own reactions reach for, and the five pools where it
is the model's only producer or only consumer — the return paths whose absence was
the error of record — together with the requirement that it resolve alone.

## ADDED Requirements

### Requirement: The glycolytic intermediates are declared as mass edges

The four registry metabolites this module's reactions touch but does not
integrate SHALL each be declared as a mass edge, with the direction mass flows
relative to this module:

| Species | Direction | Reaction |
|---|---|---|
| `M_13dpg_c` | inbound | PGK3 draws |
| `M_pep_c` | inbound | PYK3 draws |
| `M_3pg_c` | outbound | PGK3 supplies |
| `M_pyr_c` | outbound | PYK3 supplies |

Each inbound edge names a dynamic state this module does not integrate, so that
species SHALL also appear in the module's untyped inputs.

#### Scenario: The four mass edges are declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** exactly the four mass edges above are present, each naming a registry
  species

#### Scenario: Inbound metabolites are reachable by the dynamics
- **WHEN** the module's untyped inputs are enumerated
- **THEN** `M_13dpg_c` and `M_pep_c` are present, and no other registry species is

### Requirement: The module declares the return path it is the sole holder of

Where this module is Core A′'s only producer or only consumer of a pool it owns,
it SHALL declare a currency edge naming that species and the direction its own
reactions move it:

| Species | Direction | Reaction | What the edge asserts |
|---|---|---|---|
| `M_amp_c` | inbound | ADK1 draws | the only route adenosine monophosphate has back |
| `M_gmp_c` | inbound | GK1 draws | the only route guanosine monophosphate has back |
| `M_ppi_c` | inbound | PPA draws | the only route pyrophosphate has out |
| `M_pi_c` | outbound | PPA supplies | phosphate returned from pyrophosphate |
| `M_gtp_c` | outbound | PGK3, PYK3 supply | the only source of guanosine triphosphate |

These name states the module integrates, so they SHALL NOT require a matching
untyped input.

The module SHALL NOT declare an edge for a pool it owns where another module is
the counterpart on both sides — triphosphate and diphosphate traffic with
glycolysis and with the stochastic blocks is declared by those modules, and
declaring it again here would state the same crossing twice.

#### Scenario: The five currency edges are declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** exactly the five currency edges above are present

#### Scenario: A currency edge on an owned state needs no input
- **WHEN** the module is composed
- **THEN** the inbound currency edges raise no requirement that their species be
  listed as untyped inputs

#### Scenario: The declarations do not conflict with a sibling's
- **WHEN** this module's declarations are composed with those of the modules
  owning the glycolytic and transport reactions
- **THEN** no species and direction is described as both mass and currency
- **AND** resolution succeeds

#### Scenario: No other edge kind is declared
- **WHEN** the module's coupling declarations are enumerated
- **THEN** no deferred-counter, catalytic, rate-constant, volume or clamped edge
  appears

### Requirement: The module resolves, integrates and is tested alone

The module SHALL be composable, resolvable and integrable by itself. No coupling
edge SHALL name a peer module, so each resolves against whichever module owns its
species and a single-module composition does not fail for a counterpart
deliberately absent.

Resolving alone SHALL report — not fail on — states consumed with no declared
producer and registry states no module integrates. A successful resolution SHALL
NOT be read as evidence that the Core A′ boundary is closed; this module owns the
pools that make that distinction matter, since a cost with no paying state is
exactly what its reports would show and exactly what wave 3 must fail on.

#### Scenario: The module resolves alone
- **WHEN** the module is composed by itself and its coupling is resolved
- **THEN** resolution succeeds, and no edge fails for a missing counterpart

#### Scenario: Standalone resolution reports rather than fails
- **WHEN** the module is resolved alone
- **THEN** `M_13dpg_c` and `M_pep_c` are reported as unowned states
- **AND** the report distinguishes this deliberate partial composition from an
  incomplete one

#### Scenario: The module integrates alone
- **WHEN** the module is integrated by itself over a full cell cycle from the
  registry's initial conditions at published parameters, with its inbound
  metabolites held
- **THEN** a trajectory is produced
- **AND** all eight owned states remain non-negative

### Requirement: Declared production of a state this module does not own is not executed

This module declares outbound mass edges for the 3-phosphoglycerate and pyruvate
its reactions produce, but integrates neither, and the composition machinery gives
a module no way to contribute to a derivative it does not own. Those edges are a
declaration of intent that a later change executes.

A trajectory produced by this module alone, or by a composition that does not
execute the coupling, SHALL be labelled as one in which the inbound metabolites
were held rather than integrated.

#### Scenario: Held metabolites are labelled as held
- **WHEN** a trajectory is produced from a composition that does not execute the
  metabolite coupling
- **THEN** the result records that `M_13dpg_c` and `M_pep_c` were held at fixed
  values

#### Scenario: The declaration survives for the later change to consume
- **WHEN** the module's coupling declarations are enumerated
- **THEN** the outbound mass edges on `M_3pg_c` and `M_pyr_c` are present,
  whether or not any composition executes them
