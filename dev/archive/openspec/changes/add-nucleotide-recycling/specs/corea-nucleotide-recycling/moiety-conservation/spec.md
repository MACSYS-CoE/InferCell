## Purpose

The three numerical checks that say this module does the job it exists for:
adenylate and guanylate conserved over a full cell cycle rather than a burst,
phosphate closed with pyrophosphate counted as two, and the dead end of record
reproduced and then removed — so that "closes a moiety" is a measurement rather
than a claim.

## ADDED Requirements

### Requirement: Adenylate and guanylate are conserved, and the check runs over a full cycle

The sums `M_atp_c + M_adp_c + M_amp_c` and `M_gtp_c + M_gdp_c + M_gmp_c` SHALL
each be invariant under this sub-model's own dynamics. Every reaction it carries
preserves both: `ADK1` moves phosphate within the adenylate pool, `GK1` moves one
between the two pools without changing either total, and `PGK3` and `PYK3`
phosphorylate GDP to GTP.

The check SHALL be run over a full cell cycle rather than a short interval. The
failure mode it exists to catch is a slow one-way drain, and the adenylate dead
end of record looked healthy for the first 144 seconds of a 6,300-second
trajectory. A conservation check that runs for less than a cycle would have passed
on the model that was wrong.

The two sums SHALL be checked independently, so a failure names which moiety
drifted, and the tolerance SHALL be derived from the integrator's settings rather
than chosen by hand.

#### Scenario: Both moieties survive a full cycle
- **WHEN** the sub-model is integrated alone over a full cell cycle from the
  registry's initial conditions at published parameters
- **THEN** each sum equals its initial value at every reported time, within
  integrator tolerance

#### Scenario: A short run is not sufficient evidence
- **WHEN** the conservation check is configured
- **THEN** its interval is at least a full cell cycle
- **AND** the recorded reason names the 144-second-in-6,300 failure it exists to
  catch

#### Scenario: A failure names the moiety
- **WHEN** a reaction's stoichiometry is perturbed so that one pool gains a
  species the other loses
- **THEN** the check fails naming that moiety
- **AND** the other moiety's check still passes

### Requirement: Phosphate is closed, with pyrophosphate counted as two

The system SHALL check phosphate closure over this sub-model's own dynamics:
inorganic phosphate plus the phosphate carried by every phosphorylated species it
owns, with `M_ppi_c` contributing two.

Three of the five reactions close phosphate entirely within the owned states, and
the check SHALL hold exactly for a composition exercising only those. `PGK3` and
`PYK3` carry phosphate across the module boundary from species this module does
not own, so with those two active the check SHALL account for what crosses,
rather than being relaxed or skipped.

This is the check that catches a missing pyrophosphatase.

#### Scenario: Closure holds exactly without the boundary-crossing reactions
- **WHEN** the sub-model is integrated with `PGK3` and `PYK3` inactive
- **THEN** the total phosphate over the owned states is invariant within
  integrator tolerance
- **AND** pyrophosphate contributes two phosphates per molecule throughout

#### Scenario: Closure accounts for the boundary flux when all five are active
- **WHEN** the sub-model is integrated with all five reactions active and its
  inbound metabolites held
- **THEN** the total phosphate over the owned states, less the phosphate delivered
  across the declared inbound edges, is invariant within tolerance

#### Scenario: Removing pyrophosphatase breaks closure
- **WHEN** `PPA` is removed and pyrophosphate is supplied
- **THEN** the phosphate check fails, and pyrophosphate accumulates rather than
  returning to inorganic phosphate

### Requirement: Removing ADK1 reproduces the dead end of record

The system SHALL demonstrate, numerically, that this module closes the adenylate
cycle. Against a consumer draining adenosine triphosphate to monophosphate and
pyrophosphate at the rate the lumped charging step implies — roughly 553 per
second — the module SHALL keep the adenylate pool intact and the triphosphate pool
positive across a full cycle.

With `ADK1` removed under the same drain, the triphosphate and diphosphate pools
SHALL be exhausted within a few minutes, reproducing the failure the earlier
specification of this model contained.

This is what distinguishes a module that closes a moiety from one that merely
declares reactions that could.

#### Scenario: The pool survives the charging drain
- **WHEN** the sub-model is integrated over a full cycle against a consumer
  draining at the lumped charging rate
- **THEN** the adenylate sum is conserved
- **AND** `M_atp_c` remains positive throughout

#### Scenario: Without ADK1 the pool is exhausted
- **WHEN** the same integration is run with `ADK1` removed
- **THEN** the triphosphate and diphosphate pools reach zero well within the cycle
- **AND** the time at which they do is of the order of the recorded 144 seconds

#### Scenario: Without PPA pyrophosphate accumulates
- **WHEN** the same integration is run with `PPA` removed
- **THEN** pyrophosphate rises without bound rather than settling
- **AND** the magnitude it reaches is of the order recorded for the unhydrolysed
  pool
