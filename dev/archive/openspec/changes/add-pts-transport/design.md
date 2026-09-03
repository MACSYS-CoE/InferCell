## Context

See `proposal.md` — Why. The constraints that shape the approach, beyond the ones
`add-central-glycolysis` already records (frozen interface, no network on compute
nodes, six sibling branches open, a module writes only its own derivatives):

- **This module's parameters have no priors at all.** The transport source file
  carries four SBtab tables — `Compartment`, `Compound`, `Reaction`, `Quantity` —
  and no `Parameter` table. There is no `Mode` column, no geometric standard
  deviation, nothing to inherit. Eleven point values.
- **The states this module owns are proteins, not metabolites.** They sit in the
  registry's `:protein` copy-number regime at 290–831 copies, two orders of
  magnitude below the metabolite pools, which is what makes them the states a
  later change reads for surface area.
- **Two of its four owned carriers' initial conditions are not in any table.**
  The registry records all eight `:pts` states as `not_imported`, and the two
  candidate sources in the model disagree.

## Goals / Non-Goals

**Goals**

- One sub-model carrying six reactions and nine states, constructible, integrable
  and testable alone.
- The eleven asserted priors visible as asserted, since this module is where that
  distinction first has teeth.
- Four independent protein-conservation invariants checked here, numerically.
- The ptsG membrane flag discoverable by a module that did not declare it.

**Non-Goals**

- Growth. This module owns the states surface area is computed from and flags
  them; computing area, radius or volume is `add-growth-volume-coupling`'s.
- Executing the metabolite coupling. `M_g6p_c` and `M_pyr_c` production is
  declared and left unexecuted, as in `add-central-glycolysis`.
- Choosing priors for the eleven constants beyond what makes them samplable. They
  are marked asserted; how wide to assert is a question for whichever change first
  frees one.
- Carbon balance. It spans modules and belongs to wave 3 — but note that this
  module is the one that makes it *possible* to pose.

## Decisions

### D1 — Mass-action for the cascade, and the reason it differs from glycolysis

The five cascade steps take reversible second-order mass action, written exactly
as `transport_NoH2O_Zane-TB-DB.tsv:126-130` has them:

```
v = KF · A · B  −  KR · C · D
```

No saturation, no enzyme multiplier — the carriers *are* the reactants. This is
genuinely a different rate law from `add-central-glycolysis`'s modular form, and
it is the source model's own split: intracellular reactions get convenience
kinetics, transport gets mass action.

The consequence worth stating: the cascade has no Michaelis constants, so nothing
here bounds a rate as concentrations rise. With carriers conserved that is
harmless — the pools cap themselves — but it is why the conservation check is the
right invariant to lean on rather than a flux ceiling.

### D2 — Eleven constants, all asserted, and the table shape that guarantees it

| Identifier | Value | Identifier | Value |
|---|---|---|---|
| `KF_0_R_GLCpts0` | 6600 | `KR_0_R_GLCpts0` | 4000 |
| `KF_1_R_GLCpts1` | 200000 | `KR_1_R_GLCpts1` | 8000 |
| `KF_2_R_GLCpts2` | 61000 | `KR_2_R_GLCpts2` | 47000 |
| `KF_3_R_GLCpts3` | 3900 | `KR_3_R_GLCpts3` | 310 |
| `KF_4_R_GLCpts4` | 0.88 | `KR_4_R_GLCpts4` | 1.00e-05 |
| `P_R_L_LACt2r` | 5.00e-09 m/s | | |

The extract `src/organisms/coreA/data/pts_transport.tsv` deliberately carries
**no `GeometricStd` column**. `read_source_table` derives informedness from what
it finds, and "a row with no uncertainty column at all is `:asserted` — any prior
on it is this project's, not the model's". So the absence of the column is what
makes all eleven load as asserted; nothing has to be declared per row, and no
future edit can quietly promote one to `:balanced` without adding a column that
does not exist upstream.

Priors are `LogNormal(log(value), log(g))` with a single explicit `asserted_gstd`
constructor keyword rather than eleven separate numbers, so that the one thing we
are inventing is one number and is visible as one number. Every parameter is
`fixed = true` by default, freed by a `free` keyword, as in
`add-central-glycolysis` D5.

*Alternative considered:* an explicit `Informedness` column set to `asserted`.
Rejected as redundant and less robust — it states in data what the table's shape
already guarantees, and it could drift from the shape.

### D3 — Initial conditions from copy number × proteomics fraction

Totals are the published per-cell copy numbers; the split is the published
model's `proteomics_fraction`, read at startup by `MinCell_CMEODE.py:199-203`:

| Carrier | Gene | Copies | Unphos. frac | Unphos. (mM) | Phos. (mM) |
|---|---|---|---|---|---|
| ptsI | JCVISYN3A_0233 | 353 | 0.05 | 0.0008746 | 0.0166179 |
| Crr | JCVISYN3A_0234 | 314 | 0.05 | 0.0007780 | 0.0147820 |
| ptsH | JCVISYN3A_0694 | 290 | 0.05 | 0.0007185 | 0.0136521 |
| ptsG | JCVISYN3A_0779 | 831 | 0.15 | 0.0061769 | 0.0350025 |

`M_lac__L_e` starts at 0.0, from `setICs_two.py:277`.

**Not** the transport reconstruction's uniform 0.024 mM per form. That value puts
969 copies behind every carrier — 2.7×, 3.1× and 3.3× the proteomics counts for
ptsI, Crr and ptsH — so it is a placeholder, not a measurement. The extract
records it in a comment alongside the derived rows, so the discrepancy is on the
page rather than rediscovered.

Both inputs are published, so the split is `:published`, not asserted; only the
arithmetic is ours, and 20,180 particles per mM at r = 200 nm is the registry's
own cell. `dev/scripts/extract_pts_transport.jl` does the arithmetic and records
both inputs in each row's `UpstreamRow` field, so a derived concentration traces
back to the copy number and the fraction it came from.

*Alternative considered:* the commented-out block at `setICs_two.py:286-288`,
which carries the same fractions against rounded totals. It agrees to four
decimals but is disabled code; the CSVs are live model data and are the better
citation.

### D4 — Two extract files, because provenance is per file

`src/organisms/coreA/data/pts_transport.tsv` — logical name `"transport"`, 11
rows, no uncertainty column.
`src/organisms/coreA/data/pts_initial_conditions.tsv` — logical name
`"model_ics"`, 9 rows.

Two `SourceTable`s rather than one, because `ParameterSource` records one file per
table and a single extract would report `"transport"` as the source of values that
came from `setICs_two.py` and two CSVs. The identifier sets are disjoint, so
`ambiguity_report` over the pair is empty — and that emptiness is asserted, since
the wave-0 contract is explicit that an empty report is a positive result rather
than a failure to run.

The `conc_<species>` prefix is used for the nine initial conditions even though
the registry records no value for any of them, so `_check_registry_agreement`
returns early. It costs nothing now and means the check fires automatically if a
later change gives the registry a value.

### D5 — The lactate volume ratio, and the number behind the default

`d[M_lac__L_e]/dt = +v_export / R`, with `R = V_medium / V_cell` and

```
v_export = P · (lac_c − lac_e) · 3 / r_cell
         = 0.075 · (lac_c − lac_e)  s⁻¹   at P = 5e-9 m/s, r = 200 nm
```

The default `R = 1e5` is chosen against a stated criterion, not picked. At Core
A′'s ~1,106 glucose/s, lactate production is 2,212/s = 0.10961 mM/s, so steady
cytosolic lactate is 0.10961 / 0.075 ≈ 1.46 mM. One cycle exports 13.9 M lactate,
688.8 mM at initial volume:

| R | external lactate after a cycle | as a fraction of steady cytosolic |
|---|---|---|
| 1e3 | 0.689 mM | 47% |
| 1e4 | 0.0689 mM | 4.7% |
| **1e5** | **0.00689 mM** | **0.47%** |

The criterion is "below one percent over a full cycle", which 1e5 meets with
margin and 1e4 does not. `R = 1` is representable and is what the shared-volume
alternative would be; the spec requires a test that it visibly saturates, so the
default's purpose is demonstrated rather than asserted.

`R` is exposed, marked asserted, and registered in `reduction_notes` naming both
what it reproduces (the published constant external pool) and what it replaces
(a genuinely dynamic state, which the registry requires).

*Alternative considered:* clamping `M_lac__L_e` at 0 with `origin = :published`,
which is literally what the source does. Rejected because the registry lists the
species among its 32 dynamic states and the resolver rejects a clamp on a state a
module integrates — changing that is an interface bug to raise on `main`, not an
edit on this branch.

### D6 — The cell radius is a named quantity, not folded into 0.075

`3P/r` collapses to a single constant at r = 200 nm, and folding it would be
tempting and wrong: `add-growth-volume-coupling` makes the radius grow, at which
point the export rate constant changes with it. `P` and `r_cell` are carried
separately, `r_cell` defaults to the registry's cell, and the export rate is
computed from both. That is ten characters of arithmetic per step against a
silent bug in wave 2.

### D7 — The ptsG membrane flag is data on the module, discoverable by others

`membrane_protein_states(m)` returns `(:M_ptsg_c, :M_ptsg_P_c)` with the 28.0 nm²
footprint, defaulting to empty on `AbstractSubModel` — the same shape
`reduction_notes` and `coupling` take. A growth change then sweeps the
composition rather than knowing that `PtsTransport` is where to look.

The footprint is the published model's calibrated constant, not a measurement:
28.0 nm² was chosen to reproduce 54% membrane coverage for ~9,600 membrane
proteins. It is imported with that recorded, so wave 2 does not treat it as
measured. Note `getProtSA`'s docstring says 35 nm² while its code uses 28.0; the
code governs, and the discrepancy is recorded in the extract's `README.md`.

*Alternative considered:* leaving wave 2 to hard-code the two state names.
Rejected: it puts organism knowledge in a framework-level change and silently
breaks if a second membrane protein is ever added.

### D8 — Edges, directions, and what is deliberately absent

| Species | Kind | Dir | Owned | Why |
|---|---|---|---|---|
| `M_glc__D_e` | clamped | in | no (chemostat) | GLCpts4 draws, held at 40 mM, `origin = :published` |
| `M_pep_c` | mass | in | no | GLCpts0 draws |
| `M_lac__L_c` | mass | in | no | L_LACt2r draws |
| `M_pyr_c` | mass | out | no | GLCpts0 supplies |
| `M_g6p_c` | mass | out | no | GLCpts4 supplies |

`inputs` holds exactly two names, `M_pep_c` and `M_lac__L_c` — the two inbound
mass edges on states this module does not integrate. `M_glc__D_e` is deliberately
**not** in `inputs`: the resolver rejects a chemostatted input, and its error text
prescribes this exact shape, "declare a `ClampedEdge` with its `held_value` and
carry the value as a fixed parameter instead".

Unlike `add-central-glycolysis`, there are no protein-count inputs. The carriers
are states this module owns, not enzymes scaling a rate law, so translation
feeding them is an inbound mass edge on an owned state — a wave-2 crossing, and
one this module does not have to declare from its side to be correct today.

Every edge leaves `peer` unnamed, for the reason `add-central-glycolysis` D7
records.

The held value is 40 mM, the registry's. The transport reconstruction says 42.77
mM; `setICs_two.py:279` overrides it with 40, and what runs governs. Declaring
42.77 would in any case throw, since `ClampedEdge` validates a chemostat's held
value against the registry.

### D9 — Protein conservation is four checks, not one

Each cascade step transfers a phosphate between adjacent carriers, so each of the
four sums is invariant independently. Summing all four into one number would let
two errors of opposite sign cancel, so the test asserts four bounds and names the
carrier that drifts. A mutation test perturbing one step's stoichiometry confirms
the check fails and that the other three still pass, so a failure is localised
rather than global.

Tolerance is derived from the integrator's `abstol`/`reltol`, as in
`add-central-glycolysis` D10, so tightening the solver tightens the check.

## Risks / Trade-offs

- **Eleven parameters with priors we invent.** This module contributes every
  asserted prior in wave 1, and an asserted prior that is too tight is
  indistinguishable, in a posterior, from information. → All eleven are marked
  `:asserted` and enumerable; `asserted_gstd` is one visible number rather than
  eleven buried ones; and all are `fixed = true` by default, so nothing is
  sampled from an invented prior without an explicit act.

- **The volume ratio could become load-bearing.** A carbon-balance result computed
  with `R = 1e5` depends on a number we chose. → It is registered in
  `reduction_notes`, the criterion behind the default is written down, and the
  shared-volume case is exercised by a test so the alternative is visible.

- **The cascade is stiff and fast.** Forward constants span 0.88 to 200,000, and
  the fastest step turns over its pool in microseconds while lactate export runs
  at 0.075 s⁻¹ — seven orders of magnitude. → A stiff solver with explicit
  tolerances, pinned in this module's tests. The fast cascade is expected to reach
  a quasi-steady phospho-distribution almost immediately, which the conservation
  check is insensitive to by construction.

- **A standalone trajectory is not physiological.** With PEP and cytosolic lactate
  held, the cascade drives to a fixed point and export runs at a constant rate.
  → The acceptance criteria are non-negativity and the four conservation
  invariants, not that the trajectory looks biological, and the held inputs are
  labelled.

- **Wave 2 may want the carriers to grow.** Translation adds new PTS protein, at
  which point the four sums are no longer invariant and this module's headline
  check becomes wrong as stated. → The check is scoped in the spec to this
  sub-model's own dynamics. Wave 2 inheriting it needs the invariant restated as
  "conserved up to what translation adds", and this is flagged in the handoff
  rather than discovered.

## Migration Plan

Additive. One `include` in `src/InferCell.jl` after the registry, one in
`test/runtests.jl`. `membrane_protein_states` ships as a plain function over
`PtsTransport` in this module's own file, not as a protocol function with a
default on `AbstractSubModel` — that would mean editing `src/interface.jl`, which
is read-only to a wave-1 branch. Promoting it is wave 2's call; see Open
Questions. No existing sub-model, test or exported name changes, so rollback is
reverting the branch.

## Open Questions

- **Where `membrane_protein_states` is declared.** D7 wants it as a protocol
  function with an empty default, which means `src/interface.jl` — read-only here.
  The alternatives are declaring it in this module's own file as a plain function
  over `PtsTransport` only (works today; wave 2 then needs a `try`/fallback or a
  type check), or raising it on `main` as a small interface addition before this
  change lands. It changes nothing about the reactions, the states or the checks,
  and only wave 2 consumes it, so it can be settled when wave 2 starts. **Until
  then, implement it as a plain function on `PtsTransport`.**

- **How wide the asserted prior should be.** `asserted_gstd` needs a default to be
  constructible; whether that default is defensible matters only once one of the
  eleven is freed, which no change yet does. Worth deciding against a sensitivity
  run rather than by argument.
