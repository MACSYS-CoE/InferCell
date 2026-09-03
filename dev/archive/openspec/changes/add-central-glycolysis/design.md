## Context

See `proposal.md` — Why. The constraints that actually shape the approach:

- **The interface is frozen.** `src/edges.jl`, `src/resolver.jl`, `src/loader.jl`,
  `src/labels.jl`, `src/interface.jl`, `src/orchestrator.jl` and
  `src/organisms/coreA/registry.jl` are read-only on this branch. Every design
  choice below has to fit them as they are.
- **A module can only write its own derivatives.** `_build_rhs`
  (`src/orchestrator.jl:203-215`) slices `u` per sub-model, calls `dynamics` on
  the slice, and concatenates the results. There is no channel by which this
  module contributes to `d[M_atp_c]/dt`, which `add-nucleotide-recycling` owns.
- **Compute nodes have no network** and the Julia depot is periodically wiped, so
  anything the test suite needs must be in the repository and parseable with
  Base only.
- **The source model is on disk** at `/fred/oz022/tkimpson/Minimal_Cell`, at
  commit `db048ac` — the commit the registry already cites. It is a checkout, not
  a dependency: nothing at run time may require it.
- **Six sibling branches are open at the same time.** Anything this branch adds
  to a shared file is a merge conflict six times over.

## Goals / Non-Goals

**Goals**

- One sub-model type carrying ten reactions and thirteen states, constructible,
  integrable and testable with nothing else present.
- Every imported number reaching the model through `load_parameter`, so that its
  file, identifier, uncertainty and informedness survive into the composed model.
- The redox invariant checked here, numerically, on a nominal trajectory.
- The one place this port departs from what the published simulator runs stated
  in the code, not only in this document.

**Non-Goals**

- Executing the energy coupling. The outbound currency edges are declared and
  left unexecuted; `add-hook-1s-coupling` closes them.
- Carbon balance, adenylate conservation, phosphate closure. Those span modules
  and belong to wave 3, or to the sibling that owns the moiety.
- Choosing an integrator or a tolerance for the assembled Core A′ model. This
  change picks what its own tests need and no more.
- Deciding which parameters to infer. The module exposes all of them; the choice
  of the first three to six is a later change's.

## Decisions

### D1 — The rate law is `Rxns.Enzymatic`, not the SBtab `KineticLaw` string

The source file carries a `KineticLaw` column, and for `R_PGI` it reads a
`kcrg`/`keq`/`hco` modular form. **That column is not what runs.**
`defMetRxns.py:1538` builds every central-metabolism rate with
`Rxns.Enzymatic(substrates, products)` and substitutes `kcatF`, `kcatR` and the
per-species `Km`; the `KineticLaw` string is inert for these reactions. So the
port takes the Python form:

```
v = onoff · E · ( kcatF·Π(Sᵢ/KmSᵢ) − kcatR·Π(Pⱼ/KmPⱼ) )
                / ( Π(1+Sᵢ/KmSᵢ) + Π(1+Pⱼ/KmPⱼ) − 1 )
```

`onoff` is 1 throughout Core A′ and is dropped. Substrate/product term counts per
reaction: PGI (1,1), PFK (2,2), FBA (1,2), TPI (1,1), GAPD (3,2), PGK (2,2),
PGM (1,1), ENO (1,1), PYK (2,2), LDH_L (2,2) — 32 Michaelis constants in total,
which is the count the import must produce.

*Alternative considered:* the SBtab `KineticLaw` form, which would let `keq`
constrain the reverse rate by a Haldane relation instead of carrying an
independent `kcatR`. Rejected because it is not the published model's behaviour,
and Core A′ is a port.

### D2 — Michaelis constants come from the balanced `Parameter` table

Decided in `proposal.md` — "Two decisions of record". The design consequence is
that this is a **labelled deviation**, which means three things in code: the
eight differing constants are named in `reduction_notes`, the note text carries
both values, and the note is phrased so that `reduction_declarations` output
reads as a sentence a person can act on.

*Alternatives considered:* the `Quantity` table's `Value` column (what the
simulator runs, but carries no uncertainty at all, so eight of the module's
parameters would become asserted priors — and they are the loose ones, exactly
the ones inference would want to target); or `Quantity` values with balanced
`gstd` as prior width (runs what runs, keeps the widths, but leaves eight priors
whose median is not their point value, which is a subtler thing to explain in
every downstream result than a column choice).

### D3 — Kinetics arrive as a vendored, derived extract

`read_source_table` assumes one ID column per row. The balanced `Parameter` table
keys each row on `(QuantityType, Reaction id, Compound id)` and has no ID column
at all. Since `src/loader.jl` is read-only here, the table must be reshaped
before the loader sees it.

`dev/scripts/extract_central_glycolysis.jl` reads the upstream file from a
caller-supplied checkout path and writes `src/organisms/coreA/data/central_glycolysis.tsv`:

```
!!SBtab TableType='Quantity' TableName='central balanced (Core A′ glycolysis extract)'
!ID	!Mode	!GeometricStd	!UpstreamRow
kcatF_R_PGI	804.3384	1.0513	substrate catalytic rate constant|R_PGI|
kcatR_R_PGI	650	1.0513	product catalytic rate constant|R_PGI|
km_R_PGI_M_g6p_c	22.9419	1.0512	Michaelis constant|R_PGI|M_g6p_c
…
conc_M_g6p_c	3.7076	1.2785	concentration||M_g6p_c
```

Sixty-five rows: 10 `kcatF`, 10 `kcatR`, 32 `km`, 13 `conc`. `Mode` and
`GeometricStd` are the loader's default column names, so no keyword arguments are
needed at the call site. `UpstreamRow` is not read by the loader; it exists so a
renamed identifier stays traceable to the row it came from, which is what makes
the reshape auditable rather than a retyping.

The `conc_<species>` prefix is deliberate: it is what fires
`_check_registry_agreement`, so the thirteen initial conditions are checked
against `registry.jl` on every load. Two copies of one number is the trap this
project keeps hitting; here the two copies check each other.

Provenance records the file as `"central_balanced"` — the logical name the
registry already uses — so the extract and the registry rows agree on where a
value came from, and a later change vendoring the nucleotide file produces
genuine two-file ambiguities the loader's `governing` machinery can resolve.

*Alternatives considered:* vendoring the 5,000-row upstream file verbatim
(rejected: it does not fit the loader, and fixing that means editing a read-only
file); reading a configured checkout path at construction (rejected: CI and every
fresh clone would need the checkout, and the suite would exercise synthetic
fixtures rather than the real numbers).

`src/organisms/coreA/data/README.md` records the upstream filename, the commit,
what the reshape does and does not change, and the command that regenerates it.

### D4 — Priors are log-normal from the balanced mode and geometric standard deviation

Every imported value becomes `LogNormal(log(mode), log(gstd))` — the balancing
distribution's own shape, so the prior is inherited rather than asserted. Three
of the thirteen initial conditions (`M_g3p_c` and `M_lac__L_c` here) sit at
`gstd = 10`, the prior default; the loader marks them uninformed and the width is
kept rather than narrowed.

The band this yields is the reason Core A′ is the right place to test recovery:
forward `kcat` gstd runs from 1.0512 (LDH_L) to 1.7466 (FBA), against reverse
constants as loose as 33.03 (PFK). FBA's forward constant being the loosest in
the core is a cross-check on the extract: it is the number the scoping note
quotes.

### D5 — Parameters default to fixed; freeing one is explicit

Sixty-five parameters, all free, is not an inference problem anyone wants to
start with, and the scoping note is explicit that three to six is the first
target. So every `InferParameter` this module constructs carries `fixed = true`,
and the constructor takes a `free` keyword naming the ones to release:

```julia
CentralGlycolysis(; free = [:kcatF_R_FBA])
```

The module therefore exposes the whole parameter vector — nothing is hidden — but
composing it does not silently produce a 65-dimensional posterior.

*Alternative considered:* `fixed = false` for rates, as `LightMetabolism` and
`TierBMetabolism` do. Rejected on scale: those models carry five rates, this one
carries fifty-two.

### D6 — Enzyme concentration is nominal-from-copy-number, and replaceable

Each rate is scaled by an enzyme concentration held in the struct, initialised
from the published copy number at the model's initial volume (20,180 particles
per mM at r = 200 nm, the registry's cell):

| | PGI | PFK | FBA | TPI | GAPD | PGK | PGM | ENO | PYK | LDH_L |
|---|---|---|---|---|---|---|---|---|---|---|
| copies | 266 | 458 | 775 | 410 | 1355 | 411 | 322 | 998 | 551 | 1100 |
| mM | 0.013181 | 0.022696 | 0.038404 | 0.020317 | 0.067146 | 0.020367 | 0.015956 | 0.049455 | 0.027304 | 0.054509 |

Not `defaultPtnConcentration = 0.001 mM`: that is the published model's value for
reactions with **no** GPR rule, which are decoupled from the boundary by
construction. All ten reactions here have single-gene GPRs, so using it would
understate every flux by 13× to 67×.

The ten protein counts are also declared in `inputs`, so that when
`add-corea-translation` and the wave-2 hook exist they overwrite these
concentrations exactly as `createGeneExpression` does. Until then the inputs
resolve against nothing and the nominal values stand. The values are marked as
nominal, not measured — they are a stand-in for a live count, not an imported
concentration.

### D7 — Every edge leaves `peer` unnamed

The resolver fails an edge whose named peer is absent from the composition, and
its own error text says an unnamed peer "resolves against whichever module owns
the species, which is what a single-module composition under test should use".
Six of the seven counterpart modules do not exist yet, so naming peers would make
this module unresolvable on its own branch — the one thing wave 1's decomposition
depends on.

The cost is that a genuinely missing counterpart is not caught at composition.
That is by design: wave 3 owns the completeness assertion, and the wave-0 spec is
explicit that a successful resolve is not evidence of a closed boundary.

The intended counterparts are recorded in prose next to each edge —
`add-nucleotide-recycling` for the three currencies, `add-pts-transport` for the
four shared mass edges — so the intent is on the page even though the type does
not carry it.

### D8 — Directions, and which edges exist at all

Five currency edges and four mass edges, and nothing else. The rule generating
them: **declare an edge exactly where one of this module's reactions touches a
registry species another Core A′ module also touches, with the direction mass
flows relative to this module.**

| Species | Kind | Dir | Owned here | Why |
|---|---|---|---|---|
| `M_atp_c` | currency | in | no | PFK draws |
| `M_atp_c` | currency | out | no | PGK, PYK supply |
| `M_adp_c` | currency | in | no | PGK, PYK draw |
| `M_adp_c` | currency | out | no | PFK supplies |
| `M_pi_c` | currency | in | no | GAPD draws |
| `M_g6p_c` | mass | in | yes | PTS supplies |
| `M_pyr_c` | mass | in | yes | PTS supplies |
| `M_pep_c` | mass | out | yes | PTS draws |
| `M_lac__L_c` | mass | out | yes | export draws |

`M_pyr_c` is the one that needs stating: this module both produces pyruvate (PYK)
and consumes it (LDH_L) internally, so the edge does not describe its own
reactions — it describes the crossing, and the crossing is PTS feeding pyruvate
in. Direction follows mass across the boundary, not net flux inside it.

`inputs` therefore holds thirteen names: the three currencies (which
`_check_inputs_consistency` requires, because they are inbound edges on states
this module does not integrate) and the ten protein counts (which no edge kind
can name, and which the same check skips because they are not registry species).
The four inbound-and-owned mass edges need no input — `e.species in own &&
continue` in the resolver — and adding them to `inputs` would be wrong, since
`inputs` resolves against *other* modules' states.

Currency rather than mass for the three energy species because that is what the
published model does: ATP and GTP traffic is routed through a shared pool node
rather than module to module, and the edge kind exists to say so. `pool` defaults
to the species itself, which is right for both.

### D9 — Standalone integration holds the currencies, and says so

`inputs` resolves against states some module in the composition integrates, so
integrating this module alone needs something that owns `M_atp_c`, `M_adp_c` and
`M_pi_c`. The test file defines a local `HeldEnergyPool <: AbstractSubModel` that
owns exactly those three with `du = 0` — a two-line double, not a model of
anything.

It goes in `test/test_corea_central_glycolysis.jl` rather than the shared
`test/corea_test_models.jl`. Three sibling branches will want the same double,
and six branches editing one shared file is the merge conflict the wave plan
exists to avoid. If three copies do appear, promoting them to the shared file is
a cheap change on `main` afterwards.

`reduction_notes` carries the held-currency fact, so any trajectory produced this
way is labelled rather than mistaken for a closed loop.

### D10 — Redox balance is a numerical test on a nominal trajectory

`M_nad_c` and `M_nadh_c` are touched only by GAPD and LDH_L, with mirrored
stoichiometry, and by no declared edge. So their sum is invariant, and the test is
direct: integrate alone from the registry's initial conditions, assert
`|NAD⁺ + NADH − (NAD⁺₀ + NADH₀)|` stays within tolerance at every saved point.

The tolerance is the integrator's, not a fudge: the test sets `abstol`/`reltol`
explicitly and asserts against a bound derived from them, so a tightened solver
tightens the check rather than leaving it slack. A mutation test — perturbing
GAPD's NAD⁺ stoichiometry — confirms the check can fail, so that a passing
conservation test is evidence rather than a tautology.

## Risks / Trade-offs

- **Eight Michaelis constants differ from what the published simulator runs, one
  of them by 227×.** → The deviation is registered in `reduction_notes` with both
  values, so it reaches any result through `reduction_declarations`; and the
  extract keeps `UpstreamRow`, so producing the `Quantity`-column variant later is
  a script change, not a re-derivation. The scoping note gains a
  correction-of-record row so the next reader does not rediscover it.

- **The system is stiff.** Forward `kcat` spans 59.7 to 3204 and Michaelis
  constants span 0.0014 to 59.7 mM, and `M_13dpg_c` starts at 0.0098 mM against
  `M_pi_c` at 17.8 mM. A non-stiff solver will crawl or fail. → The module's tests
  pin a stiff solver and explicit tolerances; the choice for assembled Core A′
  stays open for wave 3, which is where trajectory cost gets measured anyway.

- **A nominal standalone trajectory is not a physiological one.** With ATP, ADP
  and Pi held, glycolysis runs as an open pathway: G6P is not replenished and
  lactate accumulates with nowhere to go. → The acceptance criterion is
  non-negativity and redox conservation, not that the trajectory looks
  biological. Anything stronger needs siblings that do not exist yet.

- **Nominal enzyme concentrations are a stand-in that could quietly become
  load-bearing.** A result computed before the hook exists would be reporting
  fixed enzyme levels. → They are marked nominal rather than measured, and the
  held-currency note in `reduction_notes` makes an uncoupled trajectory
  self-identifying.

- **The vendored extract could drift from upstream.** → It is generated, not
  typed; the generator is committed; the thirteen `conc_` rows are checked against
  the registry on every load, which catches drift in the part of the extract the
  registry also holds.

- **`peer` unnamed means a missing counterpart passes.** → Accepted, and stated:
  wave 3 asserts completeness. This is the wave-0 contract's explicit position,
  not a shortcut taken here.

## Migration Plan

Additive. One `include` in `src/InferCell.jl` after the registry, one `include` in
`test/runtests.jl`. No existing sub-model, test or exported name changes, so
rollback is reverting the branch.

The include sits after `registry.jl` as a reading convention; every registry call
is inside a function body, so load order is not a constraint. The extract is data
loaded at construction, not at package load, so `using InferCell` does no file
I/O.

## Open Questions

- **Which of the fifty-two kinetic parameters is the control target.** The
  scoping note wants "one glycolytic `k_cat` with a genuinely informative prior"
  as the check that recovery works. FBA's forward constant is the loosest
  (gstd 1.7466) and LDH_L's the tightest (1.0512), so the two ends are known —
  but the choice belongs with the change that runs the recovery, and it changes
  nothing here because D5 exposes all of them.

- **Whether the assembled model wants the `Quantity`-column variant as a
  sensitivity run.** Cheap to produce once the generator exists, and it would
  measure how much the column choice actually costs. Not needed to build the
  module, and it would not change these specs.
