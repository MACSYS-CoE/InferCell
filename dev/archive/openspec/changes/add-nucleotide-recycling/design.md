## Context

See `proposal.md` — Why. Beyond the constraints the sibling wave-1 changes already
record — frozen interface, no network on compute nodes, six branches open, a
module writes only its own derivatives — three things are particular to this one:

- **This module owns the pools everything else routes through.** Eight states,
  and every currency and deferred-counter edge in Core A′ terminates on one of
  them. That makes it the module whose correctness is hardest to see from inside
  itself: alone, with nothing draining ATP, the checks it most needs to pass are
  trivially satisfied.
- **Every value it imports is ambiguous.** All five reactions and all eight
  initial conditions exist in both balanced files. There is no value here that can
  be imported without a decision.
- **A full cell cycle is 6,300 s of stiff integration**, and the module's headline
  check requires exactly that, because 144 s of it looked fine on the model that
  was wrong.

## Goals / Non-Goals

**Goals**

- Five reactions, eight states, constructible and integrable alone.
- Every import carrying a declared governing file, a recorded rejected
  alternative, and an informedness derived from the source rather than asserted.
- Conservation demonstrated against a drain, not asserted in a vacuum.
- The dead end of record reproducible on demand, so the module's reason for
  existing is a test rather than a comment.

**Non-Goals**

- The lumped charging step itself. This module is built against a drain that
  imitates it; `add-lumped-trna-charging` owns the real one.
- Executing the metabolite coupling with glycolysis. Declared, unexecuted, as in
  both siblings.
- Asserting that the Core A′ boundary is closed. This module owns the states that
  make that question meaningful and still cannot answer it alone; wave 3 can.
- Deciding whether ADK1 belongs in the inference target set. That depends on the
  charging stoichiometry, which is a different change's.

## Decisions

### D1 — The same modular rate law, with repeated stoichiometry expanded

`addReactionsToModel` (`defMetRxns.py:1474`) runs one loop over the nucleotide,
central, lipid, amino-acid and cofactor modules and builds every one of them with
`Rxns.Enzymatic(substrates, products)`. So the rate law is identical to
`add-central-glycolysis` D1, and the two ODE modules differ in their parameters,
not their form.

The wrinkle here is stoichiometric coefficients above one, which Core A′'s
glycolytic reactions do not have. The Python builds one rate-law term per unit of
stoichiometry, so `2.0 M_adp_c` becomes `Prod1` and `Prod2` both bound to ADP with
the same Michaelis constant:

```
ADK1   v ∝ kcatF·(amp/Kamp)·(atp/Katp) − kcatR·(adp/Kadp)²
       den = (1+amp/Kamp)(1+atp/Katp) + (1+adp/Kadp)² − 1
PPA    v ∝ kcatF·(ppi/Kppi) − kcatR·(pi/Kpi)²
       den = (1+ppi/Kppi) + (1+pi/Kpi)² − 1
```

That squaring is why the module has 17 Michaelis constants across 5 reactions
rather than 19: ADK1 and PPA each name one product species twice.

PPA is written reversibly in the reconstruction, and the scoping note's
`PPi -> 2 Pi` shorthand is about water not being a state, not about direction. The
reverse term is not negligible — phosphate sits at 17.8 mM against a Michaelis
constant of 0.0976, so `(pi/Kpi)² ≈ 3.3e4` — but at the registry's initial
conditions the forward term still wins by 4.3×. Carrying it reversibly is both
faithful and necessary: it is what bounds pyrophosphate. See Risks for the
steady-state arithmetic.

### D2 — Two extracts, so the ambiguity is real rather than described

| File | Logical name | Contents |
|---|---|---|
| `src/organisms/coreA/data/nucleotide_recycling.tsv` | `nucleotide_balanced` | 10 kcat, 17 K_M, 8 `conc_` rows |
| `src/organisms/coreA/data/nucleotide_recycling_central.tsv` | `central_balanced` | the central file's rows for the same 35 identifiers |

Vendoring only the governing file would make every `load_parameter` call see a
single holder, and the `governing` declarations would pass without ever being
exercised. The second extract is the point: with both loaded, every one of the 35
identifiers has two holders, so an undeclared import fails and a declared one
records what it rejected. A test asserts that removing the rival extract makes the
governing declarations stop being load-bearing — the check that the check is
working.

`central_balanced` is the same logical name `add-central-glycolysis` uses, so if
both modules are composed the provenance strings agree and a later change can
merge the two extracts without renaming anything.

The `conc_<species>` prefix on the eight initial conditions fires
`_check_registry_agreement` against `registry.jl`, which records a value for all
eight. Since the registry also records a `source_file` per species, the module
additionally asserts that its governing choice equals that field — the loader
checks the value, this module checks the file, and between them the registry and
the extract cannot disagree in either dimension.

### D3 — The nucleotide file governs all five reactions, and the data says why

The scoping note flags PGK3 and PYK3. Extending that to all five is not a
judgement call; the central file's rows for these reactions are degenerate:

| Reaction | nucleotide fwd kcat | gstd | central fwd kcat | gstd |
|---|---|---|---|---|
| PGK3 | 140.8 | 1.0513 | 319.4605 | 2.11e+63 |
| PYK3 | 672.84 | 1.0513 | 1873.9829 | 2.11e+63 |
| ADK1 | 319.156 | 1.0513 | 228.0578 | 7.01e+35 |
| GK1 | 410.227 | 1.0513 | 7.5991 | 2.11e+63 |
| PPA | 646.727 | 1.0513 | 583611.6071 | 1.69e+33 |

All 17 of these reactions' Michaelis constants read exactly 0.1 mM at gstd 10 in
the central file — the prior median at prior width, 17 times over. So the choice
is mechanical: a gstd of 1e33 is balancing with no data, and the loader's own rule
classifies anything at or above prior width as a prior default. The governing
declaration therefore agrees with what informedness already says, and a test
asserts that agreement rather than restating the decision.

Note which way the errors would have gone. Taking the central file would have run
pyrophosphatase **902× too fast** and guanylate kinase **54× too slow** — and the
first of those makes a phosphate-closure check pass more easily, which is how it
would have survived review.

### D4 — Three reverse constants are uninformed even in the governing file

Every forward catalytic constant here carries gstd 1.0513, the tightest band in
Core A′. Three reverse constants do not:

| Constant | gstd | classification |
|---|---|---|
| `kcatR_R_PYK3` | 20.9517 | prior default |
| `kcatR_R_GK1` | 31.1453 | prior default |
| `kcatR_R_PPA` | 71.539 | prior default |

These load as `:prior_default` because their geometric standard deviation is at or
above the prior width, and that is correct rather than a misclassification: a gstd
of 71 means the balancing learned nothing about the constant, whatever its point
value. They are not asserted priors — the source does supply a distribution — but
they are not evidence either, and a change that later frees one should know which
of the two it is freeing. `uninformed` is retrievable per parameter, so it does.

Priors are `LogNormal(log(mode), log(gstd))`, and every parameter is
`fixed = true` by default with a `free` keyword, as in `add-central-glycolysis` D5.

### D5 — Enzyme concentrations, and the two shared with glycolysis

| Reaction | Locus | Copies | mM |
|---|---|---|---|
| ADK1 | JCVISYN3A_0651 | 213 | 0.010555 |
| PPA | JCVISYN3A_0344 | 190 | 0.009415 |
| GK1 | JCVISYN3A_0203 | 186 | 0.009217 |
| PGK3 | JCVISYN3A_0606 | 411 | 0.020367 |
| PYK3 | JCVISYN3A_0221 | 551 | 0.027304 |

The last two are the same enzymes `PGK` and `PYK` use, at the same
concentrations, because the copy number is the gene's and both reactions read it.
Two modules each carrying a nominal for one enzyme is a divergence waiting to
happen, so the module exposes its locus-to-concentration map and a test asserts
it agrees with `add-central-glycolysis`'s for the two shared loci. Until the
wave-2 hook supplies live counts, both derive from copy number at the registry's
initial volume, so agreement is a property of using the same rule — the test
records that they do.

Sharing catalysis also means the two reactions compete for one enzyme in the real
cell, which the reduction does not model: each rate law gets the full
concentration. That is inherited from the published model, which does the same, so
it is not labelled as ours — but it is worth knowing it is there.

### D6 — Conservation is checked against a drain, not in a vacuum

Alone, this module conserves adenylate and guanylate trivially: every reaction it
carries preserves both totals, so a standalone check passes even if the module
does nothing useful. That is not evidence.

So the tests build a `ChargingDrain <: AbstractSubModel` double — the smallest
thing that reproduces the failure of record — draining ATP to AMP and PPi at 553.1
per second, the rate 3,484,518 charging events over a 6,300 s cycle implies. Three
runs, all over a full cycle:

| Configuration | Expected |
|---|---|
| all five reactions | adenylate conserved, ATP positive throughout, PPi settling |
| ADK1 removed | ATP and ADP exhausted at ~141 s (78,111 particles / 553.1 s⁻¹) |
| PPA removed | PPi rising without bound |

The middle row is the module's reason for existing, expressed as a test that fails
when the module is broken. The ~144 s the scoping note quotes uses the full
adenylate pool of 79,790 particles; the module's own arithmetic gives ~141 s
because AMP is already at 0.0832 mM. Both are the same number to the precision
that matters, and the test asserts an order of magnitude rather than a value.

Phosphate closure takes two forms, because three reactions close it exactly and
two do not. With PGK3 and PYK3 inactive, `Pi + ATP·3 + ADP·2 + AMP + GTP·3 +
GDP·2 + GMP + PPi·2` is exactly invariant. With all five active, PGK3 and PYK3
import phosphate on 13DPG and PEP, so the check subtracts what crossed the two
inbound mass edges. Relaxing the tolerance instead would have hidden a real
error; accounting for the flux keeps the check sharp.

Adenylate and guanylate are checked separately so a failure names the moiety, and
every tolerance derives from the integrator's settings, as in
`add-central-glycolysis` D10.

### D7 — Edges: declare the crossing you are the only holder of

This module owns the pools, so under the sibling changes' rule — declare an edge
where your own reactions touch a species another module also touches — it would
declare a currency edge in both directions on almost everything it owns, because
its reactions are reversible. That would be seventeen edges saying almost nothing.

The rule adopted instead: **declare mass edges for every species you do not own,
and currency edges only where you are Core A′'s sole producer or sole consumer of
a pool you do own.**

| Species | Kind | Dir | Owned | What it asserts |
|---|---|---|---|---|
| `M_13dpg_c` | mass | in | no | PGK3 draws |
| `M_pep_c` | mass | in | no | PYK3 draws |
| `M_3pg_c` | mass | out | no | PGK3 supplies |
| `M_pyr_c` | mass | out | no | PYK3 supplies |
| `M_amp_c` | currency | in | yes | the only route AMP has back |
| `M_gmp_c` | currency | in | yes | the only route GMP has back |
| `M_ppi_c` | currency | in | yes | the only route PPi has out |
| `M_pi_c` | currency | out | yes | phosphate returned from PPi |
| `M_gtp_c` | currency | out | yes | the only source of GTP |

Nine edges, and the five currency ones are exactly the statements the scoping
note's dead-end analysis needed someone to make. ATP and ADP get no edge here:
glycolysis and the stochastic blocks are the counterparts on both sides, and
declaring the crossing twice would add no information.

`inputs` holds two names, `M_13dpg_c` and `M_pep_c`. The five currency edges name
owned states, so `_check_inputs_consistency` exempts them (`e.species in own &&
continue`), and adding them to `inputs` would be wrong. Every edge leaves `peer`
unnamed, for the reason `add-central-glycolysis` D7 records.

**Cross-module consistency**, checked because the resolver fails on one species
described as both mass and currency in the same direction:

| Species | glycolysis | PTS transport | this module |
|---|---|---|---|
| `M_pi_c` | currency in | — | currency out |
| `M_pep_c` | mass out | mass in | mass in |
| `M_pyr_c` | mass in | mass out | mass out |
| `M_3pg_c`, `M_13dpg_c` | — | — | mass |
| `M_atp_c`, `M_adp_c` | currency in/out | — | — |

No species is described as both kinds, in either direction, anywhere. A test
composing all three modules' declarations asserts it, so the three branches cannot
drift apart unnoticed.

### D8 — No column choice, and that is worth recording

All 17 Michaelis constants agree between the nucleotide file's `Quantity` table
and its balanced `Parameter` table. The hand-patch layer that split
`add-central-glycolysis` — 8 of 32 disagreeing, one by 227× — has no counterpart
here.

So this module inherits every prior cleanly and labels no deviation, while its
sibling labels one. The asymmetry is the source model's, not ours, and the design
records it so that a reader comparing the two modules' `reduction_notes` does not
read the difference as an oversight.

## Risks / Trade-offs

- **A full-cycle run is the module's headline test, and 6,300 s of stiff
  integration is not free.** Forward catalytic constants span 140 to 647 while
  reverse ones reach 1.4e-3, and phosphate sits at 17.8 mM against a Michaelis
  constant of 0.0024. → A stiff solver with explicit tolerances and saved points
  at the cadence the check needs rather than every step. If the run proves too
  slow for the suite, the honest fix is to mark it and run it on Slurm, not to
  shorten the interval — shortening it is precisely the error of record.

- **Pyrophosphate's steady state depends on a reverse constant that is
  uninformed.** PPA's reverse constant carries gstd 71.5, and phosphate enters its
  rate law squared. At the drain rate, forward flux balances demand at roughly
  0.39 mM pyrophosphate — comfortably bounded, and ~4× the starting pool. → The
  test asserts that pyrophosphate settles to a bounded value rather than
  asserting the value, and the sensitivity of that plateau to the reverse constant
  is recorded as something a later change can measure rather than something this
  one pretends to know.

- **Conservation passes trivially without a consumer.** This is the failure mode
  that let the original dead end through review. → D6's drain is not optional
  scaffolding; the specs make the drained runs the acceptance criterion and the
  ADK1-removed run the proof the check can fail.

- **Two modules carry a nominal for one enzyme.** → Asserted equal by a test
  across the two changes, and flagged for wave 2, where live protein counts make
  the question moot by removing both nominals.

- **The rejected extract could rot.** Nothing reads
  `nucleotide_recycling_central.tsv` except the ambiguity machinery, so an error
  in it would surface only as a governing declaration that stops failing. → A test
  asserts every one of the 35 identifiers is held by both extracts, so a truncated
  rival file fails loudly rather than quietly disarming the check.

## Migration Plan

Additive. One `include` in `src/InferCell.jl` after the registry, one in
`test/runtests.jl`. No existing sub-model, test or exported name changes, so
rollback is reverting the branch. The two extracts are read at construction, not
at package load.

If `add-central-glycolysis` lands first, `src/organisms/coreA/data/README.md`
already exists; this change appends to it rather than replacing it.

## Open Questions

- **Whether the charging step's stoichiometry changes ADK1's importance.** The
  scoping note asks whether a `trna + 2 ATP -> trna_chg + 2 ADP + 2 Pi` lumping
  would be preferable; it closes the same moieties by fiat and puts zero traffic
  through ADK1 rather than 553 per second, which changes how much ADK1's
  equilibrium constant can move the ATP/ADP ratio the 60 s rebuild reads. It is a
  cheap comparison run and it decides whether ADK1 belongs in the inference target
  set. It changes nothing here: this module implements the published reactions
  either way, and only the drain double's stoichiometry would differ.

- **Whether the full-cycle checks belong in the default test suite.** They are the
  module's headline evidence, and they are also the slowest thing wave 1 adds.
  Deciding needs the measured wall-clock, which the implementation produces.
