# Handoff

**Session date:** 2026-09-05
**Branch:** `phase4`

## Latest: phase 4 — the 60 s rebuild (2026-09-05)

Spec §11 phase 4, delivered as PR #45 from branch `phase4`. **The coupling is
now bidirectional in execution, not only in topology.** Live ODE pools re-enter
the stochastic block as recomputed rate constants, held piecewise-constant
between refreshes as the published model holds them. Five of the seven edge
kinds now execute — mass and currency (phase 1), catalytic and deferred counter
(phase 3), rate constant (this phase) — and **two remain declare-only: volume
(phase 5) and clamped.**

**The mechanism (task 4.1): the constants live in the jump block's parameter
vector**, written by the hook, exactly as phase 3's catalytic channel writes the
ODE block's. The rejected alternative is the one the drafted transcription
design chose — mutable state on the sub-model struct — and the reason to reject
it is not style: a constant held there is invisible to `remake(prob; p = θ)`
*and* to `model_free_params`, so it could never reach a posterior, the
provenance table T1, or K6's identifiability Jacobian.

**Two protocol functions, both defaulting to empty**, the rate-constant twins of
phase 1's `contributed_states`/`contributions`:

```
rebuilt_params(m)               the module's own free parameters the rebuild fills
rate_constants(p, t, m, pools)  their values, given the pools its inbound
                                RateConstantEdges name, in mM
```

`RateConstantEdge` carries no `param_slot`, and one pool feeds many constants —
the nucleotide pools set all seventeen transcription rate constants — so the
edges name the inputs and the module names the outputs. `pools` arrives as an
argument rather than through `inputs()`, which phase 3 forbids across the
boundary.

**Where it fires.** After the debit and before the SSA step, so the pools it
reads are the ones the hook just left behind — the published order. One
correctness consequence worth knowing: phase 3 gated `reset_aggregated_jumps!`
on `isempty(d.debits)` alone, and a rate-constant write rewrites what every
propensity is computed from, so the gate is now "a debit or a rebuild happened
this handshake". A rebuild fires at the handshake whose end time is a multiple
of the edge's interval, so the declared nominal values hold for the first
interval and there is no rebuild at t = 0.

**Measured, not asserted:**

| What | Measured |
|---|---|
| task 4.2 refresh count | **10** in 600 s at a 60 s edge, **20** at a 30 s edge — counted from the recorded jump parameter vector, not from the schedule |
| task 4.2 written value | bitwise equal to the module's own law at that handshake's recorded pool |
| task 4.3 piecewise-constancy | rebuilt slot bitwise unchanged across each inter-refresh window while the pool moves at >150 of 180 handshakes |
| task 4.5 elasticity | 15 (km, pool) settings against the closed form `km/(km+pool)`, worst relative error **7e-9**; gain spans **0.0025–0.667**; the `km = 0.175`, 3.6529 mM row sits at **0.0457**, inside the 0.044–0.051 band phase 10 compares against |
| task 4.6, dynamics | drain-aligned and paired per seed: 5 s at **−0.17% ± 0.16%**, 60 s at **−0.04% ± 0.17%** — *neither resolved above noise*, both far below D10's 1% |
| task 4.6, sawtooth | reported separately: 60 s at **+0.65% ± 0.02%** against a closed form of **+0.64%**; 5 s at **+0.06% ± 0.02%** against **+0.04%** |
| task 3.8 re-run | `run_handshake!` now copies the jump parameter vector per handshake: worst extrapolation **0.019 s** per 6,300 s trajectory, up from 0.014 s, still ~500× inside K1's budget |

Full tables in `dev/scripts/rebuild_channel_result.md`, which names the job and
commit that produced it.

**Three things a reader should not over-read.**

1. **The granularity measurement was got wrong first, and how it was wrong is
   the part phase 13 needs.** The first version sampled at every handshake, and
   its headline "0.65% granularity cost" was the *sawtooth of the outstanding
   debit*: between drains a coarse configuration holds up to `drain − interval`
   seconds of unpaid cost, so its pools sit high by exactly that, and `argmax`
   landed at t = 599 — the instant furthest from the last drain. In closed form
   59 s × 40 particles/s ÷ 20,180 ÷ 18.13 mM = 0.00645 against the 0.0065
   reported: three significant figures, a deterministic bookkeeping offset
   needing no simulation. Three method points follow, all now in §12's
   amendment B and in D10 itself: sample at instants that are multiples of
   every drain compared; pair the difference per seed, because coarsening the
   drain changes how often the propensity aggregation is rebuilt and so
   consumes different randomness (three configurations at one seed give three
   different paths); and state the multiplicity behind any maximum, since a
   maximum over a noisy field is a selection statistic. With all three applied,
   **the dynamical effect is not resolved above noise at either granularity**
   on this toy, and the sawtooth agrees with its closed form to within one
   standard error.
2. **The granularity numbers are the toy's.** Its ATP pool is ~404,000 particles
   against a ~40 particle/s drain; the pools D10 is really about turn over in
   109 s (adenylate) and 30 s (guanylate), which is why D10 expects 60 s to fail
   *there* and 5 s to pass. Both rows here are below 1%; that bounds the
   mechanism and nothing else.
3. **The elasticity is the toy's law, not the model's.** What the table
   establishes is that the diagnostic recovers a *known* elasticity exactly, so
   phase 10's number will be the published rate law's and not the estimator's.

**A continuous cadence is refused, by name (task 4.4).** The spec offered
either implementing it or rejecting it; rejection is the honest option, because
the outer-loop mechanism of task 3.1 holds a constant between refreshes by
construction — which is precisely what makes task 4.3's assertion possible —
and a genuinely continuous cadence needs the `JumpProblem`-over-`ODEProblem`
shape 3.1 rejected. Running it at 1 s and calling that continuous would be a
label that overstates what runs. The error points at a shorter
piecewise-constant interval instead, which is what §10 R11's cadence comparison
actually varies. A declared-but-unexecuted continuous edge is still labelled a
deviation by `reduction_declarations`, unchanged.

**Ten refusals phases 10–12 inherit.** A rebuilt name that is not one of the
module's own free parameters; a name another module in **either** block also
declares (the `param_slot` dedup trap arriving on the rebuild side — with
seventeen genes and shared gene-expression globals this is the likely real
mistake, and across the boundary the two parameter vectors are separate so the
hook would rewrite one copy and leave the other); a pool no ODE module
integrates; rebuilt parameters with no edge; an inbound edge with nothing to
rebuild; an ODE module's outbound edge no jump module consumes; a *jump*
module's **outbound** edge, the direction convention read backwards; an `:ode`
module declaring `rebuilt_params`; two intervals on one module; and an interval
that is not a whole number of handshakes. Phase 3's catalytic `param_slot`
guard had the same one-block name scan and is widened to match.

**`drain_interval` is new on the driver**, defaults to the handshake interval,
must be a whole number of them, and is a labelled reduction when coarser. With
it came `driver_declarations`, which enumerates departures carried by the
driver's *policy* rather than by any module's declarations —
`reduction_declarations` takes models and cannot see them. It labels the coarse
drain and, as a carried-along fix, a non-fractional rounding policy, which
phase 3's own docstring already called "a labelled departure" while nothing
labelled it. `REDUCTION_CATEGORIES` gained `:coarse_drain` and
`:rounding_policy`.

**The pre-merge review changed the phase materially**, and `/check-PR`'s
report is the record: two blocking findings (the granularity number above, and
a §12 amendment that was owed and unwritten) plus eight majors. The three that
matter downstream: `clipping_census`'s fraction was per *handshake*, so a
coarse drain diluted it by exactly `steps_per_drain` — at 60 s a run in which
every debit clipped would have reported 1.7% and passed K5's 5% gate, and it is
now per drain with the pending accrual surfaced; `driver_declarations` could
not reach `reduction_report`, and the one script running a coarse drain never
called it, so the committed result file carried no label at all; and two
refusals were missing, a jump module's outbound rate-constant edge (the
direction convention read backwards, and the only remaining *quiet* mistake)
and an ODE module sharing a rebuilt parameter's name.

**Suite:** 1230 passed, 0 failed (job 16221973); 1126 before the phase.

### Next steps

1. Phase 5 (growth and volume) is the natural next one: it makes `radius_nm`
   live, and the conversion already takes it as an argument rather than closing
   over it, so that is a call-site change and not an equation change. It is the
   user's call to start it.
2. The ODE track (phases 6–9) remains independent and may still fan out.
3. When writing any module: everything in phase 3's list still holds, plus —
   a rebuilt rate constant must be one of the module's *own* free parameters and
   must carry a name no other module in *either* block uses; the rebuild law
   must be a function of the pools and of the module's *other* parameters, never
   of the slot it fills, which would compound its own previous value; and the
   consuming module declares `direction = :in`, the pool's owner `:out`.
4. **Known and deliberately unfixed:** a rebuilt slot is still an entry of
   `model_free_params`, so `build_turing_model` and the ABC path sample it and
   the hook discards the draw at the first refresh. Harmless today — no hybrid
   inference path exists — but at phase 10 that would be seventeen posterior
   dimensions coming back shaped like their priors. Phases 15 to 17 own the
   fix; the consequence is documented on `rebuilt_params`.
5. Check 7 and K5 are scored at the **published 1 s drain** (§12, amendment C).
   A coarse `drain_interval` makes each debit `steps_per_drain` times larger
   against the same pool while there are that many fewer of them, so its census
   is about our choice rather than the published model's.

## Previous: phase 3 — the 1 s handshake, on a two-module toy (2026-09-05)

Spec §11 phase 3, **the kill phase**, delivered as PR #44 from branch
`phase3`. **Neither of the phase's two stated kill conditions fired, and K5 did
not fire on the toy. The project continues.**

**What now runs.** A mixed ODE/jump composition returns a `HandshakeDriver`
instead of throwing (spec §2, G1 — the one-line refusal at
`src/orchestrator.jl:27` is gone). Phase 1 made mass and currency edges
execute; this phase adds catalytic and deferred counter, so four of the seven
edge kinds now run and three — rate constant, volume, clamped — remain
declare-only (G2). Every simulated second: protein counts fill
rate-law parameter slots through a `CatalyticEdge`, the metabolic block
integrates one second, accrued costs are debited against the pools through a
`DeferredCounterEdge` under the published clamped policy, and the stochastic
block advances one second.

**Mechanism (task 3.1): an outer loop over two stepped integrators**, mirroring
the published `hookSimulation` at `delt = 1.0` s. The two rejected candidates
are in the PR. The one worth knowing: a `JumpProblem` over an `ODEProblem` — the
SciML-native hybrid, and the obvious choice — makes propensities track the ODE
pools *continuously*, which is a different model from the published one and
leaves phase 4's task 4.3 ("propensities bitwise unchanged between refreshes")
with nothing to be true about. No new dependency: `init`, `step!` and
`reset_aggregated_jumps!` come from SciMLBase and JumpProcesses, both already
direct deps, so `Project.toml` is untouched.

**The two blocks keep separate state vectors** — ODE in mM over `Float64`, jump
in particles over `Int`. That is what makes the conversion a real boundary, and
it is the fact phases 4 to 12 most need to hold onto.

**Measured, not asserted:**

| What | Measured |
|---|---|
| conversion factor | 20,180 particles/mM at 200 nm, **derived** from Avogadro and a sphere |
| check 0, fractional carry | drift **0.0** particles at 630 and at 6,300 hooks; round trip exactly zero |
| check 0, deterministic | **252.0 → 2,520.0** — exactly linear, so the spec's rejection is now evidenced |
| check 0, stochastic | RMS **11.63 → 41.12** over 60 seeds, ratio 3.53 vs √10 = 3.16 |
| task 3.4 | every count and pool **bitwise unchanged** across all 600 handshakes |
| task 3.5 | 250-particle cost on a 100-particle pool: floors at **0.0**, deficit **150.0** exactly, repaid next hook |
| task 3.6 | 400 vs 800 copies: slot filled with exactly `n/20180.39`, flux ratio exactly **2.0**, nothing else moved |
| check 7 census | **0** deficits at nominal; **0 of 200** prior draws clip |
| task 3.8 wall-clock | **1.50–1.51e-06 s** per simulated second frozen (steady across horizons), 1.57e-06 s live at 600 s rising to 2.16e-06 at 60 s → **0.014 s per 6,300 s trajectory** worst case against a 10 s budget |

**Three things a reader should not over-read.**

1. **The wall-clock is a floor, not Core A′'s cost.** One gene and two metabolite
   states against seventeen and thirty-two. It establishes that the exchange
   carries no fatal per-handshake overhead; K1 is decided by task 13.7.
2. **The clipping census bounds the mechanism, not Core A′'s counter.** The toy
   drains ~40 particles/s from ~74,000. D13/D14 put the charged-tRNA counter at
   ~553 residues/s against a pool of order 10³, which is the case K5 is really
   about, and it is re-scored in phase 14.
3. **Nothing here executes the 60 s rebuild or the volume chain.** `RateConstantEdge`
   and `VolumeEdge` are still declare-only. Phases 4 and 5 hook into the same loop.

**A new refusal phases 6–12 inherit: a state name may not be owned in both
blocks.** The blocks hold separate vectors, so one name for two states makes
every crossing that mentions it resolve to whichever block was asked first —
phase 2's aliasing bug arriving *between* blocks instead of inside one.
`_check_state_ownership` does span both blocks — it iterates every model
regardless of formalism — but skips any name the registry does not know, and
phase 10's transcripts are non-registry by task 10.2, so that is exactly the
case it misses. This surfaced because rewriting the test that pinned the old
mixed-formalism refusal showed `TranscriptionTranslation` +
`StochasticGeneExpression` had quietly *stopped* failing. Recorded as an
annotation on task 3.2 rather than a §12 amendment: it is a new rule, not a
contradicted one. Two more rules arrived with it — a module may not reach
across the boundary through `inputs()` or `written_states()`, and may not fill
a catalytic `param_slot` whose name a second ODE module also declares.

**`param_slot` now has semantics.** It was a required field on `CatalyticEdge`
that no code in `src/` read. It is resolved at build time to the declaring
module's own free parameter, and a slot that is not one is refused by name.

**Two facts about the test environment**, both cost time if met cold. The Julia
depot had been wiped, so the first run rebuilt the whole stack (~20 min). And
because caches are keyed on CPU target (`-C native`), alternating between the
login node and compute nodes invalidates them and forces a full recompile each
way — run the suite through Slurm and stay there.

**Suite:** 1126 passed, 0 failed (job 16213224); 1006 before the phase.

### Next steps

1. Phase 4, the 60 s rebuild, is the natural next one: it hooks into
   `handshake_step!` between the debit and the SSA step, where a comment marks
   the slot. It is the user's call to start it.
2. Phase 5 (growth and volume) then makes `radius_nm` live; the conversion
   already takes it as an argument rather than closing over it, so that is a
   call-site change and not an equation change.
3. The ODE track (phases 6–9) remains independent and may still fan out.
4. When writing any module: a state name may not be owned in both blocks; debits
   to ODE pools go through a `DeferredCounterEdge`, never through an affect; and
   a catalytic edge's `param_slot` must name the declaring module's own free
   parameter.

## Orientation note added (2026-09-05)

`dev/notes/2026-09-05-reading-the-coupling-figure.md` reconciles the coupling
figure against what phases 0–3 delivered, for anyone navigating by the picture
rather than by the phase list. Three points a future session should not
rediscover: `fig1r_state_graph_reduced.pdf` is superseded by
`fig1c_state_graph_corea.pdf` (tRNA charging moved to the ODE block, §12
amendment 1); the figure's HOOK, REBUILD and GROWTH devices are phases 3, 4 and
5 while its boxes are phases 6–12, so nothing merged so far builds a box; and
the inference layer of phases 15–17 has no representation in the figure at all.
Docs only — no decisions, and the spec wins on any disagreement.

## Previous: phase 2 — several jump modules compose (2026-09-04)

Spec §11 phase 2, delivered as PR #43 from branch `phase2`. **Two or more
jump modules now compose correctly, and a jump module can read and write a
declared peer's state.** The three gene-expression modules of phases 10–12 may
assume composability; the drafted transcription design
(`dev/archive/openspec/changes/add-corea-transcription/design.md`) said that
`build_problem` "already dispatches on `formalism = :jump`; the existing
stochastic gene-expression models exercise that path", which was true for one
module and false for a composition (spec §2, G3). It is true now.

**The `reactions` contract changed shape.** `reactions(m)` returns
`Reaction(rate, affect!)` values written in the module's *own* coordinates —
`rate(u, p, t, u_inputs)` and `affect!(u, u_inputs)`, the same arguments
`dynamics` receives — and `_build_jump_problem` wraps each in a
`ConstantRateJump` over views of the composed vectors using the
`SubModelContext` it previously computed and discarded. A module cannot
address a slot outside its slice whatever position it is composed at; nothing
captures `p`, so `remake(prob; p = θ)` on the ABC path is unchanged. A module
still returning `ConstantRateJump`s is refused at build time by name. The
alternative — passing an explicit index context and having authors write
`u[ctx.s[1]]` — was rejected because it relies on the discipline whose absence
is the bug. Both shipped jump models are converted and reproduce their
pre-change seeded trajectories event for event (`test/fixtures/ssa_reference.jl`,
regenerated by `dev/scripts/make_ssa_reference.slurm`).

**Peer writes are declared with `written_states(m)`** — the jump-side twin of
phase 1's `contributed_states`, and a §12 amendment (2026-09-04). The spec said
"gated on an outbound edge", but no edge kind can name a non-registry species
and phase 10 keeps the transcripts outside the registry, so decay's write to a
transcript could carry no edge. Now: `written_states ⊆ inputs` (resolver
check); `u_inputs` is a `PeerView` whose `setindex!` throws, naming module and
state, unless the write is declared, so removing the declaration makes the
write fail at its first firing; where the written state *is* a registry species
the resolver also requires a mass or currency edge on it, `:in` where the
module consumes the pool and `:out` where it produces into it, and an outbound
edge alone satisfies the inputs contract for that species (the pool is an input
only as the write channel). The converse is not required: an inbound edge may
be a pure read. An ODE module listing `written_states` is refused by the
resolver and pointed at `contributed_states`. The three declaration checks
run in `resolve_coupling`, so a module checked alone hears them; the write
refusal itself is at run time, in the orchestrator's `PeerView`. Two facts about the
fixture: it was generated on a tree whose `src/` equals the merge base
(`e03d102`), and it is valid only while the random stream is unchanged — a
JumpProcesses or Julia RNG change invalidates it legitimately, so a test 2.6
failure after a Manifest bump means "rerun the slurm script", not "regression".

**Consequences for phases 10–12.** Transcription owns the 17 transcripts (plus
counters) and writes nothing of a peer's. Translation lists the transcripts in
`inputs()` and reads them in its propensities; it writes only its own proteins
and counters. Decay lists the transcripts in both `inputs()` and
`written_states()` and decrements them through `u_inputs`; its task 12.2 now
reads "throwing once the `written_states` declaration is removed". None of
these writes is a registry species, so none needs an edge for the write; the
deferred-counter, rate-constant and clamped edges are unaffected. Debits to
ODE pools still go through the hook phase 3 builds, never through an affect.

**Done-when, measured.** Composed `BirthOwner` + `PeerBirthDeath` against a
hand-written `BirthDeath`, λ = 12, 400 replicates each at ten relaxation
times: means 12.04 (composed) and 12.12 (single), variances 13.68 and 11.71, against tolerances of 0.69 on each mean from λ and 4.9 on the variance difference. Task 2.4's peer read: over 200 replicates the mean firing count went from 19.78 at X = 20 to 39.47 at X = 40, a ratio of 1.995 against a four-standard-error tolerance of 0.155.

**Suite:** 1006 passed, 0 failed (job 16185576, after both review passes; 1002 at job 16185348 before them); 932 before the phase.

### Next steps

1. Phase 3, the kill phase: the 1 s handshake on a two-module toy, with
   wall-clock measured. It is the user's call to start it.
2. The ODE track (phases 6–9) is independent of phases 2–5 and may fan out.
3. When phases 11 and 12 are written: read peers through `inputs()`, write only
   through `written_states()`, and keep every ODE-pool debit on a
   `DeferredCounterEdge`, never in an affect.

## Previous: phase 1 — a module contributes to a state it does not own (2026-09-03)

Spec §11 phase 1, delivered as PR #42 from branch `phase-1`. **Mass and currency
edges now execute rather than merely resolve.** Two protocol functions, both
defaulting to empty: `contributed_states(m)` names registry states of *other*
modules this one adds derivative terms to, and `contributions(u, p, t, m[,
u_inputs])` returns one signed term per name as a static vector. The
orchestrator resolves each name to the owner's global index (new `contrib_idxs`
on `SubModelContext`) and adds the terms after the per-module slices are
assembled. Nothing in the five existing models changed, and the
`TranscriptionTranslation` trajectory is byte-identical to a fixture generated
before the change (`test/fixtures/txl_reference.jl`, regenerated by
`dev/scripts/make_txl_reference.slurm`).

**One amendment, in §12.** The spec said *outbound* edges. A module that draws on
a foreign pool declares an *inbound* edge, and that consumption is a derivative
term the owner never sees unless executed, so the channel and its drift check
cover mass/currency edges in **both directions** on non-owned dynamic registry
species. Approved at planning time, before implementation.

**What the resolver now enforces** (`_check_contributions_consistency`): a
contributed species must be a registry name and not chemostatted; where a module
declares coupling, the set of its contributed states must equal the set of
non-owned dynamic registry species on which it has a mass or currency edge, in
either direction. A module with contributions but no edges is *not* skipped.
Ownership of the target is the orchestrator's check, so standalone
`resolve_coupling(module)` still succeeds and reports the species as unowned
(phase 6.5 depends on that). A `:jump` module declaring contributions is an
error in the resolver; its mass and currency edges are declared and resolved but
not held to the channel, because its writes go through its reactions, which
phase 2 builds.

**Consequences for the module phases (6–9).** Every module that produces into or
draws from a pool it does not own declares the edge *and* lists the species in
`contributed_states`, returning the signed term from `contributions`. Glycolysis
contributes to ATP/ADP/Pi/NAD(H); PTS to PEP/pyruvate and the carrier states are
its own; charging to ATP, AMP, PPi. Five `CoreAStub` fixtures in
`test/test_resolver.jl` gained `contribs = [...]` for this reason.

**State container (task 1.6).** Re-taken at 32 states; numbers in
`dev/scripts/bench_rhs_containers_result.md` (Slurm job 16184994, commit b14839a):

| RHS | bytes / call | ns / call | Tsit5 solve | Rodas5P solve |
|---|---|---|---|---|
| legacy out-of-place (no contributions; Vector→SVector conversion added) | 3152 | 2087 | 3.9 ms | 18.4 ms |
| **phase-1 static out-of-place** | **0** | **26** | 0.4 ms | 7.6 ms |
| hand-written in-place `Vector` | 0 | 11 | 0.2 ms | 4.2 ms |

The in-place row is one hand-written function with no per-module dispatch, so it
bounds what a mutating framework could reach rather than what one would cost;
reaching it would rewrite every `dynamics` in mutating form and break ForwardDiff
through `remake(p=Dual)` without a `DiffCache`. The legacy row needed an explicit
result conversion: it handed each module a plain `Vector` slice, so a
broadcasting `dynamics` returned a `Vector` and OrdinaryDiffEq refused the solve;
the real models only ever worked by returning `SA[...]` by hand.

Decision: static out-of-place stays. `_build_rhs` now holds the modules as a
`Tuple` and every index as an `SVector{_,Int}`, so slicing, `reduce(vcat, …)`
and the contribution fold are all static; the composed RHS allocates nothing on
the 32-state doubles, asserted in `test/test_contributions.jl`. Two things to
know: (1) the zero-allocation claim is the framework's, not every model's —
`TranscriptionTranslation`'s own `dynamics` branches on `gene_names` at runtime
and is not inferable, so a composition containing it still allocates. Measured
in the same benchmark: 560 → 224 bytes per call, but **549 → 866 ns per call**,
because the uninferable return type now goes through dynamic dispatch inside the
tuple machinery. Nothing in Core A′ is built from that model, and the Core A′
modules return static vectors, but lifting the gene count into a type parameter
is now a real follow-up rather than a nicety. (2) At most 31 sub-models compose: `Base.Any32` matches any tuple
of 32 or more and Base stops unrolling tuple `map` there, and the guard says so.

**Pre-existing side finding, out of scope.** Multi-gene `TranscriptionTranslation`
`dynamics` allocates `Vector{eltype(u)}`, so calling it directly with Float64 `u`
and Dual `p` throws; harmless in `infer`, which promotes `u` first.

**Suite:** 932 passed, 0 failed (job 16171912); 829 before the phase.

### Next steps

1. Phase 2 (several jump modules compose) can start; it is on the framework
   track and independent of the ODE track (phases 6–9), which may now fan out.
2. When phases 6–9 are written: the `contributions` term for a currency edge is
   the *net* rate of that species across the module's reactions; declare one
   edge per (species, direction) that actually occurs.

## Previous: coupling figure redrawn with charging in the ODE block (2026-09-03)

New figure `dev/notes/figures/corea-coupling/fig1c_state_graph_corea.{pdf,png}`,
built by `make_corea.py` from fig 1r's machinery with the tRNA charging node
moved onto the ODE grid and its edges redrawn by kind. Fig 1r is untouched and
kept as the record of the wave plan's placement; spec §2 now cites fig 1c and
§12 has a dated entry. Two things to know before rebuilding: fig 1r's
byte-identity check against the shipped fig 1 raster cannot pass on the cluster
(Graphviz 2.44 renders a larger canvas), so fig 1c reports the mismatch and
continues, and fig 1r itself does not build here for the same reason. No
`--verify` step and no ImageMagick dependency. The README records the
judgement calls, including that currency traffic is drawn by fig 1's pool-node
convention rather than as a labelled edge into Nucleotide.

## Previous: executive summary added to the spec (2026-09-03)

`spec/spec.md` now opens with a §0 executive summary, eleven bullets written
for a scientist who knows whole-cell modelling but not this project. It carries
no decisions and no status; for current status it points here. The header's
"Last amended" date is set to 2026-09-03, matching the §12 entry. Nothing else
in the spec changed. Next step is unchanged: phase 0 of §11.

## Previous: the spec is amended for an unshared-parameter boundary (2026-09-03)

Four amendments to `spec/spec.md`, all logged in its §12 with evidence. **Read
that log first**; this is the summary.

**The finding.** Core A′ has **no parameter shared between the ODE block and the
stochastic block**, and that is inherited from the published model, not
introduced by our reduction. Metabolism went through SBtab balancing and gene
expression did not; the derived parameter-coupling figure draws them as separate
panels with no edge between them; and 22 of 23 upstream ODE synthetase reactions
are commented out to avoid duplication.

**Amendment 1: the lumped tRNA charging step moved to the ODE block.** This is a
**correction, not a departure**. The scoping note always put it there — its
CME-side count is 52 reactions with no charging reaction, both tRNA species are
in its dynamic-ODE-state list, and the charged pool is in its energy-interface
table. `dev/plans/reduced-syn3a-wave-plan.md` assigned it to the CME and the spec
inherited the error. The frozen registry agrees with the note, marking both
species metabolite-scale.

It also removes 3.49 million jump events per trajectory. Charging fired 553 times
a second, 353 times every other stochastic event combined, which was a
forward-cost problem before it was an inference one.

**Two costs of that move, both recorded rather than buried.** The dominant
adenylate forward channel becomes two-stage and buffered by the tRNA pool size,
which is a quantity we assert — about 81% of ATP turnover was charging, and its
flux now contains no stochastic-block quantity. And the non-smooth boundary
*relocates* rather than disappearing: translation's charged-tRNA counter debits
~553 residues per second against a pool of order 10³ particles, so it can clip on
ordinary Poisson fluctuation. Kill criterion K5 may be more likely to fire, not
less.

**Amendments 2 to 4.** The unshared-parameter finding is recorded as new decision
D13, accepted as a first-stage assumption with the extension path written up in
the survey note. The inference architecture becomes a three-block conditional
scheme (D10 rewritten), with the path-update mechanism deliberately deferred to
phase 16 — because with charging gone the stochastic block is purely first order
and may need no particles at all. And the joint-versus-cut comparison arm is
**dropped**: `src/boundary.jl` matches parameter names to decide what to pass, so
for Core A′ it passes nothing.

**One claim retracted.** An earlier draft of D10 said conditioning on the path
makes the ODE block smooth. It does not: the clip depends on the ODE pool and
hence on the ODE parameters.

**Nine defects fixed in the same pass.** The most serious: translation credited
uncharged tRNA without ever debiting the charged pool, so the composed model had
no steady state. Also the adenylate exhaustion time cannot survive a mass-action
rate law, check 8's summation theorem omitted the charging rate constant, the
charging acceptance test was circular against its own calibration, and §2's
reaction count was wrong twice. All nine are listed in §12.

**Phases 9 to 12 renumbered** so charging precedes translation: charging 12→9,
transcription 9→10, translation 10→11, decay 11→12. Translation's rate constant
reads an ODE state charging owns, so the old order had the dependency backwards.
The ODE track is now phases 6 to 9 and is **not** a clean fan-out, since phase 9
depends on phase 8.

**Also changed:** the survey note
`dev/notes/modular-bayesian-inference-heterogeneous-modules.md` gained a section
on this boundary and what the cut and iterative protocols would need if a shared
parameter ever appeared; the scoping note gained a clarifying row; the wave plan's
superseded header records where the charging error entered.

### Next steps

1. Review the amended spec, starting at §12, then D13, D14 and the rewritten D10.
2. Then phase 1, still alone: a module contributing to a state it does not own.
   It is now a prerequisite for charging as well as for the other ODE modules.
3. Open questions worth knowing: the tRNA pool size is a single asserted number
   with three consequences and no value yet, and whether protein degradation is a
   published reaction at all is unresolved (phase 11 task 9).

---

## This session: OpenSpec out, one authoritative spec in (2026-09-03)

`spec/spec.md` is now the authoritative document for Core A′ — 1,860 lines, 18
phases, 135 tasks each carrying its own verification clause. It replaces the
OpenSpec workflow entirely. `openspec/` moved to `dev/archive/openspec/`,
read-only, so the frozen wave-0 interface requirements and the four drafted
module designs stay citable.

**Nothing executable changed.** `src/` and `test/` are untouched.

### The four drafted changes were untracked

`openspec/changes/add-central-glycolysis`, `add-pts-transport`,
`add-nucleotide-recycling` and `add-corea-transcription` were never `git add`ed.
They are now at `dev/archive/openspec/changes/`, still untracked. **`git add -A`
before committing** or they are lost outright — they are not in history and
cannot be recovered. Their derived findings are absorbed into `spec/spec.md` §4
and into the scoping note's correction table, but the designs themselves are
worth keeping.

### What the spec changed about the plan

**The wave ordering is reordered, and `dev/plans/reduced-syn3a-wave-plan.md`
carries a superseded header saying so.** Four framework gaps were verified in
source this session and none of them appears in the wave plan or in the four
drafts:

| Gap | Where |
|---|---|
| Mixed ODE/jump composition is refused in one line | `src/orchestrator.jl:27`, pinned by `test/test_stochastic_ge.jl:58` |
| No execution layer for coupling — no callback, no periodic hook, no operator splitting anywhere | six of seven edge kinds are declare-only |
| Jump composition mis-indexes state and parameters — the contexts are computed and discarded, and `reactions` gets no inputs | `_build_jump_problem` |
| Inference dispatch reads only the first module of a composition | `src/inference.jl:84` |

Two of the fixes change the sub-model protocol, and the wave plan freezes those
files against module branches — so doing them after the seven modules re-opens
seven merged PRs. They move ahead of the modules. The 1 s handshake becomes
**phase 3, the kill phase**, proven on a two-module toy with wall-clock measured
before any module is built.

One drafted design's central assumption is wrong and is corrected in the spec:
`add-corea-transcription/design.md` says the jump path is already exercised by
the existing models. True for one module, false for a composition.

### Five corrections of record added to the scoping note

Absorbed from the four drafts, in the note's own inline correction table: the
Michaelis-constant column trap (8 of 32 constants differ, one by 227×); the
cross-file trap being wider than recorded (all five recycling reactions, GK1 by
54× and PPA by **902×**); the base-to-nucleotide mapping bug (makes every
transcription rate constant **1.9× too sensitive to GTP**, on the only ODE→CME
channel); the PTS phospho-split, which the model's own data files do specify; and
transcription's pyrophosphate, which the note's phosphate accounting omits.

Plus one measurement: **channel 4's elasticity is ~0.045**. Live, bidirectional,
and weak. The spec turns this into a kill criterion on the posterior rather than
on the gain, because 0.045 costs tens to hundreds of cells, which is affordable.

### Decisions taken, so they are not re-litigated

- Scope runs through synthetic recovery, coverage and per-module calibration.
- The claim is calibrated inference across the boundary; the architecture is the
  enabling result, not the claim.
- First target set is six parameters, **and it deviates from the scoping note**:
  the polymerase turnover constant is replaced by the decay constant, because
  only the product of the former and the promoter strengths enters the rate law
  and the cap never binds, so the pair is exactly degenerate. The ptsG promoter
  is added as the headline, being the only one crossing by two channels.
- ~~The protocol comparison is joint versus explicit cut.~~ **Reversed by the
  later amendment:** the cut arm is dropped too, because `src/boundary.jl` matches
  parameter names to decide what to pass and Core A′ shares none. Both the cut
  and the iterative protocol are now §7 non-goals.

### Three findings carried as unverified

Derived this session, in no note, not checked against a run. Each is marked in
§9 with the phase that settles it:

1. ~~Whether the stochastic block's likelihood is closed-form (phase 12).~~
   **RESOLVED by the later amendment**: moving charging out makes every remaining
   reaction first order, so the kernel is closed form and factorises over genes.
2. Whether the charged-tRNA channel is weaker than transcription's (now phase 11).
   The scoping note hopes it rescues the reverse direction; a first pass suggests
   the opposite.
3. Whether protein fold change over a cycle reaches two (now phase 11). A first pass
   suggests it may fall well short.

### One thing to know about CLAUDE.md

It is **gitignored** (`.gitignore:8`, the `/*.md` root pattern). It was updated
this session to point at `spec/spec.md`, but that edit is local only and will not
reach another clone or a teammate. If the spec pointer should travel, CLAUDE.md
needs force-adding or the pattern needs narrowing.

## Next steps

1. `git add -A` (**mandatory** — see above), review `spec/spec.md`, commit and
   open the PR. Phase 0's remaining items are the finding-by-finding diff (0.2)
   and the test-suite confirmation (0.6).
2. Then phase 1, alone: make a module able to contribute to a state it does not
   own. It is the prerequisite for the three ODE modules, so nothing else starts
   until it lands.
3. Then the four-way fan-out: the ODE track of phases 6 to 9 concurrently with
   the framework track of phases 2 to 5. Phase 9 depends on phase 8.

## Open questions this session did not settle

In `spec/spec.md` §9, with the three unverified findings above. The two that most
change the work: the closed-form likelihood question, and what handshake
granularity the pools require — the GTP pool turns over in half the rebuild
interval, so 60 s is expected to fail.

---

Everything below is the handoff from the wave-0 sessions, kept for context.

**Session date:** 2026-08-28
**Branch:** `refactor/split-organism-from-framework`

## This session: split the organism from the framework (2026-08-28)

`src/corea/` held two different kinds of thing under one name: the reduced-syn3A
species registry (organism data) and the edge kinds, resolver, loader and labels
(framework). Filing the framework under an organism name would have made a second
organism reach sideways into `corea/` to compose itself, so the split ran the
other way — `src/` root was already the framework layer, so the four framework
files moved up into it and only the organism kept a directory:

```
src/edges.jl  src/resolver.jl  src/loader.jl  src/labels.jl   <- from src/corea/
src/organisms/coreA/registry.jl                               <- from src/corea/
```

Four test files renamed to match (`test_corea_edges.jl` → `test_edges.jl`, and
likewise resolver/loader/labels). A test is named after the source file it
covers, so `test_corea_registry.jl` keeps its name; `corea_test_models.jl` keeps
its name because it is a shared `CoreAStub` double with no corresponding source
file, and renaming it would be churn. (It is not itself organism-specific — it
names no registry symbol. Neither split is clean for tests: all four renamed
framework tests still use real registry species names, because the framework
validates against the registry.) All nine moves are git renames, so history
follows. 829/829 tests pass, identical to the last pre-move run (job 16003213,
whose `src/` and `test/` trees are byte-identical to the merge base) — that
count identity is the check that nothing was dropped in the rename.

**Timing was the reason to do it now.** Wave 1 fans out into seven parallel
branches that all add sub-model files here. `dev/plans/reduced-syn3a-wave-plan.md`
gained a "Where the files go" section under Wave 1 fixing the destination
(`src/organisms/coreA/<module>.jl`) so seven branches do not invent seven
locations.

**What this did NOT do:** decouple the framework from Core A′. `resolver.jl`
reaches for nine registry symbols (`is_registered`, `is_chemostatted`,
`is_dynamic`, `dynamic_species`, `species_group`, `species_index`, `held_value`,
`COREA_SPECIES`, `TRANSCRIPTION_ATOL`) and `loader.jl` four (`is_registered`,
`species_entry`, the `SpeciesEntry` layout, `TRANSCRIPTION_ATOL`). `edges.jl`
and `labels.jl` name none — so the split is two clean, two coupled. Every one of
those calls is inside a function body, resolved at call time, so the include
order in `src/InferCell.jl` is a reading convention and not a load-time
constraint; the comment there says so. Parameterising resolver and loader over a
registry object is what organism #2 costs, and it is deliberately unpaid. Out of scope for the same reason: the
OpenSpec capability id `corea-interface` and `docs/api/corea-interface.md`, both
arguably misnamed now, but not worth renaming with wave 1 about to start.

---

## Previous session: archived the change (2026-08-28, after #38 merged)

PR #38 merged as `d8ecf13`. This session ran `/opsx:archive
establish-corea-interface` — resolving the decision left open in "Next steps"
below in favour of archiving. Pre-archive checks: all 4 artifacts done, 34/34
tasks ticked, and all three delta specs verified byte-identical (requirement
bodies) to `openspec/specs/corea-interface/`, so no sync ran. The change
directory moved, as a pure git rename, to
`openspec/changes/archive/2026-08-28-establish-corea-interface/` on PR #39.

**Next:** merge PR #39, then propose the seven wave-1 changes against
`openspec/specs/corea-interface/` (see step 2 under "Next steps" below —
still current).

---

Everything below is the handoff from the PR #38 sessions, kept for context.

## What this session did

Implemented OpenSpec change `establish-corea-interface` — wave 0 of the Core A′
port, per `dev/plans/reduced-syn3a-wave-plan.md`. The change was proposed and
applied in the same session, so `openspec/changes/establish-corea-interface/`
carries the proposal, three delta specs, design and tasks alongside the code.

Wave 0 is contract only: it declares and validates a boundary, and executes
nothing. Wave 2 is where the resolved edges start driving actual coupling.

### New source

| File | What it holds |
|---|---|
| `src/organisms/coreA/registry.jl` | The 32 dynamic states and 5 chemostats, named and ordered once |
| `src/edges.jl` | The seven `CouplingEdge` kinds from fig1r's legend |
| `src/resolver.jl` | `resolve_coupling`, its checks and its reports |
| `src/loader.jl` | Provenance-carrying parameter import, cross-file ambiguity report |
| `src/labels.jl` | `reduction_declarations` — enumerating what is ours, not the model's |

### Changed source

- `src/parameters.jl` — added `ParameterSource` and a seventh `provenance` field
  on `InferParameter`, behind a six-argument outer constructor so all 54
  existing construction sites (47 in `src/models/`, 7 in `test/`) work
  untouched.
- `src/interface.jl` — added `coupling`, `module_id` and `reduction_notes`, all
  with defaults. `inputs` is unchanged.
- `src/orchestrator.jl` — `_resolve_coupling` now calls `resolve_coupling`;
  `_validate_shared_params` warns when one parameter name arrives from two files.
- `src/InferCell.jl` — includes and exports for the above.

### Two decisions worth knowing

**The boundary hazards became fields, not hard-coded behaviour.** The clipped
`max(0, ·)` ATP drain and the 60 s piecewise-constant rebuild are both
declarable, and both default to what the published model does. Departing from
the published model is now something an author has to type, and
`reduction_declarations` enumerates every departure so it can reach a report.
This is how the scoping note's two open questions are handled without answering
them prematurely.

**The loader parses SBtab-shaped TSV with Base only.** `DelimitedFiles` stopped
being a stdlib at Julia 1.9, so using it would mean a new `[deps]` entry, a
`[compat]` bound, and a network fetch. Compute nodes have no network. The format
needs no quoting logic, so a small Base parser was the lower-risk choice.

## Two pre-existing problems fixed in passing

Both were broken on `main` before this branch, and both blocked this change's
own verification steps:

1. **The test suite could not run on this cluster at all.** `Pkg.test()` resolves
   a fresh test environment needing Aqua, Aqua was not in the Julia depot, and
   compute nodes have no network — so every `sbatch test/run_tests.slurm` died in
   ~20 s with `Could not resolve host: pkg.julialang.org`. Fixed by populating
   the depot from the login node, which does have network and does have Julia at
   `/apps/modules/software/Julia/1.10.5/bin/` (the wrapper needs `EBROOTJULIA`
   set and its own `bin` on `PATH`). **If the depot is wiped, this recurs** — see
   `reference_cluster_env` in auto-memory.
2. **`mkdocs build --strict` was failing.** `docs/positioning.md` was deleted in
   `390aa22` but was still referenced from `mkdocs.yml`'s nav *and* from a link
   in `docs/index.md`. Under `strict: true` either one aborts the build, so
   `.github/workflows/docs.yml` was red. Both references removed.

## Status

- `openspec validate establish-corea-interface --strict` — passes.
- `mkdocs build --strict` — passes, with the new `docs/api/corea-interface.md`
  rendering. Verified in a scratch venv, since mkdocs is not installed on the
  cluster.
- Julia test suite — **759 passed, 0 failed, 0 errored** (job 15972827, ~3m
  wall; 701 before the review-driven hardening, 269 before this change). All
  34 tasks in `tasks.md` are ticked.

Two failures were found and fixed along the way, both worth knowing about:
`SpeciesEntry`'s validating constructor originally took eight untyped positional
arguments, which is the signature Julia auto-generates, so it overwrote the
generated one and made precompilation illegal — it is an inner constructor now.
And `test/corea_test_models.jl` extended `states` after a bare `using InferCell`,
which Julia rejects; it needs an explicit `import`. That second one silently took
all five new test files with it, so the first green-looking run had in fact never
executed any of the new tests. Worth remembering as a failure mode: a passing
count that did not go *up* is not a passing run.

## Pre-merge review (2026-08-27, later session)

`/simplify` + `/check-PR` ran against PR #38: nine parallel review agents,
every MAJOR finding adversarially verified. Verdict: merge after fixes, all of
which are applied on this branch. The substance:

- **Hardened the contract** where the review found silent-pass gaps, all
  additive validation: an inbound `MassEdge` must appear in `inputs()` (the
  drift wave-1 authors would hit — only `inputs` wires a state into dynamics);
  producer-with-no-consumer is now a reported dead end (`DeadEnd` gained
  `missing_role`, replacing `consumers` with `modules`); a `ClampedEdge` on a
  chemostat must match the registry's held value; `conc_<species>` imports are
  checked against the registry row's `initial_value` (the GTP 1.6627-vs-0.1
  trap now fails at load); a `governing` declaration binds even with a single
  holder; truncated TSV rows and misspelled `Informedness` cells error instead
  of being skipped; `CurrencyEdge.pool` is registry-checked; positional edge
  construction validates (inner constructors).
- **`DeferredCounterEdge` gained `smoothing`**, required exactly when
  `clip = :smoothed` — the edge-kinds spec's "parameter controlling the
  smoothing is exposed" clause, previously unimplemented.
- **Simplified**: vocab validation consolidated on `_check_vocab` (now in
  `parameters.jl`, shared with the registry's gstd⟺informedness invariant);
  `deviates_from_published` derives from `deviation_reason`; dead code dropped
  (`_module_name`, unused accessors, write-only `SourceTable.path`, unreachable
  `state_index::Nothing` arm); six unused exports removed; corea exports moved
  into their own files so seven wave-1 branches don't all append to one block
  in `src/InferCell.jl`.
- **Spec/docs reconciled with the code**: the edge-kinds table no longer names
  a `debited_pool` field that never existed, marks which fields default to the
  published model's behaviour, and the unknown-kind failure is honestly
  resolution-time; design.md sketches and the informedness vocabulary now match
  the source.

Wave-1 advisories the review surfaced but deliberately did not act on: lower
`ResolvedEdge`s into concrete typed callback state before any hot path iterates
them; resolve chemostat `held_value`s to plain `Float64` at build time if an
RHS ever reads them; assert registry-vs-table agreement when wave 1 vendors the
real balanced tables. (A third advisory — the cost-with-no-payer throw blocking
standalone validation — was resolved by the later `claude code-review` pass
below.)

## Automated code review fixes (2026-08-27, third session)

`claude code-review` (high effort) reviewed PR #38 and reported 10 findings
plus below-cap items; each was verified against the code, design.md and the
spec deltas before acting. Seven findings plus four cleanups were fixed
directly, two were resolved as design decisions (user-approved), one was
documentation-only, and three cleanup suggestions were rejected (one would
violate the spec's "registry positions it touches" requirement; two were
YAGNI). The substance:

- **Standalone validation now works for importing modules** (design amendment):
  the cost-with-no-payer throw in `_find_dead_ends` is demoted to a report —
  the species shows in `unowned_states` and, absent a producer, as a
  `:producer` dead end. The edge-kinds spec delta's "A cost with no paying
  state" scenario was amended to match. **This gives up a check the assembled
  model needs** — in a complete composition a cost with no payer is an error,
  and now nothing fails on it — so the spec states that a successful resolve is
  not evidence of a closed boundary, and the obligation to assert completeness
  is recorded against wave 3 in `dev/plans/reduced-syn3a-wave-plan.md`. It is
  deliberately not a spec requirement here: promising behaviour this change does
  not implement is the exact failure the earlier review round flagged as
  blocking.
- **Hybrid modules unblocked**: `_check_inputs_consistency` only holds
  registry-species inputs to the typed contract, so a module with typed
  coupling can still read legacy state (mRNA, protein) through `inputs()`.
- **Clamp payloads compared edge-vs-edge**: two clamps on one species must
  agree on `held_value` even when the registry records none
  (`_check_clamp_agreement`).
- **Chemostatted inputs get a real diagnostic**: the orchestrator now points at
  the ClampedEdge + fixed-parameter pattern instead of "not owned by any
  sub-model". Wiring held values into dynamics stays a wave-2 non-goal
  (documented in design.md).
- **Loader hardening**: registry agreement uses `atol=5e-5` (the 4-decimal
  transcription's half-unit) instead of exact `==`, so full-precision wave-1
  tables won't spuriously fail; the truncation guard covers the GeometricStd
  and Informedness columns; a non-empty unparseable or NaN gstd throws instead
  of misclassifying; the default logical file name is the extension-free
  basename, matching the registry's `central_balanced`/`nucleotide_balanced`
  names (the old `basename(path)` default could never match them, and the
  registry comment claiming the loader resolves logical names to paths was
  corrected). The `conc_<species>` identifier convention is now documented in
  the registry so wave-1 tables adopt it rather than silently killing the
  agreement check.
- **Jump path parity**: `_build_jump_problem` now runs
  `_validate_shared_params`, so `:jump` compositions get the shared-parameter
  equality check and the cross-file provenance warning.
- **Cleanups**: `RateConstantEdge`'s keyword constructor no longer silently
  discards an explicit `interval` under `cadence=:continuous` (sentinel
  default); duplicate `:asserted_prior` labels deduped via `unique_params`;
  `species_in_group` reuses `_check_vocab`; two quadratic `reduce(vcat, …)`
  flattened.

`CoreAStub` gained a `form` (formalism) field for the jump-path test. Suite:
**788 passed, 0 failed** (job 15974254; 759 before this session).

After that session, `7a7d010` ran the delta→main promotion: the three specs now
also live under `openspec/specs/corea-interface/`, byte-identical to the deltas
modulo the title header. The sync that earlier drafts of this file listed as a
pending next step is therefore **done, in-branch, before merge**.

## Fourth review round (2026-08-28)

A fresh `/check-PR` (six agents, adversarial verification) plus a second
`claude code-review` pass (8 finder angles, 13 verifier passes) ran over the
whole PR. One finding was CRITICAL and design-level; everything confirmed was
fixed on this branch. Suite: **829 passed, 0 failed** (job 16003213, 1m11s;
788 before this round). `openspec validate --strict` passes.

- **The kind-agreement collision unit was redefined — the spec amendment of
  this round.** The spec (and `_check_kind_agreement`, faithfully) treated any
  two kinds on one global `(species, direction)` pair as a disagreement. But
  the published boundary puts three kinds on ATP and on GTP at once — currency
  traffic, a deferred-counter debit, a rate-constant rebuild — against only two
  directions, so by pigeonhole the flagship composition this contract was
  frozen for could not be declared without a resolver error. No test had
  composed two kinds on one species and passed, which is how three review
  rounds missed it. The amended rule: **mass and currency are mutually
  exclusive per (species, direction)** — two descriptions of one continuous
  transport — and every other kind combination coexists as distinct mechanisms.
  Both spec copies, the resolver, the docs and the tests changed together; a
  new test declares the full ATP triple and passes. The `RateConstantEdge`
  direction convention this depends on is now pinned: the module whose rate
  constants are rebuilt declares `:in`, the pool's owner `:out`.
- **Silent-pass gaps closed** (each verifier-confirmed): the converse
  inputs-drift check now covers `CurrencyEdge`, not just `MassEdge` — a
  currency consumer omitting the `inputs()` entry used to build and silently
  never receive the coupling; a `ClampedEdge` on a state another module
  integrates now throws (a value cannot be both held and evolving); a
  chemostatted species in a coupled module's `inputs()` now fails in the
  resolver with the actual remedy (drop the input, keep the ClampedEdge)
  instead of passing the contract and then always failing `build_problem` with
  circular advice — the orchestrator diagnostic remains as the backstop for
  legacy modules; `module_id` uniqueness is checked for participating modules
  (legacy duplicates still compose); `clip = :unclamped` is now a labelled
  deviation (`deviation_reason` recognised only `:smoothed`, so an unclamped
  pool — negative counts — shipped unlabelled, falsifying the "no departure
  reaches a result unlabelled" guarantee).
- **Labelling and provenance hardening**: `unique_params` prefers the
  provenance-carrying copy over a bare one, so `:asserted_prior` labelling no
  longer depends on module composition order; `ParameterSource` validates in an
  inner constructor (positional construction could bypass the informedness
  vocabulary); `ReductionLabel.category` is a validated vocabulary
  (`REDUCTION_CATEGORIES`, gaining `:unclamped_counter`), with the category
  derived per edge kind by `deviation_category` in `edges.jl` instead of a
  bare-else `isa` chain in `labels.jl`; the producer-side chemostat exemption
  is recorded in `chemostat_exemptions` like the consumer side; clamp-vs-registry
  agreement uses the shared `TRANSCRIPTION_ATOL` (5e-5) rather than exact `!=`;
  the loader warns when `role = :initial_condition` arrives outside the
  `conc_<species>` convention (the agreement check silently forfeits otherwise);
  the cross-file provenance `@warn` fires once per parameter per session.
- **Docs/specs reconciled**: `docs/api/corea-interface.md` no longer lists the
  demoted cost-with-no-payer check as a throw (the one place the bfac1a5
  demotion missed — found independently by four review agents); proposal.md and
  tasks.md 3.5/3.8 now describe the final contract; the state-registry spec's
  ownership requirement no longer calls an unowned state "an error" while its
  own scenario reports it; registry.jl names the upstream repo, commit and
  physical filenames (`Luthey-Schulten-Lab/Minimal_Cell` @ `db048ac`); the
  central fixture's header says its GTP row is *deliberately* stale rather than
  claiming to mirror the registry.
- **Simplifications applied** (from the review's Agent E): a `caught(f)` helper
  replaces the 5-line try/catch idiom at 31 sites; `RateConstantEdge`'s
  sentinel keyword constructor collapsed to a cadence-dependent default;
  `check_gradient_safety` has one return path; the loader's truncation bound is
  hoisted out of the row loop. Two suggestions deliberately not taken:
  `_owned_states`' Dict values are now read (by the clamp-ownership check), and
  the governing-file branch merge would lose the bespoke stale-table message.

## Next steps

1. Push, let CI go green, re-run `openspec validate establish-corea-interface
   --strict`, then merge PR #38. The delta→main spec sync is already done
   (`7a7d010`); do **not** run it again. After merge, either archive the change
   (`/opsx:archive establish-corea-interface`, the wave plan's assumption) or
   leave it open for amendment while wave 1 reads the contract — decide once,
   in one place.
2. Only then propose the seven wave-1 changes. A wave-1 proposal written
   against anything but `openspec/specs/corea-interface/` will invent its own
   interface instead of consuming this one. Propose all seven before applying
   any — that is the last good chance to catch an interface assumption two
   modules disagree about.

## Open questions this change did not settle

Deliberately deferred, and none of them changes the contract:

- Which clip policy Core A′ actually runs under for inference. All three are
  declarable; whether NUTS on the ODE block needs `:smoothed` is answered by
  trying it, in wave 3.
- Whether the 60 s rebuild stays piecewise-constant. Same shape — declarable,
  defaults to published, and the deviation should be measured rather than
  assumed. A wave-2 task on `add-cme-rebuild-60s`.
- The on-disk format and location of the real balanced tables. The loader takes
  a path and a logical file name; only the parser behind it would differ. Wave 1
  vendors them.

## Scope notes for the reviewer

Three things went slightly beyond the literal task list, each flagged rather
than absorbed:

- **`module_id` is a new protocol function** not named in `tasks.md`. The specs
  require coupling errors to name both modules involved; using the type name
  alone would name two instances of one type identically. It defaults to the
  type name, so nothing changes for the five existing sub-models.
- **The two `mkdocs` nav fixes** are unrelated to this change but were blocking
  its docs verification and CI.
- **`test/corea_test_models.jl`** establishes a test-double convention the suite
  did not have — no existing test file subtypes `AbstractSubModel`, because they
  all drive the five real models. The resolver cannot be tested without doubles.
