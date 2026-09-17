# Handoff

**Session date:** 2026-09-17 (phase 10 written 2026-09-10, rebased and
reviewed 2026-09-17)
**Branch:** `phase-10-transcription`

## Latest: phase 10, transcription (2026-09-10)

**What it is.** The first Core A′ module in the stochastic block: seventeen
genes, one constitutive transcription reaction each, 22 states (17 transcripts,
5 cost counters), none of them registry species. It owns no registry state,
declares eleven edges with every peer unnamed, and composes and runs alone. It
vendors its data extract under `src/organisms/coreA/data/`, the convention
phases 6, 7 and 8 landed first — the branch was written beside them and rebased
onto them afterwards.

**The rate law reproduces the spec's numbers exactly, which is the main
evidence anything is wired right.** The seventeen constants span 1.2579e-3 to
8.2909e-3 per second against the Done-when's 1.26e-3 to 8.29e-3, and the
elasticity to all four nucleotide pools spans 0.043742 to 0.050579 against
0.044 to 0.051. Neither number was tuned; both fall out of the published
formula at the balanced concentrations once the extract is right.

**Four decisions that are not derivable from the goal.**

| Decision | Why |
|---|---|
| The transcript-decay double is an acceptance criterion, not scaffolding | Phase 12 owns decay. With no consumer the counts only accumulate and every gene's time-averaged count is an order of magnitude high, so task 10.8's calibration check would be vacuous rather than weak. Same shape as phase 8.6's charging drain and phase 9.5's translation demand. It uses the published per-gene law, derivable from the extract's own length column, so it introduces no upstream data the phase had not already vendored |
| Promoter strengths are free parameters; the rate constants are the rebuilt ones | The seventeen promoter strengths are what D9 says transcript counts identify. The rate constants are derived, and the driver overwrites them at every refresh — so `rate_constants` reads the strengths from `p` and never the slots it fills, which would compound their own previous values |
| The polymerase globals are fixed | D13: the published stochastic block is a star in which ~19 scalars set all 455 rate constants, and this module inherits its hub. A free `rnaPolKcat` is a ridge through all seventeen constants at once, which is exactly the ridge D11 rejects it for |
| The four nucleotide concentrations are held on the struct *as well as* declared | A fixed parameter does not reach the composed parameter vector, so it cannot be read from `p` inside the rate law. Declaring them keeps the provenance report honest; holding them keeps a standalone run self-consistent with the pool owner before any driver exists |

**Two things the numbers said that the spec did not.** Both are §12 amendments
dated 2026-09-10.

The first is worth carrying. **Design D9's predicted steady-state band of 0.44
to 2.87 copies was computed from a constant the published model never uses.**
`MinCell_CMEODE.py` defines `rnaDegRate = 0.00578/2` at line 295 and it appears
nowhere else — not in that file, not in `MinCell_restart.py`. The reaction the
model actually adds is `DegradationRate(rnaMetID, rnasequence)` at line 444,
which is per-gene: one catalytic constant over transcript length. At the law
that runs, predicted counts span 0.332 to 2.406 against measured 0.292 to
2.178, and **sixteen of seventeen** genes agree within a factor of two rather
than all seventeen. This is the same trap as D2's cross-file ambiguity and D3's
first-minute setup constants — a number that looks corroborated because it is
written down, and is not what executes. Worth assuming the next such number is
the same until checked.

The outlier is **FBA, `JCVISYN3A_0131`, at 3.67×**, and it is not an extract
error: 775 copies matches the scoping note's own table and 0.3469 matches
`mRNA_counts.csv` line 110. Abundant protein, rare transcript — a disagreement
between two measurements, which is what an out-of-sample check is for. The
Done-when now states §3's bound, Spearman ≥ 0.7 and a factor of two for at
least fifteen of seventeen, which §3 already said and phase 10 contradicted.

The second: **the extract needs five upstream files, not the three design D8
names.** Length and base counts come from `syn3A.gb` alone, but that record
carries no AOE protein ids at all — `grep -c AOE` returns 0 against `syn2.gb`'s
454 — so the copy number runs locus → `MMSYN1_xxxx` →
`Syn3A_annotation_compilation.xlsx` → `JCVSYN2_xxxxx` → `syn2.gb`
`/protein_id` → `proteomics.xlsx`. It also carries a `!First2` column D8 omits,
because the rate law reads the first two bases' concentrations as `C₁` and `C₂`.

**A trap this phase walked into, recorded so the next one does not.** The
`!First2` column was added *after* the extract was first committed, and the
Slurm suite reads the working tree rather than `HEAD` — so the tests passed
while the committed extract, generator and README disagreed with each other.
`git status` on the data directory alone hid it, because the generator's change
was outside that path. The four commits were rebuilt. **Check `git status` with
no path filter before trusting a green suite**, and remember that a Slurm run
proves the working tree, not the commit.

**What this leaves for later phases.** Phase 11 reads the seventeen transcripts
as peer states and needs the same loci and the same extract keying. Phase 12
writes them through `written_states`, needs the `!Length` column for its own
decay law, and supersedes this phase's double for composed runs. Phase 8's
phosphate closure must see transcription's pyrophosphate — four counters' worth
— which the scoping note's PPA argument attributed to charging alone. Phase 15's
observable design inherits the proteomics circularity, which `reduction_notes`
names against all seventeen parameters.

**Channel 4's gain, for wave 3 to size its expectations against.** The
nucleotide-dependent term is about 4.5% of the rate constant's denominator and
transcript length is the rest, so the only ODE→stochastic channel in Core A′ is
live, bidirectional and weak: 0.044 to 0.051 to all four pools together.
Correcting the base mapping roughly halves the GTP half of it — 0.0079 to
0.0117 corrected against 0.0156 to 0.0196 published, a per-gene ratio spanning
1.53 to 2.45 with a median of 1.9, which is D4's factor measured rather than
argued.

**Also worth knowing.** Task 10.7's cross-check against the metabolic modules'
copy numbers is a skipped test naming phases 6 and 7 — the second `@test_skip`
in the repo, after phase 8's in `test_corea_nucleotide_recycling.jl`. And the first full-suite run
returned one failure in `test_bursty_gene_expression.jl:50`, an unseeded
500-sample Monte Carlo test whose 0.05 threshold sits at 2.48σ and so fails
about 1.3% of runs on any branch; the rerun was clean. It is pre-existing and
untouched, and it will keep doing this until someone seeds it.

**Suite.** Job 16614893 on the rebased tree, **2513 passed, 0 failed, 2
broken** in 4m40.4s, against 2265 before the phase. The pre-rebase run — job
16373862, 1708/1708 against a 1468 base — was measured on a tree without
phases 6, 7 and 8 and does not describe what merges. Both progress figures
rebuilt, with phase 10 parsed at 8/8 and phase 8 now reading as merged.

**Next.** Phases 6, 7 and 8 — the ODE-track fan-out this branch ran beside —
have since landed (#50, #49, #52), and this branch is rebased onto them; phase 9
depends on 8; phases 11 and 12 depend on this one.

---

## Previous: phase 8, nucleotide recycling (2026-09-11)

**The third Core A′ module, and the one that closes the energy moieties.** Five
reactions — PGK3 and PYK3, which make GTP a live product of glycolysis, and
ADK1, PPA and GK1, which close a moiety — owning the eight pools every other
module routes energy through. Phases 6 (#50) and 7 (#49) landed while this was
in flight and it was rebased onto them, so it inherits the vendored-extract
layout under `src/organisms/coreA/data/`, the conservation-check conventions and
the mutation-test pattern rather than creating them. With phases 6, 7 and 8
composed, `resolve_coupling` closes every species but `M_trna_c` and
`M_trna_chg_c` — which is phase 9.

**Measured, not asserted** (`dev/scripts/full_cycle_recycling_result.md`, Slurm
job 16410838, commit f6ad54f):

| What | Measured |
|---|---|
| all five reactions, full cycle | adenylate residual **5.68e-14 mM**, guanylate **7.26e-14 mM**, against single-run integrator bounds of 4.17e-8 and 2.28e-8 mM |
| ATP | never below **3.491551 mM** (initial 3.6529) |
| pyrophosphate | settles at **0.369718 mM**, moving 1.72e-15 mM over the last ten save points |
| adenylate kinase removed | ATP crosses 1% of initial at **t = 141.0 s**, against the **141.2 s** the scoping note's constant-drain arithmetic gives for the same pool |
| pyrophosphatase removed | pyrophosphate **128-fold** to 12.7902 mM, holding **72.9%** of the 35.0920 mM phosphate budget; ATP below 1% at t = 460 s |
| phosphate closure | exact **6.39e-14 mM**; flux-corrected **3.34e-13 mM** against an uncorrected drift of 0.3092 mM |
| wall-clock | **0.005 s** per 6,300 s trajectory |

**Two amendments, both in spec §12 dated 2026-09-11.**

1. **§3's exact-conservation exception gains a second gate.** §3 asked every
   conservation residual to fall fivefold when tolerances tighten tenfold. Over
   a nine-rung ladder spanning nine decades, adenylate wanders between 128 and
   1,119 ulps of its sum with no trend while the evaluation count grows 36-fold.
   It neither falls nor accumulates. But phase 6's exception, landed the day
   before, is gated on the composed right-hand side returning `sum(nᵢ·duᵢ) ===
   0.0` **bitwise**, and this composition does not — adenylate is bitwise zero
   at only 68–80 of 200 sampled states, because the terms arrive from three
   modules and do not cancel bitwise even though `nᵀf ≡ 0` exactly. So the
   exception gets a *second* gate rather than an override: `|Σ nᵢ·duᵢ| ≤ 1 ulp
   of the conserved sum` per evaluation — measured 1.000, 0.125 and 0.115
   against a mutated 1.7e16 — then the ≥6-decade ladder, bounded by each rung's
   own `tol_C` rather than a flat ulp count. **Nothing inherits it
   retroactively**: phases 6 and 7 already assert the first gate, and phase 14
   keeps the particle restatement phase 6 assigned it.

   *A first draft of this amendment was wrong and is worth remembering.* It
   claimed the fivefold rule binds only where the trajectory is restarted, and
   exempted all three moieties by argument. That reinstated a carve-out phase 7
   had already withdrawn on rebase, and the conclusion did not follow from the
   premise: if the solver preserves a linear invariant exactly, it does so after
   a restart too, so the assembled residual is set by integer rounding — half a
   particle is ~2.5e-5 mM, eight orders above the round-off floor — and does not
   scale with tolerance either. It would have handed task 14.4 an unpassable
   rule.
2. **Removing the pyrophosphatase strands the phosphate moiety; it does not make
   pyrophosphate diverge, and in a closed model it cannot.** The scoping note's
   173 mM is one pyrophosphate per charging event for a whole cycle — open-pool
   arithmetic assuming a phosphate supply Core A′ does not have, since check 4b
   says its phosphate is closed.

**Eleven edges, not the archived design's nine, and the reason is the rate law.**
The published PGK3 and PYK3 laws are reversible and read `M_3pg_c` and
`M_pyr_c` twice over — once in the reverse numerator, once in the denominator's
product bracket. An ODE module reaches a foreign state only through `inputs()`,
and `src/resolver.jl:224` holds every registry input to an *inbound* edge, with
the `written_states` escape reserved for jump modules. So a declaration carrying
those two outbound only cannot evaluate its own rate law. Each product now
carries a mass edge in both directions — phase 1's rule of one edge per
(species, direction) that actually occurs, and a true statement about a
reversible reaction — while `contributed_states` still names each once and
returns the net signed rate. Annotated in place on task 8.5, **not** a §12
amendment: the spec's own rule fixes no direction.

**A trap worth carrying: the framework omits fixed parameters from the parameter
vector entirely.** `_build_p0` and `SubModelContext.param_idxs` both come from
`model_free_params`, so a `fixed = true` constant is not in `p` and cannot be
read from it, and freeing one shifts every index after it. The archived design's
"fixed by default with a `free` keyword" is therefore not implementable as a
plain positional read. This module holds every kinetic constant on the struct
and carries a `pidx` mapping each to its slot in its *own* free list, zero where
it is held; `free = [...]` flips both together. **Phase 6 hit the same wall and
solved it the same way** — `central_glycolysis.jl` carries the identical
`pidx`-into-its-own-free-list pattern — so this is a shared convention rather
than a finding ahead of that phase.

**The five enzyme concentrations are the module's only asserted priors**, and
calling them anything else would be a friendlier label than the source supports:
the proteomics table states a copy number and no width, so the value is the
published model's and the prior is ours. §4 D7 counts thirteen asserted priors
across Core A′ — eleven transport constants plus `k_chg` and the tRNA pool —
and enzyme concentrations are not among them. **Every metabolic module that sets
a concentration from a copy number adds more, so spec task 13.6's count is owed
a reconciliation at assembly.**

**The double did not hold its pools, and that changed what the run showed.**
`dynamics` returned a zero derivative for the four glycolytic species, which is
not the same as holding them: `NucleotideRecycling` names all four in
`contributed_states`, and `_accumulate` folds a contribution into the *owner's*
`du`. So the pools moved at exactly ±v_PGK3 and ±v_PYK3 — 13DPG below 1% of
initial at **t = 0.5 s**, PEP at **t = 18.5 s** — and the GTP branch died of
substrate starvation inside the first save interval of the run whose stated
purpose is to show that PGK3 and PYK3 make GTP a live product of glycolysis.
The result file's explanation of its own inbound-flux number was wrong in the
same way: it said the 0.0505 mM crossing "is bounded because guanylate is",
while GDP ended at 79% of its initial value.

Each pool now relaxes to its registry setpoint at `k_gly` — the four lumped
reactions phase 6 supplies: GAPD makes 13DPG, PGM and ENO carry 3PG to PEP, LDH
drains pyruvate. None is in the tracked phosphate sum, so nothing enters the
closed moiety. **What bounds the branch is now guanylate, which is the honest
limit**: nothing in Core A′ consumes GTP until phase 9's translation, so PGK3
and PYK3 run until GDP is spent and the crossing is GDP₀ + GMP₀ = 0.3098 mM.
The suite asserts the pools stay within 1% of setpoint, which nothing did
before.

**What the doubles do, and one thing a first version got wrong.**
`HeldGlycolytic` owns the four glycolytic species the GTP branch reads and
rephosphorylates ADP as `ADP + Pi -> ATP`. The phosphorylation is not
decoration: nothing in this module produces ATP — the adenylate kinase returns
AMP *at the cost of an ATP* — so without a source ATP falls to zero in about
130 s in every configuration, and the kinase-removed crossing would then be
dominated by ATP draining into ADP rather than by adenylate stranding as AMP,
which is the mechanism the note's 144 s describes. An earlier version took the
phosphate from the 13DPG pool, as the published PGK does; but 13DPG does not
come from glycolysis here, so 345 mM of phosphate entered a closed moiety,
free phosphate climbed without bound, and pyrophosphate rode up to 51 mM instead
of settling at 0.371 mM. Glycolysis from G3P is `G3P + Pi + ADP -> 3PG + ATP`
once GAPD and PGK are composed, so the stand-in takes it from the free pool.

**One skipped test, naming what it waits on**: the four-module
mass-versus-currency agreement needs phase 9, which has not started. The other
skip is gone — one enzyme at one concentration across modules waited on phase 6,
phase 6 landed first, and the spec's fan-out rule makes the assertion this pull
request's. It is written, and it was worth writing: it failed on the first
attempt, because this module divided copy numbers by a transcribed 20180 where
`central_glycolysis.jl` divides by `corea_particles_per_mM()` = 20,180.39, so
`JCVISYN3A_0606` ran at 0.020367 mM in one module and 0.020366 in the other.
Both now derive the factor.

**Suite:** 2265 passed, 0 failed, 1 broken (the one remaining skip), 2266 total
(job 16410902, HEAD `68ce588`, clean tree); **2023 before the phase**, which is
phase 7's recorded count on the rebased tree. Earlier phase-8 figures of
1684/1686 against a 1468 base were measured on the pre-rebase tree — job
16374108's own header says `HEAD: 57b4e5e (tree dirty)` — and that commit is not
an ancestor of this branch. The two full-cycle testsets add tens of seconds,
almost all of it Rodas5P specialising on two new problem types — the horizon
itself is free at 3 ms a trajectory, which is why the suite asserts the
done-when at the horizon the done-when names rather than a shortened one.

**GitHub Actions is blocked on billing, not on this branch.** Every workflow
job since ~09:10 fails in three seconds with zero steps run and the annotation
"The job was not started because recent account payments have failed or your
spending limit needs to be increased." CI passed normally earlier the same day.
Until that is settled the Slurm suite is the only test signal, so quote the job
number rather than a green tick.

**A cluster fact that cost half an hour.** Precompilation caches are keyed on
CPU target and every Slurm job lands on a different node, so a cold node costs
15 to 25 minutes before a single line runs. Two jobs submitted together
precompile concurrently and both stall. Submit one at a time.

### Next steps

1. Phase 8 is **not merged**. PR #52 has completed `/check-PR`; its eight
   blocking findings are corrected. The GitHub-hosted checks did not start
   because the organization exhausted its Actions allowance, while Slurm job
   16410902 passed 2,265 tests with no failures and one intentional phase-9
   skip. Merge awaits the recorded CI waiver and the user's explicit go-ahead.
2. Phase 10 is the remaining independent module phase. Phase 9 depends on this
   one.
3. **Owed to phase 9:** the charging demand is recorded as 553.1 residues/s,
   derived from 3,484,518 residues over 6,300 s, as the *demand* the real module
   must meet — never a value to calibrate `k_chg` against and then re-check,
   which would verify arithmetic (§4 D14). A test asserts the production module
   contains no such constant, so the anti-circularity guard is mechanical rather
   than a note. Phase 8's drain double is superseded for composed runs once
   phase 9 lands and kept only for standalone ones.
4. When writing any Core A′ module: everything in phases 3, 4 and 5's lists
   still holds, plus — a fixed parameter is not in the parameter vector, so hold
   it on the struct; declare an edge per (species, direction) that actually
   occurs, which for a reversible reaction on a foreign pool is both; and a
   standalone linear invariant is preserved exactly, so its evidence is the
   mutation test rather than a tolerance sweep.

## Previous: phase 7, phosphotransferase transport and lactate export (2026-09-10)

**The second Core A′ module, and the first with a stochastic-block neighbour.**
Phase 6 (central glycolysis, #50) landed first and this branch was rebased onto
it. `PtsTransport` owns the eight carrier phospho-states and external lactate,
and it is what makes carbon balance posable at all — glycolysis makes the
lactate, this exports it. It owns the eight carrier phospho-states and external lactate, and it is
what makes carbon balance posable at all.

**Four ways the archived design had gone stale, all of them protocol rather
than chemistry.** `dev/archive/openspec/changes/add-pts-transport/` was written
against the wave-0 interface, and four things moved under it:

| The design said | What the protocol now does |
|---|---|
| `M_g6p_c` and `M_pyr_c` production is "declared and left unexecuted" | Phase 1 made contributions execute, and the resolver holds `contributed_states` to the mass edges in both directions, so all four foreign pools carry a signed term |
| Five edges: pep in, lac_c in, pyr out, g6p out, the clamp | **Ten.** Both reversible steps *read* the pool their forward term fills, and `inputs` — the only channel wiring a foreign state into `dynamics` — is held to the inbound edges. So pyruvate and glucose-6-phosphate each carry an edge in both directions |
| `membrane_protein_states` as a plain function on `PtsTransport` | The protocol function phase 5 promoted |
| (silent — it did not exist) | `extracellular_states`, which task 7.4 needs so growth does not dilute lactate that has already left the cell |

Only the second is a judgement call, and it is forced: one outbound edge alone
would declare the production and leave the rate law unable to see the pool it
draws back from.

**One thing about `fixed` that the next module phase will hit.** Every parameter
is `fixed = true` by default, so nothing is sampled from an invented prior
without an explicit act — that is design D7 and D14, and it is right. But `fixed` in this
codebase means *absent from the composed parameter vector entirely*:
`model_free_params` filters it out, so `p` never carries it and a rate law
cannot read it from there. So the fourteen scalars are carried on the struct as
well, built from the same `load_parameter` call so they cannot drift, and the
free ones are spliced back into their local slots at every call. Phase 6's extract is about 65 rows and phase 8's about 35, so both have
several dozen; they need the same thing.

`:r_cell_nm` is always free, whatever `free` is passed, because a `param_slot`
must be a free parameter to reach the vector at all. That is the trap 5b.10
documents: inference samples it and the handshake overwrites it.

**Reading `proteomics.xlsx` without `openpyxl`.** The four carrier copy numbers
— 353, 314, 290, 831 — exist upstream only in that workbook, behind the
published loader's `MMSYN1 → JCVISYN2 → AOE` chain through a *second* workbook.
No environment on this machine has `openpyxl`, and compute nodes have no
network. `zipfile` plus `ElementTree` reads it in about twenty lines: one sheet,
header on row 2, column **V** is "Absolute abundance (copy number)", which is
what the published loaders reach positionally as `iloc[:, 21]`. The four rows
are keyed by AOE id (`AOE93321.1`, `AOE93571.1`, `AOE93322.1`, `AOE93596.1`), so
the annotation workbook is not needed, and the script asserts each row's NCBI
description as well as its rounded count — a shifted sheet then fails loudly
rather than importing another protein's abundance. **Phase 10 should reuse it:
the promoter proxy is the same table's copy number over 180.**

The cross-check is `setICs_two.py:286-288`, a disabled block stating the same
four totals as concentrations. Disabled code, so not the citation — but an
independent plain-text statement of a quantity otherwise derived from a
spreadsheet, and the script asserts its own arithmetic against it to four
decimals.

**A conversion-factor divergence worth knowing about.** The archived design's
seven-decimal table used a rounded **20180** particles per mM; the code derives
**20180.3873819** from Avogadro. Four of the eight phospho values therefore
differ from D3 in the seventh decimal. The derived factor is the right one,
because the initial conditions are read back with it — building them with the
rounded one would put the carrier totals at 353.007 rather than 353, and "sums
to its published copy number exactly" would have to become a tolerance.

**The phase amended §3.** The tolerance principle asked every conservation
residual to fall at least fivefold when tolerances tighten tenfold. Measured, the
four carrier sums do not: 3.47e-17, 2.78e-17, 8.67e-18 and 1.32e-16 at
`(1e-10, 1e-8)`, eight orders under their bounds, and unchanged at
`(1e-11, 1e-9)` — the ptsI residual *rose*, to 9.02e-17. A carrier sum is a
linear invariant whose two derivative terms are exact IEEE negations, so the
right-hand side conserves it identically and there is no truncation error to
shrink. §3 had already noticed the same of the redox pair without drawing the
consequence.

**Phase 6 got there first, with a better rule.** It hit the same wall on the
redox pair, measured a ladder over eight decades, and landed an *exception* to
the tolerance principle rather than my two-class split — gated on a mechanical
criterion a phase must assert rather than argue: the composed right-hand side
must return the moiety's weighted derivative sum as literally `0.0`, bitwise.
It also fixed §6 F3, which my draft had left demanding a positive slope it had
just declared impossible.

So on rebase my §3 amendment was dropped and phase 7's check 5 now asserts
phase 6's criterion: `du[unphos] + du[phos] === 0.0` at every saved state, then
a six-rung ladder over ten decades with every rung within 100 ulps and the
largest within 100× of the smallest. **Phase 7 is the independent second
instance** — different module, different chemistry, same floor — which is worth
more as corroboration than a rival amendment would have been.

**The mutation test remains the real falsifier**, and that is the part that
matters

**The full-cycle question, answered for one case.** §9 asks whether full-cycle
checks belong in the default suite and says to decide on measured wall-clock.
Phase 7's 6,300 s lactate-export **integration costs 4.56 ms**; the **19.67 s**
around it is compiling the stiff-solver path, which the suite pays once for its
first stiff solve whatever the horizon (job 16364852). It therefore runs by
default. My first draft quoted ~19.5 s as the cost *of the full cycle* — right
magnitude, wrong meaning, and sourced from a failing run of a test version that
no longer existed. Phases 8 and 14 should expect the same shape: horizon is
nearly free on a small ODE module, and it is the handshake driver, not the
length of the integration, that will make the assembled-model checks expensive.
That settles it for a nine-state standalone ODE module and for nothing else:
checks 2, 4 and 4b run on the assembled model across a handshake driver, which
is the expensive case and is still open.

Posing that check at all needed a lactate source, because steady cytosolic
lactate is glycolysis's and glycolysis is phase 6 — so the `HeldMetabolites`
double supplies it at the published two-per-glucose rate. D6's figures come
back: 1.46836 mM cytosolic, 0.006877 mM external after a cycle, **0.468%**
against D6's stated 0.47%, and a ratio of 1e4 missing the criterion as its 4.7%
row says.

**Two failures on the way, both real and both instructive.** The first draft
asserted `steady == production/0.075`, which is wrong by exactly the external
pool — export is driven by the *difference*, so the 0.47% excess is the lactate
that has left. The second asserted the corrected form to 1e-6 and failed at
1e-5, because the system is *quasi*-steady rather than steady: external lactate
is still rising at `production/R`, cytosolic tracks it, and the export flux
falls short of production by exactly one part in `R`. The test now carries the
drift term, which is a better check than either draft.

**What phase 11 and phase 13 inherit.** The four carrier sums are conserved only
as long as nothing else makes PTS protein. Once phase 11's translation feeds the
carriers the invariant must be restated as "conserved up to what translation
adds", by subtraction rather than by widening the bound. The check as written is
scoped to this sub-model's own dynamics and says so.

**Suite.** `sbatch test/run_tests.slurm` job **16374055**, **2023/2023** in
4m58.2s on the rebased tree. Full-cycle timing is job **16364852**. Before the
rebase, on the 1468-test base, phase 7 alone was 1654/1654 (job 16364952); the
post-rebase total is against phase 6's tree and the check-5 rewrite changed
phase 7's own count, so the two deltas are not comparable and no combined one
is claimed here.

**Next.** Phases 6, 8 and 10 are the remaining fan-out, all independent of this
one. Phase 9 waits on phase 8. Three done-when clauses across the fan-out are
still not independent — task 8.4 needs phase 6, task 8.5 composes all four ODE
modules, task 10.7 checks copy numbers against the metabolic modules' — and
belong to whichever PR lands last, or to phase 13.

---

## Previous: phase 6, central glycolysis (2026-09-10)

**The first phase that writes a box on the state graph.** Every merged phase
before this one is framework — that is D0 on purpose — so `src/organisms/coreA/`
held nothing but `registry.jl`. It now holds a sub-model: ten reactions from
glucose-6-phosphate through lactate, owning the eleven glycolytic intermediates
and the redox pair, with all 65 kinetic constants and initial concentrations
imported through the loader from a vendored extract.

**1765 tests pass, up from 1468** (job 16364934, 3m32s).

**The one thing that made the spec move.** Task 6.7's redox check was written as
§3's tolerance principle prescribes — solve at `(1e-10, 1e-8)`, solve again a
decade tighter, require the residual to fall at least fivefold — and it failed on
the first Slurm run: 1.42e-14 mM against 7.99e-15 mM, a ratio of 1.78. A ladder
over eight decades then showed the residual **flat**:

| abstol / reltol | drift (mM) | ulps of the conserved sum |
|---|---|---|
| 1e-4 / 1e-2 | 1.33e-15 | 3 |
| 1e-5 / 1e-3 | 8.88e-16 | 2 |
| 1e-6 / 1e-4 | 4.00e-15 | 9 |
| 1e-8 / 1e-6 | 1.33e-15 | 3 |
| **1e-10 / 1e-8** | 1.42e-14 | 32 |
| 1e-12 / 1e-10 | 2.66e-14 | 60 |

(Six of the nine rungs measured; the suite re-runs four of them.)

`eps` at NAD⁺ + NADH = 2.2097 mM is 4.4e-16. GAPD and LDH_L are the only
reactions that touch the pair, both live in this module, and their stoichiometry
mirrors — so the composed right-hand side returns the two derivatives as
bit-for-bit negatives and the sum is conserved whatever step the solver takes.
There is no integrator error in it to shrink.

§3 check 3 already said the invariant is "exactly invariant there" while the
done-when asked for a residual that shrinks. **Both could not hold, and the
measurement says check 3 was right.**

**Read the fence, not just the exception.** Pre-merge review found the first
draft of the amendment sound in its physics and loose in its wording, and the
fixed version is the one to carry forward. The criterion is **mechanical and
must be asserted**: a check may use the exception only where a test asserts
`sum(nᵢ · duᵢ) === 0.0` bitwise on the composed right-hand side. *Crossing a
module boundary was the wrong discriminator* — it is a property of the
composition, not the moiety, and phase 8 standalone owns every adenylate
species, so tasks 7.7, 8.7 and 9.6 could each have claimed the exception on the
very configuration they are run in. The assertion carries numbers (100 ulps per
rung, 100× across the ladder), and it is scoped to **one continuous solve**,
because growth dilution rewrites every concentration at each of ~6,300
handshakes and task 14.4 would otherwise inherit a bound that breaks at
composition scope.

The residual is also not literally flat: it rises mildly with step count (3 → 32
→ 60 ulps over the last four decades), so the claim is "bounded at the floor and
not falling", never "constant". The mutation test is what keeps this from being
a weakened check: a GAPD stoichiometry of 2 rather than 1 on NAD⁺ moves the
residual by **14.2 orders of magnitude**, to 2.21 mM, the whole pool, and does
not shrink with tolerance either.

**Four decisions that are not derivable from the goal.**

| Decision | Why |
|---|---|
| The include sits after `interface.jl`, not after the registry | `<: AbstractSubModel` is resolved when the struct is *defined*, not when a method is called. The archived design's "after the registry" was written before there was a struct to place |
| The ten protein counts are declared in `inputs` and the double owns them | Task 6.4 says to declare them, but `build_problem` now refuses an input no module owns — the archived design was written when resolution was lenient. Declaring them and having the local double own them keeps the spec's text and costs one extra double. They are a forward wire for phase 11; `dynamics` reads the nominal struct values today, and `protein_sources` is a keyword so translation can re-point them under whatever names it gives its proteins |
| Enzyme concentration divides by `corea_particles_per_mM()`, not a literal 20,180 | The derived factor is 20,180.39, so five of the ten differ from the archived table's six-decimal values in the sixth decimal. A transcribed constant would also disagree with the conversion the handshake itself performs, which is the thing these concentrations have to agree with once translation supersedes them |
| The energy double is local to the test file | Three sibling branches want the same two-line double and six branches editing `corea_test_models.jl` is the conflict the fan-out exists to avoid. Promoting it once 7, 8 and 10 have landed is cheap |
| The element type reaches the rate law as a type parameter | Computing `T = promote_type(...)` as a local inside the body and using it there defeats specialisation: the composed right-hand side allocated **55,520 bytes per call** while inferring fine and returning identical values, so only an allocation probe catches it. Worth knowing for phases 7–9, which write the same shape of rate law |

**One correction of record against the archived reference.** Its task 1.1 names
`conc_M_atp_c = 3.6529` among four spot values while also fixing the
concentration count at thirteen. ATP is not one of the thirteen owned states, so
the two cannot both hold; the count is what spec §11 task 6.1 carries, and
`conc_M_g6p_c = 3.7076` is asserted in its place.
`src/organisms/coreA/data/README.md` records it, and whoever vendors the
adenylate species — phase 8, which needs the central file's rival values anyway
— picks `M_atp_c` up there.

**Four things the vendored extract falsified, all now annotated in the spec.**
§4 D1's "with its geometric standard deviation" names a column that does not
exist — the balanced table offers only `UnconstrainedGeometricStd`, so every
prior pairs a balanced median with an unconstrained width. That is free for all
13 concentrations and all 32 Michaelis constants, whose `Mode` equals the
unconstrained geometric mean exactly, and **not free for ten of the twenty
catalytic constants**, worst `R_TPI` reverse at mode 4 against 65,341.7 carrying
a 1.05 width. §4 D3's declaration understated itself: the column choice
**reverses the sign of the pathway's entry reaction** at the registry's initial
state — R_PGI at −0.3325 mM/s against +4.5963 on what the simulator runs — and
slows R_FBA by 253×. §4 D11's two controls are misdescribed (`kcat_ENO` is not
the tightest prior; `kcat_FBA` is the loosest *forward* one only), though the
target set survives. And task 6.4's "copies over 20,180 to six decimals" cannot
be met by a derived conversion factor.

**One framework gap, loud rather than fixed.** A freed initial condition is
**sampled and then ignored**: `_collect_ic_values` builds `u0` from each
parameter's stored value, and no rate law reads a `conc_` slot. The likelihood
is exactly flat in it, so the posterior equals the prior — which reads as "the
data does not constrain the initial pool" rather than as a channel that was
never wired. Phase 6 refuses to free one; phases 15–17 own the wiring, alongside
the driver-written-slot problem 5b recorded, which is the same defect in a
different channel.

**Two numbers §9's open questions now have.** Whether full-cycle checks belong in
the default suite: phase 6's six full-cycle solves cost about 1m40–2m10 of a
3m32 suite (baseline 2m01s at 1468 tests, job 16341043), which is affordable for
a 13-state ODE block; §9 is annotated with that and left open for the assembled
hybrid, which task 13.7 decides. And
`uninformed_params` on this module returns **four**, not the two prior-default
concentrations: `kcatR_R_PFK` and `kcatR_R_PYK` carry geometric standard
deviations of 33.03 and 11.34, at or above the prior width. The reverse
constants are the loose ones, which matters for D11's target set.

**What phase 6 does not do, and who owns it.** The standalone trajectory is not
physiological and is not asserted to be: G6P is never replenished, so the module
drains its own carbon into lactate and the acceptance criteria are non-negativity
and redox conservation alone. Carbon balance spans transport and export (phase
7), adenylate and phosphate closure span recycling (phase 8), and all three are
phase 14's. Task 8.4's shared-enzyme assertion needs this module and belongs to
whichever of the fan-out lands last.

**Next.** Phases 7, 8 and 10 are running concurrently on their own worktrees.
Phase 9 waits on 8.

---

## Previous: phase 5b, the inbound volume channel (2026-09-09)

**Why it exists, which is the part worth carrying.** Phases 6, 7, 8 and 10 are
mutually independent and were about to be built as four parallel PRs. R15 makes
framework code off-limits on a module branch — and phase 7 breaks that rule by
construction. Task 7.2 needs the cell radius to reach the lactate exporter's
`3P/r`, and phase 5 refuses an inbound `VolumeEdge` by name because the type
carries no parameter slot. The phase 7 author would meet a live `ArgumentError`
in a file they may not edit.

R15's escalation path is written for a need that *surfaces* during a phase. This
one was diagnosed by phase 5 and already in the spec; what was missing was a
place in the ordering. So it is now **phase 5b**, between 5 and 6, and phase 7
is a pure module phase. Recorded as a §12 amendment dated 2026-09-09.

**What the channel is.** The catalytic channel with the source removed. An
inbound `VolumeEdge` carries a `param_slot` naming a position in the declaring
module's own rate law, exactly as `CatalyticEdge` does, and the driver fills it
at step 0 of every handshake. The resolved record is `GeometryExchange`, which
is `CatalyticExchange` minus `count_idx` — and that absence *is* the channel:
the value comes from a driver field, not from any block's state.

**Four decisions that are not derivable from the goal.**

| Decision | Why |
|---|---|
| `quantity` is unit-bearing — `:radius_nm`, `:volume_litres`, `:area_nm2` | Nothing else in the codebase records a unit; neither `InferParameter` nor `SpeciesEntry` has the field. The symbol is the only place a nanometre is told from a centimetre, and a `μm`-for-`nm` slip would be absorbed into the permeability it multiplies rather than caught |
| `species` names the **declarer's own state** whose rate law reads the geometry | The geometry is a sum over every flagged state, so it belongs to no module — but `src/resolver.jl` requires every edge to name a registry species, so a pseudo-species is not available. Naming the membrane protein instead would silently become a different claim the moment task 7.6 flags a second protein |
| Written **outside** `_update_volume!` | It returns early on `isempty(d.growth)`, so a write inside it would skip exactly the fixed-cell composition phase 7 needs standalone |
| Rate laws get the geometry of the **capped** volume | `_update_volume!` caps volume and leaves radius uncapped, as `in_out.py:107` does. Harmless while the radius is only reported; wrong once it drives `3P/r`, whose whole content is that `3/r` is a sphere's surface-to-volume ratio. Labelled `:capped_rate_law_geometry` |

**One refinement the first Slurm run forced.** Deriving the rate-law geometry
from the capped volume *unconditionally* round-trips radius → volume → radius
and loses an ulp, so a fixed cell at exactly 200 nm handed its rate law
200.00000000000003. Below the cap the reported geometry is already one sphere,
so it is now passed through untouched and only the capped branch reconstructs.
The below-cap assertion is `==`, not `≈`, which is what says so.

**What 5b deliberately does not fix, and who owns it.** A `param_slot` must be a
*free* parameter, and free non-observation parameters are exactly what inference
samples — so a driver-written slot is sampled and then overwritten. Pre-existing:
`rebuilt_params` documents it at `src/interface.jl`, the catalytic channel has it
undocumented. It is sharper here, because a geometry slot is not an unknown at
all and in `3P/r` its Jacobian column is exactly proportional to the
permeability's, so a rank deficiency of one is an artefact of the declaration
rather than a finding. 5b adds the warning to `VolumeEdge` and
`driver_written_params(driver)` to enumerate the slots.

**The fix stays where `rebuilt_params` already put it: the inference phases, 15
to 17.** An earlier draft of this phase assigned it to task 13.3 — that was
wrong twice over. 13.3 is hybrid *dispatch*, its text says nothing about the
sampled set, and its test (a mixed composition dispatching or throwing) cannot
catch this, because the composition that silently samples a geometry slot is the
*homogeneous* ODE block. Core A′'s four ODE modules compose without any jump
module, so `infer` on that block runs today and would sample the radius.

**A refusal considered and rejected.** Refusing an inbound edge in a homogeneous
`build_problem`. Rejected on two grounds: all four cross-block channels are
accepted by a homogeneous build today, deliberately and documented; and task
phase 7's own done-when wants the export rate checked "at the registry's
radius", which is a homogeneous-path assertion. Instead the module declares its
geometry slot as the registry radius, and `src/interface.jl` records that the
driver writes 200.03505 nm over it — the two agree to the four significant
figures the source states the radius in and differ beyond them, deliberately.

**Also closed, as a side effect.** One writer per ODE parameter slot is now
enforced across the catalytic and geometry channels both (`_claim_slot!`). Two
catalytic edges into one slot is a hole the catalytic channel had today: both
lower, both write, the later wins, and the loser is declared, lowered, reported
and has no effect.

**Suite.** See the run recorded in spec task 5b.11. Both progress figures were
rebuilt and all five of their build-time guards re-checked in a scratch copy,
plus a sixth this phase adds: a `### Phase` heading the regex cannot read now
raises rather than being skipped. The old regex did not skip the new heading
quietly — it attributed 5b's `**PR:**` line to phase 5 and tripped the
pre-existing two-PR-lines assertion — but it would have gone quiet for a phase
whose heading carried no PR line, and the guard closes that.

**A trap for the next person who ticks a spec box.** The figures are generated
from §11 and the PNGs are byte-reproducible, so ticking a task without
re-running `make_progress.py` leaves them understating the work by exactly that
many tasks, silently. It happened once in this session and the review caught it.
`dev/notes/figures/corea-progress/README.md` now says so.

**Next.** 5b landed as PR #48. The fan-out can start: phases 6, 7, 8 and 10 as four
parallel PRs. Three of their done-when clauses are still not independent —
task 8.4 needs phase 6, task 8.5 composes all four ODE modules, task 10.7 checks
copy numbers against the metabolic modules' — and belong to whichever PR lands
last, or to phase 13. §11's parallelism paragraph now says so.

---

## Previous: two progress figures, parsed from the spec (2026-09-07)


**No executable change. No Julia, no Slurm run.** Six PRs had merged and it had
become hard to see what they delivered, for a structural reason: every merged
phase is framework, and the coupling figure draws the model. Two figures now
answer "what is built", in
[`dev/notes/figures/corea-progress/`](../dev/notes/figures/corea-progress/):

| Figure | What it claims |
|---|---|
| `fig1d_state_graph_progress` | **only what the reduction keeps**, laid out fresh, with every module and device badged by the spec §11 phase that writes it, and the edge-kind key gaining an **executes since** column |
| `fig2d_phase_roadmap` | all eighteen phases in dependency order, including the five that put nothing on a state graph |

**Fig 1d breaks the chain's one rule, deliberately.** Figs 1r and 1c draw the
reduction *onto* fig 1's canvas, which is right for "what did the reduction do"
and wrong for "what is built": half the ink is four deleted modules and thirty
dead couplings in grey, and the live machine has to be subtracted by eye. A
first version of fig 1d was that overlay plus badges, and it was unreadable. Fig
1d now keeps only the survivors — seven modules, three devices, two clamped
drawing devices, thirty couplings — on a layout chosen for them: two blocks,
three coupling devices, one cycle. **Fig 1c is not superseded**; it stays the
record of what the reduction cut, and the only figure that overlays the
published model.

**The one design decision worth carrying.** Fig 1c already spends a red ✗ and a
grey wash on *removed by the reduction*. Progress is a second, independent axis,
so it gets its own vocabulary — a green ✔ pill carrying the PR for merged, a
hollow ☐ pill carrying the phase number for not written — and its own key,
stating both axes side by side. A crossed-out box says nothing about whether
code exists, and a hollow badge says nothing about whether the module is in the
model; conflating the two is the misreading the figures exist to prevent.

**Two independent derivations, both asserted at build time, so neither figure
can state something the repo does not.**

1. *Status.* `progress_spec.phases()` parses `spec/spec.md` §11 — the headings,
   the `**PR:**` line, the `- [x]` ticks. Every badge, pill, title count and
   footer count comes from that parse. Only the phase → figure-element mapping
   is authored.
2. *The coupling set.* Laying fig 1d out fresh means authoring its edges, which
   is exactly how an arrow gets lost. So `live_couplings()` re-runs fig 1c's own
   edge stage and drops the removed modules and dead couplings its spec names,
   and `check_edges()` asserts fig 1d's table is **exactly** that set plus only
   what `ADDED_EDGES` declares with a reason.

Five guards, all exercised in a scratch copy: a coupling dropped from the table,
a coupling invented and not declared, a renamed `### Phase N —` heading, a
merged phase with an unticked task, and an `ELEMENT_PHASE` id the figure does
not draw. Flipping phase 6 to merged in a scratch spec turns its badge green.

**One coupling fig 1d draws that fig 1c does not, and it is the only one.**
`HOOK → ODE block`, the catalytic channel's ODE half. Fig 1 draws that channel
once, into Cofactor; Core A′ deletes Cofactor, so fig 1c greys the edge and puts
the real targets in a footnote. With the deleted modules gone there is no
footnote to hang it on, so fig 1d draws it where the channel runs — into the
block, once, which is how fig 1 labels it anyway.

Badge positions are *measured*, not guessed: the build renders once through
`neato -n -Tdot`, reads back the true box geometry (realigned against the hook,
because `-Tdot` renormalises the origin) and only then places the badges. A
clean rebuild reproduces both PNGs byte for byte.

**Two staleness items fixed, because the figures read from them.**

1. `spec/spec.md` phase 5 said `**PR:** _open_`; #46 merged as `dc5c288`. Now
   `#46 (merged 2026-09-05)`. A status field catching up, **not** a §12
   amendment — nothing the spec said became false.
2. `dev/notes/2026-09-05-reading-the-coupling-figure.md` listed phase 4 as
   *next* and said *four* of seven edge kinds execute. Both were true when
   written. Corrected to six merged phases and **six of seven** edge kinds —
   only `clamped` is still declare-only, no phase schedules it, and its held
   value still travels as a fixed parameter. The note now carries a dated
   update header and points at the two figures, which cannot go stale the way it
   did.

**What the figures say, in one line each.** Fig 1d: six of seven edge kinds
execute and not one box in the figure exists in code. Fig 2d: every merged phase
is framework, that is D0 on purpose, and phases 15–17 put nothing on a state
graph at all.

**One thing fig 1d's key now says that no earlier figure did.** The two solid
intra-CME arrows are not `CouplingEdge`s at all — spec §12 (2026-09-04) settled
that a jump module's peer writes are gated on `written_states` rather than on an
edge, and phase 2 is what made them work. The key carries that as an eighth row
below a divider, deliberately not counted among the seven kinds.

**Next.** Phase 6, central glycolysis — the first phase that writes a box.

## Previous: phase 5 — growth and volume (2026-09-05)

Spec §11 phase 5, on branch `phase-5-growth-and-volume`. **Six of the seven
edge kinds now execute. Only the clamped edge remains declare-and-validate,
its held value still travelling as a fixed parameter, and no phase schedules
it.** The volume chain is the sixth coupling channel and the one on neither
clock: every handshake the flagged membrane-protein counts set the surface
area, hence the radius, hence the volume, hence the particles-per-mM factor
that both conversion directions read — and a change in that factor dilutes
every ODE concentration at constant particle count, exactly as the published
model's `Rxns.partTomM` and `in_out.mMtoPart` do by calling
`calcCellVolume(pmap)` at every hook.

**Two protocol functions, both defaulting to empty.**

```
membrane_protein_states(m)      the module's own states whose counts set surface area
extracellular_states(m)         its own states referred to a volume other than the cell's
```

`membrane_protein_states(models)` is the sweep: a composition is *asked* for
its membrane proteins rather than the growth code knowing which sub-model to
ask, which is what `add-pts-transport`'s D7 wanted and what its Open Questions
deferred to "wave 2's call". Task 7.6 and task 11.5 already name it. **Either
block may declare it** — the toy owns `M_ptsg_c` as a jump state, task 7.6 has
`PtsTransport` own both phospho-forms in the ODE block — and a jump state's
count is read directly while an ODE state's is `u × factor` at the factor in
force at that hook, continuous rather than rounded.

**`extracellular_states` is one declaration more than the phase's tasks name,
and here is why.** Task 5.3 says growth dilutes "every ODE concentration".
True of the toy; false of the assembled model. Task 7.4 integrates external
lactate as a dynamic ODE state whose millimolar is a *medium* concentration
through D6's 1e5 ratio, and diluting it by the cell's volume ratio would
destroy lactate that has already left and open carbon balance. The registry has
no extracellular marker — `M_lac__L_e` is an ordinary dynamic species there —
so the exemption is a module declaration. Recorded in the PR and annotated on
task 5.3 in place, **not** as a §12 amendment: nothing the spec said became
false, the same call tasks 2.2 and 3.2 made.

**Measured, not asserted:**

| What | Measured |
|---|---|
| task 5.2, the radius | `sqrt(502831/4π)` = **200.03505 nm**, not 200.000 — the published area is **176 nm² larger** than an exact 200 nm sphere |
| task 5.2, the baseline | **479,563 nm²** = 502,831 − 831 × 28.0, derived at build time from the composition's own initial counts, never typed |
| task 5.3, dilution | at 831 → 1662 copies every concentration falls by **0.93440**, against `(200.035046/204.610919)³ = 0.93440` from the geometry; counts preserved to `rtol = 1e-14`; carried remainders bitwise unchanged; the exempt state untouched |
| task 5.4, the cap | volume stops at exactly `2 × V(0)` at 50,000 copies and stays there at 500,000, while the radius keeps rising |
| task 5.5, the arithmetic | area **502,831.0 → 526,099.0 nm² exactly**, radius **200.03505 → 204.61092 nm**, volume **1.07021×** |
| task 5.7, wall-clock | worst extrapolation **0.019 s** per 6,300 s trajectory (job 16226707) against K1's 10 s budget, unchanged from phase 4; the chain's own cost against a **matched control** is **+1.9% / +4.4% / +2.2%** at the three horizons against within-configuration spreads of **1.5–8.5%**, so it is **not resolved above the scatter** and is not claimed to be |

**Four things a reader should not over-read.**

1. **"Exactly 200.0 nm" is the one thing the spec got wrong, and it is a word
   rather than a number.** `4π(200 nm)²` is 502,654.8 nm², the published
   `setICs_two.py:342` area is 502,831, and `sqrt(502831/4π)` is 200.03505.
   Both notes say "returns r = 200.0 nm exactly"; that is four significant
   figures. The **area stays the primitive** — task 5.5's arithmetic closes only
   with 502,831 — and the assertion is made at the precision the source states
   it, with the full-precision value pinned beside it. Phase 3 did the same with
   20,180 particles/mM against a true 20180.39. Consequence worth knowing: **a
   growing composition's initial factor is 20,190.998, not 20,180.39**, a 0.05%
   difference, which is why passing `radius_nm` alongside a growing cell is
   refused rather than quietly ignored.
2. **A growing cell has a *larger* factor, not a smaller one.** Particles per
   millimolar scales with volume, so growth raises the factor and lowers every
   concentration at fixed count. Getting this backwards is easy and cost this
   session two test failures before the suite caught it.
3. **Dilution is a fourth operation.** The three existing conversion sites
   (catalytic write, debit, write-back) do not imply it: the ODE state is held
   in mM between hooks, so nothing rescales it unless the driver does. Each
   state is rewritten as *the count it held divided by the new factor*, in that
   order — going through the volume ratio instead would conserve the ratio and
   let the counts drift. Exactness is **one ulp, not bitwise**: `(u·f_old)/f_new`
   multiplied back by `f_new` need not return `u·f_old` in binary floating point.
4. **The numbers are the toy's, and the wall-clock one is not resolved.** Two
   ODE states diluted against Core A′'s thirty-two, one membrane protein against
   ptsG's two phospho-forms; the chain's cost scales with the number of states
   diluted. And on a shared node these timings scatter by several percent, which
   is the size of the effect — the benchmark now prints the within-configuration
   spread beside each row precisely so that a future reader checks it before
   quoting a gap. The first version of the measurement compared the growing toy
   against a row that starts at **zero** protein, and since the protein count
   fills the ODE rate law through a catalytic edge, that gap was a different ODE
   trajectory as much as a volume chain; the pre-merge review of PR #46 caught
   it, and the matched control is the fix.

**Growth may be reported as fractional growth or time-to-threshold, and as
nothing else — in code (task 5.6).** `reporting_constraints(driver)` returns
that as structured data so a composed model is asked rather than read, and
`doubling_time(driver)` throws, naming the ~92% reduction. `growth_report` and
`time_to_threshold` are the two admissible reportings, the latter read off the
recorded trajectory rather than extrapolated from a rate — a rate being one
algebraic step from the quantity being refused. Phase 14 task 14.9 consumes it.

**Fifteen refusals phases 6–12 inherit**, each asserted on its own message.
*From the sweep (2):* a flag on a state its declarer does not own; the same
species flagged twice, whose count would enter the area twice. *From the
lowering (8):* an **inbound** `VolumeEdge`; a flag with no `VolumeEdge`, which
would execute undeclared; an edge with nothing flagged, declared and never
executed; an edge on a species the module does not flag; a flagged state with no
edge of its own — the trap a module flagging both ptsG phospho-forms and edging
one would hit; a `:jump` module declaring `extracellular_states`, which has
counts and not concentrations; an exempt state the declarer does not own; and a
state both exempt and membrane-flagged, which would make the area depend on the
volume it sets. *From the build (5):* a growth keyword on a fixed cell and
`radius_nm` on a growing one; a non-positive initial area; a non-positive
footprint; and a derived baseline below zero.

**`REDUCTION_CATEGORIES` gained two**, both emitted by `driver_declarations`
whenever the chain is live: `:calibrated_constant` for the 28.0 nm² footprint —
the published model's *calibrated* value, chosen to reproduce 54% coverage for
~9,600 membrane proteins, and its own `getProtSA` docstring says 35 nm² while
its code uses 28.0, so the code governs and the contradiction is recorded — and
`:exogenous_growth` for the frozen non-ptsG baseline, which is **ours, not the
model's** and is why growth reaches ~1.07× rather than approaching the 2× cap.
§6 T2 already lists exogenous membrane growth as a row; the footprint is a row
T2 does not yet have, and gains one through this label.

**Suite:** 1381 passed, 0 failed (job 16226536); 1230 before the phase.

### Next steps

1. **Phase 5 is not merged.** The PR is open and awaits `/check-PR` and the
   user's go-ahead; the user asked explicitly for no merge without their say-so.
2. Phase 6 (central glycolysis) is the natural next one, and the ODE track of
   phases 6–9 may fan out — 6, 7 and 8 are mutually independent, 9 depends on 8.
   It is the user's call to start.
3. **Owed to phase 7, and deliberately not built here:** volume re-entering an
   ODE rate law, which task 7.2 wants for the lactate exporter's `3P/r`. The
   radius is live on the driver but reaches no rate law, `VolumeEdge` carries no
   parameter slot to write it into, and an inbound edge is refused by name
   rather than left to resolve and never run. Delivering it means an edge kind
   gains a field, which by R15 is a conversation on `main` and its own phase.
   Task 7.2 now says so in the spec.
4. When writing any module: everything in phases 3 and 4's lists still holds,
   plus — flag a membrane protein on the module that *owns* the state, one
   outbound `VolumeEdge` per flagged state; and declare `extracellular_states`
   for anything whose millimolar is not referred to the cell's volume, which for
   Core A′ is `M_lac__L_e` and nothing else.

## Previous: phase 4 — the 60 s rebuild (2026-09-05)

Spec §11 phase 4, delivered as PR #45 from branch `phase4`. **The coupling is
now bidirectional in execution, not only in topology.** Live ODE pools re-enter
the stochastic block as recomputed rate constants, held piecewise-constant
between refreshes as the published model holds them. Five of the seven edge
kinds now execute — mass and currency (phase 1), catalytic and deferred counter
(phase 3), rate constant (this phase) — and **two remain declare-only: volume
(phase 5) and clamped.** *(Phase 5 has since landed the volume chain; one
remains.)*

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
| task 4.6, dynamics | drain-aligned, paired per seed, one row per pool: on the debited pool **+0.02% ± 0.02%** (5 s) and **+0.01% ± 0.02%** (60 s); on the product pool **−0.17% ± 0.16%** and **−0.04% ± 0.17%** — *none resolved above noise*, all far below D10's 1% |
| task 4.6, sawtooth | reported separately, on the debited pool: 60 s at **+0.65% ± 0.02%** against a closed form of **+0.64%**; 5 s at **+0.06% ± 0.02%** against **+0.04%** |
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

**Suite:** 1230 passed, 0 failed (job 16222708, on the final tree); 1126 before the phase.

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
   *(Both have since landed; only the clamped edge remains.)*

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
