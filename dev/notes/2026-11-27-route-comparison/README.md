# KIM-195: the candidate routes for the full cell's block, compared on one frame

Linear: [KIM-195](https://linear.app/kimpson/issue/KIM-195/compare-the-candidate-routes-for-the-full-cells-block-what-each).
Feeds the route decision, KIM-187, and T1. The frame (the cases each route is costed
on) is in [`frame.md`](frame.md), committed at d465a35 before any timing in this issue ran.

Status: complete for KIM-195's six items. §5 says which routes stay in contention for KIM-187. The choice is KIM-187's.

## 1. The candidates: what each approximates, and where its error is measured

Five routes are compared. The first three are KIM-187's original options. The last two
were added on 2026-11-27 (Decision log).

| # | route | in one line |
|---|---|---|
| R1 | **mean field + NUTS** (the mixed route carried over) | the bulk factor's mean from a deterministic mean field of the whole cell; the transcript factor reads its rates from the mean-field pools; NUTS on the block |
| R2 | **data plug-in** | stage 1: the transcript factor with its rates taken from the observed panel's pools, so no ODE is solved; stage 2: a block-local bulk factor for metabolism over stage-1 draws (spec D7, variant 2) |
| R3 | **mean field + learned surrogate** for the part that does not close | R1, with the non-closing reactions' contribution to the means replaced by a learned map trained on SpeedyCell runs |
| R4 | **emulated synthetic likelihood** | R1's likelihood, with the bulk factor's offset δ and covariance S learned as functions of θ from SpeedyCell runs, in place of fixing them at θ_ref |
| R5 | **neural SBI on SpeedyCell** | a neural posterior or likelihood estimator trained on simulations of the full cell; no closure at all |

### 1.1 How each route handles each row of KIM-177's inventory

The row ids are those of `dev/notes/2026-11-26-kim177-nonclosure/README.md` §4. Key:
- **PM**: product of means, the naive closure. It is not exact.
- **CA**: conservation-aware treatment of a 0/1 site (O + F₁ = 1), as spec D14 describes.
- **exact**: the route makes no approximation for that row.
- **sim**: the row is handled inside SpeedyCell simulations, which are exact draws from
  the model.
- **δ(θ)**: the row's bias in the bulk mean is absorbed by the learned offset.
- **surr**: the learned surrogate stands in for the row.
- **obs**: the row's inputs are read from observed data.
- **—**: the row does not enter that route's likelihood.

| row | what it is (count) | R1 mean field | R2 data plug-in | R3 MF + surrogate | R4 emulated SL | R5 neural SBI |
|---|---|---|---|---|---|---|
| 1a | DnaA × a 0/1 site (99) | CA, or PM as the second arm (D14) | — (replication is not in stage 1) | surr | δ(θ) | sim |
| 1b | DnaA × a multi-copy site (2) | PM | — | surr | δ(θ) | sim |
| 1c | site × site within one cascade (177) | exact (these close) | — | exact | exact | sim |
| 2 | initiator × chromosome (2) | PM | — | surr | δ(θ) | sim |
| 3 | tRNA charging (62) | PM | obs, if charged tRNA is on the panel; otherwise fixed at θ_ref | surr | δ(θ) | sim |
| 4 | cost payments (20) | PM; reads nothing downstream, so it cannot change an observable | — | PM | PM | sim |
| 5 | gene copy number × transcription (493 genes, 486 reactions) | mean dosage from the mean replication state | **latent and not handled**: fixed at one copy, as KIM-185 timed it | surr | δ(θ) for bulk; the transcript factor still takes the mean dosage | sim |
| 6 | fork steps with count-dependent constants (491) | rates rebuilt from the mean state each minute | — | as R1 | as R1 | sim |
| 7 | `max(1, pool)` floor (29 nodes, 1,434 rules) | applied to the mean | applied to the observed pool | as R1 | as R1 | sim |
| 8 | own-protein promoter, `min(c·P/180, 180)` (29 rules) | applied to the mean protein count | obs, if the protein is observed; otherwise from a stage-1 pass | as R1 | as R1 | sim |
| 9 | rounded `CellV` and the volume cap (1 + 2 nodes) | continuous volume from the lipid means, rebuilt each minute | volume from the known cell number and the panel | as R1 | as R1 | sim |
| 10 | complex AND rule, `min` over subunits (9 nodes, 18 ODE reactions) | `min` over the subunit means | stage 2 only, as R1 | as R1 | as R1 | sim |
| 11 | hook rules: counter clamp, rounding, `CellSA` | dropped or made continuous (item 3 lists which) | — | as R1 | δ(θ) | sim |
| 12 | doubling time as a threshold crossing of `CellSA` | not in the likelihood. It is the phenotype that picks the block (D9), read from SpeedyCell. | — | — | — | — |

What the routes approximate beyond the inventory:
- **R1 and R3.** The bulk factor's δ and S are fixed at θ_ref (spec, "Assumptions").
  One draw of paths per θ.
- **R2.** It ignores that metabolism feeds back on the pools between saves. The stage-1
  rates use only the observed pools, with the panel's noise treated as the true state.
  On the reduced model (KIM-116) this shifted the promoters up from the joint: by +0.20
  joint SD on average (16 of 17 up), and by +0.65 on promoter 0607, which failed the
  "close" test. Against the truth the plug-in's promoters sit at −0.61 SD (14 of 17
  negative), the joint's offset with that upward shift laid over it, not a correction of
  it. Its D3 predictions of Q have about twice the joint's CRPS (median about 1% high).
- **R4.** The emulator's error in δ(θ) and S(θ). That error grows with distance from
  the training design.
- **R5.** The amortisation error of the network, and the summary statistics it is given.
  The posterior is exact only in the limit of infinite simulations and a well-specified
  estimator.

### 1.2 Where each route's error can be measured

| route | M0 (exact sampler, SBC) | reduced model (joint posterior) | M13 (exact check with copy number and second-order reactions) |
|---|---|---|---|
| R1 | yes. Done: passes SBC, about 80× cheaper; but on m0_v6 it sits 0.42–0.56 SD from the long exact reference, with no cause named (KIM-90). | yes, against the truth (there is no exact posterior): the joint's promoters sit −0.81 SD below the truth on KIM-92's dataset (14 of 17), an offset KIM-116 locates in the transcript factor. Over ten replicate truths the mean is −0.985 SD, with four replicates 6–11 SE below zero: gate 2 passes as frozen, but the route is not calibrated across truths (KIM-178, PR #210, merged as 50123cc; Decision log 2026-11-28). | **yes, its main test, not yet run**. Row 1a under CA and PM, row 3 PM, row 5 the mean dosage (D14). |
| R2 | yes, but M0 has no pool feedback to cut | **yes. Done (KIM-116):** fails "close" against the joint (+0.65 SD on one promoter, +0.20 on average); against the truth −0.61 SD, beside the joint's −0.81 (KIM-92); Q's CRPS about 2× the joint's | row 5 untested: it holds copy number fixed by construction |
| R3 | no: nothing on M0 needs a surrogate | partly: the surrogate's error on the reduced model's means | yes. The surrogate is trained where M13 is exact and is scored there. M17 is its own check. |
| R4 | **yes**: δ(θ) against the exact sampler. KIM-90's postmean arm is a one-point precursor: δ and S recomputed at the posterior mean closed 0.18 SD on S_0779, 0.16 on krnadeg and about 0.02 on ENO, of gaps of 0.42–0.56 SD. KIM-90 names no cause. Its one lead, the transcripts counted twice, is outside what R4 changes. Against that lead for krnadeg, the mixed fit already sits about 0.4 SD high with no transcripts (one fit, ESS 39; whether it is the same gap was not tested). This note reads that as pointing at the bulk factor, which a learned δ(θ) acts on (argued). On the reduced model, KIM-178 finds the offset already present at δ = S = 0, though δ and S at θ_ref shrink it on the worst replicate (−3.84 to −2.95 on 011) and the offset tracks max \|δ\| (r = −0.84). Its hypothesis, untested, puts the rest in the transcript factor, which R4 keeps. | yes: against the joint | yes |
| R5 | **yes, and needed**: SBC and the exact sampler on M0. This would be a new pilot, proposed in item 6 if R5 stays in. | yes: against the joint | yes |

### 1.3 What each route changes in the Claim, the non-goals and the positioning

| route | Claim and contribution sentence | non-goals | against PreCell (amortised SBI on the full cell, 28 group multipliers) |
|---|---|---|---|
| R1 | unchanged. "Inferred across the ODE–stochastic boundary" holds. Gate 3 (M13) decides whether the closure is good enough. | unchanged | distinct: a likelihood-based posterior on a targeted block, verified at each tier. PreCell has no exact check. |
| R2 | **changes**: it does not cross the boundary. §3.4 and the contribution sentence have to be reworded (spec risk). | unchanged | distinct, but weaker: a plug-in is not compositional inference |
| R3 | unchanged. **M17 moves ahead of M15**, because the surrogate becomes an enabler rather than §4.2's demonstration. | unchanged. §4.2's swap becomes part of the route. | distinct. The learned part is local and checked on M13. |
| R4 | unchanged in substance: the likelihood still spans the boundary, and its bulk mean is learned. "Synthetic likelihood" enters the methods. | unchanged | between the two: learned like PreCell, but local to δ and S, with a mechanistic transcript factor |
| R5 | **changes**: the mechanistic likelihood goes, and the contribution becomes amortised inference plus verification | **breaks** "no neural SBI comparison beyond the one surrogate swap" | **same ground** as PreCell. The paper becomes a comparison with it. |

## 2. The mean field on a realistic system: what a gradient costs (item 3)

**A 20-parameter gradient of the full cell's realistic mean field takes 359–426 s
serial** (job 21842557), depending on where the block sits. That is **about 280–335×
KIM-184's floor** of 1.277 s, measured on the first-order surrogate. It is also well
under the roughly 20 minutes that KIM-184 superseded. The cost is linear in the number
of passes: 177–230 s per pass of 10 partials, whatever B and wherever the block sits.

The code is `dev/scripts/fullcell_route_cost/` (`model.jl`, `route_cost_meanfield.jl`).
The runs of record (run 2) are at 491cb90, from the detached worktree
`.worktrees/kim195-run-491cb90`, with outdir
`/fred/oz022/tkimpson/infercell-kim195-2026-11-27/run2/`. The stage outputs are copied
into `results/`. All jobs ran on Milan nodes (AMD EPYC 7543), Julia 1.11.9, LinearSolve
3.87.0, with one BLAS thread. Run 1 (jobs 21785042–47, 70326e7, `run1/`) was the same
pipeline on the first clamp (§2.1); its gradients were 2.6–3.2× cheaper, and it is
superseded.

| stage | job | what |
|---|---|---|
| precompile | 21842554 | both environments |
| reference | 21842555 | SpeedyCell's own values at x₀ and at x₀ + 1, and its payments (`reference_speedycell.jl`, in KIM-183's environment) |
| check | 21842556 | the correctness checks below: **PASS** |
| time | 21842557 | serial: the primal, and every placement × B |
| library | 21842558 | LinearSolve 3.87's own dual path against DualSparseKLU, over the whole window |
| threaded | 21842559 | 4 threads, B = 20 and 40 |
| clamp | 21842560 | both payment forms at θ = 0 against 8 SpeedyCell cells (`clamp_speedycell.jl`) |

**Compute used:** 11.5 CPU-h for the runs of record, and 26.7 CPU-h over every job, run 1
and smokes 1–5 included (sacct CPUTimeRAW). The outdir holds:
- `sacct.tsv`: every job;
- `SHA256SUMS`: 153 files. The sum file's own sha256 begins `5459dced8c340314`.

The smokes, which are plumbing and not results:

| smoke | jobs | commit | what it found |
|---|---|---|---|
| 1 | 21718940–45 | 272a6f7 | world age at the first flux call |
| 2 | 21724625–30 | c60e6be | the finite-difference check on one small component |
| 3 | 21770013–18 | c454ad9 | dATP negative: no clamp on the counters |
| 4 | 21780916–21 | 70326e7 | clean |
| 5 | 21841153–54, 21841156–60 | 491cb90 | clean, on the new clamp and the new checks |

Jobs 21718943–45, 21724628–30 and 21770018 waited on failed dependencies and never ran.

### 2.1 The system

The system is the whole syn3A CME and ODE, as one deterministic ODE on all 2,661 CME
counts (`model.jl`'s header has the detail). It has:
- 2,002 first-order reactions, and 363 second-order reactions taken as the product of
  means (k·x_a·x_b; this is not exact);
- the 175 ODE fluxes, generated from `data/model/params/ode_reactions.csv` with the named
  constants as runtime values, and turned into counts per second at the cell volume;
- the hook's 26 cost counters, kept as states and paid into their pools at
  κ·c·p/(p + 1), with κ = 1 s⁻¹; the NMP counters are paid at κ·c. Where the pool holds
  many counts this pays the whole counter at rate κ, as the hook pays all of it each
  second. As the pool empties the payment falls to zero, so no pool goes negative: a pool
  smaller than its counter is held at a fraction of a count, where the hook pays it to 0;
- at each minute boundary, the 2,365 CME rate constants, the 79 enzyme concentrations,
  the volume and the surface areas are rebuilt from the mean state through the
  `nodes.csv` DAG. In a run, `round` is the identity (CellV only).

The Jacobian is analytic for the CME and payment terms. For the fluxes it is a
ForwardDiff Jacobian, nested under the outer duals. It has 24,549 stored non-zeros
(check stage), against KIM-184's 21,095.

**Declared simplifications.** These are cost choices, not fidelity claims.
- **The 1 s cadence.** SpeedyCell rebuilds the enzyme levels, the volume and the ODE
  state every 1 s, at the hook. This model rebuilds them every 60 s, with the CME rates.
- **The counter payment.** It is continuous here, not the hook's per-second
  min(counter, pool).
- **The hook's rounding.** The hook rounds the ODE state to counts each second, rounds
  `CellSA_Lip` and truncates `CellSA_Prot`, and rounds `CellV`. This model keeps all four
  continuous (row 11).
- **The second-order closure.** The product of means is not exact. M13 measures that
  error.

**Two changes from the plan, both Tom's calls (Decision log 2026-11-27 and 2026-11-28).**
1. **The counters are paid through a clamp.** In the plan they were paid at once with no
   clamp. That drove dATP to −2,887 over the window even at θ = 0, and broke the MET
   B = 20 solve in smoke 3 (job 21770016). The first clamp, κ·c·p/(c + p + 1), was
   described as a smooth min. It is not one: where counter and pool are comparable it pays
   about half of min(c, p). The review found this, and Tom chose κ·c·p/(p + 1) on
   2026-11-28. The `clamp` stage (job 21842560) runs both forms at θ = 0 against 8
   SpeedyCell cells (procid 1, seeds 1–8; `results/clamp.md`):
   - **The new form is much closer to SpeedyCell in the replication counters.** The dATP
     counter is 676, 877 and 12 at minutes 10, 20 and 40, against SpeedyCell's mean of
     1,364, 1,593 and 1,340. The old form's grows to 13,220, 37,520 and 87,450.
     These are snapshots of a quantity that oscillates. Averaged over minutes 21–52 (from
     `results/clamp/`), the dATP counter is 615 under the new form, 78,640 under the old,
     and 1,623 in SpeedyCell: the new form is about 2.6× low, the old about 48× high.
   - **The pools are close under the new form, not the old.** At minute 40 the dGTP pool is
     29,420 and dTTP 92,640, against SpeedyCell's 21,170 and 81,080; the old form holds
     dGTP at 42. Averaged over minutes 21–52 the dATP pool is 1,660 under the new form, 84
     under the old, and 1,322 in SpeedyCell.
   - **Replication is dATP-limited in the cells too.** Each of the 8 cells has dATP below
     10 counts at 16–23 of its minute-ends, and a hook pays its dATP pool to 0 in 28–38
     of its minutes. The mean field's pool is below 10 at 20 minute-ends, spread over
     minutes 5–51. So a pool near zero is the model's own behaviour, not an artefact.
   - **What differs is the onset.** The mean field's replication counters are non-zero
     from minute 1 (29 counts of dATP), while all 8 cells' are 0 at minutes 1 and 5 and
     first go non-zero at minutes 8–21. At minute 5 the mean field's dATP pool is 0.03
     counts, against SpeedyCell's mean of 25,090. That early start is most likely the
     closure of rows 1a and 2, which turns a 0/1 initiation event into a fraction from the
     start. It is argued, not measured here; M13 measures it.
   - **The cost.** The new form takes 16,487 steps over the window at θ = 0, against 6,415
     for the old, which accounts for the rise in §2.4's times over run 1's.
2. **The finite-difference gradient check is judged against the block's scale**, as
   |dual − FD| ≤ 1e-6 · max|g|. The rule was set after smoke 2's check (job 21724627)
   failed on one small component, so it is **post-hoc** (Decision log 2026-11-27). That
   component's central differences scatter as h shrinks: they are −1.32450, −1.32458,
   −1.32457, −1.32412 and −1.32503 for h = 1e-3 to 1e-5 (login node), against a dual of
   −1.32450. The run of record shows the same thing (job 21842556). For BOTH θ₁ (dual
   0.2471), the central difference is off by 2.4e-4 at h = 1e-4 and by 4.8e-6 at
   h = 1e-3. An error that shrinks as h grows is the mark of noise in the objective (here
   about 1e-8 to 5e-8, as 2h·|dual − FD|), not of a wrong derivative.

**Why it is not KIM-184's 1.3 s.** The primal takes 16,487 steps over the window,
against KIM-184's 174 (job 21842557).
- A dual pass takes 13,118–18,141 steps. KIM-184's took 144.
- The cost per step is about 12.2–14.4 ms for a dual pass, against KIM-184's 4.4 ms. The
  Jacobian is now state-dependent, and the fluxes' part is computed by AD.
- **The steps are not comparable one for one.** DiffEqBase's error norm for a dual state
  averages over the value and its N partials, n(1 + N) terms, so a dual pass does not
  take the primal's steps. The figure per step, and the comparison with KIM-184, carry
  that caveat.
- **Where the steps go.** The steps per minute at θ = 0 are in the clamp stage's
  `results/clamp/meanfield.csv` (job 21842560); the pools below are from a login-node
  diagnosis.
  - Minute 1 takes 866 steps under either clamp. Before the clamp, 365 of its 843 fell
    in its first second: a start-up transient from x₀ (login node, before the clamp).
  - Minutes 2–4 take about the same under both forms (225, 183 and 116 new; 225, 182 and
    127 old). Over minutes 5–52 the new form takes 101–557 a minute, against
    91–296 under the old, mostly about 93 (16,487 against 6,415 in all).
  - Only the dATP pool sits below 10 counts with its counter above it, at 20 minute-ends
    with a minimum of 0.027 counts. A pool held at a fraction of a count, under the
    absolute tolerance of 1e-10, is the likely cause. Since the cells also run
    dATP-limited (above), most of this cost is most likely the model's, not the
    closure's (argued).
- **At looser tolerances** (1e-6 absolute and 1e-6 relative), a GE B = 20 gradient takes
  118.4 s against 358.8 s, about 3.0× less, with 4,557 steps in its first pass (job
  21842557, one repeat). That is a lever for KIM-194, not a result: what the inference
  needs is not measured here.

### 2.2 Correctness checks (job 21842556, PASS)

| check | tolerance | observed |
|---|---|---|
| The rate rules, p and y₀ at x₀ against SpeedyCell's own `evaluate_rates!` and `hook_inputs` | 1e-12 relative | rates 2.2e-16 (2,358 of 2,365 bitwise); p 3.5e-16; y₀ 2.2e-16; CellV 335 = 335 |
| The ODE fluxes at (y₀, p) against SpeedyCell's `flux!` | 1e-12 relative | exact (0.0) |
| **The RHS at x₀ against an independent BigFloat sum of the table's propensities** (SpeedyCell's rates, both reaction types), plus the ODE part from SpeedyCell's fluxes | **1e-10 relative** | **4.05e-15** over 856 non-zero states; all 1,805 zero states are 0 |
| The surface areas at x₁ = x₀ + 1 against SpeedyCell's `lipid_sa!` | 1 (absolute) | 0.43: the hook rounds `CellSA_Lip`, this model does not; `CellSA_Prot` exact |
| The rate rules, p, y, the ODE fluxes and the volume at x₁, rebuilt from x₁, against SpeedyCell's at x₁ | 1e-12 relative | rates, p, fluxes and V exact; y 1.9e-16; CellV equal |
| **The RHS at x₁ against an independent BigFloat sum**, with the 26 payments built from SpeedyCell's own list of them | **1e-10 relative** | **1.1e-14** over 2,578 live states; no other state non-zero; all 26 payments present |
| The Jacobian at x₁ against ForwardDiff of the RHS, row by row, relative to each row's largest entry | 1e-8 | 4.0e-15; the 44 payment rows 2.4e-16; no entry off the sparsity pattern |
| J·v against central differences of the RHS, at x₀ + 1 | 1e-6 | 3.0e-10 |
| Dual gradient against central differences, B = 20, tight tolerances (1e-14 absolute, 1e-12 relative), over the first 300 s (5 rebuilds), 3 parameters per placement | 1e-6 · max\|g\| | ≤ 1.1e-7 |
| DualSparseKLU against LinearSolve 3.87's own dual path, one pass, first minute | 1e-8 | partials 6.0e-15 (GE) and 3.3e-15 (MET); the same steps |

**What the checks cover.** At x₀, 1,591 of the 2,365 reactions have a reactant at zero
count (299 of the 363 second-order ones), and every counter is 0. So the RHS check at x₀
exercises about a third of the network and none of the payments. The checks at x₁ cover
the rest: every reaction is live there, every counter is 1, and the surface, volume and
rebuild are compared with SpeedyCell's own code. J·v alone has little power on the
payment rows, since its tolerance is set by metabolic rows of about 3×10³; the row-by-row
Jacobian check covers them.

**The finite-difference check, per kind.** In BOTH, max|g| comes from the MET
parameters, so the scale rule alone barely tests BOTH's GE ones. The check also reports
each component against the largest |g| of its own kind (not part of the pass rule):
- GE: ≤ 1.0e-7 at h = 1e-4. MET: ≤ 8.7e-10.
- BOTH, the MET component: 6.2e-10.
- BOTH, the two GE components: θ₁ 5.5e-6 at h = 1e-4 and 1.1e-7 at h = 1e-3, the mark of
  noise in the objective at the smaller step; θ₁₀ 1.6e-6 and 1.4e-6, just over 1e-6.
  θ₁₀'s two differences fall on either side of the dual (43.24608): 43.24601 at
  h = 1e-4 and 43.24614 at h = 1e-3. At the smaller step, 2h·|dual − FD| is 1.4e-8, the
  same noise as θ₁'s. At the larger step the gap is about 9× what that noise would give,
  so it is most likely truncation, of order h². No step between the two was run.

**Not a check, an observation.** The rate rules at x₀ differ from `reactions.csv`'s
`rate` column in 1,435 of 2,365 reactions, by up to 10×. That column is the first minute's
rates from `MinCell_CMEODE.py`. SpeedyCell does not use it, and neither does this model:
both rebuild from the rules at every minute, the first included. KIM-184's surrogate used
the column.

### 2.3 The blocks

These are ranked by integrated flux at θ = 0 over the window. `results/check.md` lists
all 40 in each block.
- **GE.** The top 40 are all tRNA charging: the charging cascades of 8 amino acids,
  4 steps each, and their 8 cost payments. That is KIM-177's rows 3 and 4 (the payments
  are row 4), not promoters or translation. Of the 40, 32 are second-order.
- **MET.** These are kcatF constants, led by glycolysis: GAPD, PGM, ENO, LDH_L, PGK, PGI,
  TPI, FBA, PFK, PYK, and so on.
- **BOTH.** The top 10 or 20 of each.

The cost does not care much which: 177–230 s a pass in every placement. KIM-180's block,
picked on doubling time, will differ from these.

### 2.4 Seconds per gradient (job 21842557, serial; median of 3 after a 60 s warm-up)

| placement | B | passes | gradient (s) | per pass (s) | steps, pass 1 | run 1 (s) |
|---|---|---|---|---|---|---|
| GE | 10 | 1 | 180.8 | 180.8 | 14,316 | 69.3 |
| GE | **20** | 2 | **358.8** | 179.3 | 14,296 | 139.9 |
| GE | 40 | 4 | 703.6 | 177.4 | 14,288 | 274.6 |
| MET | 10 | 1 | 212.3 | 212.3 | 16,729 | 79.4 |
| MET | **20** | 2 | **426.2** | 214.3 | 17,129 | 158.9 |
| MET | 40 | 4 | 811.1 | 204.4 | 16,824 | 297.7 |
| BOTH | 10 | 1 | 229.8 | 229.8 | 18,141 | 72.9 |
| BOTH | **20** | 2 | **397.3** | 198.9 | 13,939 | 131.0 |
| BOTH | 40 | 4 | 755.6 | 189.5 | 13,118 | 264.5 |

The primal solve (θ = 0) takes a median 25.3 s over 16,487 steps. The run 1 column is the
first clamp (job 21785045), for comparison only.

### 2.5 DualSparseKLU against LinearSolve 3.87's own dual path (job 21842558)

One pass of 10 partials over the whole window, at B = 20:

| placement | LinearSolve 3.87 (`Rodas5P(linsolve = KLUFactorization())`, LinearSolve's own dual handling; this arm's cache type was not logged) | DualSparseKLU (median of 2) | ratio | steps (both) | partials, relative difference |
|---|---|---|---|---|---|
| GE | 246.4 s | 199.6 s | 1.23× | 14,296 | 7.9e-9 |
| MET | 304.0 s | 257.2 s | 1.18× | 17,129 | 5.1e-11 |

The partials differ by up to 7.9e-9 over the whole window, against 6.0e-15 over the first
minute (§2.2): the two paths' round-off grows over 14,000 steps. This comparison has no
tolerance of its own. The DualSparseKLU arm's cache type is logged (`check.md`:
`LinearCache` / `LinearSolveFunction`). The per-pass times here are on another node from
§2.4's (dave11 against dave7).

**On 3.87, the stock dual path is within 25% of DualSparseKLU.** KIM-184 measured 860× on
LinearSolve 3.73.0. So a version bump would retire DualSparseKLU and the three
LinearSolve internals it reaches into, at about 20% in cost. **But 3.87.0 does not load
on Julia 1.10.5,** which InferCell runs on. Its extensions are triggered by SparseArrays,
which it also lists as a hard dependency, and 1.10's loader fails with KeyError
"SparseArrays". This item therefore ran on Julia 1.11.9, which is also the version
SpeedyCell needs. Moving InferCell to Julia 1.11 is then the price of the stock path. It
is a decision for KIM-187 or KIM-194, not taken here.

**The stack.** Every number in §2 is on Julia 1.11.9 with LinearSolve 3.87.0, not on
InferCell's Julia 1.10.5 with LinearSolve 3.73.0. DualSparseKLU's figures there are not
measured. KIM-184's suggest they would be close, since its algebra does not depend on the
version, but that is argued.

**These are upper bounds on what this implementation costs**, at this model's fidelity
and the pinned tolerances. The fluxes' Jacobian is a dense
ForwardDiff Jacobian, nested under the outer duals: 146 inputs in 13 sweeps, about 2.2 MB
allocated a call, where its pattern has 493 non-zeros and 8 colours. The review's
micro-benchmarks put it at about 50–60% of a dual step (login node). So a coloured sparse
Jacobian could cut the cost per step by up to about 2×. That is inferred, not measured,
and it would leave the step counts as they are. SpeedyCell's `ODE.jl` already builds
such a Jacobian (`AD_BACKEND`, `FLUX_JAC_PREP`), which makes it a lever for KIM-194.
Every mean-field figure in §2, §4 and §5 (R1, R2's stage 2, R3 and R4) is an upper bound
on this account.

### 2.6 Threaded (job 21842559)

4 threads (`--threads=4 --gcthreads=4,1`), one pass per thread, median of 2 after a 60 s
warm-up. B = 10 has one pass, so threading cannot help it, and it is not rerun.

| placement | B | passes | gradient (s), 4 threads | serial (s), job 21842557 | speed-up, across jobs | median pass (s) |
|---|---|---|---|---|---|---|
| GE | 20 | 2 | **202.8** | 358.8 | 1.77× | 200.1 |
| GE | 40 | 4 | **246.2** | 703.6 | 2.86× | 240.6 |
| MET | 20 | 2 | **242.8** | 426.2 | 1.76× | 242.5 |
| MET | 40 | 4 | **275.0** | 811.1 | 2.95× | 263.6 |
| BOTH | 20 | 2 | **242.4** | 397.3 | 1.64× | 218.9 |
| BOTH | 40 | 4 | **267.8** | 755.6 | 2.82× | 235.5 |

**Threads contend, as in KIM-184.** The median pass on a shared node takes 200–264 s,
against 177–214 s alone. The serial and threaded jobs ran on different nodes (dave7 and
dave5), so the speed-ups mix the node with the threading.

## 3. The emulator and neural SBI: cost, argued

Every row here is **argued**, from the sources named in it. Nothing in this section was
run. The unit is one SpeedyCell cell simulated over the fit window (0–3,150 s, 52.5
minutes), on one Milan core.

**The cost of one simulated cell over the fit window: about 24 s (argued).** This is
from KIM-183 (job 21664418, `dev/scripts/fullcell_timing_cycle/fullcell_timing_cycle_result.md`
§3–4). A simulated minute takes 0.34 s at minute 1 (0.37 s on average over minutes
1–10), rising to 0.58 s by minute 62. So 52.5 minutes come to about 0.45 s × 53 ≈ 24 s.
That is consistent with 49.5 s to minute 105, since the cost per minute rises as the
cell grows. Throughput follows the same source:
- 100 cells on the 64 threads of one node, in one process, run a median 1.60× slower
  per cell than on one core (1.71× while all 64 threads were busy), at 1,912 cells per
  node-hour for full cycles (measured);
- one worker per core would run at about one-core speed (argued there).

So the simulation budgets below are CPU-h at one-core speed, and their wall-clock is
the budget divided by the cores used.

| | what is simulated | cells | CPU-h at 24 s a cell | wall-clock | source |
|---|---|---|---|---|---|
| **R5, PreCell's design on our window** | 2×10⁵ θ × 8 cells | 1.6×10⁶ | **≈ 10,700** | ≈ 7 days on 64 cores; ≈ 17 h on 640 | PreCell's budget: 2×10⁵ θ × 8 cells ≈ 28,000 core-h for full cycles, plus 80–150 A100-h of training, with the loss not saturated |
| R5, PreCell's design, full cycles | the same, to minute 120 | 1.6×10⁶ | ≈ 28,000 | as above × 2.6 | PreCell's own figure |
| R5, a sequential (non-amortised) estimator for one dataset | 10⁴–5×10⁴ θ × 8 cells | 8×10⁴–4×10⁵ | ≈ 530–2,700 | hours on a few nodes | argued: sequential NPE/NLE on a 20-dimensional block typically needs 10⁴–10⁵ simulations, and whether 8 cells per θ suffice is untested |
| R5, its M0 check | SBC against the exact sampler (M0's 50 datasets) | — | M0-sized, small | — | it would be the new pilot (§1.2) |
| **R4, an emulator of δ(θ) and S(θ), B = 20** | 200–1,000 design points × 100–800 cells | 2×10⁴–8×10⁵ | **≈ 130–5,300** | hours to a day, in parallel | argued: 10–50 design points per parameter is the usual Gaussian-process rule; 800 cells is the reduced model's reference ensemble (KIM-92), and 100 is M12's per-draw count |
| R4, with M12's draws reused | M12's 10³ draws × 100 cells (KIM-58) | 10⁵, to be simulated by M12, **which has not run** | **≈ 0 beyond M12's own cost**, if M12 saves the bulk panel at the 52 fit-window saves | — | applies to **θ_H drawn** only: M12 draws every parameter, so its draws are an emulator design over θ_H as well as the block. For θ_H fixed, M12's draws are off the design and do not help. |
| R4, then the posterior | R1's chains, with the emulated δ(θ) and S(θ) in the bulk factor | — | R1's cost (§5) | R1's | the emulator's own gradient is cheap beside the mean-field solve |

What the table shows:
- **R5's CPU-h is of R1's order, and it costs a different kind of thing.** Its budget is
  simulation, which parallelises without limit, plus GPU training. Its wall-clock is not
  binding. It needs an M0 pilot before its error is known (§1.2), and it puts the paper
  on PreCell's ground (§1.3).
- **R4 adds simulation to R1, and does not replace it.** It keeps R1's per-gradient cost
  and adds the emulator's training. Once M12 has run, reusing its draws would make that
  training nearly free, but only for θ_H drawn. For θ_H fixed, the case R1's chains are
  costed on, it is a fresh design of about 130–5,300 CPU-h.
- **None of these rows is a measurement.** The cost of 24 s per cell is the only measured
  input, and it is a published-path figure. KIM-183 notes that the `params=` path, which
  every θ draw uses, was not timed.

## 4. The per-gradient total, case by case (item 4)

The per-gradient total is the mean-field gradient of §2 plus the transcript factor's
gradient. Which transcript figure applies depends on the route (frame §2.2). All rows use
the pinned tolerances.

**R1, the mean-field route.** The transcript factor reads its rates from the mean-field
pools (spec D7). In the code, `route_transcription_rates` (`src/mixed_likelihood.jl`)
rebuilds every gene's rates from the route solution's states at the 53 switches. Under
θ-duals, every gene therefore carries all B partials, wherever the block sits. KIM-185's
naive pass puts width-12 duals through the whole factor, so it has the same cost
structure. This departs from the issue's wording, "the transcript factor's per-gene
gradient": per gene applies to R2, not to R1. Each
factor pass carries 10 partials, so B partials take ⌈B/10⌉ passes.
- **Serial.** The factor's term is KIM-185's full-width pass: 194.11 s at Dual width 12,
  measured over all 455 genes (job 21665221). It is applied here to width 10, which is
  **argued** and conservative.
- **Threaded.** That pass on 16 workers is taken as 194.11 / 6.82 = 28.5 s, using the 6.82×
  speed-up that KIM-185 measured for its per-gene pass at 16 workers (job 21665222). This
  is **argued**: the full-width pass was never threaded.
- The mean-field term is measured: job 21842557 serial and job 21842559 on 4 threads.
- **The totals mix two stacks.** The factor's terms (KIM-185) ran on Julia 1.10.5 with
  ForwardDiff 0.10.39; the mean-field terms on Julia 1.11.9 with LinearSolve 3.87.0
  (§2.5).

| placement | B | mean field, serial (s) | transcript factor, serial (s) | **total, serial (s)** | mean field, 4 threads (s) | factor, 16 workers (s) | **total, threaded (s)** |
|---|---|---|---|---|---|---|---|
| GE | 10 | 180.8 | 1 × 194.1 = 194.1 | **374.9** | 180.8 (1 pass: serial) | 1 × 28.5 = 28.5 | **209.3** |
| GE | 20 | 358.8 | 2 × 194.1 = 388.2 | **747.0** | 202.8 | 2 × 28.5 = 56.9 | **259.7** |
| GE | 40 | 703.6 | 4 × 194.1 = 776.4 | **1,480.0** | 246.2 | 4 × 28.5 = 113.8 | **360.0** |
| MET | 10 | 212.3 | 1 × 194.1 = 194.1 | **406.5** | 212.3 (1 pass: serial) | 1 × 28.5 = 28.5 | **240.8** |
| MET | 20 | 426.2 | 2 × 194.1 = 388.2 | **814.4** | 242.8 | 2 × 28.5 = 56.9 | **299.7** |
| MET | 40 | 811.1 | 4 × 194.1 = 776.4 | **1,587.5** | 275.0 | 4 × 28.5 = 113.8 | **388.9** |
| BOTH | 10 | 229.8 | 1 × 194.1 = 194.1 | **423.9** | 229.8 (1 pass: serial) | 1 × 28.5 = 28.5 | **258.3** |
| BOTH | 20 | 397.3 | 2 × 194.1 = 388.2 | **785.5** | 242.4 | 2 × 28.5 = 56.9 | **299.4** |
| BOTH | 40 | 755.6 | 4 × 194.1 = 776.4 | **1,532.0** | 267.8 | 4 × 28.5 = 113.8 | **381.7** |

**Serial, the two terms are about equal:** the mean field is 48–54% of R1's gradient
(359 of 747 s at GE B = 20). Threaded, the mean field is 68–89% of it (about 80% at
B = 20), since its 4 threads gain less than the factor's argued 16 workers. The factor's share rests on
KIM-185's naive width-12 pass on a stand-in factor, applied here by argument, so the
split is itself argued. Either term is a lever:
- **the transcript factor, argued and not measured.** The factor depends on θ only
  through a few pool trajectories: 4 NTPs at 52 windows, plus 29 own-protein counts. So
  its gradient could be taken with respect to those inputs once, and contracted with the
  mean field's partials. That would replace carrying B partials through all 455 genes;
- **the mean field:** looser tolerances (3.0×, §2.1) and a coloured flux Jacobian (up to
  about 2×, inferred, §2.5).

Both belong to KIM-194's levers if KIM-187 chooses R1.

**R2, the data plug-in.**
- **Stage 1.** The rates come from the observed pools, so each gene carries only its own 2
  partials. That is KIM-185's per-gene pass: 53.61 s serial (job 21665221) and 7.28 s on
  16 workers (job 21665222). It covers all 455 genes, so it does not depend on B. **It
  prices KIM-185's block, the promoters and mRNA decay (456 coordinates), not the frame's
  GE block.** The frame's GE block is tRNA charging and its cost payments (§2.3), which
  are not stage-1 parameters of the plug-in: row 3 of §1.1 reads charging from the panel
  or fixes it at θ_ref. Inferring them would need a block-local factor for charging, as
  stage 2 is for metabolism. That is not costed here.
- **Stage 2.** This is the block-local bulk factor for metabolism. Its per-gradient term
  is taken as R1's MET mean-field gradient, for example 426.2 s at B = 20 (job 21842557),
  with **no** transcript factor. That is an upper bound: a solve with gene expression
  frozen at stage-1 draws, or ODE-only, would be cheaper. It is not measured.

**R3, R4 and R5.**
- **R3 and R4.** Their per-gradient cost is R1's. The learned parts are cheap to evaluate
  beside the mean-field solve and the transcript factor. That is argued.
- **R5.** It has no gradient: its cost is simulation and training (§3).

## 5. The projection: CPU-h and wall-clock per posterior, and what stays in contention (item 6)

**Gradients per posterior.** 36,000 over 4 chains, or 9,000 a chain: the frame's figure
(frame §2.1), scaled from KIM-92's measured chains and fixed before any timing. KIM-185's
36,194 (35,282 leapfrog steps of KIM-116's stage-1 chains, plus 912 metric gradients at
456 coordinates) agrees to 0.5%. R2's stage 1 keeps KIM-185's own count. The figure is
**argued**: the full cell's conditioning is unknown.
- CPU-h is gradients × seconds per gradient.
- Wall-clock is one chain's 9,000 gradients in series. The 4 chains run in parallel.
- The threaded R1 rows reserve 16 cores per chain, which the transcript factor's term
  needs. So CPU-h there is 16 × the wall-clock, summed over the 4 chains.
- One-off costs are not in the table. They are small beside the chains: the dataset is
  200 cells × about 24 s, the reference ensemble 800 × 24 s (about 7 CPU-h with the
  dataset), plus a MAP (§3's cell cost; KIM-92's prep).

### 5.1 R1, the mean-field route (θ_H fixed)

| placement | B | CPU-h, serial | wall per chain, serial (days) | CPU-h, threaded (16 cores a chain) | wall per chain, threaded (days) | source |
|---|---|---|---|---|---|---|
| GE | 10 | 3,749 | 39.1 | 33,487 | 21.8 | §4 × 36,000 |
| GE | **20** | **7,470** | **77.8** | **41,551** | **27.1** | §4 × 36,000 |
| GE | 40 | 14,800 | 154.2 | 57,605 | 37.5 | §4 × 36,000 |
| MET | 10 | 4,065 | 42.3 | 38,530 | 25.1 | §4 × 36,000 |
| MET | **20** | **8,144** | **84.8** | **47,957** | **31.2** | §4 × 36,000 |
| MET | 40 | 15,875 | 165.4 | 62,222 | 40.5 | §4 × 36,000 |
| BOTH | 10 | 4,239 | 44.2 | 41,321 | 26.9 | §4 × 36,000 |
| BOTH | **20** | **7,855** | **81.8** | **47,898** | **31.2** | §4 × 36,000 |
| BOTH | 40 | 15,320 | 159.6 | 61,068 | 39.8 | §4 × 36,000 |

Every row multiplies a per-gradient total from §4 by the argued gradient count. The
mean-field term in each total is measured (jobs 21842557 and 21842559). The transcript
factor's term is measured at a neighbouring width (job 21665221) and threaded by argument.

What this says:
- **CPU-h is about 7,500–8,100 at B = 20, serial.** That is over the 5,000 CPU-h line
  that §5.2 applies to the plug-in, and about 70% of the 11,000 CPU-h that the project
  description scaled from the reduced model.
- **Wall-clock binds.**
  - A serial chain takes about 78–85 days, and a threaded one about 27–31 days. Both are
    well over Milan's 7-day job limit, so the chains must checkpoint, as KIM-92's did.
  - At B = 20 the threading buys about 2.6–2.9× in wall-clock at about 5.6–6.1× in
    CPU-h. The transcript factor's 16 workers are idle while the mean field runs.
- **The levers measured or argued here would bring it down:**
  - looser tolerances: 3.0× on the mean-field term (job 21842557, one repeat);
  - a coloured flux Jacobian: up to about 2× on the mean-field term (inferred, §2.5);
  - the transcript factor through its pool inputs (§4, argued);
  - fewer serial gradients (KIM-194 item 3).

  None of them has been combined or measured on a chain. The first two act on the
  mean-field solve, so R2's stage 2 (§5.2) inherits them too.
- **The figure is not settled.** Changing only the payment form, to the more faithful
  one, moved the mean-field term up by 2.6–3.2× (run 1 against run 2). Most of the extra
  steps come from the dATP pool held near zero, which the cells show too (§2.1), so most of
  it is most likely the model's cost (argued). But 3,659 of the 10,072 extra steps (36%)
  fall in minutes 5–20, the onset window read as closure, when 0–5 of the 8 cells have
  dATP below 10 (`results/clamp/`). A better closure might remove up to about a third of
  the extra cost. What tolerance the inference
  needs is not measured, and it could move the figure either way.

### 5.2 R2, the data plug-in (θ_H fixed)

| placement | B | CPU-h | wall per chain | source |
|---|---|---|---|---|
| GE (the frame's: tRNA charging) | 10 / 20 / 40 | not costed | — | not a stage-1 parameter (§4) |
| promoters and mRNA decay (KIM-185's block, not the frame's) | 456 | **539** (sensitivity 787–1,800) | 135 h serial; 18 h on 16 workers (1,171 CPU-h) | KIM-185 (jobs 21665221, 21665222, 21666966): stage 1 alone |
| MET | 10 / 20 / 40 | 2,123 / 4,262 / 8,111 | 22.1 / 44.4 / 84.5 days serial | stage 2 alone: 36,000 × R1's MET mean-field gradient (job 21842557). An upper bound, argued. |
| BOTH | 10 / 20 / 40 | 2,662 / 2,662 / 4,801 | stage 1, then stage 2 | KIM-185's 539 for the promoters in place of the frame's GE half, plus stage 2 at B/2 MET parameters (one pass up to 10) |

The plug-in's feasibility line is 5,000 CPU-h (gate 1). The promoter case and BOTH up to
B = 40 are under it; MET at B = 40 is over it, at 8,111. The frame's own GE case is not
costed. R2 does not cross the boundary, it holds gene copy number fixed (row 5), and on
the reduced model it failed KIM-116's "close" test (§1.2).

### 5.3 R3, R4 and R5

| route | CPU-h per posterior | wall-clock | θ_H | source |
|---|---|---|---|---|
| R3 (R1 + a surrogate for the non-closing part) | R1's (§5.1) + the surrogate's training (M17; not costed) | R1's | fixed | argued |
| R4 (R1 + an emulated δ(θ), S(θ)) | θ_H fixed: R1's + 130–5,300 for a fresh design. θ_H drawn: the design is M12's draws, ≈ 0 beyond M12 once M12 has run (it has not), but the chains for θ_H drawn are not costed (R1's likelihood needs a pseudo-marginal or a joint, frame §1.3) | R1's; the design runs in parallel beforehand | fixed, or drawn. Drawn, R1's likelihood must still take one θ_H a solve (frame §1.3). | §3, argued |
| R5 (neural SBI) | about 10,700 (PreCell's design on our window), or 530–2,700 (sequential), + 80–150 A100-h | hours to days: simulation parallelises | fixed or drawn, at the same cost | §3, argued |

**What is not costed.** The Done when asks for every route × case. These cells are left
open, each for the reason given:
- **R3's training.** The surrogate's design and training belong to M17 and are not
  costed. R3's row is a lower bound.
- **R4 and R5 at other B and placements.** Their costs are simulation, which does not
  depend on where the block sits. R4's design grows with B (10–50 points per parameter),
  so B = 10 and 40 would roughly halve and double it; that is argued. R5's PreCell-scale
  design is not sized by B.
- **R2's MET and BOTH, threaded.** Stage 2 was not timed on threads. R1's MET mean-field
  term on 4 threads (§2.6) is the nearest guide.
- **R1 and R2 with θ_H drawn.** Each evaluation needs one θ_H, so "drawn" means a
  pseudo-marginal or a joint over thousands of parameters. Frame §1.3 gives them a
  sentence, not a number, and so does this note.

### 5.4 Which routes stay in contention for KIM-187

This is a reading of the tables above. The choice is KIM-187's, and Tom's.

**The 5,000 CPU-h line is a yardstick here, not a gate.** It is gate 1's line for the
plug-in's gene-expression block (KIM-185); frame.md sets no pass line for KIM-195. At
B = 20 it is met by R2 (serially: its promoter case at 539, or 1,171 on 16 workers; its
metabolism stage at 4,262; its BOTH case at 2,662) and by R5 with a sequential estimator
(530–2,700). It is not met by R1 (7,500–8,100 serial), by R3 and R4, which cost R1's
chains and more, or by R5 at PreCell's design (about 10,700).

| route | in contention? | why |
|---|---|---|
| **R1, mean field + NUTS** | **yes**, over the line at B = 20 | **Cost.** Serial, 3,700–4,200 CPU-h at B = 10, under the 5,000 CPU-h line, and 7,500–8,100 at B = 20, over it. Threaded at B = 20 it is 42,000–48,000 CPU-h, for 27–31 days a chain. The figure moved up 2.6–3.2× with the payment form alone and is not settled (§5.1); the levers in §5.1 are KIM-194's. **Error.** It can be measured at every tier. On M0 it passes SBC but sits 0.42–0.56 SD from the long exact reference, with no cause named (KIM-90). On the reduced model its promoters sit −0.81 SD from the truth (KIM-92), and across ten truths it is not calibrated (KIM-178). M13, its main test, has not run. **Claim.** It is the only route that keeps the Claim as written with no learned component. |
| R2, data plug-in | as a fallback only | **Cost.** It is the cheapest: 539 CPU-h serial for KIM-185's promoter block (135 h a chain), or 1,171 CPU-h at 18 h a chain on 16 workers. That block is not the frame's. Its metabolism stage costs 2,100–8,100 CPU-h serial (22–85 days a chain), over the line at B = 40. The frame's own GE case is not costed. **Error.** On the reduced model it is shifted from the joint by +0.20 joint SD on average and +0.65 on one promoter, so it fails "close" (KIM-116). Against the truth it sits at −0.61 SD, beside the joint's −0.81: the same offset with an upward shift, so neither route is accurate there. On dQ, the knockdown's effect, it scores the same as the joint (CRPS 1.848e4 against 1.857e4). It is untested across truths, and its gene copy number is untested. **Claim.** It changes the Claim and the contribution sentence. |
| **R3, mean field + surrogate** | **yes**, conditional on M17 | It costs R1's chains plus training. It does not make R1 cheaper. Its case is accuracy where the closure fails, measured on M13. It keeps R1's transcript factor, so it inherits whatever of KIM-178's offset lies there. |
| **R4, emulated synthetic likelihood** | **yes** | It costs R1's chains plus a design: 130–5,300 CPU-h for θ_H fixed, or M12's draws for θ_H drawn once M12 has run. It targets the bulk factor's δ and S fixed at θ_ref. On M0 that is not shown to be the gap's cause: KIM-90 names none, recomputing δ and S at the posterior mean closed at most 0.18 SD of gaps of 0.42–0.56, and KIM-90's metabolites-only fit, on this note's reading, points at the bulk factor for krnadeg (argued; one fit, ESS 39). On the reduced model KIM-178 finds the offset already present at δ = S = 0, though δ and S shrink it on the worst replicate and it tracks max \|δ\| (r = −0.84). Whether a learned δ(θ) would remove it is untested: KIM-178's hypothesis, also untested, puts it in the transcript factor, which R4 keeps. The chains for θ_H drawn are not costed. Its error is checkable on M0. |
| R5, neural SBI | **only if Tom reopens the non-goal** | Its CPU-h is comparable or lower, and it has no closure error. But it breaks the non-goal "no neural SBI comparison beyond the one surrogate swap", puts the paper on PreCell's ground, and has no M0 check yet. |

**The M0 pilot (proposed, not filed).** R4 stays in contention, so the item calls for its M0
pilot against the exact sampler as a separate issue. The proposal:
- learn δ(θ) and S(θ) on M0 from a design of simulations;
- rerun the mixed route with them;
- compare the result with the exact sampler's posterior on M0's datasets, by the same
  shift and SBC tests M0 used.

KIM-90's postmean arm is a one-point precursor there (0.16–0.18 SD closed on krnadeg and S_0779, about 0.02 on ENO). The cost is M0-sized:
simulations at about M0's cost per cell, plus M0's 6.6 CPU-h per mixed-route posterior.
**This is offered for Tom's approval before it is created** (`Future:` or M11/M13, his call).
An SBI pilot is not proposed unless the non-goal is reopened.

## Reproduce

```bash
# from a detached worktree at a pushed commit (the launcher refuses a dirty tree)
dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir>          # precompile, reference, check, time, library, threaded, clamp
SMOKE=1 dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir>  # plumbing pass: a 300 s window, 2 repeats
dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir> time     # one stage (check writes ranking.jls first)
```

Julia 1.11.9 is at `/fred/oz022/tkimpson/julia-1.11.9`. The reference and clamp stages
run SpeedyCell in KIM-183's environment and depot, where it is developed by path. Their
inputs are the syn3a-ensemble repo at 2075809cbb825140cb31a7bec8bc3529770c8e13, read
only.

What a stranger needs, and what they will get:
- **syn3a-ensemble is private**, and its path is fixed in `model.jl`. Without access the
  model cannot be built. The mean-field stages (check, time, library, threaded and the
  mean-field half of clamp) assert that the repo is at the commit above with
  `data/model` unmodified, and record it beside the InferCell commit in their headers.
  The two SpeedyCell-side scripts record its commit but assert nothing, and nothing
  checks `code/SpeedyCell` for local edits.
- **The environments are not instantiated by the launcher.** Both depots were built on a
  login node, since the compute nodes have no network. On a new machine, run
  `Pkg.instantiate()` in `dev/scripts/fullcell_route_cost` (default depot) and in
  `dev/scripts/fullcell_timing_cycle` (the depot the launcher's `jlref` sets) first,
  with Julia 1.11.9, or edit `jl` and `jlref` in the launcher. The SpeedyCell
  environment's Manifest pins SpeedyCell by absolute path.
- **The timings reproduce only to run-to-run variation.** They depend on the node and on
  what else shares it (§2.6). The step counts did repeat: 16,487 at θ = 0 in jobs
  21842557, 21842559 and 21842560, on three nodes. Residuals at round-off level may differ by machine.
- The outdir's `dev/`, `dev-loginnode*/` and `review/` hold login-node diagnoses and the
  review's scratch. They are not runs of record, and they are outside `SHA256SUMS`.
