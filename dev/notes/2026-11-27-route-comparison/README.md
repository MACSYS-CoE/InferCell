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
| R1 | yes. Done: passes SBC, about 80× cheaper; but on m0_v6 it sits 0.42–0.56 SD from the long exact reference, with no cause named (KIM-90). | yes, against the truth (there is no exact posterior): the joint's promoters sit −0.81 SD below the truth on KIM-92's dataset (14 of 17), an offset KIM-116 locates in the transcript factor. Over ten replicate truths the mean is −0.985 SD, with four replicates 6–11 SE below zero: gate 2 passes as frozen, but the route is not calibrated across truths (KIM-178, PR #210 open; Decision log 2026-11-28). | **yes, its main test, not yet run**. Row 1a under CA and PM, row 3 PM, row 5 the mean dosage (D14). |
| R2 | yes, but M0 has no pool feedback to cut | **yes. Done (KIM-116):** fails "close" against the joint (+0.65 SD on one promoter, +0.20 on average); against the truth −0.61 SD, beside the joint's −0.81 (KIM-92); Q's CRPS about 2× the joint's | row 5 untested: it holds copy number fixed by construction |
| R3 | no: nothing on M0 needs a surrogate | partly: the surrogate's error on the reduced model's means | yes. The surrogate is trained where M13 is exact and is scored there. M17 is its own check. |
| R4 | **yes**: δ(θ) against the exact sampler. KIM-90's postmean arm is a one-point precursor: δ and S recomputed at the posterior mean closed about 0.17 SD of the 0.42–0.56 SD M0 gap. KIM-90 names no cause, and its open lead (the transcripts counted twice) is outside what R4 changes. | yes: against the joint | yes |
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

**A 20-parameter gradient of the full cell's realistic mean field takes 131–159 s
serial** (job 21785045), depending on where the block sits. That is **about 100–125×
KIM-184's floor** of 1.277 s, measured on the first-order surrogate. It is also well
under the roughly 20 minutes that KIM-184 superseded. The cost is linear in the number
of passes: 65–80 s per pass of 10 partials, whatever B and wherever the block sits.

The code is `dev/scripts/fullcell_route_cost/` (`model.jl`, `route_cost_meanfield.jl`).
The runs of record are at 70326e7, from the detached worktree
`.worktrees/kim195-run-70326e7`, with outdir
`/fred/oz022/tkimpson/infercell-kim195-2026-11-27/run1/`. The stage outputs are copied
into `results/`. All jobs ran on Milan nodes (AMD EPYC 7543), Julia 1.11.9, LinearSolve
3.87.0, with one BLAS thread.

| stage | job | what |
|---|---|---|
| precompile | 21785042 | both environments |
| reference | 21785043 | SpeedyCell's own values at x₀ (`reference_speedycell.jl`, in KIM-183's environment) |
| check | 21785044 | the correctness checks below: **PASS** |
| time | 21785045 | serial: the primal, and every placement × B |
| library | 21785046 | LinearSolve 3.87's own dual path against DualSparseKLU, over the whole window |
| threaded | 21785047 | 4 threads, B = 20 and 40 |

**Compute used:** 5.6 CPU-h for the runs of record, and 11.5 CPU-h over every job,
smokes 1–4 included (sacct CPUTimeRAW). The outdir holds:
- `sacct.tsv`: every job;
- `SHA256SUMS`: 83 files. The sum file's own sha256 begins `3a3f3313336e4b9a`.

The smokes, which are plumbing and not results:

| smoke | jobs | commit | what it found |
|---|---|---|---|
| 1 | 21718940–45 | 272a6f7 | world age at the first flux call |
| 2 | 21724625–30 | c60e6be | the finite-difference check on one small component |
| 3 | 21770013–18 | c454ad9 | dATP negative: no clamp on the counters |
| 4 | 21780916–21 | 70326e7 | clean |

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
   - **The pools agree after minute 20 under the new form, not the old.** At minute 40 the
     dGTP pool is 29,420 and dTTP 92,640, against SpeedyCell's 21,170 and 81,080. The old
     form holds dGTP at 42.
   - **Neither form tracks dATP over minutes 5–20.** The mean field's dATP pool is 0.03–0.17
     counts under the new form and 80–94 under the old, against SpeedyCell's 4,165–25,090.
     The mean field's replication counters are non-zero from minute 1 (29 counts of dATP),
     while all 8 cells' are 0 at minutes 1 and 5. So the mean field begins replication
     earlier than any of the cells and drains dATP. That is most likely the closure of
     rows 1a and 2, which turns a 0/1 initiation event into a fraction from the start. It
     is argued, not measured here; M13 measures it.
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

**Why it is not KIM-184's 1.3 s.** The primal takes 6,415 steps over the window, against
KIM-184's 174 (job 21785045).
- A dual pass takes 5,155–6,460 steps. KIM-184's took 144.
- The cost per step is about 13 ms for a dual pass, against KIM-184's 4.4 ms. The
  Jacobian is now state-dependent, and the fluxes' part is computed by AD.
- **Where the steps go is not measured on the runs of record.** A login-node diagnosis,
  on the model before the counter clamp and in minute 1 only (which starts from x₀, not
  after a rebuild), suggests two causes:
  - **a start-up transient.** 365 of minute 1's 843 steps fell in its first second;
  - **the pinned absolute tolerance of 1e-10 on counts.** The binding states included
    fractional mean-field states of the DnaA cascade (`ssDNABoundSite_*`, about 0.01
    counts) and metabolite pools moving by about 10% a second.

  The runs of record give about 100 steps for each later minute, on average (6,415 over
  53 segments, less minute 1's). Whether those later minutes carry a rebuild transient
  too is unmeasured.
- **At looser tolerances** (1e-6 absolute and 1e-6 relative), a GE B = 20 gradient takes
  49.2 s against 139.9 s, about 2.8× less (job 21785045, one repeat). That is a lever for
  KIM-194, not a result: what the inference needs is not measured here.

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
  noise in the objective at the smaller step; θ₁₀ 1.6e-6 and 1.4e-6. θ₁₀'s gap is the same
  at both steps, so it is neither noise nor truncation. Its cause is not found. It is
  1.4e-6 of the component, just over 1e-6.

**Not a check, an observation.** The rate rules at x₀ differ from `reactions.csv`'s
`rate` column in 1,435 of 2,365 reactions, by up to 9×. That column is the first minute's
rates from `MinCell_CMEODE.py`. SpeedyCell does not use it, and neither does this model:
both rebuild from the rules at every minute, the first included. KIM-184's surrogate used
the column.

### 2.3 The blocks

These are ranked by integrated flux at θ = 0 over the window. `results/check.md` lists
all 40 in each block.
- **GE.** The top 40 are all tRNA charging: the charging cascades of 8 amino acids,
  4 steps each, and their 8 cost payments. That is KIM-177's row 3, not promoters or
  translation. Of the 40, 32 are second-order.
- **MET.** These are kcatF constants, led by glycolysis: GAPD, PGM, ENO, LDH_L, PGK, PGI,
  TPI, FBA, PFK, PYK, and so on.
- **BOTH.** The top 10 or 20 of each.

The cost does not care which: 65–80 s a pass in every placement. KIM-180's block,
picked on doubling time, will differ from these.

### 2.4 Seconds per gradient (job 21785045, serial; median of 3 after a 60 s warm-up)

| placement | B | passes | gradient (s) | per pass (s) | steps, pass 1 |
|---|---|---|---|---|---|
| GE | 10 | 1 | 69.3 | 69.3 | 5,327 |
| GE | **20** | 2 | **139.9** | 69.7 | 5,325 |
| GE | 40 | 4 | 274.6 | 68.4 | 5,326 |
| MET | 10 | 1 | 79.4 | 79.4 | 6,159 |
| MET | **20** | 2 | **158.9** | 79.1 | 6,348 |
| MET | 40 | 4 | 297.7 | 74.8 | 6,460 |
| BOTH | 10 | 1 | 72.9 | 72.9 | 6,084 |
| BOTH | **20** | 2 | **131.0** | 65.4 | 5,155 |
| BOTH | 40 | 4 | 264.5 | 64.8 | 5,228 |

The primal solve (θ = 0) takes a median 10.5 s over 6,415 steps.

### 2.5 DualSparseKLU against LinearSolve 3.87's own dual path (job 21785046)

One pass of 10 partials over the whole window, at B = 20:

| placement | LinearSolve 3.87 (`Rodas5P(linsolve = KLUFactorization())`, LinearSolve's own dual handling; its cache type was not logged) | DualSparseKLU (median of 2) | ratio | steps (both) | partials, relative difference |
|---|---|---|---|---|---|
| GE | 83.3 s | 69.1 s | 1.20× | 5,325 | 1.4e-8 |
| MET | 103.6 s | 83.4 s | 1.24× | 6,348 | 1.2e-12 |

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

### 2.6 Threaded (job 21785047)

4 threads (`--threads=4 --gcthreads=4,1`), one pass per thread, median of 2 after a 60 s
warm-up. B = 10 has one pass, so threading cannot help it, and it is not rerun.

| placement | B | passes | gradient (s), 4 threads | serial (s), job 21785045 | speed-up, across jobs | median pass (s) |
|---|---|---|---|---|---|---|
| GE | 20 | 2 | **77.0** | 139.9 | 1.82× | 76.6 |
| GE | 40 | 4 | **92.4** | 274.6 | 2.97× | 90.5 |
| MET | 20 | 2 | **89.6** | 158.9 | 1.77× | 87.9 |
| MET | 40 | 4 | **108.5** | 297.7 | 2.74× | 99.2 |
| BOTH | 20 | 2 | **85.9** | 131.0 | 1.53× | 79.7 |
| BOTH | 40 | 4 | **104.1** | 264.5 | 2.54× | 94.1 |

**Threads contend, as in KIM-184.** The median pass on a shared node takes 77–99 s, against 65–79 s
alone. The serial and threaded jobs ran on different nodes, so the speed-ups mix the node
with the threading.

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
- 64 cells on one node in one threaded process run 1.6–1.7× slower per cell (1,912
  cells per node-hour for full cycles, measured);
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
- **R5 is affordable in CPU-h but costs a different kind of thing.** Its budget is
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
- The mean-field term is measured: job 21785045 serial and job 21785047 on 4 threads.

| placement | B | mean field, serial (s) | transcript factor, serial (s) | **total, serial (s)** | mean field, 4 threads (s) | factor, 16 workers (s) | **total, threaded (s)** |
|---|---|---|---|---|---|---|---|
| GE | 10 | 69.3 | 1 × 194.1 = 194.1 | **263.4** | 69.3 (1 pass: serial) | 1 × 28.5 = 28.5 | **97.8** |
| GE | 20 | 139.9 | 2 × 194.1 = 388.2 | **528.2** | 77.0 | 2 × 28.5 = 56.9 | **133.9** |
| GE | 40 | 274.6 | 4 × 194.1 = 776.4 | **1,051.0** | 92.4 | 4 × 28.5 = 113.8 | **206.3** |
| MET | 10 | 79.4 | 1 × 194.1 = 194.1 | **273.5** | 79.4 (1 pass: serial) | 1 × 28.5 = 28.5 | **107.8** |
| MET | 20 | 158.9 | 2 × 194.1 = 388.2 | **547.1** | 89.6 | 2 × 28.5 = 56.9 | **146.5** |
| MET | 40 | 297.7 | 4 × 194.1 = 776.4 | **1,074.2** | 108.5 | 4 × 28.5 = 113.8 | **222.3** |
| BOTH | 10 | 72.9 | 1 × 194.1 = 194.1 | **267.0** | 72.9 (1 pass: serial) | 1 × 28.5 = 28.5 | **101.3** |
| BOTH | 20 | 131.0 | 2 × 194.1 = 388.2 | **519.2** | 85.9 | 2 × 28.5 = 56.9 | **142.8** |
| BOTH | 40 | 264.5 | 4 × 194.1 = 776.4 | **1,040.9** | 104.1 | 4 × 28.5 = 113.8 | **217.9** |

**The transcript factor, not the ODE, dominates R1's gradient:** 388 of 528 s at GE B = 20,
serial. That gives a lever, **argued and not measured**. The factor depends on θ only
through a few pool trajectories: 4 NTPs at 52 windows, plus 29 own-protein counts. So its
gradient could be taken with respect to those inputs once, and contracted with the mean
field's partials. That would replace carrying B partials through all 455 genes. It
belongs to KIM-194's levers if KIM-187 chooses R1.

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
  is taken as R1's MET mean-field gradient, for example 158.9 s at B = 20 (job 21785045),
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
- Wall-clock is one chain's 9,048 gradients in series. The 4 chains run in parallel.
- The threaded R1 rows reserve 16 cores per chain, which the transcript factor's term
  needs. So CPU-h there is 16 × the wall-clock, summed over the 4 chains.
- One-off costs are not in the table. They are small beside the chains: the dataset is
  200 cells × about 24 s, the reference ensemble 800 × 24 s (about 6 CPU-h with the
  dataset), plus a MAP (§3's cell cost; KIM-92's prep).

### 5.1 R1, the mean-field route (θ_H fixed)

| placement | B | CPU-h, serial | wall per chain, serial (days) | CPU-h, threaded (16 cores a chain) | wall per chain, threaded (days) | source |
|---|---|---|---|---|---|---|
| GE | 10 | 2,634 | 27.4 | 15,644 | 10.2 | §4 × 36,000 |
| GE | **20** | **5,282** | **55.0** | **21,428** | **14.0** | §4 × 36,000 |
| GE | 40 | 10,510 | 109.5 | 33,002 | 21.5 | §4 × 36,000 |
| MET | 10 | 2,735 | 28.5 | 17,253 | 11.2 | §4 × 36,000 |
| MET | **20** | **5,471** | **57.0** | **23,439** | **15.3** | §4 × 36,000 |
| MET | 40 | 10,742 | 111.9 | 35,571 | 23.2 | §4 × 36,000 |
| BOTH | 10 | 2,670 | 27.8 | 16,213 | 10.6 | §4 × 36,000 |
| BOTH | **20** | **5,192** | **54.1** | **22,852** | **14.9** | §4 × 36,000 |
| BOTH | 40 | 10,409 | 108.4 | 34,872 | 22.7 | §4 × 36,000 |

Every row multiplies a per-gradient total from §4 by the argued gradient count. The
mean-field term in each total is measured (jobs 21785045 and 21785047). The transcript
factor's term is measured at a neighbouring width (job 21665221) and threaded by argument.

What this says:
- **CPU-h is about 5,200–5,500 at B = 20, serial.** That is about half the 11,000 CPU-h
  that the project description scaled from the reduced model.
- **Wall-clock binds.**
  - A serial chain takes about 55 days, and a threaded one about 14–15 days. Both are
    over Milan's 7-day job limit, so the chains must checkpoint, as KIM-92's did.
  - The threading buys about 4× in wall-clock at about 4× in CPU-h. The transcript
    factor's 16 workers are idle while the mean field runs.
- **The levers measured or argued here would bring it down:**
  - looser tolerances: 2.8× on the mean-field term (job 21785045, one repeat);
  - the transcript factor through its pool inputs (§4, argued);
  - fewer serial gradients (KIM-194 item 3).

  None of them has been combined or measured on a chain.

### 5.2 R2, the data plug-in (θ_H fixed)

| placement | B | CPU-h | wall per chain | source |
|---|---|---|---|---|
| GE (the frame's: tRNA charging) | 10 / 20 / 40 | not costed | — | not a stage-1 parameter (§4) |
| promoters and mRNA decay (KIM-185's block, not the frame's) | 456 | **539** (sensitivity 787–1,800) | 135 h serial; 18 h on 16 workers | KIM-185 (jobs 21665221, 21665222, 21666966): stage 1 alone |
| MET | 10 / 20 / 40 | 794 / 1,589 / 2,977 | 8.3 / 16.6 / 31.0 days serial | stage 2 alone: 36,000 × R1's MET mean-field gradient (job 21785045). An upper bound, argued. |
| BOTH | 10 / 20 / 40 | 1,333 / 1,333 / 2,128 | stage 1, then stage 2 | KIM-185's 539 for the promoters in place of the frame's GE half, plus stage 2 at B/2 MET parameters (one pass up to 10) |

Every costed case is under the plug-in's feasibility line of 5,000 CPU-h (gate 1). The
frame's own GE case is not costed. But R2 does not
cross the boundary, it holds gene copy number fixed (row 5), and on the reduced model it
failed KIM-116's "close" test (§1.2).

### 5.3 R3, R4 and R5

| route | CPU-h per posterior | wall-clock | θ_H | source |
|---|---|---|---|---|
| R3 (R1 + a surrogate for the non-closing part) | R1's (§5.1) + the surrogate's training (M17; not costed) | R1's | fixed | argued |
| R4 (R1 + an emulated δ(θ), S(θ)) | θ_H fixed: R1's + 130–5,300 for a fresh design. θ_H drawn: R1's + M12's draws as the design, ≈ 0 beyond M12, once M12 has run (it has not) | R1's; the design runs in parallel beforehand | fixed, or drawn. Drawn, R1's likelihood must still take one θ_H a solve (frame §1.3). | §3, argued |
| R5 (neural SBI) | about 10,700 (PreCell's design on our window), or 530–2,700 (sequential), + 80–150 A100-h | hours to days: simulation parallelises | fixed or drawn, at the same cost | §3, argued |

### 5.4 Which routes stay in contention for KIM-187

This is a reading of the tables above. The choice is KIM-187's, and Tom's.

| route | in contention? | why |
|---|---|---|
| **R1, mean field + NUTS** | **yes** | Its CPU-h is affordable: about 5,300 at B = 20. Its wall-clock binds at about 14–15 days a chain threaded, and the levers in §5.1 are what KIM-194 would work on. Its error is measured at every tier, with M13 testing what the full cell adds (§1.2). It is the only route that keeps the Claim as written with no learned component. |
| R2, data plug-in | as a fallback only | It is the cheapest (539–2,977 CPU-h), but it changes the Claim and the contribution sentence. Its known bias (+0.65 SD on the reduced model) and its untested gene copy number count against it. |
| **R3, mean field + surrogate** | **yes**, conditional on M17 | It costs R1's chains plus training. It does not make R1 cheaper. Its case is accuracy where the closure fails, measured on M13. |
| **R4, emulated synthetic likelihood** | **yes** | It costs R1's chains plus a design: 130–5,300 CPU-h for θ_H fixed, or M12's draws for θ_H drawn once M12 has run. It targets the bulk factor's δ and S fixed at θ_ref. On M0 that is not shown to be the gap's cause: KIM-90 names none, and recomputing δ and S at the posterior mean closed about 0.17 of 0.42–0.56 SD. Its error is checkable on M0. |
| R5, neural SBI | **only if Tom reopens the non-goal** | Its CPU-h is comparable or lower, and it has no closure error. But it breaks the non-goal "no neural SBI comparison beyond the one surrogate swap", puts the paper on PreCell's ground, and has no M0 check yet. |

**The M0 pilot (proposed, not filed).** R4 stays in contention, so the item calls for its M0
pilot against the exact sampler as a separate issue. The proposal:
- learn δ(θ) and S(θ) on M0 from a design of simulations;
- rerun the mixed route with them;
- compare the result with the exact sampler's posterior on M0's datasets, by the same
  shift and SBC tests M0 used.

KIM-90's postmean arm is a one-point precursor there (about 0.17 SD of the gap closed). The cost is M0-sized:
simulations at about M0's cost per cell, plus M0's 6.6 CPU-h per mixed-route posterior.
**This is offered for Tom's approval before it is created** (`Future:` or M11/M13, his call).
An SBI pilot is not proposed unless the non-goal is reopened.

## Reproduce

```bash
# from a detached worktree at a pushed commit (the launcher refuses a dirty tree)
dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir>          # precompile, reference, check, time, library, threaded
SMOKE=1 dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir>  # plumbing pass: a 300 s window, 2 repeats
dev/scripts/fullcell_route_cost/route_cost_meanfield.sh <outdir> time     # one stage (check writes ranking.jls first)
```

Julia 1.11.9 is at `/fred/oz022/tkimpson/julia-1.11.9`. The reference stage runs in
KIM-183's environment and depot, where SpeedyCell is developed by path. Its inputs are the
syn3a-ensemble repo at 2075809cbb825140cb31a7bec8bc3529770c8e13, read only.
