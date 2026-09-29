# Spec: Phase 16 — Recovery, coverage, and the reference

**Status:** draft
**Created:** 2026-09-30  ·  **Last amended:** —
**Parent spec:** [`spec/spec.md`](../spec.md) §11 phase 16. That spec's §2 (background), §7
(non-goals) and §8 (kill criteria) are inherited, not restated. Where this file and
the parent disagree on phase 16, this file wins. On anything else, the parent wins.

This sub-spec was spawned by `/implement-spec` on 2026-09-30. Phase 16 as written
was unclear rather than merely long, for three reasons:
- Task 16.2 is an open design decision that every later task depends on.
- The sampler 16.1 needs does not exist.
- The compute sizing for K1, K2 and coverage was not derived.

The phase splits in two. **16a** builds the sampler and validates it against the
reference. **16b** runs recovery on Core A′ at whatever scale 16a's measured cost
affords. Each is one PR. The parent's tasks 16.1 to 16.11 map as follows. Every
one is kept in the parent as an annotation.

| Parent | Here | | Parent | Here |
|---|---|---|---|---|
| 16.1 | 16a.9 | | 16.7 | 16b.4 |
| 16.2 | 16a.7 | | 16.8 | 16a.10 |
| 16.3 | 16a.6 | | 16.9 | 16b.5 |
| 16.4 | 16b.1 | | 16.10 | 16b.7 |
| 16.5 | 16b.2 | | 16.11 | 16b.8 |
| 16.6 | 16b.3 | | — | 16a.1–5, 16a.8, 16a.11, 16b.6 (new) |

---

## 1. The claim

**The three-block conditional scheme of D10 samples the joint posterior of D11's
six targets, and that posterior is calibrated.** The data are 200 cells of the
published model, with exact transcript counts and lognormal metabolite panels. The
sampler is validated against an exact long-run reference on a two-gene system (M0).
So a calibration result can be attributed to the method rather than to a converged
wrong answer (R4). On Core A′ it delivers F6's two load-bearing cells, the tight-prior
control, K2's verdict and coverage. Each is delivered at a scale set by a measured
cost against a fixed budget, and the scale is stated beside every number. **The
claim is only as large as the rung it was measured on.** A coverage result on M0 is
reported as M0's, never as Core A′'s.

## 2. Background

What phase 16 starts from, all on main at `49994d3`:

- **The data.** The 15.7 dataset (`dev/data/corea_dataset_15/meta.md`, SHA-256
  `5e6382da…`):
  - 200 published-model cells (clamped drain, fractional carry), seeds 30001 to
    30200, over a 6,300 s cycle.
  - Truth drawn with `Xoshiro(1507)`, Haldane-consistently.
  - 17 metabolites observed lognormally at σ = 0.1, floored at one particle, from
    `Xoshiro(1508)`.
  - 17 transcripts as exact counts.
  - 105 save points from 60 s to 6,300 s. t = 0 is excluded.
- **Identifiability.** K6 does not fire: rank 6 of 6 against a split-half floor
  (15.8, job 17704779).
- **The observation model.** `Modality`, `NoiseModel` and `observation_loglik`
  (`src/likelihoods.jl`). The `:poisson` modality only labels the transcripts, and it
  is not how they are scored (D16.7).
- **The simulator.** The hybrid driver (`src/handshake.jl`, `HandshakeDriver`,
  `handshake_step!`, `set_parameters!`) runs the SSA and the ODE together.
  - **Nothing drives it from a recorded jump path.** Phase 14b's `replay` in
    `dev/scripts/corea_validation_14b.jl` is a linear-block script helper, not this.
  - Cost: 32.2 s warm per 6,300 s cycle, 5.11 ms per handshake (13.7).
- **Priors.** Every target is `LogNormal`. Promoters and `krnadeg` are at
  `log(2.0)` width (ours, asserted). ENO is at gstd 1.170 and FBA at 1.747.
  - `draw_parameters(...; reverse = :derived)` derives each reverse constant through
    the nominal equilibrium constant (14c.1).
- **Rate laws the conditionals rest on.**
  - Transcription's constant is
    `min(RNAPOL_KCAT·S_g, ceiling) / denom_g(NTP pools)`, rebuilt every 60 s.
  - Decay's propensity is `krnadeg / n_g · m_g`.
  - The ceiling never binds at nominal (8.85 against 180, D11).
- **Samplers on disk.** Turing's `AdvancedPS`, `SSMProblems` and `Libtask` are in the
  Manifest as transitive dependencies. They are not direct dependencies.
  `src/inference.jl` runs NUTS and ABC only on uniform (ODE-only or jump-only)
  compositions, so there is no hybrid sampler.
- **The smoothed model is not the sampler here.** 14c built the "sampler model"
  (smoothed drain, `:continuous` pools) for a gradient step. 14c.5's coupled half
  failed on GTP (3.75% against 1%) because of the continuous pools. D16.1 takes the
  gradient out of phase 16, so this model is no longer on the sampling path.

## 3. Scientific validity

### Model: the target posterior

Write `c` for a cell, `X_c` for its jump path, `θ = (θ_ODE, θ_CME)` and `σ` for the
metabolite noise scale. Then:

```
p(θ, σ, X_1:C | y)  ∝  p(θ) p(σ) · Π_c  p(X_c | θ) · 1[ m(X_c, t_k) = y^tx_c,k  ∀k ]
                                        · Π_k Π_j  N( log y^met_c,j,k ; log max(x_j(θ, X_c, t_k), 1), σ )
```

- **`X_c`** is every transcription, translation, decay and translocation event with
  its time.
- **`x(θ, X_c, t)`** is the published model's ODE state, replayed deterministically
  from the path. Fractional carry is deterministic, so a path and θ fix the whole
  trajectory.
- **`p(X_c | θ)`** is the jump process's path density: `Π_events λ_e(t_e⁻) ·
  exp(−∫ Σ λ dt)`. Every propensity is zeroth or first order, and the constants are
  piecewise-constant between 60 s rebuilds. So the integral is an exact finite sum.
  **`θ_ODE` enters it**, because the rebuild reads the ODE pools.
- **The transcript likelihood is an indicator.** Transcripts are observed exactly, so
  the path must reproduce every observed count. This is D10's exact conditional, and
  it is not a Poisson modality.
- **The one-particle floor** matches how the data were made (§12 2026-09-29).

### The three blocks (D10, amended by D16.1)

| Block | Conditional | Cost per update |
|---|---|---|
| B1 | `θ_ODE, σ \| X, y`: prior, plus Σ_c [log p(X_c \| θ) + metabolite log-likelihood]. A new θ_ODE needs every cell replayed. σ alone needs no replay | C replays per θ_ODE evaluation |
| B2 | `θ_CME \| X, θ_ODE`: each promoter and `krnadeg` enters only its own propensities. With the pool trajectory fixed, the likelihood is `k^n exp(−k·A)`, where `n` counts events and `A` is exposure summed over cells. The prior is lognormal, so this is not gamma. It is a one-dimensional log-concave density, sampled exactly | Negligible: sufficient statistics only |
| B3 | `X_c \| θ, y_c`, per cell and independent across cells | The path update of D16.3 |

### Assumptions, each with its regime

| Assumption | Holds when | What breaks outside it |
|---|---|---|
| A path and θ determine the ODE trajectory exactly | Fractional carry and the clamped drain are deterministic | Stochastic rounding would add a latent. It is not the data policy. Asserted by C2 |
| Rebuild instants coincide with save points (every 60 s) | Published cadence | A window's constants would straddle observations. Asserted in 16a.1 |
| The turnover ceiling never binds | `RNAPOL_KCAT·S_g < 180` at every state visited | B2 becomes a truncated rather than a plain `k^n e^{−kA}` likelihood. Asserted per proposal, and a violation throws |
| Cells are independent given θ | Separate seeds and separate cells, as generated | Nothing here: true by construction |
| The initial transcript count follows the generator's rule | The rule the dataset was built with. Rounded mean or Poisson draw is checked in 16a.1 | If drawn, the t = 0 count is a latent that B3 must sample |
| Block 1's step needs no gradient | Two free ODE constants plus σ | Freeing many more ODE parameters would bring the gradient, and 14c's sampler model, back. §7 |

### Correctness checks, each with a tolerance

| # | Check | Asserted | Where |
|---|---|---|---|
| C1 | Path density | `E_{X~p(·\|θ0)}[p(X\|θ1)/p(X\|θ0)] = 1` within 3 SE, over 10⁴ M0 paths at two θ pairs (one CME-only shift, one ODE-only shift). A one-gene, frozen-rate path matches a hand computation to 1e-12. **Mutation:** dropping the exposure term, or reading rate constants from the previous window, fails it | 16a.2 |
| C2 | Replay fidelity | Replaying a generated cell's recorded path at its θ reproduces the latent trajectory **bitwise** at every handshake and all 105 save points, on 3 dataset cells. If the replay changes the solver's step sequence, the fallback is ≤ 1e-8 relative, and the reason is recorded. **Mutation:** removing one translation event changes GTP | 16a.1 |
| C3 | Exact transcript-only posterior | One gene, rate constants frozen (NTP pools chemostatted), transcripts observed every 60 s and nothing else. The exact posterior of `(S, krnadeg)` on a 2-D grid comes from the birth-death transition matrix, truncated at count 30 with truncation mass < 1e-12. B2 with B3's bridge must match it: the 5, 25, 50, 75 and 95% marginal quantiles each within 3 Monte Carlo SE (ESS-based) | 16a.4 |
| C4 | Block 1 on a fixed path | On M0, 3 cells, path fixed: B1's samples of `(ln kcat_ENO, ln kcat_FBA)` match a 2-D grid posterior computed by direct replay, with quantiles within 3 MC SE. `assert_haldane` holds at every accepted state | 16a.5 |
| C5 | Path update invariance | With the metabolite likelihood switched off, the path kernel must leave the bridge-conditioned prior invariant. Per-gene, per-window event-count histograms over 10⁴ sweeps match exact FFBS marginals (χ², p > 0.01, with the multiplicity stated). The particle-state independence test: mutating one particle's snapshot leaves every other unchanged | 16a.7 |
| C6 | Reference convergence | At least 4 chains from overdispersed starts. Rank-normalised R̂ < 1.01, and bulk and tail ESS ≥ 1,000 per parameter. The Monte Carlo SE of every reported 5% and 95% quantile is ≤ 0.05 posterior SD. This is what "the run's length justified rather than chosen" means | 16a.9 |
| C7 | Reference calibration | Simulation-based calibration of the M0 sampler at production settings. Rank-ECDF within 95% simultaneous bands for every free parameter, at R replications (≥ 50, K1's floor; 100 if 16a's cap allows) | 16a.9 |
| C8 | F10's own floor | The divergence between the production and reference posteriors is reported beside the reference-vs-reference divergence of its split chain halves. A production divergence inside that floor is "not resolved", not "zero" | 16a.10 |

## 4. Approach

**D16.1. Block 1 is gradient-free on the published model. This amends D10's
"gradient-based sampler".**
**Why:** Block 1 has three dimensions: `ln kcat_ENO`, `ln kcat_FBA` and `ln σ`. A slice
or adaptive Metropolis step at three dimensions costs about as many replays as a
gradient step would. It samples the data-generating model exactly: clamped drain,
fractional carry, with no smoothing and no continuous pools. That removes 14c.5's
3.75% GTP departure from phase 16 as a data-and-sampler mismatch. It also sidesteps
K5 ("gradient-based sampling is invalid as posed") rather than working around it.
- **The reverse constants are derived at every proposal**, through the nominal
  equilibrium constant, as the truth was drawn (`draw_parameters(...; reverse =
  :derived)`). Otherwise the chain samples a different model from the one that made
  the data.
- **σ is updated by its own one-dimensional slice step** on the cached replay, with
  no new replay.

**Alternatives considered:**
- NUTS or HMC on 14c's sampler model, as D10 had it. That needs automatic
  differentiation through 6,300 handshakes (`rate_constants` builds an
  `SVector{n, Float64}` and would need work). It also carries 14c.5's GTP gap into
  every posterior.
- Building both and choosing on M0. This was declined: it doubles 16a for a
  three-dimensional block.

14c's sampler model stays in `src` for when the ODE block grows (§7). Task 16a.6
(parent 16.3) still runs, as the record of why the gradient was dropped.

**D16.2. Block 2 is exact, from sufficient statistics.**
**Why:** Given the path and the replayed pools, gene `g`'s transcription events give
`S^n_g · exp(−S·A_g)`, where `A_g = Σ_c ∫ RNAPOL_KCAT/denom_g dt`. Decay gives the
same form in `krnadeg`, with `A = Σ_c Σ_g ∫ m_g/n_g dt`. The lognormal prior makes
each one a one-dimensional log-concave density. A univariate slice sampler on the log
scale samples it exactly, with no replay. **The conjugate direction is the degenerate
one** (D10, D11). With the polymerase constant fixed, that costs nothing here.
**Alternatives considered:** a gamma prior, which would give an exact conjugate draw.
It was rejected because it changes the prior the truth was drawn from.

**D16.3. The path update is conditional SMC over 60 s windows, and 16a.7 chooses
its variant by measurement.**

The parent's three options, examined against the model as built:

| Option | Status |
|---|---|
| No augmentation | **Closed.** The 60 s aggregate drain fails granularity at +288% (13.4) |
| Exact forward filtering over a truncated state space | **Not available for the joint path.** The ODE state is continuous and a function of the whole path, so no finite-state filter carries it. It survives as a component: each gene's transcript path inside a window is a birth-death bridge between two exactly observed counts, and it is sampled exactly by FFBS on a truncated count space (uniformisation, cap 30) |
| Particle Gibbs | **The route.** Per cell, a conditional SMC over the 105 windows. Each particle carries the full driver state: ODE vector, jump counts, counters, rounding carries and rate constants. The proposal is the exact transcript bridge plus translation and translocation from their conditional prior given the bridge. The weight is the window's metabolite likelihood, and the bridge's transcript density cancels against the path prior's, since the bridge is exact |

**What 16a.7 decides among,** at N ∈ {5, 10, 20, 50} particles:
- **PG without ancestor sampling.** O(N·T) replays per sweep. It is prone to path
  degeneracy in early windows.
- **PGAS.** The state at a window boundary is Markov, but the reference path's future
  continuous state changes with the ancestor. So the ancestor weight needs the
  reference's future replayed from each candidate: O(N·T²).
- **PGAS truncated at lag L.** O(N·T·L).

**The criterion:**
- Every window's segment changes in at least 10% of sweeps. That is the update rate,
  and the threshold is ⚠️ DRAFT.
- Cost per sweep, measured on M0 and on one full-scale cell.
- The choice is recorded against 13.4, D13 and 14c.7's gradient report, as the parent
  requires.

**The two installed traps (parent D10), met by construction:**
- The CSMC is hand-written over explicit state snapshots. It does not go through
  Turing's `@model` or Libtask.
- Ancestor sampling does not depend on `SSMProblems`.
- A solver object cannot be shared across particles, because each particle restores a
  copied state vector into one integrator. C5's independence test asserts this.

**Alternatives considered:**
- Single-site Metropolis-within-Gibbs on one gene's window segment. Every proposal
  needs a replay from that window to the end, about half a cycle, so a sweep costs
  about 900 cycles per cell.
- A whole-path independence proposal. Its acceptance collapses under 105 metabolite
  panels.
- AdvancedPS on `SSMProblems`. It is possible, but the driver's mutable integrator
  state needs the same snapshotting anyway, and it adds a direct dependency.

**D16.4. M0 is Core A′ cut to two genes over 600 s.** ⚠️ DRAFT: the genes are
ptsG (`JCVISYN3A_0779`, the headline, the only membrane protein) and GAPD
(`JCVISYN3A_0607`, the best-determined promoter).
- The other 15 proteins are held at their proteomics counts as fixed enzymes, which is
  a labelled M0 reduction in `reduction_report`.
- **Free:** `S_0607`, `S_0779`, `krnadeg`, `kcatF_R_ENO`, `kcatF_R_FBA` and σ. That
  is D11's set less `S_PGI`, so every block runs as it does in production.
- **Cells per M0 dataset:** ⚠️ DRAFT, 50. The number is set by 16a's cost cap, not by
  D11's precision argument, and it is stated beside every M0 result.

**Why:** It is small enough that a very long run is defensible, per D10. It still
crosses the boundary in both directions and runs all three blocks.
**Alternatives considered:** one gene. That cannot exercise B2 across genes, and C3
already covers the one-gene case exactly.

**D16.5. "Production" and "reference" are the same sampler at different settings.**
- The reference runs at N_ref ≫ N particles and chains far longer than C6 requires.
- Production runs at the settings 16a.7 picks for Core A′.
- **F10 therefore measures what production's particle count, variant (for example
  truncated AS) and chain length cost.** It does not measure a structural
  approximation, because the production sampler is exact in the limit. That is
  narrower than "exact versus approximate", and F10's caption says so.

**D16.6. Measure, then allocate, inside 50k CPU-h for the whole phase.**
- **16a is capped at 5k CPU-h.**
- **Its last task (16a.11)** does four things:
  - measures the cost of one posterior on M0, and of one sweep on a full-scale cell;
  - projects a full-scale posterior from those;
  - sets K1's bound, with T3 recording that it was set after the measurement;
  - writes a costed 16b allocation inside what remains.
- **The allocation follows a fixed priority order.** ⚠️ DRAFT, and yours to reorder:
  1. 16b.1, the tight-prior control. This is K4's stop gate.
  2. 16b.2, the six-parameter recovery.
  3. F6's two load-bearing cells: ptsG under metabolites only, and ENO under
     transcripts only.
  4. 16b.3, coverage at ≥ 50 replications.
  5. 16b.5, K2 at 1,000 cells.
  6. 16b.6, K7's two fits.
  7. F6's remaining cells.
- **What doesn't fit at full scale drops down K1's ladder:** shorter horizon, then five
  genes, then M0. If nothing affordable remains, T3 scores the item "not affordable at
  50k CPU-h". Each drop is a §12 entry, and **16a.11 is a stop point.** 16b does not
  start without your approval of the allocation.

**Why:** The projection before any measurement is about 53k CPU-h per full-scale
posterior: 2,000 sweeps × (10 particles + 5 block-1 evaluations) × 200 cells × 32.2 s.
16b as written needs at least eight posteriors plus coverage. Projections have been
wrong in this spec before (D10's 5 s drain; K1's 10 s proxy), so the allocation rests
on the measurement.

**D16.7. The data model, as the dataset was made.**
- Transcripts are scored by the path indicator, not as Poisson.
- The metabolite modality has `floor = 1.0`.
- t = 0 is excluded.
- ⚠️ DRAFT: **`M_lac__L_e` is kept.** Its scale rests on the asserted
  medium-to-cell volume ratio. That ratio is fixed at the same value in the generator
  and the sampler, so synthetic recovery is unaffected. It is flagged in T2 as an
  observable that would not transfer to Step 2 unchanged.

**D16.8. K2 gets its 1,000 cells.** The 15.7 truth is extended by 800 cells at seeds
30201 to 31000, with the same model, truth and noise seed rule. The combined set is
used only for 16b.5's transcript-only fit. It stays subject to D16.6's allocation.

**D16.9. K7 gets a task (16b.6).**
- `kcatF_R_PYK` is freed, with its reverse constant derived.
- It is fitted twice on the joint data: once as built, and once with the dropped
  reactions of the PYK gene (`JCVISYN3A_0221`) stubbed in as a competing sink.
- **The reaction list is enumerated from upstream at `db048ac` first.** The parent
  says "seven", and that is checked, not assumed.
- The sink's form is ⚠️ DRAFT: each dropped reaction takes the gene's protein pool
  with its published rate law and its substrates chemostatted, so it competes for
  enzyme and not for metabolites. A labelled reduction.
- **K7 fires on a posterior width ratio above 2.**

## 5. Data and reproducibility

- **Datasets:**
  - the 15.7 dataset;
  - M0 datasets (16a.8);
  - the 1,000-cell extension (16b.5);
  - coverage datasets (16b.3), one truth per replication, drawn by `draw_truth(...;
    purpose = :coverage)` from recorded seeds.

  Each follows 15.7's convention: truth, seeds, noise seed, T2 row, `reduction_report`
  and Julia version in a tracked `meta.md`. Observations are regenerable and ignored,
  with a SHA-256.
- **Chains** are saved with their sampler seeds, settings and commit. Every quoted
  number cites a Slurm job id.
- **The ledger.** `dev/scripts/phase16_ledger.md` records CPU-h per job against the
  50k budget and 16a's 5k cap, and it is updated with every run of record.
- **Runs of record go in a worktree at a pushed commit.** The branch is never switched
  under a running job.
- **The environment is the parent's §5.** No new direct dependency is expected
  (D16.3). If one is added, it is pinned and it is in the depot before any Slurm job,
  since compute nodes have no network.

## 6. Outputs

| Output | From | Claim |
|---|---|---|
| The path-update decision record, with the measured update rates and costs | 16a.7 | Makes D10's open choice |
| C3's exact-versus-sampled figure | 16a.4 | The exact conditional is exact |
| M0's rank-ECDF (F8 on M0) | 16a.9 | The reference is calibrated |
| **F10** | 16a.10 | C4: makes every calibration claim attributable |
| F14, and K1's bound | 16a.11 | K1 |
| T4, and the 16b.1 perturbation plot | 16b.1, 16b.2 | C3, C4 |
| F9 coverage curve, and F8 where affordable | 16b.3 | C4 |
| **F6** | 16b.4 | C3: **the figure the work is for** |
| K2's verdict sentence | 16b.5 | C3 |
| K7's width ratio | 16b.6 | Scope honesty |
| The caveats block, and **T3** | 16b.7, 16b.8 | Scope honesty, falsifiability |

Every output carries its rung (M0, reduced, or full scale) and its cell count in the
caption.

## 7. Non-goals

- **Phase 17's conditional checks** (D12 parts 1 to 3, and K3). Phase 16 runs the
  joint sampler and coverage only.
- **A gradient step in block 1**, and any further validation of 14c's sampler model.
- **Freeing parameters beyond D11's six**, except `kcatF_R_PYK` for K7 alone.
- **The polymerase ridge as a sampling target.** F13 is 15.8's.
- **Any full-scale result that D16.6's allocation did not fund.** It is scored as not
  affordable, never extrapolated.

## 8. Kill criteria

The parent's K1, K2, K4 and K7 are scored here. K6 was scored in 15.8. K3 is phase
17's. There are two local stops.

- **KR — the reference is wrong.** C3, C5 or C7 fails. Stop: nothing downstream is
  attributable, which is the point of R4.
- **KP — the path update cannot mix.** At the largest N that fits 16a's cap, some
  window's update rate stays below the D16.3 threshold for every variant. Stop and
  amend. The fallbacks are a shorter window, a guided translation proposal, or K1's
  ladder.

## 9. Open questions

- [NEEDS CLARIFICATION: is a replay cheaper than a simulation?] 13.7's 5.11 ms per
  handshake includes the SSA. Replay drops it but keeps the stiff solve. 16a.1
  measures this, and it moves every projection in D16.6.
- [NEEDS CLARIFICATION: how degenerate are the particle weights?] GTP's 200-cell
  resolution of 4.5% (14c.5) implies a per-cell coefficient of variation near 64%, so
  metabolite panels discriminate paths strongly. Whether 50 particles suffice is
  16a.7's measurement.
- [NEEDS CLARIFICATION: is the initial transcript count fixed or drawn in 15.7's
  generator?] If drawn, it is a latent in B3 (16a.1).
- [NEEDS CLARIFICATION: which of PYK's upstream reactions are the "seven"?]
  Enumerated in 16b.6.
- The update-rate threshold, M0's genes and cell count, D16.6's priority order,
  `M_lac__L_e`, and K7's sink form are the ⚠️ DRAFT items above.

## 10. Risks

| # | Risk | Early warning | Response |
|---|---|---|---|
| R16.1 | Cost exceeds the budget at every useful rung | 16a.11's projection | K1's ladder, by D16.6's written rule. T3 publishes it |
| R16.2 | Path degeneracy in early windows | Update rate in windows 1 to 10 falling as T grows, at fixed N | PGAS, or truncated AS; KP |
| R16.3 | Replay is not bitwise | C2 | Record the reason, and use the 1e-8 fallback. A non-deterministic replay would mean a hidden latent, which is a stop |
| R16.4 | A reverse constant is sampled freely and breaks Haldane | `assert_haldane` in C4 | D16.1's derived reverse constants |
| R16.5 | F10 reads zero because the reference is too short | C8's floor | Lengthen the reference. Never report a divergence below its floor as agreement |
| R16.6 | 16a overruns its 5k cap on pilots | The ledger | Stop, and report at 16a.11 with what was measured |

## 11. Task list

### Phase 16a — The sampler and the reference

**Goal:** build the three-block sampler and validate it exactly where that is possible
and against a long reference where it is not. Then measure what one posterior costs.
**Done when:**
- C1 to C8 pass.
- The path-update variant is chosen and recorded with its measurements.
- The M0 reference has converged (C6) and is calibrated (C7).
- F10 is reported against its floor.
- K1's bound is set, and a costed 16b allocation inside 50k CPU-h is approved.

**PR:** _not started_

- [ ] 16a.1 Drive the published hybrid from a recorded jump path. Record the path
  during generation, assert that rebuilds sit on save points, and find out whether
  15.7 drew initial transcripts. — verify by C2 on 3 dataset cells (bitwise at every
  handshake), by its mutation, and by the replay's wall-clock per cycle against
  13.7's 32.2 s, reported.
- [ ] 16a.2 Implement `log p(X | θ)` from the replayed rate constants — verify by C1,
  by the 1e-12 hand check, and by both mutations failing.
- [ ] 16a.3 Implement the exact per-gene transcript bridge (FFBS by uniformisation,
  cap 30) — verify by sampled bridges always hitting both endpoints, and by their
  event-count distribution matching the exact marginal (χ², p > 0.01, 10⁵ draws per
  test case, multiplicity stated).
- [ ] 16a.4 Implement block 2 from sufficient statistics, with the ceiling asserted per
  proposal — verify by C3.
- [ ] 16a.5 Implement block 1: gradient-free, with reverse constants derived and σ
  updated separately — verify by C4.
- [ ] 16a.6 (parent 16.3) Confirm a clipped drain under a fixed path still has a
  discontinuous derivative in an ODE parameter — verify by a test, as the parent's
  annotation specifies:
  - on the clipped model, a finite-difference slope in `ln kcat` that jumps across a
    clip and does not shrink with the spacing;
  - paired with 14c.5's continuity test on the smoothed model;
  - check 7's census cited as the source of the clip, not assumed.

  Recorded as the reason for D16.1.
- [ ] 16a.7 (parent 16.2) Build the conditional SMC path update over snapshots, and
  choose PG, PGAS or truncated PGAS. Verify by:
  - C5;
  - per-window update rates and cost per sweep, on M0 and on one full-scale cell, for
    each variant at N ∈ {5, 10, 20, 50};
  - the choice recorded in §4 against 13.4, D13 and 14c.7's gradient report, with the
    two installed traps addressed as D16.3 states.
- [ ] 16a.8 Build M0 and generate its datasets — verify by:
  - the build passing completeness mode;
  - `reduction_report` naming the 15 held enzymes;
  - each dataset carrying its truth, seeds and report;
  - two runs at one seed being identical.
- [ ] 16a.9 (parent 16.1) Run the M0 reference long — verify by C6, which justifies
  the length, and by C7 at production settings on R ≥ 50 datasets, with R stated.
- [ ] 16a.10 (parent 16.8) Compare production against the reference on M0 — verify
  by F10: overlaid marginals, one 2-D contour for the ptsG and ENO pair, and a
  per-marginal divergence beside C8's floor. The caption states D16.5's narrower
  meaning.
- [ ] 16a.11 Measure cost and allocate 16b. Verify by:
  - the ledger showing 16a ≤ 5k CPU-h;
  - the cost per posterior on M0, and a full-scale projection from a measured sweep;
  - K1's bound written into parent §8, with a note that it was set after the
    measurement;
  - an allocation table funding each 16b task at a named rung, summing to ≤ 50k minus
    16a's spend;
  - a §12 entry;
  - **your approval before 16b starts.**

### Phase 16b — Recovery on Core A′, at the rungs 16a.11 funds

**Goal:** the parent's recovery, coverage and shrinkage results, each at a stated
rung.
**Done when:** the parent's phase 16 condition holds at the funded rung:
- coverage within Monte Carlo error of nominal for every target;
- the tight-prior control shrinks and recovers;
- the loose-prior control is visibly wider.

In addition, K2, K4 and K7 are scored, T3 is published, and every unfunded item is
scored as such.
**PR:** _not started_

- [ ] 16b.1 (parent 16.4) Recover the tight-prior control alone — verify by:
  - the posterior concentrating on truth well inside its prior;
  - a deliberate perturbation of the truth moving the posterior with it;
  - shrinkage below 0.9.

  **If it does not recover, the machinery is wrong and the phase stops (K4).**
- [ ] 16b.2 (parent 16.5) Recover the full six-parameter set — verify by every
  marginal containing truth, and by the cost per posterior reported against 16a.11's
  projection.
- [ ] 16b.3 (parent 16.6) Compute coverage over repeated datasets — verify by:
  - 50, 80, 90, 95 and 99% intervals covering truth at those rates within binomial
    error, per parameter;
  - the replicate count and rung stated;
  - the tight-prior control inside K4's band, and K4's third trigger (equal widths for
    ENO and FBA) checked.
- [ ] 16b.4 (parent 16.7) Produce F6 — verify by all six targets under
  transcripts-only, metabolites-only and joint data, with the two load-bearing cells
  present: ptsG under metabolites only and ENO under transcripts only. Any cell 16a.11
  did not fund is shown as not run.
- [ ] 16b.5 (parent 16.9) Generate the 1,000-cell extension (D16.8) and decide K2.
  Verify by:
  - the extension's `meta.md`;
  - the prior-to-posterior divergence under transcript-only data, for every metabolic
    target, against 0.05 nat;
  - the verdict written as a sentence, with its consequence for framing.
- [ ] 16b.6 Measure K7 (D16.9). Verify by:
  - PYK's dropped sibling reactions enumerated from upstream at `db048ac`, with the
    count checked against the parent's "seven";
  - the two fits' posterior width ratio on `kcatF_R_PYK`;
  - K7 scored against twofold.
- [ ] 16b.7 (parent 16.10) Label what the result does and does not license — verify
  by the report carrying:
  - the scoping note's caveats: best-measured region, enzyme competition partly
    removed, architecture evidence rather than well-posedness;
  - K7's factor;
  - every result's rung.
- [ ] 16b.8 (parent 16.11) Fill T3 — verify by all seven criteria carrying a
  threshold, a measured value and a verdict, **published whether or not everything
  passed**, with K1's threshold marked as set after its measurement.

## 12. Amendment log

_No amendments yet._
