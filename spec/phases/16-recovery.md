# Spec: Phase 16 — Recovery, coverage, and the reference

**Status:** met at M0; full scale budget-bound (2026-10-28). 16a is ticked in §11, and 16b closed in #105, #107, #109 and #110. 16b.5 and 16b.6 are budget-bound (§12 "T3 at M0")
**Created:** 2026-09-30  ·  **Last amended:** 2026-10-29
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
| 16.1 | 16a.8 (build), 16a.9 (run) | | 16.7 | 16b.4 |
| 16.2 | 16a.7 | | 16.8 | 16a.10 |
| 16.3 | 16a.6 | | 16.9 | 16b.5 |
| 16.4 | 16b.1 | | 16.10 | 16b.7 |
| 16.5 | 16b.2 | | 16.11 | 16b.8 |
| 16.6 | 16b.3 | | — | 16a.1–5, 16a.11, 16b.6 (new) |

**Labels.** In this file, V1 to V9 and V2b are its own verification checks. C1 to C5 are
the parent's sub-claims (§6).

---

## 1. The claim

**The three-block conditional scheme of D10 samples the joint posterior of D11's
six targets, and that posterior is calibrated.** The data are 200 cells of the
published model, with exact transcript counts and bulk metabolite measurements, one per pool per save, of the population mean (amended 2026-10-05, §12). The
sampler is validated against an exact long-run reference on a two-gene system (M0).
So a calibration result can be attributed to the method rather than to a converged
wrong answer (R4). On Core A′ it delivers F6's load-bearing cell (ptsG under
metabolites only) and its measured null (ENO under transcripts only, at 15.7's
truth; amended 2026-10-05, §12), the tight-prior control, K2's verdict and coverage. Each is delivered at a scale set by a measured
cost, and the scale is stated beside every number. **The claim is only as large as
the rung it was measured on.** A coverage result on M0 is reported as M0's, never as
Core A′'s.

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
  - **Initial transcripts are fixed, not drawn.** `corea_models` builds
    `CoreATranscription` with `seed = nothing`, so each starts at
    `round(Int, mean_mrna)` (`transcription.jl:403–406`, `assembly.jl:78`). The only
    other randomness in the driver is `:stochastic` rounding, which is not the data
    policy.
- **Identifiability.** K6 does not fire: rank 6 of 6 against a split-half floor
  (15.8, job 17704779).
- **The observation model.** `Modality`, `NoiseModel` and `observation_loglik`
  (`src/likelihoods.jl`). The `:poisson` modality only labels the transcripts, and it
  is not how they are scored (D16.7).
- **The simulator.** The hybrid driver (`src/handshake.jl`, `HandshakeDriver`,
  `handshake_step!`, `set_parameters!`) runs the SSA and the ODE together.
  - **Nothing drives the full hybrid from a recorded jump path.** The nearest thing is
    `replay` (`test/corea_langevin_doubles.jl:360`), a Langevin-block test double,
    which is not this.
  - Cost: 32.2 s warm per 6,300 s cycle, 5.11 ms per handshake (13.7).
- **The order inside one handshake**, which the path update must respect
  (`handshake_step!`, `handshake.jl:1764–1908`):
  1. growth;
  2. geometry and catalytic writes;
  3. the ODE step (t − 1 → t);
  4. debit of the deferred counters;
  5. at `step % 60 == 0`, the rebuild, reading the ODE at t;
  6. the jump step over [t − 1, t].

  **So on the jump clock the rate constants change at 60m − 1 s, not at 60m.** The
  constants for [60m − 1, 60m + 59) come from pools that already reflect window m's
  jump events up to 60m − 1, but not its last second.
- **Priors.** Every target is `LogNormal`. Promoters and `krnadeg` are at
  `log(2.0)` width (ours, asserted). ENO is at gstd 1.170 and FBA at 1.747. σ's is
  `LogNormal(log 0.2, 1.0)` (`likelihoods.jl`).
  - `draw_parameters(...; reverse = :derived)` derives each reverse constant through
    the nominal equilibrium constant (14c.1).
  - `d11_models()` also frees `kcatR_R_ENO` and `kcatR_R_FBA`, which carry their own
    priors. D16.1 says why those priors are not used.
- **Rate laws the conditionals rest on.**
  - Transcription's constant is
    `min(RNAPOL_KCAT·S_g, 180) / denom_g(NTP pools)`. The ceiling binds only above
    S ≈ 153, far outside any prior.
  - Decay's propensity is `krnadeg / n_g · m_g`.
  - Translation (`k_g·m_g`) and translocation are first order.
- **Samplers on disk.** Turing's `AdvancedPS`, `SSMProblems` and `Libtask` are in the
  Manifest as transitive dependencies. They are not direct dependencies.
  `src/inference.jl` runs NUTS and ABC only on uniform (ODE-only or jump-only)
  compositions, so there is no hybrid sampler.
- **The gradient report.** The resolver's gradient obstructions are
  `obstructs_gradients` and `check_gradient_safety`. On the assembled sampler model
  they read none (14c.4, "seven before").
- **The smoothed model is not the sampler here.** 14c built the "sampler model"
  (smoothed drain, `:continuous` pools) for a gradient step. 14c.5's coupled half
  failed on GTP (3.75% against 1%) because of the continuous pools. D16.1 takes the
  gradient out of phase 16, so this model is no longer on the sampling path.

## 3. Scientific validity

### Model: the target posterior

Write `c` for a cell. The free parameters are:

```
θ_CME = (S_0607, S_0445, S_0779, krnadeg)
θ_ODE = (kcatF_R_ENO, kcatF_R_FBA)
σ_b   = the bulk metabolite noise scale (σ, the per-cell panel's, before 2026-10-05)
```

The reverse constants are deterministic functions of θ_ODE (D16.1). The target is

```
p(θ, σ_b, X_1:C | y)  ∝  p(θ_CME) p(θ_ODE) p(σ_b) · Π_c  p(X_c | θ) · 1[ m(X_c, t_k) = y^tx_c,k  ∀k ]
                                          · Π_k Π_j  N( log z_j,k ; log max( x̄_j,k , 1 ), σ_b )

x̄_j,k = (1/C) Σ_c x_j(θ, X_c, t_k)
```

**Amended 2026-10-05 (§12):**
- The metabolites are **bulk**: one measurement `z_j,k` per pool per save, of the
  C-cell mean particle count, at an assay CV `σ_b`. This replaces the per-cell panel
  `y^met_c,j,k` at σ.
- The cells averaged are the same C cells at every save, and the mean is of counts,
  not concentrations. Both are stated idealisations.
- **The cells are no longer independent given θ.** They couple through `x̄`.

- **The priors.**
  - `p(θ_ODE)` is the product of the two forward-constant priors only.
  - `p(σ_b) = LogNormal(log 0.2, 1.0)`, the prior σ had. 15.7's truth has
    σ_b = 0.10.
- **`X_c`** is every transcription, translation, decay and translocation event with
  its time.
- **`x_j(θ, X_c, t)`** is what the generator observes: the **carry-inclusive
  particle count at the current volume**, `u_j·factor + remainder_j` (the
  `emit_observables!` / `observe_latent` path, `observables.jl:250, 430`). It comes
  from replaying the published model deterministically from the path. Near the
  one-particle floor the carry is worth several σ, so leaving it out is not a
  rounding choice.
- **`p(X_c | θ)`** is the jump process's path density: `Π_events λ_e(t_e⁻) ·
  exp(−∫ Σ λ dt)`. Every propensity is zeroth or first order, and the constants are
  piecewise-constant on the jump clock's intervals [60m − 1, 60m + 59). So the
  integral is an exact finite sum. **`θ_ODE` enters it**, because both rebuilds
  (transcription from ATP and GTP, translation from charged tRNA) read the ODE pools.
- **The transcript likelihood is an indicator.** Transcripts are observed exactly, so
  the path must reproduce every observed count. This is D10's exact conditional, and
  it is not a Poisson modality.

### The three blocks (D10, amended by D16.1)

| Block | Conditional | Cost per update |
|---|---|---|
| B1 | `θ_ODE, σ_b \| X, y`: the prior, plus Σ_c log p(X_c \| θ), plus the bulk log-likelihood of `x̄`. A new θ_ODE needs every cell replayed. σ_b alone needs no replay | C replays per θ_ODE evaluation |
| B2 | `θ_CME \| X, θ_ODE`: each promoter and `krnadeg` enters only its own propensities. With the pool trajectory fixed, the likelihood is `k^n exp(−k·A)`, where `n` counts events and `A` is exposure summed over cells. The prior is lognormal, so this is not gamma. It is a one-dimensional log-concave density | Negligible: sufficient statistics only |
| B3 | `X_c \| θ, X_−c, y`, per cell, **in turn** (amended 2026-10-05): the bulk likelihood couples the cells, so each conditions on the others' current histories | The path update of D16.3, with cells sequential |

### Assumptions, each with its regime

| Assumption | Holds when | What breaks outside it |
|---|---|---|
| A path and θ determine the observed trajectory exactly | Fractional carry and the clamped drain are deterministic, and initial transcripts are fixed | Stochastic rounding would add a latent. It is not the data policy. Asserted by V2 |
| Rate constants change at 60m − 1 s on the jump clock | The handshake order above | Nothing to fix: the path update is written for it (D16.3), and V1 tests that boundary |
| A particle is its driver's full mutable state | The snapshot covers every field `handshake_step!` reads or writes (D16.3) | A missed field makes a restored particle run a different model. Asserted by V2b |
| The turnover ceiling never binds | `RNAPOL_KCAT·S_g < 180`, that is S below about 153 | B2's form becomes truncated. Asserted per proposal, and a violation throws |
| ~~Cells are independent given θ~~ | ~~Separate seeds and separate cells, as generated~~ | Amended 2026-10-05: true of the transcripts, not of the bulk metabolites, which couple the cells through `x̄`. B3 updates cells in turn |
| Block 1's step needs no gradient | Two free ODE constants plus σ | Freeing many more ODE parameters would bring the gradient, and 14c's sampler model, back. §7 |

### Verification checks, each with a tolerance

| # | Check | Asserted | Where |
|---|---|---|---|
| V1 | Path density | Three parts: the importance identity `E_{X~p(·\|θ0)}[p(X\|θ1)/p(X\|θ0)] = 1`, the hand checks, and the mutations below | 16a.2 |
| V2 | Replay fidelity | Replaying a generated cell's recorded path at its θ reproduces the carry-inclusive latent trajectory **bitwise** at every handshake and all 105 save points, on 3 dataset cells. If the replay changes the solver's step sequence, the fallback is ≤ 1e-8 relative, and the reason is recorded. **Mutation:** removing one translation event changes GTP | 16a.1 |
| V2b | Snapshot restore | Snapshot a replay at a window boundary, run another particle, restore, and continue. The result must equal the uninterrupted replay (bitwise, or 1e-8 with the reason recorded), at every window boundary of one cell. The intervening particle is constructed to differ in every snapshotted field, including a non-zero deficit, a different membrane ptsG count and a different SSA state, so no field passes vacuously. **Mutation:** leaving out `DeferredDebit.deficit` or `factor` fails it. **Amended 2026-09-30 (§12):** the SSA's cached next jump is redrawn at every handshake, so it needs no snapshot. A test asserts that corrupting it in a restored particle changes nothing, with the random stream seeded alike | 16a.1 |
| V3 | Exact transcript-only posterior | One gene, rate constants frozen (NTP pools chemostatted), transcripts observed every 60 s and nothing else. The exact posterior of `(S, krnadeg)` on a 2-D grid comes from the birth-death transition matrix. B2 with B3's bridge must match it: the 5, 25, 50, 75 and 95% marginal quantiles each within 3 Monte Carlo SE (ESS-based). **It cannot see a weight error that varies between particles, since the rates are frozen. V5 covers that** | 16a.4 |
| V4 | Block 1 on a fixed path | On M0, 3 cells, path fixed: B1's samples of `(ln kcatF_ENO, ln kcatF_FBA)` match a 2-D grid posterior computed by direct replay under the same forward-only prior. The grid density is evaluated on the θ scale, or on the ln θ scale with the Jacobian, and never as the same code path the sampler uses. Quantiles must agree within 3 MC SE. `assert_haldane` holds at every accepted state | 16a.5 |
| V5 | Path update against brute force, with the rates varying | See below the table | 16a.7 |
| V6 | Reference convergence | At least 4 chains from overdispersed starts. Rank-normalised R̂ < 1.01, and bulk and tail ESS ≥ 1,000 per parameter. The Monte Carlo SE of every reported 5% and 95% quantile is ≤ 0.05 posterior SD. This is what "the run's length justified rather than chosen" means | 16a.9 |
| V7 | Reference calibration | Simulation-based calibration of the M0 sampler at production settings. Each replicate's truth draws **all** free quantities from their priors: M0's targets by `draw_truth(...; names = M0_TARGETS)`, and σ_b from `p(σ_b)` (σ before 2026-10-05), each with a recorded seed. Rank-ECDF within 95% simultaneous bands for every free parameter, at R replications (≥ 50, K1's floor; 100 if the cap allows). **Amended 2026-10-12 (§12):** each band is at 1 − 0.05/6, so V7 as a whole is at 95%, on a complete set of R replicates | 16a.9 |
| V8 | F10's own floor | The divergence between the production and reference posteriors is reported beside the reference-vs-reference divergence of its split chain halves. *Amended 2026-10-22 (§12):* the divergence is Wasserstein-1 per marginal in reference SDs. A floor matched to production's length is added, and the reading is against that floor: "not resolved" means not above every matched null. A production divergence inside that floor is "not resolved", not "zero" | 16a.10 |
| V9 | Snapshot independence | Mutating one particle's snapshot leaves every other unchanged. This is the installed trap of solver objects shared across particles, tested directly | 16a.7 |

**V1 in detail.**
- The importance identity must hold within 3 SE, over 10⁴ M0 paths at two θ pairs:
  one CME-only shift, and one ODE-only shift that moves the pools.
- Two hand checks, each to 1e-12:
  - a one-gene path with frozen rates;
  - a path whose rates differ by at least 2× across a rebuild, with events placed on
    both sides of 60m − 1 s.
- Three mutations, each caught by the second hand check:
  - dropping the exposure term;
  - reading the rates from the previous interval;
  - putting the rebuild at 60m instead of 60m − 1.

**V5 in detail.** The kernel runs at fixed θ, since the reference is for block 3
alone. Two cases.
- **The M0 case: one gene, 3 windows.** It is small enough to enumerate by
  rejection. Acceptance is about 5 to 10%, so about 10⁵ runs of 180 s.
  - *Metabolite likelihood off.* Forward-simulate the full hybrid, rates live, and
    accept runs whose transcript counts match the data. The kernel's per-window
    event-count histograms and translation totals over 10⁴ sweeps must match the
    accepted runs' (χ², p > 0.01, multiplicity stated).
  - *Metabolite likelihood on.* Weight the accepted runs by it. The kernel's
    expectations of the same functionals must match within 3 SE. The reference's
    weighted ESS is reported and must be at least 1,000. If it is not, σ is inflated
    for this check until it is, and the value is stated.
  - **Amended 2026-10-05 (§12): under bulk metabolites, the "on" half runs on a
    population of C = 3 cells.**
    - The rejection reference simulates the three jointly. It keeps runs whose
      transcripts match all three cells' data, and weights them by the bulk
      likelihood of their mean.
    - The kernel is a sweep that updates the three cells in turn, each conditioned
      on the other two.
    - Expectations must match within 3 SE, with σ_b inflated only if the ESS needs
      it, and the value stated.
    - The toy case runs the same way.
- **An amplified case for the weight.** In M0 the transition factor varies by only
  about 0.07% between particles, too little for the χ² to see. So the same
  comparison also runs on a test composition built so the rates really differ
  between particles. Its transcription constant reads a pool that the translated
  protein drains strongly. It must meet two conditions, one for each factor of the
  weight:
  - two particles' constants differ by at least 2× within a window, which tests the
    60 s transition factor;
  - each particle's constant jumps by at least 2× across every rebuild, with
    transcription rates of order 1 per second, so the last second often carries an
    event. This tests the 1 s correction, which depends on a particle's own jump
    across 60m − 1 and not on differences between particles.

  **Mutations:** dropping the transition-probability factor, or dropping the 1 s
  correction alone, each fails this case. If a composition meeting the second
  condition cannot be built, V5 says so, and the weight unit test below is the only
  check of the 1 s correction.
- **A unit test of the weight function.** It is checked against the enumerated
  formula, with rate sets deliberately 2× apart, to 1e-12.

## 4. Approach

**D16.1. Block 1 is gradient-free on the published model, with a prior on the
forward constants only. This amends D10's "gradient-based sampler".**
**Why:** Block 1 has three dimensions: `ln kcatF_ENO`, `ln kcatF_FBA` and `ln σ`. A
slice or adaptive Metropolis step at three dimensions costs about as many replays as
a gradient step would. It samples the data-generating model exactly: clamped drain,
fractional carry, with no smoothing and no continuous pools. That removes 14c.5's
3.75% GTP departure from phase 16. It also sidesteps K5 ("gradient-based sampling is
invalid as posed") rather than working around it.
- **The reverse constants are derived at every proposal**, through the nominal
  equilibrium constant, exactly as the truth was drawn (`draw_parameters(...;
  reverse = :derived)`).
- **They carry no prior term and no Jacobian.** The prior density is
  `log p(kcatF_ENO) + log p(kcatF_FBA)`, and the reverse constants are deterministic
  functions of those. `d11_models()` frees `kcatR_R_ENO` and `kcatR_R_FBA` with priors
  of their own. **Do not reuse `inference.jl`'s one-prior-per-free-parameter vector
  or `model_free_params`.** Adding log p(kcatR) targets a different prior from the
  one the truth was drawn from, and V7 and coverage would fail for a reason that has
  nothing to do with the method.
- **σ is updated by its own one-dimensional slice step** on the cached replay, with
  no new replay.
- **Every step works on u = ln θ, and the prior term there carries the
  Jacobian.** For a `LogNormal(μ, s)` prior, the density in u is `logpdf(LogNormal,
  θ) + ln θ`, which is `Normal(u; μ, s)`. This holds for both forward constants and
  for σ, and the same rule applies in B2. Dropping the `+ ln θ` shifts each prior by
  −s² in u: 0.16 prior SD for ENO, the control K4 is scored on. The "no Jacobian" of
  the reverse constants above is a different statement. They are not sampled at all.

**Alternatives considered:**
- NUTS or HMC on 14c's sampler model, as D10 had it. That needs automatic
  differentiation through 6,300 handshakes (`rate_constants` builds an
  `SVector{n, Float64}` and would need work). It also carries 14c.5's GTP gap into
  every posterior.
- Building both and choosing on M0. This was declined: it doubles 16a for a
  three-dimensional block.

14c's sampler model stays in `src` for when the ODE block grows (§7). Task 16a.6
(parent 16.3) still runs, as the record of why the gradient was dropped.

**D16.2. Block 2 is updated from sufficient statistics, with no replay.**
**Why:** Given the path and the replayed pools, gene `g`'s transcription events give
`S^n_g · exp(−S·A_g)`, where `A_g = Σ_c ∫ RNAPOL_KCAT/denom_g dt`. Decay gives the
same form in `krnadeg`, with `A = Σ_c Σ_g ∫ m_g/n_g dt`. The lognormal prior makes
each one a one-dimensional log-concave density. It is updated by a univariate slice
step on u = ln k, with the Jacobian D16.1 states. Without it, the step behaves as if
each gene had one fewer event. That is a valid kernel, not an independent draw; adaptive
rejection sampling would give exact draws if mixing ever demands them. **The
conjugate direction is the degenerate one** (D10, D11). With the polymerase constant
fixed, that costs nothing here.
**Alternatives considered:** a gamma prior, which would give an exact conjugate draw.
It was rejected because it changes the prior the truth was drawn from.

**D16.3. The path update is conditional SMC over 60 s windows, and 16a.7 chooses
its variant by measurement.**

The parent's three options, examined against the model as built:

| Option | Status |
|---|---|
| No augmentation | **Closed.** The 60 s aggregate drain fails granularity at +288% (13.4) |
| Exact forward filtering over a truncated state space | **Not available for the joint path.** The ODE state is continuous and a function of the whole path, so no finite-state filter carries it. It survives as a component: each gene's transcript path between two exactly observed counts is a birth-death bridge, sampled exactly by FFBS on a truncated count space (uniformisation) |
| Particle Gibbs | **The route.** Per cell, a conditional SMC over the 105 observation windows W_m = (60(m − 1), 60m] |

**The proposal in window m, per particle, in the handshake order.**

For m = 1, R_0 is read from the driver's state after `set_parameters!`, not
recomputed. `set_parameters!` rebuilds only modules whose parameters changed, so a
new θ_CME must re-derive transcription's R_0 from the initial pools.

1. Over the first 59 s, where the constants R_{m−1} are fixed, propose each gene's
   transcript segment from the **60 s bridge at R_{m−1}**, ending at the observed
   count y_m. Keep its path up to 60m − 1. Translation and translocation come from
   their conditional prior given the transcripts, simulated forward.
2. Run the handshakes up to step 60m's rebuild. That gives R_m, from this particle's
   own pools.
3. Over the last second [60m − 1, 60m], propose each gene's segment from the
   **exact 1 s bridge at R_m**, from the count at 60m − 1 to y_m. Then take the jump
   step.

**The incremental weight is not 1, and not the metabolite likelihood alone.** Per
gene, write P_R(b | a; τ) for the birth-death transition probability over τ seconds
at constants R. Then:

```
w_m  =  g_m(y_met)  ·  Π_g  P_{R_{m−1}}(y_m | y_{m−1}; 60) · P_{R_m}(y_m | m_{60m−1}; 1)
                             / P_{R_{m−1}}(y_m | m_{60m−1}; 1)
```

- **`g_m(y_met)`** is window m's metabolite likelihood.
- **Why the transition probabilities stay in.** Target over proposal leaves them,
  and **they vary between particles**, because R is rebuilt from each particle's own
  pools. A first draft dropped them, and V5 is the check that they are back.
- **Cost.** The probabilities come from the same uniformisation that samples the
  bridges.
- **The count cap.** Its floor is set per dataset, as max(observed count) + 20
  (`bridge_cap`). **Amended 2026-10-01 (§12):** the floor is not always enough, so
  the cap grows in steps of 5 until the truncation mass is below 1e-12 at the
  constants in use (`adequate_cap`). Prior-drawn truths can put GAPD's mean near 35.

**A particle is the driver's full mutable state.** That is:
- the ODE state vector and the integrator's adaptive state (step size and cache);
- the jump state and its rate-constant slots;
- every counter;
- every `DeferredDebit.deficit` (the clamped drain's carried shortfall);
- the rounding remainders;
- the ODE parameter slots written by the catalytic and geometry channels;
- the growth state (`area_nm2`, `radius_nm`, `volume_litres`, `factor`);
- `n_handshakes`, which sets the drain and rebuild schedule;
- both integrators' times `t`;
- the SSA's cached propensities and pre-drawn next jump. A deep copy carries them,
  and every handshake redraws them.

**The random stream is not part of the driver** (amended 2026-09-30, §12). The SSA
draws from the task's global generator, which no copy of the driver captures. A
particle that simulates forward therefore sets that generator from its own recorded
seed before each forward step (16a.7).

**How the proposal is simulated.** Transcription and decay are suppressed in the SSA
and imposed from the bridge as scheduled events. Each imposed event applies its full
affect, including the cost counters, and refreshes the SSA's cached translation
propensities. Translation and translocation run
in the SSA forward, from the particle's own restored state and RNG stream.

V2b tests that restoring this reproduces an uninterrupted run.

**What 16a.7 decides among,** at N ∈ {5, 10, 20, 50} particles:
- **PG without ancestor sampling.** O(N·T) window-replays per sweep, about N cycles
  per cell. It is prone to path degeneracy in early windows.
- **Full PGAS.** The reference path's future continuous state changes with the
  ancestor. So each ancestor weight needs the reference's future replayed from that
  candidate, about N·T/2 cycles per cell per sweep. It is exact.
- **PGAS truncated at lag L.** About N·(1 + L) window-replays per window, which is
  N·(1 + L) cycles per cell. **It is not exact:** the dropped terms bias the kernel
  at any chain length (D16.5).

**The criterion:**
- Every window's segment changes in at least 10% of sweeps. That is the update rate,
  and the threshold is ⚠️ DRAFT.
- Cost per sweep, measured on M0 and on one full-scale cell.
- The choice is recorded against 13.4, D13 and the resolver's gradient report (14c.4:
  no obstruction), as the parent requires.

**The two installed traps (parent D10), met by construction:**
- The CSMC is hand-written over explicit snapshots. It does not go through Turing's
  `@model` or Libtask.
- Ancestor sampling does not depend on `SSMProblems`.
- One integrator is reused, and a particle is restored into it by copying state, so
  no solver object is shared across particles. V9 asserts this.

**Amended 2026-10-02 (§12): each window's proposal is annealed.** 16a.7's
measurement fired KP. The metabolite likelihood puts proposals from the
translation prior hundreds of nats below the reference within one window, so no
variant at N = 50 updates every window. The proposal for window m, per particle,
is now annealed importance sampling (AIS) inside the window:
- **The tempered targets.** π_β ∝ q · w^β over the window's path, where q is the
  proposal above and w its incremental weight. β runs over a fixed schedule
  0 = β_0 < β_1 < … < β_K = 1. ⚠️ DRAFT: geometric from β_1 = 1e-4.
- **The move at stage k** is Metropolis–Hastings at β_k. It draws s uniformly in
  the window's first 59 s, keeps the path before s, and redraws the rest from q's
  own conditional: each gene's bridge from its count at s to y_m at R_{m−1}, the
  other reactions forward from the state at s, then the last second as above. q
  cancels, so the acceptance is min(1, (w′/w)^β). The move is π_β-reversible.
- **The particle's weight** is the AIS weight, Π_k w(x_{k−1})^{β_k − β_{k−1}}, an
  unbiased estimate of ∫ q·w. It replaces w as the incremental weight. Its path is
  the chain's final state.
- **The reference particle** gets its auxiliary chain conditionally. From its path
  x_K, it draws x_{k−1} by the stage-k move from x_k, for k = K down to 1. Each move
  is its own reversal, so this is the extended target's conditional. Its weight is
  computed the same way, and its path stays the reference's.
- **Particle Gibbs only** (lag 0). Ancestor sampling would have to score the
  reference's future on the extended space, and 16a.7's measurement shows it did
  not help.
- **Cost.** About K + 1 window replays per particle per window, so about N(K + 1)
  cycles per sweep. At N = 5 and K = 50 that is about 255 cycles per cell per
  sweep, which projects to about 0.9M CPU-h per full-scale posterior on D16.6's
  assumptions. **The annealed sampler is affordable on M0, not at full scale
  under the cap.**

**Chosen 2026-10-05 (16a.7c): PG, lag 0, N = 20, under bulk metabolites** (corrected
after the #81 review from N = 10).
- **Per cell, over the windows whose true path has events,** the lowest update rate
  is 0.50 at N = 20 and 0.12 at N = 10, from 16 sweeps per cell (SE about 0.08). So
  N = 10 does not robustly clear 10%.
- **Windows with an empty true path** (transcripts 0 → 0) have posteriors that sit
  on the empty path. Their low rate is expected and is not counted.
- **Cost:** about 16.5 cycles per cell-sweep.
- **The measurement is on M0 only:** 5 of 50 cells updated, the rest held at their
  true paths, from true starts. A full-scale cell under bulk data is measured in
  16a.11. Against 13.4, D13 and 14c.4:
the path stays exact, at the published drain granularity, with no gradient
needed.

**Amended 2026-10-05 (§12): under bulk metabolites, the window weight's metabolite
term is the bulk likelihood with this cell's particle substituted.**
- **The weight.** For cell c, `g_m` becomes the bulk likelihood at
  `x̄_k = (S_−c,k + x_c,k)/C`. `S_−c,k` is the other cells' summed latent at save k,
  from their current histories, and `x_c,k` is the particle's own.
- **Order.** Cells are updated in turn within a sweep, each after the previous
  cell's update, which is a valid Gibbs scan. It costs wall-clock, not CPU-h. After
  a cell's update, `S` is refreshed from its new history.
- **Expected effect.** One cell moves `x̄` by about 1/C of its deviation, so the
  weights are gentle, and PG is expected to mix as the metabolite-off runs did
  (97 to 100% of windows changed). 16a.7c checks this, not assumes it.
- **The annealed window stays in the code,** V5-validated and off by default. It is
  the fallback if 16a.7c fails.

**Alternatives considered:**
- Single-site Metropolis-within-Gibbs on one gene's window segment. Every proposal
  needs a replay from that window to the end, about half a cycle.
- A whole-path independence proposal. Its acceptance collapses under 105 metabolite
  panels.
- AdvancedPS on `SSMProblems`. It is possible, but the driver's mutable state needs
  the same snapshotting anyway.

**D16.4. M0 is Core A′ cut to two genes over 600 s.** ⚠️ DRAFT: the genes are
ptsG (`JCVISYN3A_0779`, the headline, the only membrane protein) and GAPD
(`JCVISYN3A_0607`, the best-determined promoter).
- The other 15 proteins are held at their proteomics counts as fixed enzymes, which is
  a labelled M0 reduction in `reduction_report`.
- **`M0_TARGETS`:** `S_0607`, `S_0779`, `krnadeg`, `kcatF_R_ENO` and `kcatF_R_FBA`.
  σ is free too, but it is drawn separately and kept out of this list:
  `draw_truth` would throw on it. That is D11's set less `S_PGI`, so every block runs as it does in
  production.
- **Cells per M0 dataset:** ⚠️ DRAFT, 50. The number is set by the cap, not by D11's
  precision argument, and it is stated beside every M0 result.

**Why:** It is small enough that a very long run is defensible, per D10. It still
couples the two blocks in both directions and runs all three blocks. Amended
2026-10-05: information is shown to cross only one way for ENO at 15.7's truth.
**Alternatives considered:** one gene. That cannot exercise B2 across genes. V3 and
V5 already use one-gene cases.

**D16.5. "Production" and "reference" are the same sampler at different settings,
and what F10 measures depends on the variant 16a.7 picks.**
- The reference runs **full PGAS** (or PG, if that alone passes V6) at N_ref ≫ N, with
  chains far longer than V6 requires.
- **If production is PG or full PGAS,** it is exact in the limit. The particle count
  costs mixing, not bias. F10 then measures only the Monte Carlo and mixing error of
  production's settings.
- **If production is truncated PGAS,** it is a structural approximation whose bias
  does not vanish with chain length. F10 then measures that bias too.
- F10's caption states which case holds.
- *Amended 2026-10-18 (§12 2026-10-13):* the reference runs at N_ref = N = 20, as
  V6's 4 chains of 1,650 sweeps. Production is PG, so the particle count costs mixing,
  not bias, and V6 checks the reference's convergence directly. F10 compares a
  production-length chain with that reference on the same dataset.

**D16.6. A 50k CPU-h cap, labelled, and reopened at 16a.11** (parent §8's opening,
amended 2026-09-30).

**The projection per variant.** Take 2,000 sweeps and 200 cells at 32.2 s per cycle,
with E block-1 evaluations per sweep. E is between 5 and 15.

| Variant | Cycles per cell per sweep | CPU-h per full-scale posterior |
|---|---|---|
| PG, N = 10 | 10 + E | about 54k to 89k |
| Truncated PGAS, N = 10, L = 5 | 60 + E | about 230k to 270k |
| Full PGAS, N = 10 | about 525 + E | about 1.9M |

**So the cap cannot fund one full-scale posterior at any variant.**

*Amended 2026-10-29 (§12, O1 part 2):* the cap is reopened for track B's B8 with no
ceiling, by the user's decision. B8's spend is logged and labelled; the exact
full-scale posterior (O5) is not funded by it.

**Amended 2026-10-05 (§12, after the #81 review): wall-clock under bulk
metabolites.**
- **Why it matters:** B3 updates cells in turn (D16.3). The CPU-h above are
  unchanged, but a sweep's B3 can no longer run its cells in parallel.
- **The measured cost:** 16.5 cycles per cell-sweep at the chosen N = 20 on M0
  (16a.7c).
- **Full scale:** 200 cells × 16.5 cycles × 32.2 s ≈ 30 h per sweep on one core,
  or about 1.5 h with the 20 particles of each cell-sweep in parallel. Over 2,000
  sweeps that is about 60,000 h serially, or about 3,000 h, some four months, in
  parallel.
- **M0:** 50 cells × 16.5 × about 1.5 s ≈ 21 min per sweep, so about 29 days
  serially, or about 1.5 days in parallel.
- **The CPU-h projection** rises with N: PG at N = 20 is 20 + E cycles per cell per
  sweep, about 1.7 times the N = 10 row above.
- **So M0 is feasible in wall-clock, and a full-scale posterior needs
  particle-parallel sweeps at least** (R16.8). 16a.11 measures the real figure.
- **16a is capped at 5k CPU-h** (annotated 2026-10-06, §12: 16a.9 is recorded
  against it but not held to it). Its last task, 16a.11, does four things:
  - measures the cost of one M0 posterior, and of one sweep on a full-scale cell,
    per variant;
  - sets K1's bound, with T3 recording that it was set after the measurement;
  - proposes a costed 16b allocation;
  - **reopens the cap with you.**
- **16b does not start without your approval of the allocation.**
- **The priority order for the allocation.** ⚠️ DRAFT, and yours to reorder:
  1. 16b.1, the tight-prior control. This is K4's stop gate.
  2. 16b.2, the six-parameter recovery.
  3. F6's load-bearing cell, ptsG under metabolites only, and its measured null,
     ENO under transcripts only (amended 2026-10-05, §12).
  4. 16b.3, coverage.
  5. 16b.5, K2 at 1,000 cells.
  6. 16b.6, K7's two fits.
  7. F6's remaining cells.
- **What isn't funded at full scale drops down K1's ladder:** shorter horizon, then
  five genes, then M0. If nothing affordable remains, T3 scores the item **not run,
  budget-bound**. Every verdict the cap forces is labelled budget-bound in T3, never
  as the model's feasibility. Each drop is a §12 entry.

**Why:** 16b as written needs at least eight posteriors plus coverage. Projections
have been wrong in this spec before (D10's 5 s drain; K1's 10 s proxy), so the
allocation rests on the measurement.

**D16.7. The data model, as the dataset was made.**
- Transcripts are scored by the path indicator, not as Poisson.
- Metabolites are the carry-inclusive particle counts, with `floor = 1.0`.
- t = 0 is excluded.
- ⚠️ DRAFT: **`M_lac__L_e` is kept.** Its scale rests on the asserted
  medium-to-cell volume ratio. That ratio is fixed at the same value in the generator
  and the sampler, so synthetic recovery is unaffected. It is flagged in T2 as an
  observable that would not transfer to Step 2 unchanged.

**D16.8. K2 gets its 1,000 cells.** The 15.7 truth is extended by 800 cells at seeds
30201 to 31000, with the same model, truth and noise rule. The combined set is used
only for 16b.5's transcript-only fit. It stays subject to D16.6's allocation.

**D16.9. K7 gets a task (16b.6).**
- `kcatF_R_PYK` is freed, with its reverse constant derived and a forward-only prior.
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
  - coverage datasets (16b.3).

  Every replicate truth draws all the free quantities of its rung from their priors:
  `draw_truth(...; names = <the rung's targets>, purpose = :coverage)`, plus σ_b from
  `p(σ_b)` (σ before 2026-10-05), each seed recorded. Each dataset follows 15.7's convention: truth, seeds,
  noise seed, T2 row, `reduction_report` and Julia version in a tracked `meta.md`.
  Observations are regenerable and ignored, with a SHA-256.
- **Chains** are saved with their sampler seeds, settings and commit. Every quoted
  number cites a Slurm job id.
- **The ledger.** `dev/scripts/phase16_ledger.md` records CPU-h per job against the
  50k cap and 16a's 5k, and it is updated with every run of record.
- **Runs of record go in a worktree at a pushed commit.** The branch is never switched
  under a running job.
- **The environment is the parent's §5.** No new direct dependency is expected
  (D16.3). If one is added, it is pinned and it is in the depot before any Slurm job,
  since compute nodes have no network.

## 6. Outputs

The Claim column uses the parent's sub-claims C1 to C5.

| Output | From | Claim |
|---|---|---|
| The path-update decision record, with the measured update rates and costs per variant | 16a.7 | Makes D10's open choice |
| V3's exact-versus-sampled figure, and V5's brute-force comparison | 16a.4, 16a.7 | The conditionals are right |
| M0's rank-ECDF (F8 on M0) | 16a.9 | C4: the reference is calibrated |
| **F10** | 16a.10 | C4: makes every calibration claim attributable |
| F14, K1's bound, and the proposed allocation | 16a.11 | K1 |
| T4, and the 16b.1 perturbation plot | 16b.1, 16b.2 | C3, C4 |
| F9 coverage curve, and F8 where affordable | 16b.3 | C4 |
| **F6** | 16b.4 | C3: **the figure the work is for** |
| K2's verdict sentence | 16b.5 | C3 |
| K7's width ratio | 16b.6 | Scope honesty |
| The caveats block, and **T3** | 16b.7, 16b.8 | Scope honesty, falsifiability |

Every output carries its rung (M0, reduced, or full scale), its cell count, and a
budget-bound label where one applies.

## 7. Non-goals

- **Phase 17's conditional checks** (D12 parts 1 to 3, and K3). Phase 16 runs the
  joint sampler and coverage only.
- **A gradient step in block 1**, and any further validation of 14c's sampler model.
  That includes the end-to-end test of 14c.5's GTP mismatch, which parent §12
  retires (2026-09-30).
- **Freeing parameters beyond D11's six**, except `kcatF_R_PYK` for K7 alone.
- **The polymerase ridge as a sampling target.** F13 is 15.8's.
- **Any full-scale result the allocation did not fund.** It is scored as not run, and
  never extrapolated.

## 8. Kill criteria

How the parent's criteria are placed here:

| Criterion | Where it is scored |
|---|---|
| K1 | Its bound is set in 16a.11; its verdict is written in 16b.8 |
| K2 | 16b.5 |
| K4 | 16b.1 and 16b.3 |
| K5 | 16b.8. It **fired** in 14b. 14c answered it, and D16.1 makes phase 16's sampler not depend on it. T3 records all three |
| K6 | 15.8: it did not fire |
| K7 | 16b.6 |
| K3 | Phase 17's |

**Two local stops:**
- **KR — the reference is wrong.** V1, V2, V2b, V3, V5, V7 or V9 fails. Stop: nothing
  downstream is attributable, which is the point of R4.
- **KP — the path update cannot mix.** At the largest N that fits 16a's cap, some
  window's update rate stays below the D16.3 threshold for every variant. Stop and
  amend. The fallbacks are a shorter window, a guided translation proposal, or K1's
  ladder.
  - **Fired 2026-10-02** on PG, PGAS and truncated PGAS at N = 50
    (`dev/scripts/csmc_variants_16a7_result.md`). The amendment takes the guided
    proposal, as D16.3's annealed window proposal (§12 2026-10-02).
  - **If the annealed proposal also misses the threshold on M0 within the cap,** KP
    fires again. The remaining fallback is then to condition the path on transcripts
    alone, with the metabolites used in block 1 only. That is an approximation,
    labelled in T3.
  - **Fired again 2026-10-02** on the annealed proposal at K ≤ 100.
  - **Resolved 2026-10-05 by changing the observation model, not the sampler**
    (§12). The per-cell panel was unmeasurable, and the fallback above, a cut, would
    remove F6's ptsG cell. Under bulk metabolites, KP is re-scored by 16a.7c.

## 9. Open questions

- [NEEDS CLARIFICATION: is a replay cheaper than a simulation?] 13.7's 5.11 ms per
  handshake includes the SSA. Replay drops it but keeps the stiff solve. 16a.1
  measures this, and it moves every row of D16.6's projection.
- ~~[NEEDS CLARIFICATION: how degenerate are the particle weights?] GTP's 200-cell
  resolution of 4.5% (14c.5) implies a per-cell coefficient of variation near 64%, so
  metabolite panels discriminate paths strongly. Whether 50 particles suffice is
  16a.7's measurement.~~ — resolved 2026-10-02: completely degenerate. Fifty
  particles do not suffice. `M_pi_c`, `M_3pg_c` and `M_2pg_c` near depletion put
  proposals hundreds of nats below the reference (§12 2026-10-02).
- ~~[NEEDS CLARIFICATION: is the window proposal right at full scale?] On cell 30001,
  7 of 17 pools' ranks sit at KS p = 0.02 to 0.03, over windows of one trajectory.
  M0's ranks are uniform. Task 16a.7b decides it on independent cells.~~ —
  resolved 2026-10-02: it is. 16a.7b passes on five cells.
- [NEEDS CLARIFICATION: which of PYK's upstream reactions are the "seven"?]
  Enumerated in 16b.6.
- ~~[NEEDS CLARIFICATION: is the initial transcript count fixed or drawn in 15.7's
  generator?]~~ — resolved: fixed, at `round(Int, mean_mrna)` (§2).
- The ⚠️ DRAFT items above are still open: the update-rate threshold, M0's genes and
  cell count, D16.6's priority order, `M_lac__L_e`, and K7's sink form.

## 10. Risks

| # | Risk | Early warning | Response |
|---|---|---|---|
| R16.1 | Cost exceeds the cap at every useful rung | D16.6 already projects this. 16a.11 measures it | The cap is reopened at 16a.11. Anything else is labelled budget-bound |
| R16.2 | Path degeneracy in early windows | Update rate in windows 1 to 10 falling as T grows, at fixed N | PGAS or truncated PGAS, with its cost row; KP |
| R16.3 | The replay or a restore is not bitwise | V2 and V2b | Record the reason, and use the 1e-8 fallback. A non-deterministic replay would mean a hidden latent, which is a stop |
| R16.4 | A reverse constant gets a prior term, or is sampled freely | V4's `assert_haldane`, and V7 failing on ENO and FBA alone | D16.1's forward-only prior and derived reverse constants |
| R16.5 | F10 reads zero because the reference is too short | V8's floor | Lengthen the reference. Never report a divergence below its floor as agreement |
| R16.6 | 16a overruns its 5k cap on pilots | The ledger | Stop, and report at 16a.11 with what was measured |
| R16.7 | The annealed proposal mixes, but only at a K no rung can fund | 16a.7a's cost per sweep | M0 only; full-scale rungs drop down K1's ladder or are scored budget-bound |
| R16.8 | Updating cells in turn makes a full-scale posterior take months of wall-clock (D16.6, amended 2026-10-05) | 16a.11's measured sweep time | Run each cell-sweep's particles in parallel. A parallel scan over cells would need its own exactness argument; otherwise full-scale rungs drop down K1's ladder |

## 11. Task list

### Phase 16a — The sampler and the reference

**Goal:** build the three-block sampler and validate it exactly where that is possible
and against a long reference where it is not. Then measure what one posterior costs.
**Done when:**
- V1 to V9, and V2b, pass.
- The path-update variant is chosen and recorded with its measurements.
- The M0 reference has converged (V6) and is calibrated (V7).
- F10 is reported against its floor.
- K1's bound is set, and the reopened budget and 16b allocation are approved.

**PR:** #81 (merged 2026-10-05) carries 16a.1 to 16a.8, with 16a.7a to 16a.7c and
parent 15.7b. 16a.9a is its own PR; 16a.9b to 16a.11 follow.

- [x] 16a.1 Drive the published hybrid from a recorded jump path, and snapshot and
  restore the driver's full state (D16.3's list). Record the path during generation.
  — verify by V2 on 3 dataset cells (bitwise at every handshake), by V2b and its
  mutation, and by the replay's wall-clock per cycle reported against 13.7's 32.2 s.
  **Met** (path_replay_16a1_result.md: V2 bitwise on three cells at every handshake, V2b at all 104 boundaries; jobs 17755845 and 17755846).
- [x] 16a.2 Implement `log p(X | θ)` from the replayed rate constants, on the jump
  clock's [60m − 1, 60m + 59) intervals — verify by V1, including both hand checks
  and all three mutations. **Annotated 2026-10-01 (§12):** the importance identity
  ran at suite size on the phase 4 toy. Its 10⁴-path run on M0 follows 16a.8.
  **Met** (path_density_16a2_result.md: hand checks to 1e-12, the mutations caught, and the 10⁴-path identity on M0 within 1.22 and 1.07 SE).
- [x] 16a.3 Implement the per-gene transcript bridges (60 s and 1 s, by FFBS with
  uniformisation) and their transition probabilities, with the per-dataset cap —
  verify by:
  - sampled bridges always hitting both endpoints;
  - their event-count distribution matching the exact marginal (χ², p > 0.01, 10⁵
    draws per case, multiplicity stated);
  - the transition probabilities matching a matrix exponential to 1e-12;
  - the truncation mass asserted below 1e-12. **Amended 2026-10-01 (§12):** the
    +20 floor fails that at amplified rates, and `adequate_cap` grows it.
  **Met** (bridge_block2_16a34_result.md: exact against matrix exponentials and an augmented chain, 27 of 27; job 17758553).
- [x] 16a.4 Implement block 2 from sufficient statistics, with the ceiling asserted per
  proposal — verify by V3.
  **Met** (bridge_block2_16a34_result.md: V3 within 3 MC SE, 16 of 16; job 17758553).
- [x] 16a.5 (worked after 16a.8, since V4 runs on M0; §12 2026-10-01) Implement block 1: gradient-free, forward-only prior, reverse constants
  derived, σ updated separately — verify by V4, and by a test that the log-prior
  contains no `kcatR` term.
  **Met** (block1_16a5_result.md: V4 passes (job 17839209), and again on the bulk form (csmc_bulk_16a7c_result.md); FBA reading high is a watch item for V6 and V7).
- [x] 16a.6 (parent 16.3) Confirm a clipped drain under a fixed path still has a
  discontinuous derivative in an ODE parameter — verify by a test, as the parent's
  annotation specifies:
  - on the clipped model, a finite-difference slope in `ln kcat` that jumps across a
    clip and does not shrink with the spacing;
  - paired with 14c.5's continuity test on the smoothed model;
  - check 7's census cited as the source of the clip, not assumed.

  Recorded as the reason for D16.1.
  **Met 2026-10-13 on the toy** (`fixed_path_kink_16a6_result.md`; job 17801857 at
  `b7da3f0`, `test/test_fixed_path_kink.jl`, 9 of 9).
  - The clamp's slope jump is 1.265 at 30 steps and 1.267 at 3. It does not shrink.
  - The smoothed drain's jump falls from 0.236 to 0.024.
  - Check 7's census places the clip at handshake 3.

  **The Core A′ scan stays open, by decision of 2026-10-01.** It is not part of
  this clause. ENO and FBA cannot move the first clip, and PGK3's zero slope is
  unexplained. Block 1 is gradient-free (D16.1), so nothing downstream rests on it.
- [x] 16a.7 (parent 16.2) Build the conditional SMC path update, with D16.3's weight,
  and choose PG, PGAS or truncated PGAS. Verify by:
  - V5 in both cases, its mutations and the weight unit test, and V9;
  - per-particle streams: two restores of one snapshot, seeded alike, give identical
    forward proposals, and seeded differently they differ (§12 2026-09-30);
  - per-window update rates and cost per sweep, on M0 and on one full-scale cell, for
    each variant at N ∈ {5, 10, 20, 50};
  - the choice recorded in §4 against 13.4, D13 and 14c.4's gradient report, with the
    two installed traps addressed as D16.3 states.

  **Measured 2026-10-02; KP fired (§12 2026-10-02).** The rates and costs are in
  `dev/scripts/csmc_variants_16a7_result.md`, at N = 50, the largest. The smaller
  N were not run, since N = 50 decides KP. The choice moves to 16a.7a, and the
  full-scale check to 16a.7b.
  **Annotated 2026-10-05 (§12):** the variant is PG (lag 0) under bulk metabolites,
  at an N 16a.7c sets.
  **Closed 2026-10-13: met in part, and superseded in part.**
  - **Met:** V5 in both cases, the weight unit test, V9 and the per-particle streams
    (`test/test_csmc.jl`, in the suite of job 18461739).
  - **Superseded:** the rates at N < 50, which were not needed once N = 50 fired KP.
    The choice is made in 16a.7c (PG, lag 0, N = 20). The rates and cost on a
    full-scale cell under bulk data are measured in 16a.11, as 16a.7c says.
    The clause's "its mutations" has no definition: §3's V5 lists none (the three
    mutations there are V1's), so nothing further is ticked for it.
- [x] 16a.7a (§12 2026-10-02) Build D16.3's annealed window proposal, under particle
  Gibbs, and measure it on M0. Verify by:
  - a unit test that the reference's backward chain uses the stage moves in reverse
    and returns its own path;
  - a unit test that the AIS weight is unbiased: on the toy, its mean matches plain
    importance sampling's within 3 SE;
  - V5 in both cases with the annealed kernel, at a short schedule (K = 5), since the
    construction's correctness does not depend on K;
  - per-window update rates and cost per sweep on M0's five cells, at N ∈ {2, 5} and
    K ∈ {20, 50, 100};
  - the choice, or KP, recorded in §4.

  **Measured 2026-10-02; KP fires again, and it is open for your decision.**
  - V5 passes on both cases, and the tests pass.
  - No setting reaches the threshold: at K ≤ 100, no window changes in more than 5%
    of sweeps.
  - Inside one window, K ≥ 200 makes windows 2 and 3 exchangeable, and window 1
    resists K = 1,000. At about 42k CPU-h per M0 posterior from K = 200, that is
    above the cap.
  - `dev/scripts/csmc_variants_16a7_result.md`.

  **Superseded for production 2026-10-05 (§12)** by the bulk observation model. The
  code and its V5 pass stand, with the annealing off by default.
  **Closed 2026-10-13 as superseded.** KP's second firing was answered by the
  observation model, not by a choice here.
- [x] 16a.7b (§12 2026-10-02) Check the window proposal at full scale. Verify by the
  transcript-weighted rank of the true path's panel pools among proposals from its
  own state, on five independent 15.7-truth cells (30001 to 30005), with each pool's
  ranks uniform by KS at p > 0.01 after Bonferroni over the 17 pools.
  **Passed 2026-10-02** (jobs 17904315 and 17904316): over 525 windows, the smallest
  KS p is 0.057, and every Bonferroni p is ≥ 0.97.
- [x] 16a.7c (§12 2026-10-05) Implement bulk metabolites in the sampler, and choose
  N. Verify by:
  - M0 datasets with the bulk modality (parent 15.7b);
  - B1 scoring `x̄`, with V4 rerun on its bulk form;
  - B3 updating cells in turn with the substituted-mean weight, with a unit test
    that the weight equals the bulk likelihood of the population with that cell's
    particle in place;
  - V5 in both cases, with the bulk "on" half on C = 3 cells (§3);
  - per-window update rates and cost per sweep on M0 at N ∈ {5, 10, 20, 50};
  - **the choice of N, or KP,** recorded in §4. The variant is PG, lag 0.

  **Measured 2026-10-05; KP does not fire** (`dev/scripts/csmc_bulk_16a7c_result.md`).
  - **Per cell, over windows whose true path has events,** the lowest rate is 0.50
    at N = 20 and 0.12 at N = 10 (16 sweeps per cell). Corrected after the #81
    review, which found the first figure pooled over cells.
  - The cost is about 0.84N cycles per cell-sweep.
  - **The conditions:** M0 only; cells 1 to 5 of 50 updated, the rest held at their
    true paths; true starts.
  - **V5's bulk half passes on both cases** against the pooled reference
    (|z| ≤ 1.79). The fresh-only, fully independent checks give off χ² p ≥ 0.032
    and |z| ≤ 2.07. M0's bulk half ran at σ_b inflated 8× (0.8).

  **Met 2026-10-05.**
  - **V4 on the bulk form passes**, with all ten quantiles within 3 SE (jobs
    18061355 to 18061363). The largest deviation is FBA at its 25% point, +2.57 SE.
    The block 1 target check agrees to 1.7e-12.
  - **15.7b is done.**
  - **The window-2 watch item is resolved** by a fresh, independent reference.
  - **The first bulk V5 merge understated its reference SE** (shared runs). It now
    bootstraps over runs.
  - **Chosen: PG, lag 0, N = 20** (D16.3).
- [x] 16a.8 (worked before 16a.5; §12 2026-10-01) (parent 16.1, build) Build M0 and generate its datasets, with σ drawn per
  replicate. `generate_dataset` gains a `names` pass-through to `check_truth`, which
  defaults to `D11_TARGETS` and would refuse an M0 truth — verify by:
  - the build passing completeness mode;
  - `reduction_report` naming the 15 held enzymes;
  - each dataset carrying its truth (σ included), seeds and report;
  - two runs at one seed being identical.
  **Met** (m0_16a8_result.md: 21 of 21; job 17800819).
- [x] 16a.9 (split into 16a.9a and 16a.9b — see amendment 2026-10-06) (parent 16.1, run) Run the M0 reference long — verify by V6, which
  justifies the length, and by V7 at production settings on R ≥ 50 datasets, with R
  stated.
- [x] 16a.9a (§12 2026-10-06) Compose the three blocks into one Gibbs chain on M0,
  and pilot it. Verify by:
  - B2 being exact on the hybrid: its conditional's difference between two values
    of each promoter and of `krnadeg` equals the replay's path log-density
    difference plus the prior difference, to 1e-10;
  - a cell's base rebuilt by writing θ onto a pristine driver replaying bitwise
    as one written incrementally;
  - the threaded B3 sweep and B1 replay being bitwise equal to serial at one seed;
  - a chain resumed from a checkpoint being bitwise equal to an uninterrupted one;
  - a pilot of 4 chains from overdispersed starts on one M0 replicate (C = 50, PG,
    N = 20), reporting the cost per sweep per block, bulk and tail ESS per sweep
    and R̂ per parameter, per-window update rates, and a CPU-h and wall-clock
    projection for V6 and V7.
  **Met** (gibbs_pilot_16a9a_result.md: tests 28,658 of 28,658, job 18109288 at
  `b030512`; pilot jobs 18108335 to 18108339 at `6f8bcbc`. The slowest bulk ESS
  is 0.108 per sweep, on the ptsG promoter, and the lowest live update rate is
  0.34. V6 projects to 11.4k CPU-h and V7 at R = 50 to 57k).
- [x] 16a.9b (§12 2026-10-06; staged by §12 2026-10-07, sized by §12 2026-10-08) Run the M0
  reference long — verify by V6, which justifies the length, and by V7 at
  production settings on R ≥ 50 datasets, with R stated. Sized from 16a.9a's
  pilot.
  **Met.**
  - **V6** (v6_16a9b_result.md, merge 18420990): 4 chains at 1,650 sweeps. Worst
    R̂ 1.007, bulk ESS 1,081, quantile MCSE 0.049 SD, all on the GAPD promoter.
  - **V7** (v7_16a9b_result.md, merge 18689602 at `0ff6cff`): **R = 50**, under
    §12 2026-10-12's rule. All six parameters are within their bands at
    1 − 0.05/6, and also without the four replicates below ESS 30. F8 is in
    `dev/notes/figures/v7-16a9b/`.
  - **Caveat:** 23 of the 50 replicates have a smallest bulk ESS below L = 100,
    so the check is less sensitive than a well-mixed R = 50.
  - **Cost:** V7 used 33.9k CPU-h, and phase 16 is at about 42k of the 50k.
- [x] 16a.10 (parent 16.8) Compare production against the reference on M0 — verify
  by F10: overlaid marginals, one 2-D contour for the ptsG and ENO pair, and a
  per-marginal divergence beside V8's floor. The caption states which case of D16.5
  holds.
  **Done 2026-10-22** (`dev/scripts/f10_16a10_result.md`; job 18751635, 507.5
  CPU-h). D16.5's PG case holds.
  - **Result:** a production-length chain is not resolved from the reference on any
    marginal, under the rule of §12 2026-10-22.
  - **Biggest gap:** chain 5 ranks above all four matched nulls only on
    `kcatF_R_FBA`, by 0.020 SD. About 1.2 of six are expected to by chance.
  - **V8's split-half floor alone** would have flagged five of six, because it
    measures the noise of much longer runs.
  - **Figures:** `dev/notes/figures/f10-16a10/`.
- [x] 16a.11 Measure cost, set K1's bound, and reopen the budget. Verify by:
  **Done 2026-10-22** (§12 "16a.11", approved).
  - An M0 posterior costs 507.5 to 678 CPU-h.
  - A full-scale production posterior is projected at about 21k CPU-h.
  - K1's bound is 1,000 CPU-h per posterior. It fires at full scale (budget-bound)
    and not at M0.
  - Phase 16 closes at M0. 16b is funded at about 2.6k CPU-h.
  - the ledger showing 16a ≤ 5k CPU-h;
  - the cost per posterior on M0, and a full-scale projection per variant from
    measured sweeps;
  - K1's bound written into parent §8, with a note that it was set after the
    measurement;
  - a proposed allocation table funding each 16b task at a named rung, with any
    budget-bound item named;
  - a §12 entry;
  - **your approval of the budget and the allocation before 16b starts.**

  **Annotated 2026-10-13 (parent §12):** read
  `dev/notes/2026-10-13-inference-paths.md` before proposing the allocation. It
  weighs a faster integrator, modular inference and amortised SBI against this
  sampler, and none is decided.

### Phase 16b — Recovery on Core A′, at the rungs 16a.11 funds

**Goal:** the parent's recovery, coverage and shrinkage results, each at a stated
rung.
**Done when:** the parent's phase 16 condition holds at a funded rung:
- coverage within Monte Carlo error of nominal for every target;
- the tight-prior control shrinks and recovers;
- the loose-prior control is visibly wider.

In addition, K1, K2, K4, K5 and K7 are scored in T3. **If coverage is funded at no
rung, 16b ends incomplete.** T3 then scores coverage as budget-bound, and the
parent's phase 16 stays unticked. It does not count as done.
**PR:** #105 (merged 2026-10-28) carries 16b.1 and 16b.4's two runs; #107 (merged 2026-10-28) draws F6 at M0 and ticks 16b.4

**Funded at the M0 rung (§12 "16a.11", approved 2026-10-22).** Each task is
annotated below with its rung and funding.
- **The done-when** is read at M0. The parent's phase 16 is ticked as "met at M0;
  full scale budget-bound".
- **V7 serves as 16b.3,** with its caveats carried into T3.

- [x] 16b.1 (parent 16.4) Recover the tight-prior control alone — verify by:
  - the posterior concentrating on truth well inside its prior;
  - a deliberate perturbation of the truth moving the posterior with it;
  - shrinkage below 0.9.
  **M0, funded (~700 CPU-h).** Recovery is read off V6. The perturbation is one new
  M0 dataset with ENO's truth moved by +1 prior SD, and one production posterior on it.

  **If it does not recover, the machinery is wrong and the phase stops (K4).**
  **Done 2026-10-28 (§12 "16b.1 and 16b.4 as run"):** recovers, so K4 does not fire.
  Shrinkage 0.258, the moved truth at quantile 0.099, and a shift ratio of
  0.951 ± 0.021 (job 19034901; `dev/scripts/m0_16b_result.md`).
- [x] 16b.2 (parent 16.5) Recover the full six-parameter set — verify by:
  - every marginal containing truth;
  - at least one target shrinking below 0.5 (K4);
  - the cost per posterior reported against 16a.11's projection.
  **M0, reuse V6 (0 CPU-h).** Five of the six; `S_PGI` is not in M0.
  **Done 2026-10-28 at M0** (`dev/scripts/m0_16b_result.md`, "16b.2 and 16b.3"):
  - every truth lies inside its central 95% (quantiles 0.07 to 0.95);
  - all six shrink below 0.5 (0.060 to 0.331);
  - the 16b production posteriors cost 532.8 to 586.5 CPU-h, inside 16a.11's
    507.5 to 678.
- [x] 16b.3 (parent 16.6) Compute coverage over repeated datasets — verify by:
  - 50, 80, 90, 95 and 99% intervals covering truth at those rates within binomial
    error, per parameter, at ≥ 50 replications, with the count and rung stated;
  - **the tight-prior control at 100 replications**, inside K4's band of 90 to 99 out
    of 100 at nominal 95%. If 100 is not funded, K4's band is restated at the R that
    is run through a §12 amendment, never silently;
  - K4's third trigger (equal widths for ENO and FBA) checked.
  **M0, reuse V7 (0 CPU-h; call 1, approved).** Read off by
  `dev/scripts/v7_coverage_16a11.py`. K4's band is restated at R = 50 as 44 to 50 at
  nominal 95%.
  **Done 2026-10-28 at M0** (job 19399322; `dev/scripts/m0_16b_result.md`):
  - 29 of 30 cells are in band; the miss is krnadeg at 99%;
  - ENO covers 46 of 50 at 95%, inside 44 to 50;
  - K4's third trigger does not fire (ENO 0.267 against FBA 0.331 of prior SD);
  - call 1's four caveats are carried into T3.
- [x] 16b.4 (parent 16.7) Produce F6 — verify by all six targets under
  transcripts-only, metabolites-only and joint data, with ptsG under metabolites only
  (load-bearing) and ENO under transcripts only (the measured null) present. Any cell 16a.11
  did not fund is shown as not run. **Amended 2026-10-05 (§12):** ENO under
  transcripts only is reported against its information bound
  (`dev/scripts/eno_path_information_16a7_result.md`) as a measured null. A
  shrinkage there is a sampler error to investigate, not a result.
  **M0, funded (~1.4k CPU-h).** Transcripts-only and metabolites-only fits on V6's
  dataset; the joint fit is V6. First check that the Gibbs driver can drop a stream.
  **Held open 2026-10-28 (§12 "16b.1 and 16b.4 as run"):** both chains ran. Under
  metabolites only the transcript parameters did not mix (ESS 6–7), so ptsG's
  load-bearing cell is not measured. The user decides between reporting that arm
  as unmixed and redesigning its move, then rerunning it. F6 follows that choice.
  **Decided 2026-10-28 (§12 "16b.4 at M0: metabolites only reported as unmixed"):**
  the first option. F6 at M0 is drawn from the three fits as they stand, with
  S_0607 and ptsG under metabolites only shown as unmixed, not as a shrinkage, and
  krnadeg shown with its ESS of 39. The redesign and rerun is kept as an option to
  revisit. 16b.4 is ticked when F6 is drawn.
  **Done 2026-10-28 at M0:** F6 is drawn in `dev/scripts/m0_16b_result.md` (job
  19395555 at `f0dd5b2`), under §12 "F6 at M0: the divergence, its floor and how
  each cell reads".
  - Both named cells are present. ptsG under metabolites only is unmixed and not
    read.
  - ENO under transcripts only is at 0.0041 nat, below its floor of 0.0194, so it
    is not resolved from the prior.
- [ ] 16b.5 (parent 16.9) Generate the 1,000-cell extension (D16.8) and decide K2.
  Verify by:
  - the extension's `meta.md`;
  - the prior-to-posterior divergence under transcript-only data, for every metabolic
    target, against 0.05 nat;
  - the verdict written as a sentence, with its consequence for framing.
  **Budget-bound at 1,000 cells.** M0's transcripts-only divergence is read off 16b.4
  and labelled as M0's, not as K2's verdict.
  **M0's read-off, 2026-10-28 (job 19395555):** under transcripts only, ENO is at
  0.0041 nat and FBA at 0.0081, both below 0.05 nat and below their floors. K2's
  firing condition is met at M0, 50 cells. That is not K2's verdict.
  **Closed in T3 as budget-bound** (§12 "T3 at M0"). Left unticked: the
  1,000-cell extension was never run.
- [ ] 16b.6 Measure K7 (D16.9). Verify by:
  - PYK's dropped sibling reactions enumerated from upstream at `db048ac`, with the
    count checked against the parent's "seven";
  - the two fits' posterior width ratio on `kcatF_R_PYK`;
  - K7 scored against twofold.
  **Not run, budget-bound** (not funded at 16a.11).
  **Closed in T3 as budget-bound** (§12 "T3 at M0"). Left unticked.
- [x] 16b.7 (parent 16.10) Label what the result does and does not license — verify
  by the report carrying:
  - the scoping note's caveats: best-measured region, enzyme competition partly
    removed, architecture evidence rather than well-posedness;
  - K7's factor;
  - every result's rung and budget-bound label.
  **Done 2026-10-28 at M0:** `dev/scripts/phase16_m0_license.md` carries:
  - every result's rung and label;
  - the scoping note's three caveats;
  - K7 as not measured (budget-bound);
  - that neither direction of cross-block information is measured on the exact
    sampler at M0.
- [x] 16b.8 (parent 16.11) Fill T3 — verify by:
  - all seven criteria carrying a threshold, a measured value and a verdict,
    **published whether or not everything passed**;
  - K1's verdict written here, with its threshold marked as set after its
    measurement;
  - K5 recorded as fired in 14b and answered by 14c and D16.1;
  - every verdict the cap forced labelled budget-bound.
  **Done 2026-10-28 at M0:** `dev/scripts/phase16_t3.md` (§12 "T3 at M0").

## 12. Amendment log

### 2026-10-29 — O1 part 2 withdrawn: B8 is parked, and the cap is not reopened (the user's decision)

**Trigger:** the v0.3.0 launch decisions (`spec/spec.md` §12 2026-10-29). They
park B8 outside v0.3.0 and paper 1.

**Change:** the entry below, "O1 part 2", is superseded. The only B8 jobs run
under it were the six profiling jobs (19487970/71, 19489069/70, 19506190/91), and
no AV7 job was submitted. See `approximate-method.md` §12 2026-10-29, "B8 parked".

**Sections touched:** §12.

### 2026-10-29 — O1 part 2: the cap is reopened for track B's B8, with no ceiling (the user's decision)

**Trigger:** 16a.11 deferred O1 part 2 "until an approximate method passes 5.2". B7
passed AV6's SBC (`approximate-method.md` §12 2026-10-27). B8, AV7's ≥ 200
injections on Core A′, is projected roughly at 10k–25k CPU-h, against about 4.6k
left of D16.6's 50k.
**Change:** the user's decision (2026-10-29): B8 runs in full and is not held to
the cap. Its spend is logged in `phase16_ledger.md` as before, and labelled B8.
The gate on B8 is wall-clock (`approximate-method.md` §12 2026-10-29, "B8 as
planned"). O5, an exact Core A′ posterior, is not funded by this.
**Sections touched:** D16.6 (by reference), §12.

### 2026-10-28 — T3 at M0: K2 scored budget-bound; the forward direction unmeasured

**Trigger:** 16b.8. T3 needed two scoring choices that the results left open:
- K2 is defined at 1,000 cells, which is not funded, yet M0's read-off meets its
  firing condition at 50 cells (job 19395555);
- K2's response assumes the forward direction is measured, and on the exact
  sampler at M0 it is not (ptsG under metabolites only is unmixed).

**Change (the user's decisions, 2026-10-28, after the results):**
- **K2 is scored budget-bound,** with M0's read-off beside it: the firing condition
  is met at 50 cells, and that is not K2's verdict. K2's definition is unchanged.
- **The forward direction is recorded beside T3 as unmeasured** on the exact
  sampler at M0. It rests on track A's linearised estimate and the approximate
  route until measured. The redesigned move and rerun (§12 "16b.4 at M0") stays an
  option, to revisit when the paper's framing is set.

**Result:** T3 is in `dev/scripts/phase16_t3.md`.
- K1 does not fire at M0, and fires at full scale (budget-bound).
- K2 is budget-bound.
- K3 is phase 17's.
- K4 does not fire.
- K5 fired in 14b, and was answered by 14c and D16.1.
- K6 does not fire.
- K7 is budget-bound.

16b.8 is ticked. 16b.5 and 16b.6 stay unticked and are closed as budget-bound in
T3. The parent's phase 16 is met at M0, with full scale budget-bound (§12 "16a.11",
call 2).

*Sections:* header status, §11 16b.5, 16b.6 and 16b.8.

### 2026-10-28 — F6 at M0: the divergence, its floor and how each cell reads

**Status: recorded 2026-10-28, before any F6 divergence was computed.** The
shrinkages are the ones already in job 19378382.

**Trigger:** F6 reports "the prior-to-posterior divergence" without saying which
one. K2 names `kl_divergence` (`src/boundary.jl`), and 16b.5's M0 read-off reuses
whatever F6 uses.

**The divergence.** It is KL(posterior ‖ prior) in nats, per marginal, on the log
scale: `kl_divergence(post, prior; n_grid)`.
- `post` holds the kept draws: sweeps 101 to 656 of each 16b.4 chain, and 101 to
  1,650 of each of V6's four chains, pooled.
- `prior` holds 10⁵ draws of each parameter's prior on the log scale, from
  `Xoshiro(16_400_100)`. Every prior is log-normal, so this is
  Normal(μ, σ) from `params(prior)`.
- **The grid.** The default 200 points span the union of both supports, which
  would leave σ_b's posterior about one point per SD. `n_grid` is raised to 4,000.
  It is checked by 8,000: every cell must move by less than 2% or 0.002 nat.
- **A cross-check.** Each cell is also given the Gaussian KL, with the posterior
  fitted as a Normal. A cell where the two differ by more than 20% is reported
  with both, and the difference is named as the estimator's, not the data's.

**The floor (R16.5).** Two samples of one distribution give KL > 0. For each
cell, 20 independent samples of the prior, each the size of that cell's bulk
ESS (rounded), are scored against the reference prior sample. The 95th
percentile of those 20 is the floor. A divergence at or below its floor is
reported as "not resolved from the prior", never as zero.

**How each cell reads.**
- **Not read:** a cell with split-R̂ above 1.05. Its SD ratio and divergence are
  not reported. It is shown as unmixed, with its ESS and R̂. This covers S_0607 and
  ptsG under metabolites only (§12 "16b.4 at M0: metabolites only reported as
  unmixed").
- **Read, flagged:** a cell with bulk ESS below 100 and split-R̂ at most 1.05. It
  is reported with its ESS beside it.
- **Held:** σ_b under transcripts only is held at its start, so it has no
  posterior.
- **The joint column** is V6, as in job 19378382.

**16b.5.** ENO's and FBA's divergence under transcripts only, read against 0.05
nat, is M0's K2 read-off. It is labelled as M0's at 50 cells, not as K2's verdict.
16b.5 stays unticked.

*Sections:* §11 16b.4, 16b.5.

### 2026-10-28 — 16b.4 at M0: metabolites only reported as unmixed

**Trigger:** the choice that §12 "16b.1 and 16b.4 as run" left open. Under
metabolites only, S_0607 and ptsG did not mix (ESS 7 and 6, split-R̂ 1.13 and
1.16), so ptsG's load-bearing cell is not measured on the exact sampler.

**Change (the user's decision, 2026-10-28, after the results):** report the
metabolites-only transcript parameters as unmixed. No rerun is made now.
- F6 at M0 is drawn from the three fits as they stand: transcripts only (job
  19034902), metabolites only (job 19034903) and joint (V6).
- In the metabolites-only column, S_0607 and ptsG are shown as unmixed, with
  their ESS and split-R̂, and no shrinkage is read from them. krnadeg is shown
  with its ESS of 39.
- The metabolic parameters under metabolites only, and every cell under
  transcripts only, are reported as measured.

**Kept as an option to revisit:** redesign that arm's sampler move (for example
a joint S_g-and-path update), then rerun it, about 550 CPU-h at today's cost per
sweep. Whether the joint update fixes the mixing is not tested. The user will
revisit this depending on the story the paper tells: it is needed if the paper
claims, on the exact sampler at M0, that ptsG shrinks under metabolites only.

**It follows that:**
- at M0 the exact sampler does not measure the load-bearing cell. Until it does,
  that cell rests on track A's estimate (`approximate-method.md`, A1: 0.47 of prior
  SD against a bracket of 0.68) and on the mixed route's fit;
- B7.8 compares the stream arms on the parameters the exact fits mixed. The
  mixed route's ptsG shrinkage under metabolites only is reported alone, without
  an exact counterpart.

*Sections:* §11 16b.4.

### 2026-10-28 — 16b.1 and 16b.4 as run: 16b.1 recovers; 16b.4 is held open

**Trigger:** the three production chains finished (jobs 19034901–19034903 at
`03b4c4f`, 1,669.7 CPU-h). Their summary (job 19378382 at `626a1cc`) is in
`dev/scripts/m0_16b_result.md`.

**16b.1 recovers, by the rule fixed in advance.** ENO's shrinkage is 0.258, the
moved truth is at quantile 0.099, and the mean shift ÷ truth shift is
0.951 ± 0.021. K4 does not fire. 16b.1 is ticked.

**16b.4: what was measured.**
- **Transcripts only:** ENO's null is not flagged, by 0.00014. SD / prior SD is
  0.93445 against a flag line of 0.93431. The margin is 0.4% of one SE, far
  smaller than the uncertainty in the ESS estimate, so the verdict is not robust
  either way. FBA's is 1.007.
- **Metabolites only:** ENO, FBA and σ_b mixed (ESS 140 to 568), and match the
  joint fit.
- **Metabolites only, not measured:** S_0607 and ptsG have bulk ESS 7 and 6 and
  split-R̂ 1.13 and 1.16, so ptsG's load-bearing cell is not measured.
- **The symptom named in advance did not occur.** The named risk was the
  bootstrap proposal, checked by a live cell-window updating at a rate below 0.10.
  Every live cell-window updated at a rate of 0.37 or more (job 19376253). That
  shows the paths are resampled, not that they decorrelate from S_g, so whether
  the bootstrap move causes the slow mixing is not tested.
- **The early check was not run as fixed.** §12 2026-10-22 set `rates` at about
  sweep 20 on each chain. No such run exists. The only `rates` run is job
  19376253, on met_only after 656 sweeps. Before it, the per-sweep log lines
  ("windows changed") were read up to sweep 6.
- **Running on costs more than is left.** ESS 100 would take about 7,300 CPU-h
  more, against about 4.6k left. That assumes ESS grows in proportion to the
  draws, which an ESS of 6 does not establish.

**Change (the user's decision, 2026-10-28, after the results):** close 16b.1 now
and hold 16b.4 open. Its options are:
- report the metabolites-only transcript parameters as unmixed;
- redesign that arm's sampler move, then rerun it.

**It follows that** F6 is not drawn until that is decided. B7.8
(`approximate-method.md`) reads the same chain, so its ptsG comparison waits on the
same choice. B7.8 also needs a new floor: each arm is one chain, so the planned
chain-vs-chain W1 does not exist.

*Sections:* §11 16b.1, 16b.4.

### 2026-10-22 — 16b.1 and 16b.4 at M0: the runs and their reading rules

**Status: approved 2026-10-23, before any 16b production chain has run,** with the
ratio band of 0.5 to 1.5 and all three chains launched side by side. It fixes the
runs 16a.11 funded, and how each is read, before any of their draws exist. Only
the smoke runs (3 cells, 120 s, 2 sweeps) had been read.

**The runs.** Three production-length chains, side by side, at F10's settings:
- 50 cells, 600 s, N = 20, 656 sweeps, 5 threads, and 100 sweeps of burn-in;
- each starts from θ and σ_b drawn from the prior with `Xoshiro(16_400_000 + i)`,
  i = 1, 2, 3, and from paths drawn by `gibbs_init`;
- driver: `dev/scripts/m0_16b.{jl,sh,slurm}`, about 3 × 510 CPU-h.

| Mode | Task | Dataset | Streams |
|---|---|---|---|
| `perturbed` | 16b.1 | M0 replicate 16085, ENO's forward constant × exp(0.1568), its reverse derived | Both |
| `tx_only` | 16b.4 | V6's (replicate 16085) | Transcripts |
| `met_only` | 16b.4 | V6's | Bulk metabolites |

**What changes in `src`.**
- `m0_dataset(seed; forwards)` moves named forward constants of the drawn truth,
  each with its reverse derived (`derived_ode_values`). Every draw is the
  replicate's own, so the perturbed dataset shares V6's cells, σ_b and noise
  seeds. Only ENO moves, so the comparison is paired.
- `gibbs_setup(ds; streams)`. With the default, a chain is unchanged.
  - **Transcripts only:** the panel is empty, so no bulk term enters any block. σ_b
    is held at its start and is not inferred.
  - **Metabolites only:** block 3 bridges no gene. It simulates transcription
    forward with every other reaction and weights each window by the bulk term
    alone, which makes it a bootstrap particle Gibbs. Block 2 still reads the
    transcription counts from the path.

**16b.1's rule (K4's gate).** Read off `m0_16b.jl summary` for ENO, against V6
(sweeps 101 to 1,650). **It recovers if all three hold:**
1. **Shrinkage:** posterior SD / prior SD (0.157) is below 0.9.
2. **Coverage of the moved truth:** ln ENO = 4.3386 lies in the posterior's central
   95%.
3. **It moves with the truth:** the posterior mean's shift from V6's, divided by the
   truth's shift (+0.1568), lies in 0.5 to 1.5.
   - Under a conjugate reading at V6's shrinkage (0.27), the expected ratio is
     about 1 − 0.27² ≈ 0.93.
   - It is not exactly 1, because moving ENO changes the pools. That changes the
     rebuilt transcription rates and the paths, not only the bulk data.
   - **Why not a Monte Carlo error rule:** the means' MC SEs are about 0.005, far
     below that systematic difference, so such a rule would fail a correct
     sampler.
- V6's truth's quantile under the perturbed posterior is reported beside them.
- **If any of the three fails, K4 fires:** the phase stops and the cause is
  diagnosed before anything is rerun.

**16b.4's reading.** A measurement, not a gate.
- **F6:** the five parameters under transcripts only, metabolites only and joint
  (V6). Each with mean ± SD, shrinkage, bulk ESS, the truth and its quantile.
- **The load-bearing cell:** ptsG (`S_JCVISYN3A_0779`) under metabolites only,
  reported with its shrinkage against A1's linearised 0.47 at σ_b = 10%
  (`sensitivity_88_m0_result.md`). A1 is a bound under its own settings, not a
  prediction of this posterior.
- **The measured null:** ENO under transcripts only. The 2026-10-05 amendment
  bounds its ratio near 1 (`eno_path_information_16a7_result.md`; Core A′ at 200
  cells, at 15.7's truth).
  - A ratio below 1 − 2·SE, with SE = ratio / √(2 · bulk ESS), is flagged as a
    sampler error to investigate, not reported as a result.
  - The bound's scope (one truth, Core A′) is stated beside it.
- **For 16b.5:** M0's transcripts-only prior-to-posterior divergence for ENO and
  FBA is read off the same chain, labelled as M0's.

**The early check.** At about sweep 20, run `m0_16b.jl rates` on each chain.
- `met_only`'s bootstrap proposal is the mixing risk.
- A live cell-window below 0.10 (the pilot's threshold) is reported to you before
  the chain continues.
- Otherwise the chains run to 656. A chain that hits its time limit is resubmitted
  with the same command.

*Sections:* §11 16b.1, 16b.4, 16b.5.

### 2026-10-22 — 16a.11: the measured cost, K1's bound, and closing phase 16 at M0

**Status: approved 2026-10-22, with both calls answered "yes".** K7 (16b.6) is not
funded.

**Scope.** This is #90's task C2, "O1 and O5, at 16a.11", in the 2026-10-19
framing note (on #87's branch; parent §12 will fold it).
- **Kept separate from tracks A and B.**
  - It edits only this entry, and changes no ledger row; 16b's runs will add theirs.
  - Track A (#88, PR #92) and track B (#89, PR #93) keep their own ledger rows and
    sub-spec (`spec/phases/approximate-method.md`).
  - This entry funds none of their tasks, and claims none of their results as its
    own.
  - It reads two of them as inputs: A1 (`sensitivity_88_m0_result.md`) and B3's KG
    verdict (`jensen_gap_89_result.md`).
- **Not decided here:**
  - O2, the claim;
  - O6, the solver;
  - the folding of the 10-19 notes.
  The rest of #90 C2 owns those.

**1. The measured cost of one M0 posterior (production: PG, N = 20, 50 cells, 600 s).**

| Run | Sweeps | CPU-h | Wall | Source |
|---|---|---|---|---|
| F10's chain 5 | 656 | 507.5 | 4.2 d at 5 threads | `f10_16a10_result.md` |
| V7, mean of 50 replicates | 656 | 678 | 4.5–6+ d at 5 threads (16 resumed) | `v7_16a9b_result.md`, ledger |

- **Per sweep:** about 0.77 to 1.03 CPU-h at 5 threads. The spread comes from θ:
  draws that make longer paths cost more.
- **A production posterior holds 30 to 500 effective draws per parameter.** The
  slowest are the GAPD promoter, `krnadeg` and, in F10, ENO (`f10_16a10_result.md`).
- **A reference posterior (V6: 4 × 1,650 sweeps, at 20 threads)** cost about
  7.7k CPU-h, counting its pilot and stages (ledger: 896 + 2,308 + 4,039 + 463).

**2. The full-scale projection.** This multiplies measured parts. It is not a
measured full-scale sweep.
- **The parts:**
  - 16.5 cycles per cell-sweep at N = 20 (16a.7c, on M0 and on full-scale cell
    30001);
  - 32.2 s per warm 6,300 s cycle (13.7);
  - 200 cells, updated in turn.
- **One sweep:** about 30 CPU-h, or about 1.5 h of wall with each cell's 20
  particles in parallel (D16.6).
- **Mixing at full scale is unmeasured.** Taking M0's slowest rate, 0.17 bulk ESS
  per sweep:
  - **A production-length posterior** (100 burn-in plus 100 effective draws, about
    690 sweeps): about **21k CPU-h**, and about 6 weeks of wall.
  - **A V6-quality reference** (bulk ESS ≥ 1,000 over 4 chains, about 6,300 sweeps):
    about **190k CPU-h**.
  - **50 production-length replicates for coverage:** about **1.0M CPU-h**.
- **The solver changes the size, not the shape.** FBDF + KLU, if adopted (O6, not
  decided here), cuts each figure by up to 3.4×. The production posterior falls to
  about 6k, and coverage stays about 300k.
- **A measured full-scale sweep** would cost about 30 to 60 CPU-h and 2 to 4 h of
  wall. It is not proposed: it would firm up the 30 CPU-h, not the mixing, and the
  mixing is the unknown.

**3. K1's bound, set after the measurement** (parent §8, which T3 notes).
- **Proposed bound:** 50 calibration replications must fit inside phase 16's cap. So
  K1 fires where one production posterior costs more than 50k / 50 = **1,000 CPU-h**.
- **At M0:** 510 to 680 CPU-h, so K1 does not fire. V7 ran the 50.
- **At full scale:** about 21k CPU-h projected, so K1 fires.
- **K1's ladder** (parent §8): shorten the horizon, then five genes, then calibrate
  on the reduced reference system only. Phase 16 lands on the **last rung**: calibration
  on M0.
  - **Why the middle rungs are skipped:** a shorter horizon or five genes still means
    200 cells updated in turn. Each cuts the cost by roughly its own factor, about
    2 to 10×, which is far short of the ~20× needed.
  - T3 labels the full-scale verdict **budget-bound**.
- **The headline is not restated** to "coverage at a single truth". The ladder did
  not bottom out, because V7's coverage is at M0's rung over 50 truths.
- **Parent §8's K1 text is not edited by this entry.** On approval, a line is added
  there, citing this entry.

**4. The 16b allocation at the M0 rung.**
- **Phase 16's remaining cap:** about 7.4k CPU-h on `main`.
- **If #92 and #93 merge:** about 7.1k. Their rows, 279.2 and 4.3 CPU-h, are
  charged to phase 16's cap by those PRs, not by this entry.

| 16b task | At M0 | New CPU-h | Status |
|---|---|---|---|
| 16b.1, the tight-prior control (K4) | Recovery: V6's posterior, read off. Perturbation: one new M0 dataset with ENO's truth moved by +1 prior SD, and one production posterior on it | ~700 | **Fund** |
| 16b.2, recovery of the full set | V6's posterior. M0 has five of the six; `S_PGI` is not in M0 | 0 | Reuse |
| 16b.3, coverage | V7's 50 replicates, at 0 CPU-h (your call 1, below) | 0 | Reuse, if you agree |
| 16b.4, F6 at M0 | Transcripts-only and metabolites-only fits on V6's dataset; the joint fit is V6 | ~1.4k | **Fund**, after checking that the Gibbs driver can drop a stream (it may need a small code change) |
| 16b.5, K2 | The 1,000-cell extension is not run. At M0, K2's divergence is read off 16b.4's transcripts-only fit at 50 cells, labelled as M0's, not as K2's verdict | 0 | **Budget-bound** at 1,000 cells |
| 16b.6, K7 (PYK) | PYK is not in M0. It needs PYK freed and two fits | ~1.4k plus code | **Not run, budget-bound**, unless you fund it |
| 16b.7, 16b.8 | Labels and T3 | 0 | Writing |
| Resumes, smoke runs, suite runs | | ~500 | **Fund** |

- **The funded total is about 2.6k CPU-h,** which leaves about 4.5k, or about 3k
  with K7 funded.
- **The new posteriors run side by side,** in about 5 to 6 days of wall.
- **What else could use the remainder, none of it claimed here:**
  - track C's demo (#90 C1, about 50 to 500 CPU-h);
  - track B's next step.

**5. O1 and O5, given what tracks A and B have found.**
- **O1, part 1: close phase 16 at M0.** Proposed, as above.
- **O1, part 2: reopen the cap for paper section 5.3's injections.** Not proposed now.
  - Those injections test the approximate method. B3 found the mean-field
    description fails its Jensen-gap gate (KG) at every truth, and B4 waits on a
    rethink (`jensen_gap_89_result.md`).
  - Reopen when an approximate method passes 5.2. That is a later §12 entry.
- **O5, one exact Core A′ posterior as a full-scale yardstick:** about 21k CPU-h, or
  about 6k with FBDF + KLU. It exceeds the remaining cap either way.
  - **Recommended: not now.** Its value is as a yardstick for an approximate method
    at full scale, and there is no such method yet.
  - Revisit with O1 part 2.
- **A1's input:** at M0's size, metabolites carry a ptsG signal: 0.47 against a
  bracket of 0.68 at σ_b = 10% (`sensitivity_88_m0_result.md`, on PR #92). So 16b.4's
  metabolites-only fit has something to find. A1 is a linearised bound, not a
  posterior.

**6. Your two calls.**

**Call 1: does V7 also serve as 16b.3's coverage result?**
- **The evidence:** coverage read off V7's ranks (`dev/scripts/v7_coverage_16a11.py`;
  R = 50, L = 100), against the binomial 95% band:

  | Level | GAPD | ptsG | `krnadeg` | ENO | FBA | σ_b | Band |
  |---|---|---|---|---|---|---|---|
  | 50% | 23 | 27 | 22 | 31 | 24 | 26 | 18–32 |
  | 80% | 36 | 41 | 37 | 42 | 40 | 42 | 34–45 |
  | 90% | 42 | 47 | 45 | 45 | 42 | 48 | 41–49 |
  | 95% | 46 | 47 | 46 | 46 | 44 | 49 | 44–50 |
  | 99% | 49 | 48 | **47** | 50 | 49 | 49 | 48–50 |

  - 29 of 30 cells are in band. One or two misses are expected by chance over 30.
  - The miss, `krnadeg` at 99%, is where L = 100 is coarsest: three truths sit at an
    extreme rank.
- **For:**
  - It meets 16b.3's "≥ 50 replications, with the count and rung stated".
  - It is the same sampler at production settings.
  - It costs nothing.
- **Against, and to be stated with it:**
  - 23 of the 50 replicates have bulk ESS below 100.
  - L = 100 makes the 99% level coarse.
  - It is M0's five parameters plus σ_b, not D11's six.
  - **K4 asks for the tight-prior control at 100 replications.** V7 has 50. On
    approval, K4's band is restated at R = 50 here, as 16b.3 requires: 44 to 50 of
    50 at nominal 95%. ENO covers 46, inside the band.
- **Recommended: yes,** with those four statements carried into T3.

**Call 2: does "closed at M0" meet phase 16's done-when?**
- **The parent's condition, read at the M0 rung:**
  - coverage within Monte Carlo error for every target: V7, call 1;
  - the tight-prior control shrinks and recovers. V6's posterior SD for ENO is
    0.042 against a prior SD of 0.157, a ratio of 0.27, below K4's 0.9 and its 0.5.
    16b.1 adds the perturbation;
  - the loose-prior control is visibly wider. FBA's posterior SD is 0.184 against
    ENO's 0.042, and its ratio to its prior is 0.33 (prior SD 0.558).
- **Proposed:** yes, at the M0 rung, once 16b.1 and 16b.4 land.
  - The parent's phase 16 is ticked as **"met at M0; full scale budget-bound"**,
    never as "met".
  - 16b.7 carries that label on every result.
- **The alternative:** leave phase 16 open until a full-scale posterior is funded.
  That ties it to O5.

**On approval:**
- 16a.11 is ticked.
- Parent §8 K1 gets its bound line.
- The 16b tasks are annotated with their rung and funding.
- 16b's branch starts with 16b.1's perturbed dataset and 16b.4's stream check.

*Sections:* §11 16a.11 and 16b, D16.6, parent §8 K1 and K4.

### 2026-10-22 — F10's divergence, and a floor matched to production's length

**Status: recorded 2026-10-22, before any F10 divergence was computed.** The
production chain (job 18751635) had finished, and only its sweep log had been
read.

**Trigger:**
- V8 says "the divergence" without naming one.
- V8's floor, the reference's split halves, is too low to read a single
  production chain against. Each half holds about 540 effective draws on the GAPD
  promoter. The production chain holds about 97. So that floor sits below the
  Monte Carlo noise a production chain carries by construction, and V8 alone would
  read that noise as a difference.

**The divergence.**
- **Measure:** per marginal, the 1-D Wasserstein-1 distance between the two
  draw sets, on the log scale, divided by the reference's posterior SD.
- **Reported beside it,** as description only: the difference in means in
  reference SDs, and the ratio of SDs.
- **Production arm:** sweeps 101 to 656 of chain 5.
- **Reference:** sweeps 101 to 1,650 of V6's chains 1 to 4, as in V6's merge.

**The floors.**
- **V8, as written:** chains 1 and 2 pooled against chains 3 and 4, each on sweeps
  101 to 1,650.
- **Added, matched to production:** sweeps 101 to 656 of each V6 chain against the
  other three chains' sweeps 101 to 1,650.
  - Each V6 chain started from the prior in the same way chain 5 did (the pilot's
    `Xoshiro(16_090_000 + k)`), with the same code path and burn-in. So its first
    656 sweeps *are* a production-length chain. That gives four null draws of the
    divergence.

**The reading rule.**
- **Per parameter,** chain 5's divergence is ranked among the five: itself and the
  four null draws.
- **If chain 5 is a production-length chain like the others,** its rank is uniform
  on 1 to 5. Across six parameters, about 1.2 come out largest by chance alone.
- **"Not resolved at production length":** the divergence is not above every
  null draw.
- **Above all four:** the divergence is reported with its margin. It counts as a
  finding only if it is large against the nulls' spread, or if several parameters
  do it. The F10 result file judges this, and states its reasoning.
- **No gate.** F10 is a figure. A production divergence well beyond the nulls
  would be raised with the user before 16a.11, not silently absorbed.

*Sections:* V8, §11 16a.10.

### 2026-10-13 — F10 when the reference runs at production's N

**Status: approved 2026-10-18, at 656 sweeps.** The alternatives weighed at approval:
- 700 sweeps, rejected because 16b's chains run 656, and F10 should measure the
  length production uses;
- the zero-cost hold-out, rejected because it leaves the reference below V6;
- several fresh chains (about 2k CPU-h), left out of the remaining cap;
- a true N_ref ≫ N reference, which is unaffordable.

One production chain gives one draw of the divergence. A divergence inside V8's
floor is read as "not resolved at production's length". *Superseded 2026-10-22
(§12):* the reading is against four nulls matched to production's length, not V8's
floor alone.

**The run:** `gibbs_pilot_16a9a.jl run 5 656` from `f10_16a10.sh`. That is V6's
own setup (replicate 16085, 50 cells, N = 20) as a fifth chain, started from the
prior with `Xoshiro(16_090_005)`, at 5 threads, in its own directory so that no V6
merge reads it. `src/` changed only in docstrings between V6's commit `281d548` and
this run. Job 18751635 at `a0060de` runs it. Its first checkpoint records the same
replicate, cells, N, truth and σ_b as V6's `chain1.jls`, and a different start.

**Trigger:** V6 is met (`v6_16a9b_result.md`, merge 18420990) with PG at N = 20,
production's own setting, over 4 chains of 1,650 sweeps. That is just past V6, not
"far longer". D16.5 asked for N_ref ≫ N and chains far longer than V6 requires.
Neither holds.

**Why F10 still measures something:**
- D16.5's own case analysis says that when production is PG, the particle count
  costs mixing, not bias. F10 then measures only the Monte Carlo and mixing error
  of production's settings.
- V6 checks the reference's convergence directly, by R̂, ESS and quantile MCSE. So
  the reference's accuracy does not rest on N_ref.
- **What F10 compares is the run length, not the sampler.** A production-length
  chain is set against the long reference on the same dataset. That is the
  question 16b's single chains raise.

**Proposed change:**
- **Production arm:** one fresh chain on V6's dataset, at V7's production length
  (100 sweeps of burn-in plus at least 100 effective draws, about 656 sweeps), at 5
  threads. About 500 to 700 CPU-h and about 6 days.
  - **Corrected 2026-10-17, before approval (#86 review).** 656 sweeps does not
    give 100 effective draws. On V6's dataset the slowest bulk rate is 0.174 per
    sweep (the GAPD promoter), so 556 kept sweeps give about 97. In V7, 23 of 50
    replicates fell short of 100 at 656 sweeps (`v7_16a9b_result.md`).
  - **The proposal is now 656 sweeps, stated as such.** That is the length
    production runs, which is what F10 compares. It gives about 97 effective
    draws on the slowest parameter.
  - **The cost** is about 510 CPU-h at 0.78 CPU-h per sweep. V7's sweeps took 450
    to 1,100 s, so the wall time is 3.4 to 8.4 days, and the chain may need one
    resumption past the 6-day limit.
  - **Alternative:** 700 sweeps for at least 100 effective draws, about 550 CPU-h.
- **Reference arm:** V6's 4 chains, with V8's floor from their split halves.
- **D16.5** is annotated: N_ref = N, justified by the PG case, and F10's caption
  says so.
- **Alternative, at no cost:** hold out one V6 chain's first 656 sweeps as the
  production arm against the other three. The chains are independent, but the
  reference falls to 3 chains and below V6's ESS of 1,000.

**Recommended:** the fresh chain. It is clean, and it fits in the remaining cap.

*Sections:* D16.5, §11 16a.10.

### 2026-10-12 — V7's pass rule: α/6 per parameter, on a complete set of replicates

**Status: approved 2026-10-12, before any rank was computed.** Two of the 50
replicates (33 and 40) had finished; nothing had been merged or read.

**Trigger:** `/check-PR` on #82.
- **The rule as written is too strict.** V7 has six parameters, and each gets a
  95% simultaneous band. If V7 passes only when all six are within, a calibrated
  sampler fails up to about 26% of the time. Under KR, a failure stops the phase.
- **The automatic merge would rank a biased set.** Its `afterany` dependency
  fires when the array ends, and ranks only the replicates that finished. Which
  replicates finish depends on their speed, and their speed depends on θ. So that
  set is not a sample from the prior.

**Change.**
- **The pass rule.** Each parameter's rank-ECDF band is at level
  1 − 0.05/6 = 0.9917. **V7 passes when all six are within their bands.** By the
  union bound, a calibrated sampler then fails V7 at most 5% of the time.
  - Checked by simulation at R = 50 and L = 100, on uniform ranks, from the merge
    script's own `ecdf_bands`. Per-parameter coverage is 0.9915. With six
    independent parameters, all six are within 95.15% of the time.
- **A complete set only.**
  - The merge refuses to run while any of the R replicates is missing or short of
    its sweeps.
  - The automatic merge 18276795 was cancelled.
  - The replicates that hit the 6-day limit are resubmitted from the v7 worktree
    at `dfe570a`. Their run path is unchanged.
  - The merge runs from the commit that carries this entry.
- **Slow-mixing replicates are named, and the verdict is still on all R.**
  - A replicate whose smallest bulk ESS over its kept sweeps is below 30 is
    listed. Below 30, an AR(1) chain's thinned ranks start to widen.
  - The bands are read again without the listed replicates, as a sensitivity.
  - That reading is reported beside the verdict but does not decide it.
  - Excluding a replicate for slow mixing is not done: mixing depends on θ, so
    exclusion would bias the set as an incomplete set does.

*Sections:* §3 V7, §8 KR, §11 16a.9b.

### 2026-10-08 — 16a.9b's size: V6 at 20 threads, V7 at 5, both with 100 sweeps of burn-in

**Status: approved 2026-10-08.** This decides option 3 of §12 2026-10-07.

**Trigger:** V6 stage 1 and option 2's measurements.
- **The files:** `dev/scripts/v6_stage1_16a9b_result.md`,
  `profile_b3_16a9b_result.md` and `thread_scaling_16a9b_result.md`.
- **The jobs:** 18157561 to 18157565, 18271545, 18255035, 18256436, 18256437,
  18258719 and 18258720.

The measurements:
- **Burn-in.** With 100 sweeps dropped, R̂ is at most 1.012, and the posterior
  means agree with a 350-sweep burn-in to 0.04 SD.
- **What V6 still needs.** The quantile MCSE binds, at 5,807 kept sweeps.
- **The cost of a sweep falls with fewer threads,** and the chain is bitwise
  the same at 1, 5, 10 and 20 threads. Per sweep it is 1.22 CPU-h at 20 threads,
  1.00 at 10, 0.78 at 5 and 0.80 at 1. Block 3's parallel efficiency at 20
  threads is 40 to 51%.

**Change.**
- **V6:** the four chains continue at 20 threads to 1,552 sweeps each, with 100
  sweeps of burn-in. That is about 4.1k CPU-h and 2.2 days.
- **V7:** R = 50 replicates. Each is one chain from a prior-drawn start, at 5
  threads, with 100 sweeps of burn-in plus at least 100 effective draws. That is
  about 656 sweeps per replicate, at the slowest ESS rate of 0.180 per sweep, and
  about 25.6k CPU-h.
- **Projected phase total:** about 33.4k of the 50k.
- **D16.6 is unchanged.**

**Annotated 2026-10-09: V7 runs slower than sized, and continues.** You decided
to let both V6 and V7 run.
- **The rate.** After 55 to 130 sweeps per replicate (array 18276794), the mean
  sweep takes about 784 s at 5 threads, against the 562 s of
  `thread_scaling_16a9b_result.md`. Across replicates it ranges from about 450 to
  980 s.
- **The cost.** V7 now projects to about **35.7k CPU-h**, against 25.6k, and the
  phase to about 43k of the 50k.
- **The time limit.** About half the replicates will reach the 6-day limit, and
  each is resubmitted once from its checkpoint.
- **The cause is untested.** Every node is fully loaded, and each replicate's
  dataset differs.

**Annotated 2026-10-17: V7 is done, and 23 of 50 replicates fell short of 100
effective draws.**
- **The shortfall.** The smallest bulk ESS over the kept sweeps ranges from 15 to
  325, median 109. The 656 sweeps came from V6's one dataset. Across prior-drawn
  θ, some replicates mix more slowly than that.
- **The verdict is unaffected.** The verify clause asks for R ≥ 50 at production
  settings, not a per-replicate ESS. §12 2026-10-12's rule passes on all 50.
- **The check is less sensitive.** The loss is recorded in
  `v7_16a9b_result.md`.
- **The cost:** 33.9k CPU-h, against the 35.7k projected above.

*Sections:* §11 16a.9b.


### 2026-10-07 — 16a.9b as written costs more than the phase cap

**Status: approved 2026-10-07** (options 1 and 2 now, option 3 deferred). Option 1
runs as the pilot's 4 chains extended from 200 to 700 sweeps each. That is the
same chains from the same overdispersed starts, so the 500 new sweeps cost about
2.2k CPU-h. The code is unchanged since `6f8bcbc` outside the test extras'
compat.

**Trigger:** 16a.9a's pilot (jobs 18108335 to 18108339;
`dev/scripts/gibbs_pilot_16a9a_result.md`).
- **The rate.** The ptsG promoter mixes slowest, at bulk ESS 0.108 per sweep.
  It rests on a bulk ESS of 43, so it is rough.
- **The cost.** One sweep takes 199 s at 20 threads, or 1.11 CPU-h. B3 is 126 s
  of that, and core use is 62%.
- **V6 at N_ref = 20** needs about 10,300 sweeps over 4 chains, set by the
  quantile MCSE. That is **11.4k CPU-h and 6.0 days**.
- **V7 at R = 50** needs 1,028 sweeps per replicate, so **57k CPU-h**. R = 100 is
  114k.
- **Together, about 68k CPU-h.** The phase cap is 50k (D16.6), and phase 16 has
  spent about 1.3k. 16a.9's hours are not held to 16a's 5k (§12 2026-10-06), but
  the phase cap still applies.

**Options:**
1. **Stage V6.** Submit the 4 resumable chains for 500 sweeps each, about 2.2k
   CPU-h and 1.2 days. Then re-project from the measured rate and extend. V6's
   criteria are unchanged.
2. **Speed the sweep first.** Profile B3 and B1. Any gain scales V6 and V7 alike.
   The size of the gain is unknown until measured.
3. **Resize V7.** Each of these is a decision about what V7 certifies:
   - **Fewer effective draws per replicate.** At 50 rather than 100, R = 50 costs
     about 31k CPU-h.
   - **Fewer cells per replicate.** Cost falls roughly with C, but V7 would then
     not calibrate the sampler at M0's C = 50.
   - **Raise the phase cap,** at 16a.11 or now.
   - **Defer V7's sizing** until the refined V6 rate and any speed-up are in.

**Proposed:** options 1 and 2 now. Option 3 waits until the V6 rate is refined
and any speed-up is measured, and that choice comes back to you.

*Sections:* §11 16a.9a and 16a.9b; D16.6 (not yet changed).


### 2026-10-06 — 16a.9 splits: the Gibbs driver and a pilot first, then the reference

**Status: approved 2026-10-06.**

**Trigger:** planning 16a.9. The three blocks exist and are validated one at a
time (V3, V4, V5), but nothing composes them into a chain, and the chain's ESS per
sweep is unmeasured. A projection from unit costs puts V7 alone at about 11,000
CPU-h and a serial 200-sweep chain at about 5.5 days, both resting on a guessed
ESS per sweep.

**Change.**
- **§11:** 16a.9 is annotated in place and split. 16a.9a builds the driver, its
  exactness tests, threads in B3's particles and B1's cells, and a pilot that
  measures mixing and cost. 16a.9b runs V6 and V7 at the size the pilot sets.
- **D16.6:** 16a.9's CPU-h are recorded in the ledger but not held to 16a's 5k
  cap, by your decision. Multi-day waits are flagged before submission.
- Handoff follow-ups 1 and 3 (`_poisson_weights`' stop, `_draw_index` in the
  fallbacks) land in 16a.9a, since the pilot measures cost.

*Sections:* §11 16a.9, 16a.9a and 16a.9b; D16.6.


### 2026-10-05 — bulk metabolites, and the seam's ENO cell as a measured null

**Status: approved 2026-10-05**, with the parent's amendment of the same date,
which carries the trigger in full. The reviewed text is
`spec/amendment-draft-2026-10-05.md`.

**Trigger, in short:**
- KP fired twice. The per-cell panel is unmeasurable at any cited single-cell
  detection limit.
- A cut at the history would remove F6's ptsG cell.
- Bulk metabolites keep it: posterior over prior SD 0.10 at a 10% CV.
- A full-path information bound refutes the ENO cell: 99.9% of the prior over 200
  cells.

**Change.**
- **§1:** the data are bulk metabolites.
- **§3:** the target uses `x̄` and σ_b. B1 scores `x̄`, and B3 updates cells in turn.
  The independence assumption is struck.
- **V5:** the "on" half is on C = 3 cells.
- **V7:** draws σ_b.
- **D16.3:** the substituted-mean weight and the sequential scan; annealing is kept,
  off.
- **§8 KP:** annotated as resolved by the observation model.
- **§11:** 16a.7a is superseded for production, 16a.7c is added, and 16b.4 reports
  ENO against its bound.

**Corrected after the #81 review** (2026-10-05/06):
- **N = 10 → N = 20.** The update rate quoted first was pooled over cells. Per cell,
  over windows whose true path has events, N = 10's lowest rate is 0.12 from 16
  sweeps, which is not robust.
- **D16.6** gains the wall-clock of updating cells in turn, and **R16.8** is added.
- **§1, D16.4, D16.6's allocation and 16b.4** carry the one-way claim.
- **σ_b replaces σ** where it remained.
- **16a.1 to 16a.5 and 16a.8 are ticked** with their evidence.

*Sections:* §1, §3, §8, D16.3, D16.4, D16.6, §10 R16.8, §11 16a.1 to 16a.8, 16a.7a,
16a.7c and 16b.4.

### 2026-10-02 — KP fired; each window's proposal is annealed

**Status: approved 2026-10-02.** You chose this over revisiting the observation
floor or conditioning the path on transcripts alone.

**Trigger:** task 16a.7's variant measurement (jobs 17878617 to 17878621, 17882309
to 17882311 and merge 17884470; `dev/scripts/csmc_variants_16a7_result.md`).
- **The rates.** At N = 50 on M0's five cells over 40 sweeps, PG, ancestor sampling
  at lag 1 and 5, and full PGAS each leave a window that never changed. On
  full-scale cell 30001, PG changed none of 840 window-sweeps.
- **The cause** (jobs 17879029 and 17879030). The transcript part of the weight is
  the same across particles. The metabolite part puts proposals a mean 570 to
  1,870 nats below the reference on M0, and 40 to 400 at full scale, within one
  window. `M_pi_c`, `M_3pg_c` and `M_2pg_c` carry it: they drain toward empty, and
  their response to the path is nonlinear. Every reaction's count and exposure
  explain 50 to 59% of it on M0, and protein counts at most 16%.
- **Not a proposal error on M0.** The true path's ranks among proposals are uniform
  over 50 windows (job 17880206). At full scale they are a watch item (job
  17880205).
- **A bug, fixed at `2566e79`.** Ancestor sampling's future replay threw on a firing
  with zero propensity. That firing now gives density zero.

**Change.**
- **D16.3:** each window's proposal is AIS over π_β ∝ q · w^β, with resample-the-tail
  MH moves and a conditional backward chain for the reference, under particle
  Gibbs. Its cost, about N(K + 1) cycles per sweep, is stated, and is unfundable at
  full scale under the cap.
- **§8 KP:** annotated as fired, with the remaining fallback named.
- **§9:** the weight-degeneracy question is resolved, and a full-scale proposal
  question is added.
- **§10:** R16.7.
- **§11:** 16a.7 is annotated, and 16a.7a and 16a.7b are added.

**Rejected:**
- **More particles or ancestor sampling.** The gap is hundreds of nats inside one
  window.
- **A shorter window.** The data arrive every 60 s.
- **K1's ladder.** M0 fails in window 1.
- **A linear tilt toward the observed pools.** It explains at most about half the
  variance on M0.

*Sections:* D16.3, §8, §9, §10, §11 16a.7, 16a.7a and 16a.7b.

### 2026-10-01 — the bridge cap's +20 is a floor, and M0 is built before block 1

**Status: approved 2026-10-01.**

**Trigger:**
- **The cap.** Task 16a.3 (job 17758372 at `5d4bd8c`, then job 17758553 at
  `d9b4aa5`; `dev/scripts/bridge_block2_16a34_result.md`). At k = 1/s and
  μ = 0.2/s, over 20 s from a count of 10, max(observed) + 20 = 30 drops 7.6e-12,
  against D16.3's 1e-12. At the 15.7 truth's rates the rule drops 2.4e-20. So the
  rule holds where the data sit, and fails at the amplified rates V5 needs.
- **The order.** 16a.5's check, V4, runs block 1 on M0. V1's 10⁴-path importance
  identity needs M0 too. M0 is task 16a.8. Running V4 on full Core A′ cells instead
  costs hundreds of CPU-hours, because every slice evaluation replays three
  6,300 s cells.

**Change.**
- **D16.3:** max(observed) + 20 is the cap's floor. The cap grows in steps of 5
  until the truncation mass is below 1e-12 at the constants in use
  (`adequate_cap`). 16a.3 is annotated.
- **§11:** 16a.8 is worked before 16a.5, and V1's 10⁴-path identity runs on M0
  once it exists. No task's content changes. 16a.2, 16a.5 and 16a.8 are annotated.

*Sections:* D16.3, §11 16a.2, 16a.3, 16a.5 and 16a.8.

### 2026-09-30 — V2b's SSA clause: the random stream is not the driver's, and the cached next jump is dead at a boundary

**Status: approved 2026-09-30.**

**Trigger:** task 16a.1 (jobs 17754678 and 17754679,
`dev/scripts/path_replay_16a1_result.md`). V2b asked for a snapshot to carry the
SSA's random stream and its pre-drawn next jump, and for dropping the latter to
fail the check. Neither is possible in the driver as built:
- The SSA draws from JumpProcesses' `DEFAULT_RNG`, which is `Random.default_rng()`,
  the task's global generator. A deep copy of the driver does not copy it.
- Every handshake calls `reset_aggregated_jumps!` before the jump step, which
  redraws the next jump and writes it to `integrator.tstop`. So the cached value is
  dead at a window boundary, and leaving it out could not fail.

Replay draws nothing, so V2 and V2b hold for everything a replay uses.

**Change.**
- **V2b's SSA mutation is replaced by a test** that corrupting a restored
  particle's cached next jump changes nothing when the random stream is seeded
  alike.
- **Per-particle streams move to 16a.7.** A particle that simulates forward sets the
  global generator from its own recorded seed before each forward step. 16a.7
  verifies that two restores of one snapshot seeded alike are identical, and that
  seeded differently they differ.
- **Rejected: giving the driver its own generator.** The 15.7 dataset was generated
  under the global stream, so it would no longer regenerate identically.

**Also found in 16a.1, not an amendment.** V2b asks for the restore at every
window boundary of one cell, and the first test did one boundary. The test and the
full-cycle driver now cover every boundary.

*Sections:* §3 V2b, D16.3 (particle state), §11 16a.7.
