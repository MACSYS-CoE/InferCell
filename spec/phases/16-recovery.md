# Spec: Phase 16 — Recovery, coverage, and the reference

**Status:** in progress: 16a.1 to 16a.8, with 16a.7a to 16a.7c, merged in PR #81 (2026-10-05); 16a.9a is in flight; 16a.9b to 16a.11 and 16b remain
**Created:** 2026-09-30  ·  **Last amended:** 2026-10-06
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
| V7 | Reference calibration | Simulation-based calibration of the M0 sampler at production settings. Each replicate's truth draws **all** free quantities from their priors: M0's targets by `draw_truth(...; names = M0_TARGETS)`, and σ_b from `p(σ_b)` (σ before 2026-10-05), each with a recorded seed. Rank-ECDF within 95% simultaneous bands for every free parameter, at R replications (≥ 50, K1's floor; 100 if the cap allows) | 16a.9 |
| V8 | F10's own floor | The divergence between the production and reference posteriors is reported beside the reference-vs-reference divergence of its split chain halves. A production divergence inside that floor is "not resolved", not "zero" | 16a.10 |
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
- [ ] 16a.6 (parent 16.3) Confirm a clipped drain under a fixed path still has a
  discontinuous derivative in an ODE parameter — verify by a test, as the parent's
  annotation specifies:
  - on the clipped model, a finite-difference slope in `ln kcat` that jumps across a
    clip and does not shrink with the spacing;
  - paired with 14c.5's continuity test on the smoothed model;
  - check 7's census cited as the source of the clip, not assumed.

  Recorded as the reason for D16.1.
- [ ] 16a.7 (parent 16.2) Build the conditional SMC path update, with D16.3's weight,
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
- [ ] 16a.7a (§12 2026-10-02) Build D16.3's annealed window proposal, under particle
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
- [ ] 16a.9 (split into 16a.9a and 16a.9b — see amendment 2026-10-06) (parent 16.1, run) Run the M0 reference long — verify by V6, which
  justifies the length, and by V7 at production settings on R ≥ 50 datasets, with R
  stated.
- [ ] 16a.9a (§12 2026-10-06) Compose the three blocks into one Gibbs chain on M0,
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
- [ ] 16a.9b (§12 2026-10-06; sizing proposed in §12 2026-10-07) Run the M0
  reference long — verify by V6, which justifies the length, and by V7 at
  production settings on R ≥ 50 datasets, with R stated. Sized from 16a.9a's
  pilot.
- [ ] 16a.10 (parent 16.8) Compare production against the reference on M0 — verify
  by F10: overlaid marginals, one 2-D contour for the ptsG and ENO pair, and a
  per-marginal divergence beside V8's floor. The caption states which case of D16.5
  holds.
- [ ] 16a.11 Measure cost, set K1's bound, and reopen the budget. Verify by:
  - the ledger showing 16a ≤ 5k CPU-h;
  - the cost per posterior on M0, and a full-scale projection per variant from
    measured sweeps;
  - K1's bound written into parent §8, with a note that it was set after the
    measurement;
  - a proposed allocation table funding each 16b task at a named rung, with any
    budget-bound item named;
  - a §12 entry;
  - **your approval of the budget and the allocation before 16b starts.**

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
**PR:** _not started_

- [ ] 16b.1 (parent 16.4) Recover the tight-prior control alone — verify by:
  - the posterior concentrating on truth well inside its prior;
  - a deliberate perturbation of the truth moving the posterior with it;
  - shrinkage below 0.9.

  **If it does not recover, the machinery is wrong and the phase stops (K4).**
- [ ] 16b.2 (parent 16.5) Recover the full six-parameter set — verify by:
  - every marginal containing truth;
  - at least one target shrinking below 0.5 (K4);
  - the cost per posterior reported against 16a.11's projection.
- [ ] 16b.3 (parent 16.6) Compute coverage over repeated datasets — verify by:
  - 50, 80, 90, 95 and 99% intervals covering truth at those rates within binomial
    error, per parameter, at ≥ 50 replications, with the count and rung stated;
  - **the tight-prior control at 100 replications**, inside K4's band of 90 to 99 out
    of 100 at nominal 95%. If 100 is not funded, K4's band is restated at the R that
    is run through a §12 amendment, never silently;
  - K4's third trigger (equal widths for ENO and FBA) checked.
- [ ] 16b.4 (parent 16.7) Produce F6 — verify by all six targets under
  transcripts-only, metabolites-only and joint data, with ptsG under metabolites only
  (load-bearing) and ENO under transcripts only (the measured null) present. Any cell 16a.11
  did not fund is shown as not run. **Amended 2026-10-05 (§12):** ENO under
  transcripts only is reported against its information bound
  (`dev/scripts/eno_path_information_16a7_result.md`) as a measured null. A
  shrinkage there is a sampler error to investigate, not a result.
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
  - every result's rung and budget-bound label.
- [ ] 16b.8 (parent 16.11) Fill T3 — verify by:
  - all seven criteria carrying a threshold, a measured value and a verdict,
    **published whether or not everything passed**;
  - K1's verdict written here, with its threshold marked as set after its
    measurement;
  - K5 recorded as fired in 14b and answered by 14c and D16.1;
  - every verdict the cap forced labelled budget-bound.

## 12. Amendment log

### 2026-10-07 — 16a.9b as written costs more than the phase cap

**Status: proposed 2026-10-07, awaiting your decision.**

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
