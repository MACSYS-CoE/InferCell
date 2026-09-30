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
published model, with exact transcript counts and lognormal metabolite panels. The
sampler is validated against an exact long-run reference on a two-gene system (M0).
So a calibration result can be attributed to the method rather than to a converged
wrong answer (R4). On Core A′ it delivers F6's two load-bearing cells, the tight-prior
control, K2's verdict and coverage. Each is delivered at a scale set by a measured
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
σ     = the metabolite noise scale
```

The reverse constants are deterministic functions of θ_ODE (D16.1). The target is

```
p(θ, σ, X_1:C | y)  ∝  p(θ_CME) p(θ_ODE) p(σ) · Π_c  p(X_c | θ) · 1[ m(X_c, t_k) = y^tx_c,k  ∀k ]
                                        · Π_k Π_j  N( log y^met_c,j,k ; log max(x_j(θ, X_c, t_k), 1), σ )
```

- **The priors.**
  - `p(θ_ODE)` is the product of the two forward-constant priors only.
  - `p(σ) = LogNormal(log 0.2, 1.0)`.
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
| B1 | `θ_ODE, σ \| X, y`: the prior, plus Σ_c [log p(X_c \| θ) + metabolite log-likelihood]. A new θ_ODE needs every cell replayed. σ alone needs no replay | C replays per θ_ODE evaluation |
| B2 | `θ_CME \| X, θ_ODE`: each promoter and `krnadeg` enters only its own propensities. With the pool trajectory fixed, the likelihood is `k^n exp(−k·A)`, where `n` counts events and `A` is exposure summed over cells. The prior is lognormal, so this is not gamma. It is a one-dimensional log-concave density | Negligible: sufficient statistics only |
| B3 | `X_c \| θ, y_c`, per cell and independent across cells | The path update of D16.3 |

### Assumptions, each with its regime

| Assumption | Holds when | What breaks outside it |
|---|---|---|
| A path and θ determine the observed trajectory exactly | Fractional carry and the clamped drain are deterministic, and initial transcripts are fixed | Stochastic rounding would add a latent. It is not the data policy. Asserted by V2 |
| Rate constants change at 60m − 1 s on the jump clock | The handshake order above | Nothing to fix: the path update is written for it (D16.3), and V1 tests that boundary |
| A particle is its driver's full mutable state | The snapshot covers every field `handshake_step!` reads or writes (D16.3) | A missed field makes a restored particle run a different model. Asserted by V2b |
| The turnover ceiling never binds | `RNAPOL_KCAT·S_g < 180`, that is S below about 153 | B2's form becomes truncated. Asserted per proposal, and a violation throws |
| Cells are independent given θ | Separate seeds and separate cells, as generated | Nothing here: true by construction |
| Block 1's step needs no gradient | Two free ODE constants plus σ | Freeing many more ODE parameters would bring the gradient, and 14c's sampler model, back. §7 |

### Verification checks, each with a tolerance

| # | Check | Asserted | Where |
|---|---|---|---|
| V1 | Path density | Three parts: the importance identity `E_{X~p(·\|θ0)}[p(X\|θ1)/p(X\|θ0)] = 1`, the hand checks, and the mutations below | 16a.2 |
| V2 | Replay fidelity | Replaying a generated cell's recorded path at its θ reproduces the carry-inclusive latent trajectory **bitwise** at every handshake and all 105 save points, on 3 dataset cells. If the replay changes the solver's step sequence, the fallback is ≤ 1e-8 relative, and the reason is recorded. **Mutation:** removing one translation event changes GTP | 16a.1 |
| V2b | Snapshot restore | Snapshot a replay at a window boundary, run another particle, restore, and continue. The result must equal the uninterrupted replay (bitwise, or 1e-8 with the reason recorded), at every window boundary of one cell. The intervening particle is constructed to differ in every snapshotted field, including a non-zero deficit, a different membrane ptsG count and a different SSA state, so no field passes vacuously. **Mutation:** leaving out `DeferredDebit.deficit`, `factor`, or the SSA's next-jump state fails it | 16a.1 |
| V3 | Exact transcript-only posterior | One gene, rate constants frozen (NTP pools chemostatted), transcripts observed every 60 s and nothing else. The exact posterior of `(S, krnadeg)` on a 2-D grid comes from the birth-death transition matrix. B2 with B3's bridge must match it: the 5, 25, 50, 75 and 95% marginal quantiles each within 3 Monte Carlo SE (ESS-based). **It cannot see a weight error that varies between particles, since the rates are frozen. V5 covers that** | 16a.4 |
| V4 | Block 1 on a fixed path | On M0, 3 cells, path fixed: B1's samples of `(ln kcatF_ENO, ln kcatF_FBA)` match a 2-D grid posterior computed by direct replay under the same forward-only prior. The grid density is evaluated on the θ scale, or on the ln θ scale with the Jacobian, and never as the same code path the sampler uses. Quantiles must agree within 3 MC SE. `assert_haldane` holds at every accepted state | 16a.5 |
| V5 | Path update against brute force, with the rates varying | See below the table | 16a.7 |
| V6 | Reference convergence | At least 4 chains from overdispersed starts. Rank-normalised R̂ < 1.01, and bulk and tail ESS ≥ 1,000 per parameter. The Monte Carlo SE of every reported 5% and 95% quantile is ≤ 0.05 posterior SD. This is what "the run's length justified rather than chosen" means | 16a.9 |
| V7 | Reference calibration | Simulation-based calibration of the M0 sampler at production settings. Each replicate's truth draws **all** free quantities from their priors: M0's targets by `draw_truth(...; names = M0_TARGETS)`, and σ from `p(σ)`, each with a recorded seed. Rank-ECDF within 95% simultaneous bands for every free parameter, at R replications (≥ 50, K1's floor; 100 if the cap allows) | 16a.9 |
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
- **An amplified case for the weight.** In M0 the transition factor varies by only
  about 0.07% between particles, too little for the χ² to see. So the same
  comparison also runs on a test composition built so the rates really differ
  between particles. Its transcription constant reads a pool that the translated
  protein drains strongly, so two particles' constants differ by at least 2× within
  a window. **Mutations:** dropping the transition-probability factor, or dropping
  the 1 s correction alone, each fails this case.
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
- **The count cap.** It is set per dataset, as max(observed count) + 20, and the
  truncation mass is asserted below 1e-12 at every set of constants used.
  Prior-drawn truths can put GAPD's mean near 35.

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
- the SSA's internal state: the pre-drawn next jump time, the aggregated
  propensities and the RNG.

**How the proposal is simulated.** Transcription and decay are suppressed in the SSA
and imposed from the bridge as scheduled events. Translation and translocation run
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
crosses the boundary in both directions and runs all three blocks.
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
- **16a is capped at 5k CPU-h.** Its last task, 16a.11, does four things:
  - measures the cost of one M0 posterior, and of one sweep on a full-scale cell,
    per variant;
  - sets K1's bound, with T3 recording that it was set after the measurement;
  - proposes a costed 16b allocation;
  - **reopens the cap with you.**
- **16b does not start without your approval of the allocation.**
- **The priority order for the allocation.** ⚠️ DRAFT, and yours to reorder:
  1. 16b.1, the tight-prior control. This is K4's stop gate.
  2. 16b.2, the six-parameter recovery.
  3. F6's two load-bearing cells: ptsG under metabolites only, and ENO under
     transcripts only.
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
  `draw_truth(...; names = <the rung's targets>, purpose = :coverage)`, plus σ from
  `p(σ)`, each seed recorded. Each dataset follows 15.7's convention: truth, seeds,
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

## 9. Open questions

- [NEEDS CLARIFICATION: is a replay cheaper than a simulation?] 13.7's 5.11 ms per
  handshake includes the SSA. Replay drops it but keeps the stiff solve. 16a.1
  measures this, and it moves every row of D16.6's projection.
- [NEEDS CLARIFICATION: how degenerate are the particle weights?] GTP's 200-cell
  resolution of 4.5% (14c.5) implies a per-cell coefficient of variation near 64%, so
  metabolite panels discriminate paths strongly. Whether 50 particles suffice is
  16a.7's measurement.
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

**PR:** _not started_

- [ ] 16a.1 Drive the published hybrid from a recorded jump path, and snapshot and
  restore the driver's full state (D16.3's list). Record the path during generation.
  — verify by V2 on 3 dataset cells (bitwise at every handshake), by V2b and its
  mutation, and by the replay's wall-clock per cycle reported against 13.7's 32.2 s.
- [ ] 16a.2 Implement `log p(X | θ)` from the replayed rate constants, on the jump
  clock's [60m − 1, 60m + 59) intervals — verify by V1, including both hand checks
  and all three mutations.
- [ ] 16a.3 Implement the per-gene transcript bridges (60 s and 1 s, by FFBS with
  uniformisation) and their transition probabilities, with the per-dataset cap —
  verify by:
  - sampled bridges always hitting both endpoints;
  - their event-count distribution matching the exact marginal (χ², p > 0.01, 10⁵
    draws per case, multiplicity stated);
  - the transition probabilities matching a matrix exponential to 1e-12;
  - the truncation mass asserted below 1e-12.
- [ ] 16a.4 Implement block 2 from sufficient statistics, with the ceiling asserted per
  proposal — verify by V3.
- [ ] 16a.5 Implement block 1: gradient-free, forward-only prior, reverse constants
  derived, σ updated separately — verify by V4, and by a test that the log-prior
  contains no `kcatR` term.
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
  - per-window update rates and cost per sweep, on M0 and on one full-scale cell, for
    each variant at N ∈ {5, 10, 20, 50};
  - the choice recorded in §4 against 13.4, D13 and 14c.4's gradient report, with the
    two installed traps addressed as D16.3 states.
- [ ] 16a.8 (parent 16.1, build) Build M0 and generate its datasets, with σ drawn per
  replicate. `generate_dataset` gains a `names` pass-through to `check_truth`, which
  defaults to `D11_TARGETS` and would refuse an M0 truth — verify by:
  - the build passing completeness mode;
  - `reduction_report` naming the 15 held enzymes;
  - each dataset carrying its truth (σ included), seeds and report;
  - two runs at one seed being identical.
- [ ] 16a.9 (parent 16.1, run) Run the M0 reference long — verify by V6, which
  justifies the length, and by V7 at production settings on R ≥ 50 datasets, with R
  stated.
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
  - every result's rung and budget-bound label.
- [ ] 16b.8 (parent 16.11) Fill T3 — verify by:
  - all seven criteria carrying a threshold, a measured value and a verdict,
    **published whether or not everything passed**;
  - K1's verdict written here, with its threshold marked as set after its
    measurement;
  - K5 recorded as fired in 14b and answered by 14c and D16.1;
  - every verdict the cap forced labelled budget-bound.

## 12. Amendment log

_No amendments yet._
