# Inference paths after phase 16: exact, modular, amortised

**Date:** 2026-10-13
**Status:** open. No decision is recorded here. Return to it before 16a.11 proposes
16b's allocation (spec §12 2026-10-13).
**Builds on:** [`modular-bayesian-inference-heterogeneous-modules.md`](modular-bayesian-inference-heterogeneous-modules.md),
especially its §5.2 (amortise at the interface) and §7 (Core A′ has no shared
parameter).

## 1. Why this note

One exact posterior on the reduced syn3A model is projected to cost about 60k to
150k CPU-h. The phase 16 cap is 50k, and V7 alone takes about 36k of it. Coverage
needs at least 50 posteriors. So the current route cannot deliver phase 16's claim
at full scale. It was a question of which path to take next, and it is written
down so it is not re-derived.

## 2. The numbers, and where they come from

### One forward simulation is cheap

| Code | Model | 105 min cycle | Source |
|---|---|---|---|
| Thornburg et al.'s published code | Full syn3A | about 2.9 h (98.2 s per simulated min) | syn3a-ensemble `results/speed/speed_summary.csv` |
| syn3a-ensemble port | Full syn3A | about 45 s alone on one core (0.425 s per simulated min); about 47 s per cell on a full node run as one worker process per core | same file; 27,473 twenty-minute cells per node-hour (its README, lane 6a, F-B3) |
| InferCell | Reduced syn3A (32 ODE states, 17 genes) | 32.2 s (5.11 ms per handshake) | spec 13.7, job 17298242 |
| InferCell | M0 | 1.835 s | `dev/scripts/profile_b3_16a9b_result.md` |

"27k cells per node-hour" counts **20-minute** cells (19 simulated minutes). That is
about 4,900 full-length cycles per node-hour. It is one Distributed worker per core;
one 64-thread process gets 15,316 (about 2,700 full-length cycles), and the gap is
garbage collection (below, "Shared"). Updated 2026-10-13 from lane 6a; the first
version of this note quoted the threaded figure.

### A posterior needs millions of them

The exact sampler reruns every cell's trajectory, with every particle, at every
sweep (D16.6):

| Factor | Value |
|---|---|
| Cells | 200 |
| Cycles per cell per sweep | 16.5 measured on M0 at N = 20 (16a.7c); up to about 35 with block 1's replays |
| Sweeps | about 2,000 (D16.6's assumption; M0 mixes at 0.11 to 0.18 effective draws per sweep) |
| **Cycles per posterior** | **about 7 to 14 million** |

At 32.2 s each, that is about 60k to 125k CPU-h, and D16.6's table gives up to about
150k. **These are projections. 16a.11 measures the real figure.**

**The cost is the number of simulations, not the cost of one.** Most of the
millions are proposals the particle filter discards. Faster simulation scales the
bill down by a constant. It does not change the count.

### Coverage multiplies it

Calibration is a claim about repetition: a 90% interval covers truth 90% of the
time. Each replicate draws its own truth and dataset and needs its own posterior.
The covered count is binomial. At R = 50 and 90% nominal, a calibrated sampler
covers about 41 to 49 times, so R = 50 detects only gross miscalibration. 16b.3
asks for R ≥ 50; K4 asks for 100 on the tight-prior control. With the exact
sampler, coverage costs R posteriors. V7 is this check on M0 (R = 50, about 36k
CPU-h projected).

## 3. Speed: what syn3a-ensemble shows

**Per simulated second, the port runs the full cell (7.1 ms) at about 1.4× the
cost of InferCell's 8% reduction (5.11 ms).** The reduction bought almost no speed.
The cost is in how InferCell steps its ODE, not in how much biology it simulates.
The two have not been profiled like for like, so treat the ratio as a pointer.

InferCell's block 3 spends 97% of its time in the ODE step
(`profile_b3_16a9b_result.md`). syn3a-ensemble measured the same choices on the
full ODE (its `code/ode/README.md`, task 2b.2) and moved away from them:

| | InferCell (`src/handshake.jl:1379`) | syn3a-ensemble | Their measurement |
|---|---|---|---|
| Solver | Rodas5P | FBDF | Rodas5P 9.39 against FBDF 2.96 ms per simulated s, 3.2×, at matched tolerance, restarted every 1 s (they rebuild the problem at each hook; we keep one integrator and mark its state modified, which may treat a multistep method differently) |
| Jacobian | Dense ForwardDiff each step (39% of our time) | Sparse AD on the flux Jacobian, `J = S · ∂V/∂u`, column-coloured | 8 colours against 85; `jac!` about 5× faster |
| Linear algebra | Out-of-place `SVector`; 32×32 `lu` falls back to LAPACK and allocates (≥ 13%, plus GC) | In place; KLU with a `jac_prototype` including the diagonal | KLU over dense LU, 12% at 146 states |

**Guess, not measured:** 3 to 5× on block 3, from the solver plus less allocation.
Even then, a full-scale posterior stays above the remaining budget.

**Constraints:**
- Spec task 1.6 chose static out-of-place at 32 states. Revisiting it needs a §12
  entry.
- Our tolerances (atol 1e-10, rtol 1e-8) are not theirs. FBDF must be checked at
  ours, on our smallest pools.
- A new integrator invalidates V2's bitwise replay, and V4 and V5 need rerunning.
  V6 and V7 belong to the current integrator, so the change comes after V7.

**Shared, and diagnosed on their side:** both codes slow down per cell on a full
node. syn3a-ensemble sees 1.85×; InferCell sees about 1.5× contention plus 10 to 14%
GC. syn3a-ensemble's lane 6a (its README, F-B1 to F-B3) found the cause: Julia's
stop-the-world garbage collector. Each of its cells allocates 0.63 GB, 92% of it in
the hook's ODE, and at 64 threads the collector takes 44% of the wall time. Pinning,
per-NUMA copies and GC threads change little; one worker process per core gives
1.79× the cells per node-hour, bitwise identical. Cutting the ODE's allocation is
not done there. For us this points at processes rather than threads for particles,
and at measuring allocation per handshake alongside wall time
(`dev/scripts/bench_ode_solver.jl`). syn3a-ensemble's 6b (GPU batching
with DiffEqGPU) suits our particles too: 20 copies of one small ODE from different
states.

## 4. Path 1: modular Bayesian inference

Fit a module to its own data, with its interface inputs given. Then propagate the
posterior through the rest of the model.

**On Core A′ specifically:**
- **The gene-expression block, given the pool trajectories,** is independent
  per-gene processes with known time-varying rates. The bridges of 16a.3 and block
  2's conjugacy may give it almost exactly and cheaply. Not checked.
- **The metabolic block, given the enzyme trajectories,** is an ordinary IVP
  likelihood. Gradients come back, on the smoothed model of 14c.
- **No hidden-path inference over the coupled system.** That inference is what
  dominates the exact sampler's cost.

**What it loses:**
- **Feedback.** Downstream data cannot inform upstream parameters. The current
  load-bearing cell, ptsG under metabolites only, is exactly what a cut removes
  (§12 2026-10-05 rejected the cut for this reason).
- **Where the interface inputs come from.** Three choices: held at truth
  (synthetic only); measured (real data); or drawn from the other module's
  posterior predictive. The third is the "cut over state" that the modular note's
  §7 calls "a different and unbuilt object".
- **What it calibrates.** The cut posterior is not the Bayesian posterior. Its
  calibration is a legitimate claim, but a different one.

**Why the loss may be small here:**
- D13 means no parameter is split across two fits.
- ENO under transcripts only is a measured null (posterior SD 99.9% of prior).
- FBA is open (curvature upper-2σ about 0.72).
- If the ptsG direction also turns out weak, the cut is nearly free. The M0
  comparison in §7 measures this.

Note that `src/boundary.jl`'s `sequential_infer` and `iterative_infer` match on
parameter names, so on Core A′ they pass nothing. Path 1 needs new code, not those.

## 5. Path 2: amortised simulation-based inference

Draw θ from the prior, simulate a dataset, and train a network on the pairs. Every
later posterior is one network evaluation.

**Cost.** A 200-cell dataset is about 1.8 CPU-h on InferCell today. 10⁴ datasets
is about 18k CPU-h, and 10⁵ is about 180k. So a single posterior is only a few
times cheaper than exact. **The gain is reuse:** coverage at R = 50 or R = 10⁴
costs the same, so a sharp calibration check becomes affordable.

**The worry: many parameters.**
- The simulation count tracks how complicated the posterior is (correlations,
  ridges, multimodality), not the volume of the prior. Simulations are prior
  draws, not a grid.
- In practice, amortised SBI is dependable at tens of parameters, not hundreds.
  A whole-cell model is the latter.
- The data dimension matters as much: 200 cells of time series need summary
  statistics or a learned embedding.
- Sequential SBI spends simulations near the observed data, but it gives up
  amortisation.
- Fixing what the data cannot inform keeps the dimension down, which needs
  identifiability screening first.

**What it gives up:** exactness. The network's calibration must be shown by SBC,
which is cheap here because posteriors are free. Julia's SBI tooling is thin. The
practical route is simulations in Julia and training in Python's `sbi`.

## 6. Combining them

**Modularity is how amortised inference escapes the dimension problem.**
- **Amortise per module.** Train one network per module, conditioned on its
  interface inputs. A thousand-parameter cell becomes, say, fifty problems of
  twenty. This is the modular note's §5.2 recommendation, and it was made before
  the speed numbers existed.
- **Prefer likelihoods to posteriors per module.** A neural likelihood or ratio
  per module composes by multiplication inside an ordinary sampler. Composing
  per-module posteriors divides out priors and is fragile (modular note §3.2,
  §7 item 2).
- **Compose at two strengths:**
  - **Cut composition:** Path 1 with every module fit instant.
  - **Coupled composition:** iterate distributions over interface trajectories
    between modules. `iterative_infer` was the toy form, and it was overconfident
    (`2026-05-22-calibration-curve.md`). Whether a working form converges to the
    joint posterior is open research.
- **Or a hierarchy:** modular cuts pin most parameters, and whole-cell SBI or the
  exact sampler handles the few whose information crosses a boundary.

**Nothing built so far is wasted:**
- **The exact sampler and the M0 reference** become the yardstick every other path
  is judged against.
- **syn3a-ensemble's speed** feeds both paths: training simulations for Path 2,
  and forward propagation for Path 1.

## 7. A first experiment: three posteriors on M0

On V7's 50 M0 datasets, compare:
1. **Exact:** V7's chains.
2. **Modular:** block-wise fits with the interface given.
3. **Amortised:** NPE trained on M0 simulations.

**Measure, for each:**
- the per-marginal gap from exact, beside V8's floor;
- the ptsG and ENO pair;
- SBC under V7's own rule (bands at 1 − 0.05/6).

**Cost (rough):**
- An M0 dataset (50 cells × 1.835 s) is about 90 s.
- 10⁴ training datasets is about 250 CPU-h, and 10⁵ about 2.5k.
- The modular fits are mostly conjugate or ODE-only.
- It is small next to V7.

**What each outcome would mean:**

| Outcome | Reading |
|---|---|
| Modular matches exact and is calibrated | Feedback carries little in this model. Path 1 scales to the full cell |
| Modular is calibrated but misses ptsG's shrinkage | The cost of the cut is measured, and it is the feedback claim. Path 1 needs a coupled form |
| NPE matches exact at affordable simulation counts | Path 2 is viable at this dimension. The open question becomes dimension, and §6 |
| NPE is miscalibrated where exact is not | SBC on M0 caught it, which is what the reference is for |

## 8. Other avenues, not yet weighed

**Make the exact sampler cheaper or mix faster:**
- **Constant factors:** §3's solver work, particle-parallel sweeps (R16.8), and
  GPU batching.
- **Fewer sweeps:**
  - non-centred parameterisations, or joint θ–path moves (the Gibbs chain is slow
    because θ and the path are tightly coupled);
  - correlated pseudo-marginal MCMC, which marginalises the path instead of
    sampling it;
  - SMC².

**Approximate the stochastic block:**
- **LNA, moment closure or the chemical Langevin equation** (Schnoerr, Sanguinetti
  and Grima 2017) give a smooth likelihood and gradient-based inference over the
  whole model.
- **Caution:** transcript means are 0.3 to 2.4 copies (spec §3), where these
  approximations are poor.
- **Differentiable or smoothed SSA** (StochasticAD; the differentiable Gillespie
  algorithm) is another option.

**Exploit independent cells:**
- Transcripts are independent across cells. Compositional score-based SBI
  (Geffner, Papamakarios and Mnih 2023; Linhart et al.) trains on single-cell
  simulations and composes for N cells. Training could then cost single cells
  (about 32 s each), not 200-cell datasets.
- Bulk metabolites are a population mean and break the factorisation. Its
  likelihood term would have to be handled separately.

**Emulators and history matching:**
- Gaussian-process or neural emulators of the expensive model, or of the
  interface map.
- History matching rules out implausible regions of a large parameter space
  cheaply before any posterior is attempted. It is standard for large simulator
  calibration, and it answers the many-parameters worry directly.
- Linked emulators follow the module graph (modular note §3.7).
- The spec's "two surrogate proposals" (Beyond the horizon) belong here.

**Multi-fidelity:** M0 → reduced syn3A → full syn3A as a fidelity ladder, with the
cheap rung correcting the expensive one.

**Semi-modular inference** (Carmona and Nicholls): interpolates between cut and
full, with the interpolation weight itself chosen by a loss. It needs a cut first.

**Experimental design:** which measurements carry the most information across the
boundary. The bulk-identifiability analysis of 16a.7 is a first instance.

## 9. Decisions this leaves open

1. **Does the claim change?** From "calibrated joint inference across the
   boundary" to something like "what inference on a hybrid model can trade for
   scale, measured against an exact reference". That needs a §12 entry, or a new
   spec.
2. **Order:**
   - Does §7's M0 comparison come before 16a.11, alongside it, or after 16b?
   - Does the solver speed-up go first, since it changes every cost above?
3. **The relation to syn3a-ensemble.** The 2026-10-12 decision keeps it standalone
   until its phase 8. Copying the solver and Jacobian choices into InferCell is
   independent of that decision.
