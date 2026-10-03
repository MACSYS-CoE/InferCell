# 16a.7: the variant choice, and the KP stop

**Verdict: KP fires.** No variant reaches the update-rate threshold on M0 at
N = 50, and particle Gibbs never changed a window of the full-scale cell. The
cause is the metabolite likelihood. Proposals from the translation path's prior
land hundreds of nats below the reference within a single window. More
particles and ancestor sampling cannot close a gap that size. The sub-spec's
§8 says to stop and amend, which is where 16a.7 stands.

## The KP measurement (`dev/scripts/csmc_variants_16a7.jl`, runs at `0b23895`; cell 4's ancestor-sampling runs at `2566e79`; merge job 17884470)

The setup:
- **M0:** five cells of replicate 16085 (σ = 0.0727), at the replicate's truth.
  Each recorded path matches the dataset's latent.
- **Full-scale:** cell 30001 at 15.7's truth. The panel is observed afresh, at
  σ = 0.1 from `Xoshiro(1508)`, because 15.7's observations are not tracked.
- **Chains** start from each cell's true path, which is an exact posterior draw at
  the true θ. So the rates need no burn-in.
- **The update rate** is the fraction of sweeps in which a window's events change.
  The threshold is 10% (⚠️ DRAFT).
- **Cost** is wall time per sweep divided by a plain cycle timed in the same job.
  Each task's first sweep, which compiles, is dropped.

| Cell | Variant | N | Sweeps | Cycles per sweep | Windows below 10% | Per-window rate, windows 1 to 10 |
|---|---|---|---|---|---|---|
| M0 | PG | 50 | 40 | 51.2 | 9 of 10 | 0.00 0.00 0.00 0.00 0.05 0.00 0.00 0.03 0.20 0.07 |
| M0 | AS, lag 1 | 50 | 40 | 90.5 | 5 of 10 | 0.00 0.05 0.07 0.12 0.15 0.03 0.00 0.15 0.10 0.12 |
| M0 | AS, lag 5 | 50 | 40 | 198.2 | 6 of 10 | 0.00 0.07 0.03 0.10 0.12 0.03 0.00 0.10 0.10 0.05 |
| M0 | Full PGAS | 50 | 40 | 228.6 | 9 of 10 | 0.05 0.03 0.05 0.03 0.05 0.00 0.03 0.05 0.17 0.05 |
| Full-scale | PG | 50 | 8 | 49.1 | 105 of 105 | none of 840 window-sweeps changed |

- **Every variant has a window that never changes in 40 sweeps.** The upper 95%
  bound on such a window's rate is 7.2%.
- **The costs match D16.6's per-variant model.** PG costs about N cycles. Full
  PGAS on M0's 10 windows costs about 4.6N, against the projected 5.5N.
- **The other N** (5, 10, 20) were not run. Mixing at N = 50, the largest, decides
  KP, and the smoke runs at N = 5 changed no window either (jobs 17877524 to
  17877529).

### A bug the measurement found: an impossible firing (fixed at `2566e79`)

Ancestor sampling replays the reference's future from each candidate particle's
state. A recorded firing can have zero propensity there. The replay applied it
anyway, which left a count negative, the next propensity read below zero, and
`log` threw. That killed M0 cell 4 under all three ancestor-sampling variants
(job 17878618 task 4 and its siblings).
- **The fix:** `replay_step!` now logs −Inf for a firing with zero propensity, and
  `path_logdensity` returns −Inf once one occurs. That candidate gets no ancestor
  weight.
- **Tested** in `test/test_path_density.jl`.
- **The 16a tests pass at the fix,** 160 of 160 (job 17882299).
- Cell 4 was rerun under the fix (jobs 17882309 to 17882311).

## Why: the metabolite likelihood (`dev/scripts/csmc_collapse_16a7.jl`, jobs 17879029 and 17879030 at `7c5c155`)

200 proposals per window start from the true path's own state, so they differ
from the reference only inside the window.

**The transcript part of the weight is the same for every particle.** The
metabolite part separates them:

| Cell | Metabolite log-likelihood, proposal minus reference: mean (sd), windows 1, 2, 3 |
|---|---|
| M0 | −573 (1,360); −1,648 (2,510); −1,867 (2,489) |
| Full-scale | −40 (55); −160 (118); −401 (255) |

**Three pools carry the gap: `M_pi_c`, `M_3pg_c` and `M_2pg_c`.** All three sit at
lower glycolysis and drain toward empty:
- on M0, from 817, 172 and 34 particles at 60 s down to 485, 91 and 18 at 180 s;
- in M0's window 2, the 10th-percentile proposal leaves 3PG at about 2 particles,
  against the reference's 227;
- the proposals' spread on these pools, in σ units of log, is 16 to 35 on M0 and 4
  to 8 at full scale;
- an observation of 18 particles at σ = 0.07 pins the pool to about ±1.3 particles.

**The pools' response to the path is nonlinear,** shown by how much of the three
pools' deviations each set of regressors explains, over windows 1 to 3:

| Regressors | M0 | Full-scale |
|---|---|---|
| End-of-window protein counts, every gene (R²) | ≤ 0.16 | ≤ 0.23 |
| Translation exposure, protein-seconds (adjusted R²) | ≤ 0.07 | ≤ 0.07 |
| Every reaction's count and exposure (adjusted R²; M0 14 regressors, full-scale 98 to 104) | 0.50 to 0.59 | 0.25 to 0.37 |

- **The GAPD hypothesis is refuted.** GAPD's single correlation is 0.15 to 0.23.
- **The implication:** a proposal guided by a linear tilt toward the observed
  pools would close only part of the gap.

## Is the proposal right? (`dev/scripts/csmc_rank_16a7.jl`, jobs 17880205 and 17880206 at `a34ee7d`)

At every window, the true path's 17 panel pools are ranked among proposals from
its own state. The proposals are weighted by the transcript part of the weight
alone, so the ranks are uniform when the proposal and weight are right.

- **M0** (5 cells × 10 windows, 100 proposals each): uniform. 16 of 17 pools have
  KS p > 0.14, and `M_lac__L_e` has p = 0.03.
- **Full-scale** (105 windows, 50 proposals each): `M_pi_c` is uniform (p = 0.33).
  Seven pools sit at p = 0.02 to 0.03:
  - `M_lac__L_c`, `M_g3p_c`, `M_dhap_c`, `M_3pg_c`, `M_2pg_c`, `M_fdp_c` and `M_gmp_c`;
  - mean ranks 0.43 to 0.57;
  - the 17 tests share one trajectory and are correlated, so the asymptotic p is
    optimistic.
- **Watch item, not a pass:** possible proposal bias at full scale. V5 has only
  ever been run on M0 cut to ptsG and the toy. Several independent full-scale
  cells would decide it. It does not change KP: the collapse is as severe on M0,
  where the proposal checks out.

## Compute

About 15.3 CPU-h in all (`phase16_ledger.md`), which puts 16a at about 145 of its
5k.

# After the amendment (§12 2026-10-02)

## 16a.7b: the window proposal at full scale passes (jobs 17904315 and 17904316 at `8c66944`)

The rank check, on five independent cells (30001 to 30005) at 15.7's truth, over
525 windows with 50 proposals each:
- **Every pool's ranks are uniform.** The smallest unadjusted KS p is 0.057
  (`M_2pg_c`), and every Bonferroni p is 0.97 or more.
- **Cell 30001's earlier p = 0.02 to 0.03 was chance.**

## 16a.7a: annealed windows pass V5, and miss the threshold

**The code** is `annealed_window`, with `csmc_sweep(...; schedule)`, at `2e5f7e9`.
The 16a tests pass, 199 of 199 (job 17904432). They include:
- the backward chain keeping the reference's path and state;
- a head move keeping the path before its cut;
- the AIS weight matching plain importance sampling within 3 SE on the toy.

**V5 at K = 5 passes on both cases.** It reuses the rejection draws of `644c140`'s
run.

| Case | Jobs | Off: χ² p | On: z |
|---|---|---|---|
| M0 | 17905037 to 17905039 | ≥ 0.086 | within ±1.97 |
| Toy | 17905042 to 17905044 | ≥ 0.166 | within ±1.49 |

**The update rates miss the threshold everywhere** (jobs 17905020 to 17905026,
M0's five cells, 40 sweeps each, geometric schedules from β1 = 1e-4).

| K | N | Cycles per sweep | Windows below 10% | Highest window rate |
|---|---|---|---|---|
| 20 | 2 | 42.1 | 10 of 10 | 0.00 |
| 20 | 5 | 106.3 | 10 of 10 | 0.03 |
| 50 | 2 | 101.9 | 10 of 10 | 0.00 |
| 50 | 5 | 261.8 | 10 of 10 | 0.05 |
| 100 | 2 | 202.6 | 10 of 10 | 0.03 |
| 100 | 5 | 519.6 | 10 of 10 | 0.05 |

### Inside one annealed window (job 17905489 at `e7560ab`)

The setup: M0 cell 1, windows 1 to 3, from the true path's state, with 5 forward
chains per K and the reference's backward chain. "Replacement" is the chance that
a fresh particle displaces the reference at N = 5.

| Window | K | Acceptance, last fifth | Final log w − reference's, median | AIS log weight − reference's, median | Replacement at N = 5 |
|---|---|---|---|---|---|
| 1 | 50 | 0.14 | −15.5 | −16.2 | 0.001 |
| 1 | 200 | 0.10 | −19.2 | −18.1 | 0.000 |
| 1 | 1,000 | 0.08 | −1.6 | −8.4 | 0.003 |
| 2 | 50 | 0.10 | −15.7 | −21.8 | 0.117 |
| 2 | 200 | 0.03 | −3.7 | −1.6 | 0.601 |
| 2 | 1,000 | 0.07 | −3.0 | −2.0 | 0.771 |
| 3 | 50 | 0.18 | −25.7 | −38.7 | 0.000 |
| 3 | 200 | 0.11 | −0.6 | −8.0 | 0.496 |
| 3 | 1,000 | 0.10 | −0.4 | −3.1 | 0.345 |

- **Annealing closes the gap.** Unannealed proposals sat hundreds of nats below the
  reference. With K ≥ 200, the chains end within a few nats of it.
- **Windows 2 and 3 become exchangeable from K = 200,** at a replacement chance of
  0.35 to 0.77.
- **Window 1 does not, even at K = 1,000.** The chains reach the right log weight,
  but the AIS weight carries an 8-nat penalty built up along the schedule.
- **The late stages mix slowly,** at 3% to 14% acceptance.
- **Cost.** K = 200 at N = 5 is about 1,000 M0 cycles per sweep. That is about 42k
  CPU-h for one M0 posterior of 50 cells × 2,000 sweeps, and K = 1,000 about five
  times that. Both are above 16a's cap and near or above the phase's.

**KP fires again on M0 within the cap.** This needs your decision; see the handoff.

## Option 1: the observation floor at a cited detection limit (jobs 17933661 and 17933662 at `10cf581`, 2026-10-03)

**The cited limit.** The most sensitive single-cell metabolite detection limit
found is 0.2 amol of NAD⁺, by capillary electrophoresis with enzymatic cycling
(Lin, Trouillon, Safina and Ewing, "Chemical Analysis of Single Cells", *Anal. Chem.*
2011). That is about 1.2 × 10⁵ molecules. Single-cell mass spectrometry reports
limits of order an attomole (6 × 10⁵ molecules) and above.

**The panel against it.** In a syn3A cell the 17 panel pools hold from about 4
(`M_lac__L_e`) to about 400,000 (`M_fdp_c`) particles. Only `M_fdp_c` is above
1.2 × 10⁵; the next largest, `M_atp_c`, holds about 39,000.

**The run.** The collapse attribution reran with the panel re-observed at a floor
of 1.2 × 10⁵ particles, keeping each observation's own noise draw
(`CSMC_FLOOR=120000`).
- **The metabolite log-likelihood gap vanishes.** Proposal minus reference is
  0.00 ± 0.02 nats in every window, on M0 and at full scale, against 40 to 1,870
  nats at the one-particle floor.
- **The pools still move as before,** for example 3PG and 2PG on M0. The panel can
  no longer see them.

**What it means.**
- **At a cited detection limit, the per-cell panel carries almost no information
  about a syn3A cell's path.** The sampler would then mix, as V5's metabolite-off
  half does, but because the data say nothing, not because the sampler improved.
- **The panel at a one-particle floor and σ = 0.1 is an idealisation that no
  existing single-cell assay achieves.** It is the idealisation that makes exact
  path inference infeasible.
- **This reaches past 16a.7.** 15.6's panel choice and 15.8's rank-6 identifiability
  were measured at the one-particle floor. At this floor, ENO's and FBA's per-cell
  identifiability would rest on `M_fdp_c` alone. That has not been measured.
