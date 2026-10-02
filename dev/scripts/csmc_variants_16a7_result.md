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
