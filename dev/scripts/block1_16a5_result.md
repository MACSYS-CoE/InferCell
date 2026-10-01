# 16a.5: block 1, and V4

## Unit tests (job 17801360 at `e1841e7`; `test/test_block1.jl`, 18 of 18)

- **Reverse constants.** `derived_ode_values` reproduces the truth's derived
  reverse constants to 1e-12 and keeps every equilibrium constant.
- **The prior.** It is the forward constants' `logpdf(LogNormal, θ) + ln θ` to
  1e-12. A reverse constant passed as a forward is refused.
- **The likelihood.** With the data equal to the replay, the residual is exactly
  zero, and the path term equals the replay's own density at the full truth to
  1e-12.
- **The σ step** matches its exact 1-D conditional at five quantiles, within
  3 batch-means SE over 40,000 draws.
- **A sweep** runs and can hold σ.

## V4 on M0 (`dev/scripts/block1_16a5.jl`)

The setup:
- Replicate 16085: three M0 cells at its truth, with σ = 0.0727 held, and each
  path re-recorded and checked against the dataset's latent.
- 20 chains of 50 sweeps, each started from a grid draw, with the first 5 sweeps
  dropped.

**The exact posterior is a 2-D grid that shares no code with block 1's target.**
- **Agreement.** Block 1's log target equals the grid's at 15 points across the
  posterior, to a spread of 7.4e-13 (job 17830938).
- **ENO is narrower than the coarse grid could see.** Posterior sd is 0.019 in
  ln ENO, against a prior of 0.157. The first fine box truncated it (0.39 of the
  peak at its edge). The fix was stage 3 at ×3 width, then stage 4, extended at
  the same step to ±8 SD, with edge mass 7.0e-13.
- **FBA needed its own grid.** The main grid's FBA axis kept the coarse stage's
  2.4 points per SD, which biased FBA's quantiles by up to 0.010 in CDF units. The
  FBA-fine grid has 161 points over ±8 SD against 41 ENO points, with edge mass
  1.2e-9 (jobs 17839207 and 17839208).

**The first final test failed** (job 17830702). ln FBA's 75% point was at 3.02 SE
against the coarse FBA quantiles.
- The diagnosis, at your request, found mixing sound: lag-1 autocorrelation was
  −0.03 for ENO and 0.07 for FBA, giving about one effective draw per sweep, and
  the spread of chain means matched the SE.
- What it found instead was the grid's FBA resolution above.

**The run of record passes** (job 17839209, at `18f2c69`). All ten points are
within 3 SE:

| | p = 0.05 | 0.25 | 0.50 | 0.75 | 0.95 |
|---|---|---|---|---|---|
| ln ENO, chain fraction below the exact quantile | 0.0444 ± 0.0068 | 0.2578 ± 0.0172 | 0.5022 ± 0.0191 | 0.7500 ± 0.0154 | 0.9456 ± 0.0073 |
| ln FBA, chain fraction below the exact quantile | 0.0456 ± 0.0081 | 0.2800 ± 0.0130 | 0.5189 ± 0.0161 | 0.7811 ± 0.0137 | 0.9633 ± 0.0074 |

**Watch item.** FBA's fractions sit above their levels at four of five points: by
+2.3, +1.2, +2.3 and +1.8 SE. The chains place FBA slightly lower than the grid.
This passes, but it is not settled. With 900 draws it cannot tell chance from a
small bias. The M0 reference run (16a.9) samples FBA with far more draws, and its
V6 and V7 are where to look next.
