# 16a.3 and 16a.4: transcript bridges, and block 2 (V3)

Run of record: job 17758553 at `d9b4aa5` (`dev/scripts/phase16a_tests.slurm`).
16a.3 passed 27 of 27, and 16a.4 passed 16 of 16. The same job reran 16a.1 (16/16)
and 16a.2 (12/12).

## 16a.3: transcript bridges (`test/test_transcript_bridge.jl`)

Two rate cases:
- *typical*, the 15.7 truth's scale: k = 0.02/s, μ = 0.009/s, cap 22;
- *amplified*, V5's kind: k = 1/s, μ = 0.2/s, cap 30.

Results:
- **Transition probabilities** match a matrix exponential to 1e-12, and every row
  sums to 1 within 1e-12. That holds for both cases at τ = 1, 20 and 60 s. A
  login-node check put the gap at 6.7e-16 (typical, 60 s) and 2.6e-15 (amplified,
  20 s).
- **Bridges hit both endpoints** in all 4 × 10⁵ draws, with ordered times inside
  [0, τ) and ±1 steps.
- **Birth counts are exact.** Each case's histogram of births matches the exact
  distribution from a matrix exponential of the chain augmented with a birth
  counter, at χ² p > 0.01. The four cases were (typical, 60 s, 1→1),
  (typical, 60 s, 0→2), (amplified, 1 s, 3→4) and (amplified, 20 s, 5→12). The
  family-wise false-failure rate is about 4%.
- **The cap rule is a floor, not a guarantee.**
  - At typical rates, max(observed) + 20 drops 2.4e-20.
  - At amplified rates over 20 s from a count of 10, cap 30 drops 7.6e-12, above
    D16.3's 1e-12. The first run failed there (job 17758372, at `5d4bd8c`).
  - `adequate_cap` grows the cap in steps of 5 until the bound holds. The test
    asserts that the rule fails there, that the grown cap passes, and that at
    typical rates the rule is already enough.
- **An earlier run failed on a name clash** (job 17758130, at `aa5ef32`).
  `test/jump_test_models.jl` defines a `BirthDeath` double in Main, which shadowed
  the export. The type is now `TranscriptChain`.

## 16a.4: block 2 (`test/test_block2.jl`), V3

The setup:
- One gene, simulated by Gillespie at k = 0.05/s and μ = 0.02/s, starting at 2 and
  observed exactly every 60 s for an hour (60 windows).
- Priors: LogNormal(ln 0.04, ln 2) on k, and LogNormal(ln 0.015, ln 2) on μ.

The two sides:
- **Exact posterior.** A likelihood from the transition probabilities between
  consecutive counts, on a coarse grid over the prior, then a fine 451 × 451 grid
  from −7 to +12 posterior SD. The upper tails are heavy, because k and μ rise
  together along the ridge that fixes the mean count. Marginal CDFs are by
  trapezoid, with quantiles by linear interpolation, and the grid edges hold under
  1e-8 of the peak.
- **Gibbs.** An exact bridge per window, then block 2's `rate_update` on each rate
  from the path's birth count and time (k), or death count and occupancy (μ). That
  is 20,000 sweeps after 1,000 of burn-in.

Result:
- **At each exact 5, 25, 50, 75 and 95% quantile, on both rates,** the fraction of
  draws below it is within 3 batch-means SE (50 batches) of the nominal level.
- **Both posteriors are narrower than 0.9 of the prior's width.**
- **The Jacobian is carried.** With no data, the conditional equals
  `Normal(u; μ, s)` and `logpdf(LogNormal, θ) + ln θ`. An accepted value past the
  bound throws.

V3 freezes the rates, so it cannot see a weight error that varies between
particles. V5 (16a.7) covers that.
