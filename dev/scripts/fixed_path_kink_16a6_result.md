# 16a.6 (parent 16.3): a clipped drain under a fixed path still kinks

## The test of record (job 17801857 at `b7da3f0`; `test/test_fixed_path_kink.jl`, 9 of 9)

The setup:
- The phase 3 toy with a small ATP pool (400 particles) and a costly translation
  (40 per event), so the drain clips.
- A path recorded at kcat = 30 under the published model.
- Check 7's census finds the clip: the first drain to carry a deficit is at
  handshake 3.
- With that path held, kcat is varied from t = 0, and the ATP state after the
  drain at handshake 3 is read.
- The clip's kcat, θc = 17.275, is bisected on the clamp over continuous pools.

The slope jump across the clip, as a fraction of the paying side's sensitivity:

| Spacing either side (steps) | Clamp, continuous pools | Smoothed drain, continuous pools (14c.5's sampler model) |
|---|---|---|
| 30 | 1.265 | 0.236 |
| 3 | 1.267 | 0.024 |

- **The clamp's jump does not shrink with the spacing, so conditioning on the
  path leaves the kink.** It reads above 1 because the reference sensitivity is
  taken at θc/1.5, where it is smaller.
- The smoothed drain's jump shrinks tenfold with the spacing, as 14c.5's
  continuity test found without the path held.
- The published model, under fractional carry, has a zero finite-difference slope
  at 15 or more of 21 points across the clip.

This is the reason D16.1 takes block 1 off the gradient: on the published model
the derivative is either zero or kinked, and holding the path changes neither.

## On Core A′: open

`dev/scripts/fixed_path_kink_16a6.jl` repeats the test on 15.7's cells at its
truth (Xoshiro(1507)), varying one forward constant with its reverse derived. It
is recorded as open, by decision on 2026-10-01.
- **Cells 30001 and 30002 never clip** within the cycle (jobs 17802079, 17802347).
  The first to clip is cell 30003: `GTP_translat` on `M_gtp_c` at handshake
  3,409.
- **ENO and FBA cannot move that clip.** No value within 2^±20 of the truth flips
  it (jobs 17802347, 17814374). Block 1's two constants do not reach the GTP pool
  strongly enough.
- **PGK3, which makes GTP from GDP, flips it at 1.357× its nominal 140.8**
  (job 17814375). But the paying side's finite-difference slope is exactly zero
  both 1.5× and 2% past the clip (jobs 17814375, 17830703), and that is not
  explained. It may be the script or the pool. Diagnosing it would need the pool
  before the debit, the accrual and the state after, around the clip.
