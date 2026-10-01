# 16a.2: the path log-density (V1)

Run of record: job 17756874 at `262e662` (`dev/scripts/phase16a_tests.slurm`),
`test/test_path_density.jl`, 12 of 12 pass.

- **Hand check, one gene with frozen rates.** `ToyExpression` over 5 s with five
  firings (two transcriptions, two translations, one decay). The log-density
  matches the hand sum to 1e-12, and so do the counts and transcription's exposure.
- **Hand check, across a rebuild.** `ToyRebuiltExpression` over 120 s, with the ATP
  pool dropped to 0.1 mM before the rebuild at handshake 60. That moves the
  transcription constant by more than 2×. There are firings at 58.5 s and 59.5 s,
  on both sides of 60 − 1 s. The log-density matches a hand sum using the old
  constant on [0, 59), the new one on [59, 119) and the next on [119, 120), to 1e-12
  relative.
- **Mutations.** Three alternative values each differ from the hand sum by more than
  1e-3: no exposure term, the previous interval's constant over [59, 60), and a
  rebuild at 60 s. So an implementation with any of those errors fails the check.
- **Importance identity at suite size.** Over 400 toy paths of 180 s each, the mean
  of p(X|θ1)/p(X|θ0) is within 3 SE of 1 for a decay shift (×1.1) and for an
  ODE-only shift (the pool's kcat ×2).
  - **The ODE-only shift barely moves the toy's rates.** The toy's rebuild law has
    elasticity about 0.05 at its pool, so this case is weak.
  - **The 10⁴-path identity on M0 that V1 specifies runs when M0 exists (16a.8).**

**A first run failed the identity, from a toy inconsistency** (job 17756205, at
`31c437f`).
- `set_parameters!` recomputes a module's rebuilt constants from the initial pools
  whenever any of its parameters changes.
- The toy declared `k_tx = 2.0` against its rebuild law's 5.71 at the initial ATP.
  So a decay shift also moved transcription over [0, 59 s), and θ0 and θ1 were two
  different models: every weight was about 0.
- The fix declares the toy's constant at its law. Core A′'s constructors compute
  their initial constants from the same law as the rebuild, so it is not exposed.

## The 10⁴-path identity on M0 (job 17801682, at `e1841e7`)

`dev/scripts/path_density_16a2_m0.jl`, 100 tasks, merged by
`path_density_16a2_m0_merge.jl`.
- θ0 is a truth over `M0_TARGETS` drawn with `Xoshiro(1609)`.
- Each path is simulated at θ0 and recorded, then replayed with its density at θ0
  and at each shift.
- Paths have a median of 594 firings (167 to 1,004) over 600 s.

| Shift | Mean of p(X\|θ1)/p(X\|θ0) | SE | \|mean − 1\| / SE | Weight SD | Largest share | ESS |
|---|---|---|---|---|---|---|
| CME: GAPD's promoter ×1.1 | 1.00297 | 0.00243 | 1.22 | 0.2430 | 2.9e-4 | 9,445 |
| ODE: ENO's forward and reverse ×1.2 | 0.99981 | 0.00017 | 1.07 | 0.0175 | 1.6e-4 | 9,997 |

**Both are within 3 SE, so V1's identity holds on M0.** The ODE-only shift moves
the density through the pools alone, with a weight SD of 1.75%. That is the case
the toy could barely test. A two-path smoke run preceded it (job 17801198).
