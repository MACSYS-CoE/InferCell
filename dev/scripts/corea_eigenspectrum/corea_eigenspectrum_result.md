# KIM-181: the Hessian's eigen-spectrum with the ten glycolytic forward constants free

Linear: [KIM-181](https://linear.app/kimpson/issue/KIM-181/hessian-eigen-spectrum-on-the-reduced-model-with-a-wider-metabolic).
Feeds Fig. 3b and §3.3.

**The result.** The reduced model was refitted with all ten glycolytic forward constants free,
not only ENO and FBA. At the new MAP the data reach:

- **24 of the 29 directions** beyond the prior. A direction counts when the whitened −H
  eigenvalue exceeds 2, which is where the data's curvature exceeds the prior's.
- **5 of the 10 directions in the metabolic block.** Its marginal precision is whitened by the
  prior. With only ENO and FBA free, KIM-178's block had 2 of 2.

Widening the block adds three directions above the prior: PFK plus PYK (twice) and GAPD. It
adds five that are not: PGI, TPI, LDH_L, PGM and PGK. Nothing else moves. ENO and FBA's
eigenvalues barely change, from 11.8 and 39.9 to 12.1 and 40.2, and ranks 1 to 19 of the whitened
spectrum match KIM-178's to about 0.1%. The 21 coordinates shared with
KIM-178 shift by at most 0.18 Laplace SD.

**The count depends on the prior table.** Five of the new constants have a prior SD of 0.05
(PGI, TPI, PGK, PYK, LDH_L), against 0.69 for the promoters. A second count does not depend on
that table. It measures the likelihood's curvature alone, −H − P, at a common width of
geometric SD 2: a direction counts when the data pin it tighter than a prior of that width
would. That count is **28 of 29**, and **9 of 10** for the metabolic block.

So the data constrain almost every glycolytic constant on an absolute scale, PGI, TPI, PGM and
PGK among them. They do not constrain them beyond the tight priors the reduced model already
gives them. Each reaction's own marginal precision under the likelihood alone, 1/(L⁻¹)_ii
with L = −H − P and no prior on any coordinate, × (ln 2)², ranks them:

| PYK | GAPD | PFK | ENO | PGK | FBA | PGM | TPI | PGI | LDH_L |
|---|---|---|---|---|---|---|---|---|---|
| 550 | 200 | 175 | 153 | 133 | 60.3 | 49.4 | 20.5 | 3.10 | −0.013 |

Only LDH_L is flat in both readings. PGI's value of 3.10 means the data alone give it about a
geometric SD of 1.5; its whitened eigenvalue, 1.016, means that is no tighter than its prior. The issue's "done when" count is the first,
24 and 5. Which count Fig. 3b shows is Tom's call. The second is local curvature at this MAP,
under these priors, not a refit at geometric SD 2. Both counts were in the script at `49f345c`
before the run.

## The run

- Job 21842933, after the precompile job 21842932, on dave308. Commit `49f345c`, 1 CPU, 2 h 6
  min of wall time.
- The script is `corea_eigenspectrum.jl`.
  - It starts from KIM-92's dataset (seed 17200001) and the fitted likelihood, with δ and S
    from KIM-92's prep.
  - The model is KIM-92's, with all 20 glycolytic constants free: 10 forward and 10 reverse.
    The reverses are derived through the Haldane relation, so the 8 newly freed forward
    constants add 8 coordinates, for 29 in all.
  - Each new constant's truth is its nominal value, which is also its prior median and its
    warm start. So the run measures curvature where the truth sits at the prior's centre. It
    does not test recovery of a truth away from it.
  - The inputs are KIM-92's archive (`infercell-kim92-archive-2026-11-17`, `dataset.jls` and
    `prep.jls`) and KIM-178's MAP (`infercell-kim178-2026-11-27/kim92`). That directory is
    identical to `infercell-kim178-archive-2026-11-28/kim92`, which is checksummed.
- The MAP is warm-started from KIM-178's MAP on the same dataset, with the 8 new constants at
  their prior medians.
- The outputs are archived in `/fred/oz022/tkimpson/infercell-kim181-archive-2026-11-28/`, with
  `SHA256SUMS`. Two smoke runs (jobs 21839624 and 21840130) checked the pipeline and are not
  results.

## The checks

| check | value | criterion | |
|---|---|---|---|
| the wider model scores the same data: density at the warm start, less the narrow model's at KIM-178's MAP and the 8 new priors | Δ 1.9e-9 | < 1e-6 | holds |
| the MAP converged | log density +1.56 above the warm start; \|∇\|·sd 0.014; Newton decrement 1.9e-7 | `reference_theta` returns `:map` | holds |
| −H symmetric before symmetrising | 9.91e-9 | ≤ 1e-8 (the issue's) | holds, by 0.9% |
| H[S_0445, S_0445] against second differences of the density | −820.446 vs −820.529 | rel < 1e-2 | holds (1.0e-4) |
| H[kcatF_PYK, kcatF_PYK] | −2184.72 vs −2184.73 | rel < 1e-2 | holds (2.7e-6) |
| H[kcatF_PFK, kcatF_ENO], a new cross term | −0.953911 vs −0.953861 | rel < 1e-2 | holds (5.3e-5) |
| −H positive definite | yes; nearest the ±6 SD box at 2.65 prior SD | | |

The symmetry check passes with little margin. The three second-difference checks are the
evidence for −H's accuracy.

## The metabolic block

The marginal precision of the 10 forward constants, whitened by their prior. Each
eigenvector is close to one reaction or to one pair, so the loadings name them. The
squared loadings shown are those ≥ 0.10.

| eigenvalue | loads on | likelihood alone, × (ln 2)² |
|---|---|---|
| 40.22 | FBA 1.00 | 871.1 |
| 12.12 | ENO 0.96 | 233.7 |
| 7.177 | PFK 0.53, PYK 0.43 | 205.7 |
| 3.902 | GAPD 0.98 | 197.2 |
| 3.281 | PYK 0.53, PFK 0.44 | 122.1 |
| 1.776 | PGM 0.96 | 60.14 |
| 1.673 | PGK 0.95 | 44.48 |
| 1.107 | TPI 1.00 | 20.5 |
| 1.016 | PGI 1.00 | 3.077 |
| 1.000 | LDH_L 1.00 | −0.013 |

- **The likelihood column is the separate spectrum of the block's marginal likelihood
  precision.** Its rows are sorted on their own, so a row does not have to be the same
  direction as the eigenvalue beside it. The two spectra's eigenvectors differ, because the
  priors' widths differ.
- **The counts do not hinge on the threshold.** The whitened eigenvalues nearest 2 are 1.78
  and 3.28. The likelihood's smallest, −0.013, is indistinguishable from zero: it comes from
  subtracting a prior precision of 400, and the second-difference checks agree to about 0.08
  absolute. The 3.08 above it is signal.
- **The full spectrum's bottom 10 equal the block's marginal spectrum** to 4e-4 relative. The
  metabolic block is where the data run out.
- **The likelihood column is the block's marginal posterior precision less its prior's.** That
  is not in general the likelihood's own Schur complement on the block: the script's comment
  at line 236 overstates it. Here the two agree to 1.6e-4 relative, and both count 9 of 10,
  because the block barely couples to the other 19 coordinates. The likelihood's smallest
  eigenvalue on those 19 is 339. `_reads.jl` prints both.

The full −H spectra, whitened and likelihood-only, are in the archive's `spectrum.md`, beside
KIM-178's narrow block, along with each coordinate's MAP, prior SD and Laplace SD.

## Reproduce

- `dev/scripts/corea_eigenspectrum/corea_eigenspectrum.sh <outdir>`, run from a worktree at a
  pushed commit. The script writes `spectrum.md`.
- `julia --project=. dev/scripts/corea_eigenspectrum/corea_eigenspectrum_reads.jl <outdir or
  archive>` prints the loadings, the threshold margins, the Newton decrement and the shift from
  KIM-178's MAP. It takes a minute on the login node.
