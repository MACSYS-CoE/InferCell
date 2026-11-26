# KIM-161: the promoter shift and the ENO/FBA marginals

Output of `dev/scripts/pilot_promoter_shift.py`, run on the login node on 2026-11-26 against
`/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17`. Its inputs match that archive's
`SHA256SUMS`:

- `pilot_coordinates.tsv`: c7288568…a8ed
- `pilot_draws.tsv`: c91f7d93…8694

Both were written by `dev/scripts/corea_pilot_kim92.jl` at de24b27 (KIM-92, PR #185).

```
## 1. Promoters' standardised errors (n = 17)
mean z -0.808, SD 0.758, below truth 14 of 17
SE independent 0.243 -> mean/SE -3.33
SE with the posterior's correlations 0.264 -> mean/SE -3.06
## 2. Mean pairwise posterior correlation among the promoters: 0.012 (min -0.082, max 0.177)
## 3. Mean shift in log rate: -0.0225, factor exp(shift) - 1 = -2.22%
## 4. Skewness of the marginals (pooled draws, log scale)
  log_kcatF_R_ENO    skewness 0.046
  log_kcatF_R_FBA    skewness 0.056
```

**Reading:**

- **The shift is not an artefact of correlated errors.** Allowing for the posterior's
  correlations moves the promoters' mean standardised error from 3.33 to 3.06 SE from zero.
- **Shared noise has little support.** The posterior barely correlates the promoters: the
  mean pairwise correlation is 0.012 and the largest is 0.177. Noise that every promoter
  shared would need the data to tie them together, and the posterior shows no such tie.
- **The size of the shift is small.** On the rate scale the promoters sit about 2% low on
  average.
- **ENO and FBA are nearly symmetric**, with skewness 0.05 and 0.06. Their posterior means
  sit 2.3 and 1.4 posterior SD from θ_ref, and that offset is not skew. θ_ref is the
  bounded MAP at δ = S = 0 (`corea_pilot_kim92.jl`, `prep`), not the MAP of the fitted
  likelihood, so the two points need not coincide.
