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
- **The likelihood's own noise does not explain the shift.** The posterior barely
  correlates the promoters (mean pairwise 0.012, largest 0.177), so noise that the
  likelihood represents would scatter the promoters independently; the 3.06 SE already
  allows for it. A shared component that the likelihood does not represent (a fluctuation
  of this dataset in the pools every gene's transcription reads, or a mean-field bias)
  would leave these correlations unchanged. One dataset cannot separate those two.
  (Corrected in the PR #194 review; the first version read the 0.012 as evidence against
  shared noise of any kind.)
- **The size of the shift is small.** On the rate scale the promoters sit about 2% low on
  average.
- **ENO and FBA sit away from θ_ref for a known reason.** Their posterior means sit 2.3
  and 1.4 posterior SD from θ_ref. θ_ref is the bounded MAP at δ = S = 0
  (`corea_pilot_kim92.jl`, `prep`), not the MAP of the fitted likelihood, and ENO and FBA
  are learned from the bulk, where δ (up to 0.19 on the log scale, prep log) acts. So the
  two points need not coincide. Their marginals are also nearly symmetric (skewness 0.05
  and 0.06). Whether the optimiser stopped short is not tested; its return code was
  successful (prep log "θ_ref by map").
