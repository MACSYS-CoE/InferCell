# The archived V7 and B7 posteriors hold draws (KIM-18, 2026-10-29)

**Verdict:** both archives hold full draws, named and ordered as `M0_TARGETS`
plus σ_b. The truth for every dataset is stored beside them and agrees exactly
between the two archives. The crossed design (completion plan §3.5) needs no
refits.

**Evidence:** `dev/scripts/kim18_inspect_archives.jl` at `868d37c` (on main
`b58863d`), Slurm job 19509733, 2m17s. The log is
`/fred/oz022/tkimpson/infercell-kim18-archive-2026-10-29/kim18_inspect_19509733.log`.
Every file was deserialized on main, with no checkout of the commits that wrote
them.

## Exact: V7 (Gibbs)

`/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b/replicate<r>.jls`,
r = 1–50. Written by `dev/scripts/sbc_v7_16a9b.jl`.

- **What each file is:** `(; state, meta)`. `state` is a `gibbs_checkpoint`, and
  `state.trace` holds one entry per sweep: `(sweep, cme, u, σ, changed, seconds)`.
- **Draws:**
  - all 50 replicates reached 656 sweeps;
  - after V7's 100 burn-in sweeps, each has 556 draws × 6 parameters, one chain;
  - a draw on the log scale is `vcat(log.(t.cme), t.u, log(t.σ))`;
  - all finite.
- **Names:** `meta.names` is
  `[:S_JCVISYN3A_0607, :S_JCVISYN3A_0779, :krnadeg, :kcatF_R_ENO, :kcatF_R_FBA, :sigma_b]`.
  This is `M0_TARGETS` plus σ_b, in that order. The order comes from `g.cme` then
  `g.forwards`, so a reader should map by name rather than rely on it.
- **Truth:**
  - `meta.truth` holds it on the log scale, in `meta.names` order;
  - the seed is `meta.seed` = 16_200_000 + r;
  - for every target, the stored value is identical (difference 0) to the one
    `draw_truth(Xoshiro(seed), m0_models(); names = M0_TARGETS, purpose = :calibration)`
    regenerates.
- **What a file does not hold:** the dataset. It is rebuilt with
  `m0_dataset(seed; ncells = 50)`, or read from B7's `r<r>/prep.jls` (below).
- **Caveat:** the kept draws are autocorrelated. V7's merge log reports 23 of
  50 replicates with a smallest bulk ESS under 100. A consumer should treat the 556 draws as
  correlated, or thin them.

## Mixed: B7 (NUTS)

`/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/<fit>/chain_<c>.jls`.
Written by `dev/scripts/av6_89_b7.jl`.

- **The 116 files:** 58 fits × 2 chains.
  - r1–r50 are V7's datasets: 100 files.
  - v6 is `m0_v6`.
  - The other seven are diagnostic arms (`ks_v6`, `ks_low1`, `ks_low2`, `tx_v6`,
    `bulk_v6`, `null_v6`, `null_low1`): 14 files.
- **What each file is:** `(; fit, c, seed, start, start_tries, chain, seconds, nadapt, ndraw, commit, job)`.
- **Draws:**
  - `chain` is an MCMCChains `Chains` of 200 × 18 × 1, from 200 adaptation steps
    and 200 draws;
  - the parameters are `x[1]`–`x[6]`, already on the log scale, with NUTS
    internals in the other columns;
  - all 116 files hold 200 × 6 finite draws.
- **Names:** `x[i]` is `mixed_parameter_names(mixed_setup(ds))[i]`, which is
  `[:log_S_JCVISYN3A_0607, :log_S_JCVISYN3A_0779, :log_krnadeg, :log_kcatF_R_ENO, :log_kcatF_R_FBA, :log_sigma_b]`.
  That is the same order as V7.
- **Truth:**
  - each base fit's `prep.jls` holds the whole dataset `ds`, so the truth is
    `ds.truth`, with σ_b as `ds.sigma`;
  - for r1–r50 it is identical (difference 0, σ_b included) to V7's `meta.truth`
    for the same r;
  - the KS, stream and null arms have no `prep.jls` of their own, or one without
    `ds` (KS), and take theirs from their base fit (`fit_kind` in
    `av6_89_b7.jl`).
- **Writing commits:** the chain files were written at four commits (7a102fc
  ×86, 6526b45 ×20, 8da6045 ×6, 09be2f5 ×4). This is consistent with the run's resume history
  (`jobs_before_resume*.tsv`). None of it gets in the way of loading.

## For KIM-28 (`predict` ingesting these)

- Both formats deserialize on current main. The ingest can read them directly.
- For each r in 1–50, the exact and mixed posteriors share one dataset and one
  truth, so the crossed comparison pairs them by r.
- Draw counts per dataset:
  - exact: 556, from one correlated chain;
  - mixed: 400, from 2 chains × 200.

  Any matched-size comparison has to subsample.
