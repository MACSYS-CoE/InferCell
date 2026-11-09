# KIM-92 pre-registration: the reduced-model pilot's gates, design and scoring targets

Frozen at tag `prereg/kim-92` on 2026-11-10, before the pilot's dataset is drawn
and before any posterior on the reduced model exists. This is a snapshot of the
criteria, not a copy of the plan. The plan lives in Linear (project InferCell,
issue KIM-92; Spec D1–D4; Decision log 2026-11-09). A change to anything below
after this tag is a Decision log entry that quotes Was and Now, and says what
had been read by then.

The PPC threshold was drafted as "≤ 5% of bulk rows outside the central 99%
band" and marked as Tom's call. On 2026-11-10, before any result existed, the
KIM-92 session proposed freezing it as drafted and Tom gave the go-ahead.

## 1. The fit

- **Model.** `d11_models()`, the reduced syn3A composition, through the mixed
  route (route (ii) mean field for the bulk factor, exact transcript factors).
- **Dataset.** `corea_dataset(17_200_001; ncells = 200, horizon = COREA_CYCLE_S,
  fit_to = COREA_FIT_S, targets = COREA_TARGETS)`. That is 200 cells, run to
  6300 s, with the record cut to the 52 saves at 60–3120 s. The seed is the
  first of the injection block that `test/test_corea_mixed.jl` reserves for the
  pilot (17_200_000 + r). It is fixed here, before its truth or data are read
  at this size.
- **Free set: 21 coordinates**, `mixed_parameter_names(s)`. These are
  `log_<t>` for the 20 names of `COREA_TARGETS` (17 promoter strengths,
  krnadeg, kcatF_R_ENO, kcatF_R_FBA), then `log_sigma_b`. Every other constant
  stays at its default.
- **Setup.** `mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS)`,
  at its default tolerances: route (ii) runs at abstol 1e-10, reltol 1e-8, with
  at most 10^4 steps and a support box of ±6 prior SD. The bulk factor covers
  `M0_PANEL`, 17 pools × 52 saves = **884 rows**.
- **Reference θ, δ and S**, as in B5/B7 (spec §12 2026-10-23, "B5 as planned"):
  - θ_ref is `reference_theta` at δ = S = 0, maxiters 200.
  - The reference ensemble is 800 cells, `reference_seeds(ds, 800)`, run by
    `reference_ensemble(values, seeds; models = d11_models(), horizon = 3120.0)`.
  - δ and S come from `delta_and_S(cells, m_ref; C = 200)`, and both are held
    fixed for the chains.
- **Chains.** NUTS with target acceptance 0.8 and maximum tree depth 10
  (Turing's default).
  - 4 chains, each NADAPT = 200 adaptation steps then NDRAW = 200 kept draws.
  - Seeds are `mixed_chain_seeds(ds, 4; arm = 0)`.
  - Each chain starts at `mixed_jittered_start(rng, s, θ_ref.x; scale = 0.25)`.
- **One extension.** Each of the 4 chains continues for NDRAW = 200 more draws
  from its own final state, keeping the adapted step size and mass matrix, with
  no new adaptation. The gates are then recomputed on all 400 kept draws per
  chain. That is the only extension. Adding chains does not count as one.
- **Sizing may change, the gates may not.** If item 2's profile shows that
  this design cannot finish within the partition's limits, a smaller design
  (NADAPT, NDRAW, chain count or cell count) may replace it. That change goes
  in the Decision log before any chain starts, and §2 below does not change.
  Moving down the ladder (fewer parameters, fewer cells, a shorter window) is
  KIM-93's call, not this file's.

## 2. The gates (the kill criteria)

The pilot **fails** its kill check if, on the dataset above, either of these
holds.

1. **Convergence.** Some coordinate of the 21 has R̂ > 1.02 or bulk ESS < 100
   after at most one extension (§1).
   - R̂ is the rank-normalised split R̂, and ESS the bulk ESS, from
     MCMCDiagnosticTools over the 4 chains' kept draws (B7's `diagnostics`).
   - The gate is checked on the first 200 draws first. If it fails there, the
     chains are extended once and the gate is checked again. The result is the
     check after the extension, or the first check if no extension ran.
2. **PPC.** More than 5% of the 884 bulk rows lie outside their central 99%
   posterior predictive band. That is, **more than 44 rows outside** fails.
   - **The predictive** is the mixed route's own bulk model, the one the
     posterior was fitted with. For a posterior draw x = (log θ, log σ_b):
     log z_rep ~ N(m(θ) + δ, σ_b² I + S), where m(θ) =
     `mean_field_log_bulk` at θ, and δ and S are the fit's.
   - **The draws.** 200 posterior draws, evenly spaced through the pooled kept
     draws (all chains, after any extension). Each gets 10 replicate draws of
     log z_rep, so there are 2000 predictive samples per row.
   - **The band.** Each row's band runs from the 0.5% to the 99.5% empirical
     quantile of its 2000 samples. A row is outside when the observed log z is
     below the lower end or above the upper one.
   - This tests whether the fitted approximate likelihood reproduces its own
     data. It does not test the mean-field approximation against the hybrid.
     The approximation's protection is recovery (§3) and KIM-95's scored
     predictions (Spec, "Assumptions").
3. **Cost** has no cap (Tom, Decision log 2026-11-09). The pilot measures it,
   and KIM-93 decides whether it is affordable.

**Non-convergence from weak identifiability** is a fallback, not a stop. That
means marginals at their priors, or ridges or strong correlations in the pair
plots. Before calling a convergence failure, the result note shows each
coordinate's contraction, the pair plots of the worst-mixing coordinates, and
the posterior correlation matrix, and names the weak directions or shows there
are none. Which fallback to take is KIM-93's call.

## 3. Reported, not gated

- **Recovery.** For each of the 21 coordinates, whether the truth (`ds.truth`,
  and log σ for log σ_b) falls inside the central 90% posterior interval
  (5%–95% quantiles of the pooled kept draws). The note reports the count out
  of 21 and the standardised error (posterior mean − truth) / posterior SD.
  There is one truth, so this is not a coverage rate.
- **Contraction.** 1 − var_post / var_prior per coordinate, with the prior
  variance from `mixed_prior_sd`. This is B7's rule 8 shrinkage. For F4's
  caption, "the data reach" a parameter when its contraction is **≥ 0.5**.
  N is the count of such coordinates out of 21.
- **Cost.**
  - CPU-h by `sacct` (CPUTimeRAW) for every job and per stage: dataset and MAP,
    reference ensemble, chains, extension, PPC.
  - Wall-clock from first submission to the last job's end.
  - Item 2's profile: the log density and its 21-coordinate gradient at the
    pilot's setup and tolerances, median of 3, with time and peak RSS.

## 4. The scoring targets (scored in KIM-95, frozen here)

KIM-95 scores the pilot's posterior, if KIM-93 says go, against the true
simulator at the dataset's truth θ\*. The targets and scores are Spec D3's, on
D2's split:

- **Targets.**
  - **The later horizon.** `M0_PANEL`'s 17 pools and the 17 genes' transcripts
    at the saves in 3180–6300 s, the second half of the cycle, which the data
    never see. The truth runs use fresh seeds.
  - **Q** = `lactate_out` over one 6300 s cycle (`COREA_CYCLE_S`), as its
    population mean μ(θ) over cells (spec §12 2026-10-30, the M4 freeze).
  - **The ptsG knockdown.** Q under `Knockdown(; dose = 0.3)` on
    `S_JCVISYN3A_0779` (KIM-39), and ΔQ against baseline.
- **Scores**, for each target:
  - whether the truth falls inside the central 50% and 90% predictive
    intervals;
  - the CRPS;
  - the CRPS skill against the prior predictive, 1 − CRPS_post / CRPS_prior.
- **One truth**, so these are scores, not coverage rates. Poor scores are a
  finding, not a kill (Spec D4).

KIM-95 sets its own run sizes (cells per draw, draws, truth seeds) and records
them before its runs. It does not add, drop or redefine a target or a score
here without a Decision log entry dated after this tag.
