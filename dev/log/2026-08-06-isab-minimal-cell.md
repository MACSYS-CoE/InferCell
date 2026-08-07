# Session Log: ISAB Results Rebuilt on a True Three-Module Cell

**Date:** 2026-08-06
**Branch:** `minimal-cell-example` → PR #37
**Status:** Complete. Deck compiles, result measured, PR open.

## Why this session happened

The ISAB deck (12 Aug 2026) motivates itself with a cell-level number — slide 3
pushes a parameter range through a forward model, slide 8 asks "is doubling time
110 ± 5 min, or 60–200?" — and then delivered two abstract parameter corner
plots. The user's complaint was that the setup→results transition read badly in
practice, and that the "three-module" cell was really two, since
`StochasticGeneExpression` is the same four reactions as the ODE TX/TL block.

## The relabel that was considered and rejected

The first proposal was to keep the existing compute and simply describe the
stochastic module as a distinct subsystem. Rejected for a specific reason rather
than a general one: the boundary result is dramatic *because* the SSA block
shares all four of its parameters with the ODE block. Presented as a
between-module effect it would have claimed a transfer the experiment never
tested, and one question about shared-parameter counts would have undone it.

## What was found instead

The honest version already existed in the repo. `test/run_integ_iterative.jl`
builds metabolism + enzyme gene (ODE, state-coupled both ways) plus
`BurstyGeneExpression`, a telegraph-promoter regulator sharing only `k_tl`,
`gamma_mRNA`, `gamma_protein` — the ribosomes and degradation machinery. 11
unique parameters, 3 shared. No new framework code needed: `abc_smc` takes a
generic simulator closure and `sequential_infer` already does NUTS → KDE
boundary → ABC-SMC.

## What was learned

**The result holds, and is better than the one it replaces.** Production run
(job 15152954; 400 particles, 6 populations, 50 replicates, 1000 NUTS; 18 min):
the three shared parameters collapse onto truth at 13–187× tighter, the other
three barely move. Because only 3 of 11 are shared, the figure shows information
crossing the boundary *and stopping where the coupling stops* — something the
all-four-shared version could not show.

**The promoter rates are unidentifiable, not mis-inferred.** `k_on` (truth 0.2,
posterior ≈2.1) and `k_tx_burst` (truth 10, posterior ≈4.3) miss in *both* runs.
The summary statistics are replicate-mean trajectories, and a telegraph promoter
enters the mean only through `k_tx_burst·k_on/(k_on+k_off)`. That combination is
recovered — truth 2.857, alone 3.31 ± 2.99, conditioned 2.79 ± 0.24 — and the
boundary sharpens it 12×, which also explains why `k_tx_burst`'s own marginal
tightened 3× without being conditioned. This is a gap in
`src/summary_statistics.jl`, not in the model or the sampler.

**Smoke settings actively mislead here.** At 100 NUTS samples the differentiable
chain had not converged: it spanned 0.960–0.976 on `k_tx` with truth at 1.0 —
far too narrow and *excluding* truth. The boundary faithfully propagated that
overconfidence, producing conditioned sds of 0.001–0.004 that looked
spectacular. Inspecting the figures rather than the summary table is what caught
it. Do not quote the smoke table.

**Weighted quantiles matter for ABC figures.** Ranging corner-plot axes on
unweighted quantiles let near-zero-weight surviving particles set the limits and
squashed every posterior into the corner of its panel.

## Ruled out

- **Relabelling the existing figures** — see above.
- **Building a growth/division module** — the best story (its observable would
  be the doubling-time distribution the deck already asks about) but new
  modelling, new summary statistics and fresh synthetic data six days from the
  talk, with no fallback. Still the right move after the talk.

## Environment friction worth remembering

The Julia depot under `~/.julia` was empty (it was populated for the May runs).
Compute nodes have no route to the registry, so packages must be downloaded from
the login node with `JULIA_PKG_PRECOMPILE_AUTO=0` and precompiled inside the job
— 923 s for 392 dependencies. `module load` does not work in non-interactive
shells and the Julia binary is a wrapper needing `EBROOTJULIA`. The cluster's
system TeX lacks `beamer` entirely, so TinyTeX was installed at `~/.TinyTeX` to
verify the deck.

## Artefacts

- `test/run_talk_bursty_boundary.{jl,slurm}` — compute → CSV
- `dev/talks/ISAB/figures/plot_bursty_boundary.py` — CSV → figures (`ssa`,
  `ode`, `burst` modes)
- `dev/talks/ISAB/figures/isab_{ssa_posteriors,ode_corner,burst_rate}.png`
- `docs/superpowers/specs/2026-08-06-minimal-cell-boundary-figure-design.md` —
  design + both run logs with full tables
