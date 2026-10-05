# Conditional sampler pieces

The building blocks of phase 16's three-block sampler (`spec/phases/16-recovery.md`). Path replay and the path log-density are on the [Handshake driver](handshake.md) page.

## Transcript bridges (task 16a.3)

Within one interval of constant rate constants, a gene's transcript count is a birth-death chain, `TranscriptChain(k, μ, cap)`: transcription adds one at rate `k`, decay removes one at rate `μ·m`, and births stop at `cap`. The data fix the count exactly at every 60 s save point, so the path between two saves is a bridge of that chain.

- **Transition probabilities.** `transition_matrix(chain, τ)` and `transition_probability(chain, τ, a, b)` come from the uniformised series. They are checked against a matrix exponential to 1e-12.
- **Bridges.** `sample_bridge(rng, chain, τ, a, b)` draws one exact bridge by uniformisation (Hobolth and Stone 2009) and returns its event times and ±1 changes. `BridgeTable(chain, τ, b)` precomputes what every bridge to `b` shares, for repeated draws.
- **The cap.** `bridge_cap(max_observed)` is the largest observed count plus 20, and it is a floor rather than a guarantee. `truncation_mass(chain, τ, a)` bounds what the cap drops, and `adequate_cap(k, μ, τ, a; floor)` grows the cap until that is below 1e-12.
- **The independent check.** `bridge_birth_distribution` is the exact distribution of a bridge's birth count, from a matrix exponential of the chain augmented with a birth counter. It is what the sampled bridges are tested against.

## One-dimensional updates (task 16a.4)

Every update works on `u = ln θ`, so a `LogNormal(μ, s)` prior enters as `Normal(u; μ, s)`, which is `logpdf(LogNormal, θ) + ln θ`.

- **`slice_step(rng, logf, x)`** is one univariate slice-sampling update: stepping out, then shrinkage.
- **Block 2.** A promoter strength or the decay constant scales its reactions' propensities linearly. Its conditional given the path is `n·u − e^u·A + log Normal(u; μ, s)`, from the firing count `n` and the exposure per unit `θ`, `A` (`rate_log_conditional`). `rate_update(rng, θ, n, exposure, prior; bound)` takes one slice step on it. It throws if the accepted value passes `bound`, where the propensity stops being linear, such as a promoter's turnover ceiling.

## Block 1 (task 16a.5)

Block 1 samples the free forward constants and the metabolite σ given every cell's path. It runs on the published model and uses no gradient (D16.1).
- **`Block1(models, cells; forwards, panel)`** holds the composition, the cells, the forward constants' names, the observed panel and its priors. Each cell is a `Block1Cell`: a t = 0 snapshot with the stochastic parameters written, the recorded path, and the observed panel.
- **`block1_replay(b, u)`** replays every cell at `θF = exp.(u)`. The reverse constants are derived through the nominal equilibrium constant (`derived_ode_values`). It returns the summed path log-density and the panel's squared log residuals against `max(x, 1)`, where `x` is the carry-inclusive particle count.
- **`block1_logtarget(b, u, σ)`** is the log target. The prior is on the forward constants only, as `Normal(u; μ, s)`, and a reverse constant contributes nothing (`block1_log_prior`).
- **`block1_update!(rng, b, state)`** is one sweep: a slice step per `u` coordinate, then one on `ln σ` from its conditional given the residuals (`sigma_log_conditional`), which needs no replay.

## Block 3: the path update (task 16a.7)

`csmc_sweep(rng, base, spec, data, ref; N, lag, schedule, bulk)` is one conditional SMC sweep over a cell's 60 s observation windows (D16.3). A particle is a driver snapshot, and every proposal is built outside the SSA, then replayed:
- **Transcripts** are exact bridges between the observed counts: the window's first 59 s at its starting constants, then the last second at the constants rebuilt at 60m − 1.
- **Every other reaction** is simulated exactly given those transcripts (`propose_events`).
- **The incremental weight** is the metabolite likelihood times each gene's transition-probability factor (`window_log_weight`, D16.3's formula). Those transition probabilities differ between particles.

`lag = 0` is particle Gibbs. `lag = L` is ancestor sampling truncated at `L` windows, exact only at `L ≥ T`. Each particle draws from its own generator, seeded from the sweep's. `TranscriptMap`, `CellData`, `CSMCSpec`, `WindowEvents`, `window_events` and `join_path` are its bookkeeping, and `advance_window!` moves one particle across one window.

**Bulk metabolites** (amended 2026-10-05; the production model). Metabolites are observed as one bulk measurement per pool per save, of the population's mean particle count (`observe_bulk`). Block 1 scores the cells' mean (`Block1(...; bulk)`). In block 3, `csmc_sweep(...; bulk)` takes one `BulkObservation` per window, made by `bulk_window_observations` from the bulk record and the other cells' summed latent. The particle is scored at the mean with its own latent in place. The bulk likelihood couples the cells, so a caller updates them in turn, refreshing the others' sum from the `latent` each sweep returns. The chosen setting is particle Gibbs at N = 20.

**Annealed windows** (amended 2026-10-02; superseded for production, off by default). Without annealing, the metabolite likelihood puts fresh proposals hundreds of nats below the reference within one window, and no variant mixes. With `schedule` (for example `geometric_schedule(K)`), each particle's window is annealed importance sampling, by `annealed_window`. The tempered targets are π_β ∝ q · w^β. Each stage keeps the path before a uniform cut and redraws the rest from the proposal's own conditional (`advance_window!` with `head`), accepted at min(1, (w′/w)^β). The particle's weight is the AIS weight. The reference keeps its path, and its auxiliary chain is drawn backwards through the same moves. It runs under particle Gibbs only, and costs about N(K + 1) cycles per sweep.

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["transcript_bridge.jl", "conditional_updates.jl", "organisms/coreA/block1.jl", "csmc.jl"]
```
