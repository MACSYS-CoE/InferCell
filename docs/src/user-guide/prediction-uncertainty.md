# Prediction uncertainty over time

[`predict`](@ref) can keep time series, and you can ask it what happens to their
uncertainty when a fixed parameter is released to its prior.

**What these bands are.** The examples use a posterior fitted to M0's synthetic data and
carry it to the reduced syn3A model, where no truth exists. The result is unscored, so the
bands show what the framework's assumptions imply, and they are not calibrated predictions
for a real cell. See [the reduced syn3A model](reduced-syn3a.md).

**What they cost.** Each example runs 3,200 cells at 6,300 s, about 58 core-hours at
65 core-seconds per cell. Run them through a scheduler, with `--gcthreads=N,1` on many
threads ([`predict`](@ref) explains why), or try a short `horizon` and few draws first.

## Recording a trajectory

Pass `trajectory`, a list of ODE-block species. Every run keeps them, in particles,
every `every` seconds and at the horizon:

```julia
using InferCell, Random

pd = posterior_draws("replicate1.jls")   # a saved fit of M0; the API page lists what it reads
p = predict(d11_models(), pd; horizon = COREA_CYCLE_S, completion = :prior,
            n_parameter_draws = 40, n_held = 20, n_cells = 4, rng = Xoshiro(1),
            trajectory = [InferCell.LAC_EXPORT_METER, :M_gtp_c], every = 60.0)
```

On the hybrid, recording reads what the run already holds, so it changes no cell, no draw
and no result, and `every` must be a whole number of handshakes. On the mean field, `every`
must stay at 60, because the saves set what its minima, and so a viability floor, can see. A
species outside the ODE block is refused before any cell runs.

## Three layers of uncertainty

[`trajectory_bands`](@ref) returns empirical quantiles over time, in three layers:

```julia
b = trajectory_bands(p, InferCell.LAC_EXPORT_METER; baseline = true)
b.rows.lo[0.9], b.rows.hi[0.9]    # the inferred parameters
b.pairs.lo[0.9], b.pairs.hi[0.9]  # the inferred and held parameters
b.cells.lo[0.9], b.cells.hi[0.9]  # a new cell
```

`baseline = true` subtracts each cell's value at t = 0, so the metered export gives
the exported lactate Q(t).

The `pairs` layer is a quantile of each pair's mean over its cells, so it also holds
the Monte Carlo noise of those few cells, about σ̂²/n. [`trajectory_split`](@ref) gives
the variance split at each time with that noise removed. With `baseline = true`, its last
value for the metered export is [`variance_split`](@ref)'s result for Q. Early in a run the
noise can outweigh the spread: `flagged` marks the times where a component came out
negative, and a share there can exceed 1 or be `Inf` or `NaN`. A share is a number between 0 and 1
where `flagged` and `mc_undefined` are false and `total > 0`, with `NaN` still possible under very
sparse weights. That does not make it signal: `total` can still be cell noise, so compare it with
`mc`, the noise removed, and read a share only where `total` also exceeds it.

A cell that fails holds no trajectory, so it leaves the bands at every time and is counted
in `b.n_failed`. The bands at early times are therefore those of the cells that finished.

## Releasing a parameter to its prior

Most of a composition's kinetic constants are fixed. To see what a fixed constant's
uncertainty does, free it in the module's constructor. `predict` then draws it from
its prior with the other held parameters:

```julia
ms = d11_models()
ms[1] = CentralGlycolysis(enzymes = :translated,
                          free = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA,
                                  :kcatF_R_LDH_L, :kcatR_R_LDH_L])
p = predict(ms, pd; horizon = COREA_CYCLE_S, completion = :prior, n_parameter_draws = 40,
            n_held = 20, n_cells = 4, rng = Xoshiro(1),
            trajectory = [InferCell.LAC_EXPORT_METER])
```

The first four names are the reduced model's existing free constants, which a new
`free` list must keep. The last two release lactate dehydrogenase. Its forward
constant joins the held set, and its reverse constant is derived from the equilibrium
constant, not drawn ([`held_parameters`](@ref) says so in `skipped`).

Three things to keep in mind:

- **Free a reaction's reverse constant with its forward constant or a Michaelis constant.**
  `predict` derives `kcatR_*` and writes it, so it has to be free. Without it, the
  constructor passes and `predict` refuses when the first cell writes it, naming the
  constant.
- **A release is a sensitivity to priors, not a posterior.** The released constants are
  not fitted, and the inferred parameters were fitted with those constants at their
  nominal values, so the cut assumption is stressed harder the more you release.
- **Check how wide the prior is before you expect a wide band.** Balanced `kcatF`
  constants have a median prior geometric SD of 1.05, so releasing one barely moves a
  band. Michaelis constants have a median of 3.6. The example above shows the usage and
  does not claim that lactate dehydrogenase matters.

Released parameters can starve cells. To condition on that, run `predict` with
`completion = :viability_conditioned_prior`, a [`ViabilityRule`](@ref) as `viability`, and
`n_viability_cells`. The rule accepts or rejects a parameter pair, by the fraction of its
cells that stay viable. It does not filter cells, so a starved cell in a kept pair stays in
every band. The bands then default to the `:per_theta_I` weights, and
`weights = :joint` gives the other variant ([`viability_weights`](@ref)).
`trajectory_bands(p, s; weights = ones(size(p.mu)))` gives species `s`'s band without the
rule.
