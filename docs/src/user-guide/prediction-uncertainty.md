# Prediction uncertainty over time

[`predict`](@ref) can keep time series, and you can ask it what happens to their
uncertainty when a fixed parameter is released to its prior.

## Recording a trajectory

Pass `trajectory`, a list of ODE-block species. Every run keeps them, in particles,
every `every` seconds and at the horizon:

```julia
using InferCell, Random

pd = posterior_draws("replicate1.jls")
p = predict(d11_models(), pd; horizon = COREA_CYCLE_S, completion = :prior,
            n_parameter_draws = 40, n_held = 20, n_cells = 4, rng = Xoshiro(1),
            trajectory = [InferCell.LAC_EXPORT_METER, :M_gtp_c], every = 60.0)
```

Recording reads what the run already holds, so it changes no cell, no draw and no
result. A species outside the ODE block is refused before any cell runs. On the
hybrid, `every` must be a whole number of handshakes.

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
the variance split at each time with that noise removed, and its last value for the
metered export is [`variance_split`](@ref)'s result for Q.

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

- **A reverse constant goes with its forward constant and its Michaelis constants.**
  To release a Michaelis constant, free the reaction's `kcatR_*` too. The build refuses
  otherwise and names the constant it needs.
- **A release is a sensitivity to priors, not a posterior.** The released constants are
  not fitted, and the inferred parameters were fitted with those constants at their
  nominal values, so the cut assumption is stressed harder the more you release.
- **Check how wide the prior is before you expect a wide band.** Balanced `kcatF`
  constants have a median prior geometric SD of 1.05, so releasing one barely moves a
  band. Michaelis constants have a median of 3.6. The example above shows the usage and
  does not claim that lactate dehydrogenase matters.

Released parameters can starve cells. A [`ViabilityRule`](@ref) conditions on the cells
that stay viable, and `trajectory_bands(p, s; weights = ones(M, K))` gives the band
without it.
