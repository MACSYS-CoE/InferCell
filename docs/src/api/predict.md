```@meta
CurrentModule = InferCell
```

# Prediction

Source: [`src/predict.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict.jl), [`src/posterior_draws.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/posterior_draws.jl), [`src/predict_completion.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_completion.jl), [`src/predict_forward.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_forward.jl), [`src/predict_viability.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_viability.jl), [`src/predict_intervention.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_intervention.jl), [`src/predict_variance.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/predict_variance.jl).

## The entry point

[`predict`](@ref) computes the predictive distribution of a quantity Q, given the
data, as `dev/notes/2026-10-29-cut-predictive-distribution.md` defines it. Its
arguments keep five things apart:
- the posterior source, a [`PosteriorDraws`](@ref);
- the forward description: `forward = :hybrid` or `:mean_field`, with `delta`
  for the mean field;
- the completion policy: `completion = :none`, `:prior` or
  `:viability_conditioned_prior`, with `held`, `width`, `n_held` and `viability`
  (and, under V, `n_viability_cells`, which must be given);
- the estimand, the population mean of Q over cells;
- the intervention.

Its layout is one crossed matrix: posterior draws θ_I as rows, held draws θ_H as
columns, the held draws shared by every row. Each pair runs its own fresh cells,
since the Monte Carlo correction E[σ²]/n holds only if no two draws share cell
seeds. An intervention's arm reuses its pair's seeds, so the effect is paired
within a draw. The [`Prediction`](@ref) keeps every run, every decision and every
seed. Its `summary` separates the spread of μ̂ over the pairs from the Monte
Carlo term and from the cell-to-cell spread. Without completion, their
difference is the parameter variance. With completion it is not: every row
shares the same held draws, so about 1/K of the held parameters' share is
missing from it for K held draws, and all of it at one held draw. On the crossed
matrix, [`variance_split`](@ref) gives the total and the θ_I/θ_H split in both
orders, from weighted U-statistics in which the cells' noise and the shared
columns cancel. Without a viability rule it is unbiased. Under one, the weights
are ratio estimates of the conditioned law, and the estimate is biased; a
simulation check (KIM-66) found that bias at up to about 6% at M = 100, K = 20.

What `predict` leaves to the freeze (KIM-36, KIM-39): which Q, the viability rule,
the dose, which conditioned variant is the headline, and which order of the split
is the headline. It takes the first three as arguments, and returns both
variants' weights and the matrix the split needs.

On M0 no parameter is held, so `completion = :none`. That is the scored case. On
the reduced model the held parameters come from their priors, and the result is
unscored: no truth matches those data.

The toy models' [`posterior_predictive`](@ref) reads its posterior the same way,
and forks the rng the same way. It draws without replacement, or by weight
with replacement from an ABC posterior. It runs the toy composition once per
draw, an ODE solve or, for a jump module, a stochastic run, and returns the
trajectories.

## Reading a posterior

[`posterior_draws`](@ref) reads any posterior InferCell writes into one
[`PosteriorDraws`](@ref): the exact route's V7 replicate files, the mixed route's
B7 chain files, the `Chains` that `infer(models, ds; route)` returns, and the toy
models' `Chains` and `ABCPosterior`. Every column is mapped to its parameter by
name, and every value is held on the natural scale, ready to write into a model.
The bulk noise σ_b is kept in its own field: it is observation noise, not a cell
parameter, and a forward simulation never reads it.

## Completing the held parameters

On the reduced model, M0's posterior covers five parameters, and every other free
parameter is held: no data reaches it. [`held_parameters`](@ref) lists them for a
composition, skipping the slots the driver overwrites (the `enz_R_*` catalytic
slots), the reverse constants, and point priors, and says why each was skipped.
[`draw_held`](@ref) draws them from their priors, which `width` widens or narrows
about their medians for the sensitivity runs. Reverse constants are derived from
the whole set that runs, posterior and held values together, by
[`haldane_values`](@ref), so every reaction keeps its equilibrium constant.

## Forward descriptions and δ

[`forward_cells`](@ref) runs a composition at given parameter writes, read for
the predicted quantity Q:
- `forward = :hybrid` runs the reference composition, one stochastic cell per
  seed, and each seed fixes its cell;
- `forward = :mean_field` runs route (ii), one ODE, with no cells and no
  single-cell distribution.

Q is read from the run's metered lactate export, so every composition is
[`metered`](@ref) first: `:lactate_out` is the lactate exported in particles and
`:exported` is the fraction of the lactate formed that left the cell. A function
of the run is accepted too. Each run is kept as a [`CellRun`](@ref): its Q, its
clips, the smallest particle count of each state, and its failure, if it failed,
with the reason.

The mean field has no δ for Q of its own. The mixed route's δ covers only its 17
panel pools. [`delta_correction`](@ref) computes one from a hybrid ensemble at a
reference point, once, and keeps it as a [`DeltaCorrection`](@ref) to reuse.

## Viability

A viability rule V is data, a [`ViabilityRule`](@ref), and nothing in the code
fixes its thresholds. V is frozen in M4, before any prediction is read.
[`cell_status`](@ref) puts each run in one of five categories, which are not
interchangeable: failed, starved, below a floor, clipped, or clean. V is a
property of θ: [`assess_viability`](@ref) accepts a parameter set when enough of
its cells are viable, and otherwise keeps its count of cells per category as the
reason. The hybrid reads V on the same cells that give μ̂, so where Q differs
between viable and unviable cells, and the viable fraction is near the threshold,
the conditioned μ̂ is biased. Whether V reads cells of its own is fixed with V.
Over the crossed matrix of posterior draws × held draws,
[`viability_weights`](@ref) gives both conditioned variants of the methods note
as weights: (b) keeps θ_I's marginal as the posterior's, and (a) lets the filter
act on θ_I. Each comes with its acceptance rate and effective sample size.

## Interventions

An [`Intervention`](@ref) is applied after the parameter draw. The first is a
[`Knockdown`](@ref) of the ptsG promoter, which scales `S_JCVISYN3A_0779` by
`1 − dose` and changes nothing else ([`intervene`](@ref)). It assumes the
intervention acts only through that parameter, and that the posterior fitted on
unperturbed data carries over. [`intervention_cells`](@ref) runs both arms on the
same seeds, so each cell is paired with itself and the effect ΔQ keeps little of
the intrinsic noise. A dose of 0 gives ΔQ = 0 exactly.
[`delta_correction`](@ref) covers the intervened arm too.

## Reference

```@autodocs
Modules = [InferCell]
Pages = ["predict.jl", "posterior_draws.jl", "predict_completion.jl", "predict_intervention.jl", "predict_forward.jl", "predict_viability.jl", "predict_variance.jl"]
```
