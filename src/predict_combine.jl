# `combine_rows` (KIM-40): one `Prediction` from row blocks run as separate calls.
#
# Without completion, a prediction's rows are independent posterior draws, so its
# rows can be split over calls, each with its own `rng` and disjoint `draws`. M0's
# evidence run does this, as single-core tasks. This puts the blocks back together:
# the runs are concatenated row-wise, and every per-row and summary result is
# recomputed from them over all the rows, by the helpers `_assemble` uses. No
# block's summary is read, so none is averaged (KIM-70).

# Two blocks' δ are one correction when everything that defines and gives it
# agrees. A δ read back from a file is a copy, so `===` cannot say this.
_same_delta(x::Nothing, y) = y === nothing
_same_delta(x::DeltaCorrection, y) =
    y isa DeltaCorrection &&
    all(f -> isequal(getfield(x, f), getfield(y, f)),
        (:values, :horizon, :quantity, :modules, :fingerprint, :seeds, :hybrid_mean,
         :hybrid_se, :n_ok, :delta, :intervention, :delta_intervened))

"""
    combine_rows(blocks) -> Prediction

One [`Prediction`](@ref) from row blocks of one prediction, run as separate
[`predict`](@ref) calls, as if one call had run all their rows in block order.

**What the blocks must be:**
- without completion (`completion = :none`), so the rows are independent
  posterior draws. With completion every row shares the same held draws, and one
  call must run the crossed matrix;
- run on one posterior, with the same composition, forward, estimand, Q, horizon,
  intervention, settings (`n_parameter_draws` aside) and InferCell version. These
  are checked as a `Prediction` records them: Q by its label, the posterior by its
  `source`, and the composition by its reduction labels and skipped parameters. A
  `Prediction` keeps no Q function and no fingerprint of the composition, so two
  closures with one label, or two compositions alike in those records, are not
  told apart. Run a prediction's blocks from one script at one commit;
- on disjoint `draws`, each block with its own `rng`. Two calls with the same
  `rng` draw the same cell seeds whatever their `draws`, and no two rows may share
  cell seeds (the E[σ²]/n correction). Blocks whose seeds repeat are refused;
- for a corrected mean field, given one shared [`DeltaCorrection`](@ref) as
  `delta`, computed once with [`delta_correction`](@ref). A block that computed its
  own (`delta = :ensemble`) is refused, and so are blocks whose δ differ.

**What it computes.** The per-row results (`mu`, `sigma2`, `n_ok`, `effect`,
`effect_sigma2`) are recomputed from the concatenated runs. The `summary`,
`effect_summary` and status counts are recomputed over all the rows, with
uniform weights, as `predict` computes them. The blocks' own summaries are never
read, so they are not averaged.

**Provenance:** `fork_seed` and `threads` become one entry per block,
`block_rows` gives each block's row count, `settings.n_parameter_draws` is the
total, and `pairs` is concatenated.

Two posteriors run on the same `rng`s and block shapes draw the same cell seeds,
since the seeds come from `rng` whatever the `draws`. That gives their rows
common cell seeds.
"""
function combine_rows(blocks::AbstractVector{Prediction})
    isempty(blocks) && throw(ArgumentError("no blocks to combine"))
    a = first(blocks)
    for p in blocks
        p.completion === :none || throw(ArgumentError(
            "a block was run with completion = :$(p.completion): its rows share held draws, " *
            "so one call must run the crossed matrix"))
        p.provenance.settings.delta === :ensemble && throw(ArgumentError(
            "a block computed its own δ (delta = :ensemble); compute one DeltaCorrection " *
            "with delta_correction and pass it to every block"))
    end
    same(f, what) = all(p -> isequal(f(p), f(a)), blocks) || throw(ArgumentError(
        "the blocks differ in their $what, so they are not one prediction's rows"))
    same(p -> p.forward, "forward description")
    same(p -> p.estimand, "estimand")
    same(p -> p.quantity, "quantity")
    same(p -> p.horizon, "horizon")
    same(p -> p.intervention, "intervention")
    same(p -> merge(p.provenance.settings, (; n_parameter_draws = 0)), "settings")
    same(p -> p.provenance.source, "posterior")
    same(p -> (p.provenance.reductions, p.provenance.skipped), "composition")
    same(p -> p.provenance.version, "InferCell version")
    all(p -> _same_delta(a.delta, p.delta), blocks) || throw(ArgumentError(
        "the blocks carry different δ; pass one DeltaCorrection to every block"))

    draws = reduce(vcat, [p.draws for p in blocks])
    allunique(reduce(vcat, [unique(p.draws) for p in blocks])) || throw(ArgumentError(
        "a posterior draw is in more than one block; the blocks' draws must be disjoint"))
    seeds = cat((p.seeds for p in blocks)...; dims = 3)
    allunique(seeds) || throw(ArgumentError(
        "two rows share cell seeds, as blocks run on the same rng do; give each block " *
        "its own rng"))

    cells = reduce(vcat, [p.cells for p in blocks])
    cells_i = reduce(vcat, [p.cells_intervened for p in blocks])
    M, K = size(cells)
    mom = [_moments([c.q for c in cells[i]]) for i in CartesianIndices(cells)]
    mu, s2, nok = first.(mom), getindex.(mom, 2), last.(mom)
    u = fill(1 / (M * K), M, K)
    weights = (; uniform = u)
    conditioning = (; uniform = (; acceptance = 1.0, ess = _kish(u)))
    dc = a.delta
    shift = dc === nothing ? 0.0 : dc.delta
    summary = map(w -> _summary(mu, s2, nok, w, shift, a.forward), weights)

    effect = fill(NaN, M, K)
    es2 = fill(NaN, M, K)
    esummary = nothing
    if a.intervention !== nothing
        enok = zeros(Int, M, K)
        for i in CartesianIndices(cells)
            isempty(cells_i[i]) && continue
            effect[i], es2[i], enok[i] = _moments(_paired_delta(cells[i], cells_i[i]))
        end
        eshift = dc === nothing ? 0.0 : dc.delta_intervened - dc.delta
        esummary = map(weights) do w
            (; _summary(effect, es2, enok, w, eshift, a.forward)...,
             n_failed_intervened = _n_failed_intervened(cells, cells_i, w))
        end
    end

    pv = a.provenance
    provenance = (; fork_seed = [p.provenance.fork_seed for p in blocks], pv.version,
                  settings = merge(pv.settings, (; n_parameter_draws = M)), pv.source,
                  pv.reductions, pv.skipped,
                  pairs = reduce(vcat, [p.provenance.pairs for p in blocks]),
                  threads = [p.provenance.threads for p in blocks],
                  block_rows = [length(p.draws) for p in blocks])
    return Prediction(a.forward, a.completion, a.estimand, a.quantity, a.horizon, draws,
                      reduce(vcat, [p.theta_I for p in blocks]),
                      a.sigma_b === nothing ? nothing : reduce(vcat, [p.sigma_b for p in blocks]),
                      a.held, a.theta_H, seeds, cells, mu, s2, nok, a.intervention, cells_i,
                      effect, es2, a.viability, reduce(vcat, [p.viable for p in blocks]),
                      reduce(vcat, [p.decisions for p in blocks]),
                      cat((p.viability_seeds for p in blocks)...; dims = 3),
                      reduce(vcat, [p.viability_cells for p in blocks]), weights, conditioning,
                      dc, summary, esummary, _status_counts(cells), _status_counts(cells_i),
                      provenance)
end

export combine_rows
