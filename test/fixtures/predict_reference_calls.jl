# The calls behind test/fixtures/predict_reference.jl (KIM-68), and what is kept of
# each. Both dev/scripts/make_predict_reference.jl, run at the frozen commit, and
# test/test_predict_reference.jl, run at the checked-out tree, include this file,
# so the two run the same calls. It uses only the API as it stood at the M4 freeze
# (`93abbce`).
#
# The calls cover every path KIM-68 touches, at 60 s:
# - A: the M0 hybrid, completion = :none;
# - B: the reduced model under V with V's own cells and a knockdown, and its
#   variance splits;
# - B0: the same draws and seeds, with V reading μ̂'s cells and no knockdown;
# - C: the M0 mean field with δ from an ensemble, and a knockdown;
# - D: `forward_cells` and `intervention_cells`.

using InferCell
using InferCell: AbstractSubModel
using Random

const _PREF_H = 60.0

# A made-up posterior over M0's targets, near their nominal values: the one
# test_predict.jl uses, and the fixture stores its values.
function _pref_posterior(ms, n; seed = 27)
    rng = Xoshiro(seed)
    nominal = nominal_parameter_values(ms)
    X = [nominal[t] * exp(0.1 * randn(rng)) for _ in 1:n, t in M0_TARGETS]
    return PosteriorDraws(M0_TARGETS, X, 0.2 .+ 0.01 .* rand(rng, n), fill(1, n),
                          collect(1:n), nothing, nothing, (kind = :test,))
end

_pref_sorted(d::AbstractDict) = sort!([k => v for (k, v) in d]; by = kv -> string(kv[1]))

_pref_cell(c) = (; seed = c.seed, q = c.q, failure = c.failure, message = c.message,
                 drains = c.drains, n_clips = length(c.clips), minima = _pref_sorted(c.minima))
_pref_cells(cs::AbstractVector) = [_pref_cell(c) for c in cs]
_pref_cells(A::AbstractMatrix) = (; size = size(A), cells = [_pref_cells(cs) for cs in vec(A)])

function _pref_plain(p::Prediction)
    d = p.delta
    return (; draws = p.draws, held = p.held, theta_H = p.theta_H, seeds = vec(p.seeds),
            seeds_size = size(p.seeds), viability_seeds = vec(p.viability_seeds),
            viability_seeds_size = size(p.viability_seeds), cells = _pref_cells(p.cells),
            viability_cells = _pref_cells(p.viability_cells),
            cells_intervened = _pref_cells(p.cells_intervened), mu = p.mu, sigma2 = p.sigma2,
            n_ok = p.n_ok, effect = p.effect, effect_sigma2 = p.effect_sigma2, viable = p.viable,
            decisions = [(; x.viable, x.fraction, counts = _pref_sorted(x.counts))
                         for x in vec(p.decisions)],
            weights = p.weights, conditioning = p.conditioning, summary = p.summary,
            effect_summary = p.effect_summary, statuses = _pref_sorted(p.statuses),
            statuses_intervened = _pref_sorted(p.statuses_intervened),
            fork_seed = p.provenance.fork_seed,
            splits = p.viability === nothing && p.completion === :none ? nothing :
                     [w => variance_split(p; weights = w) for w in keys(p.weights)],
            delta = d === nothing ? nothing :
                    (; d.seeds, hybrid = _pref_cells(d.hybrid), mean_field = _pref_cell(d.mean_field),
                     d.hybrid_mean, d.hybrid_se, d.n_ok, d.delta,
                     hybrid_intervened = _pref_cells(d.hybrid_intervened),
                     mean_field_intervened = d.mean_field_intervened === nothing ? nothing :
                                             _pref_cell(d.mean_field_intervened),
                     d.delta_intervened))
end

# Each call as name => (threads -> its plain result).
function predict_reference_calls()
    ms = collect(AbstractSubModel, m0_models())
    dms = collect(AbstractSubModel, d11_models())
    pd = _pref_posterior(ms, 3)
    # At 60 s this floor rejects one of B's four pairs, and one of each row's two in
    # B0, so the intervened arm skips a pair and no row is left empty.
    V = ViabilityRule(; label = "toy: a floor some pairs miss", allow_clipped = true,
                      min_viable_fraction = 0.5, floors = [:M_gtp_c => 17_000.0])
    kd = Knockdown(; dose = 0.3)
    vkw = (; horizon = _PREF_H, completion = :viability_conditioned_prior, viability = V,
           n_held = 2, n_parameter_draws = 2, n_cells = 2, width = 2.0)
    vals = [n => pd.values[1, j] for (j, n) in enumerate(pd.names)]
    seeds = [6805, 6806, 6807]
    return [
        :A => t -> _pref_plain(predict(ms, pd; horizon = _PREF_H, n_cells = 2,
                                       rng = Xoshiro(6801), threads = t)),
        :B => t -> _pref_plain(predict(dms, pd; vkw..., intervention = kd,
                                       n_viability_cells = 2, rng = Xoshiro(6802), threads = t)),
        :B0 => t -> _pref_plain(predict(dms, pd; vkw..., n_viability_cells = 0,
                                        rng = Xoshiro(6802), threads = t)),
        :C => t -> _pref_plain(predict(ms, pd; horizon = _PREF_H, forward = :mean_field,
                                       delta = :ensemble, n_ref = 4, intervention = kd,
                                       rng = Xoshiro(6804), threads = t)),
        :D => t -> begin
            f = forward_cells(ms, haldane_values(metered(ms), vals), seeds; horizon = _PREF_H,
                              threads = t)
            i = intervention_cells(ms, vals, kd, seeds; horizon = _PREF_H, threads = t)
            (; forward = _pref_cells(f), baseline = _pref_cells(i.baseline),
             intervened = _pref_cells(i.intervened), delta = i.delta)
        end,
    ]
end

# The posterior the calls run on, kept in the fixture so a change to how it is
# drawn shows as itself.
function predict_reference_posterior()
    p = _pref_posterior(collect(AbstractSubModel, m0_models()), 3)
    return (; names = p.names, values = p.values, sigma_b = p.sigma_b)
end

# Equal, with floats to `rtol` and `atol` (NaN equal to NaN); everything else exactly.
_pref_match(a::AbstractFloat, b::AbstractFloat; tol...) =
    isequal(a, b) || isfinite(a) && isfinite(b) && isapprox(a, b; tol...)
_pref_match(a::AbstractArray, b::AbstractArray; tol...) =
    size(a) == size(b) && all(_pref_match(x, y; tol...) for (x, y) in zip(a, b))
_pref_match(a::NamedTuple, b::NamedTuple; tol...) =
    keys(a) == keys(b) && all(_pref_match(a[k], b[k]; tol...) for k in keys(a))
_pref_match(a::Tuple, b::Tuple; tol...) =
    length(a) == length(b) && all(_pref_match(x, y; tol...) for (x, y) in zip(a, b))
_pref_match(a::Pair, b::Pair; tol...) =
    _pref_match(a.first, b.first; tol...) && _pref_match(a.second, b.second; tol...)
_pref_match(a, b; tol...) = isequal(a, b)

# The first path at which two results differ, for the failure message. The default
# tolerance is the test's; the generator passes zero for an exact comparison.
function _pref_diff(a, b, path = ""; rtol = 1e-6, atol = 1e-9)
    _pref_match(a, b; rtol, atol) && return nothing
    if a isa NamedTuple && b isa NamedTuple && keys(a) == keys(b)
        for k in keys(a)
            d = _pref_diff(a[k], b[k], "$path.$k"; rtol, atol)
            d === nothing || return d
        end
    elseif (a isa AbstractArray && b isa AbstractArray && size(a) == size(b)) ||
           (a isa Tuple && b isa Tuple && length(a) == length(b))
        for (i, (x, y)) in enumerate(zip(a, b))
            d = _pref_diff(x, y, "$path[$i]"; rtol, atol)
            d === nothing || return d
        end
    elseif a isa Pair && b isa Pair
        d = _pref_diff(a.first, b.first, "$path.first"; rtol, atol)
        d === nothing || return d
        return _pref_diff(a.second, b.second, "$path.second"; rtol, atol)
    end
    return "$path: $(repr(a)) ≠ $(repr(b))"
end
