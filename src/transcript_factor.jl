"""
The exact transcript factor of track B's approximate likelihood
(spec/phases/approximate-method.md, task B1).

Each gene's transcript count is a birth-death chain observed exactly every 60 s.
The rate constants are piecewise constant, and they change at 60m − 1 s on the
jump clock, since `handshake_step!` rebuilds before the jump step. So window m
runs its first 59 s at the constants `R_{m−1}` and its last second at `R_m`, and
its factor is

```
L_g,m = [ P_{k_{g,m−1}, μ_g}(59) · P_{k_{g,m}, μ_g}(1) ] (y_{m−1}, y_m)
```

the same product the cSMC's particle weights estimate (`window_log_weight`). It
is evaluated here by propagating a row vector through the uniformised chain, so
it costs a few tridiagonal products per window and never forms a matrix.

The factor has to be differentiable, since track B samples it with NUTS. Nothing
in it branches on a rate. The cap is fixed per dataset from the largest birth
rate the prior allows, and the number of uniformisation terms is set from the
primal rate alone, so a dual number only carries through arithmetic.
"""

"""
    TranscriptFactor(y0, y; kmax, μmin, τ = 60, switch = 1, tol = 1e-12)

One cell's transcripts, ready for [`transcript_loglik`](@ref). `y0` holds each
gene's count at t = 0 and `y` its observed counts, genes × windows. Windows are
`τ` long, and the rate constants change `switch` before each window's end.

Each gene's cap is fixed here, once, as the smallest at or above
[`bridge_cap`](@ref) at which [`truncation_mass`](@ref) over `τ` from the
gene's largest count is below `tol` at birth rate `kmax` and death rate
`μmin`, the corner of the prior's support where truncation is worst. Either may
be a scalar or one value per gene.
"""
struct TranscriptFactor
    y0::Vector{Int}
    y::Matrix{Int}
    caps::Vector{Int}
    τ::Float64
    switch::Float64

    function TranscriptFactor(y0::AbstractVector{<:Integer}, y::AbstractMatrix{<:Integer};
                              kmax, μmin, τ::Real = 60.0, switch::Real = 1.0,
                              tol::Real = 1e-12)
        G = size(y, 1)
        length(y0) == G || throw(DimensionMismatch(
            "y0 has $(length(y0)) genes and y has $G"))
        all(>=(0), y0) && all(>=(0), y) ||
            throw(ArgumentError("transcript counts must be non-negative"))
        0 < switch < τ || throw(ArgumentError("switch must lie in (0, τ), got $switch"))
        kmaxs = _per_gene(kmax, G, :kmax)
        μmins = _per_gene(μmin, G, :μmin)
        caps = map(1:G) do g
            top = max(y0[g], maximum(view(y, g, :); init = 0))
            adequate_cap(kmaxs[g], μmins[g], τ, top; floor = bridge_cap(top), tol = tol)
        end
        new(Vector{Int}(y0), Matrix{Int}(y), caps, Float64(τ), Float64(switch))
    end
end

function _per_gene(x, G, name)
    xs = x isa Real ? fill(Float64(x), G) : Float64.(collect(x))
    length(xs) == G || throw(DimensionMismatch("$name has $(length(xs)) entries for $G genes"))
    return xs
end

"""
    transcript_loglik(f, k, μ) -> Real
    transcript_loglik(fs, k, μ) -> Real

The exact log-likelihood of one cell's transcripts `f`, or the sum over cells
`fs`, given the birth rates `k` and the per-transcript death rates `μ`.

`k` is genes × (windows + 1). Column j holds the constants in effect from the
start of window j, `R_{j−1}`: window m runs at `k[:, m]` until `switch` before
its end and at `k[:, m + 1]` after. `μ` has one entry per gene,
`krnadeg / len_g`, and is constant. Both may be dual numbers.
"""
function transcript_loglik(f::TranscriptFactor, k::AbstractMatrix, μ::AbstractVector)
    G, M = size(f.y)
    size(k) == (G, M + 1) || throw(DimensionMismatch(
        "k is $(size(k)); it needs one row per gene and windows + 1 = $(M + 1) columns"))
    length(μ) == G || throw(DimensionMismatch("μ has $(length(μ)) entries for $G genes"))
    T = float(promote_type(eltype(k), eltype(μ)))
    s = zero(T)
    for g in 1:G, m in 1:M
        a = m == 1 ? f.y0[g] : f.y[g, m - 1]
        v = _window_distribution(f, g, a, k[g, m], k[g, m + 1], μ[g])
        s += log(v[f.y[g, m] + 1])
    end
    return s
end

transcript_loglik(fs::AbstractVector{TranscriptFactor}, k::AbstractMatrix, μ::AbstractVector) =
    sum(f -> transcript_loglik(f, k, μ), fs)

# Gene g's count distribution at a window's end, from count `a` at its start,
# at `k0` until `switch` before the end and `k1` after: e_aᵀ P_{k0}(τ − s) P_{k1}(s).
function _window_distribution(f::TranscriptFactor, g, a, k0, k1, μ)
    T = float(promote_type(typeof(k0), typeof(k1), typeof(μ)))
    v = zeros(T, f.caps[g] + 1)
    v[a + 1] = one(T)
    v = _propagate(v, k0, μ, f.caps[g], f.τ - f.switch)
    return _propagate(v, k1, μ, f.caps[g], f.switch)
end

# vᵀ exp(Q t) by uniformisation, for the chain with births `k` below `cap` and
# deaths `μ·c`. The rate Λ = k + μ·cap bounds every exit rate, so each step's
# matrix I + Q/Λ is non-negative and the series has no cancellation. The term
# count depends only on the primal Λt, so derivatives pass through untouched.
function _propagate(v::AbstractVector, k, μ, cap::Int, t::Float64)
    Λ = k + μ * cap
    _primal(Λ) > 0 || return v
    λ = Λ * t
    nterms = _poisson_terms(_primal(λ))
    T = eltype(v)
    out = exp(-λ) .* v
    cur = v
    logw = -λ
    for n in 1:nterms
        cur = _uniform_step(cur, k, μ, cap, Λ)
        logw += log(λ) - log(n)
        out = out .+ exp(logw) .* cur
    end
    return out
end

# One uniformised step, vᵀ(I + Q/Λ), using Q's tridiagonal form: count c gains
# from c − 1 by a birth and from c + 1 by a death, and loses at its exit rate.
function _uniform_step(v::AbstractVector, k, μ, cap::Int, Λ)
    w = similar(v, promote_type(eltype(v), typeof(k), typeof(μ), typeof(Λ)))
    for c in 0:cap
        exit = (c < cap ? k : zero(k)) + μ * c
        x = v[c + 1] * (1 - exit / Λ)
        c > 0 && (x += v[c] * k / Λ)
        c < cap && (x += v[c + 2] * μ * (c + 1) / Λ)
        w[c + 1] = x
    end
    return w
end

# Enough Poisson(λ) terms that the tail beyond is far below 1e-16: the mean
# plus 12 standard deviations, plus a margin for small λ.
_poisson_terms(λ::Real) = ceil(Int, λ + 12 * sqrt(λ) + 30)

_primal(x::Real) = x
_primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))

export TranscriptFactor, transcript_loglik
