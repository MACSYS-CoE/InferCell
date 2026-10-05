"""
Exact per-gene transcript bridges (spec/phases/16-recovery.md task 16a.3).

Inside one interval of constant rate constants, one gene's transcript count is
a birth-death chain: transcription adds one at rate `k`, decay removes one at
rate `μ·m`. Translation reads the transcript and leaves it alone. Given the
counts at both ends of the interval, which the data fix exactly at every 60 s
save point, the path in between is a bridge of that chain, sampled here exactly
by uniformisation (Hobolth and Stone 2009). The transition probabilities the
particle weights need (D16.3) come from the same uniformisation.

The count space is truncated at `cap`, where births stop. The cap starts at the
dataset's largest observed count plus 20 ([`bridge_cap`](@ref)) and grows until
[`truncation_mass`](@ref), which bounds what the truncation drops, is below
1e-12 ([`adequate_cap`](@ref)).
"""

"""
    TranscriptChain(k, μ, cap)

One gene's transcript count over an interval of constant rate constants:
births at rate `k`, deaths at rate `μ·m`, on `0:cap`, with no birth out of
`cap`.
"""
struct TranscriptChain
    k::Float64
    μ::Float64
    cap::Int

    function TranscriptChain(k::Real, μ::Real, cap::Integer)
        k >= 0 && μ >= 0 || throw(ArgumentError("rates must be non-negative, got k = $k, μ = $μ"))
        cap >= 1 || throw(ArgumentError("the cap must be at least 1, got $cap"))
        new(Float64(k), Float64(μ), Int(cap))
    end
end

"""
    bridge_cap(max_observed) -> Int

The count cap for a dataset: its largest observed transcript count plus 20
(spec/phases/16-recovery.md D16.3).
"""
bridge_cap(max_observed::Integer) = Int(max_observed) + 20

"""
    adequate_cap(k, μ, τ, a; floor, tol = 1e-12) -> Int

The smallest cap at or above `floor` at which [`truncation_mass`](@ref) from
count `a` over `τ` is below `tol`. The dataset rule, [`bridge_cap`](@ref), is
the floor. It is not always enough: at k = 1/s and μ = 0.2/s over 20 s from a
count of 10, a cap of 30 drops 7.6e-12 (job 17758372). So the cap grows until
the bound holds (spec/phases/16-recovery.md D16.3).
"""
function adequate_cap(k::Real, μ::Real, τ::Real, a::Integer; floor::Integer, tol::Real = 1e-12)
    cap = max(Int(floor), a + 1)
    while truncation_mass(TranscriptChain(k, μ, cap), τ, a) >= tol
        cap += 5
        cap > 10_000 && error("no cap below 10,000 holds the truncation mass under $tol")
    end
    return cap
end

"The generator matrix of `bd`, indexed by count + 1."
function generator(bd::TranscriptChain)
    n = bd.cap + 1
    Q = zeros(n, n)
    for m in 0:bd.cap
        i = m + 1
        if m < bd.cap
            Q[i, i + 1] = bd.k
        end
        if m > 0
            Q[i, i - 1] = bd.μ * m
        end
        Q[i, i] = -sum(Q[i, :])
    end
    return Q
end

# The uniformised chain: rate Λ and one-step matrix R = I + Q/Λ.
function _uniformised(bd::TranscriptChain)
    Q = generator(bd)
    Λ = maximum(-Q[i, i] for i in axes(Q, 1))
    Λ > 0 || return 0.0, Matrix{Float64}(I, size(Q)...)
    return Λ, Matrix{Float64}(I, size(Q)...) + Q / Λ
end

# Poisson(λ) probabilities for n = 0, 1, ... until the tail beyond is below
# `tol`, computed in log space so a large λ does not underflow at n = 0.
function _poisson_weights(λ::Float64; tol = 1e-17)
    λ == 0 && return [1.0]
    w = Float64[]
    logp = -λ
    n = 0
    total = 0.0
    while true
        p = exp(logp)
        push!(w, p)
        total += p
        # Past the mode, stop once the remaining mass is negligible.
        n > λ && 1 - total < tol && break
        n += 1
        logp += log(λ) - log(n)
        n > 10λ + 200 && break
    end
    return w
end

"""
    transition_matrix(bd, τ) -> Matrix

`P_τ[a + 1, b + 1] = P(m(τ) = b | m(0) = a)`, by uniformisation.
"""
function transition_matrix(bd::TranscriptChain, τ::Real)
    τ >= 0 || throw(ArgumentError("τ must be non-negative, got $τ"))
    Λ, R = _uniformised(bd)
    w = _poisson_weights(Λ * τ)
    P = w[1] * Matrix{Float64}(I, size(R)...)
    Rn = Matrix{Float64}(I, size(R)...)
    for n in 2:length(w)
        Rn = Rn * R
        P .+= w[n] .* Rn
    end
    return P
end

"`P(m(τ) = b | m(0) = a)`, as [`transition_matrix`](@ref) gives it."
transition_probability(bd::TranscriptChain, τ::Real, a::Integer, b::Integer) =
    transition_matrix(bd, τ)[a + 1, b + 1]

"""
    truncation_mass(bd, τ, a) -> Float64

An upper bound on what the cap drops over `τ` from count `a`: the probability
that the chain, with the cap made absorbing, has reached it by `τ`. Under the
untruncated chain that is the probability of ever reaching the cap within the
interval, and every path the truncation changes does.
"""
function truncation_mass(bd::TranscriptChain, τ::Real, a::Integer)
    Q = generator(bd)
    Q[end, :] .= 0.0
    return (exp(Q * τ))[a + 1, end]
end

"""
    BridgeTable(bd, τ, b)

What every bridge ending at `b` after `τ` shares, computed once so repeated draws
cost only their own sampling: the uniformised rate and one-step matrix, the
Poisson weights of the jump count, and `v[n + 1] = R^n e_b`, the probability
that `n` uniformised steps take each count to `b`.
"""
struct BridgeTable
    bd::TranscriptChain
    τ::Float64
    b::Int
    R::Matrix{Float64}
    w::Vector{Float64}
    v::Vector{Vector{Float64}}
end

function BridgeTable(bd::TranscriptChain, τ::Real, b::Integer)
    0 <= b <= bd.cap || throw(ArgumentError("endpoint $b must lie in 0:$(bd.cap)"))
    Λ, R = _uniformised(bd)
    w = _poisson_weights(Λ * τ)
    v = Vector{Vector{Float64}}(undef, length(w))
    e = zeros(bd.cap + 1)
    e[b + 1] = 1.0
    v[1] = e
    for n in 2:length(w)
        v[n] = R * v[n - 1]
    end
    return BridgeTable(bd, Float64(τ), Int(b), R, w, v)
end

"""
    sample_bridge(rng, bd, τ, a, b) -> (times, deltas)
    sample_bridge(rng, table, a) -> (times, deltas)

One exact draw of the transcript path over `[0, τ)` given `m(0) = a` and
`m(τ) = b`: the event times, in order, and each event's change, `+1` for a
transcription and `-1` for a decay. By uniformisation:
1. the number of uniformised jumps, given both endpoints;
2. their times, as uniform order statistics;
3. the states from the start forward, each drawn in proportion to its one-step
   probability times its chance of still reaching `b`.

Virtual jumps are dropped. Throws if `b` cannot be reached from `a` in `τ`. The
second form reuses a [`BridgeTable`](@ref).
"""
sample_bridge(rng::AbstractRNG, bd::TranscriptChain, τ::Real, a::Integer, b::Integer) =
    sample_bridge(rng, BridgeTable(bd, τ, b), a)

function sample_bridge(rng::AbstractRNG, t::BridgeTable, a::Integer)
    bd, R, w, v = t.bd, t.R, t.w, t.v
    0 <= a <= bd.cap || throw(ArgumentError("endpoint $a must lie in 0:$(bd.cap)"))
    pn = [w[n] * v[n][a + 1] for n in eachindex(w)]
    total = sum(pn)
    total > 0 || throw(ArgumentError(
        "count $(t.b) cannot be reached from $a in $(t.τ) s at k = $(bd.k), μ = $(bd.μ)"))
    N = _draw_index(rng, pn, total) - 1
    times = sort!(rand(rng, N) .* t.τ)
    out_t = Float64[]
    out_d = Int[]
    probs = zeros(bd.cap + 1)
    x = a
    for i in 1:N
        left = N - i                       # steps still to take after this one
        for y in 0:bd.cap
            probs[y + 1] = R[x + 1, y + 1] * v[left + 1][y + 1]
        end
        y = _draw_index(rng, probs, sum(probs)) - 1
        if y != x
            push!(out_t, times[i])
            push!(out_d, y - x)
        end
        x = y
    end
    x == t.b || error("the bridge ended at $x rather than $(t.b)")
    return out_t, out_d
end

function _draw_index(rng::AbstractRNG, p, total)
    u = rand(rng) * total
    c = 0.0
    for i in eachindex(p)
        c += p[i]
        u < c && return i
    end
    return findlast(>(0), p)
end

"""
    bridge_birth_distribution(bd, τ, a, b; jmax) -> Vector{Float64}

The exact distribution of the number of births in a bridge from `a` to `b` over
`τ`, entries `j = 0:jmax`, from a matrix exponential of the chain augmented with
a birth counter. It is independent of [`sample_bridge`](@ref), which is what it
exists to check (spec/phases/16-recovery.md task 16a.3). The deaths are then
`j − (b − a)`.
"""
function bridge_birth_distribution(bd::TranscriptChain, τ::Real, a::Integer, b::Integer;
                                   jmax::Integer = 60)
    nm = bd.cap + 1
    idx(m, j) = j * nm + m + 1
    Q = zeros(nm * (jmax + 1), nm * (jmax + 1))
    for j in 0:jmax, m in 0:bd.cap
        i = idx(m, j)
        if m < bd.cap && j < jmax
            Q[i, idx(m + 1, j + 1)] = bd.k
        end
        m > 0 && (Q[i, idx(m - 1, j)] = bd.μ * m)
        Q[i, i] = -sum(Q[i, :])
    end
    row = exp(Q * τ)[idx(a, 0), :]
    p = [row[idx(b, j)] for j in 0:jmax]
    p[end] < 1e-12 * sum(p) || throw(ArgumentError(
        "jmax = $jmax holds $(p[end] / sum(p)) of the bridge's mass at its top; raise it"))
    return p ./ sum(p)
end

export TranscriptChain, bridge_cap, adequate_cap, transition_matrix, transition_probability,
       truncation_mass, BridgeTable, sample_bridge, bridge_birth_distribution
