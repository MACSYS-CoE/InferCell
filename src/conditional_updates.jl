"""
The one-dimensional conditional updates of the phase 16 sampler
(spec/phases/16-recovery.md D16.1 and D16.2): a univariate slice step, and block
2's update of a rate parameter from its sufficient statistics.

Every update works on `u = ln θ`. For a `LogNormal(μ, s)` prior the density of
`u` is `Normal(u; μ, s)`, which is `logpdf(LogNormal, θ) + ln θ`. Leaving out the
`+ ln θ` would shift the prior on `u` by `−s²` (D16.1).
"""

"""
    slice_step(rng, logf, x; w = 1.0, max_steps = 100) -> Float64

One univariate slice-sampling update of `x` under the unnormalised log-density
`logf` (Neal 2003): a level below `logf(x)`, an interval of width `w` stepped
out until both ends fall below the level, then shrunk until a draw lands above
it. It leaves `exp(logf)` invariant. Throws if `logf(x)` is not finite.
"""
function slice_step(rng::AbstractRNG, logf, x::Real; w::Real = 1.0, max_steps::Integer = 100)
    f0 = logf(x)
    isfinite(f0) || throw(ArgumentError("the slice step starts at x = $x, where logf is $f0"))
    level = f0 - randexp(rng)
    left = x - w * rand(rng)
    right = left + w
    j = floor(Int, max_steps * rand(rng))
    k = max_steps - 1 - j
    while j > 0 && logf(left) > level
        left -= w
        j -= 1
    end
    while k > 0 && logf(right) > level
        right += w
        k -= 1
    end
    while true
        x1 = left + rand(rng) * (right - left)
        logf(x1) > level && return x1
        x1 < x ? (left = x1) : (right = x1)
    end
end

"""
    rate_log_conditional(n, A, prior) -> Function

Block 2's log-density for `u = ln θ`, where `θ` scales the propensities of a set
of reactions linearly (spec/phases/16-recovery.md D16.2):

```
n·u − e^u·A + log Normal(u; μ, s)
```

`n` is how many times those reactions fired and `A` their exposure per unit
`θ`: the replay's integrated propensity divided by the `θ` it ran at. The prior
is the `LogNormal(μ, s)` on `θ`, carried to `u` with its Jacobian.
"""
function rate_log_conditional(n::Integer, A::Real, prior::LogNormal)
    n >= 0 && A >= 0 || throw(ArgumentError("need n ≥ 0 and A ≥ 0, got n = $n, A = $A"))
    μ, s = params(prior)
    return u -> n * u - exp(u) * A + logpdf(Normal(μ, s), u)
end

"""
    rate_update(rng, θ, n, exposure, prior; bound = Inf, w = 1.0) -> Float64

Block 2's update of one rate parameter. `n` and `exposure` are the firings and
integrated propensity of the reactions `θ` scales, from a replay run at `θ`, so
the exposure per unit `θ` is `exposure / θ`. One slice step on `u = ln θ` under
[`rate_log_conditional`](@ref).

`bound` is the value beyond which `θ` no longer enters linearly: for a promoter
strength, the turnover ceiling `TURNOVER_CEILING / RNAPOL_KCAT`. An accepted
value past it throws, since the conditional was wrong there.
"""
function rate_update(rng::AbstractRNG, θ::Real, n::Integer, exposure::Real,
                     prior::LogNormal; bound::Real = Inf, w::Real = 1.0)
    θ > 0 || throw(ArgumentError("θ must be positive, got $θ"))
    logf = rate_log_conditional(n, exposure / θ, prior)
    θ1 = exp(slice_step(rng, logf, log(θ); w))
    θ1 < bound || throw(ArgumentError(
        "the update accepted θ = $θ1, past $bound, where the propensity stops being " *
        "linear in θ and block 2's conditional does not hold"))
    return θ1
end

export slice_step, rate_log_conditional, rate_update
