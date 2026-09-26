"""
Metabolic control analysis on a composed ODE right-hand side (spec §3 check 8,
§11 task 14.8; reused by task 15.5 for output F2).

The summation theorems say that at a steady state, for a set of parameters
`αⱼ` that each multiply one reaction's whole rate and together cover every
reaction:

```
Σⱼ ∂ln Jₖ/∂ln αⱼ = 1     for every flux Jₖ          (flux control)
Σⱼ ∂ln xᵢ/∂ln αⱼ = 0     for every concentration xᵢ (concentration control)
```

Everything here is derived from the right-hand side alone, so no rate law is
restated and nothing can drift from what the model integrates:

- A **multiplier** names the parameter slots that are scaled together to scale
  one reaction. For the modular rate law that is `kcatF` and `kcatR` together,
  which is the same as scaling `E`.
- The **channel matrix** `G = ∂f/∂ln α` at `ln α = 0` has column `j` equal to
  `Nⱼ·vⱼ`: reaction `j`'s stoichiometry times its rate. Its range is the range
  of `N`, so its left null space is the conservation matrix `L`. It is taken at
  a *non-steady* state, where no column vanishes.
- The steady state solves the square system `[Qᵀf(x); L(x − x₀)] = 0`, where
  `Q` is an orthonormal basis of `range(N)`. The Jacobian of that system is
  nonsingular at a stable steady state even though `∂f/∂x` is not.
- The **implicit function theorem** then gives `dx/dln α = −A⁻¹[Qᵀ∂f/∂ln α; 0]`,
  exact to the solve.
- **Flux** `k` is read on a species `sₖ` it moves as
  `Jₖ = f_{sₖ}(x) − f_{sₖ}(x; αₖ = 0)`, which is `N_{sₖ,k}` times the reaction's
  rate. It needs every multiplier to scale its own reaction linearly and no
  other, which [`control_coefficients`](@ref) asserts. Then `∂Jₖ/∂ln αⱼ` at
  fixed `x` is `Jₖ·δₖⱼ`, and a constant multiple has the same log-derivative,
  so `∂ln Jₖ/∂ln α = δ + (∂ln Jₖ/∂x)·dx/dln α` is the flux control
  coefficient exactly, with no nested differentiation.
"""

"""
    ControlProblem(f, u0, p, names; held = Symbol[], multipliers)

One steady-state control analysis:
- `f(u, p)` is the composed right-hand side, with `u` a static vector;
- `u0` is the starting state, and the reference for every conserved total;
- `p` is the parameter vector;
- `names` are the state names in `f`'s layout;
- `held` are states fixed at their `u0` value and excluded from the solve.
  These are boundary species, such as external lactate, or accumulators that
  have no steady state;
- `multipliers` is a vector of `name => slots`: the slots of `p` scaled
  together to scale one reaction.
"""
struct ControlProblem{F}
    f::F
    u0::Vector{Float64}
    p::Vector{Float64}
    names::Vector{Symbol}
    internal::Vector{Int}
    multipliers::Vector{Pair{Symbol, Vector{Int}}}
end
function ControlProblem(f, u0, p, names::AbstractVector{Symbol};
                        held::AbstractVector{Symbol} = Symbol[],
                        multipliers)
    length(u0) == length(names) || throw(ArgumentError(
        "$(length(names)) state names for a state of length $(length(u0))"))
    for s in held
        s in names || throw(ArgumentError("held state :$s is not in the composition"))
    end
    internal = [i for (i, s) in enumerate(names) if !(s in held)]
    mults = Pair{Symbol, Vector{Int}}[Symbol(k) => collect(Int, v) for (k, v) in multipliers]
    for (k, v) in mults
        isempty(v) && throw(ArgumentError("multiplier :$k scales no slot"))
        all(i -> 1 <= i <= length(p), v) || throw(ArgumentError(
            "multiplier :$k names a slot outside the parameter vector"))
    end
    return ControlProblem(f, collect(Float64, u0), collect(Float64, p),
                          collect(Symbol, names), internal, mults)
end

_nstate(cp::ControlProblem) = length(cp.u0)

# The parameter vector with multiplier j's slots scaled by exp(a[j]).
function _scaled(cp::ControlProblem, a::AbstractVector{T}) where {T}
    q = promote_type(T, Float64).(cp.p)          # a copy, never cp.p itself
    for (j, (_, slots)) in enumerate(cp.multipliers)
        s = exp(a[j])
        for i in slots
            q[i] *= s
        end
    end
    return q
end

# The full state with the internal states replaced by x.
function _full(cp::ControlProblem, x::AbstractVector{T}) where {T}
    u = promote_type(T, Float64).(cp.u0)         # a copy, never cp.u0 itself
    u[cp.internal] .= x
    return SVector{_nstate(cp)}(u)
end

_f_int(cp::ControlProblem, x, a) = cp.f(_full(cp, x), _scaled(cp, a))[cp.internal]

"""
    channel_matrix(cp, x; a = zeros) -> Matrix

`∂f/∂ln α` over the internal states at internal state `x`: column `j` is
multiplier `j`'s reaction, `Nⱼ·vⱼ`.
"""
channel_matrix(cp::ControlProblem, x::AbstractVector;
               a::AbstractVector = zeros(length(cp.multipliers))) =
    ForwardDiff.jacobian(b -> _f_int(cp, x, b), a)

"""
    conservation_matrices(cp; rtol = 1e-9) -> (Q, L, rank)

`Q`, an orthonormal basis of the range of the stoichiometry, and `L`, whose
rows span the conserved combinations of the internal states. Both come from
the SVD of the channel matrix at `u0`. A column that is zero there, a reaction
with no rate at the starting state, would hide part of the range, so it
throws, naming the multiplier.
"""
function conservation_matrices(cp::ControlProblem; rtol::Real = 1e-9)
    x0 = cp.u0[cp.internal]
    G = channel_matrix(cp, x0)
    for (j, (k, _)) in enumerate(cp.multipliers)
        all(iszero, view(G, :, j)) && throw(ArgumentError(
            "Multiplier :$k has zero rate at the starting state, so its stoichiometry " *
            "cannot be read from the channel matrix there. Start from a state where " *
            "every reaction runs"))
    end
    # Normalise each column so the rank does not depend on how fast a reaction is.
    Gn = G ./ [maximum(abs, view(G, :, j)) for j in axes(G, 2)]'
    F = svd(Gn; full = true)
    r = count(>(rtol * F.S[1]), F.S)
    return (Q = F.U[:, 1:r], L = permutedims(F.U[:, r+1:end]), rank = r)
end

"""
    steady_state(cp; relax = 1e5, h0 = 1e-4, newton_tol = 1e-13, max_steps = 10_000)
        -> NamedTuple

Relax from `u0` by pseudo-transient continuation, then polish with Newton's
method on `[Qᵀf(x); L(x − x₀)] = 0`, which keeps every conserved total at its
`u0` value. The relaxation takes implicit-Euler steps, `(I − hJ)Δ = h·f`,
starting at `h0` and doubling the step after each one that keeps every state
non-negative, until the time reaches `relax`. Implicit Euler is L-stable, so
the block's stiffness does not bound the step, and it preserves every linear
invariant exactly. Newton then starts from a point inside the basin.

Returns:
- `x`, the internal steady state, and `u`, the full state;
- `residual`, `max|f|` over the internal states at `x`;
- `scale`, the largest reaction rate there, `max|G|`, for reading the residual;
- `iterations`, Newton's;
- `eigenvalues` of the Jacobian restricted to `range(N)`, `Qᵀ(∂f/∂x)Q`. All
  negative real parts means the steady state is stable within its conserved
  class.

Throws if either stage fails, rather than return a point that is not a
steady state.
"""
function steady_state(cp::ControlProblem; relax::Real = 1e5, h0::Real = 1e-4,
                      newton_tol::Real = 1e-13, max_steps::Integer = 10_000,
                      max_iter::Integer = 50)
    (; Q, L) = conservation_matrices(cp)
    a0 = zeros(length(cp.multipliers))
    fx(y) = _f_int(cp, y, a0)
    x0 = cp.u0[cp.internal]
    x = copy(x0)
    t, h, n = 0.0, Float64(h0), 0
    while t < relax
        n += 1
        n > max_steps && throw(ErrorException(
            "The relaxation did not reach t = $relax in $max_steps steps (t = $t, h = $h)"))
        Jx = ForwardDiff.jacobian(fx, x)
        y = x .+ (I - h * Jx) \ (h .* fx(x))
        if all(>=(0), y)
            x, t, h = y, t + h, 2h
        else
            h /= 4
            h < 1e-12 && throw(ErrorException(
                "The relaxation cannot keep the state non-negative (t = $t)"))
        end
    end
    F(y) = vcat(Q' * fx(y), L * (y - x0))
    it = 0
    for k in 1:max_iter
        it = k
        dx = ForwardDiff.jacobian(F, x) \ F(x)
        x .-= dx
        maximum(abs, dx) <= newton_tol * max(maximum(abs, x), 1.0) && break
        k == max_iter && throw(ErrorException(
            "Newton did not converge in $max_iter iterations; last step $(maximum(abs, dx))"))
    end
    G = channel_matrix(cp, x)
    Jx = ForwardDiff.jacobian(fx, x)
    return (x = x, u = collect(_full(cp, x)), residual = maximum(abs, fx(x)),
            scale = maximum(abs, G), iterations = it,
            eigenvalues = eigvals(Q' * Jx * Q), Q = Q, L = L)
end

"""
    control_coefficients(cp, ss; zero_flux = 1e-9) -> NamedTuple

The control coefficients at the steady state `ss` from [`steady_state`](@ref):
- `ccc[i, j] = ∂ln xᵢ/∂ln αⱼ` over the internal states. A state that is exactly
  zero has no log-derivative, and is left out and listed in `zero_states`;
- `fcc[k, j] = ∂ln Jₖ/∂ln αⱼ` over the multipliers' own fluxes. A flux below
  `zero_flux` times the largest is not a flux anything controls, and its row
  is `NaN` and listed in `zero_fluxes`;
- `flux[k]`, each `Jₖ` in the unit of `f`, with the species `sₖ` it is read on;
- `total_response[i, m] = ∂ln xᵢ/∂Tₘ`, the response to conserved total `m`.
  This is how a carrier total, which multiplies no rate, is reported;
- `ccc_sum` and `fcc_sum`, the row sums the summation theorems fix at 0 and 1.
"""
function control_coefficients(cp::ControlProblem, ss; zero_flux::Real = 1e-9)
    (; x, Q, L) = ss
    nm = length(cp.multipliers)
    a0 = zeros(nm)
    Jx = ForwardDiff.jacobian(y -> _f_int(cp, y, a0), x)
    G = channel_matrix(cp, x)
    A = vcat(Q' * Jx, L)
    dx = -(A \ vcat(Q' * G, zeros(size(L, 1), nm)))          # dx/dln α
    dT = A \ vcat(zeros(size(Q, 2), size(L, 1)), Matrix{Float64}(I, size(L, 1), size(L, 1)))

    nz = [i for i in eachindex(x) if x[i] != 0]
    ccc = dx[nz, :] ./ x[nz]

    # Each flux on the species its reaction moves most, at the steady state,
    # read as f minus f with that reaction switched off.
    sp = [argmax(abs.(view(G, :, k))) for k in 1:nm]
    off(k) = (b = zeros(nm); b[k] = -Inf; b)            # exp(-Inf) = 0
    Jk(y) = [(_f_int(cp, y, a0) - _f_int(cp, y, off(k)))[sp[k]] for k in 1:nm]
    flux = Jk(x)
    # The reading is valid only if each multiplier scales its own reaction
    # linearly: then switching it off removes exactly its channel column.
    for k in 1:nm
        d = _f_int(cp, x, a0) - _f_int(cp, x, off(k))
        err = maximum(abs, d - view(G, :, k))
        err <= 1e-8 * maximum(abs, view(G, :, k)) + 1e3 * eps() * maximum(abs, G) || throw(ArgumentError(
            "Multiplier :$(first(cp.multipliers[k])) does not scale one reaction " *
            "linearly: switching it off moves f by $err more than its channel"))
    end
    dJdx = ForwardDiff.jacobian(Jk, x)
    fcc = Matrix{Float64}(I, nm, nm) .+ (dJdx * dx) ./ flux
    big = maximum(abs, flux)
    zf = [k for k in 1:nm if abs(flux[k]) < zero_flux * big]
    fcc[zf, :] .= NaN

    names = cp.names[cp.internal]
    mnames = first.(cp.multipliers)
    return (ccc = ccc, ccc_states = names[nz], zero_states = names[setdiff(eachindex(x), nz)],
            fcc = fcc, multipliers = mnames, flux = flux, flux_species = names[sp],
            zero_fluxes = mnames[zf], total_response = dT[nz, :] ./ x[nz],
            ccc_sum = vec(sum(ccc; dims = 2)), fcc_sum = vec(sum(fcc; dims = 2)))
end

"""
    assert_summation(cc; tol = 1e-6) -> NamedTuple

Check 8: every concentration's control coefficients sum to zero, and every
nonzero flux's to one, within `tol`. Throws an `ArgumentError` naming the worst
state or flux, so a mutation shows which identity broke. Returns the worst
deviation of each.
"""
function assert_summation(cc; tol::Real = 1e-6)
    dc = abs.(cc.ccc_sum)
    i = argmax(dc)
    dc[i] <= tol || throw(ArgumentError(
        "Check 8 fails: the concentration control coefficients of :$(cc.ccc_states[i]) " *
        "sum to $(cc.ccc_sum[i]), not 0 (tolerance $tol)"))
    live = [k for k in eachindex(cc.fcc_sum) if !(cc.multipliers[k] in cc.zero_fluxes)]
    df = abs.(cc.fcc_sum[live] .- 1)
    k = argmax(df)
    df[k] <= tol || throw(ArgumentError(
        "Check 8 fails: the flux control coefficients of the :$(cc.multipliers[live[k]]) " *
        "flux sum to $(cc.fcc_sum[live[k]]), not 1 (tolerance $tol)"))
    return (ccc = dc[i], fcc = df[k])
end

"""
    omit_multiplier(cc, name) -> NamedTuple

`cc` with multiplier `name`'s column removed and the row sums recomputed. The
mutation of check 8: spec §3 says the flux identity fails without `k_chg`, and
this is how a test shows it does, and by how much.
"""
function omit_multiplier(cc, name::Symbol)
    j = findfirst(==(name), cc.multipliers)
    j === nothing && throw(ArgumentError("no multiplier :$name"))
    keep = [k for k in eachindex(cc.multipliers) if k != j]
    ccc = cc.ccc[:, keep]
    fcc = cc.fcc[keep, keep]
    return merge(cc, (ccc = ccc, fcc = fcc, multipliers = cc.multipliers[keep],
                      flux = cc.flux[keep], flux_species = cc.flux_species[keep],
                      ccc_sum = vec(sum(ccc; dims = 2)), fcc_sum = vec(sum(fcc; dims = 2))))
end

"""
    group_coefficients(C, multipliers, groups) -> (Matrix, Vector{Symbol})

Sum the columns of a control coefficient matrix over groups of multipliers, for
example a gene's reactions. `groups` is a vector of `name => [multiplier, …]`.
Every group member must be a multiplier. A multiplier in no group is left out.
"""
function group_coefficients(C::AbstractMatrix, multipliers::AbstractVector{Symbol}, groups)
    cols = map(groups) do (g, ms)
        idx = map(ms) do m
            j = findfirst(==(m), multipliers)
            j === nothing && throw(ArgumentError("group :$g names unknown multiplier :$m"))
            j
        end
        vec(sum(C[:, idx]; dims = 2))
    end
    return reduce(hcat, cols), Symbol[first(g) for g in groups]
end

export ControlProblem, channel_matrix, conservation_matrices, steady_state,
       control_coefficients, assert_summation, omit_multiplier, group_coefficients
