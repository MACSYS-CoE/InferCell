# Completing a posterior from priors (KIM-29): the held parameters θ_H.
#
# M0's posterior covers its five targets. On the reduced model every other free
# parameter is held: no data reaches it, so `predict` completes it from its prior
# (dev/notes/2026-10-29-cut-predictive-distribution.md §3.2). The cut is
# deliberate. D does not update θ_H, and θ_H does not feed back into θ_I.

"""
    held_parameters(models, inferred; driver = nothing) -> (; held, skipped)

The parameters a prediction completes from their priors: every free parameter
of `models` that is not among `inferred` (names, or a [`PosteriorDraws`](@ref)).
`held` is in composition order. Four kinds are skipped, and each is listed in
`skipped` as `(name, reason)`, so a prediction's provenance can say what it left
out:
- `:inferred`: the posterior has it;
- `:driver_written`: the driver overwrites it at every handshake, so a value
  set there is never read. The `enz_R_*` catalytic slots are of this kind
  ([`driver_written_params`](@ref));
- `:haldane_derived`: a reverse constant, which [`haldane_values`](@ref) derives
  from the values that run;
- `:point_prior`: a `Dirac` prior, with nothing to draw.

`driver` is the composition's built [`HandshakeDriver`](@ref); one is built over
a one-second span if none is passed.
"""
function held_parameters(models::AbstractVector{<:AbstractSubModel},
                         inferred::AbstractVector{Symbol}; driver = nothing)
    ms = collect(AbstractSubModel, models)
    d = driver === nothing ? build_problem(ms; tspan = (0.0, 1.0), complete = true) : driver
    d isa HandshakeDriver || throw(ArgumentError(
        "held_parameters needs a hybrid composition, whose driver lists the slots it writes"))
    written = Set(w.param_slot for w in driver_written_params(d))
    reverses = Set(r.reverse for r in haldane_relations(ms))
    free = unique_params(reduce(vcat, (model_free_params(parameters(m)) for m in ms)))
    held = Symbol[]
    skipped = NamedTuple{(:name, :reason), Tuple{Symbol, Symbol}}[]
    for p in free
        reason = p.name in inferred ? :inferred :
                 p.name in written ? :driver_written :
                 p.name in reverses ? :haldane_derived :
                 p.prior isa Dirac ? :point_prior : nothing
        reason === nothing ? push!(held, p.name) : push!(skipped, (name = p.name, reason))
    end
    return (; held, skipped)
end

held_parameters(models::AbstractVector{<:AbstractSubModel}, pd::PosteriorDraws; kwargs...) =
    held_parameters(models, pd.names; kwargs...)

"""
    draw_held(rng, models, held, n; width = 1.0) -> Vector{Vector{Pair{Symbol, Float64}}}

`n` independent prior draws of the parameters `held`, each as `name => value`
pairs, with priors `width` times as wide ([`draw_parameters`](@ref)). Reverse
constants are not in the draws. A prediction derives them from the whole set
that runs, posterior and held values together, with [`haldane_values`](@ref).
"""
function draw_held(rng::AbstractRNG, models::AbstractVector{<:AbstractSubModel},
                   held::AbstractVector{Symbol}, n::Integer; width::Real = 1.0)
    n >= 1 || throw(ArgumentError("n must be at least 1, got $n"))
    reverses = Set(r.reverse for r in haldane_relations(models))
    bad = [h for h in held if h in reverses]
    isempty(bad) || throw(ArgumentError(
        "$bad are reverse constants, which are derived, not drawn"))
    return [filter(p -> !(first(p) in reverses), draw_parameters(rng, models, held; width))
            for _ in 1:n]
end

export held_parameters, draw_held
