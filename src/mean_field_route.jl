"""
Route (ii): the reduced syn3A model with gene expression swapped to its mean
field and run as one ODE (spec/phases/approximate-method.md task B4, §12
2026-10-23).

The hybrid driver executes the coupling edges between two blocks at each 1 s
handshake. Here there is one block, so each edge becomes part of the
right-hand side or a callback:

- **Catalytic:** a rate law reads the mean protein count directly, converted to
  mM at the live volume. The count is not rounded.
- **Deferred counters** become continuous fluxes: the counter's derivative,
  divided by the live particles-per-mM factor, is debited from or credited to the
  pool. Nothing is clamped, and a debited pool below zero stops the solve.
- **The rate-constant rebuild** is a discrete callback at `k · interval_r − Δ`,
  where Δ is the reference's 1 s handshake. That is B1's switch time. It reads
  the ODE's own pools at that instant, where the hybrid, whose ODE runs one
  handshake ahead, reads them Δ later. The rebuilt slots are also recomputed at
  t = 0 from the initial pools and the current parameters.
- **Growth** runs continuously. The membrane count N is an extra state; the
  volume follows from it by the hybrid's own law, and every state the hybrid
  dilutes gets a `−c · dV/dt / V` term.

The edges are lowered by the hybrid's own functions, with the mean-field modules
in the place of the jump block, so the two paths refuse the same compositions
with the same messages. Each departure is declared by
[`route_declarations`](@ref).
"""

"""
    mean_field_swap(models) -> Vector{AbstractSubModel}

`models` with every `:jump` module replaced by one [`MeanField`](@ref) group
holding all of them, coupling edges kept, placed after the `:ode` modules. A
module with no reactions, such as `HeldEnzymeCounts`, joins the group and its
counts are held. A composition with no `:jump` module is returned unchanged.
"""
function mean_field_swap(models::AbstractVector{<:AbstractSubModel})
    jumps = AbstractSubModel[m for m in models if formalism(m) === :jump]
    rest = AbstractSubModel[m for m in models if formalism(m) !== :jump]
    isempty(jumps) && return collect(AbstractSubModel, models)
    return AbstractSubModel[rest..., MeanField(jumps; edges = :declared)]
end

"""
    rate_constants(p, t, m::MeanField, pools)

Each member's rebuilt constants from its own slice of `p` and its own pools, in
the order of `rebuilt_params(m)`. `pools` is in the order of the group's inbound
rate-constant edges, which is the members' edges in member order.
"""
function rate_constants(p, t, m::MeanField, pools)
    out = Any[]
    k = 0
    for (j, mem) in enumerate(m.members)
        n = count(e -> e isa RateConstantEdge && is_consumer(e), coupling(mem))
        if !isempty(rebuilt_params(mem))
            own = SVector{n}(ntuple(i -> pools[k + i], n))
            append!(out, rate_constants(view(p, m.pidx[j]), t, mem, own))
        end
        k += n
    end
    return [v for v in out]
end

"""
    MeanFieldRoute

A route (ii) build: `prob`, the `ODEProblem` with its callbacks, and the lowered
records, which share their types and field names with [`HandshakeDriver`](@ref)'s
so that [`unexecuted_edges`](@ref) reads either.

State indices are into `prob.u0`, whose last entry is the membrane count when the
cell grows (`n_idx`, else 0). Parameter indices are into `prob.p`, whose names are
`param_names`.
"""
struct MeanFieldRoute{P}
    prob::P
    models::Vector{AbstractSubModel}
    state_names::Vector{Symbol}
    param_names::Vector{Symbol}
    catalytic::Vector{CatalyticExchange}
    geometry::Vector{GeometryExchange}
    debits::Vector{DeferredDebit}
    contributions::Vector{ContributionRecord}
    rebuilds::Vector{RateConstantRebuild}
    growth::Vector{GrowthChain}
    dilute_idxs::Vector{Int}
    n_idx::Int
    interval::Float64
    rebuild_times::Vector{Float64}
    area_baseline_nm2::Float64
    footprint_nm2::Float64
    volume_cap_litres::Float64
    initial_area_nm2::Float64
    initial_volume_litres::Float64
    radius0_nm::Float64
end

# The cell's geometry at membrane count `n`: the hybrid's `_update_volume!` and
# `_geometry_value` as functions, generic in the number type so a gradient
# passes through them. With no growth chain the cell is fixed.
function _route_geometry(rt, n)
    if rt.n_idx == 0
        v = rt.initial_volume_litres
        return (; factor = particles_per_mM(v), area = rt.initial_area_nm2,
                radius = rt.radius0_nm, volume = v, capped = false)
    end
    area = surface_area_nm2(rt.area_baseline_nm2, n, rt.footprint_nm2)
    radius = radius_from_area_nm(area)
    v = cell_volume_litres(radius)
    capped = v >= rt.volume_cap_litres
    capped || return (; factor = particles_per_mM(v), area, radius, volume = v, capped)
    v = oftype(v, rt.volume_cap_litres)
    r = radius_from_volume_nm(v)
    return (; factor = particles_per_mM(v), area = 4 * π * r^2, radius = r, volume = v, capped)
end

function _geometry_quantity(geo, q::Symbol)
    q === :radius_nm && return geo.radius
    q === :volume_litres && return geo.volume
    q === :area_nm2 && return geo.area
    error("Internal: :$q is not one of $VOLUME_QUANTITIES")
end

# What the right-hand side substitutes into one module's parameters: catalytic
# counts at local slots `cat_k` from states `cat_u`, and geometric quantities at
# local slots `geo_k`.
struct _RouteSubs{C, G}
    cat_k::C
    cat_u::C
    geo_k::G
    geo_q::Vector{Symbol}
end

_substitute(p::SVector, ::_RouteSubs{SVector{0, Int}, SVector{0, Int}}, u, geo) = p
function _substitute(p::SVector{N}, s::_RouteSubs, u, geo) where {N}
    out = SVector{N, promote_type(eltype(p), eltype(u))}(p)
    for j in eachindex(s.cat_k)
        out = setindex(out, counts_to_mM(u[s.cat_u[j]], geo.factor), s.cat_k[j])
    end
    for j in eachindex(s.geo_k)
        out = setindex(out, _geometry_quantity(geo, s.geo_q[j]), s.geo_k[j])
    end
    return out
end

struct _RouteRHS{M, S, B, R}
    models::M
    slots::S
    subs::B
    route::R
    drain_counter::Vector{Int}
    drain_pool::Vector{Int}
    drain_coef::Vector{Float64}
    flagged::Vector{Int}
end

function (f::_RouteRHS)(du, u, p, t)
    rt = f.route
    geo = _route_geometry(rt, rt.n_idx == 0 ? zero(eltype(u)) : u[rt.n_idx])
    parts = map(f.models, f.slots, f.subs) do m, sl, sb
        u_local = u[sl.sidx]
        p_local = _substitute(p[sl.pidx], sb, u, geo)
        u_inputs = u[sl.iidx]
        (dynamics(u_local, p_local, t, m, u_inputs),
         contributions(u_local, p_local, t, m, u_inputs))
    end
    fill!(du, zero(eltype(du)))
    _write_parts!(du, f.slots, parts)
    # The drains, after every module's own term, so a counter's derivative is
    # its whole accrual rate.
    for k in eachindex(f.drain_pool)
        du[f.drain_pool[k]] += f.drain_coef[k] * du[f.drain_counter[k]] / geo.factor
    end
    rt.n_idx == 0 && return nothing
    # Growth. The membrane count changes by the flagged states' reaction terms
    # alone, since dilution conserves counts; the volume is A^{3/2} below the
    # cap, so dV/dt / V = 3/2 · footprint · dN/dt / A.
    ndot = zero(eltype(du))
    for i in f.flagged
        ndot += du[i]
    end
    ndot *= geo.factor
    du[rt.n_idx] = ndot
    geo.capped && return nothing
    g = 1.5 * rt.footprint_nm2 * ndot / geo.area
    for i in rt.dilute_idxs
        du[i] -= u[i] * g
    end
    return nothing
end

@inline _write_parts!(du, ::Tuple{}, ::Tuple{}) = du
@inline function _write_parts!(du, slots::Tuple, parts::Tuple)
    sl, (d, c) = first(slots), first(parts)
    for j in eachindex(sl.sidx)
        du[sl.sidx[j]] += d[j]
    end
    for j in eachindex(sl.cidx)
        du[sl.cidx[j]] += c[j]
    end
    return _write_parts!(du, Base.tail(slots), Base.tail(parts))
end

# One rebuild's constants from the pools in `u`, in the order of `r.fill_idxs`.
function _route_rebuilt_values(r, p, u, t)
    pools = [k == 0 ? r.pool_held[j] : u[k] for (j, k) in enumerate(r.pool_idxs)]
    vals = rate_constants(view(p, r.p_idxs), t, r.model, pools)
    length(vals) == length(r.fill_idxs) || error(
        "rate_constants() for $(r.declared_by) returned $(length(vals)) values " *
        "but rebuilt_params() names $(length(r.fill_idxs))")
    return vals
end

# Every rebuild, from the pools in `u` now, into `p`.
function _route_rebuild!(p, u, t, rt::MeanFieldRoute)
    for r in rt.rebuilds
        for (k, v) in zip(r.fill_idxs, _route_rebuilt_values(r, p, u, t))
            p[k] = v
        end
    end
    return nothing
end

# A debited pool this far below zero, in particles, has been overdrawn rather
# than missed by the solver's tolerance.
const _OVERDRAW_PARTICLES = 0.01

function _route_callbacks(rt::MeanFieldRoute, debited::Vector{Int})
    times = Set(rt.rebuild_times)
    rebuild = DiscreteCallback(
        (u, t, integ) -> t in times,
        integ -> (_route_rebuild!(integ.p, integ.u, integ.t, rt); u_modified!(integ, true));
        initialize = (c, u, t, integ) -> (_route_rebuild!(integ.p, u, t, rt);
                                          u_modified!(integ, true)),
        save_positions = (false, false))
    isempty(debited) && return CallbackSet(rebuild)
    overdrawn(u, i) = u[i] * _route_geometry(rt, rt.n_idx == 0 ? 0.0 : u[rt.n_idx]).factor <
                      -_OVERDRAW_PARTICLES
    guard = DiscreteCallback(
        (u, t, integ) -> any(i -> overdrawn(u, i), debited),
        integ -> begin
            i = findfirst(i -> overdrawn(integ.u, i), debited)
            error("Route (ii): the pool :$(rt.state_names[debited[i]]) went below " *
                  "zero at t = $(integ.t) s. Its drains are unclamped " *
                  "(spec/phases/approximate-method.md §12 2026-10-23), so a pool " *
                  "the costs outrun is an error here, where the hybrid would clip it")
        end;
        save_positions = (false, false))
    return CallbackSet(rebuild, guard)
end

"""
    build_mean_field_problem(models; tspan = (0.0, 100.0), interval = 1.0,
                             radius_nm = nothing, initial_surface_area_nm2 = nothing,
                             footprint_nm2 = nothing, complete = false) -> MeanFieldRoute

Route (ii) of the approximate method: `models` with its `:jump` modules swapped to
their mean field ([`mean_field_swap`](@ref)), every coupling edge executed inside
one ODE, as this file's header describes. `interval` is the reference's handshake,
which sets the rebuild's switch time (`k · interval_r − interval`). The geometry
keywords and `complete` are the hybrid's.

The plain [`build_problem`](@ref) on an all-ODE composition is unchanged, and
still refuses catalytic edges.

Solve with [`solve_mean_field`](@ref). The rebuilt, catalytic and geometric slots
of `prob.p` are written by the route, not read from it
([`mean_field_written_params`](@ref)).
"""
function build_mean_field_problem(models::AbstractVector{<:AbstractSubModel};
                                  tspan = (0.0, 100.0), interval = 1.0,
                                  radius_nm = nothing, initial_surface_area_nm2 = nothing,
                                  footprint_nm2 = nothing, complete = false)
    interval > 0 || throw(ArgumentError("The interval must be positive, got $interval"))
    for m in models
        formalism(m) === :jump && !isempty(membrane_protein_states(m)) && throw(ArgumentError(
            "Module $(module_id(m)) is a :jump module flagging membrane proteins " *
            "$(membrane_protein_states(m)). Route (ii) grows the cell from membrane " *
            "proteins held as concentrations by :ode modules, as the reduced syn3A " *
            "model's are; a flagged count inside a mean-field group is not supported"))
    end
    ms = mean_field_swap(models)
    for m in ms
        formalism(m) === :ode || throw(ArgumentError(
            "Route (ii) runs :ode modules and mean-field groups only; " *
            "$(module_id(m)) is :$(formalism(m))"))
    end
    # The mean-field groups stand where the hybrid's jump block stands, after the
    # metabolic modules, so the hybrid's lowering applies as written and the
    # metabolic indices it returns are already global.
    met = AbstractSubModel[m for m in ms if !(m isa MeanField)]
    mfs = AbstractSubModel[m for m in ms if m isa MeanField]
    ms = AbstractSubModel[met..., mfs...]
    isempty(mfs) && throw(ArgumentError(
        "Route (ii) needs a mean-field group; this composition has no :jump module " *
        "and no MeanField. Use build_problem"))

    _validate_shared_params(ms)
    resolve_coupling(ms; complete)
    contexts = _build_contexts(ms)
    _resolve_coupling(ms, contexts; validate = false)
    nmet = sum(length(states(m)) for m in met; init = 0)
    met_ctx, mf_ctx = contexts[1:length(met)], contexts[(length(met) + 1):end]

    catalytic, geometry, debits, _ = _lower_exchanges(ms, met, mfs, met_ctx)
    rebuilds = _lower_rebuilds(met, mfs, mf_ctx, Float64(interval))
    growth, dilute = _lower_growth(ms, met, mfs)
    # The lowering numbers the mean-field states from 1; they follow the
    # metabolic ones.
    catalytic = [CatalyticExchange(c.count_idx + nmet, c.param_idx, c.species,
                                   c.param_slot, c.declared_by) for c in catalytic]
    for b in debits
        b.counter_idx += nmet
    end

    # The geometry, as `_build_hybrid_problem` states it.
    if isempty(growth)
        for (name, value) in ((:initial_surface_area_nm2, initial_surface_area_nm2),
                              (:footprint_nm2, footprint_nm2))
            value === nothing || throw(ArgumentError(
                "$name = $value was passed, but no module flags a membrane protein, " *
                "so the cell does not grow and the value would never be read"))
        end
        radius0 = radius_nm === nothing ? COREA_INITIAL_RADIUS_NM : Float64(radius_nm)
        area0 = 4 * π * radius0^2
    else
        radius_nm === nothing || throw(ArgumentError(
            "This composition flags membrane proteins, so its geometry comes from " *
            "initial_surface_area_nm2; pass that, not radius_nm"))
        area0 = initial_surface_area_nm2 === nothing ?
                COREA_INITIAL_SURFACE_AREA_NM2 : Float64(initial_surface_area_nm2)
        radius0 = radius_from_area_nm(area0)
    end
    footprint = footprint_nm2 === nothing ? MEMBRANE_PROTEIN_FOOTPRINT_NM2 : Float64(footprint_nm2)
    volume0 = cell_volume_litres(radius0)
    factor0 = particles_per_mM(volume0)

    u0 = _collect_ic_values(ms)
    names = reduce(vcat, (states(m) for m in ms); init = Symbol[])
    n_idx = 0
    base = 0.0
    if !isempty(growth)
        all(g -> g.block === :ode, growth) || error(
            "Internal: a membrane protein lowered to the jump side in route (ii)")
        n0 = sum(u0[g.idx] * factor0 for g in growth)
        base = membrane_area_baseline_nm2(area0, n0, footprint)
        base >= 0 || throw(ArgumentError(
            "The initial membrane proteins come to more than the initial surface area"))
        push!(u0, n0)
        push!(names, :membrane_count)
        n_idx = length(u0)
    end

    tend = Float64(tspan[2])
    times = sort(unique(Float64[k * r.interval - interval for r in rebuilds
                                for k in 1:floor(Int, (tend + 1e-9) / r.interval)
                                if k * r.interval - interval > tspan[1]]))
    p0 = _build_p0(ms)
    rt0 = MeanFieldRoute(nothing, ms, names, _block_param_names(ms), catalytic, geometry,
                         debits, _lower_contributions(ms, AbstractSubModel[]), rebuilds,
                         growth, dilute, n_idx, Float64(interval), times, base, footprint,
                         corea_volume_cap_litres(volume0), area0, volume0, radius0)

    slots = Tuple(_Slot(_svec(ctx.state_idxs), _svec(ctx.param_idxs),
                        _svec(Int[ctx.input_map[s] for s in inputs(m)]),
                        _svec(ctx.contrib_idxs)) for (m, ctx) in zip(ms, contexts))
    subs = Tuple(map(enumerate(ms)) do (i, m)
        pidx = contexts[i].param_idxs
        local_k(c) = findfirst(==(c.param_idx), pidx)
        cs = [c for c in catalytic if c.declared_by === module_id(m)]
        gs = [g for g in geometry if g.declared_by === module_id(m)]
        _RouteSubs(_svec([local_k(c) for c in cs]), _svec([c.count_idx for c in cs]),
                   _svec([local_k(g) for g in gs]), [g.quantity for g in gs])
    end)
    live = [b for b in debits if b.pool_idx != 0]
    coef(b) = b.sign < 0 ? -1.0 : b.stoich
    f = _RouteRHS(Tuple(ms), slots, subs, rt0, [b.counter_idx for b in live],
                  [b.pool_idx for b in live], [coef(b) for b in live],
                  [g.idx for g in growth])
    debited = unique(b.pool_idx for b in live if b.sign < 0)
    prob = ODEProblem{true}(f, u0, tspan, p0; callback = _route_callbacks(rt0, debited),
                            tstops = times)
    rt = MeanFieldRoute(prob, ms, names, rt0.param_names, catalytic, geometry, debits,
                        rt0.contributions, rebuilds, growth, dilute, n_idx,
                        rt0.interval, times, base, footprint, rt0.volume_cap_litres,
                        area0, volume0, radius0)
    assert_edges_executed(ms, rt)
    return rt
end

"""
    mean_field_written_params(route) -> Vector{NamedTuple}

The parameter slots route (ii) writes rather than reads: the catalytic channel's,
the inbound geometry's and the rebuild's, as [`driver_written_params`](@ref)
lists them for the hybrid. A posterior for one of them is its prior.
"""
mean_field_written_params(rt::MeanFieldRoute) = vcat(
    [(param_slot = c.param_slot, channel = :catalytic, declared_by = c.declared_by)
     for c in rt.catalytic],
    [(param_slot = g.param_slot, channel = :volume, declared_by = g.declared_by)
     for g in rt.geometry],
    [(param_slot = s, channel = :rate_constant, declared_by = r.declared_by)
     for r in rt.rebuilds for s in r.names])

"""
    mean_field_parameters(route, values; p = route.prob.p) -> Vector

A copy of `p` with `name => value` pairs written in, promoted to the values'
type. A name the route writes is refused, as `set_parameters!` refuses it.
"""
function mean_field_parameters(rt::MeanFieldRoute, values; p = rt.prob.p)
    written = Set(w.param_slot for w in mean_field_written_params(rt))
    T = promote_type(eltype(p), (typeof(v) for (_, v) in values)...)
    out = Vector{T}(p)
    for (name, v) in values
        name in written && throw(ArgumentError(
            ":$name is a slot route (ii) writes, so a value set here would not be used"))
        i = findfirst(==(name), rt.param_names)
        i === nothing && throw(ArgumentError(
            ":$name is not a free parameter of this composition"))
        out[i] = v
    end
    return out
end

"""
    solve_mean_field(route; p = route.prob.p, saveat = [], solver = Rodas5P(),
                     abstol = 1e-10, reltol = 1e-8, maxiters = 10^5)

Solve route (ii) at parameters `p`, a copy of which the rebuild writes into, at
the hybrid's pinned tolerances. The initial state is promoted to `p`'s type, so a
dual-number `p` differentiates the solve. `maxiters` caps the solver's steps, and
a solve that reaches it returns rather than throws.
"""
function solve_mean_field(rt::MeanFieldRoute; p = rt.prob.p, saveat = Float64[],
                          solver = Rodas5P(), abstol = 1e-10, reltol = 1e-8,
                          maxiters::Integer = 10^5)
    q = copy(p)
    u0 = convert(Vector{promote_type(eltype(q), eltype(rt.prob.u0))}, rt.prob.u0)
    prob = remake(rt.prob; u0, p = q)
    return solve(prob, solver; saveat, abstol, reltol, maxiters, save_everystep = false)
end

"""
    route_declarations(route) -> Vector{ReductionLabel}

Route (ii)'s own departures from the reference composition (§12 2026-10-23),
beside the modules' [`reduction_declarations`](@ref).
"""
function route_declarations(rt::MeanFieldRoute)
    labels = [
        ReductionLabel(:formalism, :continuous_drain,
            "every deferred counter is drained continuously, unclamped, at its mean " *
            "accrual rate, where the reference debits whole particles a handshake " *
            "later and clips at zero. A pool the drains overdraw stops the solve"),
        ReductionLabel(:formalism, :rebuild_snapshot,
            "the rate constants switch at k·$(Int(round(maximum((r.interval for r in rt.rebuilds); init = 0.0)))) − " *
            "$(rt.interval) s from the pools at that instant; the reference's ODE " *
            "runs one handshake ahead and reads them $(rt.interval) s later"),
        ReductionLabel(:formalism, :unrounded_counts,
            "the catalytic channel reads mean protein counts, not whole particles"),
    ]
    if !isempty(rt.growth)
        push!(labels, ReductionLabel(:formalism, :continuous_growth,
            "the volume follows the membrane count continuously, where the reference " *
            "recomputes it at each $(rt.interval) s handshake"))
        push!(labels, ReductionLabel(:exogenous_growth, :membrane_area_baseline,
            "$(round(rt.area_baseline_nm2; digits = 1)) nm² of the initial " *
            "$(round(rt.initial_area_nm2; digits = 1)) nm² surface area is held " *
            "frozen, as in the hybrid"))
    end
    return labels
end

reduction_declarations(models::Vector{<:AbstractSubModel}, rt::MeanFieldRoute) =
    vcat(reduction_declarations(models), route_declarations(rt))

export mean_field_swap, MeanFieldRoute, build_mean_field_problem,
       mean_field_written_params, mean_field_parameters, solve_mean_field,
       route_declarations
