"""
Mean-field gene expression: the second description of track B's formalism swap
(spec/phases/approximate-method.md, task B2).

A group of `:jump` modules whose propensities are affine in the state, so the
network is first order, has mean equations that close exactly:

```
d E[x] / dt = Σ_r a_r(E[x]) · ν_r
```

where `a_r` is reaction r's propensity and `ν_r` its state change. Given the
rate constants these means are exact, not approximate. The approximation is in
what reads them: a nonlinear module fed the mean sees the Jensen gap.

[`MeanField`](@ref) builds that ODE from the jump modules themselves. It reads
each reaction's propensity closure and gets `ν_r` by applying the reaction's
own `affect!` to zeros, so the two descriptions cannot drift apart. The group
is one `:ode` module because a jump module may write a peer's non-registry
state (decay writes transcription's transcripts), and the ODE path's only
write channel, [`contributions`](@ref), takes registry species alone.
"""

"""
    MeanField(members; edges = :declared, id = :MeanField)

One `:ode` module carrying the mean equations of the `:jump` modules
`members`. It has the same states, in member order, and the same parameters.
A member's counters become cumulative expected drains.

Every member input must be another member's state, so the group is closed.
Each reaction's state change must not depend on the state, and each
propensity must be affine in it. Both are checked here, since the means are
exact only for such a network.

`edges = :declared` carries the members' coupling edges and rebuilt
parameters over unchanged. `:none` drops them, for a standalone build at
pinned rate constants. Route (ii), [`build_mean_field_problem`](@ref), executes
them; the plain [`build_problem`](@ref) does not.
"""
struct MeanField{N, K <: Tuple} <: AbstractSubModel
    members::Vector{AbstractSubModel}
    names::Vector{Symbol}
    params::Vector{InferParameter}
    pidx::Vector{Vector{Int}}
    edges::Vector{CouplingEdge}
    rebuilt::Vector{Symbol}
    id::Symbol
    kernel::K
end

function MeanField(members::AbstractVector{<:AbstractSubModel};
                   edges::Symbol = :declared, id::Symbol = :MeanField)
    edges in (:declared, :none) || throw(ArgumentError(
        "edges must be :declared or :none, got :$edges"))
    members = collect(AbstractSubModel, members)
    for mem in members
        formalism(mem) === :jump || throw(ArgumentError(
            "MeanField takes :jump modules; $(module_id(mem)) is :$(formalism(mem))"))
    end
    names = reduce(vcat, (states(mem) for mem in members); init = Symbol[])
    allunique(names) || throw(ArgumentError("MeanField's members share a state name"))
    params = reduce(vcat, (parameters(mem) for mem in members); init = InferParameter[])
    allunique(p.name for p in params) || throw(ArgumentError(
        "MeanField's members share a parameter name, so their slots would collide"))
    pos = Dict(s => i for (i, s) in enumerate(names))

    sidx, pidx, iidx = Vector{Int}[], Vector{Int}[], Vector{Int}[]
    soff, poff = 0, 0
    for mem in members
        ns, np = length(states(mem)), length(model_free_params(parameters(mem)))
        push!(sidx, collect(soff .+ (1:ns)))
        push!(pidx, collect(poff .+ (1:np)))
        soff += ns
        poff += np
        for s in inputs(mem)
            haskey(pos, s) || throw(ArgumentError(
                "$(module_id(mem)) reads :$s, which no member of the group owns. " *
                "MeanField resolves every input inside the group"))
        end
        push!(iidx, [pos[s] for s in inputs(mem)])
    end

    N = length(names)
    rates, member, cols = Any[], Int[], Vector{Float64}[]
    for (j, mem) in enumerate(members)
        p = [q.value for q in model_free_params(parameters(mem))]
        writable = Set(findall(in(written_states(mem)), inputs(mem)))
        for r in reactions(mem)
            col = zeros(N)
            Δu, Δw = _state_change(r, length(sidx[j]), length(iidx[j]), 0.0)
            Δu2, Δw2 = _state_change(r, length(sidx[j]), length(iidx[j]), 7.0)
            Δu == Δu2 && Δw == Δw2 || throw(ArgumentError(
                "A reaction of $(module_id(mem)) changes the state by an amount " *
                "that depends on the state, so its mean equation does not close"))
            all(i in writable || Δw[i] == 0 for i in eachindex(Δw)) || throw(ArgumentError(
                "A reaction of $(module_id(mem)) writes an input it does not list " *
                "in written_states"))
            col[sidx[j]] .+= Δu
            col[iidx[j]] .+= Δw
            # At the declared values and at all ones, so a zero placeholder
            # cannot hide a nonlinear term.
            all(q -> _check_affine(r, q, length(sidx[j]), length(iidx[j])), (p, ones(length(p)))) ||
                throw(ArgumentError(
                    "A propensity of $(module_id(mem)) is not affine in the state, so " *
                    "the mean equations do not close and MeanField would not be exact"))
            push!(rates, r.rate)
            push!(member, j)
            push!(cols, col)
        end
    end
    ν = isempty(cols) ? zeros(N, 0) : reduce(hcat, cols)
    nz = [findall(!iszero, view(ν, :, r)) for r in axes(ν, 2)]
    declared = edges === :declared
    es = declared ? reduce(vcat, (coupling(mem) for mem in members); init = CouplingEdge[]) :
         CouplingEdge[]
    rebuilt = declared ? reduce(vcat, (rebuilt_params(mem) for mem in members); init = Symbol[]) :
              Symbol[]
    kernel = Tuple(_kernel_term(rates[r], sidx[member[r]], pidx[member[r]], iidx[member[r]],
                                nz[r], ν[nz[r], r]) for r in eachindex(rates))
    return MeanField{N, typeof(kernel)}(members, names, params, pidx, collect(CouplingEdge, es),
                                        rebuilt, id, kernel)
end

# One reaction's term in `dynamics`, with its closure and indices concretely
# typed, so the right-hand side neither dispatches nor allocates per call.
_kernel_term(rate, s, p, i, nz, ν) =
    (; rate, s = SVector{length(s)}(s), p = SVector{length(p)}(p), i = SVector{length(i)}(i),
     nz = SVector{length(nz)}(nz), ν = SVector{length(ν)}(ν))

# A reaction's state change, by applying its affect to a base state of `x`.
function _state_change(r, ns, ni, x)
    u, w = fill(x, ns), fill(x, ni)
    r.affect!(u, w)
    return u .- x, w .- x
end

# Affine in the state: a(u₁) + a(u₂) = a(u₁ + u₂) + a(0), at two fixed points.
function _check_affine(r, p, ns, ni)
    u1, w1 = collect(1.0:ns) .* 0.7, collect(1.0:ni) .* 1.3
    u2, w2 = reverse(u1) .+ 2.0, reverse(w1) .+ 0.5
    a(u, w) = r.rate(u, p, 0.0, w)
    lhs = a(u1, w1) + a(u2, w2)
    rhs = a(u1 + u2, w1 + w2) + a(zeros(ns), zeros(ni))
    return isapprox(lhs, rhs; rtol = 1e-12, atol = 1e-300)
end

states(m::MeanField) = m.names
parameters(m::MeanField) = m.params
formalism(::MeanField) = :ode
inference_mode(::MeanField) = :differentiable
coupling(m::MeanField) = m.edges
rebuilt_params(m::MeanField) = m.rebuilt
module_id(m::MeanField) = m.id

"""
    dynamics(u, p, t, m::MeanField)

`Σ_r a_r(u) · ν_r`, each propensity read from its member's own closure at the
member's slices of the mean state and the group's parameters.
"""
function dynamics(u, p, t, m::MeanField{N}) where {N}
    T = promote_type(eltype(u), eltype(p))
    du = zeros(MVector{N, T})
    _add_terms!(du, u, p, t, m.kernel)
    return SVector{N, T}(du)
end

# The reactions in order, each adding a_r · ν_r over ν_r's nonzeros.
@inline _add_terms!(du, u, p, t, ::Tuple{}) = du
@inline function _add_terms!(du, u, p, t, ks::Tuple)
    k = first(ks)
    a = k.rate(u[k.s], p[k.p], t, u[k.i])
    for j in eachindex(k.nz)
        du[k.nz[j]] += a * k.ν[j]
    end
    return _add_terms!(du, u, p, t, Base.tail(ks))
end

function reduction_notes(m::MeanField)
    ids = join(string.(module_id.(m.members)), ", ")
    own = Any[
        :formalism => "$(ids) run as their mean equations ($(m.id)). The network is " *
            "first order, so given the rate constants the means are exact. What is " *
            "dropped is cell-to-cell variation: a nonlinear module fed these means " *
            "sees the Jensen gap (spec/phases/approximate-method.md §3).",
        :coarse_drain => "$(m.id)'s counters accrue the expected drain continuously, " *
            "not in whole events.",
    ]
    return vcat(own, reduce(vcat, (reduction_notes(mem) for mem in m.members); init = Any[]))
end

"""
    gene_expression_mean_field(; genes = read_transcription_genes(), edges = :declared, kwargs...)

The reduced syn3A model's gene expression as one [`MeanField`](@ref) module:
transcription, decay and translation, built with `kwargs` (for example the
drain's `clip`) and cut to `genes`, as `m0_models` cuts them.
"""
function gene_expression_mean_field(; genes = read_transcription_genes(),
                                    edges::Symbol = :declared, kwargs...)
    genes = collect(TranscriptionGene, genes)
    return MeanField(AbstractSubModel[CoreATranscription(; genes, kwargs...),
                                      CoreATranscriptDecay(; genes, kwargs...),
                                      CoreATranslation(; genes, kwargs...)];
                     edges, id = :GeneExpressionMeanField)
end

export MeanField, gene_expression_mean_field
