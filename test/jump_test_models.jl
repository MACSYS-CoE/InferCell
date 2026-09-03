using InferCell
using Distributions: LogNormal, Normal
using JumpProcesses: ConstantRateJump

import InferCell: states, parameters, reactions, formalism, inputs, written_states,
                  module_id

# Executable doubles for jump composition (spec §11 phase 2).
#
# The two shipped jump models both own :mRNA and :protein, so composing them
# fails on duplicate ownership before any index machinery runs. These doubles
# use disjoint, non-registry state names so that what is being tested is the
# composition itself and not the resolver's registry checks.

_jrate(k, name, id) = InferParameter(k, LogNormal(log(max(k, 1e-12)), 0.5), false, name, id, :rate)
_jic(x0, name, id) = InferParameter(Float64(x0), Normal(Float64(x0), 1.0), true, name, id, :initial_condition)

# ---------------------------------------------------------------------------
# Legacy-style doubles: reactions() returns ConstantRateJumps whose closures
# index the *global* vectors at local positions, exactly as the shipped models
# did before phase 2. They exist to witness the aliasing bug (task 2.1).
# ---------------------------------------------------------------------------

"""
    LegacyOwner(; x0 = 3, k = nothing)

Owns `:legacy_x`, fires nothing. With `k` given it also carries one free rate
parameter `k_dead` that no reaction reads, so a composition's first global
parameter slot is this module's.
"""
struct LegacyOwner <: AbstractSubModel
    params::Vector{InferParameter}
end
function LegacyOwner(; x0 = 3, k = nothing)
    params = [_jic(x0, :legacy_x0, :LegacyOwner)]
    k === nothing || pushfirst!(params, _jrate(k, :k_dead, :LegacyOwner))
    return LegacyOwner(params)
end
states(::LegacyOwner) = [:legacy_x]
parameters(m::LegacyOwner) = m.params
formalism(::LegacyOwner) = :jump
reactions(::LegacyOwner) = ConstantRateJump[]

"""
    LegacyWriter(; b = 1.0)

Owns `:legacy_y`. One zeroth-order jump at rate `b_legacy` whose affect does
`integrator.u[1] += 1` and whose rate reads `p[1]` — local positions written
against the global vectors, which is the bug.
"""
struct LegacyWriter <: AbstractSubModel
    params::Vector{InferParameter}
end
LegacyWriter(; b = 1.0) = LegacyWriter([_jrate(b, :b_legacy, :LegacyWriter),
                                        _jic(0, :legacy_y0, :LegacyWriter)])
states(::LegacyWriter) = [:legacy_y]
parameters(m::LegacyWriter) = m.params
formalism(::LegacyWriter) = :jump
reactions(::LegacyWriter) =
    [ConstantRateJump((u, p, t) -> p[1], integrator -> (integrator.u[1] += 1))]

# ---------------------------------------------------------------------------
# Phase-2 doubles: reactions() returns Reactions in local coordinates.
# ---------------------------------------------------------------------------

"""
    BirthOwner(; k = 1.0, x0 = 0)

Owns `:X` and produces it at constant rate `k_birth`. With `k = 0` it is a
frozen pool, which is what the peer-read test needs.
"""
struct BirthOwner <: AbstractSubModel
    params::Vector{InferParameter}
end
BirthOwner(; k = 1.0, x0 = 0) = BirthOwner([_jrate(k, :k_birth, :BirthOwner),
                                            _jic(x0, :X0, :BirthOwner)])
states(::BirthOwner) = [:X]
parameters(m::BirthOwner) = m.params
formalism(::BirthOwner) = :jump
reactions(::BirthOwner) = [Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1))]

"""
    PeerBirthDeath(; gamma = 0.5, b = 0.0, bname = :b_peer, writes = [:X])

Owns a firing counter `:fired` and neither owns nor produces `:X` on its own
account: it reads `:X` as a peer state, removes one at rate `gamma_peer · X`
and adds one at rate `b`, counting every firing. `bname` lets a test share the
birth parameter's name with `BirthOwner` so parameter deduplication is
exercised; `writes` is what the write gate is held to.
"""
struct PeerBirthDeath <: AbstractSubModel
    params::Vector{InferParameter}
    writes::Vector{Symbol}
end
PeerBirthDeath(; gamma = 0.5, b = 0.0, bname = :b_peer, writes = [:X]) =
    PeerBirthDeath([_jrate(gamma, :gamma_peer, :PeerBirthDeath),
                    _jrate(b, bname, :PeerBirthDeath),
                    _jic(0, :fired0, :PeerBirthDeath)],
                   collect(Symbol, writes))
states(::PeerBirthDeath) = [:fired]
parameters(m::PeerBirthDeath) = m.params
formalism(::PeerBirthDeath) = :jump
inputs(::PeerBirthDeath) = [:X]
written_states(m::PeerBirthDeath) = m.writes
reactions(::PeerBirthDeath) = [
    # X -> ∅ at gamma · X, through the peer view
    Reaction((u, p, t, uin) -> p[1] * uin[1], (u, uin) -> (uin[1] -= 1; u[1] += 1)),
    # ∅ -> X at b, through the peer view
    Reaction((u, p, t, uin) -> p[2], (u, uin) -> (uin[1] += 1; u[1] += 1)),
]

"""
    PeerReader(; c = 1.0)

Owns `:fired`; its one reaction fires at `c_read · X` and increments only its
own counter. It reads a peer and writes nothing of the peer's.
"""
struct PeerReader <: AbstractSubModel
    params::Vector{InferParameter}
end
PeerReader(; c = 1.0) = PeerReader([_jrate(c, :c_read, :PeerReader), _jic(0, :fired0, :PeerReader)])
states(::PeerReader) = [:fired]
parameters(m::PeerReader) = m.params
formalism(::PeerReader) = :jump
inputs(::PeerReader) = [:X]
reactions(::PeerReader) = [Reaction((u, p, t, uin) -> p[1] * uin[1], (u, _) -> (u[1] += 1))]

"""
    BirthDeath(; k = 1.0, gamma = 0.5, x0 = 0)

The hand-written single-module equivalent of `BirthOwner` composed with
`PeerBirthDeath`: `∅ -> X` at `k_bd`, `X -> ∅` at `gamma_bd · X`.
"""
struct BirthDeath <: AbstractSubModel
    params::Vector{InferParameter}
end
BirthDeath(; k = 1.0, gamma = 0.5, x0 = 0) =
    BirthDeath([_jrate(k, :k_bd, :BirthDeath), _jrate(gamma, :gamma_bd, :BirthDeath),
                _jic(x0, :X0, :BirthDeath)])
states(::BirthDeath) = [:X]
parameters(m::BirthDeath) = m.params
formalism(::BirthDeath) = :jump
reactions(::BirthDeath) = [
    Reaction((u, p, t, _) -> p[1], (u, _) -> (u[1] += 1)),
    Reaction((u, p, t, _) -> p[2] * u[1], (u, _) -> (u[1] -= 1)),
]

# ---------------------------------------------------------------------------
# Instrumentation for task 2.2: a vector that records every index read and
# written through it, used as the parent the orchestrator's views wrap.
# ---------------------------------------------------------------------------

struct TrackedVector{T} <: AbstractVector{T}
    data::Vector{T}
    reads::Vector{Int}
    writes::Vector{Int}
end
TrackedVector(data::Vector) = TrackedVector(data, Int[], Int[])
Base.size(v::TrackedVector) = size(v.data)
Base.IndexStyle(::Type{<:TrackedVector}) = IndexLinear()
Base.getindex(v::TrackedVector, i::Int) = (push!(v.reads, i); v.data[i])
Base.setindex!(v::TrackedVector, x, i::Int) = (push!(v.writes, i); v.data[i] = x)

# What the orchestrator's affect wrapper reads from a JumpProcesses integrator.
struct FakeIntegrator{U}
    u::U
end
