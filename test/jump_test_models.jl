using InferCell
using Distributions: LogNormal, Normal
using JumpProcesses: ConstantRateJump

import InferCell: states, parameters, reactions, formalism, inputs, module_id

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
