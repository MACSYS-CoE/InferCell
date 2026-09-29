using InferCell

# The frozen-expression variant of check 8 moved to
# `src/organisms/coreA/control_variant.jl` in task 15.5. What stays here is a
# test double only.

"""
    AssembledTwoAtp(inner = TwoAtpCharging())

Phase 9's two-ATP lumping (`test/trna_test_models.jl`) with the owner's side of
translation's two tRNA counters that `TrnaCharging` declares (spec §11 task
13.1). Without them the assembled model's completeness check refuses the
composition, since the counters would debit a pool with no owner edge. Task
14.9 reruns task 9.9's comparison on the assembled model with it.
"""
struct AssembledTwoAtp{M} <: AbstractSubModel
    inner::M
end
for f in (:states, :parameters, :module_id, :inputs, :contributed_states)
    @eval InferCell.$f(m::AssembledTwoAtp) = InferCell.$f(m.inner)
end
InferCell.coupling(m::AssembledTwoAtp) = vcat(coupling(m.inner), CouplingEdge[
    CurrencyEdge(species = :M_trna_c, direction = :in),
    CurrencyEdge(species = :M_trna_chg_c, direction = :out),
])
InferCell.dynamics(u, p, t, m::AssembledTwoAtp, ui) = dynamics(u, p, t, m.inner, ui)
InferCell.contributions(u, p, t, m::AssembledTwoAtp, ui) = contributions(u, p, t, m.inner, ui)
