using InferCell

# `using` brings these into scope but does not permit extending them; adding a
# method to another module's function needs an explicit import.
import InferCell: states, parameters, coupling, inputs, reduction_notes,
                  inference_mode, module_id, formalism, contributed_states,
                  membrane_protein_states

"""
    CoreAStub(id; st, edges, ins, contribs, params, notes, mode, form, membrane)

A minimal `AbstractSubModel` for exercising the Core A′ interface contract.

The suite has no other sub-model doubles — every existing test drives one of the
five real models — so this one exists to let the resolver be tested on
compositions that do not exist yet. It carries an explicit `id` rather than
relying on the type name, because most of these tests compose several stubs and
need to tell them apart in an error message.

It implements no dynamics: nothing here builds a problem or integrates. It does
declare `contributed_states`, so the resolver's contribution checks can be
exercised on it; the executable doubles are in `contribution_test_models.jl`.
"""
struct CoreAStub <: AbstractSubModel
    id::Symbol
    st::Vector{Symbol}
    edges::Vector{CouplingEdge}
    ins::Vector{Symbol}
    contribs::Vector{Symbol}
    params::Vector{InferParameter}
    notes::Vector{Any}
    mode::Symbol
    form::Symbol
    membrane::Vector{Symbol}
end

function CoreAStub(id::Symbol;
                   st = Symbol[],
                   edges = CouplingEdge[],
                   ins = Symbol[],
                   contribs = Symbol[],
                   params = InferParameter[],
                   notes = String[],
                   mode = :simulation,
                   form = :ode,
                   membrane = Symbol[])
    return CoreAStub(id,
                     collect(Symbol, st),
                     collect(CouplingEdge, edges),
                     collect(Symbol, ins),
                     collect(Symbol, contribs),
                     collect(InferParameter, params),
                     collect(Any, notes),
                     mode,
                     form,
                     collect(Symbol, membrane))
end

# The suite's error-capture idiom: run `f`, return the exception it throws so
# its message can be inspected, or `nothing` when it doesn't throw.
caught(f) = try f(); nothing catch e; e end

module_id(m::CoreAStub) = m.id
states(m::CoreAStub) = m.st
parameters(m::CoreAStub) = m.params
coupling(m::CoreAStub) = m.edges
inputs(m::CoreAStub) = m.ins
contributed_states(m::CoreAStub) = m.contribs
reduction_notes(m::CoreAStub) = m.notes
inference_mode(m::CoreAStub) = m.mode
formalism(m::CoreAStub) = m.form
membrane_protein_states(m::CoreAStub) = m.membrane
