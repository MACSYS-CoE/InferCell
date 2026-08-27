using InferCell

# `using` brings these into scope but does not permit extending them; adding a
# method to another module's function needs an explicit import.
import InferCell: states, parameters, coupling, inputs, reduction_notes,
                  inference_mode, module_id, formalism

"""
    CoreAStub(id; st, edges, ins, params, notes, mode)

A minimal `AbstractSubModel` for exercising the Core A′ interface contract.

The suite has no other sub-model doubles — every existing test drives one of the
five real models — so this one exists to let the resolver be tested on
compositions that do not exist yet. It carries an explicit `id` rather than
relying on the type name, because most of these tests compose several stubs and
need to tell them apart in an error message.

It implements no dynamics: nothing here builds a problem or integrates.
"""
struct CoreAStub <: AbstractSubModel
    id::Symbol
    st::Vector{Symbol}
    edges::Vector{CouplingEdge}
    ins::Vector{Symbol}
    params::Vector{InferParameter}
    notes::Vector{String}
    mode::Symbol
    form::Symbol
end

function CoreAStub(id::Symbol;
                   st = Symbol[],
                   edges = CouplingEdge[],
                   ins = Symbol[],
                   params = InferParameter[],
                   notes = String[],
                   mode = :simulation,
                   form = :ode)
    return CoreAStub(id,
                     collect(Symbol, st),
                     collect(CouplingEdge, edges),
                     collect(Symbol, ins),
                     collect(InferParameter, params),
                     collect(String, notes),
                     mode,
                     form)
end

module_id(m::CoreAStub) = m.id
states(m::CoreAStub) = m.st
parameters(m::CoreAStub) = m.params
coupling(m::CoreAStub) = m.edges
inputs(m::CoreAStub) = m.ins
reduction_notes(m::CoreAStub) = m.notes
inference_mode(m::CoreAStub) = m.mode
formalism(m::CoreAStub) = m.form
