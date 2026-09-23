using Test
using InferCell
using StaticArrays: SA, SVector

import InferCell: states, parameters, reactions, formalism, inference_mode,
                  coupling

# Spec §11 phase 11a: a catalytic edge may name a protein count the jump block
# owns, the two metabolic modules whose enzymes translation makes expose those
# enzymes as slots, and a producer counter can raise a PTS carrier. Uses the
# phase 3 and phase 5 doubles (`ToyMembranePool`, `caught`, `_jic`, `_rate`),
# so it is included after them.

"""
    ToyProteome(counts; k = 0.0, counter = nothing, credits = nothing)

A jump double owning one protein count per entry of `counts` (`name => initial
particles`), each made at rate `k`. With `counter` and `credits` set it also
owns that accrual counter, bumped by one per protein made, and declares a
producer `DeferredCounterEdge` crediting `credits` — the shape translation will
use to deliver a PTS carrier (spec §12, 2026-09-23 C).
"""
struct ToyProteome <: AbstractSubModel
    names::Vector{Symbol}
    k::Float64
    params::Vector{InferParameter}
    counter::Union{Symbol, Nothing}
    edges::Vector{CouplingEdge}
end

function ToyProteome(counts; k = 0.0, counter = nothing, credits = nothing)
    names = [first(c) for c in counts]
    # The rate is held on the struct rather than in the parameter vector: this
    # double exists to own counts, not to be inferred.
    params = [_jic(last(c), Symbol(first(c), "0"), :ToyProteome) for c in counts]
    edges = CouplingEdge[]
    if counter !== nothing
        push!(params, _jic(0, Symbol(counter, "0"), :ToyProteome))
        push!(edges, DeferredCounterEdge(species = credits, direction = :out,
                                         counter = counter))
    end
    return ToyProteome(names, Float64(k), params, counter, edges)
end

states(m::ToyProteome) = m.counter === nothing ? copy(m.names) : vcat(m.names, m.counter)
parameters(m::ToyProteome) = m.params
formalism(::ToyProteome) = :jump
inference_mode(::ToyProteome) = :simulation
coupling(m::ToyProteome) = m.edges
function reactions(m::ToyProteome)
    n = length(m.names)
    bump = m.counter !== nothing
    k = m.k
    return [Reaction((u, p, t, _) -> k,
                     (u, _) -> (u[i] += 1; bump && (u[n + 1] += 1)))
            for i in 1:n]
end

# The registry's initial concentration for a species, which is what a module's
# own initial-condition parameter holds.
_registry_conc(s) = something(species_entry(s).initial_value, 0.1)

@testset "Phase 11a: catalytic edges on jump-owned protein counts" begin

    @testset "11a.1 a catalytic edge may name a non-registry count" begin
        @test !is_registered(:P_toy)

        # No membrane flag: this testset is about the catalytic channel alone, and
        # a flag with no outbound volume edge is refused at build.
        pool = ToyMembranePool(membrane = Symbol[], edges = CouplingEdge[
            CatalyticEdge(species = :P_toy, direction = :in, param_slot = :enzyme_conc)])
        graph = resolve_coupling(pool)
        r = only(graph.edges)
        @test r.kind === :catalytic && r.species === :P_toy
        @test r.state_index == 0
        @test isempty(graph.dead_ends)

        # Every other kind on the same name is still refused, naming it.
        others = CouplingEdge[
            MassEdge(species = :P_toy, direction = :in),
            CurrencyEdge(species = :P_toy, direction = :in),
            DeferredCounterEdge(species = :P_toy, direction = :in, counter = :c_toy),
            RateConstantEdge(species = :P_toy, direction = :in),
            VolumeEdge(species = :P_toy, direction = :out),
            ClampedEdge(species = :P_toy, direction = :in, held_value = 1.0,
                        origin = :ours),
        ]
        for e in others
            err = caught(() -> resolve_coupling(CoreAStub(:Stub; edges = [e])))
            @test err isa ArgumentError
            @test occursin("P_toy", err.msg) && occursin("Only a CatalyticEdge", err.msg)
        end

        # A catalytic edge on a count no jump module owns is refused at build.
        orphan = ToyMembranePool(membrane = Symbol[], edges = CouplingEdge[
            CatalyticEdge(species = :P_missing, direction = :in, param_slot = :enzyme_conc)])
        err = caught(() -> build_problem([orphan, ToyProteome([:P_toy => 400])]))
        @test err isa ArgumentError
        @test occursin("P_missing", err.msg)

        # And it executes: the count fills the slot at every handshake.
        d = build_problem([pool, ToyProteome([:P_toy => 400]; k = 2.0)];
                          tspan = (0.0, 30.0))
        c = only(d.catalytic)
        @test c.species === :P_toy && c.param_slot === :enzyme_conc
        for _ in 1:5
            handshake_step!(d)
            # The write precedes the jump step, so the slot holds the count as it
            # stood when this handshake began, which the next step then reads.
        end
        before = d.jump.u[c.count_idx]
        handshake_step!(d)
        @test d.ode.p[c.param_idx] == counts_to_mM(before, d.factor)
        @test before > 400                     # the count really moved
    end

    @testset "11a.2 glycolysis's translated-enzyme mode" begin
        nominal = CentralGlycolysis()
        tr = CentralGlycolysis(enzymes = :translated)
        ids = [r.id for r in GLYCOLYTIC_REACTIONS]
        slots = [Symbol(:enz_, id) for id in ids]

        # The default mode is phase 6's: no slot, no catalytic edge, counts in inputs.
        @test !any(p -> p.name in slots, parameters(nominal))
        @test !any(e -> e isa CatalyticEdge, coupling(nominal))
        @test InferCell.inputs(nominal) == vcat(GLYCOLYSIS_CURRENCIES, default_protein_sources())
        @test isempty(model_free_params(parameters(nominal)))

        # Translated: ten free slots at the nominal values, asserted.
        free = model_free_params(parameters(tr))
        @test [q.name for q in free] == slots
        @test [q.value for q in free] == collect(nominal.enzyme_conc)
        @test all(q -> q.provenance.informedness === :asserted, free)
        cat = [e for e in coupling(tr) if e isa CatalyticEdge]
        @test [e.species for e in cat] == default_protein_sources()
        @test [e.param_slot for e in cat] == slots
        @test coupling(tr)[1:9] == coupling(nominal)
        @test InferCell.inputs(tr) == collect(GLYCOLYSIS_CURRENCIES)
        @test resolve_coupling(tr) isa CouplingGraph

        # The derivative is bitwise the default mode's at the registry state.
        u = SVector{13}(_registry_conc.(glycolytic_states()))
        cur = [_registry_conc(s) for s in GLYCOLYSIS_CURRENCIES]
        uin_nominal = vcat(cur, fill(1.0, 10))
        p_tr = [q.value for q in free]
        @test InferCell.dynamics(u, Float64[], 0.0, nominal, uin_nominal) ===
              InferCell.dynamics(u, p_tr, 0.0, tr, cur)
        @test InferCell.contributions(u, Float64[], 0.0, nominal, uin_nominal) ===
              InferCell.contributions(u, p_tr, 0.0, tr, cur)

        # And it reads the slot: doubling one enzyme doubles that reaction alone.
        p2 = copy(p_tr); p2[5] *= 2
        v1 = reaction_rates(u, p_tr, 0.0, tr, cur)
        v2 = reaction_rates(u, p2, 0.0, tr, cur)
        @test v2[5] == 2 * v1[5]
        @test v2[[1:4; 6:10]] == v1[[1:4; 6:10]]

        # The notes say which mode it is.
        @test any(occursin("nominal stand-ins", n) for n in InferCell.reduction_notes(nominal))
        @test !any(occursin("nominal stand-ins", n) for n in InferCell.reduction_notes(tr))
        @test any(occursin("catalytic edges overwrite", n) for n in InferCell.reduction_notes(tr))

        @test caught(() -> CentralGlycolysis(enzymes = :live)) isa ArgumentError
    end

    @testset "11a.3 recycling's translated-enzyme mode, and one count per gene" begin
        nominal = NucleotideRecycling()
        tr = NucleotideRecycling(enzymes = :translated)
        @test !any(e -> e isa CatalyticEdge, coupling(nominal))
        @test [q.name for q in model_free_params(parameters(tr))] ==
              collect(RECYCLING_ENZYME_IDS)
        cat = [e for e in coupling(tr) if e isa CatalyticEdge]
        @test [e.param_slot for e in cat] == collect(RECYCLING_ENZYME_IDS)
        @test [e.species for e in cat] ==
              [Symbol(:P_, e.locus) for e in recycling_enzymes()]
        @test caught(() -> NucleotideRecycling(enzymes = :live)) isa ArgumentError

        # PGK3 and PYK3 read the counts glycolysis's PGK and PYK read.
        glyc = Dict(e.param_slot => e.species
                    for e in coupling(CentralGlycolysis(enzymes = :translated))
                    if e isa CatalyticEdge)
        rec = Dict(e.param_slot => e.species for e in cat)
        @test rec[:enz_R_PGK3] === glyc[:enz_R_PGK]
        @test rec[:enz_R_PYK3] === glyc[:enz_R_PYK]

        # Composed: thirteen counts fill fifteen slots, each from its own count.
        loci = unique(vcat(collect(values(glyc)), collect(values(rec))))
        @test length(loci) == 13
        copies = Dict(Symbol(:P_, r.locus) => r.copies for r in GLYCOLYTIC_REACTIONS)
        for e in recycling_enzymes()
            copies[Symbol(:P_, e.locus)] = e.copies
        end
        prot = ToyProteome([s => copies[s] for s in loci])
        d = build_problem([CentralGlycolysis(enzymes = :translated), tr, prot];
                          tspan = (0.0, 10.0))
        @test length(d.catalytic) == 15
        handshake_step!(d)
        for c in d.catalytic
            @test d.ode.p[c.param_idx] == counts_to_mM(d.jump.u[c.count_idx], d.factor)
        end
        slot(name) = d.ode.p[only(c for c in d.catalytic if c.param_slot === name).param_idx]
        @test slot(:enz_R_PGK3) == slot(:enz_R_PGK)
        @test slot(:enz_R_PYK3) == slot(:enz_R_PYK)
        # At the published copy numbers the live value is the nominal one.
        @test slot(:enz_R_ADK1) ≈ recycling_enzymes()[3].concentration rtol = 1e-12
    end

    @testset "11a.4 a producer counter raises a PTS carrier, and the cell sees it" begin
        pool = ToyMembranePool(edges = CouplingEdge[
            CatalyticEdge(species = :P_toy, direction = :in, param_slot = :enzyme_conc),
            VolumeEdge(species = :M_ptsi_c, direction = :out)])
        prot = ToyProteome([:P_toy => 0]; k = 5.0, counter = :ptsi_made,
                           credits = :M_ptsi_c)
        d = build_problem([pool, prot]; tspan = (0.0, 200.0))
        i = only(k for (k, s) in enumerate(states(pool)) if s === :M_ptsi_c)
        a0 = growth_census(d).area_nm2
        carriers0 = d.ode.u[i] * d.factor
        for _ in 1:100
            before = d.ode.u[i] * d.factor
            accrued = d.jump.u[end]
            handshake_step!(d)
            # Exactly the accrual, in particles, at every drain.
            @test d.ode.u[i] * d.factor ≈ before + accrued rtol = 1e-12
        end
        # Every protein made is a carrier credited or still pending.
        made = d.jump.u[1]
        @test made > 300
        @test d.ode.u[i] * d.factor + d.jump.u[end] ≈ carriers0 + made rtol = 1e-12

        # The volume chain reads the carrier through its owner's membrane flag.
        g = growth_census(d)
        @test only(g.states).species === :M_ptsi_c
        @test g.n_membrane ≈ d.ode.u[i] * d.factor rtol = 1e-12
        @test g.area_nm2 > a0
    end
end
