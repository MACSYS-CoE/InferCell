using Test
using InferCell
using StaticArrays: SVector
using OrdinaryDiffEq: Tsit5, solve, ODEProblem
using ForwardDiff
using SciMLBase: ReturnCode

# The contribution channel: a module adds derivative terms to a registry state
# another module owns (spec §11 phase 1). The doubles are in
# contribution_test_models.jl and the CoreAStub in corea_test_models.jl.

# Everything passed as an argument and nothing captured: measured inside a
# @testset, `@allocated` would count the boxed return value of the closure.
_rhs_alloc(rhs, u, p, t) = @allocated rhs(u, p, t)

@testset "Contributions" begin
    @testset "1.3 owner's derivative is its own term plus the contribution, exactly" begin
        models = [EnergyPools(), CarbonBlock()]
        prob = build_problem(models; tspan = (0.0, 10.0))
        @test length(prob.u0) == 32
        rhs = prob.f.f
        u, p = prob.u0, prob.p
        du = rhs(u, p, 0.0)

        gamma_o, k_prod, k_cons, gamma_c = p   # global parameter order
        i_atp, i_pi, i_g6p = _gidx(models, :M_atp_c), _gidx(models, :M_pi_c), _gidx(models, :M_g6p_c)

        # The accumulator computes (own term) + (contribution) in exactly this order.
        @test du[i_atp] == (-gamma_o * u[i_atp]) + (k_prod * u[i_g6p])
        @test du[i_pi] == (-gamma_o * u[i_pi]) + (-k_cons * u[i_pi])

        # Every other owner state is untouched, and the contributor's slice is its own.
        for s in ENERGY_SPECIES
            s in (:M_atp_c, :M_pi_c) && continue
            i = _gidx(models, s)
            @test du[i] == -gamma_o * u[i]
        end
        for s in CARBON_SPECIES
            i = _gidx(models, s)
            @test du[i] == -gamma_c * u[i]
        end
    end

    @testset "done-when: the owner's state rises at the contributor's declared rate" begin
        k = 0.7
        models = [InertPool(), Source(; k = k)]
        prob = build_problem(models; tspan = (0.0, 10.0))
        i_gtp = _gidx(models, :M_gtp_c)
        @test prob.f.f(prob.u0, prob.p, 0.0)[i_gtp] == k
        sol = solve(prob, Tsit5(); saveat = 0.0:1.0:10.0)
        @test sol.retcode == ReturnCode.Success
        gtp0 = prob.u0[i_gtp]
        for (t, u) in zip(sol.t, sol.u)
            @test u[i_gtp] ≈ gtp0 + k * t rtol = 1e-12
        end
    end

    @testset "contribution length must match contributed_states" begin
        prob = build_problem([EnergyPools(), MiscountedSource()]; tspan = (0.0, 1.0))
        err = caught(() -> prob.f.f(prob.u0, prob.p, 0.0))
        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("1 term", msg) && occursin("2 species", msg)
        # A non-static container is the same protocol error.
        prob2 = build_problem([EnergyPools(), VectorSource()]; tspan = (0.0, 1.0))
        err2 = caught(() -> prob2.f.f(prob2.u0, prob2.p, 0.0))
        @test err2 isa ErrorException
        @test occursin("Vector{Float64}", sprint(showerror, err2)) && occursin("static", sprint(showerror, err2))
    end

    @testset "1.4 a contribution to a chemostatted species is rejected" begin
        m = CarbonBlock(; contribs = [:M_o2_c],
                        edges = [CurrencyEdge(species = :M_o2_c, direction = :out)], ins = Symbol[])
        err = caught(() -> resolve_coupling([EnergyPools(), m]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("M_o2_c", msg) && occursin("CarbonBlock", msg) && occursin("chemostat", msg)
        # build_problem goes through the same resolver, so it fails identically.
        @test_throws ArgumentError build_problem([EnergyPools(), m])
    end

    @testset "1.4 a contribution to a name outside the registry is rejected" begin
        m = CarbonBlock(; contribs = [:M_unobtainium_c], edges = CouplingEdge[], ins = Symbol[])
        err = caught(() -> resolve_coupling([EnergyPools(), m]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("M_unobtainium_c", msg) && occursin("CarbonBlock", msg) && occursin("registry", msg)
    end

    @testset "1.4 a contribution to a species no module owns is rejected at build" begin
        m = CarbonBlock(; ins = Symbol[], contribs = [:M_atp_c],
                        edges = [CurrencyEdge(species = :M_atp_c, direction = :out)])
        # Standalone resolution is a report, not an error: the species is merely unowned.
        graph = resolve_coupling(m)
        @test :M_atp_c in graph.unowned_states
        # Building is where an unowned target becomes an error.
        err = caught(() -> build_problem([m]))
        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("M_atp_c", msg) && occursin("CarbonBlock", msg) && occursin("own", msg)
    end

    @testset "1.4 a jump module cannot contribute yet, but may declare its edges" begin
        owner = CoreAStub(:Pools; st = [:M_atp_c], form = :jump)
        j = CoreAStub(:JumpSrc; form = :jump, contribs = [:M_atp_c],
                      edges = [MassEdge(species = :M_atp_c, direction = :out)])
        err = caught(() -> resolve_coupling([owner, j]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("JumpSrc", msg) && occursin("M_atp_c", msg) && occursin("jump", msg)
        @test_throws ArgumentError build_problem([owner, j])

        # Phase 2's shape: a jump module writing a peer's state through its
        # reactions declares the edge and no contribution, and resolves cleanly.
        j_edges_only = CoreAStub(:JumpSrc; form = :jump,
                                 edges = [MassEdge(species = :M_atp_c, direction = :out)])
        graph = resolve_coupling([owner, j_edges_only])
        @test length(graph.edges) == 1
    end

    @testset "1.5 drift: an edge on a non-owned species with no contribution throws" begin
        # Outbound edge kept, contribution deleted.
        m = CarbonBlock(; contribs = [:M_pi_c])
        err = caught(() -> resolve_coupling([EnergyPools(), m]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("M_atp_c", msg) && occursin("CarbonBlock", msg) && occursin("currency", msg)

        # Inbound edge kept, contribution deleted: consumption is a derivative
        # term too (§12 amendment of 2026-09-03).
        m2 = CarbonBlock(; contribs = [:M_atp_c])
        err2 = caught(() -> resolve_coupling([EnergyPools(), m2]))
        @test err2 isa ArgumentError
        msg2 = sprint(showerror, err2)
        @test occursin("M_pi_c", msg2) && occursin("CarbonBlock", msg2) && occursin("mass", msg2)
    end

    @testset "1.5 drift: a contribution with no edge throws" begin
        m = CarbonBlock(; contribs = [:M_atp_c, :M_pi_c, :M_gtp_c])
        err = caught(() -> resolve_coupling([EnergyPools(), m]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("M_gtp_c", msg) && occursin("CarbonBlock", msg) && occursin("edge", msg)

        # A module with contributions but no edges at all is not skipped.
        lone = CoreAStub(:Lone; contribs = [:M_atp_c])
        err2 = caught(() -> resolve_coupling([CoreAStub(:Pools; st = [:M_atp_c]), lone]))
        @test err2 isa ArgumentError
        @test occursin("Lone", sprint(showerror, err2))
    end

    @testset "1.5 an edge on a species the module itself owns needs no contribution" begin
        # The pool's owner declaring the channel from its own side (the
        # RateConstantEdge convention, and test_resolver's producer/consumer pair).
        owner = CoreAStub(:Central; st = [:M_atp_c],
                          edges = [MassEdge(species = :M_atp_c, direction = :out)])
        @test resolve_coupling(owner) isa CouplingGraph
    end

    @testset "1.5 the drift check is live through reduction_declarations too" begin
        m = CarbonBlock(; contribs = [:M_pi_c])
        @test_throws ArgumentError reduction_declarations([EnergyPools(), m])
    end

    @testset "1.5 standalone resolution of a contributor still succeeds" begin
        graph = resolve_coupling(CarbonBlock())
        @test :M_atp_c in graph.unowned_states
        @test :M_pi_c in graph.unowned_states
    end

    @testset "AD: contributions differentiate, including into a Float64 owner derivative" begin
        models = [InertPool(), Source(; k = 0.7)]
        prob = build_problem(models; tspan = (0.0, 1.0))
        i_gtp = _gidx(models, :M_gtp_c)
        # InertPool's own derivative is Float64 whatever p is; Source's contribution is Dual.
        J = ForwardDiff.jacobian(p -> prob.f.f(prob.u0, p, 0.0), prob.p)
        @test J[i_gtp, 1] == 1.0
        @test all(J[i, 1] == 0.0 for i in 1:length(prob.u0) if i != i_gtp)

        models2 = [EnergyPools(), CarbonBlock()]
        prob2 = build_problem(models2; tspan = (0.0, 1.0))
        J2 = ForwardDiff.jacobian(p -> prob2.f.f(prob2.u0, p, 0.0), prob2.p)
        i_atp, i_g6p = _gidx(models2, :M_atp_c), _gidx(models2, :M_g6p_c)
        @test J2[i_atp, 2] == prob2.u0[i_g6p]          # d(du_atp)/d(k_prod)
        @test J2[i_atp, 1] == -prob2.u0[i_atp]         # d(du_atp)/d(gamma_o)
    end

    @testset "1.6 the framework adds no allocation per right-hand-side call at 32 states" begin
        prob = build_problem([EnergyPools(), CarbonBlock()]; tspan = (0.0, 1.0))
        rhs = prob.f.f
        # Static all the way through: a Vector anywhere in the path would allocate.
        @test rhs(prob.u0, prob.p, 0.0) isa SVector{32, Float64}
        _rhs_alloc(rhs, prob.u0, prob.p, 0.0)   # warm-up: compiles for these argument types
        @test _rhs_alloc(rhs, prob.u0, prob.p, 0.0) == 0
    end

    @testset "at most 31 sub-models compose" begin
        # Base.Any32 matches any tuple of 32 or more, so 32 is already the
        # non-unrolled case.
        many = [CoreAStub(Symbol(:S, i)) for i in 1:32]
        err = caught(() -> build_problem(many))
        @test err isa ErrorException
        @test occursin("31", sprint(showerror, err))
        @test build_problem(many[1:31]) isa ODEProblem
    end
end
