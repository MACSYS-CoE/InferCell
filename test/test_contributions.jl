using Test
using InferCell
using StaticArrays: SVector
using OrdinaryDiffEq: Tsit5, solve
using ForwardDiff
using SciMLBase: ReturnCode

# The contribution channel: a module adds derivative terms to a registry state
# another module owns (spec §11 phase 1). The doubles are in
# contribution_test_models.jl and the CoreAStub in corea_test_models.jl.

# Everything passed as an argument and nothing captured: measured inside a
# @testset, `@allocated` would count the boxed return value of the closure.
_rhs_alloc(rhs, u, p, t) = @allocated rhs(u, p, t)

_gidx(models, s) = findfirst(==(s), reduce(vcat, states.(models)))

@testset "Contributions" begin
    @testset "1.3 owner's derivative is its own term plus the contribution, exactly" begin
        models = [EnergyPools(), CarbonBlock()]
        prob = build_problem(models; tspan = (0.0, 10.0))
        @test length(prob.u0) == 32
        rhs = prob.f.f
        u, p = prob.u0, prob.p
        du = rhs(u, p, 0.0)
        @test du isa SVector{32, Float64}

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

end
