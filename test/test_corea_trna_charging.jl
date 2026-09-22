using Test
using InferCell
using StaticArrays: SA

# Spec §11 phase 9 — the lumped tRNA charging step, in the ODE block.

@testset "Core A′ tRNA charging" begin

    @testset "9.1 one mass-action reaction owning the tRNA pair" begin
        m = TrnaCharging(k_chg = 0.5, trna0 = 0.04, trna_chg0 = 0.16)

        # The state set is the registry's tRNA group, in strictly increasing
        # registry order, so a reordered or extended group cannot slip through.
        @test states(m) == species_in_group(:trna)
        @test states(m) == [:M_trna_c, :M_trna_chg_c]
        @test all(diff(species_index.(states(m))) .> 0)
        @test formalism(m) === :ode

        # Hand-checked derivative. At [trna] = 0.04 mM, [ATP] = 3.0 mM and
        # k_chg = 0.5 /(mM·s), v = 0.5 · 0.04 · 3.0 = 0.06 mM/s exactly in
        # decimal; the tRNA pair moves by ∓v and ATP, AMP, PPi by −v, +v, +v.
        u = SA[0.04, 0.16]
        p = Float64[]                       # k_chg fixed: nothing in the vector
        v = 0.5 * 0.04 * 3.0
        @test charging_flux(u, p, 0.0, m, SA[3.0]) == v
        @test dynamics(u, p, 0.0, m, SA[3.0]) == SA[-v, v]
        @test contributions(u, p, 0.0, m, SA[3.0]) == SA[-v, v, v]
        @test v ≈ 0.06 rtol = 1e-15

        # Freed, k_chg is read from the parameter vector and not the struct.
        mf = TrnaCharging(k_chg = 0.5, trna0 = 0.04, trna_chg0 = 0.16, free_k_chg = true)
        @test charging_flux(u, [2.0], 0.0, mf, SA[3.0]) == 2.0 * 0.04 * 3.0

        # The initial conditions reach the composed u0 by name.
        @test [p.value for p in ic_params(parameters(m))] == [0.04, 0.16]
    end

end
