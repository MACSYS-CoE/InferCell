using Test
using InferCell
using StaticArrays: SA

# Spec §11 phase 9 — the lumped tRNA charging step, in the ODE block.

@testset "Core A′ tRNA charging" begin

    @testset "9.1 one mass-action reaction owning the tRNA pair" begin
        m = TrnaCharging(k_chg = 0.5, pool_mM = 0.2, charged_fraction = 0.8)

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
        mf = TrnaCharging(k_chg = 0.5, pool_mM = 0.2, charged_fraction = 0.8,
                          free_k_chg = true)
        @test charging_flux(u, [2.0], 0.0, mf, SA[3.0]) == 2.0 * 0.04 * 3.0

        # The initial conditions reach the composed u0 by name.
        @test [p.value for p in ic_params(parameters(m))] ≈ [0.04, 0.16] rtol = 1e-15
    end

    @testset "9.2 the pool is asserted and k_chg is derived from it" begin
        m = TrnaCharging()
        byname = Dict(p.name => p for p in parameters(m))

        # All three asserted quantities load as `:asserted`, and all three reach
        # the composed model's asserted-prior enumeration.
        for n in (:k_chg, :trna_pool_mM, :trna_charged_fraction)
            @test informedness(byname[n]) === :asserted
        end
        @test Set(p.name for p in asserted_prior_params(parameters(m))) ==
              Set([:k_chg, :trna_pool_mM, :trna_charged_fraction])

        # k_chg is computed, not typed in: recompute it here from the published
        # count, independently of the module's own constants, and check the
        # defaults against the values recorded in spec §12 (2026-09-23).
        demand_mM = (3_484_518 / 6300) / corea_particles_per_mM()
        @test demand_mM ≈ 0.027408 rtol = 1e-4
        @test byname[:k_chg].value ≈ demand_mM / (0.2 * 0.25 * 3.6529) rtol = 1e-14
        @test byname[:k_chg].value ≈ 0.15006 rtol = 1e-4
        @test charging_derivation(m).k_chg == byname[:k_chg].value

        # The initial pools: 1,009 uncharged and 4,036 charged particles.
        u0 = [p.value for p in ic_params(parameters(m))]
        @test u0 ≈ [0.05, 0.20] rtol = 1e-14
        @test u0 .* corea_particles_per_mM() ≈ [1009.02, 4036.08] rtol = 1e-5

        # Change the pool and k_chg follows: doubling the total halves it, and
        # moving the split from 0.8 to 0.9 charged doubles it.
        @test TrnaCharging(pool_mM = 0.5).k_held ≈ m.k_held / 2 rtol = 1e-14
        @test TrnaCharging(charged_fraction = 0.9).k_held ≈ 2 * m.k_held rtol = 1e-12

        # At the nominal pools the derived constant delivers the demand exactly.
        atp = species_entry(:M_atp_c).initial_value
        @test charging_flux(SA[u0...], Float64[], 0.0, m, SA[atp]) ≈ demand_mM rtol = 1e-14

        # An override is recorded in provenance, not hidden.
        o = TrnaCharging(k_chg = 0.5)
        @test occursin("override",
                       only(p for p in parameters(o) if p.name === :k_chg).provenance.identifier)

        # Only k_chg can be freed; the pool quantities are always fixed.
        @test [p.name for p in free_params(parameters(TrnaCharging(free_k_chg = true)))] == [:k_chg]
        @test isempty(free_params(parameters(m)))
        @test_throws ArgumentError TrnaCharging(charged_fraction = 1.0)
        @test_throws ArgumentError TrnaCharging(pool_mM = 0.0)
    end

end
