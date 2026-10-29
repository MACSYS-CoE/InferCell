# KIM-29: completing held parameters from their priors. The draws follow the
# prior marginals at each width, every draw that runs keeps Haldane, and the
# driver-written slots (the dead `enz_R_*`) are never drawn.

using Test
using InferCell
using Random
using Statistics: mean, std, median

@testset "KIM-29: held parameters from their priors" begin
    @testset "the default held set on the reduced model, given M0's posterior" begin
        dms = d11_models()
        h = held_parameters(dms, M0_TARGETS)
        m0_loci = Set(promoter_param.([:JCVISYN3A_0607, :JCVISYN3A_0779]))
        promoters = [n for n in h.held if startswith(String(n), "S_")]
        # d11's free parameters the driver does not write are the 17 promoters,
        # krnadeg and ENO/FBA; with M0's five inferred, 15 promoters remain.
        @test h.held == promoters && length(h.held) == 15
        @test isdisjoint(h.held, m0_loci) && promoter_param(:JCVISYN3A_0445) in h.held
        reason = Dict(s.name => s.reason for s in h.skipped)
        @test all(reason[n] === :inferred for n in M0_TARGETS)
        @test reason[:kcatR_R_ENO] === :haldane_derived && reason[:kcatR_R_FBA] === :haldane_derived
    end

    @testset "prior marginals at widths ½, 1 and 2" begin
        dms = d11_models()
        held = held_parameters(dms, M0_TARGETS).held
        nominal = nominal_parameter_values(dms)
        n = 4000
        for w in (0.5, 1.0, 2.0)
            draws = draw_held(Xoshiro(29), dms, held, n; width = w)
            @test length(draws) == n && all(d -> first.(d) == held, draws)
            for (j, name) in enumerate(held)
                x = log.([d[j].second for d in draws])
                # Each promoter's prior is LogNormal(log s, log 2), widened about
                # its median. The tolerances are about four standard errors.
                @test median(x) ≈ log(nominal[name]) atol = 4 * 1.25 * w * log(2) / sqrt(n)
                @test std(x) ≈ w * log(2) rtol = 4 / sqrt(2n)
            end
        end
        @test_throws ArgumentError draw_held(Xoshiro(1), dms, held, 0)
        @test_throws ArgumentError draw_held(Xoshiro(1), dms, held, 2; width = 0)
    end

    @testset "Haldane holds on every completed draw, and enz_R_* is never drawn" begin
        # 14c's broad composition, with every informed glycolysis and recycling
        # constant free, and ENO's catalytic slot freed too, which the driver
        # overwrites.
        informed(m) = Symbol[q.name for q in parameters(m)
                             if q.role === :rate && informedness(q) === :balanced]
        g0 = CentralGlycolysis(enzymes = :translated)
        r0 = NucleotideRecycling(enzymes = :translated)
        freed(m) = unique(vcat(informed(m), [r.reverse for r in haldane_relations(m)]))
        ms = corea_models()
        ms[1] = CentralGlycolysis(enzymes = :translated, free = vcat(freed(g0), [:enz_R_ENO]))
        ms[3] = NucleotideRecycling(enzymes = :translated, free = freed(r0))
        inferred = [:kcatF_R_ENO, :kcatF_R_FBA]
        h = held_parameters(ms, inferred)
        reason = Dict(s.name => s.reason for s in h.skipped)
        @test reason[:enz_R_ENO] === :driver_written
        @test !any(startswith(String(n), "enz_R_") for n in h.held)
        @test all(r -> !(r.reverse in h.held), haldane_relations(ms))
        # The ENO Michaelis constants are held, so ENO's reverse constant depends
        # on a posterior value and a prior value together.
        @test :km_R_ENO_M_2pg_c in h.held && :km_R_ENO_M_pep_c in h.held

        nominal = nominal_parameter_values(ms)
        rng = Xoshiro(2929)
        for draw in draw_held(rng, ms, h.held, 50; width = 2.0)
            θI = [:kcatF_R_ENO => nominal[:kcatF_R_ENO] * exp(randn(rng)),
                  :kcatF_R_FBA => nominal[:kcatF_R_FBA] * exp(randn(rng))]
            vals = haldane_values(ms, vcat(θI, draw))
            @test length(assert_haldane(ms, vals)) == 15
            v = Dict(vals)
            @test v[:kcatF_R_ENO] == θI[1].second    # the posterior's value runs as drawn
            @test !haskey(v, :enz_R_ENO)
        end
        # A reverse constant passed in is replaced by the derived one.
        v = Dict(haldane_values(ms, [:kcatF_R_ENO => 2.0, :kcatR_R_ENO => 1e9]))
        @test v[:kcatR_R_ENO] != 1e9
        @test_throws ArgumentError draw_held(rng, ms, [:kcatR_R_ENO], 1)
        @test_throws ArgumentError haldane_values(ms, [:no_such => 1.0])
    end

    @testset "width scaling refuses a prior it cannot scale" begin
        ms = d11_models()
        @test_throws ArgumentError draw_parameters(Xoshiro(1), ms, [:krnadeg]; width = -1)
        # Width 1 is the plain draw, the same stream giving the same values.
        @test draw_parameters(Xoshiro(5), ms, D11_TARGETS) ==
              draw_parameters(Xoshiro(5), ms, D11_TARGETS; width = 1.0)
    end
end
