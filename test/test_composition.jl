using InferCell
using Test
using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: ReturnCode

@testset "TX/TL + Metabolism Composition" begin
    txl = TranscriptionTranslation()
    metab = LightMetabolism()

    @testset "build_problem structure" begin
        prob = build_problem([txl, metab]; tspan=(0.0, 100.0))

        @test length(prob.u0) == 5   # 2 TX/TL + 3 metabolism
        @test length(prob.p) == 7    # 4 TX/TL + 3 unique metabolism (k_tx, k_tl deduplicated)
        @test prob.tspan == (0.0, 100.0)
    end

    @testset "Shared parameter deduplication" begin
        prob = build_problem([txl, metab])

        # Flat vector: k_tx, k_tl, gamma_mRNA, gamma_protein, k_atp, k_ntp, k_aa
        @test prob.p[1] == 1.0   # k_tx
        @test prob.p[2] == 2.0   # k_tl
        @test prob.p[3] == 0.5   # gamma_mRNA
        @test prob.p[4] == 0.1   # gamma_protein
        @test prob.p[5] == 1.0   # k_atp
        @test prob.p[6] == 0.2   # k_ntp
        @test prob.p[7] == 0.4   # k_aa
    end

    @testset "Forward simulation" begin
        prob = build_problem([txl, metab]; tspan=(0.0, 200.0))
        sol = solve(prob, Tsit5(); saveat=1.0)

        @test sol.retcode == ReturnCode.Success
        @test !any(isnan, sol.u[end])
        @test !any(isinf, sol.u[end])

        # All states should be finite and non-negative at steady state
        final = sol.u[end]
        for i in 1:5
            @test final[i] >= -0.1  # allow small numerical undershoot
        end
    end

    @testset "TX/TL dynamics unchanged by coupling" begin
        prob = build_problem([txl, metab]; tspan=(0.0, 200.0))
        sol = solve(prob, Tsit5(); saveat=1.0)

        # mRNA steady state = k_tx / gamma_mRNA = 1.0 / 0.5 = 2.0
        mRNA_ss = sol.u[end][1]
        @test isapprox(mRNA_ss, 2.0; atol=0.1)

        # protein steady state = k_tl * mRNA_ss / gamma_protein = 2.0 * 2.0 / 0.1 = 40.0
        protein_ss = sol.u[end][2]
        @test isapprox(protein_ss, 40.0; atol=1.0)
    end

    @testset "Metabolites reach positive steady state" begin
        prob = build_problem([txl, metab]; tspan=(0.0, 200.0))
        sol = solve(prob, Tsit5(); saveat=1.0)
        final = sol.u[end]

        @test final[3] > 0  # ATP > 0
        @test final[4] > 0  # NTP > 0
        @test final[5] > 0  # AA > 0
    end

    @testset "Inconsistent shared params error" begin
        bad_metab = LightMetabolism(k_tx=999.0)
        @test_throws ErrorException build_problem([txl, bad_metab])
    end

    @testset "Turing model parameter deduplication" begin
        prob = build_problem([txl, metab])
        all_model_free = unique_params(reduce(vcat,
            model_free_params.(parameters.([txl, metab]))))
        all_obs_free = unique_params(reduce(vcat,
            obs_free_params.(parameters.([txl, metab]))))

        @test length(all_model_free) == 7   # 7 unique ODE params
        @test length(all_obs_free) == 1     # 1 shared sigma_obs

        names = [p.name for p in all_model_free]
        @test :k_tx in names
        @test :k_tl in names
        @test count(==(:k_tx), names) == 1  # deduplicated
        @test count(==(:k_tl), names) == 1
    end

    @testset "Backward compatibility: single model" begin
        prob = build_problem([txl]; tspan=(0.0, 100.0))
        sol = solve(prob, Tsit5(); saveat=1.0)

        @test sol.retcode == ReturnCode.Success
        @test length(prob.u0) == 2
        @test length(prob.p) == 4
    end
end
