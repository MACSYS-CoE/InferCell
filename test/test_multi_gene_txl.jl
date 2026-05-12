using OrdinaryDiffEq: Tsit5, solve
using SciMLBase: ReturnCode

@testset "Multi-gene TranscriptionTranslation" begin
    @testset "Default no-arg constructor unchanged (single-gene backward compat)" begin
        m = TranscriptionTranslation()
        @test states(m) == [:mRNA, :protein]
        @test length(parameters(m)) == 7
        @test length(model_free_params(parameters(m))) == 4
        @test m.gene_names == [:_default]
        @test m.enzyme_slot === nothing
        @test InferCell.enzyme_protein_state(m) == :protein
    end

    @testset "Multi-gene state naming" begin
        m = TranscriptionTranslation([:enzyme, :ribosome, :reporter])
        @test states(m) == [
            :enzyme_mRNA, :enzyme_protein,
            :ribosome_mRNA, :ribosome_protein,
            :reporter_mRNA, :reporter_protein,
        ]
        @test m.gene_names == [:enzyme, :ribosome, :reporter]
        @test m.enzyme_slot == :enzyme  # auto-selected because :enzyme is present
        @test InferCell.enzyme_protein_state(m) == :enzyme_protein
    end

    @testset "Multi-gene parameter naming and counts" begin
        m = TranscriptionTranslation([:enzyme, :reporter])
        # 4 rate + 2 IC per gene × 2 genes + 1 obs = 13
        @test length(parameters(m)) == 13
        @test length(model_free_params(parameters(m))) == 8  # 4 rates × 2 genes
        @test length(obs_free_params(parameters(m))) == 1

        names = [p.name for p in parameters(m)]
        for g in (:enzyme, :reporter)
            @test Symbol("k_tx_", g) in names
            @test Symbol("k_tl_", g) in names
            @test Symbol("gamma_mRNA_", g) in names
            @test Symbol("gamma_protein_", g) in names
            @test Symbol(g, "_mRNA0") in names
            @test Symbol(g, "_protein0") in names
        end
    end

    @testset "Per-gene parameter overrides" begin
        m = TranscriptionTranslation([:enzyme, :ribosome];
            overrides=Dict(:ribosome => (k_tx=0.3, k_tl=5.0)))
        pmap = Dict(p.name => p.value for p in parameters(m))
        @test pmap[:k_tx_ribosome] == 0.3
        @test pmap[:k_tl_ribosome] == 5.0
        @test pmap[:k_tx_enzyme] == 1.0  # default
    end

    @testset "enzyme_slot validation" begin
        @test_throws ErrorException TranscriptionTranslation([:gene_a]; enzyme_slot=:not_a_gene)
        @test_throws ErrorException TranscriptionTranslation(Symbol[])
    end

    @testset "Build + forward simulation (3 genes)" begin
        m = TranscriptionTranslation([:enzyme, :ribosome, :reporter])
        prob = build_problem(m; tspan=(0.0, 100.0))
        @test length(prob.u0) == 6
        sol = solve(prob, Tsit5(); saveat=0:10.0:100.0)
        @test sol.retcode == ReturnCode.Success
        @test !any(isnan, sol.u[end])
        # All per-gene steady states should equal the single-gene defaults: mRNA=2, protein=40
        for i in 1:3
            @test isapprox(sol.u[end][2i - 1], 2.0; atol=0.2)
            @test isapprox(sol.u[end][2i],     40.0; atol=2.0)
        end
    end

    @testset "Composition with metabolism (multi-gene Block 2 + Block 1)" begin
        txl = TranscriptionTranslation([:enzyme, :ribosome, :reporter])
        metab = LightMetabolism(mRNA_source=:enzyme_mRNA, enzyme_source=:enzyme_protein)
        prob = build_problem([txl, metab]; tspan=(0.0, 200.0))
        @test length(prob.u0) == 9  # 6 from txl + 3 from metab
        sol = solve(prob, Tsit5(); saveat=0:10.0:200.0)
        @test sol.retcode == ReturnCode.Success
        # Metabolites positive and bounded
        @test sol.u[end][7] > 0          # ATP
        @test sol.u[end][7] < metab.ATP_max
        @test sol.u[end][8] > 0          # NTP
        @test sol.u[end][9] > 0          # AA
    end
end
