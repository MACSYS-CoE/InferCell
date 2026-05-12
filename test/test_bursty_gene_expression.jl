using Statistics: mean, var

@testset "BurstyGeneExpression" begin
    @testset "Construction" begin
        model = BurstyGeneExpression()
        @test formalism(model) == :jump
        @test inference_mode(model) == :simulation
        @test states(model) == [:promoter, :mRNA, :protein]
        @test length(parameters(model)) == 9
        @test length(free_params(parameters(model))) == 6
        @test length(model_free_params(parameters(model))) == 6
        @test length(obs_free_params(parameters(model))) == 0
    end

    @testset "Reactions" begin
        model = BurstyGeneExpression()
        rxns = reactions(model)
        @test length(rxns) == 6
    end

    @testset "Build JumpProblem" begin
        model = BurstyGeneExpression()
        prob = build_problem(model; tspan=(0.0, 50.0))
        @test prob isa JumpProcesses.JumpProblem
    end

    @testset "Forward simulation produces integers" begin
        model = BurstyGeneExpression()
        prob = build_problem(model; tspan=(0.0, 50.0))
        sol = solve(prob, SSAStepper(); saveat=0:5.0:50)
        @test all(x -> x == round(x), sol[1, :])  # promoter ∈ {0,1}
        @test all(x -> x == 0 || x == 1, sol[1, :])
        @test all(x -> x == round(x), sol[2, :])  # mRNA
        @test all(x -> x == round(x), sol[3, :])  # protein
    end

    @testset "Promoter ON fraction matches theory" begin
        # Stationary P(ON) = k_on / (k_on + k_off). Average across independent
        # replicates at a single late time, mirroring the steady-state mean test.
        model = BurstyGeneExpression(k_on=0.2, k_off=0.5)
        prob = build_problem(model; tspan=(0.0, 200.0))
        n_reps = 500
        on_samples = zeros(n_reps)
        for i in 1:n_reps
            sol = solve(prob, SSAStepper(); saveat=[200.0])
            on_samples[i] = sol[1, end]
        end
        expected = 0.2 / (0.2 + 0.5)
        @test abs(mean(on_samples) - expected) < 0.05
    end

    @testset "Steady-state mean matches telegraph theory" begin
        # Mean mRNA = k_tx_burst * (k_on / (k_on + k_off)) / gamma_mRNA
        # Mean protein = k_tl * mean_mRNA / gamma_protein
        k_on, k_off = 0.2, 0.5
        k_tx_burst, k_tl = 10.0, 2.0
        g_m, g_p = 0.5, 0.1
        on_frac = k_on / (k_on + k_off)
        mRNA_ss = k_tx_burst * on_frac / g_m  # ≈ 5.71
        protein_ss = k_tl * mRNA_ss / g_p     # ≈ 114.3

        model = BurstyGeneExpression(; k_on, k_off, k_tx_burst, k_tl,
                                     gamma_mRNA=g_m, gamma_protein=g_p)
        prob = build_problem(model; tspan=(0.0, 500.0))

        n_reps = 200
        final_mRNA = zeros(n_reps)
        final_protein = zeros(n_reps)
        for i in 1:n_reps
            sol = solve(prob, SSAStepper(); saveat=[500.0])
            final_mRNA[i] = sol[2, end]
            final_protein[i] = sol[3, end]
        end
        @test abs(mean(final_mRNA) - mRNA_ss) < 1.5
        @test abs(mean(final_protein) - protein_ss) < 15.0
    end

    @testset "Bursting: mRNA Fano factor > 1" begin
        # Constitutive Poisson gives Fano = 1; bursting inflates it.
        # Choose k_off >> k_on so bursts are well-separated.
        model = BurstyGeneExpression(k_on=0.05, k_off=1.0, k_tx_burst=20.0,
                                     gamma_mRNA=0.5)
        prob = build_problem(model; tspan=(0.0, 2000.0))
        sol = solve(prob, SSAStepper(); saveat=200.0:5.0:2000.0)
        mRNA_trace = sol[2, :]
        mu, var_m = mean(mRNA_trace), var(mRNA_trace)
        fano = var_m / max(mu, 1e-6)
        @test fano > 1.5
    end

    @testset "Composition with another jump model fails (mixed reactions over same state)" begin
        # Mixed formalisms are disallowed; that path is exercised in test_stochastic_ge.jl.
        # Composing two jump models that own the same state name should error in coupling.
        m1 = BurstyGeneExpression()
        m2 = StochasticGeneExpression()
        @test_throws ErrorException build_problem([m1, m2])
    end
end
