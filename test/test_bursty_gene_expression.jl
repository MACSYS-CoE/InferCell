using Statistics: mean, var, std
using Random

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
        Random.seed!(20260905)
        for i in 1:n_reps
            sol = solve(prob, SSAStepper(); saveat=[500.0])
            final_mRNA[i] = sol[2, end]
            final_protein[i] = sol[3, end]
        end
        # Four standard errors of the sample mean, computed from the sample
        # rather than picked. The bursty protein's variance is far above
        # Poisson — its Fano factor here is of order fifty — so a tolerance
        # chosen by eye is either vacuous or flaky, and the hand-picked 15.0
        # this replaced was the second: it failed CI at a deviation of 16.4
        # against a standard error of about 5.5. The claim being tested is
        # that the sample mean agrees with telegraph theory, and this is that
        # claim written down.
        @test abs(mean(final_mRNA) - mRNA_ss) < 4 * std(final_mRNA) / sqrt(n_reps)
        @test abs(mean(final_protein) - protein_ss) < 4 * std(final_protein) / sqrt(n_reps)
        @info "bursty steady state" mRNA_mean = mean(final_mRNA) mRNA_theory = mRNA_ss protein_mean = mean(final_protein) protein_theory = protein_ss protein_4se = 4 * std(final_protein) / sqrt(n_reps)
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

    @testset "Composition with another jump model fails on duplicate ownership" begin
        # Both models own :mRNA and :protein, so the composition is refused for
        # duplicate ownership — the only thing that ever made it fail. Before
        # spec phase 2 this test read as if index aliasing were caught; it was
        # not, and two jump models with disjoint states now compose correctly
        # (test_jump_composition.jl). Mixed formalisms compose as of phase 3
        # and are exercised in test_hybrid_handshake.jl.
        m1 = BurstyGeneExpression()
        m2 = StochasticGeneExpression()
        @test_throws "State :mRNA is owned by multiple sub-models" build_problem([m1, m2])
    end
end
