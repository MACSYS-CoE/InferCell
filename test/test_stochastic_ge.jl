using Statistics

@testset "StochasticGeneExpression" begin
    @testset "Construction" begin
        model = StochasticGeneExpression()
        @test formalism(model) == :jump
        @test inference_mode(model) == :simulation
        @test states(model) == [:mRNA, :protein]
        @test length(parameters(model)) == 6
        @test length(free_params(parameters(model))) == 4
        @test length(model_free_params(parameters(model))) == 4
        @test length(obs_free_params(parameters(model))) == 0
    end

    @testset "Reactions" begin
        model = StochasticGeneExpression()
        rxns = reactions(model)
        @test length(rxns) == 4
    end

    @testset "Build JumpProblem" begin
        model = StochasticGeneExpression()
        prob = build_problem(model; tspan=(0.0, 50.0))
        @test prob isa JumpProcesses.JumpProblem
    end

    @testset "Forward simulation produces integers" begin
        model = StochasticGeneExpression()
        prob = build_problem(model; tspan=(0.0, 50.0))
        sol = solve(prob, SSAStepper(); saveat=0:5.0:50)
        @test all(x -> x == round(x), sol[1, :])
        @test all(x -> x == round(x), sol[2, :])
    end

    @testset "Steady-state mean matches theory" begin
        model = StochasticGeneExpression()
        prob = build_problem(model; tspan=(0.0, 200.0))

        n_replicates = 500
        final_mRNA = zeros(n_replicates)
        final_protein = zeros(n_replicates)

        for i in 1:n_replicates
            sol = solve(prob, SSAStepper(); saveat=[200.0])
            final_mRNA[i] = sol[1, end]
            final_protein[i] = sol[2, end]
        end

        # Theoretical: mRNA_ss = k_tx/gamma_mRNA = 2.0
        #              protein_ss = k_tl * mRNA_ss / gamma_protein = 40.0
        @test abs(mean(final_mRNA) - 2.0) < 0.5
        @test abs(mean(final_protein) - 40.0) < 5.0
    end

    @testset "Formalism validation" begin
        # This used to assert that any mixed :ode/:jump composition is refused.
        # Spec §11 phase 3 replaced that refusal with the handshake driver, so
        # the test now asserts what actually fails: these two models both own
        # :mRNA and :protein, and duplicate ownership is still an error. (Phase
        # 2 rewrote the bursty model's composition test the same way, and for
        # the same reason.)
        # `caught` lives in corea_test_models.jl, which runtests.jl includes
        # after this file, so the idiom is spelled out here.
        model_ode = TranscriptionTranslation()
        model_ssa = StochasticGeneExpression()
        err = try
            build_problem([model_ode, model_ssa])
            nothing
        catch e
            e
        end
        @test err !== nothing
        @test occursin("owned by multiple sub-models", sprint(showerror, err))
    end
end
