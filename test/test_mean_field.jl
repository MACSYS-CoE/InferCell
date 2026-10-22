# spec/phases/approximate-method.md task B2: mean-field gene expression, AV3 and
# AV4.

using Random
using Statistics: mean, std
using Test
using InferCell
using ForwardDiff
using Distributions: LogNormal
using JumpProcesses: SSAStepper
using OrdinaryDiffEq: Tsit5, solve, remake

# A jump module whose propensity is quadratic, so its mean equation does not
# close: MeanField must refuse it.
struct _SquareDecayMF <: AbstractSubModel end
InferCell.states(::_SquareDecayMF) = [:x_sq]
InferCell.parameters(::_SquareDecayMF) =
    [InferParameter(0.1, LogNormal(log(0.1), 1.0), false, :k_sq, :_SquareDecayMF, :rate)]
InferCell.formalism(::_SquareDecayMF) = :jump
InferCell.reactions(::_SquareDecayMF) =
    [Reaction((u, p, t, w) -> p[1] * u[1]^2, (u, w) -> (u[1] -= 1))]

# Affine at its declared values only: the quadratic term's coefficient is zero
# there, so MeanField must also check away from them.
struct _HiddenSquareMF <: AbstractSubModel end
InferCell.states(::_HiddenSquareMF) = [:x_hs]
InferCell.parameters(::_HiddenSquareMF) =
    [InferParameter(0.1, LogNormal(log(0.1), 1.0), false, :k_hs, :_HiddenSquareMF, :rate),
     InferParameter(0.0, LogNormal(log(0.1), 1.0), false, :q_hs, :_HiddenSquareMF, :rate)]
InferCell.formalism(::_HiddenSquareMF) = :jump
InferCell.reactions(::_HiddenSquareMF) =
    [Reaction((u, p, t, w) -> p[1] * u[1] + p[2] * u[1]^2, (u, w) -> (u[1] -= 1))]

_mf_value(m, name) = only(p.value for p in parameters(m) if p.name === name)

@testset "B2: mean-field gene expression" begin
    genes = [g for g in InferCell.read_transcription_genes() if g.locus in InferCell.M0_GENES]
    jump = AbstractSubModel[CoreATranscription(; genes), CoreATranscriptDecay(; genes),
                            CoreATranslation(; genes)]
    mf = gene_expression_mean_field(; genes)
    bare = gene_expression_mean_field(; genes, edges = :none)

    @testset "the same states, parameters and edges as the jump modules" begin
        @test states(mf) == reduce(vcat, states.(jump))
        @test [p.name for p in parameters(mf)] == reduce(vcat, [[p.name for p in parameters(m)] for m in jump])
        @test coupling(mf) == reduce(vcat, coupling.(jump))
        @test rebuilt_params(mf) == reduce(vcat, rebuilt_params.(jump))
        @test formalism(mf) === :ode && inference_mode(mf) === :differentiable
        @test isempty(coupling(bare)) && isempty(rebuilt_params(bare))
        cats = [l.category for l in reduction_declarations(AbstractSubModel[mf])]
        @test :formalism in cats && :coarse_drain in cats
    end

    prob = build_problem(AbstractSubModel[bare]; tspan = (0.0, 600.0))
    saves = 60.0:60.0:600.0

    @testset "the means match the closed form at constant constants" begin
        sol = solve(prob, Tsit5(); saveat = saves, abstol = 1e-12, reltol = 1e-10)
        g = only(i for (i, x) in enumerate(genes) if x.locus === :JCVISYN3A_0607)
        loc = genes[g].locus
        k = _mf_value(mf, InferCell.rate_param(loc))
        μ = _mf_value(mf, :krnadeg) / genes[g].length
        ktl = _mf_value(mf, InferCell.translation_rate_param(loc))
        m0 = round(Int, genes[g].mean_mrna)
        im = findfirst(==(InferCell.transcript_state(loc)), states(mf))
        ip = findfirst(==(InferCell.protein_state(loc)), states(mf))
        p0 = prob.u0[ip]
        for (t, u) in zip(sol.t, sol.u)
            @test u[im] ≈ k / μ + (m0 - k / μ) * exp(-μ * t) rtol = 1e-8
            @test u[ip] - p0 ≈ ktl * (k / μ * t + (m0 - k / μ) * (1 - exp(-μ * t)) / μ) rtol = 1e-8
        end
    end

    @testset "AV4: the SSA ensemble's means at pinned constants" begin
        # A jump-only build executes no rebuild and no drain: the constants stay
        # at their declared values and the counters accumulate, which is what
        # the mean field integrates.
        Random.seed!(8903)
        jprob = build_problem(jump; tspan = (0.0, 600.0))
        n = 2000
        runs = [solve(jprob, SSAStepper(); saveat = saves) for _ in 1:n]
        msol = solve(prob, Tsit5(); saveat = saves, abstol = 1e-12, reltol = 1e-10)
        jnames = reduce(vcat, states.(jump))
        @test jnames == states(mf)
        worst = 0.0
        for t in saves, s in eachindex(jnames)
            x = [r.u[findfirst(==(t), r.t)][s] for r in runs]
            se = std(x) / sqrt(n)
            target = msol.u[findfirst(==(t), msol.t)][s]
            if se == 0
                @test mean(x) ≈ target atol = 1e-6
            else
                z = abs(mean(x) - target) / se
                worst = max(worst, z)
                @test z < 4
            end
        end
        @info "AV4: largest |SSA mean − mean field| over SE" worst
    end

    @testset "AV3: ForwardDiff through a mean-field solve matches central differences" begin
        # ptsG's protein count at 600 s, membrane and cytosolic, against every
        # free parameter on the log scale (u = 0 is the nominal point).
        p0 = copy(prob.p)
        ia = findfirst(==(:P_JCVISYN3A_0779), states(mf))
        ib = findfirst(==(:Pcyto_JCVISYN3A_0779), states(mf))
        function obs(u)
            p = p0 .* exp.(u)
            sol = solve(remake(prob; u0 = convert.(eltype(p), prob.u0), p = p), Tsit5();
                        saveat = [600.0], abstol = 1e-12, reltol = 1e-12)
            return log(sol.u[end][ia] + sol.u[end][ib])
        end
        u0 = zeros(length(p0))
        grad = ForwardDiff.gradient(obs, u0)
        h = 1e-5
        fd = map(eachindex(u0)) do i
            e = zeros(length(u0))
            e[i] = h
            (obs(u0 + e) - obs(u0 - e)) / 2h
        end
        @test all(isfinite, grad)
        @test any(!iszero, grad)
        @test grad ≈ fd rtol = 1e-6
    end

    @testset "refusals" begin
        @test_throws ArgumentError MeanField(AbstractSubModel[_SquareDecayMF()])
        @test_throws ArgumentError MeanField(AbstractSubModel[_HiddenSquareMF()])
        # Decay alone reads transcripts no member owns.
        @test_throws ArgumentError MeanField(AbstractSubModel[CoreATranscriptDecay(; genes)])
        @test_throws ArgumentError MeanField(AbstractSubModel[bare])
        @test_throws ArgumentError gene_expression_mean_field(; genes, edges = :sometimes)
    end
end
