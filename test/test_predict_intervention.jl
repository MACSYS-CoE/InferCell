# KIM-32: interventions do(a), with paired seeds. A zero-dose knockdown gives
# ΔQ ≡ 0 exactly, on both forward descriptions.

using Test
using InferCell
using InferCell: AbstractSubModel

@testset "KIM-32: interventions with paired seeds" begin
    ms = collect(AbstractSubModel, m0_models())
    nominal = nominal_parameter_values(ms)
    θ = haldane_values(ms, [n => nominal[n] for n in M0_TARGETS])
    ptsg = promoter_param(:JCVISYN3A_0779)

    @testset "a knockdown scales its parameter and nothing else" begin
        a = Knockdown(; dose = 0.25)
        @test a.param === ptsg
        v = Dict(intervene(a, ms, θ))
        @test v[ptsg] == 0.75 * Dict(θ)[ptsg]
        @test all(v[n] == x for (n, x) in θ if n !== ptsg)
        # A parameter the writes leave out starts from its nominal value.
        w = intervene(Knockdown(; param = promoter_param(:JCVISYN3A_0607), dose = 0.5), ms,
                      [ptsg => 1.0])
        @test Dict(w)[promoter_param(:JCVISYN3A_0607)] ==
              0.5 * nominal[promoter_param(:JCVISYN3A_0607)]
        @test Dict(intervene(Knockdown(; dose = 0.0), ms, θ)) == Dict(θ)
        # A knockdown of a forward constant moves its reverse constant with it, so
        # the equilibrium constant stays nominal, as it does at baseline. Every
        # intervened arm (predict's, intervention_cells', δ's) goes through here.
        e = Dict(intervene(Knockdown(; param = :kcatF_R_ENO, dose = 0.5), ms, θ))
        @test e[:kcatF_R_ENO] == 0.5 * Dict(θ)[:kcatF_R_ENO]
        @test e[:kcatR_R_ENO] ≈ 0.5 * Dict(θ)[:kcatR_R_ENO] rtol = 1e-12
        @test e[:kcatR_R_FBA] == Dict(θ)[:kcatR_R_FBA]
        @test all(r -> r.relative <= r.bound, assert_haldane(ms, collect(e)))
        # A zero dose writes exactly the baseline's writes, even for a Haldane-linked
        # parameter they leave out; and a reverse constant, which is derived, is
        # refused as a target.
        @test intervene(Knockdown(; param = :kcatF_R_PGK3, dose = 0.0), ms, θ) == θ
        @test_throws ArgumentError intervene(Knockdown(; param = :kcatR_R_ENO, dose = 0.5), ms, θ)
        @test_throws ArgumentError Knockdown(; dose = 1.5)
        @test_throws ArgumentError intervene(Knockdown(; param = :no_such, dose = 0.1), ms, θ)
    end

    @testset "zero dose: ΔQ ≡ 0 exactly, on both forward descriptions" begin
        d0 = Knockdown(; dose = 0.0)
        hy = intervention_cells(ms, θ, d0, [11, 12]; horizon = 60.0)
        @test all(==(0.0), hy.delta) && length(hy.delta) == 2
        @test [c.q for c in hy.baseline] == [c.q for c in hy.intervened]
        mf = intervention_cells(ms, θ, d0, Int[]; forward = :mean_field, horizon = 60.0)
        @test only(mf.delta) === 0.0
    end

    @testset "a real knockdown moves Q, on the same cells" begin
        r = intervention_cells(ms, θ, Knockdown(; dose = 0.5), [11, 12]; horizon = 60.0)
        @test all(isfinite, r.delta) && any(!=(0.0), r.delta)
        @test [c.seed for c in r.baseline] == [c.seed for c in r.intervened] == [11, 12]
        mf = intervention_cells(ms, θ, Knockdown(; dose = 0.5), Int[];
                                forward = :mean_field, horizon = 60.0)
        @test only(mf.delta) != 0.0
    end

    @testset "δ covers the intervened arm on the same seeds" begin
        dc = delta_correction(ms, θ, [21, 22]; horizon = 60.0,
                              intervention = Knockdown(; dose = 0.0))
        @test dc.delta_intervened == dc.delta
        @test [c.seed for c in dc.hybrid_intervened] == dc.seeds
        @test delta_correction(ms, θ, [21, 22]; horizon = 60.0).intervention === nothing
    end
end
