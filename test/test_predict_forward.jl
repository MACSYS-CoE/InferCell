# KIM-31: predict's forward descriptions. At a fixed θ, an independent hybrid
# ensemble's mean agrees with the mean field plus δ within Monte Carlo error.
# Q is read from the metered run, and a failed run is kept with its reason.

using Test
using InferCell
using InferCell: AbstractSubModel, _forward_setup, _mean_field_cell
using Random
using Statistics: mean, std

# The suite loads the toy doubles first; a run of this file alone loads them here.
isdefined(@__MODULE__, :ToyPool) ||
    foreach(include, ["corea_test_models.jl", "contribution_test_models.jl",
                      "jump_test_models.jl", "hybrid_test_models.jl"])

@testset "KIM-31: hybrid and mean-field forward, with δ" begin
    H = 120.0
    ms = collect(AbstractSubModel, m0_models())
    nominal = nominal_parameter_values(ms)
    θ = haldane_values(ms, [n => nominal[n] for n in M0_TARGETS])

    @testset "metering" begin
        mm = metered(ms)
        @test mm[2] isa MeteredPtsTransport && ms[2] isa PtsTransport
        @test metered(mm)[2] === mm[2]
        @test_throws ArgumentError metered(AbstractSubModel[ToyPool()])
    end

    @testset "δ: the hybrid mean is the mean field plus δ, within MC error" begin
        dc = delta_correction(ms, θ, 101:106; horizon = H)
        @test dc.n_ok == 6 && dc.delta == dc.hybrid_mean - dc.mean_field.q
        @test dc.quantity === :lactate_out && dc.modules == [module_id(m) for m in metered(ms)]
        fresh = forward_cells(ms, θ, 201:206; horizon = H)
        q = [c.q for c in fresh]
        @test all(c -> c.failure === :none, fresh)
        mf = only(forward_cells(ms, θ, Int[]; forward = :mean_field, horizon = H))
        @test mf.q == dc.mean_field.q                 # the mean field is deterministic
        se = sqrt(dc.hybrid_se^2 + std(q)^2 / length(q))
        @test abs(mean(q) - (mf.q + dc.delta)) <= 4se
        @info "KIM-31: hybrid mean against the mean field plus δ" mean(q) mf.q dc.delta se
        @test_throws ArgumentError delta_correction(ms, θ, [1, 1]; horizon = H)
    end

    @testset "a seed fixes a cell; Q is read from the meters" begin
        a = only(forward_cells(ms, θ, [7]; horizon = H))
        b = only(forward_cells(ms, θ, [7]; horizon = H))
        @test a.q == b.q && a.minima == b.minima
        @test a.q == a.accounts.lactate_out > 0
        @test a.drains > 0 && a.seed == 7
        @test !haskey(a.minima, LAC_EXPORT_METER) && haskey(a.minima, :M_gtp_c)
        e = only(forward_cells(ms, θ, [7]; horizon = H, quantity = :exported))
        @test e.q == a.accounts.exported && 0 < e.q <= 1
        t = only(forward_cells(ms, θ, [7]; horizon = H, quantity = run -> run.t[end]))
        @test t.q == H
        mfx = only(forward_cells(ms, θ, Int[]; forward = :mean_field, horizon = H,
                                 quantity = :exported))
        @test mfx.q == mfx.accounts.exported && isempty(mfx.clips)
        @test_throws ArgumentError forward_cells(ms, θ, [7]; horizon = H, quantity = :growth)
        @test_throws ArgumentError forward_cells(ms, θ, [7]; horizon = H, forward = :ode)
    end

    @testset "the global generator is left as it was" begin
        Random.seed!(31)
        forward_cells(ms, θ, [3]; horizon = 60.0)
        after = rand()
        Random.seed!(31)
        @test after == rand()
    end

    @testset "an overdrawn mean field is a failed run, kept with its reason" begin
        # route (ii)'s overdraw toy: a large cost on a small, unreplenished pool.
        expr = ToyRebuiltExpression(; k_tx = 2.0, mrna0 = 50, cost = 10_000)
        f = _forward_setup(AbstractSubModel[ToyPool(; kcat = 0.0), expr], :mean_field;
                           horizon = 600.0, complete = false)
        c = _mean_field_cell(f, Pair{Symbol, Float64}[], run -> run.t[end])
        @test c.failure === :overdraw && isnan(c.q) && occursin("went below zero", c.message)
        @test c.accounts === nothing
    end
end
