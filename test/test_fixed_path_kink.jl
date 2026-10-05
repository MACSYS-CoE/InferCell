# spec/phases/16-recovery.md task 16a.6 (parent 16.3): conditioning on the path
# does not make block 1 smooth. Under a fixed recorded path, a clipped drain
# still puts a kink in the ODE state as a function of an ODE parameter, because
# the pool the drain is paid from depends on that parameter. Check 7's census
# finds the clip. The smoothed drain on continuous pools (14c.5's sampler model)
# takes the kink out, and its slope jump shrinks with the spacing. On Core A′,
# with ENO's constant, this is dev/scripts/fixed_path_kink_16a6.jl.

using Random
using Test
using InferCell

@testset "16a.6: a clipped drain under a fixed path still kinks" begin
    f = corea_particles_per_mM()
    # A small ATP pool against a costly translation, so the drain clips.
    pool(θ) = ToyPool(kcat = θ, atp0 = 400 / f)
    expr(clip) = ToyExpression(k_tx = 0.5, k_tl = 2.0, cost = 40, protein0 = 200,
        edges = [DeferredCounterEdge(; species = :M_atp_c, direction = :in,
                                     counter = :atp_cost, clip,
                                     (clip === :smoothed ? (; smoothing = 1.0) : (;))...)])
    H = 60.0
    variants = (published = (:clamped_deficit_carried, :fractional_carry),
                kink = (:clamped_deficit_carried, :continuous),
                sampler = (:smoothed, :continuous))
    build(θ, v) = build_problem([pool(θ), expr(variants[v][1])]; tspan = (0.0, H),
                                rounding = variants[v][2], abstol = 1e-14, reltol = 1e-12)

    # Record a path on the published model, and find its first clip by check 7's
    # census, a drain that carried a deficit.
    θ0 = 30.0
    Random.seed!(1606)
    d = build(θ0, :published)
    Random.seed!(1606)
    record_path!(d)
    nstar = 0
    for k in 1:round(Int, H)
        n = d.n_clipped
        handshake_step!(d)
        d.n_clipped > n && nstar == 0 && (nstar = k)
    end
    @test clipping_census(d).clipped > 0
    @test nstar > 0
    path = recorded_path(d)

    # The ATP state after the drain at handshake nstar, in particles, and whether
    # that drain clipped, with the path held and θ varied from t = 0.
    function after(θ, v)
        e = build(θ, v)
        r = PathReplay(path, e)
        for _ in 1:nstar
            replay_step!(e, r)
        end
        return e.ode.u[1] * e.factor, e.debits[1].clipped
    end
    @test last(after(θ0, :published))                 # the replay reproduces the clip

    # The clip's θ, bisected on the clamp over continuous pools. kcat drains ATP,
    # so a smaller kcat leaves more to pay with.
    lo, hi = θ0, θ0
    while last(after(lo, :kink))
        hi = lo
        lo /= 2
        lo > 1e-6 || error("no kcat stops the clip")
    end
    for _ in 1:60
        mid = sqrt(lo * hi)
        last(after(mid, :kink)) ? (hi = mid) : (lo = mid)
    end
    θc = sqrt(lo * hi)

    slope(θ, v, δ) = (first(after(θ * exp(δ), v)) - first(after(θ * exp(-δ), v))) / 2δ
    dP = slope(θc / 1.5, :kink, 1e-3)                # the paying side's sensitivity
    @test dP != 0
    δ = 0.01 / abs(dP)
    jump(k, v) = abs(slope(θc * exp(k * δ), v, δ) - slope(θc * exp(-k * δ), v, δ)) / abs(dP)

    # The clamp under a fixed path: the jump does not shrink with the spacing.
    for k in (30, 3)
        @test jump(k, :kink) > 0.5
    end
    @test jump(3, :kink) > 0.8 * jump(30, :kink)
    # The published model is a step function of θ: zero slope almost everywhere.
    grid = [θc * exp(k * 37δ) for k in -10:10]
    @test count(θ -> slope(θ, :published, δ) == 0, grid) >= 15
    # 14c.5's sampler model: the jump shrinks with the spacing.
    far, near = jump(30, :sampler), jump(3, :sampler)
    @test near < far / 5
    @info "16a.6 toy: slope jump across the clip under a fixed path" nstar θc dP clamp_30 = jump(30, :kink) clamp_3 = jump(3, :kink) sampler_30 = far sampler_3 = near
end
