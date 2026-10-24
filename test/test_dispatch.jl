# spec/phases/approximate-method.md task B6: the dispatcher over the three routes,
# and its verification clause (§12 2026-10-24, "B6 as planned").

using Test
using InferCell
using InferCell: _mixed_solve, AbstractSubModel
using Random
using Turing: sample, NUTS
using MCMCChains: Chains

@testset "B6: the dispatcher" begin
    @testset "selection" begin
        obs = ObservedData([0.0, 1.0], zeros(1, 2), [:M_atp_c])
        cells = (times = [60.0], genes = [:g], seeds = [1], replicate_seed = 1,
                 transcripts = zeros(Int, 1, 1, 1), bulk = zeros(1, 1))
        hybrid = collect(AbstractSubModel, m0_models())
        toy = AbstractSubModel[ToyPool(), ToyExpression()]
        ode = AbstractSubModel[CoreAStub(:DiffA; st = [:M_atp_c], mode = :differentiable)]

        @test inference_route(ode, obs) === :nuts
        @test inference_route(ode, obs; route = :nuts) === :nuts
        @test inference_route(hybrid, cells; route = :exact) === :exact
        @test inference_route(hybrid, cells; route = :mixed) === :mixed

        msg(f) = (e = caught(f); @test e isa ArgumentError; sprint(showerror, e))
        @test occursin("does not apply here", msg(() -> inference_route(ode, obs; route = :mixed)))
        m = msg(() -> inference_route(toy, obs))
        @test occursin("no path for a hybrid composition", m) && occursin("route = :exact", m)
        m = msg(() -> inference_route(hybrid, cells))
        @test occursin("one must be named", m) && occursin(":mixed", m)
        @test occursin("unknown route :gibbs", msg(() -> inference_route(hybrid, cells; route = :gibbs)))
        @test occursin("only for a hybrid composition", msg(() -> inference_route(ode, cells; route = :exact)))
        m = msg(() -> inference_route(toy, cells; route = :exact))
        @test occursin("built for M0's gene expression", m) && occursin("has 0 and 0", m)
        short = Base.structdiff(cells, NamedTuple{(:bulk,)})
        @test occursin("lacks bulk", msg(() -> inference_route(hybrid, short; route = :mixed)))
        # infer refuses before any work, by the same rules.
        @test occursin("one must be named", msg(() -> infer(hybrid, cells)))
        @test occursin("does not apply here", msg(() -> infer(ode, obs; route = :mixed)))
    end

    # The identity checks run each route twice at smoke size, about 5 s per
    # gradient, so they are past the cost guard of §12 2026-10-24 and run with
    # the integration tests.
    if get(ENV, "INFERCELL_INTEGRATION_TESTS", "false") == "true"
        H = 120.0
        ds = m0_dataset(16085; ncells = 3, horizon = H)
        ms = collect(AbstractSubModel, m0_models())
        names = [:log_S_JCVISYN3A_0607, :log_S_JCVISYN3A_0779, :log_krnadeg,
                 :log_kcatF_R_ENO, :log_kcatF_R_FBA, :log_sigma_b]

        @testset "identity: :exact is the Gibbs chain composed by hand" begin
            N, n_sweeps, n_burn = 4, 2, 1
            ch = infer(ms, ds; route = :exact, rng = Xoshiro(6), N, n_sweeps, n_burn)

            rng = Xoshiro(rand(Xoshiro(6), UInt64))
            g = gibbs_setup(ds; models = ms, N, horizon = H)
            v = Dict(draw_truth(rng, ms; names = M0_TARGETS, purpose = :calibration).values)
            σ0 = rand(rng, SIGMA_MET_PRIOR)
            s = gibbs_init(rng, g; cme = [v[n] for n in g.cme], u = [log(v[f]) for f in g.forwards],
                           σ = σ0)
            for _ in 1:n_sweeps
                gibbs_sweep!(s, g)
            end
            byname = Dict(vcat(g.cme, g.forwards) .=> 1:(length(g.cme) + length(g.forwards)))
            want = [vcat(vcat(log.(t.cme), t.u)[[byname[n] for n in M0_TARGETS]], log(t.σ))'
                    for t in s.trace[(n_burn + 1):end]]

            @test ch.info.route === :exact
            @test Symbol.(ch.name_map.parameters) == names
            @test size(ch) == (n_sweeps - n_burn, length(names), 1)
            @test Array(ch[:, names, :])[:, :, 1] == reduce(vcat, want)
            @test ch.info.start.σ == σ0

            # The default rng. The setup reseeds the global generator per cell, so
            # without the route's own stream two seeds drew the same chain, and
            # without putting the global state back consecutive calls did.
            short = (; N = 2, n_sweeps = 1)
            Random.seed!(1)
            a = infer(ms, ds; route = :exact, short...)
            b = infer(ms, ds; route = :exact, short...)
            Random.seed!(2)
            c = infer(ms, ds; route = :exact, short...)
            @test Array(a) != Array(b) && Array(a) != Array(c) && Array(b) != Array(c)
            # The caller's global stream moves by the fork's one draw, nothing more.
            after = rand()
            Random.seed!(2)
            rand(UInt64)
            @test after == rand()
        end

        @testset "identity: :mixed is B5's pipeline composed by hand" begin
            kw = (n_ref = 8, n_adapt = 1, n_draws = 2, max_depth = 2, map_maxiters = 1)
            ch = infer(ms, ds; route = :mixed, rng = Xoshiro(7), kw...)

            s0 = mixed_setup(ds)
            ref = reference_theta(s0; maxiters = kw.map_maxiters)
            r = _mixed_solve(s0, exp.(ref.x[1:length(s0.names)]))
            cells = reference_ensemble(mixed_theta_values(s0, ref.x), reference_seeds(ds, kw.n_ref);
                                       horizon = H)
            d = delta_and_S(cells, mean_field_log_bulk(s0, r.sol); C = 3)
            s = with_reference(s0, d.delta, d.S)
            hand = sample(Xoshiro(rand(Xoshiro(7), UInt64)), mixed_model(s),
                          NUTS(kw.n_adapt, 0.8; max_depth = kw.max_depth), kw.n_draws;
                          initial_params = ref.x, progress = false)

            @test ch.info.route === :mixed
            @test Symbol.(ch.name_map.parameters) == names == mixed_parameter_names(s)
            @test ch.info.theta_ref.x == ref.x && ch.info.theta_ref.method === ref.method
            @test ch.info.delta == d.delta && ch.info.S == d.S
            @test Array(ch[:, names, :]) == Array(hand[:, [Symbol("x[$i]") for i in 1:6], :])
        end
    end
end
