# spec/phases/16-recovery.md task 16a.1, V2 at full scale: three cells of the
# 15.7 dataset (its truth, seeds 30001 to 30003, 6,300 s), each simulated with
# its path recorded, then replayed from that path. Reports, per cell:
# - whether the recorded run's 60 s latent equals generate_dataset's (recording
#   is transparent on the dataset's own configuration);
# - whether the replay equals the run bitwise at every one of 6,300 handshakes;
# - the simulation's and the replay's warm wall-clock per cycle, against 13.7's
#   32.2 s;
# - the size and cost of one mid-cycle snapshot.
#
# Usage: julia --project dev/scripts/path_replay_16a1.jl
# Regenerate: sbatch dev/scripts/path_replay_16a1.slurm

using InferCell
using Random
using Printf

const TRUTH_SEED = 1507                      # 15.7's truth
const CELLS = 30_001:30_003
const HORIZON = COREA_CYCLE_S
const N = round(Int, HORIZON)

state(d) = (t = (d.ode.t, d.jump.t), ode = copy(d.ode.u), ode_p = copy(d.ode.p),
            jump = copy(d.jump.u), jump_p = copy(d.jump.p),
            carry = copy(d.rounding.remainders),
            deficit = [b.deficit for b in d.debits],
            exchanged = [b.exchanged for b in d.debits],
            geometry = (d.factor, d.area_nm2, d.radius_nm, d.volume_litres),
            census = (d.n_handshakes, d.n_drains, d.n_clipped))

models = d11_models()
truth = draw_truth(Xoshiro(TRUTH_SEED), models; purpose = :recovery)
function cell(; horizon = HORIZON)
    d = build_problem(models; tspan = (0.0, horizon), complete = true)
    set_parameters!(d, models, truth.values)
    return d
end
ode_names = InferCell._block_names(models, :ode)
jump_names = InferCell._block_names(models, :jump)

# Warm-up, so no timing below includes compilation.
let d = cell(; horizon = 120.0)
    record_path!(d)
    for _ in 1:120; handshake_step!(d); end
    p = recorded_path(d)
    e = cell(; horizon = 120.0)
    replay!(e, p, 120)
    restore(snapshot(e))
end

println("16a.1 V2 at full scale: truth Xoshiro($TRUTH_SEED), cells $(CELLS), $N handshakes")
for seed in CELLS
    ds = generate_dataset(truth, [seed]; models, horizon = HORIZON)

    Random.seed!(seed)
    d = cell()
    Random.seed!(seed)
    record_path!(d)
    sim = Vector{Any}(undef, N)
    t_sim = @elapsed for k in 1:N
        handshake_step!(d)
        sim[k] = state(d)
    end
    path = recorded_path(d)

    # The recorded run's 60 s latent, as emit_observables! would record it.
    rows = [findfirst(==(s), vcat(ode_names, jump_names)) for s in ds.species]
    at(k) = vcat(sim[k].ode .* sim[k].geometry[1] .+ sim[k].carry, Float64.(sim[k].jump))[rows]
    latent = reduce(hcat, [at(k) for k in 60:60:N])
    same_latent = latent == ds.latent[1, :, 2:end]

    e = cell()
    r = PathReplay(path, e)
    rep = Vector{Any}(undef, N)
    t_rep = @elapsed for k in 1:N
        replay_step!(e, r)
        rep[k] = state(e)
    end
    first_bad = findfirst(k -> rep[k] != sim[k], 1:N)

    f = cell()
    replay!(f, path, N ÷ 2)
    t_snap = @elapsed s = snapshot(f)
    t_rest = @elapsed restore(s)

    @printf("seed %d: %d firings; 60 s latent equals generate_dataset: %s; replay bitwise at all %d handshakes: %s%s\n",
            seed, length(path), same_latent, N, first_bad === nothing,
            first_bad === nothing ? "" : " (first differs at handshake $first_bad)")
    @printf("  wall-clock per cycle: simulate %.2f s, replay %.2f s (ratio %.3f)\n",
            t_sim, t_rep, t_rep / t_sim)
    @printf("  mid-cycle snapshot: %.1f KiB, snapshot %.2f ms, restore %.2f ms\n",
            Base.summarysize(s) / 1024, 1e3 * t_snap, 1e3 * t_rest)
end
println("Done")
