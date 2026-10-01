# spec/phases/16-recovery.md task 16a.6 (parent 16.3) on Core A′: under a fixed
# recorded path, a clipped drain still puts a kink in the ODE state as a
# function of an ODE parameter, here ENO's forward constant (block 1's), with its
# reverse constant derived. Check 7's census finds the clip. The cell is the
# first of 15.7's cells, in seed order, at its truth, whose run clips at all. Variants, all replaying that one path:
# - :published, the clamp under fractional carry;
# - :kink, the clamp on continuous pools;
# - :sampler, 14c.5's smoothed drain on continuous pools.
#
# Usage: julia --project dev/scripts/fixed_path_kink_16a6.jl
# Regenerate: sbatch dev/scripts/fixed_path_kink_16a6.slurm

using InferCell
using Random
using Printf

const SEEDS = 30_001:30_200
const θNAME = :kcatF_R_ENO
models = d11_models()
smodels = d11_models(smoothing = COREA_SAMPLER.smoothing)
truth = draw_truth(Xoshiro(1507), models; purpose = :recovery)
tv = Dict(truth.values)

function build(v; horizon)
    ms = v === :sampler ? smodels : models
    d = build_problem(ms; tspan = (0.0, horizon), complete = true,
                      rounding = v === :published ? :fractional_carry : :continuous)
    return ms, d
end

# Record each cell's path in seed order until one clips: the first handshake at
# which any consumer drain clipped (check 7's predicate), and which drain.
SEED, nstar, kstar = 0, 0, 0
nchecked = 0
local_d = nothing
for seed in SEEDS
    global SEED, nstar, kstar, nchecked, local_d
    ms0, d0 = build(:published; horizon = COREA_CYCLE_S)
    Random.seed!(seed)
    set_parameters!(d0, ms0, truth.values)
    Random.seed!(seed)
    record_path!(d0)
    nchecked += 1
    for k in 1:round(Int, COREA_CYCLE_S)
        n = d0.n_clipped
        handshake_step!(d0)
        if d0.n_clipped > n
            global SEED, nstar = seed, k
            global kstar = findfirst(b -> b.clipped, d0.debits)
            break
        end
    end
    if nstar > 0
        local_d = d0
        break
    end
end
if nstar == 0
    println("None of the $nchecked cells of 15.7 clips within the cycle at its truth: ",
            "check 7's census is zero across the dataset, so no scan.")
    exit(0)
end
println("$nchecked cell(s) checked in seed order; the first to clip is $SEED.")
d = local_d
ms = models
path = recorded_path(d)
b = d.debits[kstar]
@printf("Cell %d at 15.7's truth: the first clip in check 7's census is %s on %s at handshake %d.\n",
        SEED, b.counter, b.species, nstar)
inames = InferCell._block_names(ms, :ode)
ipool = b.pool_idx

function after(θ, v)
    ms, e = build(v; horizon = Float64(nstar))
    writes = derived_ode_values(ms, [θNAME => θ])
    set_parameters!(e, ms, vcat([p for p in truth.values if !(first(p) in first.(writes))], writes))
    r = PathReplay(path, e)
    for _ in 1:nstar
        replay_step!(e, r)
    end
    return e.ode.u[ipool] * e.factor, e.debits[kstar].clipped
end
θ0 = tv[θNAME]
last(after(θ0, :published)) || error("the replay did not reproduce the clip")
c0 = last(after(θ0, :kink))
@printf("At the truth's %s = %.6g the clamp on continuous pools %s at that drain.\n",
        θNAME, θ0, c0 ? "also clips" : "does not clip")

# Bracket and bisect the θ at which that drain's clip flips, on :kink.
lo, hi = θ0, θ0
found = false
for f in (2.0, 0.5)
    global lo, hi, found
    a = θ0
    for _ in 1:20
        b2 = a * f
        if last(after(b2, :kink)) != c0
            lo, hi = minmax(a, b2)
            found = true
            break
        end
        a = b2
    end
    found && break
end
found || (println("No ENO constant within 2^±20 of the truth flips that drain's clip: ",
                  "ENO does not reach this pool strongly enough. Reported, not scanned."); exit(0))
clo = last(after(lo, :kink))
for _ in 1:60
    global lo, hi
    mid = sqrt(lo * hi)
    last(after(mid, :kink)) == clo ? (lo = mid) : (hi = mid)
end
θc = sqrt(lo * hi)
side = clo ? hi : lo                                  # the paying side
slope(θ, v, δ) = (first(after(θ * exp(δ), v)) - first(after(θ * exp(-δ), v))) / 2δ
dP = slope(side * (side > θc ? 1.5 : 1 / 1.5), :kink, 1e-3)
δ = 0.01 / abs(dP)
jump(k, v) = abs(slope(θc * exp(k * δ), v, δ) - slope(θc * exp(-k * δ), v, δ)) / abs(dP)
@printf("That drain's clip flips at %s = %.6g (%.4g × the truth). The pool's sensitivity on the paying side is %.5g particles per unit ln θ; step %.3g in ln θ.\n\n",
        θNAME, θc, θc / θ0, dP, δ)
println("| spacing either side (steps) | clamp, continuous pools | sampler model |")
println("|---|---|---|")
for k in (100, 30, 10, 3)
    @printf("| %d | %.4g | %.4g |\n", k, jump(k, :kink), jump(k, :sampler))
end
grid = [θc * exp(k * 37δ) for k in -10:10]
@printf("\nThe published model on %d points across the clip: zero finite-difference slope at %d.\n",
        length(grid), count(θ -> slope(θ, :published, δ) == 0, grid))
println("Done")
