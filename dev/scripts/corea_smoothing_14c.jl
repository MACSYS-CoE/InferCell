# Spec §11 tasks 14c.5 and 14c.7: the sampler model against the published one,
# on the assembled Core A′. The sampler model smooths every consumer counter at
# COREA_SMOOTHING_WIDTH and writes pools back continuously (COREA_SAMPLER); the
# published model clamps under fractional carry (§12, 2026-09-28). Two
# sections, one Slurm array:
#
# - agree: 14b's ten published-parameter seeds, each run through a full cycle
#   under both models, with every candidate observable recorded at 60 s.
#   Candidate observables are task 15.3's: every ODE state as a particle count
#   (with its carried remainder, so both models report the amount they hold),
#   every transcript count, and the cell volume. Fluxes are left out, since §4
#   D8 rules them out of the likelihood. The merge pairs the two runs per seed.
# - deriv: the continuity half on the assembled model. At the first drain where
#   `GTP_translat` clips at seed 14800, the GTP state after the debit — what the
#   next interval integrates from — is scanned in PGK3's forward constant, the
#   ODE parameter that supplies GTP, across the value at which that drain stops
#   clipping. Three variants are applied to the same saved state: the published
#   model, the clamp on continuous pools, and the sampler model. The suite
#   asserts the same on a toy (test/test_corea_smoothed_drain.jl); this is the
#   record on the real model.
#
# Usage: julia --project dev/scripts/corea_smoothing_14c.jl <task> <ntasks>
#        Task 0 also runs deriv. Merge with corea_smoothing_14c_merge.jl.
# Regenerate: sbatch dev/scripts/corea_smoothing_14c.slurm

using InferCell
using Printf
using Random
using Dates

const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, "1"))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const CYCLE = SMOKE ? 120 : round(Int, COREA_CYCLE_S)
const SEEDS = SMOKE ? [14_800] : vcat([1310, 1410], 14_800 .+ (0:7))   # 14b's
const EVERY = 60
const DIR = joinpath(@__DIR__, "smoothing_14c")
const HEADER = "# commit $(strip(read(`git rev-parse --short HEAD`, String))), job " *
               "$(get(ENV, "SLURM_JOB_ID", "none")), task $TASK of $NTASK, $(now()), " *
               "Julia $VERSION$(SMOKE ? ", SMOKE" : "")"
mkpath(DIR)

_names(ms, f) = Symbol[s for m in ms if formalism(m) === f for s in states(m)]

# ---------------------------------------------------------------------------
# agree
# ---------------------------------------------------------------------------

function observe(d, ode_names, jump_names, tx)
    vals = Pair{Symbol, Float64}[]
    for (i, s) in enumerate(ode_names)
        push!(vals, s => d.ode.u[i] * d.factor + d.rounding.remainders[i])
    end
    for (j, s) in enumerate(jump_names)
        s in tx && push!(vals, s => float(d.jump.u[j]))
    end
    push!(vals, :volume_litres => d.volume_litres)
    return vals
end

function agree_run(io, seed, drain)
    # `:smoothed` labels the sampler model, which smooths the drain and carries
    # pools continuously; `:clamped` is the published model (spec §12, 2026-09-28).
    kw = drain === :smoothed ? COREA_SAMPLER : (;)
    ms = corea_models(; smoothing = get(kw, :smoothing, nothing))
    ode_names, jump_names = _names(ms, :ode), _names(ms, :jump)
    tx = Set(transcript_state(g.locus) for g in read_transcription_genes())
    Random.seed!(seed)
    d = build_corea(; kw..., tspan = (0.0, Float64(CYCLE)))
    Random.seed!(seed)
    rows(t) = for (s, v) in observe(d, ode_names, jump_names, tx)
        println(io, join((seed, drain, t, s, repr(v)), '\t'))
    end
    rows(0)
    wall = @elapsed for k in 1:CYCLE
        handshake_step!(d)
        k % EVERY == 0 && rows(k)
    end
    c = clipping_census(d)
    println(io, "# census\t$seed\t$drain\t$(c.drains)\t$(c.clipped)\t$(@sprintf("%.1f", wall))")
    flush(io)
end

open(joinpath(DIR, "agree_task_$(TASK).tsv"), "w") do io
    println(io, HEADER)
    println(io, join(("seed", "drain", "t", "observable", "value"), '\t'))
    for (j, seed) in enumerate(SEEDS)
        (j - 1) % NTASK == TASK || continue
        for drain in (:clamped, :smoothed)
            println("$(now()) agree seed $seed $drain"); flush(stdout)
            agree_run(io, seed, drain)
        end
    end
end

# ---------------------------------------------------------------------------
# deriv
# ---------------------------------------------------------------------------

TASK == 0 && open(joinpath(DIR, "deriv.md"), "w") do io
    p(args...) = (println(io, args...); println(args...); flush(io))
    p(HEADER[3:end])
    p()
    seed = 14_800
    θname = :kcatF_R_PGK3
    ms = corea_models()
    ms[3] = NucleotideRecycling(enzymes = :translated, free = [θname])
    ode_names = unique(Symbol[q.name for m in ms if formalism(m) === :ode
                              for q in model_free_params(parameters(m))])
    iθ = findfirst(==(θname), ode_names)
    igtp = findfirst(==(:M_gtp_c), _names(ms, :ode))
    gtp_debit(d) = only(k for (k, b) in enumerate(d.debits)
                        if b.counter === :GTP_translat && b.sign < 0)

    fresh() = (Random.seed!(seed);
               d = build_problem(ms; tspan = (0.0, Float64(CYCLE)), complete = true);
               Random.seed!(seed); d)

    # The first drain at which GTP_translat clips.
    d = fresh()
    kg = gtp_debit(d)
    nstar = 0
    for k in 1:CYCLE
        handshake_step!(d)
        d.debits[kg].clipped && (nstar = k; break)
    end
    if nstar == 0
        p("`GTP_translat` did not clip within $CYCLE handshakes at seed $seed; no scan.")
        return
    end
    d0 = fresh()
    for _ in 1:(nstar - 1)
        handshake_step!(d0)
    end
    θ0 = d0.ode.p[iθ]

    # One more handshake from the saved state at θ, returning the GTP state after
    # the debit in particles — what the next interval integrates from, without
    # the carried remainder no rate law reads — and whether that debit clipped.
    # Variants: `:published` as built; `:kink`, the clamp on continuous pools;
    # `:sampler`, the smoothed drain on continuous pools.
    function after(θ, variant)
        d = deepcopy(d0)
        d.ode.p[iθ] = θ
        if variant !== :published
            d.rounding = RoundingState(:continuous; nspecies = length(d.ode.u))
        end
        if variant === :sampler
            for b in d.debits
                b.sign < 0 && b.pool_idx != 0 &&
                    (b.clip = :smoothed; b.smoothing = COREA_SMOOTHING_WIDTH)
            end
        end
        handshake_step!(d)
        return d.ode.u[igtp] * d.factor, d.debits[kg].clipped
    end
    @assert last(after(θ0, :published)) "the replay did not reproduce the clip at handshake $nstar"

    # Bracket the value at which the drain stops clipping, then bisect, on the
    # clamp over continuous pools, whose clip is a point rather than a step.
    lo, hi = θ0, θ0
    while last(after(hi, :kink))
        lo = hi
        hi *= 2
        hi > 1e6 * θ0 && error("no θ up to 1e6 × nominal stops the clip")
    end
    for _ in 1:60
        mid = sqrt(lo * hi)
        last(after(mid, :kink)) ? (lo = mid) : (hi = mid)
    end
    θc = sqrt(lo * hi)

    # Slopes in ln θ. δ is set so the pool moves about 0.01 particles per step,
    # well inside the one-particle width, from the slope on the paying side.
    slope(θ, v, δ) = (first(after(θ * exp(δ), v)) - first(after(θ * exp(-δ), v))) / 2δ
    dP = slope(θc * 1.5, :kink, 1e-3)
    δ = 0.01 / abs(dP)
    jump(k, v) = abs(slope(θc * exp(k * δ), v, δ) - slope(θc * exp(-k * δ), v, δ)) / abs(dP)
    p("## 14c.5 and 14c.7: the derivative across a clip, on the assembled model")
    p()
    p(@sprintf("Seed %d: `GTP_translat` first clips at handshake %d (t = %d s). `%s` is %.6g ",
               seed, nstar, nstar, θname, θ0),
      @sprintf("at nominal, and that drain stops clipping at %.6g (%.4g × nominal). ", θc, θc / θ0),
      @sprintf("The pool's sensitivity there, dP/d ln θ on the paying side, is %.5g particles.", dP))
    p()
    p("The slope, in ln θ, of the GTP state after the debit, either side of the clip, as a ",
      @sprintf("fraction of that sensitivity. Finite-difference step %.3g in ln θ. ", δ),
      "All three variants start from the published model's saved state.")
    p()
    p("| spacing either side (steps) | pool − accrual at the spacing (particles) | clamp, continuous pools | sampler model |")
    p("|---|---|---|---|")
    for k in (100, 30, 10, 3)
        p(@sprintf("| %d | ±%.3g | %.4g | %.4g |", k, k * δ * abs(dP), jump(k, :kink), jump(k, :sampler)))
    end
    p()
    p("A continuous derivative has a jump that falls with the spacing. The clamp's ",
      "stays at the whole sensitivity, since one side pays and the other floors at zero.")
    p()
    grid = [θc * exp(k * 37δ) for k in -10:10]
    ys = [first(after(θ, :published)) for θ in grid]
    nzero = count(θ -> slope(θ, :published, δ) == 0, grid)
    p(@sprintf("The published model on the same %d points, spaced %d steps apart across the clip: ",
               length(grid), 37),
      "the state is a whole number of particles at ", count(isinteger, ys), " of them, ",
      "and the finite-difference slope is exactly zero at $nzero. Its derivative is zero ",
      "wherever no one-particle step falls inside the stencil, so a gradient through the ",
      "handshake sees no dependence on θ at all.")
end
println("Done")
