# KIM-116: the modular arm on the pilot's dataset (Linear M9; Spec D7; Decision log
# 2026-11-19). The reduced model is inferred one block at a time. Stage 1 fits the
# promoters and krnadeg to the transcripts. Stage 2 fits ENO, FBA and σ_b to the bulk
# panel, with stage 1's draws fixed. Stage 1 comes in two variants:
#   cut     the factor cut, 1(b): the transcript factor alone over all 21 coordinates,
#           with ENO, FBA and σ_b under their priors. It still solves route (ii);
#   plugin  the data plug-in, 2: the transcript factor at rates rebuilt from the
#           observed panel (ATP, GTP; CTP and UTP are chemostats) over the 18 stage-1
#           coordinates. The states come from the solve at KIM-92's θ_ref, with the
#           panel pools replaced and the volume taken from θ_ref. It solves nothing.
# The joint posterior, KIM-92's, is the reference (archive infercell-kim92-archive-2026-11-17).
#
# Usage: julia --project dev/scripts/corea_modular_kim116.jl <stage> <out>
# Submit through dev/scripts/corea_modular_kim116.sh, from a worktree at a pushed commit.
# Stages:
#   profile <out>   item 2: each stage's log density and gradient, median of 3, with time,
#                   allocation and peak RSS, beside the joint's on the same node. It runs
#                   at the truth and the prior median, as KIM-92's profile did, with the
#                   fit's δ and S. It also reports the plug-in's inputs: its rates against
#                   route (ii)'s at the truth, and the volume factor at θ_ref against 20
#                   joint draws. Writes profile.jls and profile.md.
# SMOKE=1: the same calls on a two-cell, two-save dataset, without the archive. A
# plumbing pass, never a result.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: _mixed_solve, _state_at, _route_geometry
using ForwardDiff
using Printf
using Statistics: median, quantile
using DelimitedFiles: readdlm

const ARCHIVE = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const SEED = 17_200_001
const NREPEAT = 3                                 # item 2: median of 3
const JOINT_GRADIENT_KIM92 = 40.3                 # s, Slurm 20541310, δ = S = 0
const NFACTOR_DRAWS = 20
const STAGE2 = [:log_kcatF_R_ENO, :log_kcatF_R_FBA, :log_sigma_b]

pilot_setup(ds; kw...) = mixed_setup(ds; models = d11_models(), targets = COREA_TARGETS, kw...)

function truth_x(ds, s)
    tv = Dict(ds.truth.values)
    return vcat([log(tv[n]) for n in s.names], log(ds.sigma))
end

# The pilot's dataset, δ, S and θ_ref; a SMOKE run makes a tiny stand-in.
function inputs()
    if SMOKE
        ds = corea_dataset(SEED; ncells = 2, horizon = 120.0, targets = COREA_TARGETS)
        s = pilot_setup(ds)
        n = length(s.logz)
        return (; ds, delta = zeros(n), S = zeros(n, n), x_ref = mixed_prior_median(s),
                draws = repeat(permutedims(mixed_prior_median(s)), 2))
    end
    ds = load(ARCHIVE, "dataset")
    p = load(ARCHIVE, "prep")
    draws, header = readdlm(joinpath(ARCHIVE, "pilot_draws.tsv"), '\t', Float64; header = true)
    vec(header) == string.(mixed_parameter_names(pilot_setup(ds))) ||
        error("pilot_draws.tsv's columns are not the setup's coordinates")
    return (; ds, delta = p.delta, S = p.S, x_ref = p.ref.x, draws)
end

# The four densities, each as a function of its own free coordinates, with its start.
function stages(ds, delta, S, x_ref, x)
    both = with_reference(pilot_setup(ds), delta, S)
    txo = with_reference(pilot_setup(ds; streams = (:tx,)), delta, S)
    bulko = with_reference(pilot_setup(ds; streams = (:bulk,)), delta, S)
    names = mixed_parameter_names(both)
    t_states = @elapsed states = panel_plugin_states(txo, x_ref)
    plug = modular_stage(txo, setdiff(names, STAGE2), x; loglik = (s, y) -> plugin_loglik(s, states, y))
    two = modular_stage(bulko, STAGE2, x)
    return (; both, txo, states, t_states,
            calls = [:joint => (z -> mixed_logdensity(both, z), x),
                     :stage1_cut => (z -> mixed_logdensity(txo, z), x),
                     :stage1_plugin => (z -> stage_logdensity(plug, z), stage_x0(plug)),
                     :stage2 => (z -> stage_logdensity(two, z), stage_x0(two))])
end

# `f()` run NREPEAT times: seconds, GiB allocated, GC seconds, and peak RSS after.
function timed(f)
    map(1:NREPEAT) do _
        st = @timed f()
        (; seconds = st.time, gib = st.bytes / 2^30, gc = st.gctime, maxrss_gib = Sys.maxrss() / 2^30,
         value = st.value)
    end
end
strip_value(v) = [Base.structdiff(r, (; value = 0)) for r in v]

# The plug-in's rates against route (ii)'s at x, as |Δ log k| over genes and switches.
function rate_error(txo, states, x)
    θ = exp.(x[1:length(txo.names)])
    r = _mixed_solve(txo, θ)
    k_route = route_transcription_rates(txo, r.sol, r.p)
    k_plug = transcription_rates_from_states(txo, r.p, states)
    d = abs.(log.(k_plug) .- log.(k_route))
    return (; median = median(d), q90 = quantile(vec(d), 0.9), max = maximum(d),
            first_column_max = maximum(d[:, 1]))
end

# The particles-per-mM factor at each save, at x.
function volume_factors(s, x)
    sol = _mixed_solve(s, exp.(x[1:length(s.names)])).sol
    rt = s.route
    return [_route_geometry(rt, _state_at(sol, t)[rt.n_idx]).factor for t in s.saves]
end

function run_profile(out)
    mkpath(out)
    inp = inputs()
    s0 = pilot_setup(inp.ds)
    xs = [:truth => truth_x(inp.ds, s0), :prior_median => mixed_prior_median(s0)]
    # Compile every call on a two-cell, two-save case, so no timed call pays for it.
    t_compile = @elapsed begin
        tiny = corea_dataset(SEED; ncells = 2, horizon = 120.0, targets = COREA_TARGETS)
        st = pilot_setup(tiny)
        n = length(st.logz)
        for (_, (f, z)) in stages(tiny, zeros(n), zeros(n, n), mixed_prior_median(st),
                                  mixed_prior_median(st)).calls
            f(z)
            ForwardDiff.gradient(f, z)
        end
    end
    say(@sprintf("profile: compiled in %.0f s; %d bulk rows, %d saves to %.0f s", t_compile,
                 length(s0.logz), length(s0.saves), s0.saves[end]))
    rows = []
    local built
    for (label, x) in xs
        built = stages(inp.ds, inp.delta, inp.S, inp.x_ref, x)
        say(@sprintf("plug-in states built in %.1f s (one solve at θ_ref)", built.t_states))
        for (stage, (f, z)) in built.calls
            d = timed(() -> f(z))
            g = timed(() -> ForwardDiff.gradient(f, z))
            say(@sprintf("%s at %s, %d coordinates: density %s s, gradient %s s; log density %.6f",
                         stage, label, length(z), join((@sprintf("%.2f", r.seconds) for r in d), ", "),
                         join((@sprintf("%.2f", r.seconds) for r in g), ", "), d[1].value))
            push!(rows, (; stage, label, ncoord = length(z), density = strip_value(d),
                         gradient = strip_value(g), logdensity = d[1].value,
                         gradnorm = sqrt(sum(abs2, g[1].value))))
        end
    end
    # The plug-in's inputs: its rates at the truth, and the volume it borrows from θ_ref.
    re = rate_error(built.txo, built.states, xs[1][2])
    say(@sprintf("plug-in rates at the truth vs route (ii)'s: |Δ log k| median %.3g, 90%% %.3g, max %.3g",
                 re.median, re.q90, re.max))
    f_ref = volume_factors(s0, inp.x_ref)
    pick = round.(Int, range(1, size(inp.draws, 1); length = min(NFACTOR_DRAWS, size(inp.draws, 1))))
    dev = [maximum(abs.(log.(volume_factors(s0, inp.draws[i, :])) .- log.(f_ref))) for i in pick]
    say(@sprintf("volume factor: θ_ref vs %d joint draws, max |Δ log f| over saves: median %.3g, max %.3g",
                 length(pick), median(dev), maximum(dev)))
    save(out, "profile", (; rows, compile_seconds = t_compile, rate_error = re, factor_ref = f_ref,
                          factor_dev = dev, factor_draws = pick, t_states = built.t_states);
         stage = "profile")
    write_profile_md(out, rows, re, dev, built.t_states, s0)
end

function write_profile_md(out, rows, re, dev, t_states, s)
    io = IOBuffer()
    p(xs...) = println(io, xs...)
    p("# KIM-116 item 2: each stage's density and gradient\n")
    p("Commit `", COMMIT, "`, job ", JOB, ", on ", RAN_ON.host, " (", RAN_ON.cpu, "), ",
      Threads.nthreads(), " thread(s). SMOKE = ", SMOKE, ". The pilot's dataset (KIM-92, seed ",
      SEED, "), ", length(s.saves), " saves to ", s.saves[end], " s, ", length(s.logz),
      " bulk rows, with the fit's δ and S from the archive. Each call ran ", NREPEAT,
      " times after compiling on a two-cell case. The joint row is measured on the same node. ",
      "KIM-92's profile (Slurm 20541310, δ = S = 0) measured ", JOINT_GRADIENT_KIM92,
      " s for the joint gradient at the truth.\n")
    p("| stage | x | coordinates | density s (median) | gradient s (each) | gradient s (median) | GiB allocated, gradient (median) | peak RSS GiB |")
    p("|---|---|---|---|---|---|---|---|")
    for r in rows
        gs = [q.seconds for q in r.gradient]
        p(@sprintf("| %s | %s | %d | %.2f | %s | %.2f | %.2f | %.2f |", r.stage, r.label, r.ncoord,
                   median([q.seconds for q in r.density]), join((@sprintf("%.2f", t) for t in gs), ", "),
                   median(gs), median([q.gib for q in r.gradient]),
                   maximum(q.maxrss_gib for q in vcat(r.density, r.gradient))))
    end
    p("\nThe plug-in's states took ", @sprintf("%.1f", t_states), " s to build: one solve at θ_ref, once per fit.\n")
    p("## The plug-in's inputs\n")
    p(@sprintf("- Its rates at the truth, against route (ii)'s own: |Δ log k| over genes × switches has median %.3g, 90th percentile %.3g, max %.3g (column 1, t = 0, max %.3g). This is the panel's noise, δ, the 1 s offset and the borrowed volume together.",
               re.median, re.q90, re.max, re.first_column_max))
    p(@sprintf("- The volume factor at θ_ref against %d joint draws: max over saves of |Δ log f| has median %.3g, max %.3g.",
               length(dev), median(dev), maximum(dev)))
    write(joinpath(out, "profile.md"), String(take!(io)))
    print(read(joinpath(out, "profile.md"), String))
end

function main()
    length(ARGS) >= 2 || error("usage: corea_modular_kim116.jl <stage> <out>")
    stage, out = ARGS[1], ARGS[2]
    stage == "profile" ? run_profile(out) : error("unknown stage $stage")
end

main()
