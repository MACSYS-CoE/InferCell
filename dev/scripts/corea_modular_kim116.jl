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
# The design, gates and comparison are frozen in spec/prereg/kim-116.md at tag
# prereg/kim-116. This script implements them and changes none.
#
# Usage: julia --project dev/scripts/corea_modular_kim116.jl <stage> <out> [<args>]
# Submit through dev/scripts/corea_modular_kim116.sh, from a worktree at a pushed commit.
# <v> is a variant: `cut` or `plugin`. Its files go to <out>/<v>/.
# Stages:
#   profile <out>       item 2: each stage's log density and gradient, median of 3, with
#                       time, allocation and peak RSS, beside the joint's on the same node.
#                       It runs at the truth and the prior median, as KIM-92's profile did,
#                       with the fit's δ and S. It also reports the plug-in's inputs: its
#                       rates against route (ii)'s at the truth, and the volume factor at
#                       θ_ref against 20 joint draws. Writes profile.jls and profile.md;
#   metric <out> v      the warm metric at θ_ref over v's stage-1 coordinates, by KIM-92's
#                       rule (Decision log 2026-11-12);
#   chain <out> v c     stage-1 chain c of 4: NUTS(NADAPT, 0.8), max depth 10, from the
#                       warm metric, jittered 0.25 prior SD about θ_ref, checkpointed;
#   check <out> v k     the stage-1 gate (prereg §4) on the first k kept draws per chain;
#   extend <out> v c    the one extension: NDRAW more draws, no new adaptation;
#   stage2 <out> v a b  the inner chains a..b of K (prereg §3): ψ fitted to the bulk at
#                       the k-th evenly spaced pooled stage-1 draw. One file per chain,
#                       and a chain whose file exists is skipped;
#   ppc <out> v         KIM-92's PPC on 200 cut draws × 10 (reported, prereg §5);
#   merge <out>         the comparison with the joint (prereg §5), to modular_summary.md,
#                       with the cut draws as cut_draws_<v>.tsv in pilot_draws.tsv's columns.
# SMOKE=1: a two-cell, two-save dataset, no archive, and tiny chains. A plumbing pass,
# never a result.

include(joinpath(@__DIR__, "pilot_common.jl"))
using InferCell: _mixed_solve, _state_at, _route_geometry
using ForwardDiff
using Printf
using Random
using LinearAlgebra
using Statistics: median, quantile, std, cor
using DelimitedFiles: readdlm
using Turing
using Turing: NUTS
using MCMCChains
const DPPL = Turing.DynamicPPL
const AHMC = Turing.Inference.AHMC
const LDP = Turing.Inference.LogDensityProblems
const AMCMC = Turing.Inference.AbstractMCMC
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

envint(k, default) = (v = strip(get(ENV, k, "")); isempty(v) ? default : parse(Int, v))

# prereg §3 and §4
const VARIANTS = ("cut", "plugin")
const NCHAIN = 4
const NADAPT = envint("NADAPT", SMOKE ? 5 : 200)
const NDRAW = envint("NDRAW", SMOKE ? 5 : 200)
const MAX_DEPTH = SMOKE ? 3 : 10
const TARGET_ACCEPT = 0.8
const JITTER_SD = 0.25
const K2 = SMOKE ? 4 : 200
const NADAPT2 = SMOKE ? 5 : 100
const NKEEP2 = SMOKE ? 2 : 10
const RHAT_MAX, ESS_MIN = 1.02, 100
const CLOSE_SHIFT, CLOSE_CONTRACTION = 0.5, 0.05
const PPC_DRAWS, PPC_REPS = SMOKE ? (4, 2) : (200, 10)
const PPC_LO, PPC_HI = 0.005, 0.995
const CI_LO, CI_HI = 0.05, 0.95
const FD_STEP = 1e-3                              # × the coordinate's prior SD
const CKPT_EVERY = envint("CKPT_EVERY", SMOKE ? 2 : 5)
const SEED_BASE = 95_000_000_000                  # prereg §3

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

# ---------------------------------------------------------------- the fit's pieces

const PHI_NAMES = [Symbol(:log_, n) for n in COREA_TARGETS if !(n in (:kcatF_R_ENO, :kcatF_R_FBA))]

# The setups, θ_ref, and the plug-in's states; the same in every stage.
function fit()
    inp = inputs()
    ds = inp.ds
    txo = with_reference(pilot_setup(ds; streams = (:tx,)), inp.delta, inp.S)
    bulko = with_reference(pilot_setup(ds; streams = (:bulk,)), inp.delta, inp.S)
    both = with_reference(pilot_setup(ds), inp.delta, inp.S)
    names = mixed_parameter_names(both)
    phi = [findfirst(==(n), names) for n in PHI_NAMES]
    psi = [findfirst(==(n), names) for n in STAGE2]
    length(phi) == 18 && length(psi) == 3 && sort(vcat(phi, psi)) == collect(eachindex(names)) ||
        error("φ and ψ do not partition the coordinates")
    return merge(inp, (; txo, bulko, both, names, phi, psi, x_ref = Vector{Float64}(inp.x_ref)))
end

# Stage 1 of variant v, with its states built when it is the plug-in.
function stage1(f, v)
    v == "cut" && return modular_stage(f.txo, eachindex(f.names), f.x_ref)
    v == "plugin" || error("unknown variant $v")
    states = panel_plugin_states(f.txo, f.x_ref)
    return modular_stage(f.txo, f.phi, f.x_ref; loglik = (s, y) -> plugin_loglik(s, states, y))
end

# Where v's stage-1 φ sits among its free coordinates.
phi_in_free(f, st) = [findfirst(==(i), st.free) for i in f.phi]

vdir(out, v) = (d = joinpath(out, v); mkpath(d); d)
has(d, name) = isfile(joinpath(d, "$name.jls"))

# prereg §3's seeds, refused if they meet the dataset's cells, its reference cells,
# or any of KIM-92's chain arms.
function seed_of(f, v, offset)
    s = SEED_BASE + 64 * SEED + 1000 * (v == "cut" ? 1 : 2) + offset
    ds = f.ds
    taken = Set(vcat(ds.seeds, SMOKE ? Int[] : reference_seeds(ds, 800),
                     reduce(vcat, [mixed_chain_seeds(ds, 16; arm) for arm in 0:3])))
    s in taken && error("seed $s meets one of the dataset's")
    return s
end
chain_seed(f, v, c) = seed_of(f, v, c)
extend_seed(f, v, c) = seed_of(f, v, 10 + c)
inner_seed(f, v, k) = seed_of(f, v, 100 + k)

nuts1() = NUTS(NADAPT, TARGET_ACCEPT; max_depth = MAX_DEPTH)

# A start `JITTER_SD` prior SDs about st.x in the free coordinates, clamped to the box,
# redrawn up to 20 times until the stage density is finite (mixed_jittered_start's rule).
function stage_start(rng, st)
    lo, hi = mixed_box(st.s)
    sd = mixed_prior_sd(st.s)[st.free]
    for k in 1:20
        z = clamp.(stage_x0(st) .+ JITTER_SD .* sd .* randn(rng, length(sd)), lo[st.free], hi[st.free])
        isfinite(stage_logdensity(st, z)) && return (z = z, tries = k)
    end
    error("no finite start in 20 draws")
end

function progress_callback(label, total; every = SMOKE ? 1 : 5, offset = 0, t0 = time())
    return function (rng, model, sampler, transition, state, i; kwargs...)
        i += offset
        (i % every == 0 || i == total) || return nothing
        st = transition.stat
        say(@sprintf("%s: iteration %d of %d, %.1f h; depth %d, %d steps, step size %.3g, lp %.2f",
                     label, i, total, (time() - t0) / 3600, st.tree_depth, st.n_steps,
                     st.step_size, transition.lp))
        return nothing
    end
end

# ---------------------------------------------------------------- warm metric (KIM-92's rule)

function hessian_fd(grad, x, h)
    P = length(x)
    H = zeros(P, P)
    for j in 1:P
        e = zeros(P); e[j] = h[j]
        H[:, j] = (grad(x .+ e) .- grad(x .- e)) ./ (2h[j])
        say(@sprintf("metric: column %d of %d", j, P))
    end
    asym = maximum(abs, H .- H') / maximum(abs, H)
    return (H .+ H') ./ 2, asym
end

function run_metric(out, v)
    f = fit()
    st = stage1(f, v)
    z = stage_x0(st)
    sd = mixed_prior_sd(st.s)[st.free]
    grad(y) = ForwardDiff.gradient(w -> stage_logdensity(st, w), y)
    t = @elapsed (H, asym) = hessian_fd(grad, z, FD_STEP .* sd)
    F = eigen(Symmetric(-H))
    floor_ = minimum(1 ./ sd .^ 2)
    λ = max.(F.values, floor_)
    nfloor = count(F.values .< floor_)
    minv = diag(F.vectors * Diagonal(1 ./ λ) * F.vectors')
    all(m -> isfinite(m) && m > 0, minv) || error("the metric is not positive: $minv")
    save(vdir(out, v), "metric", (; minv, H, asym, eigenvalues = F.values, floor = floor_, nfloor,
                                  fd_step = FD_STEP, seconds = t); stage = "metric")
    say(@sprintf("metric %s: %d coordinates in %.0f s; %d eigenvalues floored; asymmetry %.2g",
                 v, length(z), t, nfloor, asym))
end

# KIM-92's warm_state (corea_pilot_kim92.jl), on the stage's model.
function warm_state(rng, model, spl, x0, minv)
    vi = DPPL.default_varinfo(rng, model, spl)
    vi = DPPL.initialize_parameters!!(vi, x0, model)
    vi = last(DPPL.evaluate!!(model, vi, DPPL.DefaultContext()))
    vi = DPPL.link(vi, model)
    theta = vi[:]
    theta ≈ x0 || error("the linked start is not x0: the prior is no longer the identity's")
    ldf = DPPL.LogDensityFunction(model, vi,
                                  DPPL.SamplingContext(rng, spl, DPPL.leafcontext(model.context));
                                  adtype = spl.alg.adtype)
    h = AHMC.Hamiltonian(AHMC.DiagEuclideanMetric(copy(minv)), Base.Fix1(LDP.logdensity, ldf),
                         Base.Fix1(LDP.logdensity_and_gradient, ldf))
    z = AHMC.phasepoint(rng, theta, h)
    ϵ = AHMC.find_good_stepsize(rng, h, theta)
    kernel = Turing.Inference.make_ahmc_kernel(spl.alg, ϵ)
    adaptor = Turing.Inference.AHMCAdaptor(spl.alg, h.metric; ϵ)
    adaptor isa AHMC.StanHMCAdaptor || error("expected Stan's windowed adaptor, got $(typeof(adaptor))")
    adaptor.pc.var .= minv
    return Turing.Inference.HMCState(vi, 0, kernel, h, z, adaptor), ϵ
end

# A saved HMC state with its Hamiltonian rebuilt on `model` (KIM-92's fresh_state).
function fresh_state(rng, model, spl, st)
    ldf = DPPL.LogDensityFunction(model, st.vi,
                                  DPPL.SamplingContext(rng, spl, DPPL.leafcontext(model.context));
                                  adtype = spl.alg.adtype)
    h = AHMC.Hamiltonian(st.hamiltonian.metric, Base.Fix1(LDP.logdensity, ldf),
                         Base.Fix1(LDP.logdensity_and_gradient, ldf))
    z = AHMC.phasepoint(h, st.z.θ, st.z.r)
    return Turing.Inference.HMCState(st.vi, st.i, st.kernel, h, z, st.adaptor)
end

# ---------------------------------------------------------------- stage-1 chains

function run_chain(out, v, c)
    f = fit()
    st = stage1(f, v)
    d = vdir(out, v)
    seed = chain_seed(f, v, c)
    model, spl = modular_model(st), DPPL.Sampler(nuts1())
    total, ck = NADAPT + NDRAW, "chain_$(c)_ckpt"
    if has(d, ck)
        raw = deserialize(joinpath(d, "$ck.jls"))
        raw.commit_full == COMMIT_FULL || error("checkpoint $ck is from $(raw.commit), this is $COMMIT")
        r = raw.result
        rng, segs, start, tries, seconds, eps0 = copy(r.rng), r.segs, r.start, r.tries, r.seconds, r.eps0
        state = fresh_state(rng, model, spl, r.state)
        done = sum(size(s, 1) for s in segs)
        say(@sprintf("%s chain %d: resumed at iteration %d of %d", v, c, done, total))
    else
        rng = Xoshiro(seed)
        s0 = stage_start(rng, st)
        start, tries = s0.z, s0.tries
        state, eps0 = warm_state(rng, model, spl, start, load(d, "metric").minv)
        segs, done, seconds = MCMCChains.Chains[], 0, 0.0
        say(@sprintf("%s chain %d: seed %d, start after %d tries, log density %.4f; initial step size %.3g",
                     v, c, seed, tries, stage_logdensity(st, start), eps0))
    end
    t0 = time() - seconds
    while done < total
        n = min(CKPT_EVERY, total - done)
        t = @elapsed seg = AMCMC.mcmcsample(rng, model, spl, n; chain_type = MCMCChains.Chains,
                                            initial_state = state, nadapts = NADAPT,
                                            discard_initial = 0, progress = false, save_state = true,
                                            callback = progress_callback("$v chain $c", total;
                                                                         offset = done, t0))
        seconds += t
        state = seg.info.samplerstate
        push!(segs, MCMCChains.setrange(seg, (done + 1):(done + n)))
        done += n
        save(d, ck, (; c, seed, start, tries, eps0, segs, state, rng = copy(rng), seconds);
             stage = "chain checkpoint")
    end
    whole = reduce(vcat, segs)
    chain = MCMCChains.setinfo(whole, merge(whole.info, (; samplerstate = state)))
    size(chain, 1) == total || error("chain $c holds $(size(chain, 1)) draws, not $total")
    save(d, "chain_$c", (; c, seed, start, start_tries = tries, chain, seconds, eps0,
                         nadapt = NADAPT, ndraw = NDRAW); stage = "chain")
    say(@sprintf("%s chain %d: %d + %d iterations in %.1f h", v, c, NADAPT, NDRAW, seconds / 3600))
end

function run_extend(out, v, c)
    f = fit()
    st = stage1(f, v)
    d = vdir(out, v)
    base = load(d, "chain_$c")
    rng = Xoshiro(extend_seed(f, v, c))
    model, spl = modular_model(st), DPPL.Sampler(nuts1())
    state = fresh_state(rng, model, spl, base.chain.info.samplerstate)
    t = @elapsed chain = AMCMC.mcmcsample(rng, model, spl, base.ndraw;
                                          chain_type = MCMCChains.Chains, initial_state = state,
                                          nadapts = 0, discard_initial = 0, progress = false,
                                          save_state = true,
                                          callback = progress_callback("$v extend $c", base.ndraw))
    ε0, ε1 = base.chain[:step_size][end], unique(chain[:step_size])
    ε1 == [ε0] || error("the extension's step size $ε1 is not the chain's adapted $ε0")
    save(d, "extend_$c", (; c, chain, seconds = t, ndraw = base.ndraw); stage = "extend")
    say(@sprintf("%s extend %d: %d more iterations in %.1f h", v, c, base.ndraw, t / 3600))
end

chain_draws(ch, P; skip = 0) = Array(ch[(skip + 1):end, [Symbol("z[$i]") for i in 1:P], 1])

# Draws × free coordinates × chains, after adaptation, with the extension when it ran.
function stage1_draws(d, P; extended = all(c -> has(d, "extend_$c"), 1:NCHAIN))
    X = map(1:NCHAIN) do c
        r = load(d, "chain_$c")
        a = chain_draws(r.chain, P; skip = r.nadapt)
        extended ? vcat(a, chain_draws(load(d, "extend_$c").chain, P)) : a
    end
    return cat(X...; dims = 3), extended
end

pooled(X) = reshape(permutedims(X, (1, 3, 2)), :, size(X, 2))

function gate(X, idx)
    A = permutedims(X, (1, 3, 2))
    rh, es = MDT.rhat(A), MDT.ess(A; kind = :bulk)
    fail = [i for i in idx if !(rh[i] <= RHAT_MAX && es[i] >= ESS_MIN)]
    return (; rh, es, fail, pass = isempty(fail))
end

function run_check(out, v, k)
    f = fit()
    st = stage1(f, v)
    d = vdir(out, v)
    X, ext = stage1_draws(d, length(st.free))
    size(X, 1) >= k || error("only $(size(X, 1)) kept draws per chain, not $k")
    gi = phi_in_free(f, st)
    g = gate(X[1:k, :, :], gi)
    io = IOBuffer()
    p(xs...) = (println(io, xs...); println(xs...))
    p("# KIM-116: the stage-1 gate, variant ", v, ", on ", k, " kept draws per chain\n")
    p("Commit `", COMMIT, "`, job ", JOB, ". Extension ", ext ? "included" : "not run",
      ". Gate (prereg/kim-116 §4): R̂ ≤ ", RHAT_MAX, " and bulk ESS ≥ ", ESS_MIN,
      " on each of the 18 φ coordinates; the cut's ψ is reported, not gated.\n")
    p("| coordinate | R̂ | bulk ESS | gated | passes |")
    p("|---|---|---|---|---|")
    for (j, i) in enumerate(st.free)
        p(@sprintf("| %s | %.4f | %.1f | %s | %s |", f.names[i], g.rh[j], g.es[j],
                   j in gi ? "yes" : "no", j in gi ? (j in g.fail ? "no" : "yes") : "—"))
    end
    p("\n**The gate ", g.pass ? "passes" : "FAILS", "** on ", k, " draws per chain: ",
      length(gi) - length(g.fail), " of ", length(gi), " φ coordinates pass.",
      g.pass || ext ? "" : " The one extension (stage `extend`) runs next.")
    write(joinpath(d, "gate_$k.md"), String(take!(io)))
end

# φ's pooled kept draws, in f.phi's order.
function phi_draws(f, v)
    st = stage1(f, v)
    X, ext = stage1_draws(joinpath(OUTDIR[], v), length(st.free))
    return pooled(X)[:, phi_in_free(f, st)], ext
end
const OUTDIR = Ref("")

# ---------------------------------------------------------------- stage 2

function run_stage2(out, v, a, b)
    f = fit()
    d = vdir(out, v)
    Φ, ext = phi_draws(f, v)
    idx = round.(Int, range(1, size(Φ, 1); length = K2))
    for k in a:min(b, K2)
        name = "stage2_$k"
        has(d, name) && (say("$v stage 2 chain $k: done, skipped"); continue)
        x = copy(f.x_ref)
        x[f.phi] = Φ[idx[k], :]
        st = modular_stage(f.bulko, f.psi, x)
        z0 = stage_x0(st)
        isfinite(stage_logdensity(st, z0)) || error("stage 2 chain $k: θ_ref's ψ is not finite here")
        rng = Xoshiro(inner_seed(f, v, k))
        t = @elapsed ch = sample(rng, modular_model(st), NUTS(NADAPT2, TARGET_ACCEPT; max_depth = MAX_DEPTH),
                                 NKEEP2; initial_params = z0, progress = false)
        z = Array(ch[:, [Symbol("z[$i]") for i in 1:length(f.psi)], 1])
        stats = (; acceptance = vec(ch[:acceptance_rate]), errors = vec(ch[:numerical_error]),
                 depth = vec(ch[:tree_depth]), lp = vec(ch[:lp]))
        save(d, name, (; k, row = idx[k], phi = Φ[idx[k], :], z, stats, seconds = t, extended = ext);
             stage = "stage2")
        say(@sprintf("%s stage 2 chain %d of %d: %.0f s; mean acceptance %.2f, %d numerical errors, max depth %d",
                     v, k, K2, t, mean(stats.acceptance), count(>(0), stats.errors), maximum(stats.depth)))
    end
end

# The cut draws, K2 × NKEEP2 rows × 21 coordinates, mixed_parameter_names' order.
function cut_draws(f, v)
    d = joinpath(OUTDIR[], v)
    rows = map(1:K2) do k
        r = load(d, "stage2_$k")
        X = repeat(permutedims(f.x_ref), size(r.z, 1))
        X[:, f.phi] .= permutedims(r.phi)
        X[:, f.psi] .= r.z
        X
    end
    return reduce(vcat, rows)
end

function run_ppc(out, v)
    f = fit()
    s = f.both
    flat = cut_draws(f, v)
    idx = round.(Int, range(1, size(flat, 1); length = PPC_DRAWS))
    n = length(s.logz)
    Z = fill(NaN, n, PPC_DRAWS * PPC_REPS)
    failed = zeros(Bool, PPC_DRAWS)
    t = @elapsed Threads.@threads for j in eachindex(idx)
        x = flat[idx[j], :]
        r = _mixed_solve(s, exp.(x[1:length(s.names)]))
        r === nothing && (failed[j] = true; continue)
        m = mean_field_log_bulk(s, r.sol)
        C = cholesky(Symmetric(s.S + exp(x[end])^2 * I)).L
        rng = Xoshiro(SEED + 116_000_000 + 1000 * (v == "cut" ? 1 : 2) + j)
        for q in 1:PPC_REPS
            Z[:, (j - 1) * PPC_REPS + q] = m .+ s.delta .+ C * randn(rng, n)
        end
    end
    any(failed) && error("route (ii) failed at $(count(failed)) PPC draws")
    lo = [quantile(Z[i, :], PPC_LO) for i in 1:n]
    hi = [quantile(Z[i, :], PPC_HI) for i in 1:n]
    outside = (s.logz .< lo) .| (s.logz .> hi)
    save(vdir(out, v), "ppc", (; idx, lo, hi, outside, seconds = t); stage = "ppc")
    say(@sprintf("%s ppc: %d draws × %d in %.0f s; %d of %d rows outside", v, PPC_DRAWS, PPC_REPS,
                 t, count(outside), n))
end

# ---------------------------------------------------------------- merge

function sacct_cpuh(out)
    file = joinpath(out, "jobs.tsv")
    isfile(file) || return Dict{String, Float64}()
    cost = Dict{String, Float64}()
    for l in readlines(file)
        name, id = split(l, '\t')
        txt = try
            read(`sacct -j $id -X -n -P --format=CPUTimeRAW`, String)
        catch
            continue
        end
        h = sum((parse(Float64, strip(x)) for x in split(strip(txt), '\n') if !isempty(strip(x))); init = 0.0) / 3600
        cost[String(name)] = get(cost, String(name), 0.0) + h
    end
    return cost
end

function summarise(flat, psd)
    μ, sd = vec(mean(flat; dims = 1)), vec(std(flat; dims = 1))
    return (; μ, sd, contraction = 1 .- sd .^ 2 ./ psd .^ 2,
            lo = [quantile(flat[:, q], CI_LO) for q in axes(flat, 2)],
            hi = [quantile(flat[:, q], CI_HI) for q in axes(flat, 2)], R = cor(flat))
end

function run_merge(out)
    f = fit()
    s, names = f.both, f.names
    P = length(names)
    psd = mixed_prior_sd(s)
    tx = truth_x(f.ds, s)
    J = summarise(f.draws, psd)
    eno, fba = findfirst(==(:log_kcatF_R_ENO), names), findfirst(==(:log_kcatF_R_FBA), names)
    cost = sacct_cpuh(out)
    io = IOBuffer()
    p(xs...) = (println(io, xs...); println(xs...))
    p("# KIM-116: the modular arm beside the joint\n")
    p("Merge at `", COMMIT, "`, job ", JOB, ". SMOKE = ", SMOKE, ". Frozen at `prereg/kim-116`. ",
      "The joint is KIM-92's ", size(f.draws, 1), " draws. \"Close\" (prereg §5): every |shift| ≤ ",
      CLOSE_SHIFT, " joint SD and every |Δ contraction| ≤ ", CLOSE_CONTRACTION, ".\n")
    for v in VARIANTS
        d = joinpath(out, v)
        all(k -> has(d, "stage2_$k"), 1:K2) || (p("## ", v, "\n\nStage 2 is not complete; skipped.\n"); continue)
        C = cut_draws(f, v)
        open(joinpath(out, "cut_draws_$v.tsv"), "w") do o
            println(o, join(String.(names), '\t'))
            foreach(i -> println(o, join((@sprintf("%.17g", x) for x in C[i, :]), '\t')), axes(C, 1))
        end
        M = summarise(C, psd)
        shift = (M.μ .- J.μ) ./ J.sd
        dcon = M.contraction .- J.contraction
        far = [names[q] for q in 1:P if abs(shift[q]) > CLOSE_SHIFT]
        wide = [names[q] for q in 1:P if abs(dcon[q]) > CLOSE_CONTRACTION]
        inside = M.lo .<= tx .<= M.hi
        st2 = [load(d, "stage2_$k").stats for k in 1:K2]
        nerr = sum(count(>(0), q.errors) for q in st2)
        nonfinite = count(q -> !all(isfinite, q.lp), st2)
        ppc = has(d, "ppc") ? load(d, "ppc") : nothing
        p("## Variant: ", v, "\n")
        p(@sprintf("%d cut draws (%d stage-1 draws × %d). Close to the joint: **%s**.", size(C, 1), K2,
                   NKEEP2, isempty(far) && isempty(wide) ? "yes" : "no"),
          isempty(far) ? "" : " Shift beyond $(CLOSE_SHIFT): " * join(string.(far), ", ") * ".",
          isempty(wide) ? "" : " Contraction beyond $(CLOSE_CONTRACTION): " * join(string.(wide), ", ") * ".")
        p(@sprintf("\nStage 2: %d inner chains, %d numerical errors in kept draws, %d with a non-finite log density; mean acceptance %.2f.",
                   K2, nerr, nonfinite, mean(mean(q.acceptance) for q in st2)))
        ppc === nothing || p(@sprintf("PPC: %d of %d rows outside the central 99%% band (KIM-92's joint: 11).",
                                      count(ppc.outside), length(ppc.outside)))
        p("\n| coordinate | joint mean ± sd | cut mean ± sd | shift (joint SD) | contraction joint | contraction cut | truth in cut 90% | cut (mean − truth)/sd |")
        p("|---|---|---|---|---|---|---|---|")
        for q in 1:P
            p(@sprintf("| %s | %.4f ± %.4f | %.4f ± %.4f | %+.2f | %.3f | %.3f | %s | %+.2f |", names[q],
                       J.μ[q], J.sd[q], M.μ[q], M.sd[q], shift[q], J.contraction[q], M.contraction[q],
                       inside[q] ? "yes" : "no", (M.μ[q] - tx[q]) / M.sd[q]))
        end
        p(@sprintf("\nTruth inside the cut's 90%% interval: %d of %d.", count(inside), P))
        p("\nThe coupling the cut drops, correlations of ENO and FBA with φ (joint → cut):\n")
        p("| φ coordinate | r(ENO) joint | r(ENO) cut | r(FBA) joint | r(FBA) cut |")
        p("|---|---|---|---|---|")
        for q in f.phi
            p(@sprintf("| %s | %+.3f | %+.3f | %+.3f | %+.3f |", names[q], J.R[q, eno], M.R[q, eno],
                       J.R[q, fba], M.R[q, fba]))
        end
        p()
        open(joinpath(out, "modular_coordinates_$v.tsv"), "w") do o
            println(o, "coordinate\tjoint_mean\tjoint_sd\tcut_mean\tcut_sd\tshift\tjoint_contraction\tcut_contraction\ttruth\tinside90")
            for q in 1:P
                println(o, join((names[q], J.μ[q], J.sd[q], M.μ[q], M.sd[q], shift[q], J.contraction[q],
                                 M.contraction[q], tx[q], Int(inside[q])), '\t'))
            end
        end
    end
    p("## Cost (sacct, CPU-h)\n")
    if isempty(cost)
        p("No sacct records.")
    else
        p("| job group | CPU-h |")
        p("|---|---|")
        for (k, h) in sort!(collect(cost); by = first)
            p(@sprintf("| %s | %.1f |", k, h))
        end
        p(@sprintf("\nTotal %.1f CPU-h, beside the joint's 491.9 (Decision log 2026-11-17).", sum(values(cost))))
    end
    write(joinpath(out, "modular_summary.md"), String(take!(io)))
end

function main(args)
    length(args) >= 2 || error("usage: corea_modular_kim116.jl <stage> <out> [<args>]")
    stage, out = args[1], args[2]
    OUTDIR[] = out
    v = length(args) >= 3 ? args[3] : ""
    stage == "profile" && return run_profile(out)
    stage == "merge" && return run_merge(out)
    v in VARIANTS || error("stage $stage needs a variant, one of $VARIANTS")
    stage == "metric" && return run_metric(out, v)
    stage == "chain" && return run_chain(out, v, parse(Int, args[4]))
    stage == "check" && return run_check(out, v, parse(Int, args[4]))
    stage == "extend" && return run_extend(out, v, parse(Int, args[4]))
    stage == "stage2" && return run_stage2(out, v, parse(Int, args[4]), parse(Int, args[5]))
    stage == "ppc" && return run_ppc(out, v)
    error("unknown stage $stage")
end

main(ARGS)
