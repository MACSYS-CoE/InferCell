# Issue #89, after B3: what makes the mean field's GTP sit below the population
# mean, and does a correction pass KG (spec/phases/approximate-method.md §8)?
# Throwaway prototype, as B3's script was; the rules each stage is read against
# are in the sub-spec's §12 entry of 2026-10-22, written before the runs.
#
# Points. A point is a parameter set on M0:
# - B3's three M0 truths: m0_v6 (V6's replicate 16085), m0_low, m0_high;
# - F10's posterior on V6's dataset (`dev/scripts/f10_16a10_draws.tsv`): f10_ref
#   at its mean, and f10_d1..f10_d8, draws at Mahalanobis quantiles (k − ½)/8;
# - low_<p>_<±> and high_<p>_<±>: m0_low and m0_high moved ±1 posterior SD
#   (F10's SD on the log scale) in one target at a time.
# A kcatF that moves takes its kcatR with it, so the equilibrium constant is kept
# (draw_truth's Haldane rule).
#
# Arms. Each runs 200 (or 800) cells of the full two-way model at a point. All but
# `full` replace one source of cell-to-cell variation by its mean-field value:
# - full:       the published model; a per-cell clip census and gene expression kept;
# - drains_own: every deferred counter is drained at its expected accrual given the
#               cell's own state over the second (Euler over 1 s), not its sampled
#               events. Removes drain shot noise, keeps drain-rate heterogeneity;
# - drains_mf:  every counter is drained at the mean field's accrual. Removes both;
# - enzymes_mf: the catalytic edges and the volume read the mean field's protein
#               counts (rounded), not the cell's. The cell's own counts are carried
#               and restored around each step;
# - popmean:    the 60 s rebuild reads the population's mean pools, live, over 200
#               cells run in lockstep (the pools one 1 s step stale). Global RNG
#               seeded once, so its cells are not B3's seeds.
#
# Usage: julia --project dev/scripts/jensen_rethink_89.jl ensemble <point> <arm> <task> <ntasks> <ncells> <seed0> <outdir>
#        julia --project dev/scripts/jensen_rethink_89.jl meanfield <point> <outdir>
#        julia --project dev/scripts/jensen_rethink_89.jl popmean <point> <outdir>
#        julia --project dev/scripts/jensen_rethink_89.jl score <outdir>
#        julia --project dev/scripts/jensen_rethink_89.jl points
# Regenerate: dev/scripts/jensen_rethink_89.sh <outdir outside the worktrees>

using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates

const SMOKE = get(ENV, "SMOKE", "0") == "1"
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const B3_ARCHIVE = get(ENV, "B3_ARCHIVE", "/fred/oz022/tkimpson/infercell-89-b3-archive-2026-10-20")
const SIGMA_B = 0.10
const NKG = 200                     # KG's reference size
const HORIZON = SMOKE ? 120.0 : Float64(M0_HORIZON_S)
const ARMS = ("full", "drains_own", "drains_mf", "enzymes_mf")
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
const LOGNAMES = [:S_JCVISYN3A_0607, :S_JCVISYN3A_0779, :krnadeg, :kcatF_R_ENO, :kcatF_R_FBA]
const PARTNER = Dict(:kcatF_R_ENO => :kcatR_R_ENO, :kcatF_R_FBA => :kcatR_R_FBA)

slurm_job() = get(ENV, "SLURM_ARRAY_JOB_ID", get(ENV, "SLURM_JOB_ID", "none"))
truth_values(seed) = draw_truth(Xoshiro(seed), m0_models(); names = M0_TARGETS,
                                purpose = :calibration).values

# F10's draws on the log scale (all five chains, as exported), and their mean and SD.
function f10_draws()
    lines = readlines(F10_DRAWS)
    head = split(lines[1], '\t')
    cols = [findfirst(==(String(n)), head) for n in LOGNAMES]
    X = [parse(Float64, split(l, '\t')[c]) for l in lines[2:end], c in cols]
    return X, vec(mean(X; dims = 1)), vec(std(X; dims = 1))
end

# Move a truth's targets to new log values, carrying each kcatR with its kcatF.
function with_logs(base, logs::AbstractVector)
    d = Dict{Symbol, Float64}(base)
    for (n, l) in zip(LOGNAMES, logs)
        new = exp(l)
        haskey(PARTNER, n) && (d[PARTNER[n]] *= new / d[n])
        d[n] = new
    end
    return [k => d[k] for (k, _) in base]
end
logs_of(vals) = [log(Dict(vals)[n]) for n in LOGNAMES]

function points()
    P = Dict{String, Any}()
    v6, low, high = truth_values(16_085), truth_values(16_200_013), truth_values(16_200_002)
    P["m0_v6"], P["m0_low"], P["m0_high"] = v6, low, high
    X, μ, sd = f10_draws()
    P["f10_ref"] = with_logs(v6, μ)
    C = cov(X)
    m = [dot(X[i, :] .- μ, C \ (X[i, :] .- μ)) for i in axes(X, 1)]
    order = sortperm(m)
    for k in 1:8
        i = order[clamp(round(Int, (k - 0.5) / 8 * length(order)), 1, length(order))]
        P["f10_d$k"] = with_logs(v6, X[i, :])
    end
    for (tag, base) in (("low", low), ("high", high)), (j, n) in enumerate(LOGNAMES), s in (+1, -1)
        l = logs_of(base); l[j] += s * sd[j]
        P["$(tag)_$(n)_$(s > 0 ? "+" : "-")"] = with_logs(base, l)
    end
    return P
end
point_values(name) = points()[name]

const CME = (CoreATranscription, CoreATranscriptDecay, CoreATranslation)

# A gene-expression module with its reactions removed and everything else kept
# (B3's `Silent`).
struct Silent{M <: AbstractSubModel} <: AbstractSubModel
    m::M
end
for f in (:states, :parameters, :formalism, :inference_mode, :coupling, :inputs,
          :written_states, :rebuilt_params, :module_id, :membrane_protein_states,
          :contributed_states)
    @eval InferCell.$f(s::Silent) = InferCell.$f(s.m)
end
InferCell.rate_constants(p, t, s::Silent, pools) = InferCell.rate_constants(p, t, s.m, pools)
InferCell.reactions(::Silent) = [InferCell.Reaction((u, p, t, w) -> 0.0, (u, w) -> nothing)]
InferCell.reduction_notes(s::Silent) =
    vcat(InferCell.reduction_notes(s.m),
         [:formalism => "$(module_id(s.m))'s reactions are replaced by the mean field (#89 B3)."])

# Where the mean field's states and parameters sit in a driver built from `ms`.
function mf_maps(ms, mf)
    jump_models = [m for m in ms if formalism(m) === :jump]
    jnames = InferCell._block_names(ms, :jump)
    pnames = InferCell._block_param_names(jump_models)
    sidx = [findfirst(==(s), jnames) for s in states(mf)]
    pidx = [findfirst(==(p.name), pnames) for p in InferCell.model_free_params(parameters(mf))]
    (any(isnothing, sidx) || any(isnothing, pidx)) &&
        error("a mean-field state or parameter is missing from the driver")
    counters = Set(e.counter for m in ms for e in coupling(m) if e isa DeferredCounterEdge)
    is_counter = [s in counters for s in states(mf)]
    is_protein = [startswith(String(s), "P") for s in states(mf)]
    return (; sidx = Int.(sidx), pidx = Int.(pidx), is_counter, is_protein)
end

cme_of(ms) = [m for m in ms if any(T -> m isa T, CME)]

function rk4(x, p, t0, h, mf; n = 4)
    f(u, t) = collect(InferCell.dynamics(u, p, t, mf))
    s = h / n
    for i in 0:(n - 1)
        t = t0 + i * s
        k1 = f(x, t); k2 = f(x .+ s / 2 .* k1, t + s / 2)
        k3 = f(x .+ s / 2 .* k2, t + s / 2); k4 = f(x .+ s .* k3, t + s)
        x = x .+ s / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
    end
    return x
end

panel_rows(species) = [findfirst(==(o), species) for o in M0_PANEL]
mf_file(out, pt) = joinpath(out, "meanfield", "$pt.jls")

# B3's mean-field cell, which also keeps its mean state after every 1 s step
# (`xs[:, k + 1]` after step k), for the drains_mf and enzymes_mf arms.
function run_meanfield(pt, out)
    mkpath(joinpath(out, "meanfield"))
    full = m0_models()
    vals = point_values(pt)
    cme = cme_of(full)
    ms = AbstractSubModel[any(T -> m isa T, CME) ? Silent(m) : m for m in full]
    mf = MeanField(cme; edges = :none)
    d = build_problem(ms; tspan = (0.0, HORIZON), complete = true)
    set_parameters!(d, ms, vals)
    mp = mf_maps(ms, mf)
    x = Float64.(d.jump.u[mp.sidx])
    carry = zeros(length(x))
    rows = panel_rows(InferCell._block_names(ms, :ode))
    times, counts, xs = Float64[], Vector{Float64}[], Vector{Float64}[copy(x)]
    record!() = (push!(times, d.ode.t);
                 push!(counts, collect(d.ode.u[rows]) .* d.factor .+ d.rounding.remainders[rows]))
    record!()
    nsteps = round(Int, HORIZON / d.interval)
    for k in 1:nsteps
        handshake_step!(d)
        p = d.jump.p[mp.pidx]
        x1 = rk4(x, p, d.jump.t - d.interval, d.interval, mf)
        for i in eachindex(x)
            if mp.is_counter[i]
                carry[i] += x1[i] - x[i]
                w = floor(carry[i])
                d.jump.u[mp.sidx[i]] += Int(w)
                carry[i] -= w
            else
                d.jump.u[mp.sidx[i]] = round(Int, x1[i])
            end
        end
        x = x1
        push!(xs, copy(x))
        k % 60 == 0 && record!()
    end
    serialize(mf_file(out, pt),
              (; point = pt, values = vals, commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               species = M0_PANEL, times, counts = reduce(hcat, counts),
               mfstates = states(mf), xs = reduce(hcat, xs), clipped = clipping_census(d).clipped))
    println("$(now()) mean field $pt done")
end

# One cell of the full model at `vals` under `arm`. Returns its panel record at
# 60 s saves, the gene-expression states there, its volume, and its clips.
function run_cell(ms, vals, seed, arm, mfr, mf, mp)
    Random.seed!(seed)
    d = build_problem(ms; tspan = (0.0, HORIZON), complete = true)
    set_parameters!(d, ms, vals)
    Random.seed!(seed)
    ode_names = InferCell._block_names(ms, :ode)
    jnames = InferCell._block_names(ms, :jump)
    rows = panel_rows(ode_names)
    gerows = [mp.sidx[i] for i in eachindex(mp.sidx) if !mp.is_counter[i]]
    cidx = [mp.sidx[i] for i in eachindex(mp.sidx) if mp.is_counter[i]]
    cpos = findall(mp.is_counter)
    prot = [mp.sidx[i] for i in eachindex(mp.sidx) if mp.is_protein[i]]
    ppos = findall(mp.is_protein)
    rec = ClipRecord(d)
    times, counts, ge, vol = Float64[], Vector{Float64}[], Vector{Float64}[], Float64[]
    function record!()
        push!(times, d.ode.t)
        push!(counts, collect(d.ode.u[rows]) .* d.factor .+ d.rounding.remainders[rows])
        push!(ge, Float64.(d.jump.u[gerows])); push!(vol, d.volume_litres)
    end
    record!()
    carry = zeros(length(cidx))
    pending = zeros(length(cidx))          # drains_own: the expected accrual to drain next
    nsteps = round(Int, HORIZON / d.interval)
    for k in 1:nsteps
        if arm == "drains_mf"
            # Step k drains what accrued over step k − 1: the mean field's accrual then.
            inc = k == 1 ? zeros(length(cidx)) : mfr.xs[cpos, k] .- mfr.xs[cpos, k - 1]
            carry .+= inc
            w = floor.(carry); carry .-= w
            for (j, c) in enumerate(cidx); d.jump.u[c] = Int(w[j]); end
        elseif arm == "drains_own"
            carry .+= pending
            w = floor.(carry); carry .-= w
            for (j, c) in enumerate(cidx); d.jump.u[c] = Int(w[j]); end
        end
        xstart = Float64.(d.jump.u[mp.sidx])
        if arm == "enzymes_mf"
            own = copy(d.jump.u[prot])
            set = [round(Int, mfr.xs[ppos[j], k]) for j in eachindex(prot)]
            for (j, c) in enumerate(prot); d.jump.u[c] = set[j]; end
            handshake_step!(d)
            for (j, c) in enumerate(prot); d.jump.u[c] = own[j] + (d.jump.u[c] - set[j]); end
        else
            handshake_step!(d)
        end
        if arm == "drains_own"
            # The expected accrual over the jump step just taken, from the state it
            # started in (the exchange zeroes the counters) at the constants in force.
            xstart[cpos] .= 0
            p = d.jump.p[mp.pidx]
            a = collect(InferCell.dynamics(xstart, p, d.jump.t - d.interval, mf))
            pending = a[cpos] .* d.interval
        end
        record_clips!(rec, d)
        k % 60 == 0 && record!()
    end
    census = clipping_census(d)
    return (; times, counts = reduce(hcat, counts), ge = reduce(hcat, ge), vol,
            clipped = census.clipped, drains = census.drains, clips = clip_summary(rec),
            gestates = jnames[gerows])
end

ens_dir(out, pt, arm, seed0) = joinpath(out, "ens_$(pt)_$(arm)_s$(seed0)")

function run_ensemble(pt, arm, task, ntask, ncells, seed0, out)
    arm in ARMS || error("arm must be one of $ARMS")
    dir = ens_dir(out, pt, arm, seed0); mkpath(dir)
    ms = m0_models()
    vals = point_values(pt)
    mf = MeanField(cme_of(ms); edges = :none)
    mp = mf_maps(ms, mf)
    mfr = arm in ("drains_mf", "enzymes_mf") ? deserialize(mf_file(out, pt)) : nothing
    mfr === nothing || mfr.values == vals || error("$pt: the mean field was run at other values")
    mfr === nothing || mfr.mfstates == states(mf) || error("$pt: the mean field's states differ")
    cells = Any[]
    for k in 1:(SMOKE ? 4 : ncells)
        (k - 1) % ntask == task || continue
        seed = seed0 + k
        t0 = time()
        c = run_cell(ms, vals, seed, arm, mfr, mf, mp)
        push!(cells, (; seed, c...))
        println("$(now()) $pt $arm seed $seed ", @sprintf("%.1f s", time() - t0),
                " clipped $(c.clipped)/$(c.drains)"); flush(stdout)
    end
    serialize(joinpath(dir, "task_$(task).jls"),
              (; point = pt, arm, values = vals, commit = COMMIT, smoke = SMOKE, job = slurm_job(),
               species = M0_PANEL, cells))
end

# The live population-mean arm: every cell's rebuild reads the mean over the
# population of each pool it rebuilds from (in mM), taken before the step.
function run_popmean(pt, out)
    dir = joinpath(out, "ens_$(pt)_popmean"); mkpath(dir)
    ms = m0_models()
    vals = point_values(pt)
    n = SMOKE ? 4 : NKG
    ds = Any[]
    for k in 1:n
        Random.seed!(40_000 + k)
        d = build_problem(ms; tspan = (0.0, HORIZON), complete = true)
        set_parameters!(d, ms, vals)
        push!(ds, d)
    end
    ode_names = InferCell._block_names(ms, :ode)
    rows = panel_rows(ode_names)
    fed = unique(Symbol[s for r in ds[1].rebuilds for (s, k) in zip(r.species, r.pool_idxs) if k != 0])
    fidx = Dict(s => findfirst(==(s), ode_names) for s in fed)
    pm = Dict{Symbol, Float64}()
    f = (s, t) -> get(pm, s, nothing)
    for d in ds; override_rebuild_pools!(d, f); end
    println("$(now()) popmean $pt: $n cells, rebuild pools $(fed)"); flush(stdout)
    Random.seed!(89_000)
    cells = [(; times = Float64[], counts = Vector{Float64}[]) for _ in 1:n]
    rec!(i) = (push!(cells[i].times, ds[i].ode.t);
               push!(cells[i].counts, collect(ds[i].ode.u[rows]) .* ds[i].factor .+ ds[i].rounding.remainders[rows]))
    for i in 1:n; rec!(i); end
    nsteps = round(Int, HORIZON / ds[1].interval)
    for k in 1:nsteps
        for s in fed
            pm[s] = mean(d.ode.u[fidx[s]] for d in ds)
        end
        for (i, d) in enumerate(ds)
            handshake_step!(d)
            k % 60 == 0 && rec!(i)
        end
        k % 60 == 0 && (println("$(now()) popmean $pt step $k"); flush(stdout))
    end
    out_cells = [(; seed = 0, times = c.times, counts = reduce(hcat, c.counts)) for c in cells]
    serialize(joinpath(dir, "task_0.jls"),
              (; point = pt, arm = "popmean", values = vals, commit = COMMIT, smoke = SMOKE,
               job = slurm_job(), species = M0_PANEL, cells = out_cells, fed))
end

# --- scoring -----------------------------------------------------------------

function load_cells(dir)
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(dir)))
    isempty(files) && return nothing
    ps = [deserialize(joinpath(dir, f)) for f in files]
    cells = sort(vcat((p.cells for p in ps)...); by = c -> c.seed)
    prov = (; commit = join(unique(p.commit for p in ps), ","), job = join(unique(p.job for p in ps), ","),
            values = first(ps).values, smoke = all(p.smoke for p in ps))
    return cells, prov
end

logbulk(cs, j) = [log(max(mean(c.counts[i, j] for c in cs), 1.0)) for i in eachindex(M0_PANEL)]
logbulk_mat(cs) = reduce(hcat, [logbulk(cs, j) for j in 2:size(first(cs).counts, 2)])
kg_shift(M, L) = vec(sum((M .- L) .^ 2; dims = 2)) ./ (2SIGMA_B^2)
kg_noise(cs) = (h = length(cs) ÷ 2;
                vec(sum((logbulk_mat(cs[1:h]) .- logbulk_mat(cs[(h + 1):(2h)])) .^ 2; dims = 2)) ./ 2 ./ (2SIGMA_B^2))
mf_log(mfr) = log.(max.(mfr.counts[:, 2:end], 1.0))

function verdict(shift, noise)
    bad = [M0_PANEL[i] for i in eachindex(shift) if shift[i] > 3noise[i]]
    badf = [M0_PANEL[i] for i in eachindex(shift) if shift[i] > 3max(noise[i], 0.1)]
    return (; total = sum(shift), noise = sum(noise), bad, badf,
            pass = sum(shift) <= sum(noise) && isempty(bad),
            passf = sum(shift) <= sum(noise) && isempty(badf))
end

short(s) = replace(String(s), "M_" => "", "_c" => "")
const GTP = findfirst(==(:M_gtp_c), M0_PANEL)

function score(out)
    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    p("# Issue #89 B3 rethink: scores", SMOKE ? " — SMOKE RUN, NOT A RESULT" : "")
    p()
    p("Scored at $COMMIT (job $(slurm_job())) over `$out`. σ_b = $(Int(100SIGMA_B))%. ",
      "KG as agreed; `pass*` is the information-only variant with each pool's noise floored at 0.1 nats.")
    p()
    P = points()
    refs = Dict{Tuple{String, Int}, Any}()
    for pt in keys(P), seed0 in (40_000, 50_000)
        dir = ens_dir(out, pt, "full", seed0)
        isdir(dir) || continue
        r = load_cells(dir)
        r === nothing || (refs[(pt, seed0)] = r)
    end
    summary = Dict{String, Any}()

    # Check 1a: the full arm reproduces B3's cells at the same seeds.
    p("## Check: the full arm against B3's archive")
    p()
    for pt in ("m0_v6", "m0_low", "m0_high")
        (SMOKE || !haskey(refs, (pt, 40_000))) && continue
        cells = refs[(pt, 40_000)][1]
        b3dir = joinpath(B3_ARCHIVE, "ensemble_$pt")
        isdir(b3dir) || continue
        b3 = vcat((deserialize(joinpath(b3dir, f)).cells for f in readdir(b3dir) if endswith(f, ".jls"))...)
        b3s = vcat((deserialize(joinpath(b3dir, f)).seeds for f in readdir(b3dir) if endswith(f, ".jls"))...)
        diffs = Float64[]
        for c in cells
            j = findfirst(==(c.seed), b3s)
            j === nothing && continue
            push!(diffs, maximum(abs.(c.counts .- b3[j]) ./ max.(abs.(b3[j]), 1.0)))
        end
        p(@sprintf("- %s: %d cells matched by seed, worst relative difference %.3g.", pt, length(diffs),
                   isempty(diffs) ? NaN : maximum(diffs)))
    end
    p()

    # Stage 1: the arms at B3's truths, against the full arm's first 200 cells.
    p("## Stage 1: which variation makes the gap")
    p()
    p("GTP offset: the mean over saves of log(x̄_arm / x̄_full). `share` is the arm's offset ",
      "over the mean field's. An arm that removes the variation causing the gap sits where ",
      "the mean field does (share ≈ 1); one that removes an irrelevant variation stays with ",
      "the full model (share ≈ 0).")
    p()
    p("| Truth | Arm | cells | total | noise | KG | GTP offset | MF GTP offset | share of MF offset | clipped drains (cells with any) |")
    p("|---|---|---|---|---|---|---|---|---|---|")
    for pt in ("m0_v6", "m0_low", "m0_high")
        haskey(refs, (pt, 40_000)) || continue
        ref = refs[(pt, 40_000)][1][1:min(end, NKG)]
        L = logbulk_mat(ref)
        noise = kg_noise(ref)
        isfile(mf_file(out, pt)) || continue
        mfr = deserialize(mf_file(out, pt))
        Mmf = mf_log(mfr)
        mfoff = mean(Mmf[GTP, :] .- L[GTP, :])
        v = verdict(kg_shift(Mmf, L), noise)
        p(@sprintf("| %s | mean field | – | %.3g | %.3g | %s | %+.3f | %+.3f | 1 | – |", pt, v.total, v.noise,
                   v.pass ? "pass" : "FAIL", mfoff, mfoff))
        summary["$(pt)/mf"] = v
        clipcells = count(c -> c.clipped > 0, ref)
        clipped = sum(c.clipped for c in ref); drains = sum(c.drains for c in ref)
        p(@sprintf("| %s | full (reference) | %d | – | %.3g | – | 0 | %+.3f | 0 | %d of %d (%d) |", pt, length(ref),
                   v.noise, mfoff, clipped, drains, clipcells))
        for arm in ("drains_own", "drains_mf", "enzymes_mf", "popmean")
            dir = arm == "popmean" ? joinpath(out, "ens_$(pt)_popmean") : ens_dir(out, pt, arm, 40_000)
            isdir(dir) || continue
            r = load_cells(dir); r === nothing && continue
            cs = r[1][1:min(end, NKG)]
            M = logbulk_mat(cs)
            off = mean(M[GTP, :] .- L[GTP, :])
            va = verdict(kg_shift(M, L), noise)
            cl = haskey(first(cs), :clipped) ? @sprintf("%d of %d (%d)", sum(c.clipped for c in cs),
                                                         sum(c.drains for c in cs), count(c -> c.clipped > 0, cs)) : "–"
            p(@sprintf("| %s | %s | %d | %.3g | %.3g | %s | %+.3f | %+.3f | %.2f | %s |", pt, arm, length(cs),
                       va.total, va.noise, va.pass ? "pass" : "FAIL", off, mfoff, off / mfoff, cl))
            summary["$(pt)/$(arm)"] = (; va..., off, mfoff)
        end
    end
    p()
    p("Per-counter clips in the full arm (counter: cells that clipped, drains clipped, largest deficit):")
    p()
    for pt in ("m0_v6", "m0_low", "m0_high")
        haskey(refs, (pt, 40_000)) || continue
        ref = refs[(pt, 40_000)][1][1:min(end, NKG)]
        agg = Dict{Symbol, Tuple{Int, Int, Float64}}()
        for c in ref, r in c.clips
            a = get(agg, r.counter, (0, 0, 0.0))
            agg[r.counter] = (a[1] + 1, a[2] + r.drains, max(a[3], r.max_deficit))
        end
        p("- $pt: ", isempty(agg) ? "none" :
          join([@sprintf("%s %d, %d, %.0f", k, v...) for (k, v) in sort(collect(agg); by = x -> -x[2][2])], "; "))
    end
    p()

    # Stage 2A: an offset measured at a reference point, carried to targets.
    p("## Stage 2A: an offset at a posterior's width")
    p()
    p("δ = log x̄_ref − log x̄_MF(ref), per pool and save, from the reference point's full ",
      "ensemble (all its cells, seeds 40001…). Each target is scored against its own ",
      "200-cell full ensemble at seeds 50001…, independent of δ's.")
    p()
    # The last group is information only: δ from the dataset's truth (200 cells)
    # rather than the posterior mean, which is 1–3 posterior SDs away in places.
    groups = [("f10_ref", ["f10_d$k" for k in 1:8]),
              ("m0_low", sort([k for k in keys(P) if startswith(k, "low_")])),
              ("m0_high", sort([k for k in keys(P) if startswith(k, "high_")])),
              ("m0_v6", ["f10_d$k" for k in 1:8])]
    for (src, targets) in groups
        haskey(refs, (src, 40_000)) && isfile(mf_file(out, src)) || continue
        rc = refs[(src, 40_000)][1]
        δ = logbulk_mat(rc) .- mf_log(deserialize(mf_file(out, src)))
        p("### δ from $src ($(length(rc)) cells)")
        p()
        p("| Target | uncorrected total | corrected total | noise | KG uncorrected | KG corrected | pass* corrected | pools > 3× noise (corrected) | GTP shift corrected (noise) |")
        p("|---|---|---|---|---|---|---|---|---|")
        npass, npassf, ntot = 0, 0, 0
        for t in targets
            haskey(refs, (t, 50_000)) && isfile(mf_file(out, t)) || continue
            cs = refs[(t, 50_000)][1][1:min(end, NKG)]
            L = logbulk_mat(cs); noise = kg_noise(cs)
            M = mf_log(deserialize(mf_file(out, t)))
            vu = verdict(kg_shift(M, L), noise)
            sc = kg_shift(M .+ δ, L)
            vc = verdict(sc, noise)
            ntot += 1; npass += vc.pass; npassf += vc.passf
            p(@sprintf("| %s | %.3g | %.3g | %.3g | %s | %s | %s | %s | %.3g (%.3g) |", t, vu.total, vc.total, vc.noise,
                       vu.pass ? "pass" : "FAIL", vc.pass ? "pass" : "FAIL", vc.passf ? "pass" : "FAIL",
                       isempty(vc.bad) ? "–" : join(short.(vc.bad), ", "), sc[GTP], noise[GTP]))
            summary["2A/$src/$t"] = (; vu, vc)
        end
        p()
        p("**$src: corrected KG passes at $npass of $ntot targets ($npassf of $ntot with the floor).**")
        p()
    end
    write(joinpath(out, "scores.md"), take!(io))
    serialize(joinpath(out, "scores.jls"), summary)
    println("Wrote ", joinpath(out, "scores.md"))
end

function main()
    stage = ARGS[1]
    if stage == "ensemble"
        run_ensemble(ARGS[2], ARGS[3], parse(Int, ARGS[4]), parse(Int, ARGS[5]),
                     parse(Int, ARGS[6]), parse(Int, ARGS[7]), ARGS[8])
    elseif stage == "meanfield"
        run_meanfield(ARGS[2], ARGS[3])
    elseif stage == "popmean"
        run_popmean(ARGS[2], ARGS[3])
    elseif stage == "score"
        score(ARGS[2])
    elseif stage == "points"
        for (k, v) in sort(collect(points()); by = first)
            println(k, "  ", join([@sprintf("%s=%.4g", n, x) for (n, x) in v], " "))
        end
    else
        error("stage must be ensemble, meanfield, popmean, score or points")
    end
end

main()
