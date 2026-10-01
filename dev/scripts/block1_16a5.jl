# spec/phases/16-recovery.md task 16a.5, check V4: block 1 on M0 with three cells'
# paths fixed. Its samples of (ln kcatF_ENO, ln kcatF_FBA), at the truth's σ,
# must match a 2-D grid posterior computed by direct replay, at the 5, 25, 50, 75
# and 95% quantiles within 3 MC SE.
#
# The grid shares no code with block 1's target:
# - its prior is `logpdf(LogNormal, θ) + ln θ`, on the θ scale with the Jacobian;
# - its likelihood is phase 15's `observation_loglik`;
# - the replay and the reverse-constant rule are shared, and V1, V2 and 14c.1
#   check those.
#
# Stages, each a Slurm array or job, in order:
#   grid 1 <row>   coarse 21 × 21 over the prior ± 4.5 SD, one row per task
#   merge 1        the coarse marginals' moments set the fine box
#   grid 2 <row>   fine 41 × 41 over −6 to +6 posterior SD, one row per task
#   merge 2        marginal CDFs, exact quantiles, and 20 chain starts drawn from it
#   chain <task>   one chain of SWEEPS sweeps from its start, σ held
#   grid s, merge s, for s ≥ 3: only if merge s − 1 found mass at a box edge (above
#                  1e-8 of the peak). The ENO posterior is narrower than the coarse
#                  step, so the coarse moments misjudge it. If that axis's step is
#                  coarser than a quarter of its SD, it is widened 3× about its mean
#                  at 81 points. Otherwise it is extended at the same step to ±8 SD,
#                  and only the new rows are computed: merge s − 1 writes their
#                  indices to rows<s>.txt, and the rest are read from grid s − 1.
#   final          pool the chains, less BURN sweeps each, and test against the
#                  last merge's exact quantiles
#
# Usage: julia --project dev/scripts/block1_16a5.jl <stage> [<n>] [<task>]
# Regenerate: dev/scripts/block1_16a5.sh submits the stages with dependencies.

using InferCell
using Random
using Serialization
using Statistics
using Printf
using Distributions: LogNormal, logpdf, params

const STAGE = ARGS[1]
const DIR = joinpath(@__DIR__, "block1_16a5")
const REPLICATE = 16085
const NCELLS = 3
const NCHAINS = 20
const SWEEPS = parse(Int, get(ENV, "SWEEPS", "50"))
const BURN = 5          # chains start from a grid draw; drop their first sweeps
const FW = [:kcatF_R_ENO, :kcatF_R_FBA]
const ODE_NAMES = Set([:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA])
mkpath(DIR)

ms = m0_models()
N = round(Int, M0_HORIZON_S)

# The replicate, and each cell's path re-recorded from its seed. Recording is
# transparent (V2), so the recorded run must reproduce the dataset's latent.
function setup()
    ds = m0_dataset(REPLICATE; ncells = NCELLS)
    names = InferCell._block_names(ms, :ode)
    rows = [findfirst(==(s), names) for s in M0_PANEL]
    lrows = [findfirst(==(s), ds.species) for s in M0_PANEL]
    paths = map(enumerate(ds.seeds)) do (c, s)
        Random.seed!(s)
        d = build_m0()
        set_parameters!(d, ms, ds.truth.values)
        Random.seed!(s)
        record_path!(d)
        for k in 1:(N ÷ 60)
            for _ in 1:60
                handshake_step!(d)
            end
            x = [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in rows]
            x == ds.latent[c, lrows, k] || error("cell $s: the recorded run is not the dataset's")
        end
        recorded_path(d)
    end
    cme = [p for p in ds.truth.values if !(first(p) in ODE_NAMES)]
    base = build_m0()
    set_parameters!(base, ms, cme)
    obs = [ds.observed[c, 1:length(M0_PANEL), :] for c in 1:NCELLS]
    return (; ds, paths, cme, base, obs)
end
S = setup()
priors = [InferCell._parameter_by_name(ms, f).prior for f in FW]

# The grid's own log posterior at (u_ENO, u_FBA).
noise = NoiseModel(Modality(:metabolite, :lognormal, M0_PANEL; floor = 1.0))
function grid_logpost(u)
    lp = sum(logpdf(p, exp(x)) + x for (p, x) in zip(priors, u))
    writes = derived_ode_values(ms, [f => exp(x) for (f, x) in zip(FW, u)])
    names = InferCell._block_names(ms, :ode)
    rows = [findfirst(==(s), names) for s in M0_PANEL]
    for c in 1:NCELLS
        d = build_m0()
        set_parameters!(d, ms, vcat(S.cme, writes))
        r = PathReplay(S.paths[c], d; density = true)
        pred = Vector{Vector{Float64}}()
        for k in 1:(N ÷ 60)
            for _ in 1:60
                replay_step!(d, r)
            end
            push!(pred, [d.ode.u[i] * d.factor + d.rounding.remainders[i] for i in rows])
        end
        data = ObservedData(S.ds.times, S.obs[c], M0_PANEL)
        lp += path_logdensity(r) + observation_loglik(noise, data, i -> pred[i], [S.ds.sigma])
    end
    return lp
end

prior_box() = [(params(p)[1] - 4.5params(p)[2], params(p)[1] + 4.5params(p)[2]) for p in priors]
axes_for(stage) = stage == 1 ?
    [range(lo, hi; length = 21) for (lo, hi) in prior_box()] :
    deserialize(joinpath(DIR, "box$(stage).jls"))
# A row of grid `stage`: its own file, or the previous stage's row at the same
# ln ENO when the axis was extended rather than rebuilt.
function grid_row(stage, r)
    f = joinpath(DIR, "grid$(stage)_row$(r).jls")
    isfile(f) && return deserialize(f).vals
    off = deserialize(joinpath(DIR, "offset$(stage).jls"))
    return grid_row(stage - 1, r - off)
end

if STAGE == "grid"
    stage, row = parse(Int, ARGS[2]), parse(Int, ARGS[3])
    ax = axes_for(stage)
    vals = [grid_logpost([ax[1][row + 1], b]) for b in ax[2]]
    serialize(joinpath(DIR, "grid$(stage)_row$(row).jls"), (row = row, vals = vals))
    println("grid $stage row $row done")
elseif STAGE == "merge"
    stage = parse(Int, ARGS[2])
    ax = axes_for(stage)
    lp = reduce(vcat, [grid_row(stage, r)' for r in 0:(length(ax[1]) - 1)])
    w = exp.(lp .- maximum(lp))
    mE, mF = vec(sum(w; dims = 2)), vec(sum(w; dims = 1))
    mom(u, m) = (c = sum(u .* m) / sum(m); (c, sqrt(sum((u .- c) .^ 2 .* m) / sum(m))))
    (cE, sE), (cF, sF) = mom(ax[1], mE), mom(ax[2], mF)
    if stage == 1
        box = [range(cE - 6sE, cE + 6sE; length = 41), range(cF - 6sF, cF + 6sF; length = 41)]
        serialize(joinpath(DIR, "box2.jls"), box)
        @printf("coarse: ln ENO %.4f ± %.4f, ln FBA %.4f ± %.4f\n", cE, sE, cF, sF)
    else
        edge = maximum([mE[1], mE[end]]) / maximum(mE), maximum([mF[1], mF[end]]) / maximum(mF)
        if maximum(edge) > 1e-8
            edge[2] > 1e-8 && error("FBA's axis truncates too; extend this script to widen it")
            a = ax[1]
            if step(a) > sE / 4
                box = [range(cE - 3 * (last(a) - first(a)) / 2, cE + 3 * (last(a) - first(a)) / 2;
                             length = 81), ax[2]]
                rows = collect(0:80)
            else
                klo = max(0, ceil(Int, (first(a) - (cE - 8sE)) / step(a)))
                khi = max(0, ceil(Int, ((cE + 8sE) - last(a)) / step(a)))
                n = length(a) + klo + khi
                box = [range(first(a) - klo * step(a); step = step(a), length = n), ax[2]]
                serialize(joinpath(DIR, "offset$(stage + 1).jls"), klo)
                rows = vcat(collect(0:(klo - 1)), collect((klo + length(a)):(n - 1)))
            end
            open(io -> println(io, join(rows, ",")), joinpath(DIR, "rows$(stage + 1).txt"), "w")
            serialize(joinpath(DIR, "box$(stage + 1).jls"), box)
            @printf("stage %d: mass at the box edge (%.1e, %.1e); grid %d widens it to ln ENO %.4f to %.4f, ln FBA %.4f to %.4f\n",
                    stage, edge..., stage + 1, first(box[1]), last(box[1]), first(box[2]), last(box[2]))
            exit(0)
        end
        tcdf(u, f) = (c = [0.0; cumsum([(f[i] + f[i + 1]) / 2 * step(u) for i in 1:length(u) - 1])];
                      c ./ c[end])
        quant(cdf, u, p) = (i = findfirst(>=(p), cdf);
                            u[i - 1] + (p - cdf[i - 1]) / (cdf[i] - cdf[i - 1]) * step(u))
        ps = [0.05, 0.25, 0.5, 0.75, 0.95]
        q = [[quant(tcdf(ax[1], mE), ax[1], p) for p in ps],
             [quant(tcdf(ax[2], mF), ax[2], p) for p in ps]]
        # Chain starts: grid cells drawn by mass, jittered within the cell.
        rng = Xoshiro(16086)
        cells = vec(CartesianIndices(w))
        mass = vec(w) ./ sum(w)
        starts = map(1:NCHAINS) do _
            i = cells[findfirst(>=(rand(rng)), cumsum(mass))]
            [ax[1][i[1]] + (rand(rng) - 0.5) * step(ax[1]), ax[2][i[2]] + (rand(rng) - 0.5) * step(ax[2])]
        end
        isfile(joinpath(DIR, "exact.jls")) ||
            serialize(joinpath(DIR, "exact.jls"), (quantiles = q, ps = ps, starts = starts,
                                                   moments = ((cE, sE), (cF, sF)), edge = edge))
        serialize(joinpath(DIR, "exact_final.jls"), (quantiles = q, ps = ps, stage = stage,
                                                     moments = ((cE, sE), (cF, sF)), edge = edge))
        @printf("fine: ln ENO %.4f ± %.4f, ln FBA %.4f ± %.4f; edge mass %.1e, %.1e\n",
                cE, sE, cF, sF, edge...)
    end
elseif STAGE == "chain"
    task = parse(Int, ARGS[2])
    ex = deserialize(joinpath(DIR, "exact.jls"))
    b = Block1(ms, [Block1Cell(snapshot(S.base), S.paths[c], S.obs[c]) for c in 1:NCELLS];
               forwards = FW, panel = M0_PANEL)
    st = Block1State(b, ex.starts[task + 1], S.ds.sigma)
    rng = Xoshiro(160_860 + task)
    draws = zeros(SWEEPS, 2)
    t = @elapsed for s in 1:SWEEPS
        block1_update!(rng, b, st; update_sigma = false)
        draws[s, :] = st.u
    end
    serialize(joinpath(DIR, "chain_$(task).jls"), (task = task, draws = draws, seconds = t))
    @printf("chain %d: %d sweeps in %.0f s\n", task, SWEEPS, t)
elseif STAGE == "final"
    ex = deserialize(joinpath(DIR, "exact_final.jls"))
    @printf("exact quantiles from grid stage %d; edge mass %.1e, %.1e\n", ex.stage, ex.edge...)
    chains = [deserialize(joinpath(DIR, "chain_$(k).jls")).draws[(BURN + 1):end, :]
              for k in 0:(NCHAINS - 1)]
    @printf("V4 on M0: replicate %d, %d cells, σ = %.4f held; %d chains × %d sweeps\n",
            REPLICATE, NCELLS, S.ds.sigma, NCHAINS, size(chains[1], 1))
    ok = true
    for (j, name) in enumerate(("ln kcatF_ENO", "ln kcatF_FBA"))
        for (p, q) in zip(ex.ps, ex.quantiles[j])
            per_chain = [mean(c[:, j] .<= q) for c in chains]
            m, se = mean(per_chain), std(per_chain) / sqrt(length(per_chain))
            pass = abs(m - p) <= 3se
            global ok &= pass
            @printf("  %s at p = %.2f: exact quantile %.4f, chain fraction below %.4f ± %.4f (%s)\n",
                    name, p, q, m, se, pass ? "within 3 SE" : "OUTSIDE 3 SE")
        end
    end
    println(ok ? "V4 passes" : "V4 FAILS")
end
