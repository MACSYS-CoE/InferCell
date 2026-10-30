# KIM-35: does lactate export move with both the inferred θ_I and the held θ_H?
# A pilot on the reduced model (`d11_models`), on θ held out of the evidence run:
# θ_I from B7's `v6` fit (m0_v6, not one of the 50 evidence datasets), θ_H from
# its prior, every rng seed from M4's reserved block 18_435_xxx. Nothing here
# reads an evidence posterior.
#
# Usage (dev/scripts/sensitivity_pilot_kim35.sh submits every stage):
#   julia --project dev/scripts/sensitivity_pilot_kim35.jl <stage> <outdir>
# Stages, each at H = 6300 (one published cycle) and H = 600 (M0's horizon):
#   arm_I_<H>   θ_I: 50 v6 posterior draws; θ_H at nominal (`completion = :none`).
#               The rows are independent draws, so var − MC is Var_I.
#   arm_H_<H>   θ_I fixed at v6's posterior mean of log θ (one draw); θ_H: 50 prior
#               draws at width 1 (`completion = :prior`, M = 1). The columns are
#               independent draws, so var − MC is Var_H.
#   local_<H>   the mean field at that reference point, each of the 20 parameters
#               moved by ±0.2 in ln θ: d ln Q / d ln θ, for both Q
#   merge       the rule below, for both Q candidates and both horizons
# Each arm runs 4 cells per pair, threaded. Both Q candidates (`lactate_out`,
# `exported`) are read from the same cells, through `carbon_accounts`.
#
# The rule, declared before the run (M4 plan, approved 2026-10-30): a (Q, H)
# passes if each group's parameter variance, net of the MC term, is
#   (i) above 2 bootstrap SEs (2,000 resamples of the arm's pairs), and
#   (ii) at least 5% of Var_I + Var_H.
# If no candidate passes at the horizon the freeze picks, the note recommends a
# change of Q and M4 stops for the user. The unfiltered cell statuses are
# recorded too, so V's thresholds see held-out acceptance rates.

using InferCell
using InferCell: AbstractSubModel
using Random
using Serialization
using Statistics: mean, std, var
using Printf
using Dates

const STAGE = ARGS[1]
const OUT = ARGS[2]
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const JOB = get(ENV, "SLURM_JOB_ID", "none")
const V6 = ["/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/chain_$c.jls" for c in 1:2]
const QS = (:lactate_out, :exported)
const SMOKE = get(ENV, "SMOKE", "0") == "1"      # a plumbing pass at tiny sizes, not a result
const N_DRAWS, N_CELLS, B = SMOKE ? (3, 2, 50) : (50, 4, 2_000)
# A smoke pass runs each stage's horizon H as H / 10 (630 s and 60 s).
horizon(H) = SMOKE ? H / 10 : H

say(x...) = (println(now(), " ", x...); flush(stdout))
save(name, x) = serialize(joinpath(OUT, "$name.jls"),
                          (; result = x, commit = COMMIT, job = JOB, stage = STAGE,
                           threads = Threads.nthreads(), when = now()))
load(name) = deserialize(joinpath(OUT, "$name.jls")).result

function reference(pd)
    ref = InferCell._theta_ref(pd)
    one = PosteriorDraws(pd.names, reshape(last.(ref), 1, :), nothing, [1], [1], nothing,
                         nothing, (kind = :theta_ref, of = pd.source))
    return ref, one
end

# μ̂, σ̂² and the cell count of a quantity read from each pair's finished cells.
function moments(p::Prediction, q::Symbol)
    xs = [[getproperty(c.accounts, q) for c in cs if c.failure === :none] for cs in vec(p.cells)]
    return [isempty(x) ? NaN : mean(x) for x in xs], [length(x) >= 2 ? var(x) : NaN for x in xs],
           length.(xs)
end

# The spread over the pairs, net of the MC term; the pairs are independent draws.
net(mu, s2, n) = var(mu) - mean(s2 ./ n)

function numbers(p, q; rng)
    mu, s2, n = moments(p, q)
    ok = isfinite.(mu) .& isfinite.(s2)
    mu, s2, n = mu[ok], s2[ok], n[ok]
    boot = map(1:B) do _
        i = rand(rng, 1:length(mu), length(mu))
        net(mu[i], s2[i], n[i])
    end
    return (; q, mean = mean(mu), var = var(mu), mc = mean(s2 ./ n), var_theta = net(mu, s2, n),
            se = std(boot), sigma2 = mean(s2), pairs = length(mu), dropped = count(!, ok))
end

# Acceptance under candidate rules, read on the arm's own cells (not V's: V is
# not fixed here). Clipping allowed, 14c's starvation threshold, no floors.
function acceptance(p)
    map((0.5, 1.0)) do f
        V = ViabilityRule(; label = "KIM-35 candidate", allow_clipped = true,
                          min_viable_fraction = f)
        (; min_viable_fraction = f,
         acceptance = mean(assess_viability(cs, V).viable for cs in p.cells))
    end
end

function arm(which, H)
    pd = posterior_draws(V6)
    dms = d11_models()
    ref, one = reference(pd)
    seed = 18_435_000 + (which === :I ? 100 : 200) + (H == 6300 ? 1 : 2)
    p = which === :I ?
        predict(dms, pd; horizon = horizon(H), n_parameter_draws = N_DRAWS, n_cells = N_CELLS,
                rng = Xoshiro(seed)) :
        predict(dms, one; horizon = horizon(H), completion = :prior, n_held = N_DRAWS,
                n_cells = N_CELLS, rng = Xoshiro(seed))
    say("arm $which at $H s: $(size(p.mu)) pairs, statuses $(p.statuses)")
    rng = Xoshiro(seed + 10_000)
    nums = Dict(q => numbers(p, q; rng) for q in QS)
    for q in QS
        say("  ", nums[q])
    end
    acc = acceptance(p)
    say("  candidate acceptance: ", acc)
    save("arm_$(which)_$(Int(H))", (; numbers = nums, statuses = p.statuses, acceptance = acc,
                                    held = p.held, reference = ref, prediction = p))
end

function local_sens(H)
    pd = posterior_draws(V6)
    ms = metered(d11_models())
    ref, _ = reference(pd)
    held = held_parameters(ms, pd).held
    nominal = nominal_parameter_values(ms)
    base = vcat(ref, [h => nominal[h] for h in held])
    f = InferCell._forward_setup(ms, :mean_field; horizon = horizon(H))
    solve(vals) = InferCell._restoring_global_rng(() ->
        InferCell._mean_field_cell(f, haldane_values(ms, vals), :lactate_out))
    q(c) = c.failure === :none ? (; lactate_out = c.accounts.lactate_out,
                                  exported = c.accounts.exported) : nothing
    c0 = q(solve(base))
    say("reference at $H s: ", c0)
    # The log-scale SD each parameter is drawn with: the v6 posterior's for θ_I,
    # the promoters' prior σ = log 2 for θ_H (transcription.jl).
    sd = Dict(n => std(log.(pd.values[:, j])) for (j, n) in enumerate(pd.names))
    all(h -> startswith(String(h), "S_"), held) || error("a held parameter is not a promoter")
    rows = map(SMOKE ? vcat(pd.names[1:2], held[1:2]) : vcat(pd.names, held)) do name
        moved(s) = [n => (n === name ? v * exp(s) : v) for (n, v) in base]
        up, dn = q(solve(moved(0.2))), q(solve(moved(-0.2)))
        e = up === nothing || dn === nothing ? (; lactate_out = NaN, exported = NaN) :
            (; lactate_out = (log(up.lactate_out) - log(dn.lactate_out)) / 0.4,
             exported = (log(up.exported) - log(dn.exported)) / 0.4)
        s = get(sd, name, log(2.0))
        r = (; name, group = name in pd.names ? :inferred : :held, log_sd = s,
             elasticity = e)
        say("  ", r)
        r
    end
    save("local_$(Int(H))", (; reference = c0, rows))
end

function merge()
    io = open(joinpath(OUT, "sensitivity.md"), "w")
    out(x...) = (println(io, x...); println(x...))
    out("# KIM-35 sensitivity at $COMMIT, job $JOB, $(now())\n")
    for H in (6300, 600)
        a, h = load("arm_I_$H"), load("arm_H_$H")
        out("## H = $H s\n")
        out("| Q | Var_I ± SE | Var_H ± SE | MC (I, H) | share I | share H | pass |")
        out("|---|---|---|---|---|---|---|")
        for q in QS
            i, k = a.numbers[q], h.numbers[q]
            tot = max(i.var_theta, 0) + max(k.var_theta, 0)
            ok_i = i.var_theta > 2i.se && i.var_theta >= 0.05tot
            ok_h = k.var_theta > 2k.se && k.var_theta >= 0.05tot
            out(@sprintf("| %s | %.4g ± %.2g | %.4g ± %.2g | %.3g, %.3g | %.3f | %.3f | %s |",
                         q, i.var_theta, i.se, k.var_theta, k.se, i.mc, k.mc,
                         max(i.var_theta, 0) / tot, max(k.var_theta, 0) / tot,
                         ok_i && ok_h ? "yes" : "no (I $(ok_i), H $(ok_h))"))
        end
        out("\nStatuses, arm I: $(a.statuses); arm H: $(h.statuses).")
        out("Candidate acceptance (clipping allowed, starved at 0.5, no floors), arm I: ",
            a.acceptance, "; arm H: ", h.acceptance, ".\n")
        l = load("local_$H")
        out("Mean-field d ln Q / d ln θ at the reference (Q = $(l.reference)):\n")
        out("| parameter | group | log SD | lactate_out | exported | (ε·SD)², lactate_out |")
        out("|---|---|---|---|---|---|")
        for r in l.rows
            out(@sprintf("| %s | %s | %.3f | %.4f | %.4f | %.3g |", r.name, r.group, r.log_sd,
                         r.elasticity.lactate_out, r.elasticity.exported,
                         (r.elasticity.lactate_out * r.log_sd)^2))
        end
        out("")
    end
    close(io)
end

function main()
    mkpath(OUT)
    say("KIM-35 sensitivity pilot, stage $STAGE, at $COMMIT, job $JOB, ",
        "$(Threads.nthreads()) thread(s)")
    m = match(r"^(arm_I|arm_H|local)_(\d+)$", STAGE)
    if STAGE == "merge"
        merge()
    elseif m !== nothing
        H = parse(Float64, m[2])
        m[1] == "local" ? local_sens(H) : arm(m[1] == "arm_I" ? :I : :H, H)
    else
        error("unknown stage $STAGE")
    end
    say("done")
end

main()
