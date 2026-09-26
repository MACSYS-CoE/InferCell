# Spec §11 phase 14b: the derivative, ensemble and census checks on the
# assembled Core A′, at full scale. Check 7's census has its own driver,
# `corea_census.jl`.
#
# One section per invocation, chosen by SECTION, so the slow ones run side by
# side as Slurm jobs:
#
#   SECTION=floor   14.2  check 1b's particle-floor report over the pinned cycle,
#                         the split of spec §3 (amended 2026-09-26), and the
#                         biased-ENO mutant's trajectory for the band mutation.
#   SECTION=band    14.2  check 1b's Langevin ensemble. An array: task TASK of
#                         NTASK runs its share of the trajectories and writes
#                         them raw; SECTION=merge-band reads them all.
#   SECTION=check8  14.8  check 8 on the frozen variant, with its mutation.
#   SECTION=check6  14.9  check 6 and task 9.9's rerun on the assembled model.
#   SECTION=f5      14b.5 the two external comparisons, over F5_SEEDS cycles.
#
# Each section writes `validation_14b/<section>.md`. SECTION=report
# concatenates them into corea_validation_14b_result.md, which also carries
# F3's continuation table (spec §11 task 14.10).
#
# Usage: SECTION=<s> julia --project dev/scripts/corea_validation_14b.jl [task ntask]
# Regenerate: see corea_validation_14b.slurm.

using InferCell
using Printf
using Dates
using Random
using Statistics
using LinearAlgebra
using Serialization

const ROOT = joinpath(@__DIR__, "..", "..")
include(joinpath(ROOT, "test", "nucleotide_test_models.jl"))
include(joinpath(ROOT, "test", "trna_test_models.jl"))
include(joinpath(ROOT, "test", "corea_validation_doubles.jl"))
include(joinpath(ROOT, "test", "corea_control_doubles.jl"))
include(joinpath(ROOT, "test", "corea_langevin_doubles.jl"))

const SECTION = get(ENV, "SECTION", "")
const TASK = parse(Int, get(ARGS, 1, get(ENV, "SLURM_ARRAY_TASK_ID", "0")))
const NTASK = parse(Int, get(ARGS, 2, get(ENV, "NTASK", "1")))
const SMOKE = get(ENV, "SMOKE", "0") == "1"
const CYCLE = SMOKE ? 120 : round(Int, COREA_CYCLE_S)
const SEED = 1410
const H = 0.01                                   # §12 2026-09-26 E
const H_COARSE = 0.02                            # the halving comparison
const N_TRAJ = SMOKE ? 4 : 100
const N_COARSE = SMOKE ? 2 : 20
const F5_SEEDS = SMOKE ? [SEED] : SEED .+ (0:19)
const TRNA_SEEDS = SMOKE ? [SEED] : SEED .+ (0:4)
const OUTDIR = joinpath(@__DIR__, "validation_14b")
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test dev/scripts/corea_validation_14b.jl`, String)))
const JOB = get(ENV, "SLURM_JOB_ID", "none")

mkpath(OUTDIR)
out = IOBuffer()
p(args...) = (println(out, args...); println(args...); flush(stdout))
fmt(x) = @sprintf("%.3g", x)
provenance() = p("Section `$SECTION` at commit $COMMIT", DIRTY ? " (tree dirty)" : "",
                 ", Slurm job $JOB, on $(gethostname()), $(now()), Julia $VERSION.",
                 SMOKE ? " **SMOKE RUN, NOT A RESULT.**" : "")
save(name) = (write(joinpath(OUTDIR, "$name.md"), take!(out)); println("wrote $name.md"))

ode_pnames(ms) = unique(Symbol[q.name for m in ms if formalism(m) === :ode
                                for q in model_free_params(parameters(m))])
function seeded_build(ms; seed = SEED, n = CYCLE)
    Random.seed!(seed)
    d = build_problem(ms; tspan = (0.0, Float64(n)), complete = true)
    Random.seed!(seed)
    return d
end
const EXTERNAL = [:M_lac__L_e]
# The mutant of the band check: enolase at half its catalytic constants, which
# feeds 2PG and PEP directly.
eno_half_models() = (ms = corea_models();
                     ms[1] = CentralGlycolysis(enzymes = :translated,
                                               free = [:kcatF_R_ENO, :kcatR_R_ENO]); ms)
function eno_half_driver(ms)
    d = seeded_build(ms)
    pn = ode_pnames(ms)
    for s in (:kcatF_R_ENO, :kcatR_R_ENO)
        d.ode.p[findfirst(==(s), pn)] *= 0.5
    end
    return d
end

# Every section runs inside one function, so a loop's assignments update the
# section's variables rather than making new locals (Julia's soft scope).
function main()
# ===========================================================================
if SECTION == "floor"
    p("## Check 1b: the particle floor over the pinned cycle")
    p()
    provenance()
    p()
    ms = corea_models()
    d = seeded_build(ms)
    wall = @elapsed vr = validation_run!(ms, d, CYCLE; moieties = corea_moieties(carbon = false))
    rep = particle_floor(vr)
    split = langevin_pools(rep; external = EXTERNAL)
    p(@sprintf("Seed %d, one cycle of %d handshakes (%.0f s). Every ODE state, in particles at the volume it was held at, smallest minimum first. Flagged: minimum below %d.",
               SEED, CYCLE, wall, PARTICLE_FLOOR))
    p()
    p("| state | min | at (s) | median | time under 1 particle | flagged | check 1b |")
    p("|---|---|---|---|---|---|---|")
    for r in rep
        role = r.species in EXTERNAL ? "external, not assessed" :
               r.species in split.excluded ? "**excluded**: median under $(CONTINUUM_MEDIAN)" :
               r.species in split.cross_check ? "Langevin cross-check" :
               r.flagged ? "flagged" : ""
        p(@sprintf("| %s | %.3g | %.0f | %.4g | %.1f%% | %s | %s |", r.species, r.min_particles,
                   r.t, r.median_particles, 100r.below_one, r.flagged ? "yes" : "no", role))
    end
    p()
    p("$(length(flagged_states(rep))) of $(length(rep)) states fall below $(PARTICLE_FLOOR) particles. ",
      "Excluded outright, with a median under $(CONTINUUM_MEDIAN) particles: ",
      join(("`$s`" for s in split.excluded), ", "), ". ",
      "Cross-checked against the Langevin ensemble: ", join(("`$s`" for s in split.cross_check), ", "), ".")
    p()
    serialize(joinpath(OUTDIR, "floor.jls"), (cross_check = split.cross_check, excluded = split.excluded))

    # The band mutation's reference: the same seed with enolase halved.
    me = eno_half_models()
    de = eno_half_driver(me)
    ref = Vector{Vector{Float64}}(undef, CYCLE)
    for k in 1:CYCLE
        handshake_step!(de)
        ref[k] = collect(Float64, de.ode.u)
    end
    serialize(joinpath(OUTDIR, "eno_half.jls"), ref)
    p("The band mutation's trajectory, enolase's two catalytic constants halved at seed $SEED, is recorded for `merge-band`.")
    save("floor")

# ===========================================================================
elseif SECTION == "band"
    pools = deserialize(joinpath(OUTDIR, "floor.jls")).cross_check
    ms = corea_models()
    lb = LangevinBlock(ms, ode_pnames(ms))
    idx = [findfirst(==(s), lb.names) for s in pools]
    jobs = vcat([(H, r) for r in 1:N_TRAJ], [(H_COARSE, 10_000 + r) for r in 1:N_COARSE])
    mine = [j for (i, j) in enumerate(jobs) if (i - 1) % NTASK == TASK]
    paths = Dict{Float64, FrozenPath}()
    for h in unique(first.(mine))
        t = @elapsed paths[h] = record_frozen_path(lb, ms, seeded_build(ms), CYCLE; h)
        println("path at h = $h recorded in $(round(t; digits = 1)) s"); flush(stdout)
    end
    res = []
    for (h, r) in mine
        t = @elapsed (tr, c) = replay(lb, paths[h], MersenneTwister(r))
        push!(res, (h = h, r = r, clamped = c, x = [tr[k][i] for k in eachindex(tr), i in idx]))
        println("trajectory $r at h = $h: $(round(t; digits = 1)) s, $c steps clamped"); flush(stdout)
    end
    # The reference and its path's integrator check, once.
    extra = TASK == 0 ? Dict(h => (u_ref = [fp.u_ref[k][idx] for k in eachindex(fp.t)],
                                   t = fp.t, Ω = fp.Ω,
                                   a_rel = Dict(s => maximum(abs(fp.a[k][i]) / max(abs(fp.u_ref[k][i]), 1e-300)
                                                             for k in eachindex(fp.a))
                                                for (s, i) in ((s, findfirst(==(s), lb.names))
                                                               for s in (:M_g6p_c, :M_fdp_c, :M_g3p_c, :M_13dpg_c,
                                                                         :M_3pg_c, :M_2pg_c, :M_pep_c, :M_pyr_c,
                                                                         :M_nadh_c, :M_lac__L_c))))
                             for (h, fp) in paths) : nothing
    serialize(joinpath(OUTDIR, "band_task_$TASK.jls"), (pools = pools, res = res, extra = extra))
    println("done")

# ===========================================================================
elseif SECTION == "merge-band"
    fl = deserialize(joinpath(OUTDIR, "floor.jls"))
    pools = fl.cross_check
    parts = [deserialize(joinpath(OUTDIR, f)) for f in readdir(OUTDIR) if occursin(r"^band_task_\d+\.jls$", f)]
    res = reduce(vcat, [pt.res for pt in parts])
    extra = only(pt.extra for pt in parts if pt.extra !== nothing)
    eno = deserialize(joinpath(OUTDIR, "eno_half.jls"))
    p("## Check 1b: the chemical-Langevin cross-check")
    p()
    provenance()
    p()
    p("The ensemble is the test-local double in `test/corea_langevin_doubles.jl`: the ODE block's ",
      "reactions as a chemical Langevin equation, each reaction's forward and reverse channel ",
      "separately, integrated by local linearization, with the jump path frozen at seed $SEED's ",
      "(spec §12, 2026-09-26 E). The band is the 5th to 95th percentile at every handshake.")
    p()
    hs = sort(unique(r.h for r in res))
    ex = extra[H]
    p("**The integrator against the pinned Rodas5P over one interval**, at h = $H: the largest ",
      "`|aₖ|` relative to the state on species no counter touches, which is this integrator's own error.")
    p()
    p("| species | " * join(("h = $h" for h in sort(collect(keys(extra)))), " | ") * " |")
    p("|---|" * repeat("---|", length(extra)))
    for s in sort(collect(keys(ex.a_rel)))
        p("| $s | " * join((fmt(extra[h].a_rel[s]) for h in sort(collect(keys(extra)))), " | ") * " |")
    end
    p()
    bands = Dict{Float64, Any}()
    for h in hs
        X = [r.x for r in res if r.h == h]
        n = length(X)
        K, J = size(first(X))
        lo = [quantile([x[k, j] for x in X], 0.05) for k in 1:K, j in 1:J]
        hi = [quantile([x[k, j] for x in X], 0.95) for k in 1:K, j in 1:J]
        bands[h] = (n = n, lo = lo, hi = hi, clamped = sum(r.clamped for r in res if r.h == h))
    end
    ref = reduce(vcat, permutedims.(extra[H].u_ref))
    excursion(lo, hi, x) = [maximum(max(lo[k, j] - x[k, j], x[k, j] - hi[k, j], 0.0) / abs(x[k, j])
                                    for k in axes(x, 1)) for j in axes(x, 2)]
    worst_t(lo, hi, x, j) = extra[H].t[argmax([max(lo[k, j] - x[k, j], x[k, j] - hi[k, j], 0.0) / abs(x[k, j])
                                              for k in axes(x, 1)])]
    b = bands[H]
    exc = excursion(b.lo, b.hi, ref)
    p("| pool | trajectories | reference inside the band | largest excursion outside it | at (s) | band width at the end (particles) |")
    p("|---|---|---|---|---|---|")
    for (j, s) in enumerate(pools)
        inside = count(k -> b.lo[k, j] <= ref[k, j] <= b.hi[k, j], axes(ref, 1)) / size(ref, 1)
        p(@sprintf("| %s | %d | %.1f%% of handshakes | %s | %.0f | %.3g |", s, b.n, 100inside,
                   fmt(exc[j]), worst_t(b.lo, b.hi, ref, j), (b.hi[end, j] - b.lo[end, j]) * extra[H].Ω[end]))
    end
    p()
    p("The excursion is the smallest relative observation noise below which phase 15 must exclude the pool ",
      "(spec §3, amended 2026-09-26 B). Steps a trajectory clamped at zero: $(b.clamped) over ",
      "$(b.n) trajectories.")
    p()
    if haskey(bands, H_COARSE)
        c = bands[H_COARSE]
        ec = excursion(c.lo, c.hi, ref)
        p("**The halving comparison.** Excursion at h = $H_COARSE ($(c.n) trajectories) against h = $H: ",
          join((@sprintf("`%s` %s against %s", s, fmt(ec[j]), fmt(exc[j])) for (j, s) in enumerate(pools)), "; "), ".")
        p()
    end
    # The mutation: the enolase-halved trajectory against the nominal band.
    names = Symbol[s for m in corea_models() if formalism(m) === :ode for s in states(m)]
    X = reduce(vcat, permutedims.([u[[findfirst(==(s), names) for s in pools]] for u in eno]))
    em = excursion(b.lo, b.hi, X)
    p("**The mutation.** Enolase's catalytic constants halved, same seed, against the nominal band: ",
      join((@sprintf("`%s` %s", s, fmt(em[j])) for (j, s) in enumerate(pools)), ", "),
      ". The pool it moves most is `$(pools[argmax(em)])`.")
    serialize(joinpath(OUTDIR, "band_summary.jls"), (pools = pools, excursion = exc, mutant = em))
    save("band")

# ===========================================================================
elseif SECTION == "check8"
    p("## Check 8: the summation theorems on the frozen-expression variant")
    p()
    provenance()
    p()
    cq = frozen_control_problem()
    cm = conservation_matrices(cq)
    t = @elapsed ss = steady_state(cq)
    cc = control_coefficients(cq, ss)
    dev = assert_summation(cc)
    names = cq.names[cq.internal]
    p("The variant: the four ODE modules at their nominal enzyme counts, with `FrozenDemand` consuming ",
      "charged tRNA and GTP at the published demand at the nominal pools. External lactate, the ",
      "demand's accumulator and the held protein counts are held. $(length(cq.multipliers)) ",
      "multipliers, one per reaction (spec §3, amended 2026-09-26).")
    p()
    p(@sprintf("Steady state in %.1f s: max |f| %.2e mM/s against a largest reaction rate of %.3g mM/s, after %d Newton iteration(s). Eigenvalues of the Jacobian within the conserved class: %.3g to %.3g /s, all negative.",
               t, ss.residual, ss.scale, ss.iterations, minimum(real, ss.eigenvalues), maximum(real, ss.eigenvalues)))
    p()
    p("**Conserved combinations: $(size(cm.L, 1)).** The eight of 14a's moieties the variant keeps ",
      "(redox, adenylate, guanylate, phosphate and the four carriers) and the tRNA pair each lie in ",
      "their span. The tenth is not one 14a asserts. With glucose and lactate both at the boundary it ",
      "is a closed combination over the glycolytic intermediates, NAD/NADH, the carriers and the ",
      "nucleotides, found numerically.")
    p()
    p(@sprintf("**Summation identities.** Worst concentration sum %.2e (should be 0), worst flux sum %.2e from 1. Tolerance 1e-6.",
               dev.ccc, dev.fcc))
    p("Zero-flux reactions, left out of the flux sums: ", join(("`$z`" for z in cc.zero_fluxes), ", "), ".")
    p()
    p("| multiplier | flux (mM/s, on) | Σⱼ C^J − 1 |")
    p("|---|---|---|")
    for (k, J, s, sm) in zip(cc.multipliers, cc.flux, cc.flux_species, cc.fcc_sum)
        p(@sprintf("| %s | %.4g (`%s`) | %s |", k, J, s, isnan(sm) ? "zero flux" : @sprintf("%.1e", sm - 1)))
    end
    p()
    Cg, genes = group_coefficients(cc.ccc, cc.multipliers, frozen_gene_groups())
    p("Gene-level concentration control coefficients of the smallest steady pools, ",
      "summing a gene's reactions (PGK and PYK each two):")
    p()
    small = sortperm(ss.x[[findfirst(==(s), names) for s in cc.ccc_states]])[1:5]
    p("| gene | " * join(("`$(cc.ccc_states[i])`" for i in small), " | ") * " |")
    p("|---|" * repeat("---|", length(small)))
    for (g, gname) in enumerate(genes)
        p("| $gname | " * join((fmt(Cg[i, g]) for i in small), " | ") * " |")
    end
    p()
    mut = omit_multiplier(cc, :R_charging)
    live = [k for k in eachindex(mut.multipliers) if !(mut.multipliers[k] in mut.zero_fluxes)]
    msg = try
        assert_summation(mut); "passed (it should not)"
    catch e
        sprint(showerror, e)
    end
    p(@sprintf("**The mutation.** Without `k_chg`'s multiplier the worst concentration sum is %.3g and the worst flux sum is %.3g from 1. The check reports: %s",
               maximum(abs, mut.ccc_sum), maximum(abs.(mut.fcc_sum[live] .- 1)), msg))
    serialize(joinpath(OUTDIR, "check8.jls"), (ccc = dev.ccc, fcc = dev.fcc,
              mut_ccc = maximum(abs, mut.ccc_sum), mut_fcc = maximum(abs.(mut.fcc_sum[live] .- 1)), msg = msg))
    save("check8")

# ===========================================================================
elseif SECTION == "check6"
    p("## Check 6: the nominal trajectory, and task 9.9 on the assembled model")
    p()
    provenance()
    p()
    ms = corea_models()
    d = seeded_build(ms)
    rec = run_handshake!(d, CYCLE)
    g = growth_report(d)
    p(@sprintf("Seed %d, one full cycle. Fractional growth **%.4f** (volume over its initial value).", SEED, g.fractional))
    for x in (1.02, 1.04, 1.06)
        tt = time_to_threshold(rec, x)
        p(@sprintf("- time to %.2f×: %s", x, tt === nothing ? "not reached" : @sprintf("%.0f s", tt)))
    end
    rc = only(reporting_constraints(d))
    p("- `doubling_time` is refused in code: ", rc.verdict === :refused ? "yes" : "NO",
      ", reporting instead ", join(("`$q`" for q in rc.instead), " or "), ".")
    p()

    # --- task 9.9: the two lumpings at matched cycle-mean charging flux ---
    p("### Task 9.9: the charging stoichiometry, matched on cycle-mean flux")
    p()
    counted(charging) = (ms = corea_models(); ms[4] = ChargingCounter(charging); ms)
    ci(ms, s) = findfirst(==(s), Symbol[x for m in ms if formalism(m) === :ode for x in states(m)])
    rmod = NucleotideRecycling(enzymes = :translated)
    function lumping_run(charging, seed)
        ms = counted(charging)
        d = seeded_build(ms; seed)
        names = Symbol[x for m in ms if formalism(m) === :ode for x in states(m)]
        pn = ode_pnames(ms)
        rs = [findfirst(==(s), names) for s in states(rmod)]
        ri = [findfirst(==(s), names) for s in inputs(rmod)]
        rp = [findfirst(==(q.name), pn) for q in model_free_params(parameters(rmod))]
        ratio = Float64[]; amp = Float64[]; ppi = Float64[]; adk = Float64[]
        for k in 1:CYCLE
            handshake_step!(d)
            u = d.ode.u
            if k % 60 == 0
                push!(ratio, u[ci(ms, :M_atp_c)] / u[ci(ms, :M_adp_c)])
                push!(amp, u[ci(ms, :M_amp_c)]); push!(ppi, u[ci(ms, :M_ppi_c)])
                v = recycling_fluxes(u[rs], d.ode.p[rp], d.ode.t, rmod, u[ri])
                push!(adk, v[3] * d.factor)            # R_ADK1, particles/s
            end
        end
        chg = d.ode.u[ci(ms, :chg_residues_mM)] * d.factor / d.ode.t
        return (flux = chg, ratio_mean = mean(ratio), ratio_end = ratio[end], amp = mean(amp),
                ppi = mean(ppi), adk = mean(adk))
    end
    pub = [lumping_run(TrnaCharging(), s) for s in TRNA_SEEDS]
    target = mean(r.flux for r in pub)
    two(s) = [lumping_run(AssembledTwoAtp(TwoAtpCharging(k2_scale = s)), sd) for sd in TRNA_SEEDS]
    # Secant on the two-ATP scale, matching the seed-mean cycle-mean flux.
    s0, s1 = 1.0, 1.044
    r0, r1 = two(s0), two(s1)
    f0, f1 = mean(r.flux for r in r0) - target, mean(r.flux for r in r1) - target
    for _ in 1:6
        abs(f1) <= 1e-4 * target && break
        s0, s1 = s1, s1 - f1 * (s1 - s0) / (f1 - f0)
        r0, f0 = r1, f1
        r1 = two(s1)
        f1 = mean(r.flux for r in r1) - target
    end
    p(@sprintf("Over %d paired seeds, the published AMP + PPi lumping charges %.2f residues/s on the cycle mean. The two-ATP lumping matches it at `k2_scale` = %.5f (%.2e relative).",
               length(TRNA_SEEDS), target, s1, f1 / target))
    p()
    p("| lumping | charging flux (/s) | ATP/ADP, cycle mean | ATP/ADP at the end | AMP (mM) | PPi (mM) | ADK1 net flux (/s) |")
    p("|---|---|---|---|---|---|---|")
    for (name, rs) in (("published, AMP + PPi", pub), ("two ATP → two ADP", r1))
        p(@sprintf("| %s | %.2f | %.4f | %.4f | %.4g | %.3g | %.1f |", name, mean(r.flux for r in rs),
                   mean(r.ratio_mean for r in rs), mean(r.ratio_end for r in rs),
                   mean(r.amp for r in rs), mean(r.ppi for r in rs), mean(r.adk for r in rs)))
    end
    p()
    Δ = mean(r.ratio_mean for r in r1) / mean(r.ratio_mean for r in pub) - 1
    p(@sprintf("ATP/ADP differs by **%.2f%%** between the lumpings on the cycle mean, with live glycolysis. Phase 9's double pinned ADP and measured 0.035%%.", 100Δ))
    serialize(joinpath(OUTDIR, "check6.jls"), (fractional = g.fractional, ratio_delta = Δ, k2_scale = s1))
    save("check6")

# ===========================================================================
elseif SECTION == "f5"
    p("## F5: the two external comparisons (task 14b.5)")
    p()
    provenance()
    p()
    ms = corea_models()
    tlm = only(m for m in ms if m isa CoreATranslation)
    genes = read_transcription_genes()
    jn = Symbol[x for m in ms if formalism(m) === :jump for x in states(m)]
    tx = [findfirst(==(transcript_state(g.locus)), jn) for g in genes]
    toff = findfirst(==(first(states(tlm))), jn) - 1
    ng = length(tlm.genes)
    ip = findfirst(g -> g.locus === TL_PTSG, tlm.genes)
    burn = SMOKE ? 10 : 600
    mrna = zeros(length(F5_SEEDS), length(genes))
    folds = zeros(length(F5_SEEDS), ng)
    for (si, seed) in enumerate(F5_SEEDS)
        d = seeded_build(ms; seed)
        handshake_step!(d)
        u0 = copy(d.jump.u)
        acc = zeros(length(genes)); n = 0
        for k in 2:CYCLE
            handshake_step!(d)
            if k > burn
                acc .+= d.jump.u[tx]; n += 1
            end
        end
        mrna[si, :] = acc ./ n
        p1 = collect(Float64, d.jump.u[toff .+ (1:ng)])
        p1[ip] += d.jump.u[toff + ng + 1] - u0[toff + ng + 1]
        folds[si, :] = p1 ./ u0[toff .+ (1:ng)]
        println("seed $seed done"); flush(stdout)
    end
    pred = vec(mean(mrna; dims = 1))
    meas = [g.mean_mrna for g in genes]
    tc = transcript_comparison(pred, meas)
    fl = vec(mean(folds; dims = 1))
    lens = [only(g for g in genes if g.locus === tg.locus).length for tg in tlm.genes]
    fc = fold_change_report(fl, lens)
    p("Assembled model, $(length(F5_SEEDS)) seeds ($(first(F5_SEEDS))–$(last(F5_SEEDS))), one cycle each. ",
      "Transcripts are time-averaged per gene after $(burn) s and over seeds. Fold change is each ",
      "gene's protein count at the end over its count after the first handshake, with ptsG's ",
      "cytoplasmic pool added back (phase 11's definition), averaged over seeds.")
    p()
    p("| comparison | threshold (§3) | measured | verdict |")
    p("|---|---|---|---|")
    p(@sprintf("| transcripts: Spearman | ≥ 0.7 | %.3f | %s |", tc.rho, tc.rho >= 0.7 ? "pass" : "**miss**"))
    p(@sprintf("| transcripts: within 2× | ≥ 15 of 17 | %d of %d | %s |", tc.within, tc.n, tc.within >= 15 ? "pass" : "**miss**"))
    p(@sprintf("| fold change: median | [1.7, 2.3] | %.3f | %s |", fc.median, 1.7 <= fc.median <= 2.3 ? "pass" : "**miss**"))
    p(@sprintf("| fold change: range | none < 1.5 or > 3.0 | %.3f to %.3f (%d below, %d above) | %s |",
               fc.min, fc.max, fc.below, fc.above, fc.below + fc.above == 0 ? "pass" : "**miss**"))
    p(@sprintf("| fold change against length | negative slope | %.3f | %s |", fc.slope, fc.slope < 0 ? "pass" : "**miss**"))
    p()
    p("| gene | predicted mRNA | measured | ratio | fold change | length (nt) |")
    p("|---|---|---|---|---|---|")
    for (i, g) in enumerate(genes)
        j = findfirst(tg -> tg.locus === g.locus, tlm.genes)
        p(@sprintf("| %s | %.3f | %.3f | %.2f | %s | %d |", g.locus, pred[i], meas[i], pred[i] / meas[i],
                   j === nothing ? "—" : @sprintf("%.3f", fl[j]), g.length))
    end
    p()
    p("Per-seed fold-change medians: ", join((@sprintf("%.3f", median(folds[s, :])) for s in axes(folds, 1)), ", "), ".")
    serialize(joinpath(OUTDIR, "f5.jls"), (tc = tc, fc = fc))
    save("f5")

# ===========================================================================
elseif SECTION == "report"
    parts = ["floor", "band", "check8", "check6", "f5"]
    io = IOBuffer()
    println(io, "# Phase 14b: the derivative, ensemble and census checks on the assembled Core A′")
    println(io)
    println(io, "Assembled by `dev/scripts/corea_validation_14b.jl` (SECTION=report) from the section files ",
            "below, each with its own provenance line. Check 7's census is in `corea_census_result.md`. ",
            "Regenerate with `dev/scripts/corea_validation_14b.slurm`.")
    println(io)
    for s in parts
        f = joinpath(OUTDIR, "$s.md")
        isfile(f) ? (write(io, read(f, String)); println(io)) : println(io, "*Section `$s` missing.*\n")
    end
    # F3's continuation (spec §11 task 14.10): one mutation per check, 1b, 7 and 8.
    println(io, "## F3's mutation table, continued: checks 1b, 7 and 8")
    println(io)
    println(io, "Checks 0 to 5 are in `corea_validation_result.md` (phase 14a).")
    println(io)
    println(io, "| check | mutation | what it produced | named quantity |")
    println(io, "|---|---|---|---|")
    if isfile(joinpath(OUTDIR, "band_summary.jls"))
        b = deserialize(joinpath(OUTDIR, "band_summary.jls"))
        j = argmax(b.mutant)
        println(io, @sprintf("| 1b | enolase's catalytic constants halved | excursion %.3g outside the nominal band, against %.3g nominal | `%s` |",
                             b.mutant[j], b.excursion[j], b.pools[j]))
    end
    println(io, "| 7 | tRNA pool at 0.05 mM, a fifth of the asserted | translation's counter clips (asserted in `test/test_corea_validation_14b.jl`, 600 handshakes) | `tRNA_translat` |")
    if isfile(joinpath(OUTDIR, "check8.jls"))
        c = deserialize(joinpath(OUTDIR, "check8.jls"))
        println(io, @sprintf("| 8 | `k_chg`'s multiplier omitted | worst concentration sum %.3g, worst flux sum %.3g from 1, against 1e-6 | %s |",
                             c.mut_ccc, c.mut_fcc, replace(c.msg, "|" => "/")))
    end
    write(joinpath(@__DIR__, "corea_validation_14b_result.md"), take!(io))
    println("wrote corea_validation_14b_result.md")
else
    error("SECTION must be one of floor, band, merge-band, check8, check6, f5, report; got \"$SECTION\"")
end
end

main()
