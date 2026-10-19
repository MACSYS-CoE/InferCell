# Issue #88, tasks A1 and A2: merge `sensitivity_88.jl`'s task files into
#
# 1. the linearised information bound: posterior over prior SD per target under
#    metabolites only, transcripts only and joint, beside its split-half noise
#    bracket, at bulk assay CV σ_b. This is `bulk_identifiability_16a7.jl`'s
#    calculation (D9) on any of the two models. Each row is an observable's
#    C-cell ensemble mean at one save point; a metabolite row's resolution is
#    max(SE over C cells, σ_b · max(|mean|, floor)), a transcript row's 15.8's
#    (σ_b = 1%). Rows whose SE exceeds σ_b · mean are counted: there the σ_b
#    sweep is flat, because cell-to-cell spread, not the assay, sets the noise.
# 2. R per pool, the variance of the bulk log-likelihood over the C cells' paths
#    (inference-framing note §7): R ≈ Σ_i Var(log x̄_i) / σ_b², with x̄ the C-cell
#    mean of a panel pool at a save point, floored at one particle as the
#    likelihood floors it (src/csmc.jl). Var(log x̄) is estimated by resampling C
#    of the nominal cells, not from CV²/C, which fails at the one-particle floor.
#    The cells are unconditioned on the transcripts, so this R is an upper bound
#    on the note's R, which conditions on D_tx.
# 3. For `corea` only, 15.8's identifiability check redone (every observable,
#    every save, 1% floor, 200 cells), against 15.8's singular values.
#
# Save point t = 0, the fixed initial condition, is dropped from 1 and 2.
#
# Usage: julia --project dev/scripts/sensitivity_88_merge.jl <dir of task files>
# Writes <dir>/merged.md and <dir>/merged.jls, and prints the markdown.

using InferCell
using Printf
using Random
using Serialization
using LinearAlgebra
using Statistics: mean, std, var, median
using Distributions: LogNormal, params

const DIR = ARGS[1]
const NRESAMPLE = 4000
const D158 = [923, 400, 254, 215, 184, 119]   # 15.8's singular values, floor 72.3
const SHORT = Dict(promoter_param(:JCVISYN3A_0607) => "S_0607 (GAPD)",
                   promoter_param(:JCVISYN3A_0445) => "S_0445 (PGI)",
                   promoter_param(:JCVISYN3A_0779) => "S_0779 (ptsG)",
                   :krnadeg => "krnadeg", :kcatF_R_ENO => "ENO", :kcatF_R_FBA => "FBA")

function load(dir)
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(dir)))
    isempty(files) && error("no task files in $dir")
    parts = [deserialize(joinpath(dir, f)) for f in files]
    p1 = first(parts)
    all(p -> p.species == p1.species && p.times == p1.times && p.columns == p1.columns &&
             p.delta == p1.delta && p.model == p1.model, parts) || error("tasks disagree on layout")
    all(p -> p.commit == p1.commit && p.job == p1.job && p.smoke == p1.smoke, parts) ||
        error("the task files come from more than one run")
    return p1, merge((p.runs for p in parts)...), length(files)
end

function prior_sds(ms, cols)
    map(cols) do c
        ps = [p for m in ms for p in parameters(m) if p.name === c]
        length(ps) == 1 || error("no unique parameter $c")
        pr = only(ps).prior
        pr isa LogNormal || error("$c's prior is $(typeof(pr)), not lognormal")
        params(pr)[2]
    end
end

function main()
    p1, runs, nfiles = load(DIR)
    model = p1.model
    seeds = sort(unique(k[1] for k in keys(runs)))
    cols = p1.columns
    nsp, nt = length(p1.species), length(p1.times)
    sp(r) = p1.species[(r - 1) % nsp + 1]
    ti(r) = (r - 1) ÷ nsp + 1
    istx(o) = startswith(string(o), "mRNA_")
    floor_of(o) = o === :volume_litres ? 0.0 : istx(o) ? 1.0 : 500.0
    floors = [floor_of(sp(r)) for r in 1:(nsp * nt)]
    mat(ss, c, s) = permutedims(reduce(hcat, [vec(runs[(k, c, s)]) for k in ss]))
    nominal = mat(seeds, :nominal, 0)
    h = length(seeds) ÷ 2
    A, B = seeds[1:h], seeds[(h + 1):(2h)]

    ms = model === :m0 ? m0_models() : d11_models()
    prior_sd = prior_sds(ms, cols)
    ratio(J) = sqrt.(diag(inv(J' * J + Diagonal(1 ./ prior_sd .^ 2)))) ./ prior_sd

    live = [r for r in 1:(nsp * nt) if p1.times[ti(r)] > 0]
    metrows = [r for r in live if sp(r) in Set(M0_PANEL)]
    txrows = [r for r in live if istx(sp(r))]
    ncells = model === :m0 ? 50 : 200
    cvs = model === :m0 ? [0.05, 0.10, 0.20] : [0.01, 0.05, 0.10, 0.20]

    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    p("# Issue #88 ", model === :m0 ? "A1: M0" : "A2: Core A′", " sensitivity merge",
      p1.smoke ? " — SMOKE RUN, NOT A RESULT" : "")
    p()
    p("From $nfiles task files (job $(p1.job) at $(p1.commit)), merged at ",
      strip(read(`git rev-parse --short HEAD`, String)), ": $(length(seeds)) paired seeds, ",
      "δ = $(p1.delta), horizon $(p1.horizon) s, $(length(unique(ti.(live)))) save points after t = 0. ",
      "Split halves of $h seeds.")
    p()
    p("Prior SD in ln θ: ", join([@sprintf("%s %.3g", get(SHORT, c, string(c)), s)
                                  for (c, s) in zip(cols, prior_sd)], ", "), ".")
    p()

    # 1. The information bound.
    jac(ss, rel, fl) = ensemble_jacobian([mat(ss, c, +1) for c in cols], [mat(ss, c, -1) for c in cols],
                                         nominal; delta = p1.delta, ncells, rel, floors = fl)
    names = [get(SHORT, c, string(c)) for c in cols]
    m = vec(mean(nominal; dims = 1)); s = vec(std(nominal; dims = 1))
    J1, NA1, NB1 = jac(seeds, 0.01, floors), jac(A, 0.01, floors), jac(B, 0.01, floors)
    bound = Dict{Tuple{Float64, String}, Any}()
    p("## Information bound ($ncells cells)")
    p()
    for (cv, fl, tag) in vcat([(cv, floors, "") for cv in cvs],
                              [(0.10, [istx(sp(r)) || sp(r) === :volume_litres ? floors[r] : 1.0
                                       for r in eachindex(floors)], ", metabolite floor 1 particle")])
        Jc, NAc, NBc = (jac(seeds, cv, fl), jac(A, cv, fl), jac(B, cv, fl))
        J = vcat(Jc[metrows, :], J1[txrows, :])
        N = vcat((NAc[metrows, :] .- NBc[metrows, :]) ./ 2, (NA1[txrows, :] .- NB1[txrows, :]) ./ 2)
        nm = length(metrows)
        sedom = count(r -> s[r] / sqrt(ncells) > cv * max(abs(m[r]), fl[r]), metrows)
        p(@sprintf("### σ_b = %.0f%%%s", 100cv, tag))
        p()
        p("$sedom of $nm metabolite rows have their resolution set by the $ncells-cell SE, not σ_b.")
        p()
        p("| Design | Rows | Rank (floor) | ", join(names, " | "), " |")
        p("|---|---|---|", repeat("---|", length(cols)))
        for (name, rows) in [("metabolites only", 1:nm), ("transcripts only", (nm + 1):size(J, 1)),
                             ("joint", 1:size(J, 1))]
            Jd, Nd = J[rows, :], N[rows, :]
            id = check_identifiability(Jd, Nd)
            rJ, rN = ratio(Jd), ratio(Nd)
            bound[(cv, name * tag)] = (ratio = rJ, bracket = rN, rank = id.rank, floor = id.floor)
            p(@sprintf("| %s | %d | %d of %d (%.3g) | ", name, length(rows), id.rank, length(cols), id.floor),
              join([@sprintf("%.2f (%.2f)", a, b) for (a, b) in zip(rJ, rN)], " | "), " |")
        end
        p()
    end
    p("Each cell is posterior over prior SD, linearised, with the split-half noise matrix's ratio in ",
      "brackets. A target is informed only if its ratio is below its bracket. Save points are treated ",
      "as independent, which overstates the information.")
    p()
    if model === :m0
        k = findfirst(==(promoter_param(:JCVISYN3A_0779)), cols)
        b = bound[(0.10, "metabolites only")]
        p(@sprintf("**A1 gate (ptsG, metabolites only, σ_b = 10%%):** ratio %.2f, bracket %.2f — %s.",
                   b.ratio[k], b.bracket[k],
                   b.ratio[k] < b.bracket[k] ? "informed: the gate passes" :
                   "not informed: the gate fails; choose a larger M0 dataset"))
        p()
    end

    # 2. R per pool.
    pools = [o for o in M0_PANEL if o in p1.species]
    tlive = [t for t in 1:nt if p1.times[t] > 0]
    nom = Dict(k => runs[(k, :nominal, 0)] for k in seeds)
    rng = Xoshiro(88)
    Rtab = Dict{Tuple{Symbol, Int}, Vector{Float64}}()   # pool, C => Var(log x̄) per save
    cvcell = Dict{Symbol, Float64}()
    below1 = Dict{Symbol, Float64}()
    for o in pools
        i = findfirst(==(o), p1.species)
        X = [nom[k][i, t] for k in seeds, t in tlive]     # cells × saves
        cvcell[o] = median(std(X; dims = 1) ./ max.(mean(X; dims = 1), 1e-12))
        below1[o] = mean(X .< 1)
        for C in (model === :m0 ? (50,) : (200, 50))
            Rtab[(o, C)] = map(axes(X, 2)) do t
                var([log(max(mean(X[rand(rng, 1:length(seeds), C), t]), 1.0)) for _ in 1:NRESAMPLE])
            end
        end
    end
    Cs = model === :m0 ? (50,) : (200, 50)
    sbs = (0.05, 0.10, 0.20)
    p("## R per pool")
    p()
    p("Σ over $(length(tlive)) save points of Var(log x̄_C) / σ_b², x̄_C the mean of C cells resampled ",
      "($NRESAMPLE draws) from $(length(seeds)) nominal cells, floored at 1 particle. Per-cell CV is the ",
      "median over save points. \"< 1\" is the fraction of cell-saves below one particle.")
    p()
    p("| Pool | per-cell CV | < 1 | ", join(["R, C = $C, σ_b = $(Int(100sb))%" for C in Cs for sb in sbs], " | "), " |")
    p("|---|---|---|", repeat("---|", length(Cs) * length(sbs)))
    for o in sort(pools; by = o -> -sum(Rtab[(o, first(Cs))]))
        p("| `$o` | ", @sprintf("%.3g | %.2f | ", cvcell[o], below1[o]),
          join([@sprintf("%.3g", sum(Rtab[(o, C)]) / sb^2) for C in Cs for sb in sbs], " | "), " |")
    end
    p("| **total** | | | ", join([@sprintf("**%.3g**", sum(sum(Rtab[(o, C)]) for o in pools) / sb^2)
                                   for C in Cs for sb in sbs], " | "), " |")
    p()
    p("n = $(length(pools)) pools × $(length(tlive)) saves = $(length(pools) * length(tlive)). Regimes ",
      "(note §7): R ≪ 1, one path draw per θ suffices; 1 ≲ R ≲ n, Gaussian convolution; pools near ",
      "the one-particle floor are non-Gaussian.")
    p()

    # 3. 15.8's check, redone (Core A′ only).
    regress = nothing
    if model === :corea
        regress = check_identifiability(J1, (NA1 .- NB1) ./ 2)
        p("## 15.8's identifiability check, redone")
        p()
        p("Every observable and save, 1% floor, 200 cells: singular values ",
          join((@sprintf("%.3g", x) for x in regress.singular_values), ", "),
          @sprintf(", floor %.3g, rank %d of 6. ", regress.floor, regress.rank),
          "15.8 (job 17704779 at 7835c20): ", join(D158, ", "), ", floor 72.3, rank 6 of 6.")
        p()
    end

    write(joinpath(DIR, "merged.md"), take!(io))
    serialize(joinpath(DIR, "merged.jls"),
              (; model, commit = p1.commit, job = p1.job, nseeds = length(seeds), cols, prior_sd,
               bound, Rtab, cvcell, below1, pools, times = p1.times[tlive], regress, smoke = p1.smoke))
    println("Wrote ", joinpath(DIR, "merged.md"))
end

main()
