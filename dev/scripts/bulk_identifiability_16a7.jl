# spec/phases/16-recovery.md task 16a.7, after KP fired twice (handoff
# 2026-10-03): is F6 still possible when metabolites are observed as
# population averages with a realistic assay noise?
#
# This reuses 15.8's sensitivity run (job 17704779 at 7835c20; the sampler
# model, D11's six freed, 1470 paired seeds), and runs no new simulation. Each
# row of 15.8's Jacobian is an observable's 200-cell ensemble mean at one save
# point, in units of its resolution. 15.8's resolution is the standard error of
# that mean, floored at 1% of max(|mean|, floor). That is a population-average
# design whose metabolite precision is set by cell-to-cell spread alone.
#
# Here a metabolite row's floor is a bulk assay's coefficient of variation,
# CV ∈ {1% (15.8's), 5%, 10%, 20%}. A transcript row keeps 15.8's resolution, an
# ensemble mean of exact counts, which understates the per-cell counts the data
# carry. Save point t = 0 is dropped. Per design and CV:
# - the noise-gated rank of J (15.8's test: singular values above the largest
#   singular value of the split-half noise matrix);
# - F6's quantity, posterior over prior SD per target, linearised: the Fisher
#   information JᵀJ (unit noise per row) plus each target's lognormal prior
#   precision in ln θ;
# - the same SD ratio with the noise matrix in place of J. A target whose
#   ratio is no lower under J than under the noise has not shrunk.
#
# Rows at different save points are treated as independent. That overstates
# the information, most for the cell-to-cell part of a row's resolution, since
# the same cells are followed. The "sparse" design keeps every 15th save point
# (15 min) to show how much rests on that.
#
# Usage: julia --project dev/scripts/bulk_identifiability_16a7.jl <dir of 15.8's task files>
# Regenerate: sbatch dev/scripts/bulk_identifiability_16a7.slurm <dir>

using InferCell
using Printf
using Serialization
using LinearAlgebra
using Statistics: mean
using Distributions: LogNormal, params

const DIR = ARGS[1]
const NCELLS = 200
const CVS = [0.01, 0.05, 0.10, 0.20]

function main()
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(DIR)))
    isempty(files) && error("no task files in $DIR")
    parts = [deserialize(joinpath(DIR, f)) for f in files]
    p1 = first(parts)
    all(p -> p.species == p1.species && p.times == p1.times && p.columns == p1.columns &&
             p.commit == p1.commit && p.job == p1.job, parts) || error("tasks disagree")
    runs = merge((p.runs for p in parts)...)
    seeds = sort(unique(k[1] for k in keys(runs)))
    cols = p1.columns[1:6]
    nsp, nt = length(p1.species), length(p1.times)
    sp(r) = p1.species[(r - 1) % nsp + 1]
    ti(r) = (r - 1) ÷ nsp + 1
    floor_of(o) = o === :volume_litres ? 0.0 : startswith(string(o), "mRNA_") ? 1.0 : 500.0
    floors = [floor_of(sp(r)) for r in 1:(nsp * nt)]
    mat(ss, c, s) = permutedims(reduce(hcat, [vec(runs[(k, c, s)]) for k in ss]))
    nominal = mat(seeds, :nominal, 0)
    jac(ss, rel) = ensemble_jacobian([mat(ss, c, +1) for c in cols], [mat(ss, c, -1) for c in cols],
                                     nominal; delta = p1.delta, ncells = NCELLS, rel, floors)
    h = length(seeds) ÷ 2
    A, B = seeds[1:h], seeds[(h + 1):(2h)]

    met = Set(M0_PANEL)
    tx = Set(s for s in p1.species if startswith(string(s), "mRNA_"))
    live = [r for r in 1:(nsp * nt) if p1.times[ti(r)] > 0]
    sparse = Set(t for (i, t) in enumerate(sort(unique(ti(r) for r in live))) if i % 15 == 0)
    metrows = [r for r in live if sp(r) in met]
    txrows = [r for r in live if sp(r) in tx]

    # Each target's prior SD in ln θ.
    ms = d11_models()
    prior_sd = map(cols) do c
        ps = [p for m in ms for p in parameters(m) if p.name === c]
        length(ps) == 1 || error("no unique parameter $c")
        pr = only(ps).prior
        pr isa LogNormal || error("$c's prior is $(typeof(pr)), not lognormal")
        params(pr)[2]
    end

    ratio(J) = sqrt.(diag(inv(J' * J + Diagonal(1 ./ prior_sd .^ 2)))) ./ prior_sd

    println("# Bulk identifiability: is F6 possible with population-average metabolites?\n")
    @printf("From %d task files of 15.8's run (job %s at %s): %d seeds, δ = %g, %d save points after t = 0.\n",
            length(files), p1.job, p1.commit, length(seeds), p1.delta, length(sort(unique(ti(r) for r in live))))
    println("Prior SD in ln θ: ", join([@sprintf("%s %.3g", c, s) for (c, s) in zip(cols, prior_sd)], ", "), "\n")
    short = ["S_0607", "S_0445", "S_0779 (ptsG)", "krnadeg", "ENO", "FBA"]
    J1, NA1, NB1 = jac(seeds, 0.01), jac(A, 0.01), jac(B, 0.01)
    for cv in CVS
        Jc, NAc, NBc = cv == 0.01 ? (J1, NA1, NB1) : (jac(seeds, cv), jac(A, cv), jac(B, cv))
        # Metabolite rows at the bulk CV; transcript rows at 15.8's resolution.
        J = vcat(Jc[metrows, :], J1[txrows, :])
        N = vcat((NAc[metrows, :] .- NBc[metrows, :]) ./ 2, (NA1[txrows, :] .- NB1[txrows, :]) ./ 2)
        nm = length(metrows)
        @printf("## Metabolite assay CV %.0f%%\n\n", 100cv)
        println("| Design | Rows | Rank (floor) | ", join(short, " | "), " |")
        println("|---|---|---|", repeat("---|", 6))
        for (name, rows) in [("metabolites only", 1:nm), ("transcripts only", (nm + 1):size(J, 1)),
                             ("joint", 1:size(J, 1)),
                             ("metabolites only, sparse", [i for (i, r) in enumerate(metrows) if ti(r) in sparse]),
                             ("transcripts only, sparse", [nm + i for (i, r) in enumerate(txrows) if ti(r) in sparse])]
            Jd, Nd = J[rows, :], N[rows, :]
            id = check_identifiability(Jd, Nd)
            rJ, rN = ratio(Jd), ratio(Nd)
            cells = [@sprintf("%.2f (%.2f)", a, b) for (a, b) in zip(rJ, rN)]
            @printf("| %s | %d | %d of 6 (%.3g) | %s |\n", name, length(rows), id.rank, id.floor,
                    join(cells, " | "))
        end
        println("\nEach cell is posterior over prior SD, linearised, with the same ratio from the noise matrix in brackets. ",
                "Lower is more shrinkage; a cell no lower than its bracket has not shrunk.\n")
    end
end

main()
