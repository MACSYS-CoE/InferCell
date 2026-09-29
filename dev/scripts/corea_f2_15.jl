# Spec §11 task 15.5: output F2, the concentration control coefficient matrix,
# with check 8 as its guard, and the metabolite panel chosen from its largest
# rows (spec §4 D9). Writes dev/scripts/corea_f2_15_result.md and the matrix as
# dev/scripts/corea_f2_15_matrix.tsv.
#
# Usage: julia --project dev/scripts/corea_f2_15.jl
# Regenerate: sbatch dev/scripts/corea_f2_15.slurm

using InferCell
using Printf
using Dates

const OUT = joinpath(@__DIR__, "corea_f2_15_result.md")
const TSV = joinpath(@__DIR__, "corea_f2_15_matrix.tsv")
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src dev/scripts/corea_f2_15.jl`, String)))

io = IOBuffer()
p(args...) = (println(io, args...); println(args...))
p("# Output F2: concentration control coefficients (spec §11 task 15.5)")
p()
p("Commit $COMMIT$(DIRTY ? " (tree dirty)" : ""), job $(get(ENV, "SLURM_JOB_ID", "none")), ",
  "$(now()), Julia $VERSION. Regenerate with `sbatch dev/scripts/corea_f2_15.slurm`.")
p()
t = @elapsed f2 = concentration_control()
p(@sprintf("The frozen-expression variant's steady state, %d metabolite pools against %d genes, in %.1f s. ",
           length(f2.metabolites), length(f2.genes), t),
  @sprintf("Check 8 guards it: the worst summation deviation is %.2e (tolerance 1e-6), ", maximum(f2.summation)),
  @sprintf("and every coefficient matches re-solved steady states to %.2e (concentration) ", f2.derivatives.ccc),
  @sprintf("and %.2e (flux), against 1e-5.", f2.derivatives.fcc))
p()
p("Each entry is ∂ln x/∂ln E for a gene's enzyme, summing its reactions (PGK and PYK each two). ",
  "The four phosphotransferase genes set carrier totals, which multiply no rate, and are not columns. ",
  "The carrier states are protein counts and are not rows (spec §4 D8).")
p()
fmt(x) = abs(x) < 5e-4 ? "0" : @sprintf("%.3f", x)
p("| metabolite | steady (mM) | " * join(("$g" for g in f2.genes), " | ") * " |")
p("|---|---|" * repeat("---|", length(f2.genes)))
order = sortperm([maximum(abs, f2.C[i, :]) for i in eachindex(f2.metabolites)]; rev = true)
for i in order
    p("| `$(f2.metabolites[i])` | $(@sprintf("%.4g", f2.steady[i])) | " *
      join((fmt(c) for c in f2.C[i, :]), " | ") * " |")
end
p()
open(TSV, "w") do f
    println(f, join(vcat("metabolite", string.(f2.genes)), '\t'))
    for i in eachindex(f2.metabolites)
        println(f, join(vcat(string(f2.metabolites[i]), repr.(f2.C[i, :])), '\t'))
    end
end
panel = metabolite_panel(f2)
p("## The metabolite panel")
p()
p("Ranked by the row's largest |coefficient| over the genes. The panel is every pool at or above ",
  "0.1, where a 10% change in some enzyme moves it by at least 1%.")
p()
p("| rank | metabolite | largest \\|C\\| | gene | in panel |")
p("|---|---|---|---|---|")
for (k, r) in enumerate(panel.ranking)
    p(@sprintf("| %d | `%s` | %.3f | %s | %s |", k, r.metabolite, r.strength, r.gene,
               r.metabolite in panel.panel ? "yes" : "no"))
end
p()
p("Panel ($(length(panel.panel))): ", join(("`$m`" for m in panel.panel), ", "), ".")
write(OUT, take!(io))
println("Wrote $OUT")
