# Merge tasks 15.6 and 15.8's sensitivity run (corea_sensitivity_15.jl) into
# dev/scripts/corea_sensitivity_15_result.md (§12 2026-09-29, phase 15 planning).
#
# - 15.8: J, the six-target Jacobian in 200-cell resolution units, and its
#   split-half noise floor. The set is identifiable when all six singular values
#   exceed the floor (K6). The seven-column set adds the polymerase direction.
# - 15.6: per candidate observable and target, the largest |J| over save points,
#   against the same entry of the split-half noise matrix.
# - The metabolite panel: F2's pools less check 1b's exclusions, kept where some
#   target moves the pool by at least one resolution unit per unit ln θ, and by
#   more than three times its noise at that entry.
#
# Usage: julia --project dev/scripts/corea_sensitivity_15_merge.jl

using InferCell
using Printf
using Serialization
using Statistics: mean

const DIR = joinpath(@__DIR__, "sensitivity_15")
const OUT = joinpath(@__DIR__, "corea_sensitivity_15_result.md")
const F2 = joinpath(@__DIR__, "corea_f2_15_matrix.tsv")
const NCELLS = 200
const EXCLUDED = [:M_ppi_c, :M_13dpg_c, :M_nadh_c, :M_pep_c]   # check 1b (14b)

files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(DIR)))
isempty(files) && error("no task files in $DIR")
parts = [deserialize(joinpath(DIR, f)) for f in files]
p1 = first(parts)
all(p -> p.species == p1.species && p.times == p1.times && p.columns == p1.columns &&
         p.delta == p1.delta, parts) || error("tasks disagree on layout")
runs = merge((p.runs for p in parts)...)
seeds = sort(unique(k[1] for k in keys(runs)))
cols = p1.columns
nsp, nt = length(p1.species), length(p1.times)
rowname(r) = (p1.species[(r - 1) % nsp + 1], p1.times[(r - 1) ÷ nsp + 1])
# 14c.5's floors: 500 particles for a pool, one copy for a transcript, none for the volume.
floor_of(o) = o === :volume_litres ? 0.0 : startswith(string(o), "mRNA_") ? 1.0 : 500.0
const FLOORS = [floor_of(rowname(r)[1]) for r in 1:(nsp * nt)]

mat(ss, c, s) = permutedims(reduce(hcat, [vec(runs[(k, c, s)]) for k in ss]))
function jac(ss, cs)
    nominal = mat(seeds, :nominal, 0)    # the resolution is taken on every seed
    ensemble_jacobian([mat(ss, c, +1) for c in cs], [mat(ss, c, -1) for c in cs],
                      nominal; delta = p1.delta, ncells = NCELLS, floors = FLOORS)
end
h = length(seeds) ÷ 2
A, B = seeds[1:h], seeds[(h + 1):(2h)]
J6, N6 = jac(seeds, cols[1:6]), (jac(A, cols[1:6]) .- jac(B, cols[1:6])) ./ 2
J7, N7 = jac(seeds, cols), (jac(A, cols) .- jac(B, cols)) ./ 2
id6, id7 = check_identifiability(J6, N6), check_identifiability(J7, N7)

io = IOBuffer()
p(args...) = (println(io, args...); println(args...))
p("# Tasks 15.6 and 15.8: the ensemble sensitivity", p1.smoke ? " — SMOKE RUN, NOT A RESULT" : "")
p()
p("Merged by `dev/scripts/corea_sensitivity_15_merge.jl` from $(length(files)) task files ",
  "(runs at commit $(p1.commit), job $(p1.job); merged at ",
  "$(strip(read(`git rev-parse --short HEAD`, String)))). $(length(seeds)) paired seeds, ",
  "each through 15 configurations on the sampler model with D11's six freed: nominal, and ",
  "ln θ ± $(p1.delta) on each column. Rows are $(nsp) candidate observables × $(nt) save ",
  "points, each in units of its resolution: the $(NCELLS)-cell standard error of its nominal ",
  "mean, floored at 1% of max(|mean|, floor) as in 14c.5's gate. The ",
  "noise matrix is (J_A − J_B)/2 over two halves of $h seeds.")
p()
p("## 15.8: identifiability (K6)")
p()
fmtv(v) = join((@sprintf("%.3g", x) for x in v), ", ")
p("**D11's six** ($(join(("`$c`" for c in cols[1:6]), ", "))): singular values ",
  fmtv(id6.singular_values), @sprintf(". Noise floor %.3g. ", id6.floor),
  "Rank $(id6.rank) of 6 — **", id6.full_rank ? "identifiable, K6 does not fire" :
  "not identifiable, K6 fires", "**. ",
  @sprintf("Condition number %.3g, which the floor caps at about %.3g, so it is not held to 1e6.",
           id6.condition, first(id6.singular_values) / id6.floor))
p()
p("**With the polymerase direction** (all 17 promoters scaled together): singular values ",
  fmtv(id7.singular_values), @sprintf(". Noise floor %.3g. Rank %d of 7. ", id7.floor, id7.rank),
  @sprintf("The seventh singular value is %.3g, %.3g of the six-set's smallest. ",
           last(id7.singular_values), last(id7.singular_values) / last(id6.singular_values)),
  "With 14 promoters fixed the scale is anchored, so a partial ridge was expected (§12).")
p()

# 15.6: the influence audit, per observable and column, the largest |J| over
# save points, with the noise at that entry.
p("## 15.6: influence audit")
p()
p("Largest |∂ mean / ∂ ln θ| over the save points, in resolution units, per observable and ",
  "target. An entry is marked with * when it exceeds three times the noise at the same ",
  "entry. Observables are ordered by their largest entry.")
p()
audit = map(1:nsp) do i
    rows = [r for r in 1:size(J7, 1) if (r - 1) % nsp + 1 == i]
    vals = map(eachindex(cols)) do c
        k = argmax(abs.(J7[rows, c]))
        (v = abs(J7[rows[k], c]), sig = abs(J7[rows[k], c]) > 3abs(N7[rows[k], c]),
         t = rowname(rows[k])[2])
    end
    (o = p1.species[i], vals = vals, top = maximum(x -> x.v, vals))
end
sort!(audit; by = a -> -a.top)
short(c) = c === POLYMERASE_DIRECTION ? "polymerase" : replace(string(c), "S_JCVISYN3A_" => "S_")
p("| observable | " * join((short(c) for c in cols), " | ") * " |")
p("|---|" * repeat("---|", length(cols)))
for a in audit
    p("| `$(a.o)` | " * join((@sprintf("%.2f%s", x.v, x.sig ? "*" : "") for x in a.vals), " | ") * " |")
end
p()

# The metabolite panel.
f2rank = String[]
if isfile(F2)
    lines = readlines(F2)[2:end]
    rows = [(Symbol(first(split(l, '\t'))), maximum(abs, parse.(Float64, split(l, '\t')[2:end])))
            for l in lines]
    sort!(rows; by = r -> -r[2])
    f2rank = [string(r[1]) for r in rows]
end
proteins = Set(species_in_group(:pts))
mets = [a for a in audit if is_registered(a.o) && !(a.o in proteins) &&
        !startswith(string(a.o), "mRNA_") && a.o !== :volume_litres]
panel = [a.o for a in mets if !(a.o in EXCLUDED) && any(x -> x.v >= 1 && x.sig, a.vals)]
p("## The metabolite panel (15.5 into 15.7)")
p()
p("Kept: a pool check 1b does not exclude ($(join(("`$e`" for e in EXCLUDED), ", "))) with some ",
  "target moving it by at least one resolution unit per unit ln θ, beyond three times its noise.")
p()
p("| metabolite | F2 rank | largest entry | kept |")
p("|---|---|---|---|")
for a in mets
    k = findfirst(==(string(a.o)), f2rank)
    p("| `$(a.o)` | $(k === nothing ? "–" : k) | $(@sprintf("%.2f", a.top)) | ",
      a.o in EXCLUDED ? "excluded (check 1b)" : a.o in panel ? "yes" : "no", " |")
end
p()
p("Panel ($(length(panel))): ", join(("`$m`" for m in panel), ", "), ".")
write(OUT, take!(io))
serialize(joinpath(DIR, "merged.jls"), (J6 = J6, N6 = N6, J7 = J7, N7 = N7, species = p1.species,
                                        times = p1.times, columns = cols, panel = panel))
println("Wrote $OUT")
