# KIM-195: SpeedyCell's own values at procid 1's initial state, the reference that
# route_cost_meanfield.jl's `check` stage compares its generic code against.
#
#   julia --project=dev/scripts/fullcell_timing_cycle \
#       dev/scripts/fullcell_route_cost/reference_speedycell.jl <outdir>
#
# Runs in KIM-183's environment (SpeedyCell from syn3a-ensemble by path; Julia 1.11.9).
# It writes <outdir>/reference/:
#   rates.csv     the CME rate constants from the rate rules at the initial counts
#                 (SpeedyCell.CME.evaluate_rates!), one per reaction
#   p.csv, y0.csv the hook's enzyme concentrations and ODE initial state
#                 (SpeedyCell.ODE.hook_inputs)
#   flux.csv      the 175 ODE fluxes at (y0, p), published rates (SpeedyCell.ODE.flux!)
#   flux_params.csv  the same through the parameterised rates at the default constants
#                 (SpeedyCell.ODE.flux_params!)
#   x1.csv, rates1.csv, p1.csv, y1.csv, flux1.csv, V1.txt
#                 the same at x1 = x0 .+ 1, where every reaction and every counter is
#                 live: x1's surface areas from SpeedyCell's lipid_sa! (round and trunc)
#                 and its CellV from cell_volume, then the rates, the hook inputs, the
#                 fluxes and the volume V (L) there
#   pay.csv       the hook's payments as SpeedyCell's WriteBack lists them: kind,
#                 counter, pool, then the species it adds to (1-based species indices)
#   meta.txt      host, Julia version, the syn3a-ensemble commit
# Every number is written with repr, so it round-trips exactly.

using SpeedyCell
const CME = SpeedyCell.CME
const ODE = SpeedyCell.ODE
using InteractiveUtils: versioninfo

out = joinpath(abspath(ARGS[1]), "reference")
mkpath(out)
model = "/fred/oz022/tkimpson/syn3a-ensemble/data/model"
names = String[]; x0 = Int[]
for l in Iterators.drop(eachline(joinpath(model, "cme", "species.csv")), 1)
    i, nm, c = split(l, ',')
    push!(names, nm); push!(x0, parse(Int, c))
end
@assert names == ODE.CME_NAMES

rules = CME.RateRules(model, names)
rates = CME.evaluate_rates!(zeros(CME.nrules(rules)), rules, x0)
p, y0 = ODE.hook_inputs(x0)
V = zeros(ODE.N_REACTIONS); ODE.flux!(V, y0, p)
Vp = zeros(ODE.N_REACTIONS); ODE.flux_params!(Vp, y0, ODE.RateParams(p, ODE.K_DEFAULT))

wr(f, v) = open(joinpath(out, f), "w") do io
    for x in v
        println(io, repr(x))
    end
end
# x1: every count one above x0, with the surface areas and CellV the hook writes.
wb = SpeedyCell.Cell.WriteBack(names)
x1 = x0 .+ 1
SpeedyCell.Cell.lipid_sa!(x1, wb)
V1, cellv1 = SpeedyCell.Cell.cell_volume(x1[wb.cellsa]); x1[wb.cellv] = cellv1
rates1 = CME.evaluate_rates!(zeros(CME.nrules(rules)), rules, x1)
p1, y1 = ODE.hook_inputs(x1)
F1 = zeros(ODE.N_REACTIONS); ODE.flux!(F1, y1, p1)
pay = String[]
for (c, translat) in wb.hydro
    pool, made = translat ? (wb.gtp, wb.gdp) : (wb.atp, wb.adp)
    push!(pay, "hydro,$c,$pool,$made;$(wb.pi)")
end
for list in (wb.ntp, wb.dntp), (c, m) in list
    push!(pay, "ntp,$c,$m,$(wb.ppi)")
end
for (c, m) in wb.nmp
    push!(pay, "nmp,$c,0,$m")
end
for (c, ch, un) in wb.aatrna
    push!(pay, "aatrna,$c,$ch,$un")
end

wr("rates.csv", rates); wr("p.csv", p); wr("y0.csv", y0); wr("flux.csv", V); wr("flux_params.csv", Vp)
wr("x1.csv", x1); wr("rates1.csv", rates1); wr("p1.csv", p1); wr("y1.csv", y1); wr("flux1.csv", F1)
wr("V1.txt", [V1])
open(io -> foreach(l -> println(io, l), pay), joinpath(out, "pay.csv"), "w")
open(joinpath(out, "meta.txt"), "w") do io
    println(io, "host ", gethostname(), " job ", get(ENV, "SLURM_JOB_ID", "none"))
    println(io, "syn3a-ensemble ", readchomp(`git -C /fred/oz022/tkimpson/syn3a-ensemble rev-parse HEAD`))
    println(io, "reactions ", length(rates), " enzymes ", length(p), " ode species ", length(y0), " fluxes ", length(V))
    println(io, "flux_params == flux bitwise: ", Vp == V)
    versioninfo(io)
end
println("wrote ", out)
