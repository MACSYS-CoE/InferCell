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
wr("rates.csv", rates); wr("p.csv", p); wr("y0.csv", y0); wr("flux.csv", V); wr("flux_params.csv", Vp)
open(joinpath(out, "meta.txt"), "w") do io
    println(io, "host ", gethostname(), " job ", get(ENV, "SLURM_JOB_ID", "none"))
    println(io, "syn3a-ensemble ", readchomp(`git -C /fred/oz022/tkimpson/syn3a-ensemble rev-parse HEAD`))
    println(io, "reactions ", length(rates), " enzymes ", length(p), " ode species ", length(y0), " fluxes ", length(V))
    println(io, "flux_params == flux bitwise: ", Vp == V)
    versioninfo(io)
end
println("wrote ", out)
