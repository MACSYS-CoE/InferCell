# KIM-56: a parameter vector per run, for syn3a-ensemble's SpeedyCell (`Params`, lane 10a).
#
# `params_from_vector(model, v)` takes one multiplier per named kinetic constant, in
# `param_names()` order (805 ODE constants, then 33 CME constants; see
# results/param_names.csv), and returns the `Params` that `run_cell(...; params=P)` takes.
# v = ones(838) is the published model on the parameterised path. Validation (finite,
# ≥ 0) is Params's.
#
# Loading this file checks that the SpeedyCell in use is syn3a-ensemble's at the commit
# in SYN3A_PIN, with code/SpeedyCell and data/model clean, and that `param_names()` is the
# committed results/param_names.csv, name for name and in order. The CME names carry
# source line numbers, so another commit can reorder them with the count still 838; a
# vector indexed by the committed table would then scale the wrong constants.

import SpeedyCell
using SpeedyCell.Cell

const SYN3A = "/fred/oz022/tkimpson/syn3a-ensemble"
const SYN3A_PIN = readchomp(joinpath(@__DIR__, "SYN3A_PIN"))

function check_pin(; root=SYN3A, pin=SYN3A_PIN)
    src = joinpath(root, "code", "SpeedyCell")
    startswith(pathof(SpeedyCell), src) || error("SpeedyCell is loaded from $(pathof(SpeedyCell)), not $src")
    head = readchomp(`git -C $root rev-parse HEAD`)
    head == pin || error("syn3a-ensemble is at $head, not the pin $pin")
    st = readchomp(`git -C $root status --porcelain -- code/SpeedyCell data/model`)
    isempty(st) || error("syn3a-ensemble code/SpeedyCell or data/model is not clean:\n$st")
    return nothing
end

function check_names(csv=joinpath(@__DIR__, "results", "param_names.csv"); names=param_names())
    lines = collect(Iterators.drop(eachline(csv), 1))
    any(l -> occursin('"', l), lines) && error("$csv has a quoted name; this check splits on commas")
    committed = [split(l, ',')[2] for l in lines]
    committed == names && return nothing
    length(committed) == length(names) || error("param_names() has $(length(names)) names, $csv $(length(committed))")
    i = findfirst(k -> committed[k] != names[k], eachindex(names))
    error("param_names() differs from $csv at index $i: $(names[i]) against $(committed[i])")
end

check_pin()
check_names()

function params_from_vector(m::CellModel, v::AbstractVector{<:Real}; names::AbstractVector{String}=param_names())
    length(v) == length(names) ||
        throw(DimensionMismatch("$(length(v)) multipliers for $(length(names)) named constants"))
    return Params(m; multipliers=Dict{String,Float64}(n => Float64(x) for (n, x) in zip(names, v)))
end
