# KIM-116: check that item 8's D3 draws are item 7's posterior. The predict prep
# (corea_modular_kim116_predict.jl) assembles the plug-in's cut draws from stage2_*.jls
# itself; this checks its log_draws, and their coordinate names, against the merge's
# run/cut_draws_plugin.tsv, bit for bit.
#
# Usage: julia --project dev/scripts/corea_modular_kim116_d3check.jl <kim116-archive>
# (first run inline as Slurm 21628002, at e129fcc, with this code)

include(joinpath(@__DIR__, "corea_modular_kim116_predict.jl"))
using DelimitedFiles: readdlm

function d3check(A)
    pp = load(joinpath(A, "predict"), "prep")
    tsv = joinpath(A, "run", "cut_draws_plugin.tsv")
    T = readdlm(tsv, Float64; skipstart = 1)
    h = split(readline(tsv), '\t')
    println("names equal: ", Symbol.(h) == pp.names)
    println("size ", size(pp.log_draws), " vs ", size(T))
    println("bitwise equal: ", pp.log_draws == T)
    println("max abs diff: ", maximum(abs.(pp.log_draws .- T)))
    Symbol.(h) == pp.names && pp.log_draws == T || error("the D3 draws are not the merge's plug-in draws")
end

d3check(ARGS[1])
