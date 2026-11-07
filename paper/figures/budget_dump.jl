# Writes the reduced model's prediction, per parameter pair, from the KIM-43 evidence archive's
# w<i>.jls files, to paper/figures/data/budget_mu.tsv: width, the θ_I row, the θ_H column, μ̂
# (the mean of Q over 4 hybrid cells) and whether the viability rule kept the pair.
# make_figures.py's fig_budget reads the output. It uses no InferCell types.
#
# Run from the repository root, with no project:
#   julia paper/figures/budget_dump.jl [EVIDENCE_ARCHIVE] > paper/figures/data/budget_mu.tsv
using Serialization

function main(archive)
    println("width\trow\tcolumn\tmu\tviable")
    for i in 0:2
        x = deserialize(joinpath(archive, "runs/budget-kim42/w$i.jls")).result
        for j in CartesianIndices(x.mu)
            println(join((x.width, j[1], j[2], x.mu[j], x.viable[j]), '\t'))
        end
    end
end

main(isempty(ARGS) ? "/fred/oz022/tkimpson/infercell-predict-archive-2026-11-06" : ARGS[1])
