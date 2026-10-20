# Issue #89 B3, option 1's test: a reference-θ offset δ taken at one truth and applied to
# another, scored as KG. Reads the B3 archive and A3's live60 ensembles; no new runs.
# Usage: julia --project dev/scripts/jensen_gap_89_delta.jl
using InferCell, Serialization, Statistics, Printf
D = "/fred/oz022/tkimpson/infercell-89-b3-archive-2026-10-20"
A3 = "/fred/oz022/tkimpson/infercell-88-a3-archive-2026-10-19"
loadp(dir) = [deserialize(joinpath(dir, f)) for f in sort(readdir(dir)) if endswith(f, ".jls")]
dirs = Dict("m0_v6" => "$D/ensemble_m0_v6", "m0_low" => "$D/ensemble_m0_low", "m0_high" => "$D/ensemble_m0_high",
            "corea1" => "$A3/ensemble_live60_truth1", "corea2" => "$A3/ensemble_live60_truth2")
lb(cells, j) = [log(max(mean(c[i, j] for c in cells), 1.0)) for i in eachindex(M0_PANEL)]
σ = 0.10
data = Dict()
for (t, dir) in dirs
    ps = loadp(dir); sp = first(ps).species; times = first(ps).times
    rows = [findfirst(==(o), sp) for o in M0_PANEL]
    cells = [c[rows, :] for p in ps for c in p.cells]
    mf = deserialize("$D/meanfield_$t.jls")
    common = [τ for τ in mf.times if τ > 0 && τ in times]
    jr = [findfirst(==(τ), times) for τ in common]; jm = [findfirst(==(τ), mf.times) for τ in common]
    h = length(cells) ÷ 2
    L = reduce(hcat, [lb(cells, j) for j in jr])
    H = reduce(hcat, [lb(cells[1:h], j) for j in jr]) .- reduce(hcat, [lb(cells[(h+1):2h], j) for j in jr])
    M = log.(max.(mf.counts[:, jm], 1.0))
    data[t] = (; L, M, noise = vec(sum(H .^ 2; dims = 2)) ./ 2 ./ (2σ^2))
end
gtp = findfirst(==(:M_gtp_c), M0_PANEL)
println("| Offset from → applied to | total | noise | pools > 3× noise | GTP shift (noise) |")
println("|---|---|---|---|---|")
for group in (["m0_v6", "m0_low", "m0_high"], ["corea1", "corea2"]), a in group, b in group
    a == b && continue
    δ = data[a].L .- data[a].M
    sh = vec(sum((data[b].M .+ δ .- data[b].L) .^ 2; dims = 2)) ./ (2σ^2)
    n = data[b].noise
    bad = [String(M0_PANEL[i]) for i in eachindex(sh) if sh[i] > 3n[i]]
    @printf("| %s → %s | %.3g | %.3g | %d%s | %.3g (%.3g) |\n", a, b, sum(sh), sum(n), length(bad),
            isempty(bad) ? "" : " (" * join(replace.(bad, "M_" => ""), ", ") * ")", sh[gtp], n[gtp])
end
for b in keys(data)
    sh = vec(sum((data[b].M .- data[b].L) .^ 2; dims = 2)) ./ (2σ^2)
    @printf("uncorrected %s: total %.3g, GTP %.3g\n", b, sum(sh), sh[gtp])
end
