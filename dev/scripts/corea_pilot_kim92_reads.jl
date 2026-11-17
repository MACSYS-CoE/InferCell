# KIM-92: the result note's derived numbers, read from a finished pilot outdir.
# Read-only and light (no model evaluation), so it runs on a login node:
#   julia --project dev/scripts/corea_pilot_kim92_reads.jl <outdir>
# It reads chain_c.jls, extend_c.jls, metric.jls, prep.jls and pilot_coordinates.tsv, as
# written by dev/scripts/corea_pilot_kim92.jl at de24b27, and prints:
#   1. per chain, the kept draws' mean and SD, and each half's mean, of the coordinates
#      that failed the 200-draw gate (log_S_0213 and log_kcatF_R_ENO);
#   2. per chain and phase (adaptation 1–200, kept 201–400, extension), the mean leapfrog
#      steps per iteration, numerical errors and tree-depth counts, plus the step size;
#   3. the conditioning number of −H under the unit metric, κ(−H), and under the warm
#      diagonal metric, κ(D½ (−H) D½) with D = diag(minv);
#   4. the Laplace SD at θ_ref, √diag((−H)⁻¹), against the posterior SD of the 1,600
#      pooled kept draws (4 × (200 + 200));
#   5. the standardised errors, (posterior mean − truth)/sd, over the 17 promoters, and
#      θ_ref (the MAP) against the truth and the posterior mean, in posterior SDs;
#   6. the gradient's cost: each chain's sampling seconds over its total leapfrog steps.
# Wall times of single iterations are not stored in the chain files. The note takes them
# from the progress lines in logs/kim92_chain<c>_<job>.log, one every 5 iterations.
using InferCell, Turing, MCMCChains, Serialization, Printf, Statistics, LinearAlgebra
using DelimitedFiles

const OUT = abspath(ARGS[1])
const NADAPT = 200
ld(name) = deserialize(joinpath(OUT, "$name.jls"))

counts(v) = join(["$k=>$(count(==(k), v))" for k in sort(unique(v))], ", ")

function phase(label, ch)
    d, n = Int.(vec(ch[:tree_depth])), Int.(vec(ch[:n_steps]))
    e = Int.(vec(ch[:numerical_error]))
    @printf("  %-12s %3d iterations, mean steps %5.1f, total %5d, numerical errors %2d, depths %s\n",
            label, length(d), mean(n), sum(n), sum(e), counts(d))
    return sum(n)
end

function main()
    A = [ld("chain_$c") for c in 1:4]
    E = [ld("extend_$c") for c in 1:4]
    pn = names(A[1].result.chain, :parameters)
    T, _ = readdlm(joinpath(OUT, "pilot_coordinates.tsv"), '\t'; header = true)
    name = String.(T[:, 1])
    pm, ps, tr = Float64.(T[:, 6]), Float64.(T[:, 7]), Float64.(T[:, 8])
    idx(s) = findfirst(==(s), name)

    println("## 1. Coordinates that failed at 200 draws (kept draws 201–400)")
    for s in ("log_S_JCVISYN3A_0213", "log_kcatF_R_ENO")
        i = idx(s)
        println(s)
        for c in 1:4
            v = vec(Array(A[c].result.chain[(NADAPT + 1):end, pn[i:i], :]))
            h = length(v) ÷ 2
            @printf("  chain %d: mean %.4f, sd %.4f, halves %.4f / %.4f\n", c, mean(v), std(v),
                    mean(v[1:h]), mean(v[(h + 1):end]))
        end
    end

    println("\n## 2. Sampler, by chain and phase")
    for c in 1:4
        ch, ex = A[c].result.chain, E[c].result.chain
        @printf("chain %d (job %s; extension job %s), step size kept %.4f, extension %.4f\n",
                c, string(A[c].job), string(E[c].job), ch[:step_size][end], ex[:step_size][end])
        na = phase("adaptation", ch[1:NADAPT])
        nk = phase("kept", ch[(NADAPT + 1):end])
        nx = phase("extension", ex)
        @printf("  6. seconds per leapfrog step: chain %.1f, extension %.1f\n",
                A[c].result.seconds / (na + nk), E[c].result.seconds / nx)
    end

    println("\n## 3. Conditioning")
    m = ld("metric").result
    negH = Symmetric(-m.H)
    ev = eigvals(negH)
    D = Diagonal(sqrt.(m.minv))
    ew = eigvals(Symmetric(D * negH * D))
    @printf("unit metric: κ(−H) = %.3g (√κ %.3g); warm diagonal: κ = %.3g (√κ %.3g)\n",
            ev[end] / ev[1], sqrt(ev[end] / ev[1]), ew[end] / ew[1], sqrt(ew[end] / ew[1]))

    println("\n## 4. Laplace SD at θ_ref against the posterior SD")
    X = vcat([vcat(Array(A[c].result.chain[(NADAPT + 1):end, pn, :]),
                   Array(E[c].result.chain[:, pn, :])) for c in 1:4]...)
    s = vec(std(X; dims = 1))
    lap = sqrt.(diag(inv(negH)))
    r = lap ./ s
    @printf("ratio range %.3f to %.3f over %d pooled draws\n", minimum(r), maximum(r), size(X, 1))
    for i in eachindex(name)
        @printf("  %-24s Laplace %.4f, posterior %.4f, ratio %.3f\n", name[i], lap[i], s[i], r[i])
    end

    println("\n## 5. Standardised errors and θ_ref")
    z = (pm .- tr) ./ ps
    prom = findall(startswith.(name, "log_S_"))
    @printf("promoters: %d, below the truth %d, mean z %.3f, sd %.3f, SE if independent %.3f\n",
            length(prom), count(<(0), z[prom]), mean(z[prom]), std(z[prom]), 1 / sqrt(length(prom)))
    @printf("all %d coordinates: Σz² = %.2f\n", length(z), sum(abs2, z))
    x = ld("prep").result.ref.x
    zmap, dmap = (x .- tr) ./ ps, (pm .- x) ./ ps
    @printf("θ_ref: mean (θ_ref − truth)/sd over the promoters %.3f\n", mean(zmap[prom]))
    for i in eachindex(name)
        @printf("  %-24s θ_ref %.4f, mean %.4f, truth %.4f, (θ_ref − truth)/sd %+.2f, (mean − θ_ref)/sd %+.2f\n",
                name[i], x[i], pm[i], tr[i], zmap[i], dmap[i])
    end
end

main()
