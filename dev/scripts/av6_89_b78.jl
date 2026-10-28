# spec/phases/approximate-method.md B7.8, under §12 "B7.8 as planned": the mixed
# route's stream arms (B7's tx_v6 and bulk_v6) against 16b.4's exact single-stream
# fits on m0_v6. Per marginal: W1 in exact posterior SD, beside a floor (the larger
# of the exact chain's split-half W1 and the mixed fit's chain-vs-chain W1).
#
#   julia --project dev/scripts/av6_89_b78.jl <b7_run_dir> <exact_tx_dir> <exact_met_dir>
#
# <b7_run_dir> holds av6_{tx,bulk}_v6_draws.tsv (B7's merge, job 19292232);
# the exact dirs hold 16b's chain.jls (jobs 19034902 and 19034903).

using InferCell
using Serialization
using Statistics
using Printf
using Distributions: params
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

const BURN = 100
const NMIX = 200            # draws per mixed chain; the TSV holds chain 1, then chain 2

# W1 between two samples: ∫ |F_a − F_b| (av6_89_b7.jl's estimator).
function w1(a::AbstractVector, b::AbstractVector)
    xs = sort(vcat(a, b))
    sa, sb = sort(a), sort(b)
    s = 0.0
    for i in 1:(length(xs) - 1)
        Fa = searchsortedlast(sa, xs[i]) / length(sa)
        Fb = searchsortedlast(sb, xs[i]) / length(sb)
        s += abs(Fa - Fb) * (xs[i + 1] - xs[i])
    end
    return s
end

function read_tsv(path)
    ls = readlines(path)
    names = replace.(split(ls[1], '\t'), "log_" => "")
    return String.(names), reduce(vcat, [permutedims(parse.(Float64, split(l, '\t'))) for l in ls[2:end]])
end

# The exact chain's kept draws (sweeps BURN+1 on), sweeps × parameters, log scale.
function exact(dir)
    r = deserialize(joinpath(dir, "chain.jls"))
    x = reduce(vcat, [permutedims(vcat(log.(t.cme), t.u, log(t.σ))) for t in r.state.trace[(BURN + 1):end]])
    return string.(r.meta.names), x
end

function main(b7dir, txdir, metdir)
    ms = m0_models()
    for (lab, arm, dir) in (("transcripts only", "tx_v6", txdir), ("metabolites only", "bulk_v6", metdir))
        mnames, M = read_tsv(joinpath(b7dir, "av6_$(arm)_draws.tsv"))
        enames, E = exact(dir)
        size(M, 1) == 2NMIX || error("$arm: expected $(2NMIX) draws, got $(size(M, 1))")
        x3 = reshape(E, size(E, 1), 1, size(E, 2))
        ess, rhat = MDT.ess(x3; kind = :bulk), MDT.rhat(x3)
        h = size(E, 1) ÷ 2
        println("\n### $lab: $arm (mixed, 2 × $NMIX draws) against the exact fit ($(size(E, 1)) draws)\n")
        println("| Parameter | Exact mean ± SD | Mixed mean ± SD | W1 / exact SD | Exact split-half | Mixed chain-vs-chain | Floor | Exact ESS | Exact split-R̂ | Reads |")
        println("|---|---|---|---|---|---|---|---|---|---|")
        for (q, n) in enumerate(enames)
            m = M[:, findfirst(==(n), mnames)]
            e = E[:, q]
            if lab == "transcripts only" && n == "sigma_b"
                @printf("| %s | held at its start | %.3f ± %.3f | — | — | — | — | — | — | not compared |\n", n, mean(m), std(m))
                continue
            end
            if rhat[q] > 1.05
                @printf("| %s | unmixed | %.3f ± %.3f | — | — | — | — | %.0f | %.3f | not read |\n", n, mean(m), std(m), ess[q], rhat[q])
                continue
            end
            sd = std(e)
            d = w1(m, e) / sd
            sh = w1(e[1:h], e[(h + 1):(2h)]) / sd
            cc = w1(m[1:NMIX], m[(NMIX + 1):end]) / sd
            fl = max(sh, cc)
            reads = d > fl ? @sprintf("resolved, by %.3f", d - fl) : "not resolved"
            ess[q] < 100 && (reads *= ", exact ESS below 100")
            @printf("| %s | %.3f ± %.3f | %.3f ± %.3f | %.3f | %.3f | %.3f | %.3f | %.0f | %.3f | %s |\n",
                    n, mean(e), sd, mean(m), std(m), d, sh, cc, fl, ess[q], rhat[q], reads)
        end
    end
    # ptsG under metabolites only: the mixed fit's shrinkage, with no exact counterpart.
    mnames, M = read_tsv(joinpath(b7dir, "av6_bulk_v6_draws.tsv"))
    m = M[:, findfirst(==("S_JCVISYN3A_0779"), mnames)]
    psd = params(InferCell._parameter_by_name(ms, :S_JCVISYN3A_0779).prior)[2]
    @printf("\nptsG under metabolites only, mixed (bulk_v6): SD / prior SD %.3f, 1 − var_post / var_prior %.3f. The exact fit is unmixed and not read.\n",
            std(m) / psd, 1 - var(m) / psd^2)
end

main(ARGS[1], ARGS[2], ARGS[3])
