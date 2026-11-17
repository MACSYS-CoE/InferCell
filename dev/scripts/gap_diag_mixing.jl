# KIM-90 item 2: is the m0_v6 gap larger than either side's Monte Carlo error?
#
# Reads only archives: V6's four exact chains (F10's TSV, sweeps 101–1,650) and
# m0_v6's two mixed chains (B7's archive, 200 draws each after 200 adaptation).
# Every distance is in V6 posterior SDs, as AV6's was.
#
# Per parameter:
#   gap        W1(mixed, V6), and the mean shift, as AV6 reported them;
#   V6 noise   split-half W1 (sweeps 101–875 against 876–1,650, chains pooled);
#              chain-vs-rest W1 at full length, and at F10's matched length
#              (one chain's sweeps 101–656 against the other three, AV6's nulls);
#              MCSE of the mean (bulk ESS);
#   mixed noise chain 1 against chain 2, split-half within the pooled 400, MCSE
#              of the mean;
#   z          the mean shift over the two MCSEs combined in quadrature.
#
# Usage: julia --project dev/scripts/gap_diag_mixing.jl <out.md>

using InferCell
using Serialization
using Statistics
using Printf
using MCMCChains
using MCMCChains: MCMCDiagnosticTools
const MDT = MCMCDiagnosticTools

const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
const B7_RUN = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6"
const PNAMES = ["S_JCVISYN3A_0607", "S_JCVISYN3A_0779", "krnadeg", "kcatF_R_ENO",
                "kcatF_R_FBA", "sigma_b"]
const V6_FIRST, V6_SPLIT, V6_LAST, F10_LEN = 101, 875, 1650, 656

# AV6's W1 (dev/scripts/av6_89_b7.jl), unchanged.
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

# MCSE of the mean from bulk ESS; x is draws × chains (MDT wants draws × chains × params).
mcse_mean(x::AbstractMatrix) = std(vec(x)) / sqrt(only(MDT.ess(reshape(x, size(x)..., 1);
                                                             kind = :bulk)))

function v6_draws()
    lines = readlines(F10_DRAWS)
    head = split(lines[1], '\t')
    rows = [parse.(Float64, split(l, '\t')) for l in lines[2:end]]
    cols = [findfirst(==(n), head) for n in PNAMES]
    any(isnothing, cols) && error("F10's TSV lacks a column: $head")
    # chain × sweep × parameter, sweeps V6_FIRST:V6_LAST
    out = fill(NaN, 4, V6_LAST - V6_FIRST + 1, length(PNAMES))
    for r in rows
        c, sw = Int(r[1]), Int(r[2])
        (c in 1:4 && V6_FIRST <= sw <= V6_LAST) || continue
        out[c, sw - V6_FIRST + 1, :] = r[cols]
    end
    any(isnan, out) && error("V6's draws are incomplete")
    return out
end

function mixed_draws()
    chs = [deserialize(joinpath(B7_RUN, "chain_$c.jls")) for c in 1:2]
    # draws × parameter × chain
    return cat((Array(c.chain[:, [Symbol("x[$i]") for i in 1:6], 1]) for c in chs)...; dims = 3),
           [(c.seed, c.commit, c.job) for c in chs]
end

function main(out)
    v6 = v6_draws()
    mx, meta = mixed_draws()
    io = IOBuffer()
    say(s) = (println(s); println(io, s))
    say("# KIM-90 item 2: Monte Carlo error on both sides of the m0_v6 gap\n")
    say("V6: 4 chains × sweeps $(V6_FIRST)–$(V6_LAST) from `f10_16a10_draws.tsv` " *
        "(V6 jobs 18276511–14 and 18420986–89). m0_v6: 2 chains × 200 draws from " *
        "`$B7_RUN` (chains $(join(["seed $(m[1]), job $(m[3]), commit $(m[2])" for m in meta], "; "))).")
    say("Units: V6 posterior SDs. Gap = W1(mixed, V6); shift = mean(mixed) − mean(V6).\n")
    say("| Parameter (log) | gap W1 | shift | V6 split-half W1 | V6 chain-vs-rest W1 (full) | " *
        "AV6 matched nulls (656) | V6 MCSE(mean) | mixed chain 1 vs 2 W1 | mixed split-half W1 | " *
        "mixed MCSE(mean) | z (shift / MCSE) |")
    say("|---|---|---|---|---|---|---|---|---|---|---|")
    nsplit = V6_SPLIT - V6_FIRST + 1
    for (p, n) in enumerate(PNAMES)
        ref = vec(v6[:, :, p])
        sd = std(ref)
        mix = vec(mx[:, p, :])
        gap = w1(mix, ref) / sd
        shift = (mean(mix) - mean(ref)) / sd
        split = w1(vec(v6[:, 1:nsplit, p]), vec(v6[:, nsplit+1:end, p])) / sd
        rest(k, len) = w1(vec(v6[k, 1:len, p]), vec(v6[setdiff(1:4, k), :, p])) / sd
        cvr = [rest(k, size(v6, 2)) for k in 1:4]
        nulls = [rest(k, F10_LEN - V6_FIRST + 1) for k in 1:4]
        v6mcse = mcse_mean(permutedims(v6[:, :, p])) / sd
        m12 = w1(mx[:, p, 1], mx[:, p, 2]) / sd
        h = size(mx, 1) ÷ 2
        msplit = w1(vec(mx[1:h, p, :]), vec(mx[h+1:end, p, :])) / sd
        mxmcse = mcse_mean(mx[:, p, :]) / sd
        z = shift / hypot(v6mcse, mxmcse)
        f(v) = join((@sprintf("%.3f", x) for x in v), ", ")
        say(@sprintf("| log_%s | %.3f | %+.3f | %.3f | %s | %s | %.3f | %.3f | %.3f | %.3f | %+.1f |",
                     n, gap, shift, split, f(cvr), f(nulls), v6mcse, m12, msplit, mxmcse, z))
    end
    say("\nThe gap is beyond Monte Carlo error on a parameter when it exceeds every V6 and " *
        "mixed noise column and |z| is large; z counts only the error of the means.")
    write(out, take!(io))
end

abspath(PROGRAM_FILE) == (@__FILE__) && main(ARGS[1])
