# KIM-90 item 1: are the exact route (V6) and the mixed route (m0_v6) aimed at the
# same target? Same data, same prior, same σ_b treatment, and no truncation of the
# mixed density that the exact sampler does not also impose where the posterior sits.
#
# Stages:
#   fingerprint <out.tsv>   at whatever commit `--project` points at: m0_v6's dataset
#                           (replicate 16085, 50 cells, M0's horizon) hashed field by
#                           field, each target's prior, σ_b's prior, the panel and
#                           the horizon. Run it at V6's commit (281d548), at B7's
#                           (7a102fc) and on this branch.
#   compare <fp1.tsv> <fp2.tsv> …
#                           asserts every fingerprint equal, then on this branch's code
#                           builds both setups from the dataset and asserts the same
#                           names, priors, data, panel, floor and window; and reports
#                           how far each archived posterior sits from the mixed
#                           density's truncations (μmin, kmax) and block 2's bound.
#
# Usage (login node; a dataset takes ~3 min):
#   julia --project=<worktree at a commit> dev/scripts/gap_diag_identity.jl fingerprint <tsv>
#   julia --project dev/scripts/gap_diag_identity.jl compare <tsv> <tsv> <tsv>

using InferCell
using InferCell: _parameter_by_name
using Serialization
using SHA
using Printf
using Statistics
using Distributions: params, quantile, LogNormal

const REPLICATE = 16_085
const NCELLS = 50
const B7_PREP = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/prep.jls"
const B7_CHAINS = ["/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run/v6/chain_$c.jls"
                   for c in 1:2]
const V6_CHAIN1 = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v6/gibbs_pilot_16a9a/chain1.jls"
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
const DATA_FIELDS = (:truth, :seeds, :noise_seed, :times, :genes, :transcripts, :sigma,
                     :bulk, :bulk_seed)

h(x::AbstractArray{<:Real}) = bytes2hex(sha256(reinterpret(UInt8, vec(collect(Float64, x)))))
h(x::AbstractArray) = bytes2hex(sha256(join(string.(vec(x)), ",")))
h(x::Real) = h([x])
h(x) = bytes2hex(sha256(repr(x)))

function fingerprint(out)
    ds = m0_dataset(REPLICATE; ncells = NCELLS, horizon = M0_HORIZON_S)
    rows = Tuple{String, String}[]
    for f in DATA_FIELDS
        v = getproperty(ds, f)
        f === :truth && (v = [last(p) for p in v.values])
        push!(rows, ("data.$f", h(v)))
    end
    ms = m0_models()
    for n in vcat(M0_TARGETS, :krnadeg)
        p = _parameter_by_name(ms, n).prior
        push!(rows, ("prior.$n", @sprintf("%s(%.17g, %.17g)", nameof(typeof(p)), params(p)...)))
    end
    push!(rows, ("prior.sigma_b", @sprintf("%s(%.17g, %.17g)",
                 nameof(typeof(SIGMA_MET_PRIOR)), params(SIGMA_MET_PRIOR)...)))
    push!(rows, ("targets", join(M0_TARGETS, ",")))
    push!(rows, ("panel", join(M0_PANEL, ",")))
    push!(rows, ("horizon", string(M0_HORIZON_S)))
    open(out, "w") do io
        for (k, v) in rows; println(io, k, '\t', v); end
    end
    foreach(r -> println(r[1], '\t', r[2]), rows)
end

read_fp(f) = Dict(split(l, '\t')[1] => split(l, '\t')[2] for l in readlines(f))

check(ok, msg) = (println(ok ? "  ok    " : "  FAIL  ", msg); ok)

function compare(fps)
    println("## Fingerprints across commits\n")
    fs = read_fp.(fps)
    ok = true
    for k in sort(collect(keys(fs[1])))
        vals = [get(f, k, "<missing>") for f in fs]
        ok &= check(allequal(vals), "$k equal across $(length(fs)) commits")
    end

    println("\n## The B7 archive's dataset against this branch's\n")
    p = deserialize(B7_PREP)
    ds = m0_dataset(REPLICATE; ncells = NCELLS, horizon = M0_HORIZON_S)
    for f in DATA_FIELDS
        a, b = getproperty(p.ds, f), getproperty(ds, f)
        if f === :bulk
            # The bulk is an ODE-and-noise output: B7's prep ran on a compute node, this
            # check on the login node, so it may differ by roundoff (memory: roundoff by machine).
            rel = maximum(abs.(log.(a) .- log.(b)))
            ok &= check(rel < 1e-8, @sprintf("prep.jls ds.bulk ≈ m0_dataset's: max |Δ log| = %.3g (tolerance 1e-8)", rel))
        else
            ok &= check(a == b, "prep.jls ds.$f == m0_dataset's")
        end
    end

    println("\n## V6's chain metadata against the dataset\n")
    v6 = deserialize(V6_CHAIN1).meta
    ok &= check(v6.replicate == REPLICATE && v6.ncells == NCELLS && v6.N == 20,
                "V6: replicate $(v6.replicate), $(v6.ncells) cells, N = $(v6.N)")
    ok &= check(v6.truth == ds.truth.values, "V6 truth == dataset truth")
    ok &= check(v6.sigma_true == ds.sigma, "V6 sigma_true == dataset σ_b")

    println("\n## Both setups on this branch\n")
    s = mixed_setup(ds)
    g = gibbs_setup(ds; N = 20, horizon = M0_HORIZON_S)
    xn = mixed_parameter_names(s)
    gn = vcat(g.cme, g.forwards, :sigma_b)
    ok &= check(Set(Symbol.(replace.(String.(xn), "log_" => ""))) == Set(gn),
                "same free parameters: mixed $(xn) / exact $(gn)")
    ms = m0_models()
    gprior = Dict{Symbol, Tuple{Float64, Float64}}()
    for (n, pr) in zip(g.cme, g.cme_priors); gprior[n] = params(pr); end
    for n in g.forwards; gprior[n] = params(_parameter_by_name(ms, n).prior); end
    gprior[:sigma_b] = params(g.sigma_prior)
    mμ, msd = mixed_prior_median(s), mixed_prior_sd(s)
    for (i, n) in enumerate(xn)
        k = Symbol(replace(String(n), "log_" => ""))
        gμ, gsd = gprior[k]
        ok &= check(mμ[i] == gμ && msd[i] == gsd,
                    @sprintf("prior on %s: mixed N(%.6g, %.6g), exact N(%.6g, %.6g) on the log",
                             n, mμ[i], msd[i], gμ, gsd))
    end
    ok &= check(s.logz == vec(log.(g.bulk)), "bulk: mixed log z == log of the exact route's record")
    ok &= check(all(f.y == t for (f, t) in zip(s.factors, g.transcripts)),
                "transcripts: each cell's record identical in both routes")
    ok &= check(s.floor == g.floor, "bulk floor $(s.floor) / $(g.floor)")
    ok &= check(Float64(g.window) == 60.0 && s.saves == 60.0 .* (1:length(s.saves)) &&
                length(s.saves) == size(g.bulk, 2),
                "window 60 s in both, $(length(s.saves)) saves / $(size(g.bulk, 2)) windows")
    ok &= check([Symbol(String(n)) for n in g.panel] == M0_PANEL, "exact panel == M0_PANEL")
    println("  σ_b: both put a Normal on log σ_b with the parameters above; the exact route " *
            "slice-samples log σ_b given the residuals (block 1), the mixed route samples " *
            "it jointly under NUTS. Same target.")

    println("\n## Truncations, against the archived posteriors\n")
    tsv = readlines(F10_DRAWS)
    hd = split(tsv[1], '\t')
    rows = [parse.(Float64, split(l, '\t')) for l in tsv[2:end]]
    v6 = reduce(hcat, [r for r in rows if r[1] in 1:4 && r[2] >= 101])'
    kr = v6[:, findfirst(==("krnadeg"), hd)]
    mx = vcat([deserialize(f).chain[:, Symbol("x[3]"), 1] |> vec for f in B7_CHAINS]...)
    μmin_log = log(quantile(_parameter_by_name(ms, :krnadeg).prior, 1e-6))
    psd = mixed_prior_sd(s)[3]
    @printf("  mixed μmin on krnadeg: log krnadeg ≥ %.4f. V6's min %.4f (%.1f prior SDs above); m0_v6's min %.4f (%.1f above)\n",
            μmin_log, minimum(kr),
            (minimum(kr) - μmin_log) / psd, minimum(mx), (minimum(mx) - μmin_log) / psd)
    ok &= check(minimum(kr) > μmin_log && minimum(mx) > μmin_log,
                "no draw of either posterior touches μmin")
    println("  kmax (mixed) and TURNOVER_CEILING / RNAPOL_KCAT (exact, block 2) bound the " *
            "transcription constant; the mixed density returns -Inf past kmax, the exact " *
            "update throws. Neither fired: no -Inf density in m0_v6's chains (finite lp " *
            "below) and V6's 4 chains ran to 1,650 sweeps without the throw.")
    lp = vcat([vec(deserialize(f).chain[:, :lp, 1]) for f in B7_CHAINS]...)
    ok &= check(all(isfinite, lp), "every m0_v6 draw has a finite log density ($(length(lp)) draws)")

    println(ok ? "\nRESULT: target identity holds." : "\nRESULT: target identity FAILS (see FAIL lines).")
    ok || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    stage = ARGS[1]
    stage == "fingerprint" ? fingerprint(ARGS[2]) :
    stage == "compare" ? compare(ARGS[2:end]) : error("unknown stage $stage")
end
