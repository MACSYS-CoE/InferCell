# KIM-18: confirm the archived V7 (exact) and B7 (mixed) posteriors hold draws,
# that their coordinates line up with M0_TARGETS plus σ_b, and that each dataset's
# true θ can be recovered for scoring (completion plan §3.5).
#
# Reads only; writes nothing but stdout. Run on main, which is what `predict`
# (KIM-28) will load the files with.
# Usage: julia --project dev/scripts/kim18_inspect_archives.jl

using InferCell
using Random
using Serialization
using Statistics
using Printf
using MCMCChains

const V7_DIR = "/fred/oz022/tkimpson/infercell-16a9b-archive-2026-10-17/v7/sbc_v7_16a9b"
const B7_DIR = "/fred/oz022/tkimpson/infercell-89-b7-archive-2026-10-27/run"
const V7_BURN, V7_SWEEPS = 100, 656
const DATA_SEED0 = 16_200_000
const TARGET_NAMES = vcat([Symbol(:log_, n) for n in M0_TARGETS], :log_sigma_b)

# The truth's leading draw, as m0_dataset makes it, without simulating cells.
truth_targets(seed) = Dict(draw_truth(Xoshiro(seed), m0_models(); names = M0_TARGETS,
                                      purpose = :calibration).values)

function inspect_v7()
    println("## V7 (exact, Gibbs): $V7_DIR\n")
    files = sort(filter(f -> occursin(r"^replicate\d+\.jls$", f), readdir(V7_DIR));
                 by = f -> parse(Int, match(r"\d+", f).match))
    println("files: $(length(files))")
    x1 = deserialize(joinpath(V7_DIR, files[1]))
    println("top-level fields: ", keys(x1))
    println("state fields: ", keys(x1.state))
    println("meta fields: ", keys(x1.meta))
    println("trace element type: ", typeof(x1.state.trace[1]))
    println("trace element fields: ", propertynames(x1.state.trace[1]))
    truths = Dict{Int,Vector{Float64}}()
    worst = 0.0
    for f in files
        x = deserialize(joinpath(V7_DIR, f))
        r = x.meta.replicate
        vnames = [n === :sigma_b ? :log_sigma_b : Symbol(:log_, n) for n in x.meta.names]
        sort(vnames) == sort(TARGET_NAMES) || error("$f: names $(vnames)")
        kept = x.state.trace[(V7_BURN + 1):V7_SWEEPS]
        X = reduce(hcat, [vcat(log.(t.cme), t.u, log(t.σ)) for t in kept])'
        all(isfinite, X) || error("$f: non-finite draws")
        # The truth in the file against the seed's.
        tv = truth_targets(x.meta.seed)
        for (i, n) in enumerate(x.meta.names)
            n === :sigma_b && continue
            worst = max(worst, abs(x.meta.truth[i] - log(tv[n])))
        end
        truths[r] = x.meta.truth[[findfirst(==(n), vnames) for n in TARGET_NAMES]]
        @printf("replicate %2d: seed %d, %d sweeps, kept %d × %d, sd(log) = %s\n", r, x.meta.seed,
                length(x.state.trace), size(X)..., round.(vec(std(X; dims = 1)); sigdigits = 2))
    end
    println("\nname order in the file: ", x1.meta.names)
    @printf("largest |file truth − seed truth| over targets (log): %.3g\n", worst)
    return truths
end

function inspect_b7(v7truths)
    println("\n## B7 (mixed, NUTS): $B7_DIR\n")
    dirs = sort(filter(d -> isdir(joinpath(B7_DIR, d)) &&
                            isfile(joinpath(B7_DIR, d, "chain_1.jls")), readdir(B7_DIR)))
    nfiles = sum(count(f -> occursin(r"^chain_\d+\.jls$", f), readdir(joinpath(B7_DIR, d)))
                 for d in dirs)
    println("fit directories with chains: $(length(dirs)); chain files: $nfiles")
    c1 = deserialize(joinpath(B7_DIR, "r1", "chain_1.jls"))
    println("chain file fields: ", keys(c1))
    println("chain: ", size(c1.chain), " (iterations, parameters, chains)")
    println("chain names: ", names(c1.chain))
    p1 = deserialize(joinpath(B7_DIR, "r1", "prep.jls"))
    println("prep fields: ", keys(p1))
    pnames = mixed_parameter_names(mixed_setup(p1.ds))
    println("x[i] ↔ ", pnames)
    pnames == TARGET_NAMES || println("NOTE: order differs from M0_TARGETS + σ_b")
    worst = 0.0
    for d in dirs
        for c in 1:2
            ch = deserialize(joinpath(B7_DIR, d, "chain_$c.jls"))
            P = count(nm -> startswith(String(nm), "x["), names(ch.chain))
            X = Array(ch.chain[:, [Symbol("x[$i]") for i in 1:P], 1])
            all(isfinite, X) || error("$d/chain_$c: non-finite draws")
            @printf("%-9s chain %d: %d × %d draws, %d + %d steps, commit %s, sd = %s\n", d, c,
                    size(X)..., ch.nadapt, ch.ndraw, ch.commit,
                    round.(vec(std(X; dims = 1)); sigdigits = 2))
        end
        isfile(joinpath(B7_DIR, d, "prep.jls")) || continue
        p = deserialize(joinpath(B7_DIR, d, "prep.jls"))
        hasproperty(p, :ds) || continue
        tv = Dict(p.ds.truth.values)
        t = vcat([log(tv[n]) for n in M0_TARGETS], log(p.ds.sigma))     # TARGET_NAMES' order
        m = match(r"^r(\d+)$", d)
        if m !== nothing
            worst = max(worst, maximum(abs.(t .- v7truths[parse(Int, m[1])])))
        end
    end
    @printf("\nlargest |B7 prep truth − V7 file truth| over r1–r50 (log, incl. σ_b): %.3g\n", worst)
end

function main()
    t = inspect_v7()
    inspect_b7(t)
end

main()
