# Spec §11 task 14c.5: which of the sampler model's two departures causes the
# coupled half's GTP gap (§12 2026-09-29, the coupled-half entry). A 2×2 on one
# seed, every handshake to the pair's decoupling: A the published model (clamp,
# fractional carry), B the clamp with continuous pools, C the smoothed drain
# with fractional carry, D the sampler model. Records the GTP, GDP, GMP, ATP,
# ADP and Pi pools in particles, each GTP debit's deficit and clip flag, and the
# jump state. Analyse with corea_smoothing_14c5_decomp_merge.jl.
#
# Run of record: job 17695360, at commit 2bbb02c, seeds 18227, 16184, 19958 and
# 16391 (the 5,000-seed run's four worst coupled differences) to handshakes
# 4515, 6110, 5477 and 3315. Output in corea_smoothing_14c5_decomp_result.md.
# Usage: julia --project dev/scripts/corea_smoothing_14c5_decomp.jl <seed> <handshakes> <out.jls>
# Regenerate: sbatch dev/scripts/corea_smoothing_14c5_decomp.slurm

using InferCell, Random, Serialization
seed = parse(Int, ARGS[1]); T = parse(Int, ARGS[2]); OUT = ARGS[3]
variants = (A = (;), B = (; rounding = :continuous), C = (; smoothing = 1.0),
            D = (; smoothing = 1.0, rounding = :continuous))
res = Dict{Symbol, Any}()
for (name, kw) in pairs(variants)
    Random.seed!(seed)
    ms = corea_models(; smoothing = get(kw, :smoothing, nothing))
    d = build_corea(; tspan = (0.0, 6300.0), kw...)
    names = InferCell._block_names(ms, :ode)
    ix = Dict(s => findfirst(==(s), names) for s in (:M_gtp_c, :M_gdp_c, :M_gmp_c, :M_atp_c, :M_adp_c, :M_pi_c))
    gtpdeb = [k for (k, b) in enumerate(d.debits) if b.species === :M_gtp_c && b.sign < 0]
    labels = [d.debits[k].counter for k in gtpdeb]
    Random.seed!(seed)
    st = Dict(s => zeros(T) for s in keys(ix))
    def = zeros(T, length(gtpdeb)); clip = falses(T, length(gtpdeb))
    jumps = Vector{Vector{Int}}(undef, T)
    for k in 1:T
        handshake_step!(d)
        for (s, i) in ix
            st[s][k] = d.ode.u[i] * d.factor + d.rounding.remainders[i]
        end
        for (j, kk) in enumerate(gtpdeb)
            def[k, j] = d.debits[kk].deficit; clip[k, j] = d.debits[kk].clipped
        end
        jumps[k] = collect(Int, d.jump.u)
    end
    res[name] = (; st, def, clip, labels, jumps)
    println("done $name"); flush(stdout)
end
serialize(OUT, (seed = seed, T = T, res = res))
