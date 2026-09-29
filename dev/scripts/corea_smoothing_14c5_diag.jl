# Spec §11 task 14c.5: why the coupled half fails (§12 2026-09-29, the
# coupled-half entry). Seed 18227 carries the 5,000-seed run's largest coupled
# difference, 3.75% on `M_gtp_c` at t = 4500 s (job 17671709). This reruns that
# pair for 4,500 handshakes and prints the GTP pool under both models at each
# 60 s save point where they differ by more than two particles, with each
# model's GTP_translat deficit and its clipped drains in the preceding minute.
#
# Run of record: job 17679585, at commit 87b88f1; output in
# corea_smoothing_14c5_diag_result.md.
# Usage: julia --project dev/scripts/corea_smoothing_14c5_diag.jl

using InferCell, Random, Printf
seed = 18227; T = 4500
function run(v)
    Random.seed!(seed)
    kw = v === :sampler ? COREA_SAMPLER : (;)
    ms = corea_models(; smoothing = get(kw, :smoothing, nothing))
    d = build_corea(; tspan = (0.0, 6300.0), kw...)
    names = InferCell._block_names(ms, :ode); ig = findfirst(==(:M_gtp_c), names)
    kg = only(k for (k, b) in enumerate(d.debits) if b.counter === :GTP_translat && b.sign < 0)
    Random.seed!(seed)
    gtp = zeros(T); clip = falses(T); def = zeros(T); jumps = Vector{Vector{Int}}(undef, T)
    for k in 1:T
        handshake_step!(d)
        gtp[k] = d.ode.u[ig] * d.factor + d.rounding.remainders[ig]
        clip[k] = d.debits[kg].clipped; def[k] = d.debits[kg].deficit
        jumps[k] = collect(Int, d.jump.u)
    end
    return (; gtp, clip, def, jumps)
end
p = run(:published); s = run(:sampler)
dec = findfirst(k -> p.jumps[k] != s.jumps[k], 1:T)
println("first jump mismatch: ", dec === nothing ? "none by $T" : dec)
println("published clipped drains: ", count(p.clip), "; sampler: ", count(s.clip))
diff = s.gtp .- p.gtp
for k in (60:60:T)
    abs(diff[k]) > 2 || continue
    @printf("t=%4d  pub=%9.2f  smp=%9.2f  diff=%7.2f  rel=%6.2f%%  pubdef=%8.1f smpdef=%8.1f  clips in last 60 s: pub %2d smp %2d\n",
            k, p.gtp[k], s.gtp[k], diff[k], 100diff[k] / max(abs(p.gtp[k]), 500), p.def[k], s.def[k],
            count(p.clip[k-59:k]), count(s.clip[k-59:k]))
end
