# Tier B: refit under Block-1 metabolism mismatch.
#
# Generate synthetic data from `TierBMetabolism` (Hill response + ATP leak) and
# refit using the v0.0.1 `LightMetabolism` (Michaelis-Menten, no leak). The
# scientific claim Tier B exists to surface: shared rate parameters that have
# meaning in both models should still fall in the fitted CI; mismatched ones
# (the rates absorbing the missing physics) should drift in a documented
# direction.

using InferCell
using OrdinaryDiffEq: Tsit5, solve
using Random
using Statistics: quantile, mean
using Dates

println("=== Tier B: Block-1 metabolism mismatch refit ===")
println("Start: $(Dates.now())")
flush(stdout)

Random.seed!(2026)

# --- Tier B generator ---
println("\n1. Generating Tier B synthetic data...")
flush(stdout)
txl_gen = TranscriptionTranslation()
gen = TierBMetabolism()  # n_hill=2, k_atp_leak=0.2
prob_gen = build_problem([txl_gen, gen]; tspan=(0.0, 50.0))
sol_gen = solve(prob_gen, Tsit5())
data = observe(sol_gen, 0.0:2.5:50.0, [txl_gen, gen]; sigma=0.3)
println("   data shape: $(size(data.observations))")
flush(stdout)

# --- Refit with the v0.0.1 model (deliberate mismatch) ---
println("\n2. Refitting with v0.0.1 LightMetabolism...")
flush(stdout)
txl_fit = TranscriptionTranslation()
fit = LightMetabolism()
t0 = time()
chain = infer([txl_fit, fit], data; n_samples=500)
println("   Inference complete in $(round(time() - t0; digits=1))s")
flush(stdout)

# --- Report posterior recovery vs ground-truth ---
truth = Dict(
    :k_tx => 1.0, :k_tl => 2.0, :gamma_mRNA => 0.5, :gamma_protein => 0.1,
    :k_atp => 1.0, :k_ntp => 0.2, :k_aa => 0.4, :sigma_obs => 0.3,
)

println("\n3. Posterior summary (90% CI):")
println("   parameter        truth   median   CI_low   CI_high   in_CI")
println("   " * "-"^60)
for (name, val) in truth
    samples = vec(chain[string(name)].data)
    med = quantile(samples, 0.5)
    lo, hi = quantile(samples, [0.05, 0.95])
    flag = lo <= val <= hi ? "yes" : "NO"
    println("   $(rpad(string(name), 16)) $(rpad(string(val), 7)) " *
            "$(rpad(round(med; digits=3), 8)) " *
            "$(rpad(round(lo;  digits=3), 8)) " *
            "$(rpad(round(hi;  digits=3), 9)) $flag")
end

# Heuristic narrative: shared TX/TL rates and AA/NTP rates (no direct dependence
# on the missing enzyme physics) should land in CI. `k_atp` will absorb the
# Hill + leak mismatch — expect it shifted from truth=1.0.
println("\n=== Tier B refit done at $(Dates.now()) ===")
