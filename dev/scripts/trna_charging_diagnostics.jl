# Spec §11 tasks 9.8 and 9.9: what the asserted tRNA pool size does, and whether
# the lumped charging step's stoichiometry matters.
#
# 9.8. The pool size is one number with three consequences, and all three are
#      ours (§4 D14): the charging flux, the buffer that check 7's
#      charged-tRNA counter draws on, and the lag it puts on the adenylate
#      forward channel. It is scanned at the default charged fraction, 0.8, in
#      two ways. With k_chg and the demand double held at their defaults, the
#      pool sets the flux. With both re-derived for each pool, as D14 does,
#      the flux is pinned and the pool sets only the buffer and the lag. Check
#      1b's 500-particle floor bounds the pool from below. Check 7's census
#      needs phase 11's translation, so the buffer here is the analytic seconds
#      of demand, labelled as such, and phase 11's task 11.7 supplies the census.
#
# 9.9. The published AMP + PPi products against a 2 ATP -> 2 ADP + 2 Pi
#      lumping, matched on steady flux rather than event count (D14). The
#      comparison records the adenylate kinase traffic and the ATP/ADP ratio
#      the 60 s rebuild reads.
#
# Usage: julia --project dev/scripts/trna_charging_diagnostics.jl
# Regenerate: sbatch dev/scripts/trna_charging_diagnostics.slurm

using InferCell
using Printf
using Dates
using OrdinaryDiffEq: remake

include(joinpath(@__DIR__, "..", "..", "test", "nucleotide_test_models.jl"))
include(joinpath(@__DIR__, "..", "..", "test", "trna_test_models.jl"))

# Provenance is captured before the measurements, so a HEAD that moves during a
# long job cannot be misattributed.
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test dev/scripts`, String)))

const N = corea_particles_per_mM()
const DEMAND = TL_DEMAND_PER_S                      # residues/s, from the doubles
const F0 = CHARGING_POOL_DEFAULTS.charged_fraction
const POOLS = [0.025, 0.05, 0.1, 0.125, 0.25, 0.5, 1.0]
const FLOOR = 500                                  # check 1b, particles
const LONGEST_PROTEIN = 746                        # residues; spec task 11.1

"Steady charging flux over the last 600 s, in residues per second."
steady_flux(ms, sol) = window_flux(ms, sol, CYCLE_S - 600)

# The default demand double's constant, which the fixed-constant scan holds.
const K_TL0 = TranslationDemand().k_tl

# ---------------------------------------------------------------------------
# 9.8 — the pool-size scan
# ---------------------------------------------------------------------------

rows = []
for pool in POOLS
    # (a) Fixed constants: the defaults' k_chg and k_tl, and a different pool.
    msa = charging_models(charging = TrnaCharging(pool_mM = pool,
                                                  k_chg = derive_k_chg()),
                          demand = TranslationDemand(k_tl = K_TL0))
    sa = recycling_solve(msa)
    # (b) Re-derived: k_chg by D14, and k_tl so the double takes the demand
    #     at this pool's nominal charged pool.
    k_tl = TL_DEMAND_MM_PER_S / (F0 * pool)
    mb = TrnaCharging(pool_mM = pool)
    msb = charging_models(charging = mb, demand = TranslationDemand(k_tl = k_tl))
    sb = recycling_solve(msb)
    U, C = layout_index(msb, :M_trna_c), layout_index(msb, :M_trna_chg_c)
    ATP = layout_index(msb, :M_atp_c)
    unc_min = minimum(u[U] for u in sb.u) * N
    chg_min = minimum(u[C] for u in sb.u) * N
    atp = sb.u[end][ATP]
    k = charging_derivation(mb).k_chg
    # The low-pass filter the pool puts on the adenylate forward channel:
    # linearised, the charged pool relaxes at k_chg·[ATP] + k_tl.
    tau = 1 / (k * atp + k_tl)
    push!(rows, (pool = pool, particles = pool * N,
                 flux_fixed = steady_flux(msa, sa),
                 atp_fixed = sa.u[end][ATP],
                 flux_derived = steady_flux(msb, sb), atp_derived = atp,
                 k_chg = k, unc_min = unc_min, chg_min = chg_min,
                 buffer_s = sb.u[end][C] * N / DEMAND,
                 longest_covered = sb.u[end][C] * N / LONGEST_PROTEIN,
                 tau = tau, floor_ok = min(unc_min, chg_min) >= FLOOR))
end

# The lag, measured once at the default pool rather than only linearised: from
# the steady state, step the demand up 10% and time the charging flux's rise to
# 63% of its value 30 s later. ATP responds too, on a much slower scale, which
# is why the reference is 30 s and not the end of the cycle.
ms0 = charging_models()
s0 = recycling_solve(ms0; horizon = 1200.0)
step_ms = charging_models(demand = TranslationDemand(k_tl = 1.1 * K_TL0))
step_prob = remake(build_problem(step_ms; tspan = (0.0, 30.0)); u0 = s0.u[end])
step_sol = solve(step_prob, Rodas5P(); abstol = ABSTOL_R, reltol = RELTOL_R, saveat = 0.01)
m0 = TrnaCharging()
Us, As = layout_index(step_ms, :M_trna_c), layout_index(step_ms, :M_atp_c)
vflux(u) = charging_flux(SA[u[Us], 0.0], Float64[], 0.0, m0, SA[u[As]]) * N
v_start, v_30 = vflux(step_sol.u[1]), vflux(step_sol.u[end])
i63 = findfirst(u -> (vflux(u) - v_start) >= 0.632 * (v_30 - v_start), step_sol.u)
tau_measured = step_sol.t[i63]
tau_linear = only(r.tau for r in rows if r.pool == CHARGING_POOL_DEFAULTS.pool_mM)

# ---------------------------------------------------------------------------
# 9.6 — the tRNA ladder, every rung that was run
# ---------------------------------------------------------------------------
#
# The suite asserts six decades ending at the pinned tolerances. The residual
# is roundoff and differs between machines, so the numbers of record are these,
# from a Slurm run, over all nine rungs including the two beyond the asserted
# window. Reporting every rung run is what §3's amendment requires.

const LADDER = [(1e-4, 1e-2), (1e-5, 1e-3), (1e-6, 1e-4), (1e-7, 1e-5), (1e-8, 1e-6),
                (1e-9, 1e-7), (1e-10, 1e-8), (1e-11, 1e-9), (1e-12, 1e-10)]
lms = charging_models()
LU, LC = layout_index(lms, :M_trna_c), layout_index(lms, :M_trna_chg_c)
ltr(u) = u[LU] + u[LC]
ladder = map(LADDER) do (a, r)
    l = recycling_solve(lms; abstol = a, reltol = r)
    d = max_drift(l, ltr)
    (a = a, r = r, drift = d, ulps = d / eps(ltr(l.u[1])),
     orders = log10(tolerance_bound(l, ((LU, 1), (LC, 1)); abstol = a, reltol = r) / d))
end
lprob = build_problem(lms; tspan = (0.0, CYCLE_S))
lsol = recycling_solve(lms)
bitwise = count(u -> (d = lprob.f(u, lprob.p, 0.0); d[LU] + d[LC] === 0.0), lsol.u)
asserted = ladder[1:7]
spread_asserted = maximum(x.drift for x in asserted) / minimum(x.drift for x in asserted)

# ---------------------------------------------------------------------------
# 9.9 — the stoichiometry comparison, matched on steady flux
# ---------------------------------------------------------------------------

function lumping_run(charging)
    ms = charging_models(charging = charging)
    sol = recycling_solve(ms)
    rec = ms[1]
    ix(s) = layout_index(ms, s)
    # ADK1's net rate at each save point: the recycling module evaluated on its
    # own states and the held glycolytic inputs, which is exactly what the
    # composed right-hand side evaluates.
    adk1(u) = recycling_fluxes(SA[(u[ix(s)] for s in RECYCLING_STATES)...], Float64[],
                               0.0, rec, SA[(u[ix(s)] for s in RECYCLING_INPUTS)...])[3]
    late = sol.u[findfirst(t -> t >= CYCLE_S - 600, sol.t):end]
    u = sol.u[end]
    return (flux = steady_flux(ms, sol),
            adk1 = sum(adk1, late) / length(late) * N,
            atp = u[ix(:M_atp_c)], adp = u[ix(:M_adp_c)], amp = u[ix(:M_amp_c)],
            ppi = u[ix(:M_ppi_c)], ratio = u[ix(:M_atp_c)] / u[ix(:M_adp_c)])
end
pub = lumping_run(TrnaCharging())
two_nominal = lumping_run(TwoAtpCharging())

# Matched at the steady state, not at nominal ATP: rescale k2 by secant until the
# two-ATP form's steady flux equals the published form's. k2 derived at nominal
# ATP undershoots, because ATP settles below nominal and this law is second
# order in it (spec §12, 2026-09-23 review correction).
function match_k2(target; tol = 1e-6)
    s0, s1 = 1.0, 1.02
    f0 = lumping_run(TwoAtpCharging(k2_scale = s0)).flux - target
    f1 = lumping_run(TwoAtpCharging(k2_scale = s1)).flux - target
    for _ in 1:20
        abs(f1) < tol * target && break
        s0, s1, f0 = s1, s1 - f1 * (s1 - s0) / (f1 - f0), f1
        f1 = lumping_run(TwoAtpCharging(k2_scale = s1)).flux - target
    end
    return s1
end
k2_scale = match_k2(pub.flux)
two = lumping_run(TwoAtpCharging(k2_scale = k2_scale))

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

jobid = get(ENV, "SLURM_JOB_ID", "not under Slurm")
io = IOBuffer()
println(io, "# tRNA charging diagnostics (spec tasks 9.8, 9.9)")
println(io)
println(io, "Generated by `dev/scripts/trna_charging_diagnostics.jl` at commit $COMMIT",
        DIRTY ? " (WARNING: tree dirty)" : "",
        ", Slurm job $jobid, on $(gethostname()), $(Dates.now()), Julia $(VERSION).")
println(io)
println(io, "Regenerate with `sbatch dev/scripts/trna_charging_diagnostics.slurm`.")
println(io)
println(io, "Composition: `NucleotideRecycling` + `HeldGlycolytic` + `TrnaCharging` (with an")
println(io, @sprintf("exact cumulative counter) + `TranslationDemand`, 6,300 s, Rodas5P, abstol %.0e,", ABSTOL_R))
println(io, @sprintf("reltol %.0e. Demand %.2f residues/s, from 3,484,518 over 6,300 s. %.2f particles/mM.",
                     RELTOL_R, DEMAND, N))
println(io)
println(io, "## 9.8 — the pool size's three consequences")
println(io)
println(io, @sprintf("Charged fraction %.1f throughout. \"Fixed\" holds k_chg = %.5f /(mM·s) and the",
                     F0, derive_k_chg()))
println(io, "demand double's k_tl at their defaults. \"Derived\" re-derives both for each")
println(io, "pool, as D14 does. The floor is check 1b's 500 particles, applied to each tRNA")
println(io, "state's minimum over the cycle. The buffer is seconds of demand in the charged")
println(io, "pool at the end of the cycle. It is analytic: check 7's census needs phase 11.")
println(io)
println(io, "| pool (mM) | particles | flux, fixed (/s) | flux, derived (/s) | k_chg derived | min uncharged | min charged | 1b floor | buffer (s) | longest proteins covered | lag τ (s) |")
println(io, "|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows
    println(io, @sprintf("| %.3f | %.0f | %.1f | %.2f | %.4f | %.0f | %.0f | %s | %.2f | %.1f | %.2f |",
                         r.pool, r.particles, r.flux_fixed, r.flux_derived, r.k_chg,
                         r.unc_min, r.chg_min, r.floor_ok ? "clears" : "**fails**",
                         r.buffer_s, r.longest_covered, r.tau))
end
println(io)
println(io, @sprintf("Lag at the default pool: linearised τ = %.3f s; measured, by a 10%% demand step", tau_linear))
println(io, @sprintf("from steady state, the charging flux reaches 63%% of its 30 s response in **%.2f s**.", tau_measured))
println(io, @sprintf("Against the 1 s handshake that is %.2fx, and against the 60 s rebuild %.4fx. τ",
                     tau_measured, tau_measured / 60))
println(io, "scales linearly with the pool in the table's derived column.")
println(io)
println(io, "## 9.6 — the tRNA ladder, every rung run")
println(io)
println(io, @sprintf("The summed tRNA derivative is bitwise 0.0 at %d of %d save points (first gate).",
                     bitwise, length(lsol.u)))
println(io, "Drift is roundoff and machine-dependent; these are this job's values. The suite")
println(io, "asserts the first seven rungs (six decades, ending at the pinned tolerances).")
println(io)
println(io, "| abstol | reltol | drift (mM) | ulps of the sum | orders below tol_C | asserted |")
println(io, "|---|---|---|---|---|---|")
for (i, x) in enumerate(ladder)
    println(io, @sprintf("| %.0e | %.0e | %.3e | %.0f | %.2f | %s |", x.a, x.r, x.drift,
                         x.ulps, x.orders, i <= 7 ? "yes" : "reported only"))
end
println(io)
println(io, @sprintf("Spread over the asserted rungs: %.1fx (bound 100x).", spread_asserted))
println(io)
println(io, "## 9.9 — AMP + PPi against 2 ATP -> 2 ADP + 2 Pi, matched on flux")
println(io)
println(io, "Both forms are derived to deliver the demand at nominal ATP. ATP settles below")
println(io, "nominal, and the two-ATP law is second order in it, so at that derivation its")
println(io, @sprintf("steady flux undershoots; k2 is then rescaled by %.5f so the steady fluxes match.", k2_scale))
println(io, "Values over the last 600 s of the cycle, and the state at its end. The")
println(io, "composition uses phase 8's glycolytic double, not live glycolysis.")
println(io)
println(io, "| lumping | charging flux (/s) | ADK1 net (/s) | ATP (mM) | ADP (mM) | AMP (mM) | PPi (mM) | ATP/ADP |")
println(io, "|---|---|---|---|---|---|---|---|")
for (name, r) in (("AMP + PPi (published)", pub),
                  ("2 ATP, k2 at nominal ATP (not matched)", two_nominal),
                  ("2 ATP, k2 matched at steady state", two))
    println(io, @sprintf("| %s | %.2f | %.2f | %.4f | %.4f | %.4f | %.4f | %.4f |",
                         name, r.flux, r.adk1, r.atp, r.adp, r.amp, r.ppi, r.ratio))
end
println(io)
println(io, @sprintf("Matched, the fluxes agree to %.1e and the ATP/ADP ratio differs by %.3f%%",
                     abs(two.flux / pub.flux - 1), 100 * (two.ratio / pub.ratio - 1)))
println(io, @sprintf("between the forms (%.2f%% unmatched).", 100 * (two_nominal.ratio / pub.ratio - 1)))

out = joinpath(@__DIR__, "trna_charging_diagnostics_result.md")
write(out, String(take!(io)))
println("wrote ", out)
print(read(out, String))
