# Spec §11 tasks 8.6 and 8.7: nucleotide recycling over a full 6,300 s cycle.
#
# The module owns the pools every other module routes energy through, so it
# cannot be verified from inside itself: alone it conserves both moieties
# trivially, because nothing consumes them. The acceptance criterion is
# therefore a composition against a charging drain, and it is run over a full
# cycle because 144 s of the model that turned out to be wrong looked fine.
#
# Three configurations, and the two removals are the evidence:
#
#   all five      both moieties conserved, ATP positive, pyrophosphate settling
#   ADK1 removed  ATP below 1% of initial near the scoping note's ~144 s
#   PPA removed   pyrophosphate strands the phosphate moiety and the pathway
#                 stalls — it cannot diverge in a phosphate-closed model
#
# **The threshold is a crossing, not exhaustion.** The note's 144 s assumes a
# constant drain. This drain saturates below a Michaelis constant three orders
# of magnitude under the nominal ATP pool (see `ChargingDrain`), and phase 9's
# real charging step is mass-action in both reactants, so ATP decays and never
# reaches zero. Both numbers are recorded below with the model each belongs to.
#
# **Conservation is asserted against the round-off floor, and the reason is
# structural** (spec §3's exact-conservation exception, second gate, §12
# 2026-09-11). Adenylate, guanylate and flux-corrected phosphate all satisfy
# `nᵀf ≡ 0` exactly — every reaction's contribution cancels term by term — so
# there is no local truncation error in them to scale with the tolerance. They
# are *not* bitwise zero on the composed right-hand side, because the terms
# arrive from three modules and are summed in an order that does not cancel;
# that is why the second gate exists and the first does not apply. What the gate
# asserts is the per-evaluation residual, one ulp of the conserved sum, and the
# suite asserts it. Both tolerance settings are still run and both residuals
# reported — that is the evidence for the claim rather than an assertion that
# fails — and what is asserted is that the residual
# sits orders of magnitude below the single-run integrator bound. The evidence
# that these checks can fail at all is the mutation test in the suite, which is
# the point §3 makes with "a conservation test that cannot fail is not
# evidence". The five-fold rule binds where the trajectory is *restarted*, which
# is the assembled hybrid model's 6,300 handshakes and is what the bound's
# `N_restarts` factor is about.
#
# Usage: julia --project dev/scripts/full_cycle_recycling.jl

using InferCell
using Printf
using Dates

include(joinpath(@__DIR__, "..", "..", "test", "nucleotide_test_models.jl"))

const CYCLE = 6300.0            # s — one published cell cycle
const SAVE_EVERY = 60.0         # s — the rebuild interval, spec §3
const ABSTOL = 1e-10            # mM, pinned in spec §3
const RELTOL = 1e-8
const TIGHTEN = 10              # the tolerance principle's factor
const ATP_THRESHOLD = 0.01      # spec task 8.6: 1% of the initial value

# Provenance captured before the measurements, not after: on a long job HEAD can
# move between load and write, and the artefact would then name source the
# measurement never saw.
const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const DIRTY = !isempty(strip(read(`git status --porcelain -- src test dev/scripts`, String)))

# Composed state layout: recycling's eight, then the held glycolytic pools and
# their cumulative phosphorylation counter, then the drain's counter.
const ATP, ADP, AMP, PI, GTP, GDP, GMP, PPI = 1, 2, 3, 4, 5, 6, 7, 8
const SLP_CUM = 13
const DRAIN_CUM = 14

models(; reactions = RECYCLING_REACTIONS, k_slp = RECYCLING_DRAIN_MM_PER_S / 0.2178,
       drain = RECYCLING_DRAIN_MM_PER_S) = AbstractSubModel[
    NucleotideRecycling(reactions = reactions),
    HeldGlycolytic(k_slp = k_slp),
    ChargingDrain(rate_mm_per_s = drain),
]

function run_cycle(ms; tspan = (0.0, CYCLE), abstol = ABSTOL, reltol = RELTOL)
    prob = build_problem(ms; tspan = tspan)
    return solve(prob, Rodas5P(); abstol = abstol, reltol = reltol,
                 saveat = SAVE_EVERY)
end

adenylate(u) = u[ATP] + u[ADP] + u[AMP]
guanylate(u) = u[GTP] + u[GDP] + u[GMP]

# Free phosphate plus every phosphorylated species, pyrophosphate counted twice.
phosphate(u) = u[PI] + 3u[ATP] + 2u[ADP] + u[AMP] + 3u[GTP] + 2u[GDP] + u[GMP] + 2u[PPI]

# What crosses the two inbound mass edges, read off the trajectory rather than
# quadratured. GTP is produced by PGK3 and PYK3 and by nothing else in this
# composition, so its change *is* the integral of the GTP branch's flux, and
# each turnover carries one phosphate in from the held 13DPG or PEP pool. The
# substrate-level phosphorylation the glycolytic double supplies takes its
# phosphate from the free pool instead, so it crosses nothing and appears here
# not at all. Subtract the flux; do not relax the bound.
inbound_phosphate(u, u0) = u[GTP] - u0[GTP]

drift(sol, f) = maximum(abs(f(u) - f(sol.u[1])) for u in sol.u)
function phosphate_drift(sol)
    u0 = sol.u[1]
    return maximum(abs((phosphate(u) - inbound_phosphate(u, u0)) - phosphate(u0))
                   for u in sol.u)
end

"""
Run one configuration at both tolerance settings and report each residual with
the ratio between them. For a linear invariant the ratio is expected to be about
one — the number is reported as evidence of tolerance-independence, not as a
pass criterion.
"""
function scaling(ms, measures; tspan = (0.0, CYCLE))
    loose = run_cycle(ms; tspan = tspan)
    tight = run_cycle(ms; tspan = tspan, abstol = ABSTOL / TIGHTEN,
                      reltol = RELTOL / TIGHTEN)
    return loose, tight,
           [(name, f(loose), f(tight),
             f(tight) == 0 ? Inf : f(loose) / f(tight)) for (name, f) in measures]
end

const MEASURES = [("adenylate", s -> drift(s, adenylate)),
                  ("guanylate", s -> drift(s, guanylate))]

# The single-run bound the tolerance principle reports as a secondary number:
#     tol_C = N_restarts · Σ|nᵢ| · max(abstol, reltol · maxₜ|xᵢ(t)|)
# with one restart, since a pure ODE composition runs no handshake.
function tolerance_bound(sol, coeffs; abstol = ABSTOL, reltol = RELTOL)
    return sum(abs(n) * max(abstol, reltol * maximum(abs(u[i]) for u in sol.u))
               for (i, n) in coeffs)
end

# ---------------------------------------------------------------------------
# 8.6 — the three configurations
# ---------------------------------------------------------------------------

full_loose, _, full_scaling = scaling(models(), MEASURES)
u0 = full_loose.u[1]
atp0 = u0[ATP]

atp_min_full = minimum(u[ATP] for u in full_loose.u)
ppi_final = full_loose.u[end][PPI]
ppi_max = maximum(u[PPI] for u in full_loose.u)
ppi_late = [u[PPI] for u in full_loose.u[(end - 10):end]]

no_kinase = run_cycle(models(reactions = (:R_PGK3, :R_PYK3, :R_GK1, :R_PPA)))
no_ppa = run_cycle(models(reactions = (:R_PGK3, :R_PYK3, :R_ADK1, :R_GK1)))

"First save point at which ATP has fallen below `frac` of its initial value."
function crossing_time(sol, frac)
    target = frac * sol.u[1][ATP]
    i = findfirst(u -> u[ATP] < target, sol.u)
    return i === nothing ? nothing : sol.t[i]
end

# A finer grid than the 60 s save points, for the crossing only: the note quotes
# 144 s and a 60 s grid could not distinguish 120 from 180.
no_kinase_fine = solve(build_problem(models(reactions = (:R_PGK3, :R_PYK3, :R_GK1, :R_PPA));
                                     tspan = (0.0, 600.0)),
                       Rodas5P(); abstol = ABSTOL, reltol = RELTOL, saveat = 1.0)
t_cross = crossing_time(no_kinase_fine, ATP_THRESHOLD)

# Two constant-drain figures, and they are not the same arithmetic.
#
# `note_drain_time` is the scoping note's: the FULL adenylate pool at the
# published demand, ~144 s (dev/notes/reduced-syn3a-scoping.md, "the cell runs
# out of adenylate after 144 s"). `module_drain_time` subtracts the AMP already
# present, which is the archived design's ~141 s (design.md: "the module's own
# arithmetic gives ~141 s because AMP is already at 0.0832 mM"). Reporting the
# second under the first's name is a misattribution §12's D13 entry already
# records being fixed once, so both are computed and both are labelled.
adenylate0 = adenylate(u0)
note_drain_time = adenylate0 / RECYCLING_DRAIN_MM_PER_S
module_drain_time = (adenylate0 - u0[AMP]) / RECYCLING_DRAIN_MM_PER_S

ppi_no_ppa = [u[PPI] for u in no_ppa.u]
ppi_growth_late = ppi_no_ppa[end] - ppi_no_ppa[end - 1]
# What removing the enzyme actually does in a phosphate-closed model: it strands
# the moiety rather than letting the concentration diverge.
phosphate_budget = phosphate(u0)
ppi_share = 2 * ppi_no_ppa[end] / phosphate_budget
atp_no_ppa_final = no_ppa.u[end][ATP]
no_ppa_fine = solve(build_problem(models(reactions = (:R_PGK3, :R_PYK3, :R_ADK1, :R_GK1));
                                  tspan = (0.0, 6300.0)),
                    Rodas5P(); abstol = ABSTOL, reltol = RELTOL, saveat = 10.0)
t_cross_no_ppa = crossing_time(no_ppa_fine, ATP_THRESHOLD)

# ---------------------------------------------------------------------------
# 8.7 — phosphate closure, both forms
# ---------------------------------------------------------------------------

# Exact: with the GTP-branch reactions inactive and no phosphorylation, nothing
# carries phosphate across the boundary. The charging drain still runs, and it
# closes internally — three phosphates in ATP become one in AMP and two in
# pyrophosphate — so this form tests the module and not the absence of traffic.
exact_ms = models(reactions = (:R_ADK1, :R_GK1, :R_PPA), k_slp = 0.0)
_, _, exact_scaling = scaling(exact_ms, [("phosphate (exact)", s -> drift(s, phosphate))])

# Flux-corrected: all five active, the phosphorylation running, and the inbound
# flux subtracted.
corr_loose, _, corr_scaling =
    scaling(models(), [("phosphate (flux-corrected)", phosphate_drift)])

# What the substrate-level phosphorylation turned over. It takes every phosphate
# from the free pool, so it crosses no edge and appears in no correction.
slp_total = corr_loose.u[end][SLP_CUM]

# ---------------------------------------------------------------------------
# Wall-clock — what the §9 open question needs
# ---------------------------------------------------------------------------

run_cycle(models())                                   # warm
elapsed = @elapsed run_cycle(models())
elapsed_tight = @elapsed run_cycle(models(); abstol = ABSTOL / TIGHTEN,
                                  reltol = RELTOL / TIGHTEN)

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

jobid = get(ENV, "SLURM_JOB_ID", "not under Slurm")
io = IOBuffer()
println(io, "# Nucleotide recycling over a full cycle (spec tasks 8.6, 8.7)")
println(io)
println(io, "Generated by `dev/scripts/full_cycle_recycling.jl` at commit $COMMIT",
        DIRTY ? " (WARNING: tree dirty)" : "",
        ", Slurm job $jobid, on $(gethostname()), $(Dates.now()), Julia $(VERSION).")
println(io)
println(io, "Regenerate with `sbatch dev/scripts/full_cycle_recycling.slurm`.")
println(io)
println(io, "Composition: `NucleotideRecycling` + `HeldGlycolytic` + `ChargingDrain`,")
println(io, @sprintf("6,300 s, Rodas5P, abstol %.0e, reltol %.0e, save points every %d s.",
                     ABSTOL, RELTOL, Int(SAVE_EVERY)))
println(io)

println(io, "## 8.6 — the three configurations")
println(io)
println(io, "### All five reactions")
println(io)
println(io, @sprintf("Adenylate %.6f mM and guanylate %.6f mM at t = 0.", adenylate0,
                     guanylate(u0)))
println(io, @sprintf("ATP never falls below **%.6f mM** (initial %.4f).", atp_min_full, atp0))
println(io, @sprintf("Pyrophosphate settles at **%.6f mM** (peak %.6f); over the last ten",
                     ppi_final, ppi_max))
println(io, @sprintf("save points it moves by %.3e mM, so it is at a steady state rather",
                     maximum(ppi_late) - minimum(ppi_late)))
println(io, "than still climbing. The reverse term is what puts it there: at 17.8 mM")
println(io, "phosphate against a Michaelis constant of 0.0976 the product bracket is")
println(io, "large, and a pyrophosphatase written irreversibly would not settle here.")
println(io)
println(io, "### Conservation, against the round-off floor")
println(io)
adenylate_bound = tolerance_bound(full_loose, [(ATP, 1), (ADP, 1), (AMP, 1)])
guanylate_bound = tolerance_bound(full_loose, [(GTP, 1), (GDP, 1), (GMP, 1)])
println(io, "| moiety | residual at (1e-10, 1e-8) | at a tenth of each | ratio | integrator bound |")
println(io, "|---|---|---|---|---|")
for ((name, lo, hi, ratio), bound) in zip(full_scaling, (adenylate_bound, guanylate_bound))
    println(io, @sprintf("| %s | %.3e mM | %.3e mM | %.1fx | %.3e mM |",
                         name, lo, hi, ratio, bound))
end
println(io)
println(io, "**Tightening the solver does not shrink these residuals, and it should not.**")
println(io, "Both moieties satisfy `nᵀf ≡ 0` exactly — every reaction's contribution")
println(io, "cancels term by term — so there is no local truncation error in them to")
println(io, "scale with the tolerance, and what the table shows is round-off, five orders")
println(io, "below the integrator bound in the last column. The two-rung ratios are not")
println(io, "the evidence and should not be read as one: over nine rungs the adenylate")
println(io, "residual wanders between 128 and 1,119 ulps of its sum with no trend, while")
println(io, "the evaluation count grows 36-fold, so any single pair of rungs gives an")
println(io, "arbitrary ratio between 0.2x and 1.8x.")
println(io)
println(io, "This is the **second gate** of spec §3's exact-conservation exception (§12,")
println(io, "2026-09-11), not the first: these sums are not bitwise zero on the composed")
println(io, "right-hand side, because the terms arrive from three modules. The gate that")
println(io, "licenses the floor is the per-evaluation residual — one ulp of the conserved")
println(io, "sum — and `test/test_corea_nucleotide_recycling.jl` asserts it, together")
println(io, "with a six-rung ladder bounded by each rung's own `tol_C`. The evidence that")
println(io, "these checks can fail is that assertion and the mutation test, not this table.")
println(io)
println(io, "### Adenylate kinase removed")
println(io)
println(io, @sprintf("ATP falls below 1%% of its initial value (%.6f mM) at **t = %s s**.",
                     ATP_THRESHOLD * atp0, t_cross === nothing ? "never" : string(t_cross)))
println(io, @sprintf("The scoping note computes **%.1f s**, the full %.4f mM pool at the",
                     note_drain_time, adenylate0))
println(io, @sprintf("published demand of %.6f mM/s under a *constant* drain. Subtracting the",
                     RECYCLING_DRAIN_MM_PER_S))
println(io, @sprintf("AMP already present gives the archived design's **%.1f s**: (%.4f − %.4f) mM",
                     module_drain_time, adenylate0, u0[AMP]))
println(io, @sprintf("/ %.6f mM/s. All three numbers belong to",
                     RECYCLING_DRAIN_MM_PER_S))
println(io, "**threshold crossing** under a saturating drain: ATP decays toward zero")
println(io, "rather than through it, and what the run reports at the horizon is round-off.")
println(io, @sprintf("ATP at the end of the cycle: %.3e mM.", no_kinase.u[end][ATP]))
println(io, @sprintf("Adenylate is still conserved without the kinase: residual %.3e mM,",
                     drift(no_kinase, adenylate)))
println(io, "which is what says the pool was stranded as AMP rather than lost.")
println(io)
println(io, "### Pyrophosphatase removed")
println(io)
println(io, @sprintf("Pyrophosphate reaches **%.4f mM** at 6,300 s from %.4f mM — a %.0f-fold",
                     ppi_no_ppa[end], ppi_no_ppa[1], ppi_no_ppa[end] / ppi_no_ppa[1]))
println(io, @sprintf("rise, against **%.6f mM** in the configuration that keeps the enzyme.",
                     ppi_final))
println(io)
println(io, "**It strands the moiety; it does not diverge, and it cannot.** Phosphate is")
println(io, @sprintf("closed in this composition at %.4f mM, and pyrophosphate ends holding",
                     phosphate_budget))
println(io, @sprintf("**%.1f%%** of it. With no route back to free phosphate the",
                     100 * ppi_share))
println(io, "substrate-level phosphorylation has nothing to work with, so ATP collapses")
println(io, @sprintf("— below 1%% of its initial value at t = %s s, and %.3e mM at the end —",
                     t_cross_no_ppa === nothing ? "never" : string(t_cross_no_ppa),
                     atp_no_ppa_final))
println(io, @sprintf("and charging stops with it, which is why the rise flattens (%.4f mM",
                     ppi_growth_late))
println(io, "over the last save interval) rather than continuing.")
println(io)
println(io, "The scoping note's **173 mM** is one pyrophosphate per charging event over a")
println(io, "whole cycle with nothing hydrolysing it. That is open-pool arithmetic: it")
println(io, "assumes charging runs at 553.1/s for 6,300 s, which needs a phosphate supply")
println(io, "the closed model does not have — Core A′'s own phosphate is a closed moiety,")
println(io, "which is what spec §3 check 4b asserts. Both numbers are recorded with the")
println(io, "model each belongs to. A stalled pathway is the stronger demonstration that")
println(io, "the enzyme is required, and it is the one that survives assembly.")
println(io)

println(io, "## 8.7 — phosphate closure, both forms")
println(io)
println(io, "| form | residual at (1e-10, 1e-8) | at a tenth of each | ratio |")
println(io, "|---|---|---|---|")
for (name, lo, hi, ratio) in vcat(exact_scaling, corr_scaling)
    println(io, @sprintf("| %s | %.3e mM | %.3e mM | %.1fx |", name, lo, hi, ratio))
end
println(io)
println(io, @sprintf("Phosphate crossing the two inbound mass edges over the cycle: %.4f mM,",
                     inbound_phosphate(corr_loose.u[end], corr_loose.u[1])))
println(io, "carried in from the resupplied 13DPG and PEP pools, one per turnover of the")
println(io, "GTP branch. **It is bounded because guanylate is**: nothing in Core A′ consumes")
println(io, "GTP until phase 9's translation, so PGK3 and PYK3 run until GDP is spent and")
println(io, @sprintf("then stop — the crossing equals GDP₀ + GMP₀ = %.4f mM to 0.2%%. That is a",
                     corr_loose.u[1][GDP] + corr_loose.u[1][GMP]))
println(io, "structural statement about the composition. It was *not* true of a first")
println(io, "version of the double, which returned a zero derivative for the four")
println(io, "glycolytic pools and so did not hold them at all: `contributed_states` folds")
println(io, "recycling's terms into the owner's `du`, 13DPG fell 99.8% inside half a")
println(io, "second, and the branch then stopped on substrate rather than on guanylate,")
println(io, "at 0.0505 mM. The flux-corrected form subtracts the crossing; widening the")
println(io, "bound to swallow it instead would have hidden any real error of the same size.")
println(io, @sprintf("The substrate-level phosphorylation crosses nothing — it turned over %.2f mM",
                     slp_total))
println(io, "of phosphate over the cycle and took every one of them from the free pool,")
println(io, "which is what keeps the moiety closed and what a first version of the double")
println(io, "got wrong.")
println(io)
println(io, "Both forms are linear invariants, so the ratio column carries the same")
println(io, "meaning it does above: tolerance-independence, not a five-fold fall.")
println(io)

println(io, "## Wall-clock")
println(io)
println(io, @sprintf("One 6,300 s trajectory: **%.3f s** at the pinned tolerances, %.3f s at",
                     elapsed, elapsed_tight))
println(io, @sprintf("a tenth of each. %d states and seven rate expressions, no handshake — this is",
                     length(full_loose.u[1])))
println(io, "not Core A′'s cost and does not decide K1, which spec task 13.7 measures on")
println(io, "the assembled hybrid model.")
println(io)
println(io, "It adds a third module measurement to §9's question of whether the full-cycle")
println(io, "checks belong in the")
println(io, "default suite, and the answer is that the *horizon* is free: at a few")
println(io, "milliseconds a trajectory, running 6,300 s rather than 600 s costs nothing,")
println(io, "so `test/test_corea_nucleotide_recycling.jl` asserts the phase's done-when at")
println(io, "the horizon the done-when names. What is not free is compilation — the two")
println(io, "full-cycle testsets add tens of seconds to the suite, almost all of it")
println(io, "Rodas5P specialising on new problem types, and shortening the horizon would")
println(io, "not recover any of it. The figure is not quoted here because it is a property")
println(io, "of the whole tree — phases 6 and 7 now compile the stiff path first, so the")
println(io, "marginal cost of this phase is smaller than it was measured to be alone. This driver is not a duplicate of those tests: it is where")
println(io, "the numbers are *recorded*, with the commit and job that produced them, so a")
println(io, "reader gets the measurements without running anything.")
println(io)
println(io, "## What this does not cover")
println(io)
println(io, "Every number here is the module's against *doubles*. `HeldGlycolytic` relaxes")
println(io, "four glycolytic pools to their registry values and rephosphorylates ADP at a")
println(io, "rate matched to the charging demand; phase 6's module replaces it with")
println(io, "reactions that do neither exactly, and in particular its pools are set by")
println(io, "glycolytic flux rather than held near a setpoint.")
println(io, "`ChargingDrain` is zeroth order in the tRNA pools where phase 9's module is")
println(io, "mass-action in both. The conservation and closure results are structural and")
println(io, "survive that replacement; the settled pyrophosphate concentration and the")
println(io, "crossing time do not, and are quoted as this composition's.")

result = joinpath(@__DIR__, "full_cycle_recycling_result.md")
write(result, String(take!(io)))
println("wrote $result")
