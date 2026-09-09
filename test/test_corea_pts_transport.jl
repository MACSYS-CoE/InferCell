using Test
using InferCell
using Distributions
using StaticArrays: SVector
using OrdinaryDiffEq: Rodas5P, solve

import InferCell: states, parameters, dynamics

# Spec §11 phase 7 — the phosphotransferase cascade and lactate export, the
# first real Core A′ sub-model. `caught` comes from corea_test_models.jl and
# `_gidx` from contribution_test_models.jl, both included first by runtests.jl.

# The published copy numbers, and the carrier pairs whose sums they are.
const PTS_CARRIERS = ((:M_ptsi_c, :M_ptsi_P_c, 353),
                      (:M_ptsh_c, :M_ptsh_P_c, 290),
                      (:M_crr_c, :M_crr_P_c, 314),
                      (:M_ptsg_c, :M_ptsg_P_c, 831))

"""
    HeldMetabolites(; leak = 0.0, lac_source = 0.0)

The four foreign pools `PtsTransport` reads and writes, owned with **no
dynamics of their own**, so that the only thing moving them is the cascade.

"Held" means held by this double, not held constant: phase 1 made the
contribution channel execute, so `PtsTransport`'s terms reach these pools and
they move. That is the point — a pool nothing could change would make the
contribution channel untestable here.

This composition is **not** physiological and is not asserted to be. The
carriers relax to their own equilibrium in milliseconds and pyruvate is drawn
down doing it, because the published initial conditions put 95% of every
carrier in the phosphorylated form and the pools that would replenish
phosphoenolpyruvate live in phase 6. What the composition is for is
non-negativity and the four conservation invariants, neither of which depends
on the trajectory being biological.

`lac_source` stands in for glycolysis, which produces two lactate per glucose
and is the only reason cytosolic lactate has a steady state at all. At the
published 0.10961 mM/s it holds cytosolic lactate at 0.10961/0.075 = 1.46 mM,
which is the denominator spec §4 D6's one-percent criterion is stated against.

It carries one free parameter of its own, so `PtsTransport`'s free parameters
start at a non-zero offset in the composed vector. That is deliberate: the
module maps its own rate law's slots onto the *local* parameter slice, and an
off-by-one between local and global would be invisible in a composition where
the two coincide.
"""
struct HeldMetabolites <: AbstractSubModel
    params::Vector{InferParameter}
end

const HELD_SPECIES = (:M_g6p_c, :M_pep_c, :M_pyr_c, :M_lac__L_c)

function HeldMetabolites(; leak = 0.0, lac_source = 0.0)
    params = vcat(
        [InferParameter(leak, LogNormal(0.0, 1.0), false, :hold_leak,
                        :HeldMetabolites, :rate),
         InferParameter(lac_source, LogNormal(0.0, 1.0), false, :lac_source,
                        :HeldMetabolites, :rate)],
        [InferParameter(something(species_entry(s).initial_value, 0.1),
                        LogNormal(0.0, 1.0), true, Symbol(s, "0"),
                        :HeldMetabolites, :initial_condition)
         for s in HELD_SPECIES])
    return HeldMetabolites(params)
end

states(::HeldMetabolites) = collect(HELD_SPECIES)
parameters(m::HeldMetabolites) = m.params
# Species order is (g6p, pep, pyr, lac_c); only cytosolic lactate has a source.
dynamics(u, p, t, ::HeldMetabolites) =
    SVector(-p[1] * u[1], -p[1] * u[2], -p[1] * u[3], -p[1] * u[4] + p[2])

"""
    PtsCarrierLeak(inner)

`PtsTransport` with one cascade step's stoichiometry perturbed: GLCpts1 is
written as if it *created* phospho-HPr rather than transferring the phosphate
from phospho-EI, so the ptsH pair is no longer conserved while the other three
still are. This is the mutation that shows the conservation check can fail, and
fails locally.
"""
struct PtsCarrierLeak{N} <: AbstractSubModel
    inner::PtsTransport{N}
end

for f in (:states, :parameters, :inputs, :contributed_states, :coupling,
          :membrane_protein_states, :extracellular_states, :reduction_notes,
          :formalism, :inference_mode)
    @eval InferCell.$f(m::PtsCarrierLeak) = InferCell.$f(m.inner)
end
InferCell.module_id(::PtsCarrierLeak) = :PtsCarrierLeak
InferCell.contributions(u, p, t, m::PtsCarrierLeak, ui) =
    InferCell.contributions(u, p, t, m.inner, ui)

function dynamics(u, p, t, m::PtsCarrierLeak, u_inputs)
    du = dynamics(u, p, t, m.inner, u_inputs)
    q = InferCell._live_scalars(m.inner, p)
    _, v1, _, _, _, _ = InferCell._pts_rates(q, u, u_inputs)
    return Base.setindex(du, du[5] + v1, 5)   # M_ptsh_P_c created, not transferred
end

# Largest deviation of a carrier's two-form sum from its initial value, over
# every save point of `sol`.
function carrier_residual(sol, models, unphos, phos)
    i, j = _gidx(models, unphos), _gidx(models, phos)
    initial = sol.u[1][i] + sol.u[1][j]
    return maximum(abs((u[i] + u[j]) - initial) for u in sol.u)
end

# The single-run bound of spec §3: N_restarts * sum|n_i| * max(abstol,
# reltol * max_t|x_i|). The stoichiometric coefficients of a carrier sum are
# both one, and a standalone solve has no handshake restarts, so N_restarts = 1.
function carrier_bound(sol, models, unphos, phos; abstol, reltol)
    i, j = _gidx(models, unphos), _gidx(models, phos)
    peak = maximum(max(abs(u[i]), abs(u[j])) for u in sol.u)
    return 2 * max(abstol, reltol * peak)
end

function integrate_pts(models; tspan = (0.0, 600.0), abstol = 1e-10,
                       reltol = 1e-8, saveat = 60.0)
    prob = build_problem(models; tspan = tspan)
    return solve(prob, Rodas5P(); abstol = abstol, reltol = reltol,
                 saveat = saveat)
end

@testset "Phase 7 — phosphotransferase transport and lactate export" begin

    m = PtsTransport()
    models = [HeldMetabolites(), m]

    @testset "7.1 the two extracts, and what their shape guarantees" begin
        rates = read_source_table(joinpath(COREA_DATA_DIR, "pts_transport.tsv");
                                  file = "transport")
        ics = read_source_table(joinpath(COREA_DATA_DIR,
                                         "pts_initial_conditions.tsv");
                                file = "model_ics")

        @test length(rates.values) == 11
        @test length(ics.values) == 9

        # No GeometricStd column at all: that absence, not a per-row
        # declaration, is what makes every constant load as asserted.
        header = first(l for l in eachline(joinpath(COREA_DATA_DIR,
                                                    "pts_transport.tsv"))
                       if startswith(l, "!") && !startswith(l, "!!"))
        @test !occursin("GeometricStd", header)
        @test !occursin("Informedness", header)

        @test rates.values["KF_1_R_GLCpts1"] == 200000.0
        @test rates.values["KR_4_R_GLCpts4"] == 1.00e-05
        @test rates.values["P_R_L_LACt2r"] == 5.00e-09
        @test rates.values["KF_4_R_GLCpts4"] == 0.88
        @test all(rates.informedness[id] === :asserted for id in PTS_RATE_IDS)

        # Disjoint identifier sets, so the report is empty -- and an empty
        # report is a positive result rather than a failure to run.
        @test isempty(ambiguity_report([rates, ics]))

        # Every one of the eleven reaches the composed model's asserted-prior
        # enumeration, and a genuinely balanced parameter does not. The balanced
        # parameter stands in for phase 6's glycolytic imports, which do not
        # exist yet on this branch.
        balanced = InferParameter(319.5, LogNormal(log(319.5), log(1.051)), true,
                                  :kcat_fwd_PGK3, :Glycolysis, :rate,
                                  ParameterSource("central_balanced";
                                                  informedness = :balanced))
        enumerated = [l.subject for l in
                      reduction_declarations([CoreAStub(:Glycolysis;
                                                        params = [balanced]), m])
                      if l.category === :asserted_prior]
        for name in PTS_SCALARS
            @test name in enumerated
        end
        @test !(:kcat_fwd_PGK3 in enumerated)

        # External glucose is in that list, and belongs there: 40 mM is a point
        # value with no quantified uncertainty, which is what the registry
        # records for it too. That the *prior* would be ours is a different
        # claim from the *clamp* being ours, and the clamp is the published
        # model's — no :clamp deviation is registered for it. Task 7.5 asserts
        # that half.
        @test informedness(only(filter(q -> q.name === :glc_e_mM,
                                       parameters(m)))) === :asserted
        @test species_entry(:M_glc__D_e).informedness === :asserted

        @test unique(source_file.(rate_params(parameters(m)))) ⊇ ["transport"]
        @test unique(source_file.(ic_params(parameters(m)))) == ["model_ics"]
    end

    @testset "7.2 the rate laws, hand-checked" begin
        @test length(states(m)) == 9
        @test Set(states(m)) == Set(vcat(species_in_group(:pts), [:M_lac__L_e]))
        @test issorted(species_index.(states(m)))

        u = SVector(0.0, 0.0011, 0.0161, 0.0009, 0.0131, 0.0012, 0.0141,
                    0.0062, 0.0349)
        p = [q.value for q in model_free_params(parameters(m))]
        u_inputs = SVector(0.0409, 1.46, 3.3660, 3.7076)   # pep, lac_c, pyr, g6p
        du = dynamics(u, p, 0.0, m, u_inputs)

        kf = (6600.0, 200000.0, 61000.0, 3900.0, 0.88)
        kr = (4000.0, 8000.0, 47000.0, 310.0, 1.00e-05)
        lac_e, ptsi, ptsi_P, ptsh, ptsh_P, crr, crr_P, ptsg, ptsg_P = u
        pep, lac_c, pyr, g6p = u_inputs
        v0 = kf[1] * ptsi * pep - kr[1] * ptsi_P * pyr
        v1 = kf[2] * ptsh * ptsi_P - kr[2] * ptsh_P * ptsi
        v2 = kf[3] * ptsh_P * crr - kr[3] * ptsh * crr_P
        v3 = kf[4] * ptsg * crr_P - kr[4] * crr * ptsg_P
        v4 = kf[5] * ptsg_P * 40.0 - kr[5] * ptsg * g6p
        v_export = 5.00e-09 * (lac_c - lac_e) * 3 / (200.0 * 1e-9)

        @test du[1] ≈ v_export / 1e5
        @test du[2] ≈ -v0 + v1
        @test du[3] ≈ v0 - v1
        @test du[4] ≈ -v1 + v2
        @test du[5] ≈ v1 - v2
        @test du[6] ≈ -v2 + v3
        @test du[7] ≈ v2 - v3
        @test du[8] ≈ -v3 + v4
        @test du[9] ≈ v3 - v4

        # The four foreign pools get terms through contributions, never through
        # this module's own derivative -- it returns nine entries, one per state
        # it owns, and none of them is a pool it does not integrate.
        @test length(du) == 9
        @test contributions(u, p, 0.0, m, u_inputs) ≈ SVector(-v0, -v_export, v0, v4)

        # 3P/r at the registry's radius is 0.075 s^-1, and 3/r is a sphere's
        # surface-to-volume ratio, so halving the radius doubles it.
        @test v_export ≈ 0.075 * (lac_c - lac_e)
        half = PtsTransport(radius_nm = 100.0)
        p_half = [q.value for q in model_free_params(parameters(half))]
        @test dynamics(u, p_half, 0.0, half, u_inputs)[1] ≈ 2 * du[1]

        # Zero at equal pools, and negative when external exceeds cytosolic.
        equal_pools = SVector(lac_c, ptsi, ptsi_P, ptsh, ptsh_P, crr, crr_P,
                              ptsg, ptsg_P)
        @test dynamics(equal_pools, p, 0.0, m, u_inputs)[1] == 0.0
        flooded = SVector(2 * lac_c, ptsi, ptsi_P, ptsh, ptsh_P, crr, crr_P,
                          ptsg, ptsg_P)
        @test dynamics(flooded, p, 0.0, m, u_inputs)[1] < 0

        # Mass action, not saturation: with pyruvate absent GLCpts0 is its
        # forward term alone, and doubling both its substrates quadruples it.
        q = InferCell._live_scalars(m, p)
        no_reverse = SVector(pep, lac_c, 0.0, g6p)
        base = InferCell._pts_rates(q, u, no_reverse)[1]
        @test base ≈ kf[1] * ptsi * pep
        @test InferCell._pts_rates(q, Base.setindex(u, 2 * ptsi, 2),
                                   SVector(2 * pep, lac_c, 0.0, g6p))[1] ≈ 4 * base

        # No cascade step touches an adenylate or phosphate species: the PTS
        # spends phosphoenolpyruvate, not ATP, which is why the net yield is
        # two ATP per glucose rather than three.
        touched = Set(vcat(states(m), inputs(m), [:M_glc__D_e]))
        for s in (:M_atp_c, :M_adp_c, :M_amp_c, :M_pi_c)
            @test !(s in touched)
        end
    end

    @testset "7.3 initial conditions reproduce the published copy numbers" begin
        factor = corea_particles_per_mM()
        ic = Dict(p.name => p for p in ic_params(parameters(m)))
        for (unphos, phos, copies) in PTS_CARRIERS
            total = ic[Symbol(unphos, "0")].value + ic[Symbol(phos, "0")].value
            @test round(Int, total * factor) == copies
            # Not merely close after rounding: the extract is built with the
            # same factor it is read back with, so the round trip is exact to
            # a part in 1e12 rather than to the nearest particle.
            @test total * factor ≈ copies rtol = 1e-12
        end
        @test ic[:M_lac__L_e0].value == 0.0

        # Both published inputs travel with the value: the copy number and the
        # proteomics fraction are in the extract, and the provenance names the
        # file they were read from rather than calling the split ours.
        rows = Dict(split(l, '\t')[1] => split(l, '\t')
                    for l in eachline(joinpath(COREA_DATA_DIR,
                                               "pts_initial_conditions.tsv"))
                    if startswith(l, "conc_"))
        @test rows["conc_M_ptsi_c"][4] == "353"
        @test rows["conc_M_ptsi_c"][5] == "0.05"
        @test occursin("protein_metabolites_frac.csv", rows["conc_M_ptsi_c"][6])
        @test rows["conc_M_ptsg_P_c"][4] == "831"
        @test rows["conc_M_ptsg_P_c"][5] == "0.85"
        @test occursin("membrane_protein_metabolites.csv",
                       rows["conc_M_ptsg_P_c"][6])
    end

    @testset "7.4 external lactate and the medium-to-cell volume ratio" begin
        # R = 1 is the shared-volume alternative: total lactate is then
        # conserved across the two pools and export saturates visibly.
        shared = PtsTransport(volume_ratio = 1.0)
        u = SVector(0.0, 0.0011, 0.0161, 0.0009, 0.0131, 0.0012, 0.0141,
                    0.0062, 0.0349)
        u_inputs = SVector(0.0409, 1.46, 3.3660, 3.7076)
        p_shared = [q.value for q in model_free_params(parameters(shared))]
        gain = dynamics(u, p_shared, 0.0, shared, u_inputs)[1]
        loss = -contributions(u, p_shared, 0.0, shared, u_inputs)[2]
        @test gain ≈ loss                      # nothing created, nothing lost
        p = [q.value for q in model_free_params(parameters(m))]
        @test dynamics(u, p, 0.0, m, u_inputs)[1] ≈ gain / 1e5

        # Saturation: at R = 1 the external pool reaches the cytosolic one and
        # the rate collapses; at the default it never gets near.
        equalised = Base.setindex(u, u_inputs[2], 1)
        @test dynamics(equalised, p_shared, 0.0, shared, u_inputs)[1] == 0.0

        @test extracellular_states(m) == [:M_lac__L_e]
        @test isempty(intersect(extracellular_states(m),
                                membrane_protein_states(m)))

        notes = reduction_declarations(m)
        ratio_notes = filter(l -> l.category === :lumping, notes)
        @test length(ratio_notes) == 1
        @test occursin("100000", only(ratio_notes).description)
        @test occursin("medium-to-cell", only(ratio_notes).description)
    end

    @testset "7.4 the one-percent criterion, over a full cycle" begin
        # The criterion behind the default ratio (spec §4 D6): external lactate
        # after a full 6,300 s cycle must stay below one percent of *steady*
        # cytosolic lactate. Steady cytosolic lactate needs a source, and the
        # source is glycolysis, which is phase 6 -- so the double supplies it at
        # the published rate, two lactate per glucose at ~1,106 glucose/s.
        #
        # Slow, and gated for it: this is the only 6,300 s integration in the
        # suite. Spec §9 asks whether full-cycle checks belong in the default
        # suite and says to decide on the measured cost; the measurement is in
        # the phase 7 handoff note.
        if get(ENV, "INFERCELL_FULL_CYCLE_TESTS", "false") == "true"
            production = 2 * 1106 / corea_particles_per_mM()   # mM/s
            @test production ≈ 0.10961 rtol = 1e-3
            source = [HeldMetabolites(lac_source = production), m]
            sol = integrate_pts(source; tspan = (0.0, 6300.0), saveat = 100.0)
            @test sol.retcode == InferCell.ReturnCode.Success

            i_c = _gidx(source, :M_lac__L_c)
            i_e = _gidx(source, :M_lac__L_e)
            steady = sol.u[end][i_c]
            external = sol.u[end][i_e]
            @test steady ≈ production / 0.075 rtol = 1e-3     # ~1.46 mM
            @test external / steady < 0.01                    # D6's criterion
            @test external ≈ 0.00689 rtol = 0.05              # D6's own figure

            # And the criterion is discriminating: a ratio of 1e4 misses it.
            missed = [HeldMetabolites(lac_source = production),
                      PtsTransport(volume_ratio = 1e4)]
            sol4 = integrate_pts(missed; tspan = (0.0, 6300.0), saveat = 100.0)
            @test sol4.u[end][_gidx(missed, :M_lac__L_e)] /
                  sol4.u[end][_gidx(missed, :M_lac__L_c)] > 0.01
        end
    end

    @testset "7.5 the boundary declaration" begin
        edges = coupling(m)
        kinds = [(e.species, edge_kind(e), e.direction) for e in edges]
        @test (:M_glc__D_e, :clamped, :in) in kinds
        @test (:M_pep_c, :mass, :in) in kinds
        @test (:M_lac__L_c, :mass, :in) in kinds
        @test (:M_pyr_c, :mass, :in) in kinds
        @test (:M_pyr_c, :mass, :out) in kinds
        @test (:M_g6p_c, :mass, :in) in kinds
        @test (:M_g6p_c, :mass, :out) in kinds
        @test (:M_lac__L_e, :volume, :in) in kinds
        @test (:M_ptsg_c, :volume, :out) in kinds
        @test (:M_ptsg_P_c, :volume, :out) in kinds
        @test length(edges) == 10
        @test all(e -> e.peer === nothing, edges)
        @test !any(e -> edge_kind(e) in (:currency, :deferred_counter,
                                         :rate_constant), edges)

        # Every foreign pool the rate law reads is an input, and external
        # glucose deliberately is not: the resolver refuses a chemostatted
        # input and prescribes the clamp instead.
        @test inputs(m) == [:M_pep_c, :M_lac__L_c, :M_pyr_c, :M_g6p_c]
        @test !(:M_glc__D_e in inputs(m))
        @test contributed_states(m) == inputs(m)

        # The transport reconstruction's 42.77 mM throws against the registry's
        # 40, which setICs_two.py:279 is the authority for.
        err = caught(() -> resolve_coupling(
            CoreAStub(:Rival; edges = [ClampedEdge(species = :M_glc__D_e,
                                                   direction = :in,
                                                   held_value = 42.77,
                                                   origin = :published)])))
        @test err isa ArgumentError
        @test occursin("42.77", sprint(showerror, err))

        # Standalone resolution reports rather than fails: this is a deliberate
        # partial composition, and the pools it draws from are simply unowned.
        graph = resolve_coupling(m)
        @test :M_pep_c in graph.unowned_states
        @test :M_lac__L_c in graph.unowned_states

        # The glucose clamp is the published model's, so it must not reach the
        # what-is-ours enumeration; the volume ratio must.
        labels = reduction_declarations(m)
        @test !any(l -> l.category === :clamp, labels)
        @test !any(l -> l.subject === :M_glc__D_e, labels)
        @test any(l -> l.category === :lumping, labels)

        # Drop an input the rate law needs, or add the chemostat, and the
        # resolver says which species and why.
        dropped = caught(() -> resolve_coupling(
            CoreAStub(:Dropped; st = collect(states(m)), edges = coupling(m),
                      ins = [:M_lac__L_c, :M_pyr_c, :M_g6p_c],
                      contribs = collect(inputs(m)),
                      membrane = [:M_ptsg_c, :M_ptsg_P_c])))
        @test dropped isa ArgumentError
        @test occursin("M_pep_c", sprint(showerror, dropped))

        added = caught(() -> resolve_coupling(
            CoreAStub(:Added; st = collect(states(m)), edges = coupling(m),
                      ins = vcat(collect(inputs(m)), [:M_glc__D_e]),
                      contribs = collect(inputs(m)),
                      membrane = [:M_ptsg_c, :M_ptsg_P_c])))
        @test added isa ArgumentError
        @test occursin("chemostat", sprint(showerror, added))
    end

    @testset "7.6 ptsG is flagged through the protocol function" begin
        # Found by sweeping the composition, without naming PtsTransport's type.
        @test membrane_protein_states(models) == [:M_ptsg_c, :M_ptsg_P_c]
        @test length(membrane_protein_states(models)) == 2
        for s in states(m)
            s in (:M_ptsg_c, :M_ptsg_P_c) && continue
            @test !(s in membrane_protein_states(models))
        end
        # Each flagged state carries the outbound edge the flag owes.
        outbound = Set(e.species for e in coupling(m)
                       if e isa VolumeEdge && e.direction === :out)
        @test outbound == Set([:M_ptsg_c, :M_ptsg_P_c])
        # The footprint is the published model's calibrated constant, and the
        # driver is what declares it -- this module supplies the states.
        @test InferCell.MEMBRANE_PROTEIN_FOOTPRINT_NM2 == 28.0
    end

    @testset "7.7 four carrier sums, four independent bounds" begin
        sol = integrate_pts(models)
        @test sol.retcode == InferCell.ReturnCode.Success

        # Non-negativity for all nine owned states at every save point.
        for s in states(m)
            i = _gidx(models, s)
            @test minimum(u[i] for u in sol.u) >= -1e-10
        end

        for (unphos, phos, _) in PTS_CARRIERS
            residual = carrier_residual(sol, models, unphos, phos)
            bound = carrier_bound(sol, models, unphos, phos;
                                  abstol = 1e-10, reltol = 1e-8)
            @test residual <= bound
        end

        # The tolerance principle: the bound is derived from the integrator, so
        # tightening the solver tenfold tightens the check, and the residual
        # stays under the tightened bound. A fixed threshold could be passed by
        # loosening it; this cannot.
        tight = integrate_pts(models; abstol = 1e-11, reltol = 1e-9)
        for (unphos, phos, _) in PTS_CARRIERS
            loose_bound = carrier_bound(sol, models, unphos, phos;
                                        abstol = 1e-10, reltol = 1e-8)
            tight_bound = carrier_bound(tight, models, unphos, phos;
                                        abstol = 1e-11, reltol = 1e-9)
            @test tight_bound <= loose_bound / 5
            @test carrier_residual(tight, models, unphos, phos) <= tight_bound
        end

        # The check can fail, and fails locally: one step written so a carrier
        # is created rather than transferred breaks that carrier's sum and
        # leaves the other three intact.
        mutant = [HeldMetabolites(), PtsCarrierLeak(m)]
        bad = integrate_pts(mutant; tspan = (0.0, 60.0), saveat = 6.0)
        broken = carrier_residual(bad, mutant, :M_ptsh_c, :M_ptsh_P_c)
        @test broken > carrier_bound(bad, mutant, :M_ptsh_c, :M_ptsh_P_c;
                                     abstol = 1e-10, reltol = 1e-8)
        for (unphos, phos, _) in PTS_CARRIERS
            unphos === :M_ptsh_c && continue
            @test carrier_residual(bad, mutant, unphos, phos) <=
                  carrier_bound(bad, mutant, unphos, phos;
                                abstol = 1e-10, reltol = 1e-8)
        end
    end

    @testset "7.x the export rationale is queryable, not only documented" begin
        rationale = lactate_export_rationale(m)
        @test occursin("691", rationale)
        @test occursin("345", rationale)
        @test occursin("error", rationale)
    end
end
