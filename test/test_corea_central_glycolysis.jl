using Test
using InferCell
using StaticArrays
using Distributions: LogNormal, params as dist_params
using SciMLBase: ReturnCode

import InferCell: states, parameters, dynamics, inputs, formalism

# Spec §11 phase 6 — central glycolysis, the ten reactions from
# glucose-6-phosphate through lactate.
#
# `caught` comes from `corea_test_models.jl`, included first by `runtests.jl`.
# The energy double is local to this file rather than shared: three sibling
# module branches want the same two-line double, and six branches editing one
# doubles file is the merge conflict the fan-out exists to avoid. Promoting it
# to `corea_test_models.jl` once they have all landed is a cheap change on main.

# ---------------------------------------------------------------------------
# The double
# ---------------------------------------------------------------------------

"""
    HeldEnergyPool()

Owns the three energy currencies glycolysis reads but does not integrate, and
the ten protein counts it declares as inputs, all with zero derivative.

It is not a model of anything. `inputs` resolves against states some module in
the composition integrates, so integrating glycolysis alone needs something that
owns them; holding them fixed is what makes the standalone trajectory
interpretable — and what the module's third reduction note says out loud.
"""
struct HeldEnergyPool <: AbstractSubModel
    params::Vector{InferParameter}
    st::Vector{Symbol}
end
function HeldEnergyPool()
    st = vcat(GLYCOLYSIS_CURRENCIES, default_protein_sources())
    vals = vcat([species_entry(s).initial_value for s in GLYCOLYSIS_CURRENCIES],
                [Float64(r.copies) for r in GLYCOLYTIC_REACTIONS])
    ps = [InferParameter(v, LogNormal(log(v), 1.0), true, Symbol(s, "0"),
                         :HeldEnergyPool, :initial_condition)
          for (s, v) in zip(st, vals)]
    return HeldEnergyPool(ps, st)
end
states(m::HeldEnergyPool) = m.st
parameters(m::HeldEnergyPool) = m.params
dynamics(u, p, t, ::HeldEnergyPool) = zero(u)

# ---------------------------------------------------------------------------
# Helpers written independently of the module
# ---------------------------------------------------------------------------

# Parameter value by name, straight off the InferParameter vector.
_pv(m, name::Symbol) = parameters(m)[findfirst(p -> p.name === name, parameters(m))].value

"""
    naive_rates(m, conc) -> Vector{Float64}

The ten net fluxes, recomputed from `GLYCOLYTIC_REACTIONS` and the module's own
imported constants with Dict arithmetic and no static vectors.

Written from the formula in spec §4 D1 rather than from `src`, so that agreement
with `reaction_rates` is two independent readings of the published rate law and
not one implementation checked against itself.
"""
function naive_rates(m, conc::Dict{Symbol, Float64};
                     reactions = GLYCOLYTIC_REACTIONS)
    ec = m.enzyme_conc
    out = Float64[]
    for (k, r) in enumerate(reactions)
        num_s = 1.0; den_s = 1.0
        for (s, n) in r.substrates
            x = conc[s] / _pv(m, Symbol("km_$(r.id)_$s"))
            num_s *= x^n; den_s *= (1 + x)^n
        end
        num_p = 1.0; den_p = 1.0
        for (s, n) in r.products
            x = conc[s] / _pv(m, Symbol("km_$(r.id)_$s"))
            num_p *= x^n; den_p *= (1 + x)^n
        end
        kf = _pv(m, Symbol("kcatF_$(r.id)"))
        kr = _pv(m, Symbol("kcatR_$(r.id)"))
        push!(out, ec[k] * (kf * num_s - kr * num_p) / (den_s + den_p - 1))
    end
    return out
end

# The registry's initial concentrations for everything a rate law reads.
function registry_concentrations()
    d = Dict{Symbol, Float64}()
    for s in vcat(glycolytic_states(), GLYCOLYSIS_CURRENCIES)
        d[s] = species_entry(s).initial_value
    end
    return d
end

_u0(conc) = SVector{13, Float64}(conc[s] for s in glycolytic_states())
_uin(conc) = SVector{13, Float64}(vcat([conc[s] for s in GLYCOLYSIS_CURRENCIES],
                                       [Float64(r.copies) for r in GLYCOLYTIC_REACTIONS]))

"""
    redox_drift(sol, nad, nadh) -> Float64

The largest departure of NAD⁺ + NADH from its initial value over the saved
points. GAPD and LDH_L are the only reactions that touch the pair and their
stoichiometry mirrors, so the sum is invariant and any drift is the integrator's
— or a broken stoichiometry's.
"""
function redox_drift(sol, nad::Int, nadh::Int)
    total0 = sol.u[1][nad] + sol.u[1][nadh]
    return maximum(abs(u[nad] + u[nadh] - total0) for u in sol.u)
end

"""
    assert_redox_conserved(sol, nad, nadh; bound)

Throw naming the conserved pool and the drift where NAD⁺ + NADH has moved
further than the integrator can account for. The message is the assertion: a
mutation test that could only report `false` would not show *which* invariant
broke.
"""
function assert_redox_conserved(sol, nad::Int, nadh::Int; bound::Float64)
    drift = redox_drift(sol, nad, nadh)
    drift <= bound || throw(ArgumentError(
        "The conserved pool M_nad_c + M_nadh_c is not invariant: it drifts by " *
        "$drift mM over the trajectory, against an integrator bound of $bound mM. " *
        "GAPD and LDH_L are the only reactions that touch it and their " *
        "stoichiometry must mirror"))
    return drift
end

# tol_C of spec §3: the single-run bound on a conserved sum, from the
# integrator's own tolerances. Reported as a secondary number; the assertion is
# the scaling, not this.
function conservation_bound(sol, idxs, abstol, reltol)
    return sum(max(abstol, reltol * maximum(abs(u[i]) for u in sol.u)) for i in idxs)
end

# One reaction's stoichiometry mutated, everything else identical.
function mutate_stoichiometry(id::Symbol, species::Symbol, coefficient::Int)
    return [r.id === id ?
            (id = r.id, enzyme = r.enzyme, locus = r.locus, copies = r.copies,
             substrates = [s === species ? (s => coefficient) : (s => n)
                           for (s, n) in r.substrates],
             products = [s === species ? (s => coefficient) : (s => n)
                         for (s, n) in r.products]) : r
            for r in GLYCOLYTIC_REACTIONS]
end

const GLYCOLYSIS_ABSTOL = 1e-10
const GLYCOLYSIS_RELTOL = 1e-8
const GLYCOLYSIS_CYCLE = 6300.0

@testset "Phase 6 — central glycolysis" begin
    m = CentralGlycolysis()
    owned = glycolytic_states()
    nad = findfirst(==(:M_nad_c), owned)
    nadh = findfirst(==(:M_nadh_c), owned)

    @testset "6.1 the vendored extract" begin
        rows = [split(l, '\t') for l in readlines(CENTRAL_GLYCOLYSIS_TABLE)
                if !startswith(l, "!") && !startswith(l, "%") && !isempty(strip(l))]
        ids = [String(r[1]) for r in rows]
        @test length(rows) == 65
        @test count(startswith("kcatF_"), ids) == 10
        @test count(startswith("kcatR_"), ids) == 10
        @test count(startswith("km_"), ids) == 32
        @test count(startswith("conc_"), ids) == 13
        @test length(unique(ids)) == 65

        # Byte for byte against the upstream file, not merely numerically equal:
        # a reshape that reprinted 650 as 650.0 would still pass an isapprox.
        cell(id, col) = String(rows[findfirst(==(id), ids)][col])
        @test cell("kcatF_R_PGI", 2) == "804.3384"
        @test cell("kcatF_R_PGI", 3) == "1.0513"
        @test cell("kcatF_R_FBA", 2) == "59.7"
        @test cell("kcatF_R_FBA", 3) == "1.7466"
        @test cell("km_R_PGI_M_g6p_c", 2) == "22.9419"
        # The archived reference names conc_M_atp_c here while also fixing the
        # count at thirteen; ATP is not an owned state, so the two cannot both
        # hold. src/organisms/coreA/data/README.md records the substitution.
        @test cell("conc_M_g6p_c", 2) == "3.7076"
        @test cell("conc_M_g6p_c", 3) == "1.2785"
        @test cell("kcatR_R_PGI", 2) == "650"

        # Every row carries the triple it came from, so a renamed identifier
        # stays traceable.
        @test all(r -> count(==('|'), String(r[4])) == 2, rows)
    end

    @testset "6.2 the states and the rate law" begin
        @test length(states(m)) == 13
        @test Set(states(m)) ==
              union(Set(species_in_group(:glycolytic)), Set(species_in_group(:redox)))
        idx = species_index.(states(m))
        @test all(idx[i] < idx[i + 1] for i in 1:12)
        @test formalism(m) === :ode

        # One rate-law term per unit of stoichiometry, and the total is the
        # number of Michaelis constants the extract holds.
        terms = [(substrate_terms(r), product_terms(r)) for r in m.rates]
        @test terms == [(1, 1), (2, 2), (1, 2), (1, 1), (3, 2),
                        (2, 2), (1, 1), (1, 1), (2, 2), (2, 2)]
        @test sum(sum, terms) == 32

        conc = registry_concentrations()
        v = reaction_rates(_u0(conc), SVector{0, Float64}(), 0.0, m, _uin(conc))
        @test v ≈ naive_rates(m, conc) rtol = 1e-12
        @test all(isfinite, v)

        # Equal forward and reverse terms give exactly zero net rate. PGI has
        # one substrate and one product, so the balancing F6P follows in closed
        # form: kcatF·(G6P/KmS) = kcatR·(P/KmP).
        balanced = copy(conc)
        balanced[:M_f6p_c] = conc[:M_g6p_c] * _pv(m, :kcatF_R_PGI) *
                             _pv(m, :km_R_PGI_M_f6p_c) /
                             (_pv(m, :kcatR_R_PGI) * _pv(m, :km_R_PGI_M_g6p_c))
        vb = reaction_rates(_u0(balanced), SVector{0, Float64}(), 0.0, m, _uin(balanced))
        forward = _pv(m, :kcatF_R_PGI) * conc[:M_g6p_c] / _pv(m, :km_R_PGI_M_g6p_c)
        @test abs(vb[1]) < 1e-12 * forward * m.enzyme_conc[1]

        # Products far above substrates run a reaction backwards. Done one
        # reaction at a time: scaling every product at once would also scale
        # some other reaction's substrates — F6P is PGI's product and PFK's
        # substrate — and the factors would cancel rather than reverse it.
        for (k, r) in enumerate(GLYCOLYTIC_REACTIONS)
            swamped = copy(conc)
            for (s, _) in r.products
                swamped[s] = conc[s] * 1e8
            end
            vk = reaction_rates(_u0(swamped), SVector{0, Float64}(), 0.0, m, _uin(swamped))
            @test vk[k] < 0
        end

        # The derivative is the stoichiometry times the fluxes and nothing else:
        # hand-checked at the registry's initial conditions on the one state only
        # two reactions touch.
        du = dynamics(_u0(conc), SVector{0, Float64}(), 0.0, m, _uin(conc))
        @test length(du) == 13
        nv = naive_rates(m, conc)
        @test du[findfirst(==(:M_2pg_c), owned)] ≈ nv[7] - nv[8] rtol = 1e-12   # PGM in, ENO out
        @test du[nad] ≈ nv[10] - nv[5] rtol = 1e-12                             # LDH_L in, GAPD out
        @test du[nadh] ≈ nv[5] - nv[10] rtol = 1e-12
        @test du[nad] ≈ -du[nadh] rtol = 1e-14
    end

    @testset "6.3 every value arrives through the loader" begin
        ps = parameters(m)
        @test length(ps) == 65
        @test all(p -> source_file(p) == "central_balanced", ps)
        @test all(p -> provenance_of(p).identifier !== nothing, ps)
        @test length(rate_params(ps)) == 52
        @test length(ic_params(ps)) == 13

        # One file, so nothing was chosen over an alternative yet. The
        # nucleotide file's rival values arrive with phase 8.
        @test isempty(governing_choices(ps))

        # Fixed by default: composing the module does not silently produce a
        # 65-dimensional posterior.
        @test all(p -> p.fixed, ps)
        @test isempty(model_free_params(ps))

        free = CentralGlycolysis(; free = [:kcatF_R_FBA])
        @test length(model_free_params(parameters(free))) == 1
        @test only(model_free_params(parameters(free))).name === :kcatF_R_FBA
        # And it reaches the rate law from `p` rather than from the struct.
        conc = registry_concentrations()
        base = reaction_rates(_u0(conc), SVector{0, Float64}(), 0.0, m, _uin(conc))
        doubled = reaction_rates(_u0(conc), SA[2 * _pv(m, :kcatF_R_FBA)], 0.0,
                                 free, _uin(conc))
        @test doubled[3] > base[3]
        @test doubled[1] ≈ base[1] rtol = 1e-14

        # The prior is the balancing distribution's own shape.
        pgi = ps[findfirst(p -> p.name === :kcatF_R_PGI, ps)]
        @test dist_params(pgi.prior) == (log(804.3384), log(1.0513))
        @test informedness(pgi) === :balanced

        # The two states nothing informed keep the prior default's width.
        for s in (:M_g3p_c, :M_lac__L_c)
            p = ps[findfirst(q -> q.name === Symbol(s, "0"), ps)]
            @test informedness(p) === :prior_default
            @test dist_params(p.prior)[2] ≈ log(10.0)
            @test species_entry(s).informedness === :prior_default
        end
        # Four parameters carry no informed width, not two: the two
        # concentrations above, and the reverse constants of PFK and PYK, whose
        # geometric standard deviations are 33.03 and 11.34 — at or above the
        # prior default. Naming them is the point; a bare count of two would
        # have been wrong and would have looked right.
        @test Set(p.name for p in uninformed_params(ps)) ==
              Set([:M_g3p_c0, :M_lac__L_c0, :kcatR_R_PFK, :kcatR_R_PYK])

        # The registry-agreement check is live: mutating one concentration row
        # in a temporary copy makes the load throw naming the species, the
        # imported value and the registry's.
        mktempdir() do dir
            path = joinpath(dir, "central_glycolysis.tsv")
            text = read(CENTRAL_GLYCOLYSIS_TABLE, String)
            write(path, replace(text, "conc_M_pyr_c\t3.366\t" => "conc_M_pyr_c\t9.999\t"))
            err = caught(() -> CentralGlycolysis(; table_path = path))
            @test err isa ArgumentError
            msg = sprint(showerror, err)
            @test occursin("M_pyr_c", msg)
            @test occursin("9.999", msg)
            @test occursin("3.366", msg)
        end
    end

    @testset "6.4 enzyme concentrations, nominal from copy number" begin
        factor = corea_particles_per_mM()
        copies = [r.copies for r in GLYCOLYTIC_REACTIONS]
        @test copies == [266, 458, 775, 410, 1355, 411, 322, 998, 551, 1100]
        @test m.enzyme_conc == SVector{10, Float64}(c / factor for c in copies)

        # The scoping note's rounded 20,180 agrees to five decimals. It differs
        # in the sixth for PFK, GAPD, PGK, ENO and LDH_L, because the factor is
        # derived from Avogadro at the registry's radius (20,180.39) rather than
        # transcribed — which is what makes it agree with the conversion the
        # handshake performs.
        @test all(isapprox(m.enzyme_conc[i], copies[i] / 20180; atol = 1e-5) for i in 1:10)
        @test factor ≈ 20180.39 rtol = 1e-6

        # Not the published no-rule default. All ten reactions have single-gene
        # rules, so 0.001 mM would understate every flux by 13× to 67×.
        @test all(!=(0.001), m.enzyme_conc)
        @test minimum(m.enzyme_conc) > 0.013

        # Overriding one scales exactly the rates that enzyme catalyses.
        conc = registry_concentrations()
        base = reaction_rates(_u0(conc), SVector{0, Float64}(), 0.0, m, _uin(conc))
        ec = collect(m.enzyme_conc)
        ec[8] *= 3                                     # ENO
        m3 = CentralGlycolysis(; enzyme_conc = ec)
        scaled = reaction_rates(_u0(conc), SVector{0, Float64}(), 0.0, m3, _uin(conc))
        @test scaled[8] ≈ 3 * base[8] rtol = 1e-14
        @test all(scaled[i] ≈ base[i] for i in 1:10 if i != 8)

        # The ten protein counts are declared so translation can supersede them.
        @test inputs(m)[1:3] == GLYCOLYSIS_CURRENCIES
        @test inputs(m)[4:13] == default_protein_sources()
        @test length(default_protein_sources()) == 10
        @test !any(is_registered, default_protein_sources())
        renamed = CentralGlycolysis(; protein_sources = [Symbol(:X, i) for i in 1:10])
        @test inputs(renamed)[4] === :X1
    end

    @testset "6.5 the boundary" begin
        edges = coupling(m)
        @test length(edges) == 9
        @test count(e -> e isa CurrencyEdge, edges) == 5
        @test count(e -> e isa MassEdge, edges) == 4
        @test all(e -> e.peer === nothing, edges)

        kinds = [(e.species, edge_kind(e), e.direction) for e in edges]
        @test kinds == [(:M_atp_c, :currency, :in), (:M_atp_c, :currency, :out),
                        (:M_adp_c, :currency, :in), (:M_adp_c, :currency, :out),
                        (:M_pi_c, :currency, :in),
                        (:M_g6p_c, :mass, :in), (:M_pyr_c, :mass, :in),
                        (:M_pep_c, :mass, :out), (:M_lac__L_c, :mass, :out)]

        # No edge names a redox species. NAD+/NADH are touched only by GAPD and
        # LDH_L, both inside this module, which is why check 3 is per-module.
        @test !any(e -> e.species in (:M_nad_c, :M_nadh_c), edges)
        @test all(e -> is_registered(e.species), edges)

        @test contributed_states(m) == GLYCOLYSIS_CURRENCIES

        # Standalone resolution reports rather than fails, and reports the three
        # currencies among the states no module here integrates.
        graph = resolve_coupling([m])
        @test length(graph.edges) == 9
        @test issubset(GLYCOLYSIS_CURRENCIES, graph.unowned_states)
        @test length(graph.unowned_states) == n_dynamic_states() - 13
        @test !any(s -> s in graph.unowned_states, states(m))

        # The outbound currency edges execute (spec §11 phase 6's amendment
        # from D0): each moves exactly the mass its rate law says. ATP is drawn
        # by PFK and supplied by PGK and PYK, and the net is what reaches the
        # owner's derivative.
        conc = registry_concentrations()
        c = contributions(_u0(conc), SVector{0, Float64}(), 0.0, m, _uin(conc))
        nv = naive_rates(m, conc)
        @test length(c) == 3
        @test c[1] ≈ nv[6] + nv[9] - nv[2] rtol = 1e-12    # ATP: PGK + PYK − PFK
        @test c[2] ≈ nv[2] - nv[6] - nv[9] rtol = 1e-12    # ADP: the mirror
        @test c[3] ≈ -nv[5] rtol = 1e-12                   # Pi: GAPD draws
        @test c[1] ≈ -c[2] rtol = 1e-14

        # A composition that owns the currencies wires the contribution into
        # their derivative, so the pool sees glycolysis draw and supply.
        models = [m, HeldEnergyPool()]
        prob = build_problem(models; tspan = (0.0, 1.0))
        du = prob.f(prob.u0, prob.p, 0.0)
        atp = 13 + 1
        @test du[atp] ≈ c[1] rtol = 1e-12
        @test du[13 + 3] ≈ c[3] rtol = 1e-12
        # The double integrates nothing of its own, so what the pool's
        # derivative holds is the contribution and only the contribution.
        @test du[13 + 4] == 0.0
    end

    @testset "6.6 what is ours" begin
        labels = reduction_declarations([m, HeldEnergyPool()])
        ours = filter(l -> l.subject === :CentralGlycolysis, labels)
        @test length(ours) == 4
        text = join((l.description for l in ours), " ")

        @test occursin("NOX", text)
        @test occursin("LDH_L", text) && occursin("gene-protein-reaction", text) &&
              occursin("prior-dominated", text)
        @test !any(r -> r.id === :R_NOX, GLYCOLYTIC_REACTIONS)

        # All eight differing Michaelis constants, with both values.
        @test length(KM_COLUMN_DISAGREEMENTS) == 8
        for (id, ran, imported) in KM_COLUMN_DISAGREEMENTS
            @test occursin(id, text)
            @test occursin(string(ran), text)
            @test occursin(string(imported), text)
            # The imported half is what the module actually runs on.
            @test _pv(m, Symbol(id)) ≈ imported rtol = 1e-12
        end
        @test occursin("227", text)     # the FBA G3P constant's factor

        @test occursin("M_atp_c", text) && occursin("not a closed energy loop", text)
        @test occursin("nominal", text)

        # Each reads as a sentence rather than a slug.
        @test all(l -> length(split(l.description)) > 8, ours)
        @test all(l -> l.category === :lumping, ours)
    end

    @testset "6.7 integration and redox balance" begin
        models = [m, HeldEnergyPool()]
        prob = build_problem(models; tspan = (0.0, GLYCOLYSIS_CYCLE))
        saveat = 0.0:60.0:GLYCOLYSIS_CYCLE

        sol = solve(prob, Rodas5P(); abstol = GLYCOLYSIS_ABSTOL,
                    reltol = GLYCOLYSIS_RELTOL, saveat = saveat)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.u) == length(saveat)

        # Non-negativity, each state against its own integrator bound rather
        # than against zero: the solver is entitled to its tolerance. One
        # assertion over the whole trajectory, naming the state and time that
        # came closest, rather than 1,378 assertions that would report only
        # that one of them failed.
        bounds = [max(GLYCOLYSIS_ABSTOL,
                      GLYCOLYSIS_RELTOL * maximum(abs(u[i]) for u in sol.u)) for i in 1:13]
        mins = [minimum(u[i] for u in sol.u) for i in 1:13]
        margins = mins .+ bounds
        @test all(>(0), margins)
        @test all(>=(0), mins)          # in fact none of the thirteen goes negative at all
        tightest = argmin(margins)
        @info "phase 6 non-negativity" state=owned[tightest] minimum=mins[tightest] bound=bounds[tightest]

        # The redox pair is conserved, and this is the invariant spec §3's
        # tolerance principle carves its one exception for (amended 2026-09-10;
        # see §12). GAPD and LDH_L are the only reactions that touch the pair,
        # both live here, and their stoichiometry mirrors — so the composed
        # right-hand side returns du[NAD+] and du[NADH] as bit-for-bit
        # negatives and the sum is conserved whatever step the solver takes.
        # There is no integrator error in it to shrink. What is asserted is
        # therefore **flatness at the floating-point floor** across a ladder of
        # tolerances, which is a stronger statement than a fivefold fall and
        # not a weaker one.
        ladder = [(1e-4, 1e-2), (1e-8, 1e-6),
                  (GLYCOLYSIS_ABSTOL, GLYCOLYSIS_RELTOL), (1e-12, 1e-10)]
        drifts = Float64[]
        for (a, r) in ladder
            l = (a, r) == (GLYCOLYSIS_ABSTOL, GLYCOLYSIS_RELTOL) ? sol :
                solve(prob, Rodas5P(); abstol = a, reltol = r, saveat = saveat)
            @test l.retcode == ReturnCode.Success
            push!(drifts, redox_drift(l, nad, nadh))
        end

        drift_loose = drifts[end - 1]                       # at the pinned pair
        total0 = sol.u[1][nad] + sol.u[1][nadh]
        floor_ulps = drifts ./ eps(total0)
        tol_C = conservation_bound(sol, (nad, nadh), GLYCOLYSIS_ABSTOL, GLYCOLYSIS_RELTOL)

        # Every point on the ladder is a few tens of ulps of the conserved sum,
        # and eight decades of tolerance move it by less than one decade.
        @test all(<(100), floor_ulps)
        @test maximum(drifts) / minimum(drifts) < 100
        @test drift_loose <= tol_C

        @info "phase 6 redox residual, flat at the floating-point floor" ladder drifts floor_ulps eps_of_sum=eps(total0) tol_C

        # The check can fail, and says what broke. GAPD's NAD+ stoichiometry no
        # longer mirrors LDH_L's, so the sum is no longer invariant.
        mutated = CentralGlycolysis(; reactions = mutate_stoichiometry(:R_GAPD, :M_nad_c, 2))
        mprob = build_problem([mutated, HeldEnergyPool()];
                              tspan = (0.0, GLYCOLYSIS_CYCLE))
        msol = solve(mprob, Rodas5P(); abstol = GLYCOLYSIS_ABSTOL,
                     reltol = GLYCOLYSIS_RELTOL, saveat = saveat)
        err = caught(() -> assert_redox_conserved(msol, nad, nadh; bound = tol_C))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("M_nad_c + M_nadh_c", msg)
        @test occursin("drifts by", msg)
        mutated_drift = redox_drift(msol, nad, nadh)
        @test mutated_drift > 1e-6                      # millimolar, not roundoff
        @test mutated_drift > 1e6 * max(drift_loose, eps())
        # …and the leak does not shrink when the solver is tightened either,
        # which is the half of the tolerance principle that still bites: a
        # structural error is invisible to the integrator's tolerance.
        mtight = solve(mprob, Rodas5P(); abstol = GLYCOLYSIS_ABSTOL / 100,
                       reltol = GLYCOLYSIS_RELTOL / 100, saveat = saveat)
        @test redox_drift(mtight, nad, nadh) > 0.2 * mutated_drift
        @info "phase 6 redox mutation" mutated_drift orders=log10(mutated_drift / drift_loose)

        # …and the unmutated check is not vacuously true of a trajectory that
        # never moved: the states this module owns do change over the cycle.
        @test assert_redox_conserved(sol, nad, nadh; bound = tol_C) == drift_loose
        moved = [maximum(abs(u[i] - sol.u[1][i]) for u in sol.u) for i in 1:13]
        @test count(>(1e-6), moved) >= 10
    end

    @testset "Done when: thirteen states integrate over a cycle, redox conserved" begin
        models = [m, HeldEnergyPool()]
        prob = build_problem(models; tspan = (0.0, GLYCOLYSIS_CYCLE))
        sol = solve(prob, Rodas5P(); abstol = GLYCOLYSIS_ABSTOL,
                    reltol = GLYCOLYSIS_RELTOL, saveat = 0.0:60.0:GLYCOLYSIS_CYCLE)
        drift = redox_drift(sol, nad, nadh)
        lac = findfirst(==(:M_lac__L_c), owned)
        g6p = findfirst(==(:M_g6p_c), owned)

        @test sol.retcode == ReturnCode.Success
        @test all(u -> all(>=(-1e-8), u[1:13]), sol.u)
        @test drift < 1e-8
        # Standalone this is not physiological and is not asserted to be: G6P is
        # never replenished, so the module drains its own carbon into lactate.
        # Carbon balance spans transport and export and belongs to phase 14.
        @test sol.u[end][lac] > sol.u[1][lac]
        @test sol.u[end][g6p] < sol.u[1][g6p]

        @info "phase 6 done-when" states=13 reactions=length(m.rates) params=length(parameters(m)) edges=length(coupling(m)) redox_drift=drift lactate=sol.u[end][lac] g6p=sol.u[end][g6p]
    end
end
