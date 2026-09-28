# Spec §11 task 14c.1: prior draws keep each reaction's equilibrium constant.
#
# Verify by every draw reproducing each reaction's published equilibrium
# constant to roundoff, by a mutation that draws the reverse constant
# independently failing, and by a draw of D11's `kcat_ENO` or `kcat_FBA` alone
# scaling that reaction's reverse constant by the same factor. "Published" is
# the constant the nominal Mode values imply (§12, 2026-09-28, planning A).

using Test
using InferCell
using Random

# Each reaction's species, substrates then products, one entry per unit of
# stoichiometry, in the order the rate law reads them.
const _GLY_SPECIES = Dict(r.id => ([s for (s, n) in r.substrates for _ in 1:n],
                                   [s for (s, n) in r.products for _ in 1:n])
                          for r in GLYCOLYTIC_REACTIONS)
const _REC_SPECIES = Dict(
    :R_PGK3 => ([:M_13dpg_c, :M_gdp_c], [:M_3pg_c, :M_gtp_c]),
    :R_PYK3 => ([:M_gdp_c, :M_pep_c], [:M_gtp_c, :M_pyr_c]),
    :R_ADK1 => ([:M_amp_c, :M_atp_c], [:M_adp_c, :M_adp_c]),
    :R_GK1 => ([:M_atp_c, :M_gmp_c], [:M_adp_c, :M_gdp_c]),
    :R_PPA => ([:M_ppi_c], [:M_pi_c, :M_pi_c]))

# Every species at 1 mM except the reaction's first product, set so the
# mass-action ratio Π P / Π S equals `keq`. If the relation describes the rate
# law, that reaction's flux vanishes there.
function _at_equilibrium(species, keq)
    subs, prods = species
    c = Dict{Symbol, Float64}()
    first_p = prods[1]
    n = count(==(first_p), prods)
    c[first_p] = keq^(1 / n)
    return s -> get(c, s, 1.0)
end

function _glycolysis_flux(m, rel, keq)
    conc = _at_equilibrium(_GLY_SPECIES[rel.reaction], keq)
    u = [conc(s) for s in glycolytic_states()]
    ui = [conc(s) for s in (:M_atp_c, :M_adp_c, :M_pi_c)]
    k = findfirst(r -> r.id === rel.reaction, GLYCOLYTIC_REACTIONS)
    return reaction_rates(u, Float64[], 0.0, m, ui)[k], k
end

function _recycling_flux(m, rel, keq)
    conc = _at_equilibrium(_REC_SPECIES[rel.reaction], keq)
    u = [conc(s) for s in RECYCLING_STATES]
    ui = [conc(s) for s in (:M_13dpg_c, :M_pep_c, :M_3pg_c, :M_pyr_c)]
    k = findfirst(==(rel.reaction), RECYCLING_REACTIONS)
    return recycling_fluxes(u, Float64[], 0.0, m, ui)[k], k
end

@testset "Haldane-consistent prior draws (spec §11 task 14c.1)" begin
    gly = CentralGlycolysis()
    rec = NucleotideRecycling()
    nominal = nominal_parameter_values(AbstractSubModel[gly, rec])

    @testset "the relations describe the rate laws that run" begin
        rels = haldane_relations(AbstractSubModel[gly, rec])
        @test length(rels) == 15
        @test [r.reaction for r in rels] ==
              vcat([r.id for r in GLYCOLYTIC_REACTIONS], collect(RECYCLING_REACTIONS))
        # Setting the mass-action ratio to the relation's constant zeroes the
        # flux: the relation is the rate law's, not a transcription of it.
        for rel in haldane_relations(gly)
            keq = equilibrium_constant(rel, nominal)
            v, _ = _glycolysis_flux(gly, rel, keq)
            v2, _ = _glycolysis_flux(gly, rel, 2keq)
            @test abs(v) <= 1e-12 * abs(v2)
        end
        for rel in haldane_relations(rec)
            keq = equilibrium_constant(rel, nominal)
            v, _ = _recycling_flux(rec, rel, keq)
            v2, _ = _recycling_flux(rec, rel, 2keq)
            @test abs(v) <= 1e-12 * abs(v2)
        end
        # And the tie can fail: ADK1 with its squared product written as a
        # single power puts the equilibrium somewhere the rate law does not.
        adk = only(r for r in haldane_relations(rec) if r.reaction === :R_ADK1)
        wrong = HaldaneRelation(adk.reaction, adk.forward, adk.reverse,
                                adk.substrates, [first(only(adk.products)) => 1])
        v, _ = _recycling_flux(rec, adk, equilibrium_constant(wrong, nominal))
        v2, _ = _recycling_flux(rec, adk, 2equilibrium_constant(wrong, nominal))
        @test abs(v) > 1e-3 * abs(v2)
        # A disabled reaction has no relation.
        @test length(haldane_relations(NucleotideRecycling(
            reactions = filter(!=(:R_ADK1), RECYCLING_REACTIONS)))) == 4
    end

    @testset "the constants kept are the Mode-implied ones (§12 planning A)" begin
        keq = Dict(r.reaction => equilibrium_constant(r, nominal)
                   for r in haldane_relations(AbstractSubModel[gly, rec]))
        # Upstream's equilibrium-constant rows are 3.0982 and 0.0012; the
        # Modes imply these.
        @test keq[:R_ENO] ≈ 3.0976 rtol = 1e-4
        @test keq[:R_TPI] ≈ 18.858 rtol = 1e-4
        @test keq[:R_PGK] ≈ 7.545e5 rtol = 1e-3
    end

    # The broad census's set: every informed kinetic constant of the two
    # balanced modules (spec §11 phase 14b, decided 2026-09-26).
    informed(m) = Symbol[q.name for q in parameters(m)
                         if q.role === :rate && informedness(q) === :balanced]
    broad = vcat(informed(gly), informed(rec))
    ms = AbstractSubModel[gly, rec]

    @testset "every derived draw keeps every equilibrium constant" begin
        worst = 0.0
        for seed in 1:200
            draw = draw_parameters(MersenneTwister(seed), ms, broad)
            rows = assert_haldane(ms, draw)
            @test length(rows) == 15
            worst = max(worst, maximum(r.relative / r.bound for r in rows))
            # Every reverse constant is derived, none drawn: the draw names each
            # reaction's reverse constant exactly once.
            revs = [r.reverse for r in haldane_relations(ms)]
            @test count(p -> first(p) in revs, draw) == 15
        end
        @test worst <= 1
        @info "14c.1: worst Haldane residual over 200 broad draws, as a fraction of its roundoff bound" worst
    end

    @testset "the mutation: independent reverse draws fail, naming a reaction" begin
        draw = draw_parameters(MersenneTwister(1), ms, broad; reverse = :independent)
        err = try
            assert_haldane(ms, draw); nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("equilibrium constant", err.msg)
        @test occursin("Haldane", err.msg)
        # And it misses by orders, not by a hair.
        vals = merge(nominal, Dict(draw))
        misses = [abs(equilibrium_constant(r, vals) / equilibrium_constant(r, nominal) - 1)
                  for r in haldane_relations(ms)]
        @test maximum(misses) > 1e6 * 4 * 12 * eps()
    end

    @testset "a forward constant drawn alone scales its reverse by the same factor" begin
        for (f, r) in ((:kcatF_R_ENO, :kcatR_R_ENO), (:kcatF_R_FBA, :kcatR_R_FBA))
            for seed in 1:20
                draw = Dict(draw_parameters(MersenneTwister(seed), ms, [f]))
                @test keys(draw) == Set([f, r])
                @test draw[r] / nominal[r] ≈ draw[f] / nominal[f] rtol = 8eps()
                @test draw[f] != nominal[f]
            end
        end
        # Naming the reverse constant too changes nothing: it is derived, not drawn.
        a = draw_parameters(MersenneTwister(3), ms, [:kcatF_R_ENO])
        b = draw_parameters(MersenneTwister(3), ms, [:kcatF_R_ENO, :kcatR_R_ENO])
        @test Dict(a) == Dict(b)
    end

    @testset "argument checks" begin
        @test_throws ArgumentError draw_parameters(MersenneTwister(1), ms, [:no_such])
        @test_throws ArgumentError draw_parameters(MersenneTwister(1), ms,
                                                   [:kcatF_R_ENO]; reverse = :other)
        @test_throws ArgumentError draw_parameters(MersenneTwister(1), ms,
                                                   [:kcatF_R_ENO, :kcatF_R_ENO])
    end

    @testset "D11's six targets, drawn and written into a driver" begin
        dms = d11_models()
        @test length(D11_TARGETS) == 6
        draw = draw_parameters(MersenneTwister(14_900), dms, D11_TARGETS)
        @test Set(first.(draw)) == Set(vcat(D11_TARGETS, [:kcatR_R_ENO, :kcatR_R_FBA]))
        assert_haldane(dms, draw)

        d = build_problem(dms; tspan = (0.0, 120.0), complete = true)
        ode_names = InferCell._block_param_names(filter(m -> formalism(m) === :ode, dms))
        jump_names = InferCell._block_param_names(filter(m -> formalism(m) === :jump, dms))
        # The layout `set_parameters!` assumes is the one the builder made: every
        # slot the driver does not overwrite holds its name's nominal value. Most
        # values are distinct, so a permuted layout would fail here.
        nominal = nominal_parameter_values(dms)
        written = Set(w.param_slot for w in driver_written_params(d))
        for (names, p) in ((ode_names, d.ode.p), (jump_names, d.jump.p))
            @test length(names) == length(p)
            @test all(p[i] == nominal[n] for (i, n) in enumerate(names) if !(n in written))
        end
        @test length(unique(d.ode.p)) > length(d.ode.p) ÷ 2
        k_gapd = findfirst(==(Symbol("k_tx_JCVISYN3A_0607")), jump_names)
        s_gapd = findfirst(==(Symbol("S_JCVISYN3A_0607")), jump_names)
        p_ode0 = copy(d.ode.p)
        p_jump0 = copy(d.jump.p)

        # A nominal draw leaves the driver exactly as built.
        set_parameters!(d, dms, [n => nominal_parameter_values(dms)[n]
                                 for n in first.(draw)])
        @test d.ode.p == p_ode0 && d.jump.p == p_jump0

        set_parameters!(d, dms, draw)
        dd = Dict(draw)
        for (n, v) in draw
            i = findfirst(==(n), ode_names)
            i === nothing ? (@test d.jump.p[findfirst(==(n), jump_names)] == v) :
                            (@test d.ode.p[i] == v)
        end
        # The rebuilt rate constant follows the promoter from t = 0, and the
        # genes whose promoters did not move keep theirs.
        @test d.jump.p[k_gapd] / p_jump0[k_gapd] ≈ dd[jump_names[s_gapd]] / p_jump0[s_gapd] rtol = 1e-12
        untouched = [i for (i, n) in enumerate(jump_names)
                     if startswith(String(n), "k_tx_") &&
                        !(Symbol("S_", chopprefix(String(n), "k_tx_")) in keys(dd))]
        @test length(untouched) == 14
        @test d.jump.p[untouched] ≈ p_jump0[untouched] rtol = 1e-12

        @test_throws ArgumentError set_parameters!(d, dms, [:enz_R_ENO => 1.0])
        @test_throws ArgumentError set_parameters!(d, dms, [:kcatF_R_PGI => 1.0])
        handshake_step!(d)
        @test_throws ArgumentError set_parameters!(d, dms, draw)
    end
end
