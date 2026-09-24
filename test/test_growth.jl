using Test
using InferCell
using StaticArrays
using Random

# Spec §11 phase 5 — growth and volume.
#
# Doubles are in `hybrid_test_models.jl` and `corea_test_models.jl`, both
# included first by `runtests.jl`. `caught` comes from the latter.

@testset "Phase 5 — growth and volume" begin

    @testset "5.1 the membrane flag is data on the module, found by sweeping" begin
        # The flag is declared by whichever module owns the state, and the
        # composition is *swept* for it. Nothing here names a module type, and
        # the sweep is run at several positions so the answer cannot depend on
        # where the declaring module sits.
        pts = CoreAStub(:PtsLike; st = [:M_ptsg_c, :M_ptsg_P_c],
                        membrane = [:M_ptsg_c, :M_ptsg_P_c])
        glyc = CoreAStub(:GlycolysisLike; st = [:M_g6p_c, :M_f6p_c])
        nucl = CoreAStub(:RecyclingLike; st = [:M_atp_c, :M_amp_c])

        for composition in ([pts, glyc, nucl], [glyc, pts, nucl], [glyc, nucl, pts])
            @test membrane_protein_states(composition) == [:M_ptsg_c, :M_ptsg_P_c]
        end

        # The default really is empty, so a module that flags nothing needs no
        # method — and a composition of such modules reports none.
        @test isempty(membrane_protein_states(glyc))
        @test isempty(membrane_protein_states(ToyPool()))
        @test isempty(membrane_protein_states([glyc, nucl]))
    end

    @testset "5.1 the refusals, each naming what is wrong" begin
        # A flag on a state the declaring module does not own. The whole point
        # of the sweep is that the growth code need not know where to look, so
        # a flag that is not on the owner would make the answer depend on which
        # module happened to be asked.
        stray = CoreAStub(:Stray; st = [:M_g6p_c], membrane = [:M_ptsg_c])
        err = caught(() -> membrane_protein_states([stray]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("Stray", msg)
        @test occursin("M_ptsg_c", msg)
        @test occursin("does not own", msg)

        # The same species flagged twice: its count would enter the area twice.
        a = CoreAStub(:First; st = [:M_ptsg_c], membrane = [:M_ptsg_c])
        b = CoreAStub(:Second; st = [:M_ptsg_c], membrane = [:M_ptsg_c])
        err = caught(() -> membrane_protein_states([a, b]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("First", msg) && occursin("Second", msg)
        @test occursin("twice", msg)
    end

    @testset "5.2 the surface-area, radius and volume chain" begin
        # The published initial area is the primitive and the radius is derived
        # from it, which is the only way the growth arithmetic of task 5.5
        # closes. It is *not* an exact 200 nm sphere: 4pi(200 nm)^2 is
        # 502,654.8 nm^2, so the published area is 176 nm^2 larger and the
        # radius comes out at 200.035 nm. The published 200 nm is reproduced
        # to the four significant figures the source states it in, and the
        # full-precision value is asserted alongside so the gap is on the
        # record rather than hidden inside a tolerance.
        r0 = radius_from_area_nm(COREA_INITIAL_SURFACE_AREA_NM2)
        @test COREA_INITIAL_SURFACE_AREA_NM2 == 502831.0
        @test round(r0; digits = 1) == 200.0
        @test r0 ≈ 200.035046 rtol = 1e-9
        @test 4 * π * 200.0^2 ≈ 502654.8 atol = 0.1
        @test COREA_INITIAL_SURFACE_AREA_NM2 - 4 * π * 200.0^2 ≈ 176.2 atol = 0.1

        # The chain is derived, not transcribed: a sphere of the returned
        # radius has the area it was built from.
        @test 4 * π * r0^2 ≈ COREA_INITIAL_SURFACE_AREA_NM2 rtol = 1e-12

        # The frozen baseline. Deriving it is what keeps 831 copies at 28 nm^2
        # a *part* of the published area rather than all of it; getting this
        # wrong is what would make doubling ptsG double the whole cell.
        base = membrane_area_baseline_nm2(COREA_INITIAL_SURFACE_AREA_NM2, 831)
        @test base == 479563.0
        @test surface_area_nm2(base, 831) == COREA_INITIAL_SURFACE_AREA_NM2
        @test MEMBRANE_PROTEIN_FOOTPRINT_NM2 == 28.0
        @test 831 * MEMBRANE_PROTEIN_FOOTPRINT_NM2 == 23268.0

        # The cap is twice the initial volume, applied to the volume as the
        # published code does. Its hard-coded 6.70e-17 is twice its own rounded
        # 3.35e-17; ours is twice the computed initial volume.
        v0 = cell_volume_litres(r0)
        @test v0 ≈ 3.3528e-17 rtol = 1e-4
        @test corea_volume_cap_litres(v0) == 2 * v0
        @test corea_volume_cap_litres(v0) ≈ 6.70e-17 rtol = 1e-2
    end

    # The conversion factor a growing composition starts at. Derived from the
    # published *area*, so 20,191 particles per mM rather than the 20,180 an
    # exactly-200 nm cell gives — a 0.05% difference, and the reason the driver
    # refuses to be handed both a radius and an area.
    F0 = corea_particles_per_mM(radius_from_area_nm(COREA_INITIAL_SURFACE_AREA_NM2))

    # No deferred counter, so nothing writes a pool back through the rounding
    # policy and the only thing that can move an ODE state is the dilution.
    _no_counter() = ToyGrowingExpression(edges = CouplingEdge[
        VolumeEdge(species = :M_ptsg_c, direction = :out)])

    @testset "5.1 a composition flagging none is untouched" begin
        d = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 200.0))
        @test isempty(d.growth)
        @test isempty(d.dilute_idxs)
        @test d.factor == corea_particles_per_mM()          # the pre-phase-5 constant
        @test d.radius_nm == COREA_INITIAL_RADIUS_NM
        @test all(l -> l.category === :driver_policy, driver_declarations(d))

        # In trajectory form: with nothing flagged the cell never moves, so the
        # factor and the geometry are bitwise constant for the whole run.
        Random.seed!(5001)
        rec = run_handshake!(d, 100)
        @test all(g -> g == rec.growth[1], rec.growth)
        @test rec.growth[1].factor == corea_particles_per_mM()
        @test rec.growth[1].fractional == 1.0
        @test growth_census(d).n_membrane == 0.0

        # And the protocol defaults really are empty, so a module that flags
        # nothing needs no method.
        @test isempty(membrane_protein_states(ToyExpression()))
        @test isempty(extracellular_states(ToyPool()))
    end

    @testset "5.1 the driver-side refusals, each naming what is wrong" begin
        # A flag with no edge: the channel would execute undeclared.
        err = caught(() -> build_problem(
            [ToyPool(), ToyGrowingExpression(edges = CouplingEdge[
                DeferredCounterEdge(species = :M_atp_c, direction = :in,
                                    counter = :atp_cost)])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        # "outbound" since phase 5b: an inbound edge does not discharge a flag's
        # obligation, so the message has to say which direction is missing.
        @test occursin("no outbound VolumeEdge", sprint(showerror, err))

        # An edge with nothing flagged: declared and never executed, the trap
        # phase 4 spent a task closing on the rate-constant channel.
        err = caught(() -> build_problem(
            [ToyPool(), ToyGrowingExpression(membrane = Symbol[])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("never", msg) && occursin("VolumeEdge", msg)

        # The direction convention read backwards. Phase 5 refused an inbound
        # edge here, at build time, because there was no half of the channel to
        # execute it. Phase 5b built that half, so a bare inbound edge is now
        # refused *earlier* — by the constructor, for carrying no `param_slot`
        # and no `quantity` — and `test_volume_inbound.jl` asserts it there.
        # What this file still owns is that a bare inbound edge never reaches
        # the flag pairing at all.
        err = caught(() -> VolumeEdge(species = :M_ptsg_c, direction = :in))
        @test err isa ArgumentError
        @test occursin("param_slot", sprint(showerror, err))

        # An edge on a species this module does not flag: its count would reach
        # no surface area.
        err = caught(() -> build_problem(
            [ToyPool(), ToyGrowingExpression(edges = CouplingEdge[
                VolumeEdge(species = :M_atp_c, direction = :out)])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("does not flag", sprint(showerror, err))

        # A flag on a state this module does not own, caught by the sweep.
        err = caught(() -> build_problem(
            [ToyPool(), ToyGrowingExpression(membrane = [:M_ptsi_c])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("does not own", sprint(showerror, err))

        # Stating the geometry twice. The published area is not a 200 nm sphere,
        # so a radius alongside it is not a harmless restatement.
        err = caught(() -> build_problem([ToyPool(), _no_counter()];
                                         tspan = (0.0, 10.0), radius_nm = 400.0))
        @test err isa ArgumentError
        @test occursin("initial_surface_area_nm2", sprint(showerror, err))

        # A jump module has counts, not concentrations, so it has nothing to
        # exempt from dilution.
        err = caught(() -> build_problem(
            [ToyPool(), ToyGrowingExpression(extracellular = [:M_ptsg_c])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("nothing here to exempt", sprint(showerror, err))

        # An exemption on a state the declaring module does not own.
        err = caught(() -> build_problem(
            [ToyMembranePool(extracellular = [:M_g6p_c]), _no_counter()];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("does not own", sprint(showerror, err))

        # Exempt *and* flagged: the area would depend on the volume it sets.
        err = caught(() -> build_problem(
            [ToyMembranePool(extracellular = [:M_ptsi_c]), _no_counter()];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("both extracellular and a membrane", sprint(showerror, err))
    end

    @testset "5.2 the chain is live on the driver, and labelled" begin
        d = build_problem([ToyPool(kcat = 0.0), _no_counter()]; tspan = (0.0, 10.0))
        c = growth_census(d)

        # The composition starts at exactly the published area, whatever count
        # it happens to hold, because the baseline is derived from that count.
        @test c.n_membrane == 831.0
        @test c.area_nm2 == COREA_INITIAL_SURFACE_AREA_NM2
        @test c.baseline_nm2 == 479563.0
        @test c.footprint_nm2 == 28.0
        @test c.radius_nm ≈ 200.035046 rtol = 1e-9
        @test round(c.radius_nm; digits = 1) == 200.0
        @test d.factor ≈ 20191.0 atol = 0.5
        @test !c.capped

        # The two labels the growth chain carries, both retrievable rather than
        # living in a docstring.
        labels = driver_declarations(d)
        cats = [l.category for l in labels]
        @test :calibrated_constant in cats
        @test :exogenous_growth in cats
        report = reduction_report([ToyPool(kcat = 0.0), _no_counter()], d)
        @test occursin("35 nm²", report)
        @test occursin("calibrated", report)
        @test occursin("479563.0 nm²", report)
    end

    @testset "5.3 growth dilutes every ODE concentration, at constant count" begin
        # One membrane protein in each block: the jump side's count is whole
        # particles, the ODE side's is a concentration read at the live factor.
        # M_lac__L_e is declared extracellular, so it must not dilute.
        models = [ToyMembranePool(), _no_counter()]
        d = build_problem(models; tspan = (0.0, 20.0))
        @test [g.block for g in d.growth] == [:ode, :jump]
        @test [g.species for g in d.growth] == [:M_ptsi_c, :M_ptsg_c]
        @test d.dilute_idxs == [1, 2, 3]                 # every ODE state but lac_e

        # Nothing has happened yet, so the first handshake must not move the
        # cell: the baseline was derived from exactly these counts.
        Random.seed!(5002)
        handshake_step!(d)
        @test d.area_nm2 ≈ COREA_INITIAL_SURFACE_AREA_NM2 rtol = 1e-12
        @test d.factor ≈ F0 rtol = 1e-12

        counts_before = [d.ode.u[i] * d.factor for i in 1:4]
        conc_before = collect(Float64, d.ode.u)
        rem_before = copy(d.rounding.remainders)
        f_before = d.factor

        # Double the membrane protein the jump block owns.
        d.jump.u[2] = 1662
        handshake_step!(d)

        # A millimolar of a bigger cell is more particles, so the factor rises
        # and every concentration at fixed count falls by the volume ratio.
        ratio = f_before / d.factor
        @test d.factor > f_before
        @test ratio ≈ (200.035046 / 204.610919)^3 rtol = 1e-6
        @test ratio < 1

        for i in d.dilute_idxs
            # Diluted by the volume ratio, and the count is what survives. It is
            # a one-ulp statement rather than a bitwise one: (u*f_old)/f_new
            # multiplied back by f_new need not return u*f_old in binary
            # floating point.
            @test d.ode.u[i] ≈ conc_before[i] * ratio rtol = 1e-12
            @test d.ode.u[i] * d.factor ≈ counts_before[i] rtol = 1e-14
        end
        # The exempt state is referred to the medium, so it is untouched — in
        # concentration, which means its count *does* move, and that is the
        # point: it is not a count in this cell.
        @test d.ode.u[4] == conc_before[4]

        # The carried remainders are in particles and the rescale preserves
        # counts, so they are untouched by it.
        @test d.rounding.remainders == rem_before

        # Both conversion directions read the live volume: the catalytic slot is
        # the protein count at the *new* factor, not the old one.
        @test d.ode.p[3] == counts_to_mM(1662, d.factor)
        @test d.ode.p[3] != counts_to_mM(1662, f_before)
    end

    @testset "5.4 growth stops at exactly twice the initial volume" begin
        d = build_problem([ToyMembranePool(), _no_counter()]; tspan = (0.0, 20.0))
        v0 = d.initial_volume_litres
        cap = d.volume_cap_litres
        @test cap == 2 * v0

        # Twice the volume is 2^(2/3) times the area, so the cap needs about
        # 11,400 membrane proteins — three orders above the ~1.07x Core A' can
        # actually reach, which is why it has to be driven there.
        Random.seed!(5003)
        d.jump.u[2] = 50_000
        handshake_step!(d)
        @test d.volume_litres == cap
        @test growth_census(d).capped
        @test growth_census(d).fractional == 2.0

        # Still capped when the count keeps rising, and the *radius* is left
        # uncapped, exactly as in_out.py:107 leaves it.
        r_capped = d.radius_nm
        d.jump.u[2] = 500_000
        handshake_step!(d)
        @test d.volume_litres == cap
        @test d.radius_nm > r_capped
        @test d.factor == particles_per_mM(cap)
    end

    @testset "5.5 the scoping note's ptsG arithmetic" begin
        d = build_problem([ToyPool(kcat = 0.0), _no_counter()]; tspan = (0.0, 20.0))
        Random.seed!(5004)
        handshake_step!(d)

        at831 = growth_census(d)
        expected831 = toy_growth(831, 831)
        @test at831.area_nm2 ≈ expected831.area_nm2 rtol = 1e-12
        @test at831.radius_nm ≈ expected831.radius_nm rtol = 1e-12
        @test at831.volume_litres ≈ expected831.volume_litres rtol = 1e-12

        d.jump.u[2] = 1662
        rec = run_handshake!(d, 1)
        at1662 = rec.growth[end]
        expected1662 = toy_growth(1662, 831)

        # Against the closed form, then against the note's own three numbers.
        @test at1662.area_nm2 ≈ expected1662.area_nm2 rtol = 1e-12
        @test at1662.radius_nm ≈ expected1662.radius_nm rtol = 1e-12
        @test at1662.fractional ≈ expected1662.fractional rtol = 1e-12

        @test at831.area_nm2 == 502831.0
        @test at1662.area_nm2 == 526099.0
        @test round(at831.radius_nm; digits = 1) == 200.0
        @test round(at1662.radius_nm; digits = 1) == 204.6
        @test round(at1662.fractional; digits = 2) == 1.07

        @info "task 5.5 the ptsG arithmetic" area_831 = at831.area_nm2 area_1662 = at1662.area_nm2 radius_831 = at831.radius_nm radius_1662 = at1662.radius_nm fractional = at1662.fractional
    end

    @testset "5.6 the doubling-time comparison is refused in code" begin
        d = build_problem([ToyPool(kcat = 0.0), _no_counter()]; tspan = (0.0, 20.0))

        # Retrievable from the composed model, not only from prose.
        cs = reporting_constraints(d)
        c = only(filter(x -> x.quantity === :doubling_time, cs))
        @test c.verdict === :refused
        @test occursin("92%", c.reason)
        @test c.instead == [:fractional_growth, :time_to_threshold]

        # And refused when asked for.
        err = caught(() -> doubling_time(d))
        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("refused", msg)
        @test occursin("fractional_growth", msg) && occursin("time_to_threshold", msg)

        # The two admissible reportings.
        Random.seed!(5005)
        d.jump.u[2] = 1662
        rec = run_handshake!(d, 3)
        g = growth_report(d)
        @test g.growing
        @test g.fractional ≈ 1.0702 atol = 1e-3
        @test !haskey(pairs(g), :doubling_time)
        @test time_to_threshold(rec, 1.05) == rec.t[1]
        @test time_to_threshold(rec, 1.5) === nothing

        # 2.0 is the cap ratio, so `fractional` saturates there: the query would
        # return the handshake growth *stopped* at, which is a doubling time
        # wearing another name. Refused for the same reason.
        err = caught(() -> time_to_threshold(rec, 2.0))
        @test err isa ErrorException
        @test occursin("doubling time wearing another name", sprint(showerror, err))
        @test all(g -> !g.capped, rec.growth)

        # A refusal that does not depend on this composition happening to grow.
        fixed = build_problem([ToyPool(), ToyExpression()]; tspan = (0.0, 10.0))
        @test !growth_report(fixed).growing
        @test any(x -> x.quantity === :doubling_time, reporting_constraints(fixed))
    end

    @testset "5.1 a flagged state with no edge of its own is refused" begin
        # Distinct from "an edge on a species the module does not flag": here
        # the module flags two states and edges only one, so the second state's
        # count would reach the surface area with nothing declaring that it
        # does. The trap phase 7 walks into if it flags both ptsG phospho-forms
        # and declares one edge.
        err = caught(() -> build_problem(
            [ToyPool(), ToyGrowingExpression(membrane = [:M_ptsg_c, :atp_cost])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("atp_cost", msg)
        @test occursin("One outbound edge per flagged", msg)
    end

    @testset "5.3 the ODE-side count is particles, asserted independently" begin
        # The build-time baseline and every hook read counts through the same
        # function, so an error in the ODE branch of that conversion cancels
        # against the baseline and is invisible in the area. Assert the count
        # itself: 353 ptsI (ODE, held as a concentration) plus 831 ptsG (jump,
        # already particles).
        d = build_problem([ToyMembranePool(), _no_counter()]; tspan = (0.0, 20.0))
        c = growth_census(d)
        @test c.n_membrane ≈ 1184.0 rtol = 1e-12
        @test c.area_nm2 ≈ COREA_INITIAL_SURFACE_AREA_NM2 rtol = 1e-12
        @test c.baseline_nm2 ≈ 502831.0 - 1184.0 * 28.0 rtol = 1e-12

        # Were the ODE branch to return mM rather than particles, the count
        # would be 831 + 0.0175 and the baseline would absorb it silently — so
        # the count is the assertion that catches it, not the area.
        @test c.n_membrane - 831.0 ≈ 353.0 rtol = 1e-12
    end

    @testset "5.2 the frozen fraction is reported against the *initial* area" begin
        # The label is generated beside a result, which by then is many
        # handshakes old. Reporting the live area would state a frozen fraction
        # that falls as the cell grows — 89.5% instead of 95.4% on the done-when
        # composition — and it is the label that travels into §6 T2.
        models = [ToyPool(kcat = 0.0), _no_counter()]
        d = build_problem(models; tspan = (0.0, 20.0))
        Random.seed!(5007)
        d.jump.u[2] = 5000
        handshake_step!(d)
        @test d.area_nm2 > COREA_INITIAL_SURFACE_AREA_NM2      # the cell grew
        report = reduction_report(models, d)
        @test occursin("479563.0 nm² of the initial 502831.0 nm²", report)
        @test !occursin("of the initial $(round(d.area_nm2; digits = 1))", report)
    end

    @testset "5.2 the geometry is stated once, and has to be sane" begin
        # Non-positive area, and a footprint that would freeze the chain while
        # still labelling it.
        err = caught(() -> build_problem([ToyPool(kcat = 0.0), _no_counter()];
                                         tspan = (0.0, 10.0),
                                         initial_surface_area_nm2 = 0.0))
        @test err isa ArgumentError
        @test occursin("must be positive", sprint(showerror, err))

        err = caught(() -> build_problem([ToyPool(kcat = 0.0), _no_counter()];
                                         tspan = (0.0, 10.0), footprint_nm2 = 0.0))
        @test err isa ArgumentError
        @test occursin("declared, labelled and frozen", sprint(showerror, err))

        # A baseline that is not a baseline: 831 proteins at 28 nm² exceed the
        # area they are supposed to be part of, so the law stops being affine.
        err = caught(() -> build_problem([ToyPool(kcat = 0.0), _no_counter()];
                                         tspan = (0.0, 10.0),
                                         initial_surface_area_nm2 = 1000.0))
        @test err isa ArgumentError
        @test occursin("super-linear", sprint(showerror, err))

        # And the mirror of the radius refusal: a growth keyword on a cell that
        # does not grow would be accepted and never read.
        for kw in (:initial_surface_area_nm2, :footprint_nm2)
            err = caught(() -> build_problem([ToyPool(), ToyExpression()];
                                             tspan = (0.0, 10.0),
                                             (kw => 600000.0,)...))
            @test err isa ArgumentError
            @test occursin("never be read", sprint(showerror, err)) ||
                  occursin("never read", sprint(showerror, err))
        end
    end

    @testset "Done when: the volume channel runs beside the other three" begin
        # The catalytic channel, the deferred debit and growth, all live, on the
        # composition phase 3 built. Only the 60 s rebuild is absent, and phase
        # 4's own done-when covers that.
        Random.seed!(5006)
        d = build_problem([ToyPool(kcat = 0.0), ToyGrowingExpression(k_tx = 2.0)];
                          tspan = (0.0, 300.0))
        rec = run_handshake!(d, 300)

        @test d.n_handshakes == 300
        @test d.jump.u[2] > 831                       # translation added membrane protein
        @test rec.growth[end].area_nm2 > rec.growth[1].area_nm2
        @test rec.growth[end].fractional > 1.0
        @test rec.growth[end].factor > rec.growth[1].factor   # a bigger cell, so more particles per mM
        @test rec.census.clipped == 0

        # The area at each handshake is the published law at the count the jump
        # block left at the previous one — the same one-handshake lag the
        # catalytic channel carries, asserted the same way. Radius, volume and
        # factor follow from it, so the whole chain is checked and not only its
        # first link.
        @test all(i -> isapprox(rec.growth[i].area_nm2,
                                toy_growth(rec.jump[i - 1][2], 831).area_nm2;
                                rtol = 1e-12), 2:300)
        last = toy_growth(rec.jump[299][2], 831)
        @test rec.growth[end].radius_nm ≈ last.radius_nm rtol = 1e-12
        @test rec.growth[end].volume_litres ≈ last.volume_litres rtol = 1e-12
        @test rec.growth[end].factor ≈ particles_per_mM(last.volume_litres) rtol = 1e-12
        @test rec.growth[end].fractional ≈ last.fractional rtol = 1e-12

        @info "phase 5 done-when" handshakes = d.n_handshakes ptsg = d.jump.u[2] area_nm2 = d.area_nm2 radius_nm = d.radius_nm fractional = growth_report(d).fractional
    end

end
