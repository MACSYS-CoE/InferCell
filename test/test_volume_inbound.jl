using Test
using InferCell
using StaticArrays
using Random

# Spec §11 phase 5b — the inbound half of the volume channel: the cell's
# geometry re-entering an ODE rate law as a parameter, which task 7.2 needs for
# the lactate exporter's 3P/r.
#
# Doubles are in `hybrid_test_models.jl` and `corea_test_models.jl`, both
# included first by `runtests.jl`. `caught` comes from the latter.

@testset "Phase 5b — the inbound volume channel" begin

    # The growth partner, with no deferred counter, so nothing writes a pool
    # back through the rounding policy: the phase-5 idiom.
    _grower() = ToyGrowingExpression(edges = CouplingEdge[
        VolumeEdge(species = :M_ptsg_c, direction = :out)])

    # A reader that does *not* also flag, so the growth arithmetic is a function
    # of the jump block's count alone and `toy_rate_law_geometry(n, 831)`
    # applies directly.
    function _reader(; quantity = :radius_nm, param_slot = :r_cell_nm, kw...)
        return ToyExportingPool(; membrane = Symbol[],
                                edges = CouplingEdge[
                                    CatalyticEdge(species = :M_ptsg_c,
                                                  direction = :in,
                                                  param_slot = :enzyme_conc),
                                    VolumeEdge(species = :M_lac__L_e,
                                               direction = :in,
                                               param_slot = param_slot,
                                               quantity = quantity)],
                                kw...)
    end

    # Where the geometry lands in the *composed* ODE parameter vector. Computed
    # over the whole composition by `_build_p0`'s own rule — first-seen name,
    # ODE modules in order — rather than looked up on one module, because a
    # module-local index is only the global one when that module is first. Every
    # assertion below indexes `d.ode.p` with this, so getting it from a single
    # module would make the multi-module testset silently assert nothing.
    function _slot(models, name = :r_cell_nm)
        seen = Symbol[]
        for m in models
            formalism(m) === :ode || continue
            for q in model_free_params(parameters(m))
                q.name in seen || push!(seen, q.name)
            end
        end
        return findfirst(==(name), seen)
    end

    @testset "5b.1 the constructor refusals, each naming what is wrong" begin
        # An inbound edge with no slot: there would be nothing to write, which
        # is the "resolves and never executes" failure phase 4 closed.
        err = caught(() -> VolumeEdge(species = :M_lac__L_e, direction = :in,
                                      quantity = :radius_nm))
        @test err isa ArgumentError
        @test occursin("param_slot", sprint(showerror, err))

        # An inbound edge with no quantity. Nothing else in the codebase records
        # a unit, so this is the only place a nanometre is told from a
        # centimetre.
        err = caught(() -> VolumeEdge(species = :M_lac__L_e, direction = :in,
                                      param_slot = :r_cell_nm))
        @test err isa ArgumentError
        @test occursin("quantity", sprint(showerror, err))

        # The unit-free names are refused as hard as a nonsense one: :radius is
        # exactly the typo the vocabulary exists to catch.
        for bad in (:radius, :volume, :area, :radius_um, :mass)
            err = caught(() -> VolumeEdge(species = :M_lac__L_e, direction = :in,
                                          param_slot = :r_cell_nm, quantity = bad))
            @test err isa ArgumentError
            @test occursin("quantity", sprint(showerror, err))
        end

        # An outbound edge carrying either field would imply a write that never
        # happens.
        err = caught(() -> VolumeEdge(species = :M_ptsg_c, direction = :out,
                                      param_slot = :r_cell_nm))
        @test err isa ArgumentError
        @test occursin("no `param_slot`", sprint(showerror, err))

        err = caught(() -> VolumeEdge(species = :M_ptsg_c, direction = :out,
                                      quantity = :radius_nm))
        @test err isa ArgumentError
        @test occursin("no `quantity`", sprint(showerror, err))

        # The outbound form phase 5 wrote still constructs, unchanged, and
        # carries neither field.
        e = VolumeEdge(species = :M_ptsg_c, direction = :out)
        @test e.param_slot === nothing
        @test e.quantity === nothing

        # An inbound edge is still a volume edge in every other respect: it
        # carries no mass, so it cannot be counted in a conservation check.
        i = VolumeEdge(species = :M_lac__L_e, direction = :in,
                       param_slot = :r_cell_nm, quantity = :radius_nm)
        @test edge_kind(i) === :volume
        @test !carries_mass(i)
        @test is_consumer(i) && !is_producer(i)
        @test caught(() -> mass_contribution(i)) !== nothing
    end

    @testset "5b.2 an inbound edge names the declarer's own state" begin
        # The geometry is a sum over every flagged state in the composition, so
        # it belongs to no module and names no molecule. The species therefore
        # names *this* module's state whose rate law reads it — anything else
        # would describe a dependence that does not exist, and would quietly
        # become a different claim once a second module flags a second membrane
        # protein, which task 7.6 plans.
        stray = ToyExportingPool(edges = CouplingEdge[
            VolumeEdge(species = :M_g6p_c, direction = :in,
                       param_slot = :r_cell_nm, quantity = :radius_nm)])
        err = caught(() -> resolve_coupling([stray]))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("does not own", msg)
        @test occursin("M_g6p_c", msg)

        # Caught on the single module, without a composition — the standalone
        # property the resolver's docstring claims.
        ok = _reader()
        @test caught(() -> resolve_coupling([ok])) === nothing
    end

    @testset "5b.3 the lowering refusals, each naming what is wrong" begin
        # A slot that is not one of the module's own free parameters. A *fixed*
        # parameter is the same failure: it never reaches the composed vector,
        # so there is no slot to write.
        err = caught(() -> build_problem(
            [_reader(param_slot = :not_a_parameter), _grower()];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("not one of its own free parameters", sprint(showerror, err))

        # Two writers, one slot: both lower, both write, the later wins, and one
        # channel is declared, lowered, reported and has no effect.
        err = caught(() -> build_problem(
            [ToyExportingPool(membrane = Symbol[], edges = CouplingEdge[
                 CatalyticEdge(species = :M_ptsg_c, direction = :in,
                               param_slot = :enzyme_conc),
                 VolumeEdge(species = :M_lac__L_e, direction = :in,
                            param_slot = :enzyme_conc, quantity = :radius_nm)]),
             _grower()];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("One writer per slot", sprint(showerror, err))

        # A jump module cannot receive the geometry: the write targets the ODE
        # parameter vector, and its parameters live in the propensity vector.
        err = caught(() -> build_problem(
            [_reader(),
             ToyGrowingExpression(edges = CouplingEdge[
                 VolumeEdge(species = :M_ptsg_c, direction = :out),
                 VolumeEdge(species = :M_ptsg_c, direction = :in,
                            param_slot = :k_tx_toy, quantity = :radius_nm)])];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        msg = sprint(showerror, err)
        @test occursin("jump module", msg)
        @test occursin("RateConstantEdge", msg)   # names the remedy it does have
    end

    @testset "5b.4 a module may both flag a protein and read the geometry" begin
        # The load-bearing case for the producer filter in `_lower_growth`.
        # Undirected, its `e.species in names` check would tell this module that
        # its lactate edge names an unflagged species — a correct declaration
        # refused with a message about the wrong problem.
        d = build_problem([ToyExportingPool(), _grower()]; tspan = (0.0, 10.0))
        @test length(d.geometry) == 1
        @test only(d.geometry).species === :M_lac__L_e
        @test only(d.geometry).quantity === :radius_nm

        # Both flagged states still reach the surface area: the inbound edge
        # satisfied no flag's obligation.
        @test Set(g.species for g in d.growth) == Set([:M_ptsi_c, :M_ptsg_c])

        # And a flag whose only edge is inbound is still refused — the inbound
        # edge must not stand in for the outbound one.
        err = caught(() -> build_problem(
            [ToyExportingPool(edges = CouplingEdge[
                 VolumeEdge(species = :M_lac__L_e, direction = :in,
                            param_slot = :r_cell_nm, quantity = :radius_nm)]),
             _grower()];
            tspan = (0.0, 10.0)))
        @test err isa ArgumentError
        @test occursin("no outbound VolumeEdge", sprint(showerror, err))
    end

    @testset "5b.5 the slot holds this hook's geometry, at every handshake" begin
        Random.seed!(5101)
        m = _reader(perm = 0.0)
        models = [m, _grower()]
        d = build_problem(models; tspan = (0.0, 20.0))
        slot = _slot(models)

        # Before any handshake the slot is the module's declared value; after
        # one it is the driver's. Those differ in the fourth significant figure
        # — 200.0 against radius_from_area_nm(502831) = 200.03505 — so this
        # assertion cannot pass by the declared value happening to be right.
        @test d.ode.p[slot] == COREA_INITIAL_RADIUS_NM
        handshake_step!(d)
        @test d.ode.p[slot] != COREA_INITIAL_RADIUS_NM
        # `≈` rather than `==`: the closed form goes through the volume and back,
        # which the driver only does above the cap, and the round trip costs an
        # ulp. Asserting to 1e-12 says the delivered value is the sphere's
        # radius however it was arrived at.
        @test d.ode.p[slot] ≈ toy_rate_law_geometry(831, 831).radius_nm rtol = 1e-12

        # A count that moves moves the slot, and by this hook's value rather
        # than the last one's. A one-handshake lag at published parameters is a
        # ~1e-5 relative error that no tolerance test would notice, so the count
        # is driven hard enough that a lag is a whole percent, and the
        # comparison against the record is exact.
        d.jump.u[2] = 20_000
        rec = run_handshake!(d, 1)
        expected = toy_rate_law_geometry(20_000, 831).radius_nm
        @test rec.ode_p[end][slot] ≈ expected rtol = 1e-12
        @test abs(rec.ode_p[end][slot] - toy_rate_law_geometry(831, 831).radius_nm) /
              expected > 0.01
        @test rec.ode_p[end][slot] === d.ode.p[slot]

        # Read off the recorded trajectory rather than the schedule: every
        # handshake's slot is that handshake's geometry.
        m2 = _reader(perm = 0.0)
        d2 = build_problem([m2, _grower()]; tspan = (0.0, 20.0))
        slot2 = _slot([m2, _grower()])
        Random.seed!(5102)
        counts = [900, 1200, 1200, 5000, 831]
        vals = Float64[]
        for n in counts
            d2.jump.u[2] = n
            handshake_step!(d2)
            push!(vals, d2.ode.p[slot2])
        end
        @test all(i -> vals[i] ≈ toy_rate_law_geometry(counts[i], 831).radius_nm,
                  eachindex(counts))
    end

    @testset "5b.6 all three quantities, each in its own unit" begin
        Random.seed!(5103)
        for q in VOLUME_QUANTITIES
            m = _reader(perm = 0.0, quantity = q)
            models = [m, _grower()]
            d = build_problem(models; tspan = (0.0, 10.0))
            slot = _slot(models)
            d.jump.u[2] = 1662
            handshake_step!(d)
            expected = getproperty(toy_rate_law_geometry(1662, 831), q)
            @test d.ode.p[slot] ≈ expected rtol = 1e-12
        end

        # The three differ by orders of magnitude — ~2e2 nm, ~3.6e-17 L,
        # ~5e5 nm² — so a mix-up between them cannot survive any tolerance, and
        # asserting that is what makes the unit-bearing vocabulary load-bearing
        # rather than a naming style.
        g = toy_rate_law_geometry(1662, 831)
        @test g.radius_nm > 1e2 && g.radius_nm < 1e3
        @test g.volume_litres < 1e-15
        @test g.area_nm2 > 1e5

        # Dimensional consistency, written longhand and independently of the
        # implementation: 3/r in cm must equal A/V in cm. This is the assertion
        # that fails on any unit slip, in the vocabulary or in a rate law.
        r_cm = g.radius_nm * 1e-7
        a_cm2 = g.area_nm2 * 1e-14
        v_cm3 = g.volume_litres * 1e3
        @test 3 / r_cm ≈ a_cm2 / v_cm3 rtol = 1e-12
    end

    @testset "5b.7 a fixed cell is written, behind a second ODE module" begin
        # Two things at once, because they share a composition.
        #
        # A fixed cell: nothing flags a membrane protein, so `_update_volume!`
        # returns early. `ToyFixedExporter` declares `r_cell = 1.0`, so a write
        # skipped by living inside that early return is a 200x error, not a
        # subtle one.
        #
        # And the geometry slot sits on the *second* ODE module, behind
        # `ToyPool`'s three free parameters. That is the only place
        # `ode_contexts[i].param_idxs[j]` is exercised at `i > 1`: a mapping that
        # confused the module-local index for the global one would write into
        # `ToyPool`'s rate law instead — in range, type-correct, and silent.
        Random.seed!(5104)
        m = ToyFixedExporter(perm = 0.0)
        models = [ToyPool(kcat = 0.0), m, ToyExpression()]
        d = build_problem(models; tspan = (0.0, 10.0))
        slot = _slot(models)

        @test isempty(d.growth)          # the cell does not grow
        @test length(d.geometry) == 1    # and still reads the geometry
        # The offset is real: local index 2 on its own module, global 5 here.
        @test slot == 5
        @test _slot([m]) == 2

        @test d.ode.p[slot] == 1.0
        handshake_step!(d)
        @test d.ode.p[slot] == COREA_INITIAL_RADIUS_NM
        # ToyPool's own slots are untouched by the geometry write.
        @test d.ode.p[1] == 0.0          # kcat_toy, as declared
        @test d.ode.p[2] == 1.0          # km_toy, as declared

        # It stays there, rather than being written once and drifting.
        rec = run_handshake!(d, 5)
        @test all(p -> p[slot] == COREA_INITIAL_RADIUS_NM, rec.ode_p)

        # A fixed cell given an explicit radius reads that one.
        m2 = ToyFixedExporter(perm = 0.0)
        models2 = [ToyPool(kcat = 0.0), m2, ToyExpression()]
        d2 = build_problem(models2; tspan = (0.0, 10.0), radius_nm = 250.0)
        handshake_step!(d2)
        @test d2.ode.p[_slot(models2)] == 250.0
    end

    @testset "5b.8 above the cap the rate law still sees one sphere" begin
        # `_update_volume!` caps the volume and leaves the radius uncapped, as
        # in_out.py:107 does. That is harmless while the radius is only
        # reported and wrong once it drives 3P/r, whose whole content is that
        # 3/r is a sphere's surface-to-volume ratio. So above the cap the radius
        # and area a rate law receives are reconstructed from the capped volume
        # and the three quantities again describe one sphere; below it nothing
        # is reconstructed.
        Random.seed!(5105)
        m = _reader(perm = 0.0)
        models = [m, _grower()]
        d = build_problem(models; tspan = (0.0, 20.0))
        slot = _slot(models)

        # Below the cap, the rate law's radius and the reported one agree.
        handshake_step!(d)
        @test !growth_census(d).capped
        # Exactly, not approximately: below the cap the reported geometry is
        # already one sphere, so it is passed through rather than reconstructed.
        @test d.ode.p[slot] == d.radius_nm

        # Above it they must not, and the rate law's is the smaller.
        d.jump.u[2] = 50_000
        handshake_step!(d)
        @test growth_census(d).capped
        @test d.ode.p[slot] < d.radius_nm
        @test cell_volume_litres(d.ode.p[slot]) ≈ d.volume_litres rtol = 1e-12

        # The reported radius is still the published uncapped one, so nothing
        # phase 5 asserted has changed.
        @test d.radius_nm ≈ toy_growth(50_000, 831).radius_nm rtol = 1e-12

        # And the departure is declared rather than left in a comment.
        labels = driver_declarations(d)
        @test any(l -> l.category === :capped_rate_law_geometry, labels)
    end

    @testset "5b.9 the channel is visible from outside" begin
        Random.seed!(5106)
        m = _reader(perm = 0.0)
        models = [m, _grower()]
        d = build_problem(models; tspan = (0.0, 10.0))
        handshake_step!(d)

        c = only(growth_census(d).consumers)
        @test c.species === :M_lac__L_e
        @test c.quantity === :radius_nm
        @test c.param_slot === :r_cell_nm
        @test c.declared_by === :ToyExportingPool
        @test c.value == d.ode.p[_slot(models)]

        # The row is read back out of the parameter vector, not recomputed: a
        # census that recomputed would report the right geometry even if the
        # write never landed, which is the one thing it exists to witness.
        d.ode.p[_slot(models)] = -1.0
        @test only(growth_census(d).consumers).value == -1.0

        # A composition with no consumer says so, rather than omitting the row.
        plain = build_problem([ToyPool(kcat = 0.0), _grower()]; tspan = (0.0, 10.0))
        @test isempty(growth_census(plain).consumers)
        @test !isempty(d.geometry)

        # Every driver-written slot is enumerable, across all three channels —
        # the first step to excluding them from the sampled set, which is task
        # 13.3's. They are free parameters, so inference samples them today and
        # the driver overwrites the draw; a posterior for one is its prior.
        written = driver_written_params(d)
        @test (param_slot = :r_cell_nm, channel = :volume,
               declared_by = :ToyExportingPool) in written
        @test any(w -> w.channel === :catalytic, written)

        # The third channel, on phase 4's own composition. Task 13.3 reads this
        # to decide what not to sample, so a rebuilt slot silently missing from
        # it is a posterior that gets read as an identifiability result.
        reb = build_problem([toy_slow_pool(), ToyRebuiltExpression()];
                            tspan = (0.0, 120.0))
        @test (param_slot = :k_tx_rb, channel = :rate_constant,
               declared_by = :ToyRebuiltExpression) in driver_written_params(reb)
    end

    @testset "Done when: an inbound edge resolves, executes, and drives a flux" begin
        # The whole channel, on a rate law shaped like task 7.2's: the geometry
        # reaches 3P/r, and the flux it produces tracks the closed form as the
        # cell grows. Without the channel the flux would be constant.
        Random.seed!(5107)
        m = _reader(perm = 1e-4)
        models = [m, _grower()]
        d = build_problem(models; tspan = (0.0, 60.0))
        slot = _slot(models)
        lac_c, lac_e = 4, 5

        d.jump.u[2] = 831
        handshake_step!(d)
        r_small, cyto = d.ode.p[slot], d.ode.u[lac_c]

        d.jump.u[2] = 30_000
        handshake_step!(d)
        r_big = d.ode.p[slot]

        # The cell grew, so 3P/r fell, so the export slowed — compared at one
        # common cytosolic concentration, so it is the geometry and not the pool.
        @test r_big > r_small
        @test toy_export(r_big, 1e-4, cyto) < toy_export(r_small, 1e-4, cyto)

        # Lactate really moved, and in the export direction: out of the cytosol,
        # into the medium. The medium pool is exempt from dilution, so its rise
        # is the export and nothing else.
        @test d.ode.u[lac_c] < 2.0
        @test d.ode.u[lac_e] > 0.0

        @info "phase 5b done-when" handshakes = d.n_handshakes r_small = r_small r_big = r_big cytosolic = d.ode.u[lac_c] exported = d.ode.u[lac_e] consumers = length(d.geometry)
    end
end
