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
        # radius comes out at 200.0349 nm. The published 200 nm is reproduced
        # to the four significant figures the source states it in, and the
        # full-precision value is asserted alongside so the gap is on the
        # record rather than hidden inside a tolerance.
        r0 = radius_from_area_nm(COREA_INITIAL_SURFACE_AREA_NM2)
        @test COREA_INITIAL_SURFACE_AREA_NM2 == 502831.0
        @test round(r0; digits = 1) == 200.0
        @test r0 ≈ 200.0349 atol = 1e-4
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

end
