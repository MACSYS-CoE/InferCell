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

end
