using Test
using InferCell
using JumpProcesses: SSAStepper
using StaticArrays: SVector

# Jump composition (spec §11 phase 2). Doubles in jump_test_models.jl.

@testset "Jump composition" begin
    # Task 2.1: the witness. Before phase 2 these two doubles composed and B's
    # firings landed on A's state (its affect did `integrator.u[1] += 1` against
    # the global vector) while its rate read A's parameter. Both facts were
    # asserted in the first commit of this phase; the fix turns them into a
    # named refusal at build time, so the aliasing cannot recur silently.
    @testset "2.1 a module still returning ConstantRateJumps is refused by name" begin
        a = LegacyOwner(x0 = 3, k = 0.0)
        b = LegacyWriter(b = 1.0)
        err = caught(() -> build_problem([a, b]; tspan = (0.0, 50.0)))
        @test err isa ArgumentError
        @test occursin("LegacyWriter", err.msg)
        @test occursin("ConstantRateJump", err.msg)
        @test occursin("local coordinates", err.msg)
        # An empty legacy vector has nothing to alias and is not refused.
        @test build_problem(LegacyOwner(x0 = 3)) isa JumpProcesses.JumpProblem
    end

    # Task 2.2: every read and write a module's reactions make falls inside its
    # own slice, or its declared inputs, whatever position it is composed at.
    # The parent vectors record each index the orchestrator's views touch.
    @testset "2.2 reactions address only the module's own slice and its inputs" begin
        # A three-module layout: [BirthOwner :X | PeerBirthDeath :fired | BirthDeath :X2]
        # is more than the doubles allow (two own :X), so lay the slots out by
        # hand: the owner at global 1, the peer's counter at 2, and a third
        # module's two states at 3:4 so neither tested module sits first.
        owner = BirthOwner(k = 1.0)
        peer = PeerBirthDeath(gamma = 0.5, b = 0.2)
        slot_owner = InferCell._JumpSlot(SVector(4), SVector(3), SVector{0, Int}())
        slot_peer = InferCell._JumpSlot(SVector(2), SVector(1, 2), SVector(4))
        u = TrackedVector([10, 0, 7, 5])
        p = TrackedVector([0.5, 0.2, 1.0])
        for (m, sl, allowed_u, allowed_p) in ((owner, slot_owner, Set([4]), Set([3])),
                                              (peer, slot_peer, Set([2, 4]), Set([1, 2])))
            for r in reactions(m)
                j = InferCell._global_jump(r, sl)
                empty!(u.reads); empty!(u.writes); empty!(p.reads)
                j.rate(u, p, 0.0)
                j.affect!(FakeIntegrator(u))
                @test !isempty(u.writes)
                @test issubset(u.reads, allowed_u)
                @test issubset(u.writes, allowed_u)
                @test issubset(p.reads, allowed_p)
            end
        end
        # And the values seen are the slot's, not position one's.
        r_birth = reactions(owner)[1]
        @test InferCell._global_jump(r_birth, slot_owner).rate([10, 0, 7, 5], [0.5, 0.2, 1.0], 0.0) == 1.0
        r_death = reactions(peer)[1]
        @test InferCell._global_jump(r_death, slot_peer).rate([10, 0, 7, 5], [0.5, 0.2, 1.0], 0.0) == 0.5 * 5
    end

end
