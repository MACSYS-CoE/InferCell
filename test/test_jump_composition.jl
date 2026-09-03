using Test
using InferCell
using JumpProcesses: SSAStepper
using StaticArrays: SVector
using Random
using Statistics: mean, var

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
        slot_owner = InferCell._JumpSlot(SVector(4), SVector(3), SVector{0, Int}(),
                                         SVector{0, Bool}(), SVector{0, Symbol}(), :BirthOwner)
        slot_peer = InferCell._JumpSlot(SVector(2), SVector(1, 2), SVector(4),
                                        SVector(true), SVector(:X), :PeerBirthDeath)
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


    # Task 2.3: the composed layout — state counts concatenate in module order,
    # each module's initial counts sit at its own offsets, and a parameter name
    # shared by two modules occupies one global slot that both slices read.
    @testset "2.3 the composed layout has one slice per module" begin
        owner = BirthOwner(k = 0.7, x0 = 12)
        peer = PeerBirthDeath(gamma = 0.3, b = 0.7, bname = :k_birth)   # shares :k_birth
        models = [owner, peer]
        prob = build_problem(models; tspan = (0.0, 1.0))
        @test length(prob.prob.u0) == sum(length(states(m)) for m in models) == 2
        @test prob.prob.u0 == [12, 0]
        @test prob.prob.p == [0.7, 0.3]           # k_birth once, then gamma_peer
        ctxs = InferCell._build_contexts(models)
        @test ctxs[1].state_idxs == 1:1 && ctxs[2].state_idxs == 2:2
        @test ctxs[1].param_idxs == [1]
        @test ctxs[2].param_idxs == [2, 1]        # gamma_peer, then the shared k_birth
        # The peer's local p is [gamma, k_birth] = [0.3, 0.7] whichever slot each sits in.
        sl = InferCell._JumpSlot(SVector(2), SVector(2, 1), SVector(1),
                                 SVector(true), SVector(:X), :PeerBirthDeath)
        death, birth = reactions(peer)
        @test InferCell._global_jump(death, sl).rate([12, 0], [0.7, 0.3], 0.0) == 0.3 * 12
        @test InferCell._global_jump(birth, sl).rate([12, 0], [0.7, 0.3], 0.0) == 0.7
        # Reversing the composition order moves the slices, not the values.
        prob_r = build_problem([peer, owner]; tspan = (0.0, 1.0))
        @test prob_r.prob.u0 == [0, 12]
        @test prob_r.prob.p == [0.3, 0.7]
    end

    # Task 2.4: a propensity proportional to a peer's count. The owner is
    # frozen (birth rate zero), the reader fires at c · X and counts its own
    # firings, so the count over [0, T] is Poisson(c · X0 · T). Doubling X0
    # must double the mean; the tolerance is four standard errors of the ratio
    # of two Poisson sample means, not a round number.
    @testset "2.4 a jump module reads a declared peer's state in a propensity" begin
        c, T, n = 0.1, 10.0, 200
        function mean_fired(X0)
            prob = build_problem([BirthOwner(k = 0.0, x0 = X0), PeerReader(c = c)]; tspan = (0.0, T))
            Random.seed!(7)
            fired = zeros(n)
            peer_untouched = true
            for i in 1:n
                sol = solve(prob, SSAStepper(); saveat = [T])
                peer_untouched &= sol[1, end] == X0   # the reader never writes its peer
                fired[i] = sol[2, end]
            end
            @test peer_untouched
            return mean(fired)
        end
        m1, m2 = mean_fired(20), mean_fired(40)
        se_ratio = 2 * sqrt(1 / (n * 20) + 1 / (n * 40))   # ≈ 0.039
        @test abs(m1 - c * 20 * T) < 4 * sqrt(c * 20 * T / n)
        @test abs(m2 / m1 - 2) < 4 * se_ratio
        @info "task 2.4 peer read" mean_at_20 = m1 mean_at_40 = m2 ratio = m2 / m1 four_se = 4 * se_ratio
    end

    # Task 2.5: a declared peer write. The owner is frozen at 50 and the
    # decay-shaped peer removes one X per firing through the peer view, so X
    # runs to zero and the peer's counter reaches exactly 50 (the chance a
    # single copy survives 1000 s at rate 0.2 is e^-200). Remove the
    # declaration and the very same write throws at its first firing.
    @testset "2.5 a jump module writes a declared peer's state" begin
        owner = BirthOwner(k = 0.0, x0 = 50)
        decay = PeerBirthDeath(gamma = 0.2, b = 0.0, writes = [:X])
        prob = build_problem([owner, decay]; tspan = (0.0, 1000.0))
        Random.seed!(11)
        sol = solve(prob, SSAStepper(); saveat = [1000.0])
        @test sol[1, end] == 0
        @test sol[2, end] == 50
        @test all(x -> x >= 0, sol[1, :])   # never fires at zero copies: rate is gamma · X

        undeclared = PeerBirthDeath(gamma = 0.2, b = 0.0, writes = Symbol[])
        prob_u = build_problem([owner, undeclared]; tspan = (0.0, 1000.0))
        err = caught(() -> solve(prob_u, SSAStepper(); saveat = [1000.0]))
        @test err isa ArgumentError
        @test occursin("PeerBirthDeath", err.msg)
        @test occursin(":X", err.msg)
        @test occursin("written_states", err.msg)
    end

    @testset "2.5 the declaration is held to inputs() and, for registry species, to an edge" begin
        # A write to a state not read through inputs() has no view to go through.
        err = caught(() -> build_problem([BirthOwner(), PeerBirthDeath(writes = [:fired])]))
        @test err isa ArgumentError
        @test occursin("PeerBirthDeath", err.msg) && occursin(":fired", err.msg)
        # An ODE module has no reactions to write through; it contributes instead.
        err = caught(() -> build_problem([OdeWithWrites()]))
        @test err isa ArgumentError
        @test occursin("OdeWithWrites", err.msg) && occursin("contributed_states", err.msg)
        # A registry species written by a jump module needs a mass or currency
        # edge declaring the crossing; a transcript, being outside the registry,
        # is gated by the declaration alone (spec §12, 2026-09-04).
        err = caught(() -> resolve_coupling([RegistryJumpWriter(edges = CouplingEdge[])]))
        @test err isa ArgumentError
        @test occursin("RegistryJumpWriter", err.msg) && occursin(":M_atp_c", err.msg)
        @test occursin("written_states", err.msg)
        @test resolve_coupling([RegistryJumpWriter()]) isa CouplingGraph
    end
end
