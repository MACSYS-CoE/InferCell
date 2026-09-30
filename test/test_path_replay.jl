# spec/phases/16-recovery.md task 16a.1: the published hybrid driven from a
# recorded jump path (V2), and snapshots that restore a driver exactly (V2b).
# The full-cycle V2 on three 15.7 dataset cells is dev/scripts/path_replay_16a1.jl.

using Random
using Test
using InferCell

# Everything a handshake can change, so "bitwise" means every field.
function _replay_state(d)
    return (t = (d.ode.t, d.jump.t), ode = copy(d.ode.u), ode_p = copy(d.ode.p),
            jump = copy(d.jump.u), jump_p = copy(d.jump.p),
            carry = copy(d.rounding.remainders),
            deficit = [b.deficit for b in d.debits],
            exchanged = [b.exchanged for b in d.debits],
            geometry = (d.factor, d.area_nm2, d.radius_nm, d.volume_litres),
            census = (d.n_handshakes, d.n_drains, d.n_clipped))
end

function _cell(truth; horizon)
    dm = d11_models()
    d = build_problem(dm; tspan = (0.0, horizon), complete = true)
    set_parameters!(d, dm, truth.values)
    return d
end

# Run the SSA for `n` handshakes, recording, and keep every handshake's state.
function _simulate(truth, seed; horizon, n, record = true)
    Random.seed!(seed)
    d = _cell(truth; horizon)
    Random.seed!(seed)
    record && record_path!(d)
    states = [begin handshake_step!(d); _replay_state(d) end for _ in 1:n]
    return d, states
end

function _replay(truth, path; horizon, n)
    d = _cell(truth; horizon)
    r = PathReplay(path, d)
    return d, [begin replay_step!(d, r); _replay_state(d) end for _ in 1:n]
end

@testset "16a.1: path replay and snapshots" begin
    truth = draw_truth(Xoshiro(1601), d11_models(); purpose = :recovery)
    H, N, SEED = 300.0, 300, 16011

    d_sim, sim = _simulate(truth, SEED; horizon = H, n = N)
    path = recorded_path(d_sim)

    @testset "recording is transparent" begin
        _, unrecorded = _simulate(truth, SEED; horizon = H, n = N, record = false)
        @test unrecorded == sim
        @test length(path) > 0
        @test all(l -> occursin("_", String(l)), path.labels)
    end

    @testset "V2: the replay reproduces the run bitwise at every handshake" begin
        _, rep = _replay(truth, path; horizon = H, n = N)
        @test rep == sim
        # Mutation: one translation event removed changes GTP.
        tl = findall(k -> startswith(String(path.labels[k]), "CoreATranslation_") &&
                          k != findlast(l -> startswith(String(l), "CoreATranslation_"),
                                        path.labels), path.reactions)
        @test !isempty(tl)
        drop = tl[cld(length(tl), 2)]
        keep = setdiff(eachindex(path.times), drop)
        mutant = JumpPath(path.times[keep], path.reactions[keep], path.labels)
        _, rep_m = _replay(truth, mutant; horizon = H, n = N)
        igtp = findfirst(==(:M_gtp_c), InferCell._block_names(d11_models(), :ode))
        @test any(k -> rep_m[k].ode[igtp] != sim[k].ode[igtp], 1:N)
        # A path from another composition is refused.
        wrong = JumpPath(path.times, path.reactions, reverse(path.labels))
        @test_throws ArgumentError PathReplay(wrong, _cell(truth; horizon = H))
    end

    @testset "V2b: a restored snapshot continues exactly" begin
        B = 120                                   # a window boundary
        d = _cell(truth; horizon = H)
        replay!(d, path, B)
        s = snapshot(d)
        # The intervening particle differs in every field a snapshot must carry.
        x = restore(s)
        record_path!(x)
        Random.seed!(99)
        for _ in 1:30
            handshake_step!(x)
        end
        for b in x.debits
            b.deficit = 7.0
        end
        x.rounding.remainders .+= 0.25
        x.factor *= 1.01
        x.jump.u .+= 1
        # Restored, the snapshot continues exactly as the uninterrupted replay.
        y = restore(s)
        r = PathReplay(path, y)
        cont = [begin replay_step!(y, r); _replay_state(y) end for _ in (B + 1):N]
        @test cont == sim[(B + 1):N]
        # Mutations: a restore that carried the other particle's deficit, or its
        # conversion factor, does not.
        for corrupt! in (z -> foreach(b -> b.deficit = 7.0, z.debits),
                         z -> (z.factor = x.factor))
            z = restore(s)
            corrupt!(z)
            rz = PathReplay(path, z)
            @test [begin replay_step!(z, rz); _replay_state(z) end for _ in (B + 1):N] !=
                  sim[(B + 1):N]
        end
    end
end
