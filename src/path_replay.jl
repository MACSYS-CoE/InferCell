"""
Driving the hybrid from a recorded jump path, and snapshotting a driver
(spec/phases/16-recovery.md task 16a.1).

The conditional sampler of phase 16 holds a cell's stochastic path fixed and
asks what the ODE block does under it. A path and the parameters fix the whole
trajectory: fractional carry and the clamped drain are deterministic, and the
initial transcripts are fixed. So the ODE side is a replay: every handshake runs
exactly as [`handshake_step!`](@ref) does, except that the stochastic block's
interval is advanced by applying the recorded firings in order rather than by
the SSA.
"""

"""
    JumpPath

A recorded stochastic path: every firing's time and its index in the composed
jump list, in firing order, with the composition's jump labels
(`<module>_<local index>`) so a replay on a different composition is refused.
Built by [`recorded_path`](@ref).
"""
struct JumpPath
    times::Vector{Float64}
    reactions::Vector{Int}
    labels::Vector{Symbol}

    function JumpPath(times, reactions, labels)
        length(times) == length(reactions) || throw(DimensionMismatch(
            "$(length(times)) firing times for $(length(reactions)) reaction indices"))
        issorted(times) || throw(ArgumentError("a path's firing times must be in order"))
        all(k -> 1 <= k <= length(labels), reactions) || throw(ArgumentError(
            "a path names a reaction outside its $(length(labels)) labels"))
        new(collect(Float64, times), collect(Int, reactions), collect(Symbol, labels))
    end
end

Base.length(p::JumpPath) = length(p.times)

"""
    record_path!(driver) -> driver

Start recording every firing of the driver's stochastic block, discarding
anything recorded before. Recording draws no random numbers, so a recorded run
is the same run as an unrecorded one. Read the record with
[`recorded_path`](@ref).
"""
function record_path!(d::HandshakeDriver)
    e = d.events
    isempty(e.labels) && error(
        "This driver's stochastic block was built without an event log, so its " *
        "firings cannot be recorded")
    empty!(e.times)
    empty!(e.reactions)
    e.recording = true
    return d
end

"""
    recorded_path(driver) -> JumpPath

The firings recorded since [`record_path!`](@ref), as a [`JumpPath`](@ref).
"""
recorded_path(d::HandshakeDriver) = JumpPath(copy(d.events.times), copy(d.events.reactions),
                                             copy(d.events.labels))

"""
    PathReplay(path, driver)

A cursor into `path`, positioned at the first firing at or after the driver's
current time, so a replay can start from a mid-run snapshot. Throws if the path
was recorded on a composition whose jumps are labelled differently.
"""
mutable struct PathReplay
    path::JumpPath
    next::Int
end

function PathReplay(path::JumpPath, d::HandshakeDriver)
    path.labels == d.events.labels || throw(ArgumentError(
        "This path was recorded on a composition whose $(length(path.labels)) jumps " *
        "are labelled differently from this driver's $(length(d.events.labels)); " *
        "replaying it would apply each firing to the wrong reaction"))
    return PathReplay(path, searchsortedfirst(path.times, d.jump.t))
end

"""
    replay_step!(driver, replay) -> driver

One handshake in which the stochastic block is advanced by the recorded path.
Steps 0 to 3b are [`handshake_step!`](@ref)'s own. Then every recorded firing
before the end of the interval is applied, in order, through the same affect the
SSA calls, and the block's clock moves to the end of the interval. A firing
exactly at the boundary belongs to the next interval, as it does in the SSA,
which stops at the boundary before firing it.

The SSA is never called, so its cached propensities go stale. A driver that has
been replayed must not be stepped by [`handshake_step!`](@ref) afterwards.
"""
function replay_step!(d::HandshakeDriver, r::PathReplay)
    _, drained, clipped = _exchange!(d)
    target = d.jump.t + d.interval
    p = r.path
    affects = d.events.affects
    while r.next <= length(p.times) && p.times[r.next] < target
        affects[p.reactions[r.next]](d.jump)
        r.next += 1
    end
    d.jump.tprev = d.jump.t
    d.jump.t = target
    _close_handshake!(d, drained, clipped)
    return d
end

"""
    replay!(driver, path, n_steps) -> PathReplay

Run `n_steps` handshakes of [`replay_step!`](@ref) from the driver's current
time. Returns the cursor, so a caller can continue the same replay.
"""
function replay!(d::HandshakeDriver, path::JumpPath, n_steps::Integer)
    r = PathReplay(path, d)
    for _ in 1:n_steps
        replay_step!(d, r)
    end
    return r
end

"""
    snapshot(driver) -> HandshakeDriver
    restore(snapshot) -> HandshakeDriver

A particle is the driver's full mutable state (spec/phases/16-recovery.md D16.3).
A snapshot is a deep copy of the driver, and `restore` deep-copies it again, so
one snapshot can seed any number of particles and no particle shares a solver
object, a counter, a deficit, a remainder or the growth state with another.

**The random stream is not in the snapshot.** The SSA draws from the task's
global generator, which a deep copy does not copy. A replay draws nothing, so
this does not touch it. A particle that simulates forward must seed its own
stream.
"""
snapshot(d::HandshakeDriver) = deepcopy(d)
restore(s::HandshakeDriver) = deepcopy(s)
@doc (@doc snapshot) restore

export JumpPath, record_path!, recorded_path, PathReplay, replay_step!, replay!,
       snapshot, restore
