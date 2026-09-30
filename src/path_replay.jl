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
    PathReplay(path, driver; density = false)

A cursor into `path`, positioned at the first firing at or after the driver's
current time, so a replay can start from a mid-run snapshot. Throws if the path
was recorded on a composition whose jumps are labelled differently.

With `density = true` the replay also accumulates the path's log-density under
the replayed rate constants (spec/phases/16-recovery.md task 16a.2): see
[`path_logdensity`](@ref). `counts[j]` is how many times reaction `j` fired, and
`exposure[j]` is its integrated propensity `∫ λ_j dt`, which are the sufficient
statistics block 2 reads.
"""
mutable struct PathReplay
    path::JumpPath
    next::Int
    density::Bool
    logsum::Float64
    counts::Vector{Int}
    exposure::Vector{Float64}
end

function PathReplay(path::JumpPath, d::HandshakeDriver; density::Bool = false)
    path.labels == d.events.labels || throw(ArgumentError(
        "This path was recorded on a composition whose $(length(path.labels)) jumps " *
        "are labelled differently from this driver's $(length(d.events.labels)); " *
        "replaying it would apply each firing to the wrong reaction"))
    n = length(path.labels)
    return PathReplay(path, searchsortedfirst(path.times, d.jump.t), density, 0.0,
                      zeros(Int, n), zeros(n))
end

"""
    path_logdensity(replay) -> Float64

The log-density of the replayed stretch of the path under the rate constants the
replay rebuilt:

```
log p(X | θ) = Σ_firings log λ_k(t⁻) − Σ_j ∫ λ_j dt
```

Every propensity is constant between firings, because the rate constants change
only at a rebuild. On the jump clock a rebuild lands at the start of the
interval of the handshake that runs it, 60m − 1 s, since `handshake_step!`
rebuilds before the jump step. So the integral is an exact sum over the segments
between firings and interval ends. θ enters through the rebuilt constants, which
read the replayed ODE pools, and through any propensity parameter directly.
Requires a replay built with `density = true`.
"""
function path_logdensity(r::PathReplay)
    r.density || throw(ArgumentError(
        "this replay was built without density = true, so it accumulated nothing"))
    return r.logsum - sum(r.exposure)
end

# Add every propensity's integral over [t0, t1) at the current state.
function _accrue_exposure!(r::PathReplay, d::HandshakeDriver, t0, t1)
    dt = t1 - t0
    dt > 0 || return nothing
    u, p = d.jump.u, d.jump.p
    rates = d.events.rates
    for j in eachindex(rates)
        r.exposure[j] += rates[j](u, p, t0) * dt
    end
    return nothing
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
    t = d.jump.t
    while r.next <= length(p.times) && p.times[r.next] < target
        k = p.reactions[r.next]
        if r.density
            τ = p.times[r.next]
            _accrue_exposure!(r, d, t, τ)
            r.logsum += log(d.events.rates[k](d.jump.u, d.jump.p, τ))
            r.counts[k] += 1
            t = τ
        end
        affects[k](d.jump)
        r.next += 1
    end
    r.density && _accrue_exposure!(r, d, t, target)
    d.jump.tprev = d.jump.t
    d.jump.t = target
    _close_handshake!(d, drained, clipped)
    return d
end

"""
    replay!(driver, path, n_steps; density = false) -> PathReplay

Run `n_steps` handshakes of [`replay_step!`](@ref) from the driver's current
time. Returns the cursor, so a caller can continue the same replay or read
[`path_logdensity`](@ref) from it.
"""
function replay!(d::HandshakeDriver, path::JumpPath, n_steps::Integer;
                 density::Bool = false)
    r = PathReplay(path, d; density)
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
       path_logdensity, snapshot, restore
