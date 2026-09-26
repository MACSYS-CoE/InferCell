"""
The balance checks on the assembled Core A′ (spec §3, §11 phase 14a).

Every conservation check in the module phases ran on one ODE solve. The
assembled model is a hybrid that restarts its ODE block at every one of 6,300
handshakes and rewrites pools from whole particles, so spec §3 requires each
inherited check to be **restated in particles and to carry the `N_restarts`
factor**. This file is that restatement, written once, so a module-local check
and its composition-scope version are one implementation with one bound.

**What a closure counts.** A moiety is a weighted sum over the whole cell, not
over the ODE block:

```
W = Σᵢ wᵢ·(uᵢ·factor + remᵢ)          ODE pools, in particles, with the carried remainder
  + Σⱼ τⱼ·nⱼ                           jump states: transcripts by base content, proteins
  + Σ_rows κ(counter)·sign·exchanged   what a registry chemostat supplied or absorbed
  + Σ_c e_c·(pending_c + deficit_c)    accrual the hook has not yet applied
```

Every term but the first is an exact integer, read from state the driver
already holds, so nothing is integrated by quadrature. The transcript weights
come from the **extract**, not from the module, for the same reason task 11.8's
residue closure does: a module whose counters charge the wrong base count
cannot then close its own check. `e_c` is what applying one unit of counter
`c` would add to `W` across all its rows. Adding it for the accrual still in
flight makes `W` exact at every handshake, not only at a drain.
"""

# ---------------------------------------------------------------------------
# The flux meter: carbon balance needs the glucose that crossed in
# ---------------------------------------------------------------------------

"The ODE state integrating glucose uptake, in mM of the cell."
const GLC_UPTAKE_METER = :glc_uptake_meter

"The ODE state integrating lactate export, in mM of the cell."
const LAC_EXPORT_METER = :lac_export_meter

"""
    MeteredPtsTransport(inner = PtsTransport())

`PtsTransport` with two accumulator states appended: [`GLC_UPTAKE_METER`](@ref)
integrates GLCpts4's net rate and [`LAC_EXPORT_METER`](@ref) the exporter's.
Spec §3's check 2 is "glucose in against lactate out", and glucose enters
through a clamped external pool, so no state records how much has crossed. The
meter is that record.

The meters are non-registry states that nothing reads. They are intracellular,
so growth dilutes them at constant particle count, and `uᵢ·factor` is then the
particles that crossed, at whatever volume each crossed at. External lactate
cannot serve as the export record: it is referred to the medium at `1/R` and
is not diluted, so its particle count depends on which volume the medium is
taken to be.

**What it changes.** Nothing in any rate law. The adaptive step controller
sees two more states, so a metered trajectory differs from an unmetered one at
the level of the integrator tolerance. It is a measurement instrument, used by
the phase 14a checks, and `build_corea()` does not compose it unless asked.
"""
struct MeteredPtsTransport{P <: PtsTransport} <: AbstractSubModel
    inner::P
end
MeteredPtsTransport() = MeteredPtsTransport(PtsTransport())

for f in (:parameters, :inputs, :contributed_states, :coupling,
          :membrane_protein_states, :extracellular_states, :reduction_notes,
          :formalism, :inference_mode)
    @eval $f(m::MeteredPtsTransport) = $f(m.inner)
end
module_id(::MeteredPtsTransport) = :PtsTransport
states(m::MeteredPtsTransport) = vcat(states(m.inner), [GLC_UPTAKE_METER, LAC_EXPORT_METER])

const _N_PTS_STATES = 9
@inline _pts_own(u) = SVector{_N_PTS_STATES}(ntuple(i -> u[i], _N_PTS_STATES))

function dynamics(u, p, t, m::MeteredPtsTransport, u_inputs)
    own = _pts_own(u)
    du = dynamics(own, p, t, m.inner, u_inputs)
    _, _, _, _, v4, v_export = _pts_rates(_live_scalars(m.inner, p), own, u_inputs)
    return vcat(du, SVector(v4, v_export))
end

contributions(u, p, t, m::MeteredPtsTransport, u_inputs) =
    contributions(_pts_own(u), p, t, m.inner, u_inputs)

# ---------------------------------------------------------------------------
# Moieties
# ---------------------------------------------------------------------------

"""
    Moiety(name, ode; jump = [], chemostat = [])

One conserved quantity across the whole cell:
- `ode`, weights on ODE states;
- `jump`, weights on jump states;
- `chemostat`, the weight per unit a counter moves through a registry
  chemostat. That is what the counter **carries**, not the chemostat species.
  `CTP_mRNA` supplies CTP (three phosphates) and `CMP_mRNAdeg` returns CMP
  (one) into the same `M_ctp_c` chemostat.

Every term is in particles.
"""
struct Moiety
    name::Symbol
    ode::Vector{Pair{Symbol, Float64}}
    jump::Vector{Pair{Symbol, Float64}}
    chemostat::Vector{Pair{Symbol, Float64}}
end
Moiety(name::Symbol, ode; jump = Pair{Symbol, Float64}[],
       chemostat = Pair{Symbol, Float64}[]) =
    Moiety(name, _pairs(ode), _pairs(jump), _pairs(chemostat))
_pairs(v) = Pair{Symbol, Float64}[Symbol(k) => Float64(w) for (k, w) in v]

# The carbon in each glycolytic intermediate. The two meters enter with the sign
# that makes the sum conserved: glucose that crossed in is carbon the cell did
# not start with, and lactate that crossed out is carbon it no longer holds.
const _CARBON = (M_g6p_c = 6, M_f6p_c = 6, M_fdp_c = 6, M_dhap_c = 3,
                 M_g3p_c = 3, M_13dpg_c = 3, M_3pg_c = 3, M_2pg_c = 3,
                 M_pep_c = 3, M_pyr_c = 3, M_lac__L_c = 3)

# Every phosphorylated species, pyrophosphate twice. The phospho-carriers carry
# one each. NAD⁺ and NADH are left out: they carry a pyrophosphate, but GAPD
# and LDH only move a hydride between them, and redox is checked on its own.
const _PHOSPHATE = (M_atp_c = 3, M_adp_c = 2, M_amp_c = 1, M_pi_c = 1, M_ppi_c = 2,
                    M_gtp_c = 3, M_gdp_c = 2, M_gmp_c = 1,
                    M_g6p_c = 1, M_f6p_c = 1, M_fdp_c = 2, M_dhap_c = 1,
                    M_g3p_c = 1, M_13dpg_c = 2, M_3pg_c = 1, M_2pg_c = 1,
                    M_pep_c = 1,
                    M_ptsi_P_c = 1, M_ptsh_P_c = 1, M_crr_P_c = 1, M_ptsg_P_c = 1)

# What each chemostat row carries, in phosphates: a supplied NTP carries three
# and a returned NMP one (spec §11 task 14.6, annotated 2026-09-24).
const _CHEMOSTAT_PHOSPHATE = (CTP_mRNA = 3, UTP_mRNA = 3, CMP_mRNAdeg = 1, UMP_mRNAdeg = 1)

# The four carriers, each with the counter translation credits it through.
const _CARRIERS = ((name = :ptsI, unphos = :M_ptsi_c, phos = :M_ptsi_P_c, counter = :ptsI_translat),
                   (name = :ptsH, unphos = :M_ptsh_c, phos = :M_ptsh_P_c, counter = :ptsH_translat),
                   (name = :Crr,  unphos = :M_crr_c,  phos = :M_crr_P_c,  counter = :Crr_translat),
                   (name = :ptsG, unphos = :M_ptsg_c, phos = :M_ptsg_P_c, counter = :ptsG_transloc))

"""
    corea_moieties(; genes = read_transcription_genes(), carbon = true) -> Vector{Moiety}

The conserved quantities of the balance checks (spec §3), in check order:

| Moiety | Check | Jump side |
|---|---|---|
| `:carbon` | 2 | none; the two meters close it |
| `:redox` | 3 | none |
| `:adenylate`, `:guanylate` | 4 | each transcript's A or G count |
| `:phosphate` | 4b | each transcript's length, and CTP/UTP chemostat rows |
| `:carrier_ptsI` … `:carrier_ptsG` | 5 | minus the protein count translation credits from |

Carbon needs [`MeteredPtsTransport`](@ref) in the composition. Pass
`carbon = false` for a run without it.

A carrier's protein count enters with weight −1. Translation credits one new
carrier per protein made, so conservation "up to what translation adds" is a
subtraction, not a wider bound.
"""
function corea_moieties(; genes = read_transcription_genes(), carbon::Bool = true)
    tx(f) = [transcript_state(g.locus) => f(g) for g in genes]
    out = Moiety[]
    carbon && push!(out, Moiety(:carbon,
        vcat(collect(pairs(_CARBON)), [GLC_UPTAKE_METER => -6, LAC_EXPORT_METER => 3])))
    push!(out, Moiety(:redox, (:M_nad_c => 1, :M_nadh_c => 1)))
    push!(out, Moiety(:adenylate, (:M_atp_c => 1, :M_adp_c => 1, :M_amp_c => 1);
                      jump = tx(g -> g.counts.A)))
    push!(out, Moiety(:guanylate, (:M_gtp_c => 1, :M_gdp_c => 1, :M_gmp_c => 1);
                      jump = tx(g -> g.counts.G)))
    push!(out, Moiety(:phosphate, collect(pairs(_PHOSPHATE));
                      jump = tx(g -> g.length),
                      chemostat = collect(pairs(_CHEMOSTAT_PHOSPHATE))))
    for c in _CARRIERS
        spec = only(r for r in TRANSLATION_COUNTERS if r.counter === c.counter)
        push!(out, Moiety(Symbol(:carrier_, c.name), (c.unphos => 1, c.phos => 1);
                          jump = [protein_state(spec.locus) => -1]))
    end
    return out
end

# ---------------------------------------------------------------------------
# A moiety resolved against one driver
# ---------------------------------------------------------------------------

# Indices resolved once, so the per-handshake total touches no symbol.
struct _MoietyLayout
    moiety::Moiety
    ode::Vector{Tuple{Int, Float64}}         # (ODE index, weight)
    jump::Vector{Tuple{Int, Float64}}        # (jump index, weight)
    chemostat::Vector{Tuple{Int, Float64}}   # (debit index, carried weight)
    inflight::Vector{Tuple{Int, Vector{Int}, Float64}}  # (counter slot, consumer debits, e_c)
end

# The composed state names of one block, in the builders' layout.
_block_names(models, f) = Symbol[s for m in models if formalism(m) === f for s in states(m)]

function _layout(m::Moiety, ode_names, jump_names, d::HandshakeDriver)
    find(names, s, block) = something(findfirst(==(s), names), 0) |> i ->
        i > 0 ? i : throw(ArgumentError(
            "Moiety :$(m.name) weights :$s, which is not a state of this " *
            "composition's $block block. Compose the module that owns it" *
            (s in (GLC_UPTAKE_METER, LAC_EXPORT_METER) ?
             " — the carbon meters need MeteredPtsTransport, `build_corea(metered = true)`" : "")))
    ode = [(find(ode_names, s, "ODE"), w) for (s, w) in m.ode]
    jump = [(find(jump_names, s, "jump"), w) for (s, w) in m.jump]
    wode = Dict(m.ode)
    wchem = Dict(m.chemostat)
    chem = [(k, get(wchem, b.counter, 0.0))
            for (k, b) in enumerate(d.debits) if b.pool_idx == 0]
    inflight = Tuple{Int, Vector{Int}, Float64}[]
    for g in d.counters
        e = 0.0
        consumers = Int[]
        for k in g.debit_idxs
            b = d.debits[k]
            w = b.pool_idx == 0 ? get(wchem, b.counter, 0.0) :
                get(wode, ode_names[b.pool_idx], 0.0)
            e += w * b.sign * (b.sign < 0 ? 1.0 : b.stoich)
            b.sign < 0 && push!(consumers, k)
        end
        # With two consumers each pays the whole accrual, so pending + deficit
        # would count the accrual twice. No Core A′ counter has two, and the
        # driver refuses products on one (spec §12, 2026-09-24 E).
        length(consumers) <= 1 || throw(ArgumentError(
            "Counter :$(d.debits[first(consumers)].counter) has " *
            "$(length(consumers)) consumer rows, so its in-flight accrual has no " *
            "single meaning for moiety :$(m.name)"))
        e == 0.0 || push!(inflight, (g.counter_idx, consumers, e))
    end
    return _MoietyLayout(m, ode, jump, chem, inflight)
end

function _total(L::_MoietyLayout, d::HandshakeDriver)
    s = 0.0
    u = d.ode.u
    rem = d.rounding.remainders
    for (i, w) in L.ode
        s += w * (u[i] * d.factor + rem[i])
    end
    for (j, w) in L.jump
        s += w * d.jump.u[j]
    end
    for (k, w) in L.chemostat
        b = d.debits[k]
        s += w * b.sign * b.exchanged
    end
    for (c, consumers, e) in L.inflight
        flight = float(d.jump.u[c])
        for k in consumers
            flight += d.debits[k].deficit
        end
        s += e * flight
    end
    return s
end

# ---------------------------------------------------------------------------
# One recorded run
# ---------------------------------------------------------------------------

"""
    ValidationRun

What [`validation_run!`](@ref) records over a trajectory: everything the phase
14a checks read. The fields are:
- `t`, the handshake times;
- `ode`, the ODE state after each handshake, including the initial state, and
  `factor`, the particles per mM it was held at;
- `jump_min`, each jump state's minimum;
- `totals`, each moiety's whole-cell total at every handshake, with index 1
  the initial state;
- `max_remainder`, the largest carried remainder seen;
- `samples`, `(t, u, p)` every `sample_every` handshakes, for the right-hand
  side gates;
- `factor_max`;
- the solver's `abstol` and `reltol`;
- the rounding `policy`.
"""
struct ValidationRun
    ode_names::Vector{Symbol}
    jump_names::Vector{Symbol}
    moieties::Vector{Moiety}
    t::Vector{Float64}
    ode::Vector{Vector{Float64}}
    factor::Vector{Float64}
    jump_min::Vector{Int}
    totals::Dict{Symbol, Vector{Float64}}
    max_remainder::Float64
    samples::Vector{NamedTuple{(:t, :u, :p), Tuple{Float64, Vector{Float64}, Vector{Float64}}}}
    factor_max::Float64
    abstol::Float64
    reltol::Float64
    policy::Symbol
    n_restarts::Int
end

"""
    validation_run!(models, d, n; moieties = corea_moieties(...), sample_every = 60)
        -> ValidationRun

Advance the driver `n` handshakes, recording what the balance checks read.
`models` must be the list `d` was built from, since the driver keeps indices
and not names. Each moiety's total is taken at every handshake, so a closure is
asserted along the whole trajectory, not only at its end.
"""
function validation_run!(models::Vector{<:AbstractSubModel}, d::HandshakeDriver,
                         n::Integer;
                         moieties = corea_moieties(carbon = GLC_UPTAKE_METER in
                                                   _block_names(models, :ode)),
                         sample_every::Integer = 60)
    ode_names = _block_names(models, :ode)
    jump_names = _block_names(models, :jump)
    length(ode_names) == length(d.ode.u) && length(jump_names) == length(d.jump.u) ||
        throw(ArgumentError(
            "These models name $(length(ode_names)) ODE and $(length(jump_names)) jump " *
            "states, but the driver holds $(length(d.ode.u)) and $(length(d.jump.u)). " *
            "Pass the model list the driver was built from"))
    layouts = [_layout(m, ode_names, jump_names, d) for m in moieties]
    totals = Dict(L.moiety.name => [_total(L, d)] for L in layouts)
    t = [d.ode.t]
    ode = [collect(Float64, d.ode.u)]
    factor = [d.factor]
    jmin = collect(Int, d.jump.u)
    maxrem = maximum(abs, d.rounding.remainders; init = 0.0)
    samples = [(t = d.ode.t, u = collect(Float64, d.ode.u), p = collect(Float64, d.ode.p))]
    for k in 1:n
        handshake_step!(d)
        push!(t, d.ode.t)
        push!(ode, collect(Float64, d.ode.u))
        push!(factor, d.factor)
        for (j, x) in enumerate(d.jump.u)
            x < jmin[j] && (jmin[j] = x)
        end
        maxrem = max(maxrem, maximum(abs, d.rounding.remainders; init = 0.0))
        for L in layouts
            push!(totals[L.moiety.name], _total(L, d))
        end
        k % sample_every == 0 &&
            push!(samples, (t = d.ode.t, u = collect(Float64, d.ode.u),
                            p = collect(Float64, d.ode.p)))
    end
    s = solver_settings(d)
    return ValidationRun(ode_names, jump_names, collect(Moiety, moieties), t, ode, factor,
                         jmin, totals, maxrem, samples, maximum(factor), s.abstol, s.reltol,
                         d.rounding.policy, n)
end

# ---------------------------------------------------------------------------
# The checks
# ---------------------------------------------------------------------------

"""
    closure_residual(run, name) -> Vector{Float64}

A moiety's whole-cell total minus its initial value, in particles, at every
handshake.
"""
function closure_residual(run::ValidationRun, name::Symbol)
    haskey(run.totals, name) || throw(ArgumentError(
        "No moiety :$name in this run; it recorded $(sort(collect(keys(run.totals))))"))
    tot = run.totals[name]
    return tot .- first(tot)
end

"""
    state_bound(us, i, abstol, reltol) -> Float64
    state_bound(run::ValidationRun, i) -> Float64

The integrator's bound on state `i` over a trajectory `us`, a vector of state
vectors: `max(abstol, reltol · maxₜ|xᵢ|)`, in the state's own unit. It is the
single definition behind check 1's non-negativity and behind `tol_C`, used
module-locally (phase 6 passes `sol.u`) and at composition scope alike.
"""
state_bound(us::AbstractVector{<:AbstractVector}, i::Integer, abstol, reltol) =
    max(abstol, reltol * maximum(abs(u[i]) for u in us))
state_bound(run::ValidationRun, i::Integer) = state_bound(run.ode, i, run.abstol, run.reltol)

"""
    conservation_bound(us, weights; abstol, reltol, n_restarts = 1, factor = 1.0) -> Float64

Spec §3's single-run bound on the conserved sum `Σ wᵢ·xᵢ`, where `weights` is
a list of `(state index, wᵢ)` pairs:

```
tol_C = N_restarts · Σᵢ |wᵢ| · max(abstol, reltol · maxₜ|xᵢ|)
```

A module-local solve has no restarts and is in mM. At composition scope the
bound carries `N_restarts` and is converted to particles by `factor`, as §3
requires of an inherited check. One implementation, with two scopes.
"""
conservation_bound(us::AbstractVector{<:AbstractVector}, weights; abstol, reltol,
                   n_restarts::Integer = 1, factor::Real = 1.0) =
    n_restarts * factor * sum(abs(w) * state_bound(us, i, abstol, reltol) for (i, w) in weights)

"""
    moiety_drift(us, weights) -> Float64

The largest departure of `Σ wᵢ·xᵢ` from its initial value over a trajectory:
the module-local residual, in the states' own unit.
"""
function moiety_drift(us::AbstractVector{<:AbstractVector}, weights)
    total(u) = sum(w * u[i] for (i, w) in weights)
    t0 = total(first(us))
    return maximum(abs(total(u) - t0) for u in us)
end

"""
    assert_conserved(label, drift, bound; unit = "mM") -> drift
    assert_conserved(run, name; bound) -> Float64

Throw an `ArgumentError` naming the conserved quantity if `drift` exceeds
`bound`. Otherwise return the drift. The message is the assertion: a mutation
test has to see *which* quantity broke. The run form checks a moiety's
whole-cell closure at every handshake, in particles.
"""
function assert_conserved(label::AbstractString, drift::Real, bound::Real;
                          unit::AbstractString = "mM", at::AbstractString = "")
    drift <= bound && return drift
    throw(ArgumentError(
        "The conserved quantity $label is not invariant: it drifts by $drift " *
        "$unit$at, against an integrator bound of $bound $unit"))
end

_weights(run::ValidationRun, m::Moiety) =
    [(findfirst(==(s), run.ode_names), w) for (s, w) in m.ode]
_moiety(run::ValidationRun, name::Symbol) = only(x for x in run.moieties if x.name === name)

"""
    moiety_bound(run, name) -> Float64

[`conservation_bound`](@ref) at composition scope, in particles: the run's
handshake count as `N_restarts` and its largest particles-per-mM as `factor`.
Only the ODE weights count toward it. The jump, chemostat and in-flight terms
are exact integers and add nothing.
"""
moiety_bound(run::ValidationRun, name::Symbol) =
    conservation_bound(run.ode, _weights(run, _moiety(run, name));
                       abstol = run.abstol, reltol = run.reltol,
                       n_restarts = run.n_restarts, factor = run.factor_max)

function assert_conserved(run::ValidationRun, name::Symbol; bound::Real)
    r = closure_residual(run, name)
    k = argmax(abs.(r))
    return assert_conserved("the $name moiety, across the assembled cell,", abs(r[k]),
                            bound; unit = "particles", at = " at t = $(run.t[k]) s")
end

"""
    rhs_gate(d, run, name) -> NamedTuple

Which branch of spec §3's exact-conservation exception the composed right-hand
side passes for one moiety's ODE weights. The composed `f` is evaluated at every
sampled state:
- `:bitwise` if `Σ wᵢ·duᵢ === 0.0` at every sample;
- `:one_ulp` if it is at most one ulp of the conserved sum `Σ wᵢ·uᵢ`;
- `:n_ulps` if it is at most `n` ulps, where `n` is the number of weighted
  terms (amended 2026-09-25: a sum of 21 terms arriving from several modules
  rounds more than a sum of three);
- `:none` otherwise, in which case the fivefold fall applies.

Also returned:
- `n`, the number of weighted terms;
- `max_ulps`, the largest per-evaluation residual in ulps of the sum;
- `max_gamma`, the largest residual over the floating-point summation bound
  `γₙ = n·ε·Σ|wᵢ·duᵢ|`, where `n` is the number of weighted terms. At most 1
  means the residual is within what rounding in the sum alone can produce.

The criterion is mechanical: which gate a check uses is read off this, not
argued (spec §3).
"""
function rhs_gate(d::HandshakeDriver, run::ValidationRun, name::Symbol)
    idx = _weights(run, _moiety(run, name))
    f = d.ode.f
    bitwise = true
    worst = 0.0
    gamma = 0.0
    for smp in run.samples
        u = SVector{length(smp.u)}(smp.u)
        du = f(u, smp.p, smp.t)
        r = sum(w * du[i] for (i, w) in idx)
        total = sum(w * u[i] for (i, w) in idx)
        scale = length(idx) * eps() * sum(abs(w * du[i]) for (i, w) in idx)
        bitwise &= r === 0.0 || r === -0.0
        worst = max(worst, abs(r) / eps(abs(total)))
        scale > 0 && (gamma = max(gamma, abs(r) / scale))
    end
    n = length(idx)
    gate = bitwise ? :bitwise : worst <= 1.0 ? :one_ulp : worst <= n ? :n_ulps : :none
    return (gate = gate, n = n, max_ulps = worst, max_gamma = gamma)
end

"""
    first_negative(run) -> Union{Nothing, NamedTuple}

Check 1: the first ODE state, in time, that falls below the negative of its own
integrator bound (`state_bound`). Returns `(species, t, value, bound)`, or
`nothing`. A jump state below zero is reported the same way, with `bound = 0`.
Naming the first violation is the point. A global failure would say nothing
about where the model went wrong.
"""
function first_negative(run::ValidationRun)
    bounds = [state_bound(run, i) for i in eachindex(run.ode_names)]
    for (k, u) in enumerate(run.ode), i in eachindex(u)
        u[i] < -bounds[i] &&
            return (species = run.ode_names[i], t = run.t[k], value = u[i], bound = bounds[i])
    end
    j = findfirst(<(0), run.jump_min)
    j === nothing || return (species = run.jump_names[j], t = NaN,
                             value = float(run.jump_min[j]), bound = 0.0)
    return nothing
end

"""
    assert_nonnegative(run)

Throw an `ArgumentError` naming the first state and time [`first_negative`](@ref)
finds. Otherwise return `nothing`.
"""
function assert_nonnegative(run::ValidationRun)
    v = first_negative(run)
    v === nothing && return nothing
    throw(ArgumentError(
        "Check 1 fails: :$(v.species) is $(v.value) at t = $(v.t) s, below the " *
        "negative of its integrator bound, $(-v.bound)"))
end

"""
    carbon_accounts(run) -> NamedTuple

Check 2's scalars over the run, in particles:
- `glucose_in` and `lactate_out`, the two meters;
- `lactate_formed`, which is lactate out plus the rise in cytosolic lactate;
- `glucose_consumed`, which is glucose in minus the rise in the glycolytic
  intermediates, in glucose equivalents (carbon over six);
- `homolactic`, `lactate_formed / glucose_consumed`. It is exactly 2 whenever
  carbon closes, because it is the carbon closure at the run's two endpoints
  written as a ratio. It is not independent evidence that the carbon went to
  lactate, since a rise in pyruvate is subtracted from the glucose consumed;
- `exported`, the fraction of the lactate formed that left the cell. This is
  the independent quantity: it is what fails when export is removed.
"""
function carbon_accounts(run::ValidationRun)
    i(s) = findfirst(==(s), run.ode_names)
    # Particles at the factor each state was held at. The carbon species are
    # never written by a debit, so they carry no rounding remainder.
    p(k, s) = run.ode[k][i(s)] * run.factor[k]
    first_, last_ = 1, length(run.ode)
    inter(k) = sum(w * p(k, s) for (s, w) in pairs(_CARBON) if s !== :M_lac__L_c) / 6
    gin = p(last_, GLC_UPTAKE_METER) - p(first_, GLC_UPTAKE_METER)
    lout = p(last_, LAC_EXPORT_METER) - p(first_, LAC_EXPORT_METER)
    formed = lout + p(last_, :M_lac__L_c) - p(first_, :M_lac__L_c)
    consumed = gin - (inter(last_) - inter(first_))
    return (glucose_in = gin, lactate_out = lout, lactate_formed = formed,
            glucose_consumed = consumed, homolactic = formed / consumed,
            exported = lout / formed)
end

# ---------------------------------------------------------------------------
# Check 1b: the particle floor (spec §3, §11 task 14.2)
# ---------------------------------------------------------------------------

"The particle count below which spec §3 check 1b flags a continuous state."
const PARTICLE_FLOOR = 500

"""
    particle_floor(run; floor = PARTICLE_FLOOR) -> Vector{NamedTuple}

Check 1b's report: every ODE state over the run, in particles at the factor it
was held at, sorted from the smallest minimum. Each row carries:
- `min_particles` and the time `t` it was reached;
- `median_particles` over the run;
- `below_one`, the fraction of handshakes spent under one particle;
- `flagged`, whether the minimum is below `floor`.

Every state is listed, as §3 asks. The carried rounding remainder is below one
particle and is not added. External states are referred to the medium rather
than the cell, so [`langevin_pools`](@ref) leaves them out.
"""
function particle_floor(run::ValidationRun; floor::Real = PARTICLE_FLOOR)
    rows = map(eachindex(run.ode_names)) do i
        x = [run.ode[k][i] * run.factor[k] for k in eachindex(run.ode)]
        k = argmin(x)
        (species = run.ode_names[i], min_particles = x[k], t = run.t[k],
         median_particles = median(x), below_one = count(<(1), x) / length(x),
         flagged = x[k] < floor)
    end
    return sort(rows; by = r -> r.min_particles)
end

"""
    flagged_states(report) -> Vector{Symbol}

The states [`particle_floor`](@ref) puts below the floor, smallest first.
"""
flagged_states(report) = Symbol[r.species for r in report if r.flagged]

"The cycle median, in particles, below which check 1b excludes a pool outright."
const CONTINUUM_MEDIAN = 10

"""
    langevin_pools(report; n = 3, min_median = CONTINUUM_MEDIAN, external = Symbol[])
        -> (cross_check, excluded)

How check 1b splits the flagged pools (spec §3, amended 2026-09-26):
- `excluded` are the flagged pools whose cycle median is under `min_median`
  particles. Neither the ODE nor a Langevin diffusion describes a pool of about
  one particle, so no ensemble is needed to exclude them.
- `cross_check` are the `n` remaining flagged pools with the smallest medians,
  which the chemical-Langevin ensemble is run on. By median, not minimum: a
  pool of thousands that a clipped debit empties for one handshake has a
  minimum of zero and is not a small pool.

`external` states are left out of both.
"""
function langevin_pools(report; n::Integer = 3, min_median::Real = CONTINUUM_MEDIAN,
                        external = Symbol[])
    rows = [r for r in report if r.flagged && !(r.species in external)]
    excluded = Symbol[r.species for r in rows if r.median_particles < min_median]
    rest = sort([r for r in rows if r.median_particles >= min_median]; by = r -> r.median_particles)
    return (cross_check = Symbol[r.species for r in rest[1:min(n, end)]], excluded = excluded)
end

# ---------------------------------------------------------------------------
# Check 7: which counter clips, and when (spec §3, §11 task 14.8)
# ---------------------------------------------------------------------------

"""
    ClipRecord(d)

Check 7's per-counter record. [`clipping_census`](@ref) counts drains at which
*any* counter carried a deficit. K5 needs to know which counter, how often and
from when, so [`record_clips!`](@ref) is called after every
`handshake_step!` and reads the carried deficits at each drain.
"""
mutable struct ClipRecord
    clipped::Dict{Symbol, Int}          # drains at which the counter carried a deficit
    first_clip::Dict{Symbol, Float64}   # the first such drain's time
    max_deficit::Dict{Symbol, Float64}  # the largest deficit carried, in particles
    counters::Vector{Symbol}            # every consumer counter on a live pool
    drains::Int
    last_drain::Int                     # the driver's drain count when last read
end
function ClipRecord(d::HandshakeDriver)
    cs = unique(Symbol[b.counter for b in d.debits if b.sign < 0 && b.pool_idx != 0])
    return ClipRecord(Dict{Symbol, Int}(), Dict{Symbol, Float64}(),
                      Dict{Symbol, Float64}(), cs, 0, d.n_drains)
end

"""
    record_clips!(rec, d) -> rec

Read the driver after one `handshake_step!`. On a step that drained, every
consumer row carrying a deficit is counted against its counter. A deficit is
nonzero exactly when the pool paid less than was asked, since a paid debit
leaves `accrued − accrued`, which is zero.
"""
function record_clips!(rec::ClipRecord, d::HandshakeDriver)
    d.n_drains == rec.last_drain && return rec
    rec.last_drain = d.n_drains
    rec.drains += 1
    for b in d.debits
        (b.sign < 0 && b.pool_idx != 0 && b.deficit > 0) || continue
        rec.clipped[b.counter] = get(rec.clipped, b.counter, 0) + 1
        haskey(rec.first_clip, b.counter) || (rec.first_clip[b.counter] = d.ode.t)
        rec.max_deficit[b.counter] = max(get(rec.max_deficit, b.counter, 0.0), b.deficit)
    end
    return rec
end

"""
    clip_summary(rec) -> Vector{NamedTuple}

One row per counter that clipped: `counter`, `drains` clipped, `first` clip
time, and `max_deficit` in particles. Sorted by drains clipped, most first. An
empty vector is check 7's zero.
"""
clip_summary(rec::ClipRecord) =
    sort([(counter = c, drains = n, first = rec.first_clip[c], max_deficit = rec.max_deficit[c])
          for (c, n) in rec.clipped]; by = r -> -r.drains)

# ---------------------------------------------------------------------------
# Output F5: the two external comparisons (spec §3, §11 task 14b.5)
# ---------------------------------------------------------------------------

# Ranks with ties averaged, so a Spearman correlation needs no dependency.
function _ranks(x::AbstractVector)
    o = sortperm(x)
    r = similar(x, Float64)
    i = 1
    while i <= length(o)
        j = i
        while j < length(o) && x[o[j+1]] == x[o[i]]
            j += 1
        end
        r[o[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return r
end

"""
    spearman(x, y) -> Float64

Spearman's rank correlation, with tied ranks averaged.
"""
function spearman(x::AbstractVector, y::AbstractVector)
    length(x) == length(y) || throw(ArgumentError("lengths differ"))
    rx, ry = _ranks(x), _ranks(y)
    rx .-= mean(rx)
    ry .-= mean(ry)
    return sum(rx .* ry) / sqrt(sum(abs2, rx) * sum(abs2, ry))
end

"""
    transcript_comparison(predicted, measured; min_rho = 0.7, factor = 2.0, min_within = 15)
        -> NamedTuple

Spec §3's first external comparison: predicted transcript steady states
against measured mean counts, gene by gene. It passes when Spearman's
correlation is at least `min_rho` and at least `min_within` genes agree
within `factor`. Returns `rho`, `within`, `n`, the per-gene `ratio` and `pass`.
"""
function transcript_comparison(predicted::AbstractVector, measured::AbstractVector;
                               min_rho::Real = 0.7, factor::Real = 2.0,
                               min_within::Integer = 15)
    ratio = predicted ./ measured
    within = count(r -> 1 / factor <= r <= factor, ratio)
    rho = spearman(predicted, measured)
    return (rho = rho, within = within, n = length(ratio), ratio = ratio,
            pass = rho >= min_rho && within >= min_within)
end

"""
    fold_change_report(folds, lengths; band = (1.7, 2.3), bounds = (1.5, 3.0)) -> NamedTuple

Spec §3's second external comparison: protein fold change over one cycle. It
passes when the median is inside `band`, no gene is outside `bounds`, and the
least-squares slope of log fold change against log length is negative, as in
the published histogram, where long genes underproduce. Returns `median`,
`min`, `max`, `below` and `above` (genes outside `bounds`), `slope` and `pass`.
"""
function fold_change_report(folds::AbstractVector, lengths::AbstractVector;
                            band = (1.7, 2.3), bounds = (1.5, 3.0))
    x = log.(lengths)
    y = log.(folds)
    slope = sum((x .- mean(x)) .* (y .- mean(y))) / sum(abs2, x .- mean(x))
    med = median(folds)
    below = count(<(bounds[1]), folds)
    above = count(>(bounds[2]), folds)
    return (median = med, min = minimum(folds), max = maximum(folds), below = below,
            above = above, slope = slope,
            pass = band[1] <= med <= band[2] && below == 0 && above == 0 && slope < 0)
end

export GLC_UPTAKE_METER, LAC_EXPORT_METER, MeteredPtsTransport, Moiety, corea_moieties,
       ValidationRun, validation_run!, closure_residual, state_bound, conservation_bound,
       moiety_drift, assert_conserved, moiety_bound, rhs_gate, first_negative,
       assert_nonnegative, carbon_accounts
export PARTICLE_FLOOR, particle_floor, flagged_states, CONTINUUM_MEDIAN, langevin_pools
export ClipRecord, record_clips!, clip_summary
export spearman, transcript_comparison, fold_change_report
