"""
Core A′'s phosphotransferase transport module (spec §11 phase 7).

Five mass-action steps that import glucose against phosphoenolpyruvate, plus
passive lactate export. It owns the eight carrier phospho-states and external
lactate, and it is the module that makes carbon balance (spec §3, check 2)
posable at all.

The cascade's rate law is genuinely different from glycolysis's modular form,
and that is the source model's own split: intracellular reactions get
convenience kinetics, transport gets mass action. The carriers *are* the
reactants, so there is no enzyme multiplier and no saturation — nothing here
bounds a rate as a concentration rises. With the carriers conserved that is
harmless, because the pools cap themselves, and it is why the four conservation
sums rather than a flux ceiling are the invariant this module leans on.
"""

# Where this module's vendored extracts live. `src/organisms/coreA/data/README.md`
# records the upstream files, the commit and the command that regenerates them.
const COREA_DATA_DIR = joinpath(@__DIR__, "data")

"""
The scalars that reach [`dynamics`](@ref), in `parameters` order.

Every one is carried on the struct as well as in `parameters`, because a
parameter with `fixed = true` is absent from the composed parameter vector
entirely — [`model_free_params`](@ref) filters it out, so `p` never carries it
and a rate law cannot read it from there. The struct copy and the
[`InferParameter`](@ref) are built from the same `load_parameter` call, so they
cannot drift, and a test asserts they agree.
"""
const PTS_SCALARS = (:kf_glcpts0, :kr_glcpts0,
                     :kf_glcpts1, :kr_glcpts1,
                     :kf_glcpts2, :kr_glcpts2,
                     :kf_glcpts3, :kr_glcpts3,
                     :kf_glcpts4, :kr_glcpts4,
                     :p_lact2r,
                     :r_cell_nm, :lac_volume_ratio, :glc_e_mM)

"""
The eleven upstream identifiers, in [`PTS_SCALARS`](@ref) order, as
`src/organisms/coreA/data/pts_transport.tsv` names them.
"""
const N_PTS_SCALARS = length(PTS_SCALARS)

const PTS_RATE_IDS = ("KF_0_R_GLCpts0", "KR_0_R_GLCpts0",
                      "KF_1_R_GLCpts1", "KR_1_R_GLCpts1",
                      "KF_2_R_GLCpts2", "KR_2_R_GLCpts2",
                      "KF_3_R_GLCpts3", "KR_3_R_GLCpts3",
                      "KF_4_R_GLCpts4", "KR_4_R_GLCpts4",
                      "P_R_L_LACt2r")

"""
The nine states this module integrates, in registry order.

External lactate comes **first**: the registry groups it under `:other` at
position 22, ahead of the eight `:pts` proteins at 25 to 32, and `states` must
be strictly increasing in `species_index` so the composed state vector and the
registry agree on order.
"""
pts_transport_states() = vcat([:M_lac__L_e], species_in_group(:pts))

# The four foreign pools. Every one is both read and written: the cascade's
# reverse terms read pyruvate and glucose-6-phosphate, which the forward terms
# produce, and that is why each of those two carries an edge in both directions
# rather than only the outbound one its net production would suggest.
const PTS_FOREIGN = (:M_pep_c, :M_lac__L_c, :M_pyr_c, :M_g6p_c)

"""
    PtsTransport(; data_dir, asserted_gstd, free, volume_ratio, radius_nm)

The phosphotransferase cascade and lactate export.

- `asserted_gstd = 2.0` — the geometric standard deviation of every asserted
  prior. The upstream transport file has no `Parameter` table and no
  uncertainty column, so these eleven priors are **this project's, not the
  model's** (spec §4 D7). One number rather than eleven, so that the thing we
  are inventing is visible as one thing. What width is *defensible* is worth
  deciding against a sensitivity run rather than by argument, and matters only
  once one of the eleven is freed.
- `free` — names from [`PTS_SCALARS`](@ref) to expose as free parameters.
  Everything is `fixed = true` by default, so nothing is sampled from an
  invented prior without an explicit act. `:r_cell_nm` is always free whatever
  is passed; see below.
- `volume_ratio = 1e5` — the medium-to-cell volume ratio `R` of spec §4 D6.
- `radius_nm = COREA_INITIAL_RADIUS_NM` — the geometry the export law uses when
  no handshake is running.

**The radius is an inbound volume channel, not a constant.** `3P/r` collapses
to a single number at 200 nm, and folding it would be wrong, because phase 5
makes the radius grow and the export rate constant changes with it. So this
module declares an inbound [`VolumeEdge`](@ref) on `M_lac__L_e` — its own state
whose rate law reads the geometry — with `quantity = :radius_nm`, and the
driver fills `:r_cell_nm` at every handshake. Two consequences:

1. The slot is in **nanometres** while the permeability is in metres per
   second, so the rate law divides by `r_cell_nm * 1e-9`. The unit-bearing
   `quantity` symbol exists precisely so that this conversion is stated
   somewhere a reader can check.
2. A standalone homogeneous build executes no handshake, so the slot keeps the
   value declared here. It is the registry's initial radius, 200 nm, so a
   standalone run is the published rate law at the published geometry. Note the
   driver writes `radius_from_area_nm(502831) = 200.03505` nm over it rather
   than a round 200.0 — the two agree to the four significant figures the
   source states the radius in, and differ beyond them, deliberately.

`:r_cell_nm` is therefore always free: `param_slot` must name one of the
module's own free parameters, because only those reach the composed parameter
vector. Inference therefore samples it, and **only under the handshake driver
does anything overwrite it.** Under a driver its posterior is its prior and
must not be read as an identifiability result — the trap
[`rebuilt_params`](@ref) and [`VolumeEdge`](@ref) both document, whose fix
belongs to spec phases 15 to 17. In a *homogeneous* ODE composition nothing
overwrites it at all, and since only the ratio `P/r` enters the export law, a
run that frees `:p_lact2r` samples a flat ridge in the two. Pin `:r_cell_nm`
before inferring on a homogeneous build; the geometry is not an unknown there
either.
"""
struct PtsTransport{N} <: AbstractSubModel
    params::Vector{InferParameter}
    q::SVector{N_PTS_SCALARS, Float64}
    free_slots::SVector{N, Int}
end

function PtsTransport(; data_dir::AbstractString = COREA_DATA_DIR,
                        asserted_gstd::Real = 2.0,
                        free = Symbol[],
                        volume_ratio::Real = 1e5,
                        radius_nm::Real = COREA_INITIAL_RADIUS_NM)
    asserted_gstd > 1 || throw(ArgumentError(
        "asserted_gstd must exceed 1; got $asserted_gstd. It is a geometric " *
        "standard deviation, so 1 is a point mass and below 1 is not a width"))
    radius_nm > 0 || throw(ArgumentError(
        "radius_nm must be positive; got $radius_nm. The export law divides by " *
        "it, so a non-positive radius gives an infinite or undefined rate at " *
        "the first right-hand side call rather than a message here"))
    volume_ratio > 0 || throw(ArgumentError(
        "volume_ratio must be positive; got $volume_ratio. External lactate " *
        "accumulates at 1/R, so a non-positive ratio gives an infinite or " *
        "undefined derivative at the first right-hand side call"))
    unknown = setdiff(Symbol.(free), PTS_SCALARS)
    isempty(unknown) || throw(ArgumentError(
        "PtsTransport cannot free $(join(string.(":", unknown), ", ")): not one " *
        "of $(PTS_SCALARS). The nine initial conditions are held fixed and are " *
        "not freeable here"))

    rates = read_source_table(joinpath(data_dir, "pts_transport.tsv");
                              file = "transport")
    ics = read_source_table(joinpath(data_dir, "pts_initial_conditions.tsv");
                            file = "model_ics")
    tables = [rates, ics]

    freed = Set(Symbol.(free))
    push!(freed, :r_cell_nm)   # required: a param_slot must be a free parameter

    # `governing` is passed on every import below even though the two extracts'
    # identifier sets are disjoint. A declared governing file that is not the
    # holder is refused, so this turns a future overlap between them into a
    # load-time error rather than a silent choice.
    params = InferParameter[]
    for (name, id) in zip(PTS_SCALARS, PTS_RATE_IDS)
        push!(params, load_parameter(tables, id;
                                     name = name,
                                     module_id = :PtsTransport,
                                     prior = _asserted_prior(rates.values[id],
                                                             asserted_gstd),
                                     role = :rate,
                                     fixed = !(name in freed),
                                     governing = "transport",
                                     table = "Quantity"))
    end

    # Three scalars with no upstream row of their own. The radius and the
    # volume ratio are ours; the glucose clamp is the published model's value,
    # carried as a fixed parameter because the resolver refuses a chemostatted
    # input and prescribes exactly this shape.
    glucose_mM = held_value(:M_glc__D_e)
    push!(params, _asserted_scalar(:r_cell_nm, radius_nm, asserted_gstd,
                                   !(:r_cell_nm in freed)))
    push!(params, _asserted_scalar(:lac_volume_ratio, volume_ratio, asserted_gstd,
                                   !(:lac_volume_ratio in freed)))
    # `:asserted`, which is what the registry records for this species: the
    # medium concentration is a point value with no quantified uncertainty
    # anywhere upstream, so any prior on it would be ours. That is a different
    # claim from the *clamp* being ours, which it is not — the published model
    # clamps every external species, the edge carries `origin = :published`,
    # and no `:clamp` deviation is registered for it.
    push!(params, InferParameter(Float64(glucose_mM),
                                 Dirac(Float64(glucose_mM)), true,
                                 :glc_e_mM, :PtsTransport, :rate,
                                 ParameterSource("registry";
                                                 identifier = "M_glc__D_e",
                                                 informedness = :asserted)))

    for species in pts_transport_states()
        id = "conc_$(species)"
        value = ics.values[id]
        push!(params, load_parameter(tables, id;
                                     name = Symbol(species, "0"),
                                     module_id = :PtsTransport,
                                     prior = _asserted_prior(value, asserted_gstd),
                                     role = :initial_condition,
                                     fixed = true,
                                     governing = "model_ics",
                                     table = "Quantity"))
    end

    q = SVector{N_PTS_SCALARS, Float64}(ntuple(N_PTS_SCALARS) do i
        name = PTS_SCALARS[i]
        params[findfirst(p -> p.name === name, params)].value
    end)

    # p_local is `model_free_params(parameters(m))` in declaration order, so
    # slot j of that vector is the j-th free scalar here. Initial conditions are
    # always fixed, so they never enter the mapping.
    live = [findfirst(==(p.name), PTS_SCALARS) for p in model_free_params(params)]
    any(isnothing, live) && error(
        "Internal: a free PtsTransport parameter is not one of PTS_SCALARS, so " *
        "the parameter vector and the rate law would disagree on position")
    slots = SVector{length(live), Int}(Tuple(live))

    return PtsTransport{length(live)}(params, q, slots)
end

# A log-normal at the imported value, or — for external lactate, whose published
# initial condition is exactly 0.0 — a point mass at zero, since log(0) is not a
# median. Both are held fixed, so neither is ever drawn from.
_asserted_prior(value::Real, gstd::Real) =
    value > 0 ? LogNormal(log(value), log(gstd)) : Dirac(0.0)

_asserted_scalar(name::Symbol, value::Real, gstd::Real, fixed::Bool) =
    InferParameter(Float64(value), _asserted_prior(value, gstd), fixed,
                   name, :PtsTransport, :rate,
                   ParameterSource("asserted_by_this_project";
                                   identifier = String(name),
                                   informedness = :asserted))

states(::PtsTransport) = pts_transport_states()
parameters(m::PtsTransport) = m.params
inputs(::PtsTransport) = collect(PTS_FOREIGN)
contributed_states(::PtsTransport) = collect(PTS_FOREIGN)
membrane_protein_states(::PtsTransport) = [:M_ptsg_c, :M_ptsg_P_c]
extracellular_states(::PtsTransport) = [:M_lac__L_e]

"""
    coupling(m::PtsTransport)

Ten edges. Every `peer` is left unnamed, so the module resolves against
whichever composition it is placed in rather than naming a sibling it does not
depend on.

| Species | Kind | Direction | Why |
|---|---|---|---|
| `M_glc__D_e` | clamped | in | GLCpts4 draws it, held at the registry's 40 mM, `origin = :published` |
| `M_pep_c` | mass | in | GLCpts0 draws it |
| `M_lac__L_c` | mass | in | the exporter draws it |
| `M_pyr_c` | mass | in and out | GLCpts0 produces it forward and draws it back |
| `M_g6p_c` | mass | in and out | GLCpts4 produces it forward and draws it back |
| `M_ptsg_c`, `M_ptsg_P_c` | volume | out | ptsG's count sets the membrane surface area |
| `M_lac__L_e` | volume | in | the export law reads the cell radius |

**Pyruvate and glucose-6-phosphate carry an edge in each direction**, which the
archived design's five-edge table does not. Both are reversible steps whose
reverse term reads the pool the forward term fills, and `inputs` is the only
channel that wires a foreign state into `dynamics` — and `inputs` is held to the
*inbound* edges. One outbound edge alone would declare the production and leave
the rate law unable to see the pool it draws back from.

The glucose clamp is `origin = :published`: the published model clamps every
external species, so it is not a departure and must not appear in
[`reduction_declarations`](@ref). The transport reconstruction's 42.77 mM is
superseded by `setICs_two.py:279`'s 40, and declaring 42.77 would in any case
throw, since a clamp on a chemostat is validated against the registry.
"""
function coupling(::PtsTransport)
    edges = CouplingEdge[
        ClampedEdge(species = :M_glc__D_e, direction = :in,
                    held_value = held_value(:M_glc__D_e), origin = :published),
        MassEdge(species = :M_pep_c, direction = :in),
        MassEdge(species = :M_lac__L_c, direction = :in),
        MassEdge(species = :M_pyr_c, direction = :in),
        MassEdge(species = :M_pyr_c, direction = :out),
        MassEdge(species = :M_g6p_c, direction = :in),
        MassEdge(species = :M_g6p_c, direction = :out),
        VolumeEdge(species = :M_lac__L_e, direction = :in,
                   param_slot = :r_cell_nm, quantity = :radius_nm),
    ]
    for s in (:M_ptsg_c, :M_ptsg_P_c)
        push!(edges, VolumeEdge(species = s, direction = :out))
    end
    return edges
end

"""
    reduction_notes(m::PtsTransport)

The medium-to-cell volume ratio, which is ours.

The published model clamps every external species, so its export drains into a
pool pinned at zero and never saturates. The registry instead lists external
lactate among its dynamic states and the resolver rejects a clamp on a state a
module integrates, so integrating it in a single shared volume would make export
saturate as the pools equilibrate — the 691 mM would not leave, it would move
outside. The ratio reproduces the published unsaturated efflux through a
genuinely dynamic state.
"""
reduction_notes(m::PtsTransport) = [
    "external lactate accumulates at 1/R of the export rate with R declared " *
    "at $(_scalar(m, :lac_volume_ratio)), a medium-to-cell volume ratio asserted by " *
    "this project: it reproduces the published model's constant external pool, " *
    "which clamps M_lac__L_e and so never saturates the exporter, through the " *
    "dynamic state the registry requires. Chosen against a stated criterion — " *
    "external lactate must stay below one percent of steady cytosolic lactate " *
    "over a full cycle, which 1e5 meets at 0.47% and 1e4 misses at 4.7%",
]

"""
    lactate_export_rationale(m::PtsTransport) -> String

Why lactate export is not removable, recorded where the code is rather than only
in the planning artefacts.

An earlier Core A specification deleted this reaction. That was an error, and it
is recorded here so it is not repeated: glycolysis makes two lactate per
glucose, so at Core A′'s ~1,106 glucose/s a cycle produces ~13.9 M lactate —
about 691 mM at initial volume and about 345 mM even at doubled volume, forty
times the entire measured phosphate pool. Without export the core is not a
sustained pathway but one that poisons itself in minutes, and carbon balance is
open, so spec §3's check 2 cannot even be posed.
"""
lactate_export_rationale(::PtsTransport) =
    "Lactate export is required, not optional. Glycolysis makes two lactate " *
    "per glucose; at ~1,106 glucose/s a full cycle produces ~13.9 M lactate, " *
    "which is ~691 mM at initial volume and ~345 mM at doubled volume — about " *
    "forty times the measured phosphate pool. An earlier Core A specification " *
    "deleted this reaction in error; without it carbon balance is open and " *
    "spec §3 check 2 cannot be posed."

_scalar(m::PtsTransport, name::Symbol) =
    m.q[findfirst(==(name), PTS_SCALARS)]

# The live scalars: the struct's values with every free one replaced by its
# entry in the composed parameter vector. `Base.setindex` on an `SVector`
# returns a new `SVector`, so this allocates nothing and stays type-stable.
#
# **The element type comes from `p`, not from the struct.** Every
# gradient-based path carries `ForwardDiff.Dual` in `p` — `infer` defaults to
# `ForwardDiffSensitivity` and `check_identifiability` differentiates the
# forward map — and splicing a `Dual` into an `SVector{N, Float64}` would try
# to convert it to `Float64` and throw a `MethodError`. Promoting first is what
# keeps this module inside the AD contract `test/test_contributions.jl` asserts
# for every other composition.
@inline function _live_scalars(m::PtsTransport, p)
    T = promote_type(Float64, eltype(p))
    q = SVector{N_PTS_SCALARS, T}(m.q)
    @inbounds for j in eachindex(m.free_slots)
        q = Base.setindex(q, convert(T, p[j]), m.free_slots[j])
    end
    return q
end

# The six reaction rates, in one place so `dynamics` and `contributions` cannot
# disagree about them. `q` is the live scalar vector from `_live_scalars`, `u`
# this module's own slice in `states` order, and `u_inputs` the four foreign
# pools in `inputs` order.
@inline function _pts_rates(q, u, u_inputs)
    kf0, kr0, kf1, kr1, kf2, kr2, kf3, kr3, kf4, kr4, perm, r_nm, _, glc_e = q
    lac_e, ptsi, ptsi_P, ptsh, ptsh_P, crr, crr_P, ptsg, ptsg_P = u
    pep, lac_c, pyr, g6p = u_inputs

    v0 = kf0 * ptsi * pep - kr0 * ptsi_P * pyr
    v1 = kf1 * ptsh * ptsi_P - kr1 * ptsh_P * ptsi
    v2 = kf2 * ptsh_P * crr - kr2 * ptsh * crr_P
    v3 = kf3 * ptsg * crr_P - kr3 * crr * ptsg_P
    v4 = kf4 * ptsg_P * glc_e - kr4 * ptsg * g6p
    # 3/r is a sphere's surface-to-volume ratio, so the radius must reach this
    # in metres to match the permeability's m/s; the slot carries nanometres.
    v_export = perm * (lac_c - lac_e) * 3 / (r_nm * 1e-9)

    return SVector(v0, v1, v2, v3, v4, v_export)
end

"""
    dynamics(u, p, t, m::PtsTransport, u_inputs)

Each cascade step transfers a phosphate between adjacent carriers, so each
carrier's two forms appear with opposite signs and every one of the four sums is
invariant independently — which is what spec §3's check 5 asserts, as four
bounds rather than one.
"""
function dynamics(u, p, t, m::PtsTransport, u_inputs)
    q = _live_scalars(m, p)
    v0, v1, v2, v3, v4, v_export = _pts_rates(q, u, u_inputs)
    ratio = q[13]
    return SVector(
        v_export / ratio,   # M_lac__L_e, referred to the medium's volume
        -v0 + v1,           # M_ptsi_c
        v0 - v1,            # M_ptsi_P_c
        -v1 + v2,           # M_ptsh_c
        v1 - v2,            # M_ptsh_P_c
        -v2 + v3,           # M_crr_c
        v2 - v3,            # M_crr_P_c
        -v3 + v4,           # M_ptsg_c
        v3 - v4,            # M_ptsg_P_c
    )
end

"""
    contributions(u, p, t, m::PtsTransport, u_inputs)

The signed terms this module adds to the four pools it does not own, in
[`contributed_states`](@ref) order. Pyruvate and glucose-6-phosphate each get
one net term rather than two, because a reversible step's forward and reverse
directions are one flux with a sign.
"""
function contributions(u, p, t, m::PtsTransport, u_inputs)
    v0, _, _, _, v4, v_export = _pts_rates(_live_scalars(m, p), u, u_inputs)
    return SVector(-v0,         # M_pep_c, consumed by GLCpts0
                   -v_export,   # M_lac__L_c, leaving the cell
                   v0,          # M_pyr_c, produced by GLCpts0
                   v4)          # M_g6p_c, produced by GLCpts4
end

export PtsTransport, PTS_SCALARS, PTS_RATE_IDS, COREA_DATA_DIR,
       pts_transport_states, lactate_export_rationale
