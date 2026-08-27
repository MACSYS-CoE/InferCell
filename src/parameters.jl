"""
How well a value is known, shared by [`ParameterSource`](@ref) and the Core A′
species registry so that a kinetic constant and an initial condition describe
their provenance the same way.

- `:balanced` — read from a balancing distribution's `Mode` column.
- `:prior_default` — the prior median at prior width. Nothing informed it.
- `:asserted` — the source carries a point value with no quantified
  uncertainty, so any prior on it is this project's rather than the model's.
- `:not_imported` — no value has been taken from the source model yet.
"""
const INFORMEDNESS = (:balanced, :prior_default, :asserted, :not_imported)

# Shared vocabulary-membership check. `context` names the thing being
# validated (an edge kind, a species, a struct) in the error message.
function _check_vocab(context, field::Symbol, value, allowed)
    value in allowed || throw(ArgumentError(
        "$context field `$field` is :$value, which is not valid. " *
        "Must be one of $allowed"))
    return value
end

# The prior width that marks a value as uninformed: a geometric standard
# deviation at or above this is the balancing prior's default, and a row
# sitting there was informed by nothing. Shared by the species registry and
# the loader so the rule lives once.
const PRIOR_DEFAULT_GSTD = 10.0

"""
    ParameterSource(file; table, identifier, informedness, alternatives)

Where an imported value came from.

Carried per value rather than resolved once by hand, because the same
cross-file ambiguity has already produced two errors of record: PGK3 and PYK3
are parameterised in both the central and the nucleotide balanced file with
different modes, and GTP's initial concentration reads 1.6627 mM in the
nucleotide file against the 0.1 mM prior default in the central one.

`alternatives` lists the same identifier's value in every *other* file it
appears in, so a value that was chosen in the presence of an alternative says so
in the assembled model rather than only in the module's source.
"""
struct ParameterSource
    file::String
    table::Union{String, Nothing}
    identifier::Union{String, Nothing}
    informedness::Symbol
    alternatives::Vector{Pair{String, Float64}}
end

function ParameterSource(file::AbstractString;
                         table = nothing,
                         identifier = nothing,
                         informedness::Symbol = :balanced,
                         alternatives = Pair{String, Float64}[])
    _check_vocab(:ParameterSource, :informedness, informedness, INFORMEDNESS)
    return ParameterSource(String(file),
                           table === nothing ? nothing : String(table),
                           identifier === nothing ? nothing : String(identifier),
                           informedness,
                           Pair{String, Float64}[String(k) => Float64(v)
                                                 for (k, v) in alternatives])
end

"""
    was_chosen_over_alternative(s::ParameterSource) -> Bool

Whether this value's source file was chosen in the presence of another file
carrying the same identifier.
"""
was_chosen_over_alternative(s::ParameterSource) = !isempty(s.alternatives)

"""
    InferParameter(value, prior, fixed, name, module_id, role[, provenance])

A single parameter declared by a sub-model: its nominal value, prior, identity,
and role. `role` is one of `:rate`, `:initial_condition`, `:observation`; the
inference layer uses `role` and `fixed` to decide which parameters are sampled.

`provenance` is an optional [`ParameterSource`](@ref) recording the file the
value was imported from. It defaults to `nothing`, so a parameter constructed
directly — as every sub-model predating the Core A′ contract does — reports its
provenance as absent rather than as an incorrect source.
"""
struct InferParameter
    value::Float64
    prior::Distribution
    fixed::Bool
    name::Symbol
    module_id::Symbol
    role::Symbol  # :rate, :initial_condition, :observation
    provenance::Union{ParameterSource, Nothing}
end

InferParameter(value, prior, fixed, name, module_id, role) =
    InferParameter(value, prior, fixed, name, module_id, role, nothing)

"""
    free_params(params)

Return the subset of `params` that are not fixed (i.e. inferred).
"""
free_params(params::Vector{InferParameter}) = filter(p -> !p.fixed, params)

"""
    rate_params(params)

Return parameters whose `role` is `:rate` (kinetic rate constants).
"""
rate_params(params::Vector{InferParameter}) = filter(p -> p.role == :rate, params)

"""
    ic_params(params)

Return parameters whose `role` is `:initial_condition` (initial state values).
"""
ic_params(params::Vector{InferParameter}) = filter(p -> p.role == :initial_condition, params)

"""
    obs_params(params)

Return parameters whose `role` is `:observation` (e.g. measurement noise σ).
"""
obs_params(params::Vector{InferParameter}) = filter(p -> p.role == :observation, params)

"""
    model_free_params(params)

Free parameters of the dynamical model — everything except observation
parameters. These are the entries that appear in the ODE/SSA parameter vector.
"""
model_free_params(params::Vector{InferParameter}) =
    filter(p -> p.role != :observation, free_params(params))

"""
    obs_free_params(params)

Free observation parameters (e.g. measurement noise σ). These appear in the
likelihood but not in the dynamics RHS.
"""
obs_free_params(params::Vector{InferParameter}) =
    filter(p -> p.role == :observation, free_params(params))

"""
    ode_free_params(params)

Backward-compatibility alias for [`model_free_params`](@ref).
"""
const ode_free_params = model_free_params

"""
    provenance_of(p::InferParameter) -> Union{ParameterSource, Nothing}

The source record of a parameter, or `nothing` where it was constructed
directly rather than imported.
"""
provenance_of(p::InferParameter) = p.provenance

"""
    source_file(p::InferParameter) -> Union{String, Nothing}

The file a parameter's value was imported from, or `nothing`.
"""
source_file(p::InferParameter) =
    p.provenance === nothing ? nothing : p.provenance.file

"""
    informedness(p::InferParameter) -> Symbol

How well a parameter's value is known — one of [`INFORMEDNESS`](@ref).
A parameter with no provenance reports `:not_imported`.
"""
informedness(p::InferParameter) =
    p.provenance === nothing ? :not_imported : p.provenance.informedness

"""
    asserted_prior_params(params)

Parameters whose prior this project asserted rather than inherited from the
source model — values that carry no quantified uncertainty where they came
from, such as the ten PTS mass-action constants.

The balanced glycolytic parameters do not appear here: their priors are the
model's own.
"""
asserted_prior_params(params::Vector{InferParameter}) =
    filter(p -> informedness(p) === :asserted, params)

"""
    uninformed_params(params)

Parameters sitting at the prior median at prior width. Nothing informed them,
and they should not be read as measurements.
"""
uninformed_params(params::Vector{InferParameter}) =
    filter(p -> informedness(p) === :prior_default, params)

"""
    provenance_conflicts(params) -> Vector{Tuple{Symbol, Vector{String}}}

Parameter names that arrive from more than one source file, with the files
involved.

This is the cross-file trap in its most easily-missed form: the values agree, so
no consistency check fires, and the first-seen source silently wins
deduplication. Reporting it is the point — see the two errors of record in
[`ParameterSource`](@ref).
"""
function provenance_conflicts(params::Vector{InferParameter})
    files = Dict{Symbol, Vector{String}}()
    for p in params
        f = source_file(p)
        f === nothing && continue
        bucket = get!(files, p.name, String[])
        f in bucket || push!(bucket, f)
    end
    return sort!([(name, sort(fs)) for (name, fs) in files if length(fs) > 1];
                 by = first)
end

"""
    unique_params(params)

Deduplicate `params` by `name`, preserving first-occurrence order. Used when
merging parameter lists from several sub-models that share named parameters.
"""
function unique_params(params::Vector{InferParameter})
    seen = Set{Symbol}()
    result = InferParameter[]
    for p in params
        if !(p.name in seen)
            push!(seen, p.name)
            push!(result, p)
        end
    end
    return result
end
