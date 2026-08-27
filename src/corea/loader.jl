"""
Provenance-carrying parameter loading.

Reads a balanced table from a caller-supplied path and returns
[`InferParameter`](@ref)s that know which file they came from.

The tables are not vendored in this repository — wave 0 delivers the loader and
wave 1 imports through it — so the loader is path-agnostic: a caller supplies
the path and the logical name to record as provenance.

Parsing is deliberately Base-only. The format is SBtab-shaped tab-separated
text, which needs no quoting logic, and a parsing dependency would have to be
fetched over a network the compute nodes do not have.
"""

"""
    SourceTable

One parsed source file: its logical name (what provenance records), its path on
disk, the identifier-to-value map it contributed, and how well each value is
known.
"""
struct SourceTable
    file::String
    path::String
    values::Dict{String, Float64}
    informedness::Dict{String, Symbol}
end

# The prior median and prior width that mark a value as uninformed: 0.1 mM at a
# geometric standard deviation of 10 is the balancing prior's default, and a row
# sitting there was informed by nothing.
const PRIOR_DEFAULT_GSTD = 10.0

"""
    read_source_table(path; file=basename(path), id_column="ID",
                      value_column="Mode", gstd_column="GeometricStd")

Parse an SBtab-shaped tab-separated table.

Lines beginning `!!` (the SBtab document declaration), `%` (comments), and blank
lines are skipped. The header row's column names may carry SBtab's leading `!`,
which is stripped.

`value_column` defaults to `Mode` because that is the column the published
simulator reads. Reading a balancing-distribution column instead is the error
that produced a fictitious capacity bottleneck in an earlier version of the
scoping note.

How well each value is known is read from an `Informedness` column when the
table has one, and otherwise derived: a row whose geometric standard deviation
sits at the prior width is `:prior_default`, a row with no uncertainty column at
all is `:asserted` — any prior on it is this project's, not the model's — and
anything else is `:balanced`.
"""
function read_source_table(path::AbstractString;
                           file::AbstractString = basename(path),
                           id_column::AbstractString = "ID",
                           value_column::AbstractString = "Mode",
                           gstd_column::AbstractString = "GeometricStd")
    isfile(path) || throw(ArgumentError(
        "Source table not found at $path. The balanced tables are not vendored " *
        "in this repository; supply the path to a checkout of the source model"))

    lines = filter(readlines(path)) do line
        s = strip(line)
        !isempty(s) && !startswith(s, "!!") && !startswith(s, "%")
    end
    isempty(lines) && throw(ArgumentError("Source table $path has no rows"))

    header = [lstrip(strip(c), '!') for c in split(lines[1], '\t')]
    id_idx = _column_index(header, id_column, path)
    value_idx = _column_index(header, value_column, path)
    gstd_idx = findfirst(==(gstd_column), header)
    informedness_idx = findfirst(==("Informedness"), header)

    values = Dict{String, Float64}()
    informedness = Dict{String, Symbol}()

    for line in lines[2:end]
        cells = split(line, '\t')
        length(cells) < max(id_idx, value_idx) && continue

        id = String(strip(cells[id_idx]))
        isempty(id) && continue

        raw = strip(cells[value_idx])
        value = tryparse(Float64, raw)
        value === nothing && throw(ArgumentError(
            "Row :$id of $path has an unparseable $value_column value \"$raw\""))

        haskey(values, id) && throw(ArgumentError(
            "Identifier $id appears twice in $path. A source table must name " *
            "each value once"))

        values[id] = value
        informedness[id] = _row_informedness(cells, gstd_idx, informedness_idx)
    end

    return SourceTable(String(file), String(path), values, informedness)
end

function _column_index(header::Vector{<:AbstractString}, name::AbstractString,
                       path::AbstractString)
    idx = findfirst(==(name), header)
    idx === nothing && throw(ArgumentError(
        "Source table $path has no column \"$name\". Its columns are " *
        "$(join(header, ", "))"))
    return idx
end

function _row_informedness(cells, gstd_idx, informedness_idx)
    if informedness_idx !== nothing && length(cells) >= informedness_idx
        declared = Symbol(strip(cells[informedness_idx]))
        declared in INFORMEDNESS && return declared
    end

    # No uncertainty column: the source carries a point value, so any prior on
    # it is ours rather than the model's.
    gstd_idx === nothing && return :asserted
    length(cells) < gstd_idx && return :asserted

    gstd = tryparse(Float64, strip(cells[gstd_idx]))
    gstd === nothing && return :asserted
    return gstd >= PRIOR_DEFAULT_GSTD ? :prior_default : :balanced
end

"""
    AmbiguousValue

An identifier that appears in more than one source file, with each file's value
and whether the values agree.

A disagreement is the PGK3/PYK3 case — 319.5 in the central file against 140.8
in the nucleotide one. An agreement is still worth reporting, because which file
governs is then decided by composition order rather than by anyone's intent.
"""
struct AmbiguousValue
    identifier::String
    values::Vector{Pair{String, Float64}}
    agrees::Bool
end

"""
    ambiguity_report(tables) -> Vector{AmbiguousValue}

Every identifier present in more than one of `tables`, with each file's value.

Obtainable without running a simulation, and empty when every value comes from
exactly one file — an emptiness that is a positive result rather than a failure
to run.
"""
function ambiguity_report(tables::Vector{SourceTable})
    seen = Dict{String, Vector{Pair{String, Float64}}}()
    for t in tables
        for (id, v) in t.values
            push!(get!(seen, id, Pair{String, Float64}[]), t.file => v)
        end
    end

    report = AmbiguousValue[]
    for (id, entries) in seen
        length(entries) > 1 || continue
        sorted = sort(entries; by = first)
        agrees = all(p -> last(p) == last(sorted[1]), sorted)
        push!(report, AmbiguousValue(id, sorted, agrees))
    end
    return sort!(report; by = a -> a.identifier)
end

"""
    disagreements(tables) -> Vector{AmbiguousValue}

The subset of [`ambiguity_report`](@ref) whose values actually differ across
files.
"""
disagreements(tables::Vector{SourceTable}) =
    filter(a -> !a.agrees, ambiguity_report(tables))

"""
    load_parameter(tables, identifier; name, module_id, prior, role=:rate,
                   fixed=false, governing=nothing, table=nothing)

Import one value as a provenance-carrying [`InferParameter`](@ref).

When `identifier` appears in more than one of `tables`, `governing` must name
the file that wins. Without it the load fails, naming the identifier, both files
and both values — an unresolved ambiguity is an error at load time, never a
silent choice. With it, the parameter's provenance records both the file chosen
and the ones rejected.
"""
function load_parameter(tables::Vector{SourceTable}, identifier::AbstractString;
                        name::Symbol,
                        module_id::Symbol,
                        prior::Distribution,
                        role::Symbol = :rate,
                        fixed::Bool = false,
                        governing::Union{AbstractString, Nothing} = nothing,
                        table::Union{AbstractString, Nothing} = nothing)
    id = String(identifier)
    holders = [t for t in tables if haskey(t.values, id)]

    isempty(holders) && throw(ArgumentError(
        "Identifier $id is not in any of the supplied tables " *
        "($(join([t.file for t in tables], ", ")))"))

    chosen = if length(holders) == 1
        only(holders)
    else
        governing === nothing && throw(ArgumentError(
            "Identifier $id is present in more than one source file and no " *
            "governing file was declared: " *
            join(["$(t.file) = $(t.values[id])" for t in holders], ", ") *
            ". Declare which file governs; resolving this by hand once is how " *
            "it bites a third time"))

        match = findfirst(t -> t.file == governing, holders)
        match === nothing && throw(ArgumentError(
            "Governing file \"$governing\" does not hold identifier $id. " *
            "It is held by $(join([t.file for t in holders], ", "))"))
        holders[match]
    end

    alternatives = [t.file => t.values[id] for t in holders if t.file != chosen.file]

    provenance = ParameterSource(chosen.file;
                                 table = table,
                                 identifier = id,
                                 informedness = get(chosen.informedness, id, :not_imported),
                                 alternatives = alternatives)

    return InferParameter(chosen.values[id], prior, fixed, name, module_id,
                          role, provenance)
end

"""
    governing_choices(params) -> Vector{Tuple{Symbol, String, Vector{String}}}

Every parameter whose source file was chosen in the presence of an alternative,
as `(name, chosen file, rejected files)`.

This is what makes a cross-file decision visible in the assembled model rather
than only in the module's source.
"""
function governing_choices(params::Vector{InferParameter})
    out = Tuple{Symbol, String, Vector{String}}[]
    for p in params
        s = p.provenance
        s === nothing && continue
        was_chosen_over_alternative(s) || continue
        push!(out, (p.name, s.file, sort([first(a) for a in s.alternatives])))
    end
    return sort!(out; by = first)
end
