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

One parsed source file: its logical name (what provenance records), the
identifier-to-value map it contributed, and how well each value is known.
"""
struct SourceTable
    file::String
    values::Dict{String, Float64}
    informedness::Dict{String, Symbol}
end

"""
    read_source_table(path; file=<extension-free basename>, id_column="ID",
                      value_column="Mode", gstd_column="GeometricStd")

Parse an SBtab-shaped tab-separated table.

`file` is the logical name provenance records — what `governing` declarations
and the registry's `source_file` entries match against, by string equality. It
defaults to the path's extension-free basename, so reading
`.../central_balanced.tsv` records `"central_balanced"`, matching the
registry's logical source names.

Lines beginning `!!` (the SBtab document declaration), `%` (comments), and blank
lines are skipped. The header row's column names may carry SBtab's leading `!`,
which is stripped.

`value_column` defaults to `Mode` because that is the column the published
simulator reads. Reading a balancing-distribution column instead is the error
that produced a fictitious capacity bottleneck in an earlier version of the
scoping note.

How well each value is known is read from an `Informedness` column when the
table has one, and otherwise derived: a row whose geometric standard deviation
sits at or above the prior width ([`PRIOR_DEFAULT_GSTD`](@ref), 10) is
`:prior_default`, a row with no uncertainty column at all is `:asserted` — any
prior on it is this project's, not the model's — and anything else is
`:balanced`.
"""
function read_source_table(path::AbstractString;
                           file::AbstractString = first(splitext(basename(path))),
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
        # A truncated row is corruption, not something to skip: silently
        # dropping it could collapse a two-file ambiguity to one holder and
        # bypass the governs machinery entirely — or, cut short before the
        # uncertainty columns, silently re-derive an informedness the row
        # declared. `split` keeps empty cells, so a short row means missing
        # tabs, not a legitimately empty final column.
        required = max(id_idx, value_idx,
                       something(gstd_idx, 0), something(informedness_idx, 0))
        length(cells) < required && throw(ArgumentError(
            "Source table $path has a truncated row (\"$(first(line, 60))\"): " *
            "$(length(cells)) cell(s) where the widest header column is " *
            "number $required"))

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
        informedness[id] = _row_informedness(cells, gstd_idx, informedness_idx,
                                             id, path)
    end

    return SourceTable(String(file), values, informedness)
end

function _column_index(header::Vector{<:AbstractString}, name::AbstractString,
                       path::AbstractString)
    idx = findfirst(==(name), header)
    idx === nothing && throw(ArgumentError(
        "Source table $path has no column \"$name\". Its columns are " *
        "$(join(header, ", "))"))
    return idx
end

function _row_informedness(cells, gstd_idx, informedness_idx, id, path)
    if informedness_idx !== nothing
        raw = strip(cells[informedness_idx])
        if !isempty(raw)
            declared = Symbol(raw)
            # A declared informedness must be a real one: silently re-deriving
            # over a misspelling would let the declaration rot unnoticed.
            declared in INFORMEDNESS || throw(ArgumentError(
                "Row :$id of $path declares informedness \"$raw\", which is not " *
                "one of $INFORMEDNESS"))
            return declared
        end
    end

    # No uncertainty column, or an empty cell in it: the source carries a point
    # value, so any prior on it is ours rather than the model's.
    gstd_idx === nothing && return :asserted
    raw_gstd = strip(cells[gstd_idx])
    isempty(raw_gstd) && return :asserted

    # A non-empty cell must be a real width. Falling back to :asserted over
    # "N/A" would misfile a corrupt row, and NaN compares false against the
    # prior width, which would report an unquantified value as informed.
    gstd = tryparse(Float64, raw_gstd)
    (gstd === nothing || isnan(gstd)) && throw(ArgumentError(
        "Row :$id of $path has an unparseable geometric standard deviation " *
        "\"$raw_gstd\". An empty cell means the source asserts a point value; " *
        "anything else must be a number, or the row's informedness would be " *
        "silently misclassified"))
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
        push!(report, AmbiguousValue(id, sorted, allequal(last(p) for p in sorted)))
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

A `governing` declaration is also checked when only one table holds the
identifier: a declared governing file that is not the holder means the table
set and the declaration disagree — a stale or truncated table — and failing is
better than silently importing from the wrong file with clean provenance.

An identifier of the form `conc_<species>` naming a registry species is checked
against the registry: the imported value must match the registry row's
`initial_value` where one is recorded, so the registry and the loader cannot
carry two versions of the same concentration.
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
        # The declared governing file is still binding: a single holder that is
        # not the declared governor means the governing table is stale or
        # missing, and importing from the other file would be exactly the
        # silent cross-file substitution this loader exists to prevent.
        if governing !== nothing && only(holders).file != governing
            throw(ArgumentError(
                "Identifier $id declares governing file \"$governing\" but is " *
                "held only by $(only(holders).file). Either the governing table " *
                "is stale or the declaration is wrong; refusing to import " *
                "$(only(holders).values[id]) from a file the declaration rejects"))
        end
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

    _check_registry_agreement(id, chosen)

    alternatives = [t.file => t.values[id] for t in holders if t.file != chosen.file]

    provenance = ParameterSource(chosen.file;
                                 table = table,
                                 identifier = id,
                                 informedness = get(chosen.informedness, id, :not_imported),
                                 alternatives = alternatives)

    return InferParameter(chosen.values[id], prior, fixed, name, module_id,
                          role, provenance)
end

# The registry transcribes initial concentrations by hand; the loader imports
# the same tables live. Two channels for one number is the cross-file trap one
# layer up, so where an identifier names a registry species with a recorded
# value, the two must agree.
function _check_registry_agreement(id::String, chosen::SourceTable)
    startswith(id, "conc_") || return nothing
    species = Symbol(chopprefix(id, "conc_"))
    is_registered(species) || return nothing

    entry = species_entry(species)
    entry.initial_value === nothing && return nothing

    value = chosen.values[id]
    # The registry transcribes values to four decimal places while the source
    # tables carry full precision, so agreement means equal up to that
    # transcription — half a unit in the fourth decimal — not bit-for-bit.
    isapprox(value, entry.initial_value; atol = 5e-5) || throw(ArgumentError(
        "Identifier $id imports $value from $(chosen.file), but the registry " *
        "records :$species at $(entry.initial_value)" *
        (entry.source_file === nothing ? "" : " from $(entry.source_file)") *
        ". The registry row and the imported value are two copies of one " *
        "number; import from the governing file or correct the registry"))
    return nothing
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

export SourceTable, AmbiguousValue, read_source_table, ambiguity_report,
       disagreements, load_parameter, governing_choices
