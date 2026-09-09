# Derive the Core A′ central-glycolysis parameter extract for spec/spec.md §11
# task 6.1.
#
# The upstream balanced table keys every row on the triple (QuantityType,
# Reaction id, Compound id) and has no single-identifier column at all, which is
# what `read_source_table` assumes. So the table is reshaped before the loader
# sees it, rather than vendored verbatim: `src/loader.jl` is framework and is not
# amended on a module branch (spec §10 R15), and a 5,000-row upstream file would
# not fit the loader in any case (spec §5).
#
# The reshape renames and selects. It changes no value: every `Mode` and
# `GeometricStd` cell is copied as the *string* the upstream file holds, never
# parsed and reprinted, so `product catalytic rate constant R_PGI` stays `650`
# rather than becoming `650.0`. `!UpstreamRow` carries the triple the row came
# from, so a renamed identifier stays traceable to its source line and the
# reshape is auditable rather than a retyping.
#
# Two columns are deliberately not read. `!UnconstrainedGeometricMean` is the
# unconstrained estimate rather than the balanced one — spec §4 D1 records
# reading it as an error of record — and the `Quantity` table earlier in the
# same document is the hand-patch layer the published simulator actually runs.
# Eight of the 32 Michaelis constants differ between the two, and taking the
# balanced column is a declared departure (spec §4 D3), registered by the module
# rather than hidden here.
#
# Usage:
#   julia dev/scripts/extract_central_glycolysis.jl <Minimal_Cell checkout> \
#         [output path]
#
# Base only, no project: the script must run on a login node with no depot.
# Re-running it over an unchanged checkout reproduces the file byte for byte,
# which is the round-trip requirement of spec §5.

const UPSTREAM_REL = joinpath("CME_ODE", "model_data",
                              "Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv")
const UPSTREAM_COMMIT = "db048ac"

# The ten reactions from glucose-6-phosphate through lactate, in pathway order.
# NOX is not among them: LDH_L already regenerates NAD⁺, it carries no
# gene-protein-reaction rule, and its catalytic constant was prior-dominated.
# The module registers that as a reduction declaration.
const REACTIONS = ["R_PGI", "R_PFK", "R_FBA", "R_TPI", "R_GAPD",
                   "R_PGK", "R_PGM", "R_ENO", "R_PYK", "R_LDH_L"]

# The thirteen owned states, in `COREA_SPECIES` order: eleven glycolytic
# intermediates then the redox pair. The `conc_` prefix is what fires the
# loader's registry-agreement check, so these thirteen check themselves against
# `src/organisms/coreA/registry.jl` on every load.
const SPECIES = ["M_g6p_c", "M_f6p_c", "M_fdp_c", "M_dhap_c", "M_g3p_c",
                 "M_13dpg_c", "M_3pg_c", "M_2pg_c", "M_pep_c", "M_pyr_c",
                 "M_lac__L_c", "M_nad_c", "M_nadh_c"]

const N_EXPECTED = (kcatF = 10, kcatR = 10, km = 32, conc = 13)

"""
    parameter_table(path) -> (rows, header)

The rows of the `Parameter` table alone. The upstream file is one SBtab document
holding six tables; splitting on `!!SBtab` and selecting by `TableName` is what
keeps the `Quantity` table — same document, same file, different meaning — from
being read by mistake.
"""
function parameter_table(path::AbstractString)
    lines = readlines(path)
    starts = findall(l -> startswith(l, "!!SBtab"), lines)
    isempty(starts) && error("$path holds no SBtab table declaration")
    i = findfirst(k -> occursin("TableName='Parameter'", lines[k]), starts)
    i === nothing && error("$path holds no table named 'Parameter'")
    first_row = starts[i] + 2                      # the declaration, then the header
    next = findfirst(>(starts[i]), starts[(i + 1):end])
    last_row = next === nothing ? length(lines) : starts[i + next] - 1
    header = split(rstrip(lines[starts[i] + 1], '\t'), '\t')
    rows = [split(l, '\t') for l in lines[first_row:last_row] if !isempty(strip(l, ['\t', ' ']))]
    return rows, header
end

"""
    column(header, name) -> Int

Position of a named column, by its header text rather than by a hard-coded
number. `Mode` is column four and the geometric standard deviation column seven,
but a file that moved them would otherwise be read silently at the wrong offset.
"""
function column(header::Vector{<:AbstractString}, name::AbstractString)
    i = findfirst(==(name), header)
    i === nothing && error("the Parameter table has no column $name; it has $(join(header, ", "))")
    return i
end

function main()
    length(ARGS) in (1, 2) ||
        error("usage: julia dev/scripts/extract_central_glycolysis.jl <Minimal_Cell checkout> [output path]")
    checkout = ARGS[1]
    out = length(ARGS) == 2 ? ARGS[2] :
          joinpath(@__DIR__, "..", "..", "src", "organisms", "coreA", "data",
                   "central_glycolysis.tsv")

    path = joinpath(checkout, UPSTREAM_REL)
    isfile(path) ||
        error("$path does not exist. Pass the root of a Luthey-Schulten-Lab/Minimal_Cell checkout at commit $UPSTREAM_COMMIT.")

    rows, header = parameter_table(path)
    c_type = column(header, "!QuantityType")
    c_rxn = column(header, "!Reaction:SBML:reaction:id")
    c_cpd = column(header, "!Compound:SBML:species:id")
    c_mode = column(header, "!Mode")
    c_gstd = column(header, "!UnconstrainedGeometricStd")

    # (quantity type, reaction, compound) -> (mode, gstd), as strings.
    cells = Dict{Tuple{String, String, String}, Tuple{String, String}}()
    # Michaelis rows, per reaction, in upstream order: the order is the source's
    # substrates-then-products order, which is the order the rate law's terms
    # are built in, so preserving it keeps the extract readable beside the model.
    km_of = Dict(r => Tuple{String, String, String}[] for r in REACTIONS)
    for row in rows
        length(row) < c_gstd && continue
        qt, rxn, cpd = row[c_type], row[c_rxn], row[c_cpd]
        mode, gstd = row[c_mode], row[c_gstd]
        key = (qt, rxn, cpd)
        haskey(cells, key) && error("the Parameter table holds $key twice")
        cells[key] = (mode, gstd)
        qt == "Michaelis constant" && rxn in REACTIONS && push!(km_of[rxn], (cpd, mode, gstd))
    end

    take(qt, rxn, cpd) = get(cells, (qt, rxn, cpd)) do
        error("the Parameter table has no $qt row for ($rxn, $cpd)")
    end

    out_rows = Tuple{String, String, String, String}[]
    emit(id, mode, gstd, upstream) = push!(out_rows, (id, mode, gstd, upstream))

    for r in REACTIONS
        mode, gstd = take("substrate catalytic rate constant", r, "")
        emit("kcatF_$r", mode, gstd, "substrate catalytic rate constant|$r|")
    end
    for r in REACTIONS
        mode, gstd = take("product catalytic rate constant", r, "")
        emit("kcatR_$r", mode, gstd, "product catalytic rate constant|$r|")
    end
    for r in REACTIONS, (cpd, mode, gstd) in km_of[r]
        emit("km_$(r)_$cpd", mode, gstd, "Michaelis constant|$r|$cpd")
    end
    for s in SPECIES
        mode, gstd = take("concentration", "", s)
        emit("conc_$s", mode, gstd, "concentration||$s")
    end

    counts = (kcatF = count(r -> startswith(r[1], "kcatF_"), out_rows),
              kcatR = count(r -> startswith(r[1], "kcatR_"), out_rows),
              km = count(r -> startswith(r[1], "km_"), out_rows),
              conc = count(r -> startswith(r[1], "conc_"), out_rows))
    counts == N_EXPECTED ||
        error("row classes are $counts, expected $N_EXPECTED — the upstream table changed shape")

    mkpath(dirname(out))
    open(out, "w") do io
        println(io, "!!SBtab TableType='Quantity' TableName='central balanced (Core A′ glycolysis extract)'")
        println(io, "!ID\t!Mode\t!GeometricStd\t!UpstreamRow")
        println(io, "% Generated by dev/scripts/extract_central_glycolysis.jl. Do not edit by hand.")
        println(io, "% Derived from $UPSTREAM_REL at Luthey-Schulten-Lab/Minimal_Cell $UPSTREAM_COMMIT,")
        println(io, "% table TableName='Parameter'. Identifiers are renamed; no value is changed.")
        println(io, "% See src/organisms/coreA/data/README.md for the regeneration command.")
        for (id, mode, gstd, upstream) in out_rows
            println(io, "$id\t$mode\t$gstd\t$upstream")
        end
    end
    println("wrote $(length(out_rows)) rows to $out: ", counts)
end

main()
