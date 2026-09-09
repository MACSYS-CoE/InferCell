#!/usr/bin/env julia
#
# Reshape the nucleotide-recycling rows of two published balanced SBtab files
# into the single-identifier form `src/loader.jl` reads.
#
#     julia dev/scripts/extract_nucleotide_recycling.jl <Minimal_Cell checkout>
#
# Writes two extracts, and the second one is the point. Every value this module
# imports exists in *both* balanced files, so vendoring only the governing file
# would leave every `load_parameter` call with a single holder and the
# `governing` declarations would never be exercised. With both loaded each of
# the 35 identifiers has two holders, an undeclared import fails, and a declared
# one records what it rejected (spec §4 D2).
#
# The reshape renames and reorders. It changes no value: `Mode` and
# `UnconstrainedGeometricStd` are copied as the strings the upstream file holds,
# so a float round-trip cannot lose a digit, and the `UpstreamRow` column keeps
# each renamed identifier traceable to the row it came from.
#
# No dependencies beyond Base, deliberately: this runs on a login node, where
# loading the project's package stack would invalidate the compute nodes' caches.

const UPSTREAM = Dict(
    "nucleotide_balanced" => "CME_ODE/model_data/Nucleotide_Kinetic_Parameters.tsv",
    "central_balanced" =>
        "CME_ODE/model_data/Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv",
)

const OUTPUTS = Dict(
    "nucleotide_balanced" => "src/organisms/coreA/data/nucleotide_recycling.tsv",
    "central_balanced" => "src/organisms/coreA/data/nucleotide_recycling_central.tsv",
)

# Canonical order: the two GTP-regenerating reactions, then the three that close
# a moiety, which is the order `dev/notes/reduced-syn3a-scoping.md` introduces
# them in.
const REACTIONS = ["R_PGK3", "R_PYK3", "R_ADK1", "R_GK1", "R_PPA"]

# Michaelis constants, per reaction, in the upstream file's own row order —
# substrates then products, which is the order the rate law's numerator reads
# them in. Listed rather than discovered so a silently dropped row is a length
# mismatch and not a shorter extract.
const MICHAELIS = [
    "R_PGK3" => ["M_13dpg_c", "M_gdp_c", "M_3pg_c", "M_gtp_c"],
    "R_PYK3" => ["M_gdp_c", "M_pep_c", "M_gtp_c", "M_pyr_c"],
    "R_ADK1" => ["M_amp_c", "M_atp_c", "M_adp_c"],
    "R_GK1" => ["M_atp_c", "M_gmp_c", "M_adp_c", "M_gdp_c"],
    "R_PPA" => ["M_ppi_c", "M_pi_c"],
]

# The eight species this module owns, in registry order.
const CONCENTRATIONS = ["M_atp_c", "M_adp_c", "M_amp_c", "M_pi_c",
                        "M_gtp_c", "M_gdp_c", "M_gmp_c", "M_ppi_c"]

const N_ROWS = 2 * length(REACTIONS) + sum(length(last(p)) for p in MICHAELIS) +
               length(CONCENTRATIONS)

"""
    balanced_rows(path) -> Dict{Tuple{String,String,String}, Tuple{String,String}}

The balanced `Parameter` table keyed on `(QuantityType, reaction, compound)`,
which is the three-part key that has no ID column and so does not fit the
loader. Values are the `Mode` and `UnconstrainedGeometricStd` cells, verbatim.

The upstream document holds six SBtab tables; only the one declaring
`TableName='Parameter'` carries balanced modes with uncertainties. The earlier
`TableName='Quantity'` table is the hand-patch layer the simulator reads for
some modules, and picking it up here would silently swap one provenance story
for another.
"""
function balanced_rows(path::AbstractString)
    rows = Dict{Tuple{String, String, String}, Tuple{String, String}}()
    in_table = false
    header = String[]
    for line in eachline(path)
        if startswith(line, "!!!")
            continue
        elseif startswith(line, "!!SBtab")
            in_table = occursin("TableName='Parameter'", line)
            header = String[]
            continue
        end
        in_table || continue
        cells = split(line, '\t')
        if startswith(line, "!")
            header = String[strip(c, '!') for c in cells]
            continue
        end
        isempty(header) && continue
        col(name) = begin
            i = findfirst(==(name), header)
            i === nothing ? "" : (i <= length(cells) ? String(strip(cells[i])) : "")
        end
        qty = col("QuantityType")
        isempty(qty) && continue
        key = (qty, col("Reaction:SBML:reaction:id"), col("Compound:SBML:species:id"))
        haskey(rows, key) && error("duplicate balanced row $key in $path")
        rows[key] = (col("Mode"), col("UnconstrainedGeometricStd"))
    end
    return rows
end

"""
    extract_rows(rows) -> Vector{NTuple{4,String}}

The 35 `(identifier, mode, gstd, upstream)` rows, in the canonical order above.
Every lookup is by an explicitly listed key, so an upstream file missing one of
them fails here rather than producing a shorter extract that would quietly
disarm the ambiguity check.
"""
function extract_rows(rows)
    out = NTuple{4, String}[]
    fetch(qty, rxn, cmp) = begin
        key = (qty, rxn, cmp)
        haskey(rows, key) || error("upstream file holds no row $key")
        mode, gstd = rows[key]
        (mode, gstd, "$qty|$rxn|$cmp")
    end
    for r in REACTIONS
        push!(out, ("kcatF_$r", fetch("substrate catalytic rate constant", r, "")...))
        push!(out, ("kcatR_$r", fetch("product catalytic rate constant", r, "")...))
    end
    for (r, species) in MICHAELIS, s in species
        push!(out, ("km_$(r)_$(s)", fetch("Michaelis constant", r, s)...))
    end
    for s in CONCENTRATIONS
        push!(out, ("conc_$s", fetch("concentration", "", s)...))
    end
    length(out) == N_ROWS || error("expected $N_ROWS rows, built $(length(out))")
    return out
end

function write_extract(path, logical_name, upstream_file, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "!!SBtab TableType='Quantity' TableName='$logical_name " *
                    "(Core A′ nucleotide-recycling extract)'")
        println(io, "!ID\t!Mode\t!GeometricStd\t!UpstreamRow")
        println(io, "% Derived from $upstream_file at Luthey-Schulten-Lab/Minimal_Cell")
        println(io, "% commit db048ac, balanced Parameter table. Regenerate with the")
        println(io, "% command in src/organisms/coreA/data/README.md. The reshape renames")
        println(io, "% identifiers and reorders rows; it changes no value.")
        for (id, mode, gstd, upstream) in rows
            println(io, "$id\t$mode\t$gstd\t$upstream")
        end
    end
    return path
end

function main(checkout)
    isdir(checkout) || error("no such Minimal_Cell checkout: $checkout")
    for (logical, rel) in UPSTREAM
        src = joinpath(checkout, rel)
        isfile(src) || error("upstream file missing: $src")
        rows = extract_rows(balanced_rows(src))
        out = write_extract(OUTPUTS[logical], logical, basename(rel), rows)
        println("wrote $out — $(length(rows)) rows from $(basename(rel))")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 ||
        error("usage: julia dev/scripts/extract_nucleotide_recycling.jl <Minimal_Cell checkout>")
    main(ARGS[1])
end
