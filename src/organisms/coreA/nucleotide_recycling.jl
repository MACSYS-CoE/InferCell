# Core A′ nucleotide recycling (spec §11 phase 8).
#
# Five reactions: the two GDP-phosphorylating analogues of PGK and PYK, which
# make GTP a live product of glycolysis, and the three that close a conserved
# moiety. Without the last three the core is a set of dead ends rather than a
# sustained pathway — the lumped charging step converts ATP to AMP and nothing
# else in Core A′ touches AMP, so 3,484,518 charging events against an adenylate
# pool of ~79,800 particles exhausts it in ~144 s of a 6,300 s cycle.
#
#     PGK3   13DPG + GDP <=> 3PG + GTP
#     PYK3   PEP   + GDP <=> pyruvate + GTP
#     ADK1   AMP   + ATP <=> 2 ADP
#     GK1    ATP   + GMP <=> ADP + GDP
#     PPA    PPi        <=> 2 Pi
#
# This module owns the eight pools every other module routes energy through,
# which is why it cannot be verified from inside itself: standalone it conserves
# both moieties trivially. Its acceptance test is a full cycle against an
# external charging drain (spec task 8.6).

"""
Directory holding the vendored parameter extracts. See its `README.md` for the
upstream files, the commit, and the command that regenerates them.
"""
const COREA_DATA_DIR = joinpath(@__DIR__, "data")

"The five reactions, in the order the scoping note introduces them."
const RECYCLING_REACTIONS = (:R_PGK3, :R_PYK3, :R_ADK1, :R_GK1, :R_PPA)

"""
The eight states this module integrates: four adenylate, three guanylate and
pyrophosphate. `M_pi_c` is in the registry's `:adenylate` group and `M_ppi_c` in
its `:other` group alongside phase 7's external lactate, so this set is listed
rather than derived from `species_in_group`.
"""
const RECYCLING_STATES = [:M_atp_c, :M_adp_c, :M_amp_c, :M_pi_c,
                          :M_gtp_c, :M_gdp_c, :M_gmp_c, :M_ppi_c]

"""
The four glycolytic species this module reads and writes but does not own.

The two products are here because the published rate laws are *reversible* and
read them twice over — once in the reverse numerator, once in the denominator's
product bracket:

    R_PGK3 ... - (keq^(-hco/2) * M_3pg_c * M_gtp_c) ...
           / ( (1+13dpg/Km)(1+gdp/Km) + (1+M_3pg_c/Km)(1+gtp/Km) - 1 )

so a declaration carrying `M_3pg_c` outbound only could not evaluate its own
rate law. Each therefore carries a mass edge in *both* directions, which is
phase 1's rule of one edge per (species, direction) that actually occurs, and
is a true statement about a reversible reaction. `contributed_states` still
names each species once and returns the net signed rate. This supersedes the
nine-edge table in `dev/archive/openspec/changes/add-nucleotide-recycling/`
design.md D7 and its task 3.2; see spec §11 task 8.5.
"""
const RECYCLING_INPUTS = [:M_13dpg_c, :M_pep_c, :M_3pg_c, :M_pyr_c]

# Kinetic parameter slots, in `parameters()` order. Ten catalytic constants,
# seventeen Michaelis constants and five enzyme concentrations. Seventeen and
# not nineteen because ADK1 and PPA each name one product species twice under
# one constant, which is what the published builder does with a stoichiometric
# coefficient above one.
const RECYCLING_KM_SPECIES = (
    :R_PGK3 => (:M_13dpg_c, :M_gdp_c, :M_3pg_c, :M_gtp_c),
    :R_PYK3 => (:M_gdp_c, :M_pep_c, :M_gtp_c, :M_pyr_c),
    :R_ADK1 => (:M_amp_c, :M_atp_c, :M_adp_c),
    :R_GK1 => (:M_atp_c, :M_gmp_c, :M_adp_c, :M_gdp_c),
    :R_PPA => (:M_ppi_c, :M_pi_c),
)

const RECYCLING_KINETIC_IDS = let ids = Symbol[]
    for r in RECYCLING_REACTIONS
        push!(ids, Symbol("kcatF_", r), Symbol("kcatR_", r))
    end
    for (r, species) in RECYCLING_KM_SPECIES, s in species
        push!(ids, Symbol("km_", r, "_", s))
    end
    for r in RECYCLING_REACTIONS
        push!(ids, Symbol("enz_", r))
    end
    Tuple(ids)
end

const N_RECYCLING_KINETIC = length(RECYCLING_KINETIC_IDS)   # 32

# Copy number per gene at the registry's initial volume, and the locus each
# comes from (design.md D5). PGK3 and PYK3 add reactions but no genes: they are
# catalysed by the same enzymes as glycolysis's PGK and PYK, so those two rows
# are nominals a composition with phase 6 must agree on.
const RECYCLING_ENZYME_COPIES = (
    :R_PGK3 => ("JCVISYN3A_0606", 411, true),
    :R_PYK3 => ("JCVISYN3A_0221", 551, true),
    :R_ADK1 => ("JCVISYN3A_0651", 213, false),
    :R_GK1 => ("JCVISYN3A_0203", 186, false),
    :R_PPA => ("JCVISYN3A_0344", 190, false),
)

"Particles per millimolar at the registry's 200 nm cell radius."
const RECYCLING_PARTICLES_PER_MM = 20180

"""
    NucleotideRecycling(; reactions, free, tables, gstd_tables)

The five-reaction nucleotide-recycling module (spec §11 phase 8).

`reactions` selects which of [`RECYCLING_REACTIONS`](@ref) are active, so a
configuration with the adenylate kinase or the pyrophosphatase removed is the
same module rather than a second one — which is what lets the removal runs of
spec task 8.6 assert against the *same* rate laws.

Every imported constant is `fixed` by default and held on the struct. Naming one
in `free` makes it an entry of `model_free_params`, and the module then reads it
from the parameter vector instead, so a sampler moves it and a fixed one costs
no posterior dimension. Nothing is inferred through a prior nobody chose.
"""
struct NucleotideRecycling <: AbstractSubModel
    params::Vector{InferParameter}
    held::NTuple{N_RECYCLING_KINETIC, Float64}
    pidx::NTuple{N_RECYCLING_KINETIC, Int}
    active::NTuple{5, Bool}
    edges::Vector{CouplingEdge}
    contribs::Vector{Symbol}
    ins::Vector{Symbol}
end

# The extracts hold a mode and a width, and `SourceTable` keeps only the mode —
# `read_source_table` reads one value column and folds the width into an
# informedness. So each extract is read twice, the second time with the width
# column *as* the value column, which is the only way to build a prior of the
# right shape without a second copy of the number in this file. `src/loader.jl`
# is frozen against a module branch (spec §10 R15), so widening `SourceTable` is
# not this phase's to do.
_recycling_tables() = [
    read_source_table(joinpath(COREA_DATA_DIR, "nucleotide_recycling.tsv");
                      file = "nucleotide_balanced"),
    read_source_table(joinpath(COREA_DATA_DIR, "nucleotide_recycling_central.tsv");
                      file = "central_balanced"),
]

_recycling_gstds() = Dict(
    "nucleotide_balanced" => read_source_table(
        joinpath(COREA_DATA_DIR, "nucleotide_recycling.tsv");
        file = "nucleotide_balanced", value_column = "GeometricStd").values,
    "central_balanced" => read_source_table(
        joinpath(COREA_DATA_DIR, "nucleotide_recycling_central.tsv");
        file = "central_balanced", value_column = "GeometricStd").values,
)

"""
    recycling_governing_file(identifier) -> String

Which balanced file governs `identifier`. The nucleotide file governs all five
reactions: the central file's rival values are not a second estimate but an
artefact of balancing with no data, and its geometric standard deviations of
1e33 to 1e63 against 1.0513 say so mechanically (spec §4 D2).

The initial conditions are not uniform, and assuming they were would be wrong
for four of the eight. Each defers to the file `registry.jl` records as that
species' source, so the registry and the extract cannot disagree about *where* a
value came from any more than the loader's own check lets them disagree about
what it is.
"""
function recycling_governing_file(identifier::AbstractString)
    if startswith(identifier, "conc_")
        species = Symbol(identifier[6:end])
        src = species_entry(species).source_file
        src === nothing && throw(ArgumentError(
            "The registry records no source file for :$species, so no governing " *
            "declaration can be derived for $identifier"))
        return src
    end
    return "nucleotide_balanced"
end

function NucleotideRecycling(; reactions = RECYCLING_REACTIONS,
                             free::AbstractVector{Symbol} = Symbol[],
                             tables = _recycling_tables(),
                             gstds = _recycling_gstds(),
                             enzyme_gstd = 1.2)
    for r in reactions
        r in RECYCLING_REACTIONS || throw(ArgumentError(
            "NucleotideRecycling knows no reaction :$r; it carries " *
            "$(join(string.(":", RECYCLING_REACTIONS), ", "))"))
    end
    for name in free
        name in RECYCLING_KINETIC_IDS || throw(ArgumentError(
            "NucleotideRecycling cannot free :$name, which is not one of its " *
            "$(N_RECYCLING_KINETIC) kinetic parameters. Initial conditions are " *
            "held by the registry and are not freed here"))
    end

    params = InferParameter[]
    held = Float64[]

    # Kinetics: the ten catalytic and seventeen Michaelis constants, imported
    # against both extracts so every one names the file it rejected.
    for id in RECYCLING_KINETIC_IDS[1:(N_RECYCLING_KINETIC - 5)]
        identifier = String(id)
        governing = recycling_governing_file(identifier)
        gstd = gstds[governing][identifier]
        p = load_parameter(tables, identifier;
                           name = id, module_id = :NucleotideRecycling,
                           prior = LogNormal(log(_table_value(tables, governing, identifier)),
                                             log(gstd)),
                           role = :rate, fixed = !(id in free),
                           governing = governing)
        push!(params, p)
        push!(held, p.value)
    end

    # Enzyme concentrations: published copy number at the registry's volume.
    #
    # The *value* is the published model's. The *width* is not: the proteomics
    # table states a count and no uncertainty, so any prior here is this
    # project's and the informedness is `:asserted` rather than `:balanced`.
    # That is a deliberate classification and it has a consequence worth
    # knowing — §4 D7 counts thirteen asserted priors across Core A′, all of
    # them transport constants and the two charging quantities, and enzyme
    # concentrations are not among them. Every metabolic module that sets an
    # enzyme concentration from a copy number adds more, so spec §11 task 13.6's
    # count is owed a reconciliation at assembly. Reported rather than filed
    # under a friendlier label.
    for r in RECYCLING_REACTIONS
        locus, copies, _ = _enzyme_row(r)
        value = copies / RECYCLING_PARTICLES_PER_MM
        name = Symbol("enz_", r)
        push!(params, InferParameter(value, LogNormal(log(value), log(enzyme_gstd)),
                                     !(name in free), name, :NucleotideRecycling,
                                     :rate,
                                     ParameterSource("published_proteomics";
                                                     table = "proteomics count",
                                                     identifier = locus,
                                                     informedness = :asserted)))
        push!(held, value)
    end

    # Initial conditions, each from the file the registry names as its source.
    for s in RECYCLING_STATES
        identifier = "conc_$s"
        governing = recycling_governing_file(identifier)
        gstd = gstds[governing][identifier]
        push!(params,
              load_parameter(tables, identifier;
                             name = Symbol(s, "0"), module_id = :NucleotideRecycling,
                             prior = LogNormal(log(_table_value(tables, governing, identifier)),
                                               log(gstd)),
                             role = :initial_condition, fixed = true,
                             governing = governing))
    end

    # Free kinetic parameters, in `model_free_params` order, are what the
    # composed vector holds; a fixed one has slot 0 and is read from `held`.
    freelist = [p.name for p in model_free_params(params)]
    pidx = Tuple(something(findfirst(==(id), freelist), 0)
                 for id in RECYCLING_KINETIC_IDS)

    return NucleotideRecycling(params,
                               Tuple(held), pidx,
                               Tuple(r in reactions for r in RECYCLING_REACTIONS),
                               _recycling_edges(), copy(RECYCLING_INPUTS),
                               copy(RECYCLING_INPUTS))
end

function _table_value(tables, file, identifier)
    i = findfirst(t -> t.file == file, tables)
    i === nothing && throw(ArgumentError(
        "No source table named \"$file\" among $(join((t.file for t in tables), ", "))"))
    return tables[i].values[identifier]
end

_enzyme_row(r::Symbol) = last(RECYCLING_ENZYME_COPIES[findfirst(p -> first(p) === r,
                                                               RECYCLING_ENZYME_COPIES)])

"""
    recycling_enzymes() -> Vector{NamedTuple}

Per-reaction locus, copy number, concentration and whether the gene is shared
with central glycolysis. Exposed so a composition can assert that one enzyme is
not run at two concentrations: PGK3 and PYK3 are catalysed by the same genes as
PGK and PYK, so two modules carry a nominal for each (spec §11 task 8.4).
"""
recycling_enzymes() = [(reaction = r, locus = locus, copies = copies,
                        concentration = copies / RECYCLING_PARTICLES_PER_MM,
                        shared_with_glycolysis = shared)
                       for (r, (locus, copies, shared)) in RECYCLING_ENZYME_COPIES]

# The boundary. Mass edges for the four species this module does not own, in
# every direction mass actually crosses; currency edges where this module is the
# *principal* producer or consumer of a pool it does own. The original rule said
# "sole", and D13 broke it: the lumped charging step of phase 9 also produces
# AMP and pyrophosphate, so recycling is no longer the only holder of those
# crossings (spec §11 task 8.5).
#
# ATP and ADP deliberately carry no edge. This module is neither their principal
# producer nor their principal consumer — glycolysis and the charging step are —
# so claiming the crossing here would put two principals on one pool.
function _recycling_edges()
    return CouplingEdge[
        MassEdge(species = :M_13dpg_c, direction = :in),
        MassEdge(species = :M_pep_c, direction = :in),
        MassEdge(species = :M_3pg_c, direction = :out),
        MassEdge(species = :M_3pg_c, direction = :in),
        MassEdge(species = :M_pyr_c, direction = :out),
        MassEdge(species = :M_pyr_c, direction = :in),
        CurrencyEdge(species = :M_amp_c, direction = :in),
        CurrencyEdge(species = :M_gmp_c, direction = :in),
        CurrencyEdge(species = :M_ppi_c, direction = :in),
        CurrencyEdge(species = :M_pi_c, direction = :out),
        CurrencyEdge(species = :M_gtp_c, direction = :out),
    ]
end

states(::NucleotideRecycling) = RECYCLING_STATES
parameters(m::NucleotideRecycling) = m.params
coupling(m::NucleotideRecycling) = m.edges
inputs(m::NucleotideRecycling) = m.ins
contributed_states(m::NucleotideRecycling) = m.contribs
formalism(::NucleotideRecycling) = :ode

# A kinetic constant: from the parameter vector where it is free, from the
# struct where it is fixed. Converting the held value to the vector's element
# type keeps both branches one type, so the right-hand side stays inferable
# under ForwardDiff.
@inline function _k(p, m::NucleotideRecycling, i::Int)
    @inbounds slot = m.pidx[i]
    return slot == 0 ? convert(eltype(p), @inbounds(m.held[i])) : @inbounds(p[slot])
end

"""
    _modular(kcatF, kcatR, E, subs, prods)

The published modular rate law (`Rxns.Enzymatic`, spec §3):

    v = E · ( kcatF·Π(Sᵢ/KmSᵢ) − kcatR·Π(Pⱼ/KmPⱼ) )
          / ( Π(1+Sᵢ/KmSᵢ) + Π(1+Pⱼ/KmPⱼ) − 1 )

`subs` and `prods` are saturation ratios with **one entry per unit of
stoichiometry**, so a coefficient of two appears as the same ratio twice and the
term is squared. That is what the published builder does, and it is why five
reactions carry seventeen Michaelis constants rather than nineteen.
"""
@inline function _modular(kcatF, kcatR, E, subs::NTuple{NS}, prods::NTuple{NP}) where {NS, NP}
    den = prod(map(x -> 1 + x, subs)) + prod(map(x -> 1 + x, prods)) - 1
    return E * (kcatF * prod(subs) - kcatR * prod(prods)) / den
end

"""
    recycling_fluxes(u, p, t, m, u_inputs) -> SVector{5}

The five reaction rates in [`RECYCLING_REACTIONS`](@ref) order, in mM/s. An
inactive reaction contributes exactly zero rather than being absent, so every
configuration has the same shape and the same rate laws.
"""
@inline function recycling_fluxes(u, p, t, m::NucleotideRecycling, u_inputs)
    atp, adp, amp, pin, gtp, gdp, gmp, ppi = u
    dpg, pep, pg3, pyr = u_inputs
    k(i) = _k(p, m, i)
    a = m.active
    return SA[
        a[1] * _modular(k(1), k(2), k(28),
                        (dpg / k(11), gdp / k(12)), (pg3 / k(13), gtp / k(14))),
        a[2] * _modular(k(3), k(4), k(29),
                        (gdp / k(15), pep / k(16)), (gtp / k(17), pyr / k(18))),
        a[3] * _modular(k(5), k(6), k(30),
                        (amp / k(19), atp / k(20)), (adp / k(21), adp / k(21))),
        a[4] * _modular(k(7), k(8), k(31),
                        (atp / k(22), gmp / k(23)), (adp / k(24), gdp / k(25))),
        a[5] * _modular(k(9), k(10), k(32),
                        (ppi / k(26),), (pin / k(27), pin / k(27))),
    ]
end

function dynamics(u, p, t, m::NucleotideRecycling, u_inputs)
    v = recycling_fluxes(u, p, t, m, u_inputs)
    pgk3, pyk3, adk1, gk1, ppa = v
    return SA[
        -adk1 - gk1,                 # ATP   consumed by ADK1 and GK1
        2adk1 + gk1,                 # ADP   two per kinase turnover, one per GK1
        -adk1,                       # AMP   the only route AMP has back
        2ppa,                        # Pi    two per pyrophosphate hydrolysed
        pgk3 + pyk3,                 # GTP   the only source of GTP in Core A′
        -pgk3 - pyk3 + gk1,          # GDP
        -gk1,                        # GMP   the only route GMP has back
        -ppa,                        # PPi
    ]
end

# One signed net term per entry of `contributed_states`, in that order. The two
# substrates are drawn down and the two products supplied; under the reversible
# law each term changes sign with the flux, which is what the second edge on
# each product declares.
function contributions(u, p, t, m::NucleotideRecycling, u_inputs)
    v = recycling_fluxes(u, p, t, m, u_inputs)
    return SA[-v[1], -v[2], v[1], v[2]]
end

export NucleotideRecycling, RECYCLING_REACTIONS, RECYCLING_STATES, RECYCLING_INPUTS,
       RECYCLING_KINETIC_IDS, recycling_enzymes, recycling_fluxes,
       recycling_governing_file, COREA_DATA_DIR
