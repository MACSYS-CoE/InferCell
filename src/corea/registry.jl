"""
Core A′ species registry.

The single place the reduced-syn3A species are named and ordered. Every Core A′
sub-model indexes into this registry; none of them defines its own species names.
The contents are the state list of `dev/notes/reduced-syn3a-scoping.md` — 32
dynamic states and 5 chemostats — and the per-group cardinalities are asserted by
the test suite so that the note and the code cannot drift apart silently.

Names are the published model's BiGG-style identifiers (`M_g6p_c`, `M_13dpg_c`,
…) rather than the display names of the note's tables. They map one-to-one onto
the source model, and every one of them is a valid Julia identifier, which the
display names (`13DPG`, `3PG`) are not.

H2O and H+ are deliberately absent: Core A′ derives from the `NoH2O` model with
hydrogen-ion accounting removed, which is also why inorganic pyrophosphatase is
carried as `PPi -> 2 Pi`.
"""

const SPECIES_GROUPS = (:glycolytic, :adenylate, :guanylate, :redox, :other,
                        :trna, :pts, :chemostat)

const SPECIES_TREATMENTS = (:dynamic, :chemostatted)

"""
Copy-number regimes, from the scoping note's three-regime table. The regime is
what makes a formalism defensible for a species, so it is registry data rather
than a modelling choice made later: mRNA sits at 0–2 copies and requires a
discrete method, protein at 266–1355, metabolites at ~2,000–74,000.

`:mrna` is unused by the current registry — Core A′'s mRNA lives in the
stochastic block, not in this state list — and is defined here because the
transcription and decay modules will index the same vocabulary.
"""
const SPECIES_REGIMES = (:mrna, :protein, :metabolite)

"""
    SpeciesEntry

One row of the Core A′ species registry.

`initial_value` is in mM and `gstd` is the geometric standard deviation of the
balancing distribution it came from; both are `nothing` where nothing has been
imported. `source_file` records which file the value was read from — the
central-versus-nucleotide distinction that has already produced two errors of
record, so it is carried per entry rather than resolved once by hand.
"""
struct SpeciesEntry
    name::Symbol
    group::Symbol
    treatment::Symbol
    regime::Symbol
    initial_value::Union{Float64, Nothing}
    gstd::Union{Float64, Nothing}
    source_file::Union{String, Nothing}
    informedness::Symbol

    # An inner constructor, so that validation cannot be bypassed and so that
    # this does not collide with the implicit one — an eight-argument outer
    # constructor would have the same signature and overwrite it.
    function SpeciesEntry(name, group, treatment, regime, initial_value, gstd,
                          source_file, informedness)
        group in SPECIES_GROUPS || throw(ArgumentError(
            "Species :$name has invalid group :$group. Must be one of $SPECIES_GROUPS"))
        treatment in SPECIES_TREATMENTS || throw(ArgumentError(
            "Species :$name has invalid treatment :$treatment. Must be one of $SPECIES_TREATMENTS"))
        regime in SPECIES_REGIMES || throw(ArgumentError(
            "Species :$name has invalid regime :$regime. Must be one of $SPECIES_REGIMES"))
        informedness in INFORMEDNESS || throw(ArgumentError(
            "Species :$name has invalid informedness :$informedness. Must be one of $INFORMEDNESS"))
        return new(Symbol(name), group, treatment, regime,
                   initial_value === nothing ? nothing : Float64(initial_value),
                   gstd === nothing ? nothing : Float64(gstd),
                   source_file, informedness)
    end
end

# Source files are the balanced SBtab tables of the published model, referred to
# by the names the scoping note gives them. The tables are not vendored here —
# wave 0 delivers the loader, wave 1 imports through it — so these are logical
# source identifiers that the loader resolves to paths, not filenames on disk.
const CENTRAL_FILE = "central_balanced"
const NUCLEOTIDE_FILE = "nucleotide_balanced"

"""
    COREA_SPECIES

The registry: 32 dynamic states then 5 chemostats, in the scoping note's group
order. Position in this vector is the canonical ordering every module indexes
against, so entries are appended within their group rather than reordered.
"""
const COREA_SPECIES = SpeciesEntry[
    # --- Glycolytic intermediates (11) -------------------------------------
    # Initial conditions from the central balanced file's `Mode` column.
    SpeciesEntry(:M_g6p_c,      :glycolytic, :dynamic, :metabolite,  3.7076, 1.28, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_f6p_c,      :glycolytic, :dynamic, :metabolite,  0.8538, 1.92, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_fdp_c,      :glycolytic, :dynamic, :metabolite,  7.6037, 1.14, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_dhap_c,     :glycolytic, :dynamic, :metabolite,  0.6445, 2.10, CENTRAL_FILE, :balanced),
    # G3P is one of the three states nothing informed: prior median, prior width.
    SpeciesEntry(:M_g3p_c,      :glycolytic, :dynamic, :metabolite,  0.1,   10.0,  CENTRAL_FILE, :prior_default),
    SpeciesEntry(:M_13dpg_c,    :glycolytic, :dynamic, :metabolite,  0.0098, 6.02, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_3pg_c,      :glycolytic, :dynamic, :metabolite,  1.1015, 1.78, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_2pg_c,      :glycolytic, :dynamic, :metabolite,  0.0272, 4.88, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_pep_c,      :glycolytic, :dynamic, :metabolite,  0.0409, 4.46, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_pyr_c,      :glycolytic, :dynamic, :metabolite,  3.3660, 1.31, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_lac__L_c,   :glycolytic, :dynamic, :metabolite,  0.1,   10.0,  CENTRAL_FILE, :prior_default),

    # --- Adenylate (4) ------------------------------------------------------
    # Pi is grouped with the adenylate species because the scoping note's state
    # list groups it there. Phosphate closure spans this grouping — it counts Pi
    # plus every phosphorylated species — and belongs to add-nucleotide-recycling.
    SpeciesEntry(:M_atp_c,      :adenylate,  :dynamic, :metabolite,  3.6529, 1.28, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_adp_c,      :adenylate,  :dynamic, :metabolite,  0.2178, 2.90, CENTRAL_FILE, :balanced),
    # AMP is a nucleotide-module species: the nucleotide file governs, not the
    # central file's 0.1 mM prior default.
    SpeciesEntry(:M_amp_c,      :adenylate,  :dynamic, :metabolite,  0.0832, 3.76, NUCLEOTIDE_FILE, :balanced),
    SpeciesEntry(:M_pi_c,       :adenylate,  :dynamic, :metabolite, 17.8185, 1.06, CENTRAL_FILE, :balanced),

    # --- Guanylate (3) ------------------------------------------------------
    # All three from the nucleotide file. Reading these from the central file is
    # the cross-file trap: it shows 0.1 mM at gstd 10 because they are out of
    # that module's scope, which understates the guanylate pool by an order of
    # magnitude (~39,800 particles against the ~4,000 that 0.1 mM implies).
    SpeciesEntry(:M_gtp_c,      :guanylate,  :dynamic, :metabolite,  1.6627, 1.57, NUCLEOTIDE_FILE, :balanced),
    SpeciesEntry(:M_gdp_c,      :guanylate,  :dynamic, :metabolite,  0.2981, 2.65, NUCLEOTIDE_FILE, :balanced),
    SpeciesEntry(:M_gmp_c,      :guanylate,  :dynamic, :metabolite,  0.0117, 5.82, NUCLEOTIDE_FILE, :balanced),

    # --- Redox (2) ----------------------------------------------------------
    SpeciesEntry(:M_nad_c,      :redox,      :dynamic, :metabolite,  2.1844, 1.45, CENTRAL_FILE, :balanced),
    SpeciesEntry(:M_nadh_c,     :redox,      :dynamic, :metabolite,  0.0253, 4.96, CENTRAL_FILE, :balanced),

    # --- Other (2) ----------------------------------------------------------
    # PPi is the second of the three uninformed states. It barely matters once
    # PPA is present: the pool turns over ~1,700 times per cycle, so any starting
    # value is forgotten within seconds.
    SpeciesEntry(:M_ppi_c,      :other,      :dynamic, :metabolite,  0.1,   10.0,  CENTRAL_FILE, :prior_default),
    # External lactate: an export product with no imported initial condition.
    SpeciesEntry(:M_lac__L_e,   :other,      :dynamic, :metabolite,  nothing, nothing, nothing, :not_imported),

    # --- tRNA (2) -----------------------------------------------------------
    # One effective charged-tRNA pool, replacing the 20 per-amino-acid chains.
    # This lumping is ours, not the published model's; the module that carries it
    # registers the fact via `reduction_notes`, so `reduction_declarations`
    # enumerates it.
    SpeciesEntry(:M_trna_c,     :trna,       :dynamic, :metabolite, nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_trna_chg_c, :trna,       :dynamic, :metabolite, nothing, nothing, nothing, :not_imported),

    # --- PTS phospho-states (8) --------------------------------------------
    # Total copies per cell are known (ptsI 353, Crr 314, ptsH 290, ptsG 831) but
    # the split across phosphorylated and unphosphorylated forms is not, so no
    # initial condition is imported for either form. PTS protein conservation —
    # each pair conserved across its two forms — is add-pts-transport's check.
    SpeciesEntry(:M_ptsi_c,     :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_ptsi_P_c,   :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_ptsh_c,     :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_ptsh_P_c,   :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_crr_c,      :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_crr_P_c,    :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_ptsg_c,     :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_ptsg_P_c,   :pts,        :dynamic, :protein,    nothing, nothing, nothing, :not_imported),

    # --- Chemostats (5) -----------------------------------------------------
    # Held constant rather than integrated. CTP, UTP and the amino-acid pool are
    # chemostatted by this reduction, not by the published model; the medium
    # concentrations of external glucose and O2 are the model's own.
    SpeciesEntry(:M_glc__D_e,   :chemostat,  :chemostatted, :metabolite, 40.0, nothing, nothing, :asserted),
    SpeciesEntry(:M_ctp_c,      :chemostat,  :chemostatted, :metabolite, nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_utp_c,      :chemostat,  :chemostatted, :metabolite, nothing, nothing, nothing, :not_imported),
    SpeciesEntry(:M_aa_pool_c,  :chemostat,  :chemostatted, :metabolite,  0.1, nothing, nothing, :asserted),
    SpeciesEntry(:M_o2_c,       :chemostat,  :chemostatted, :metabolite, nothing, nothing, nothing, :not_imported),
]

"""
    COREA_SPECIES_INDEX

Name-to-position map over [`COREA_SPECIES`](@ref), built once at load. Position
is the canonical ordering, so this is the only lookup a module needs.
"""
const COREA_SPECIES_INDEX = Dict{Symbol, Int}(
    e.name => i for (i, e) in enumerate(COREA_SPECIES))

# A duplicated name would make the index silently lose an entry, which is
# exactly the class of error the registry exists to prevent. Catch it at load.
length(COREA_SPECIES_INDEX) == length(COREA_SPECIES) || error(
    "Duplicate species name in COREA_SPECIES: the registry has " *
    "$(length(COREA_SPECIES)) entries but only $(length(COREA_SPECIES_INDEX)) distinct names")

"""
    species_index(name::Symbol) -> Int

Canonical registry position of `name`. Throws naming the unrecognised symbol
rather than returning a default, so a typo in a module's coupling declaration
fails at composition instead of silently indexing the wrong state.
"""
function species_index(name::Symbol)
    idx = get(COREA_SPECIES_INDEX, name, nothing)
    idx === nothing && throw(ArgumentError(
        "Unknown Core A′ species :$name. It is not in the registry of " *
        "$(length(COREA_SPECIES)) species; see src/corea/registry.jl"))
    return idx
end

"""
    species_entry(name::Symbol) -> SpeciesEntry
    species_entry(i::Integer)   -> SpeciesEntry

The registry row for a species, by name or by canonical position.
"""
species_entry(name::Symbol) = COREA_SPECIES[species_index(name)]

function species_entry(i::Integer)
    checkbounds(Bool, COREA_SPECIES, i) || throw(BoundsError(COREA_SPECIES, i))
    return COREA_SPECIES[i]
end

"""
    is_registered(name::Symbol) -> Bool

Whether `name` is in the registry. The non-throwing companion to
[`species_index`](@ref), for callers that want to report several unknown species
at once rather than failing on the first.
"""
is_registered(name::Symbol) = haskey(COREA_SPECIES_INDEX, name)

"""
    species_group(name::Symbol) -> Symbol

The conserved-moiety group a species belongs to. This is what lets the dead-end
report name a moiety: AMP with no route back and GMP with no route back are the
same class of error, and reporting them as adenylate and guanylate says so.
"""
species_group(name::Symbol) = species_entry(name).group

"""
    is_chemostatted(name::Symbol) -> Bool

Whether the registry holds this species at a fixed concentration. A module that
declares dynamics for a chemostatted species is rejected at composition.
"""
is_chemostatted(name::Symbol) = species_entry(name).treatment === :chemostatted

"""
    is_dynamic(name::Symbol) -> Bool

Whether this species is integrated rather than held.
"""
is_dynamic(name::Symbol) = species_entry(name).treatment === :dynamic

"""
    dynamic_species() -> Vector{Symbol}
    chemostat_species() -> Vector{Symbol}

The registry's dynamic states and its chemostats, in canonical order.
"""
dynamic_species() = [e.name for e in COREA_SPECIES if e.treatment === :dynamic]
chemostat_species() = [e.name for e in COREA_SPECIES if e.treatment === :chemostatted]

"""
    species_in_group(group::Symbol) -> Vector{Symbol}

Registry members of one group, in canonical order.
"""
function species_in_group(group::Symbol)
    group in SPECIES_GROUPS || throw(ArgumentError(
        "Unknown species group :$group. Must be one of $SPECIES_GROUPS"))
    return [e.name for e in COREA_SPECIES if e.group === group]
end

"""
    held_value(name::Symbol) -> Union{Float64, Nothing}

The concentration a chemostatted species is held at, in mM. Throws for a dynamic
species, since asking a dynamic state for its held value is a category error
rather than a missing lookup. Returns `nothing` where the species is chemostatted
but no value has been imported.
"""
function held_value(name::Symbol)
    e = species_entry(name)
    e.treatment === :chemostatted || throw(ArgumentError(
        "Species :$name is dynamic, not chemostatted, so it has no held value. " *
        "Use species_entry(:$name).initial_value for its initial condition"))
    return e.initial_value
end

"""
    is_informed(name::Symbol) -> Bool

Whether this species' value came from a balancing distribution. False for the
three states nothing informed — G3P, PPi and cytosolic lactate, which sit at the
0.1 mM prior median with a geometric standard deviation of 10 — and for anything
not yet imported. Worth knowing before treating a starting point as trusted.
"""
is_informed(name::Symbol) = species_entry(name).informedness === :balanced

"""
    uninformed_species() -> Vector{Symbol}

Every registry species whose value is not a balanced one, in canonical order.
"""
uninformed_species() = [e.name for e in COREA_SPECIES if e.informedness !== :balanced]

"""
    n_dynamic_states() -> Int
    n_chemostats() -> Int

Registry cardinalities. The scoping note fixes these at 32 and 5; the test suite
asserts them, so a change to the note that is not carried into the registry fails
the suite rather than accumulating.
"""
n_dynamic_states() = count(e -> e.treatment === :dynamic, COREA_SPECIES)
n_chemostats() = count(e -> e.treatment === :chemostatted, COREA_SPECIES)
