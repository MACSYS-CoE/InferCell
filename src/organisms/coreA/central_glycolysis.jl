"""
Central glycolysis: the ten reactions from glucose-6-phosphate through lactate.

Spec §11 phase 6. The module owns the eleven glycolytic intermediates and the
redox pair — thirteen of the registry's 32 dynamic states — and reads the three
energy currencies it does not integrate.

Every kinetic constant and initial concentration is imported through
[`load_parameter`](@ref) from `src/organisms/coreA/data/central_glycolysis.tsv`,
so each carries its file, identifier, uncertainty and informedness into the
composed model. Nothing is typed into this file except the ten published copy
numbers and the reaction stoichiometry.
"""

# ---------------------------------------------------------------------------
# The reaction network
# ---------------------------------------------------------------------------

"""
The ten reactions, in pathway order, each with its substrates and products as
`species => stoichiometric coefficient`.

Participant order within a reaction is the upstream balanced table's, which is
also the order its Michaelis constants appear in — so `km_R_GAPD_M_g3p_c` and
this table cannot drift apart unnoticed.

**NOX is not here.** The published reaction list carries a terminal oxidase
regenerating NAD⁺ from NADH; Core A′ drops it, and the module registers that as
a reduction declaration rather than leaving it as a silent absence.

The coefficients are data rather than literals in a rate law because the redox
check's mutation test perturbs one of them: a conservation test that cannot be
made to fail is not evidence (spec §3).
"""
const GLYCOLYTIC_REACTIONS = [
    (id = :R_PGI, enzyme = :PGI, locus = :JCVISYN3A_0445, copies = 266,
     substrates = [:M_g6p_c => 1], products = [:M_f6p_c => 1]),
    (id = :R_PFK, enzyme = :PFK, locus = :JCVISYN3A_0220, copies = 458,
     substrates = [:M_atp_c => 1, :M_f6p_c => 1],
     products = [:M_adp_c => 1, :M_fdp_c => 1]),
    (id = :R_FBA, enzyme = :FBA, locus = :JCVISYN3A_0131, copies = 775,
     substrates = [:M_fdp_c => 1],
     products = [:M_dhap_c => 1, :M_g3p_c => 1]),
    (id = :R_TPI, enzyme = :TPI, locus = :JCVISYN3A_0727, copies = 410,
     substrates = [:M_dhap_c => 1], products = [:M_g3p_c => 1]),
    (id = :R_GAPD, enzyme = :GAPD, locus = :JCVISYN3A_0607, copies = 1355,
     substrates = [:M_g3p_c => 1, :M_nad_c => 1, :M_pi_c => 1],
     products = [:M_13dpg_c => 1, :M_nadh_c => 1]),
    (id = :R_PGK, enzyme = :PGK, locus = :JCVISYN3A_0606, copies = 411,
     substrates = [:M_13dpg_c => 1, :M_adp_c => 1],
     products = [:M_3pg_c => 1, :M_atp_c => 1]),
    (id = :R_PGM, enzyme = :PGM, locus = :JCVISYN3A_0729, copies = 322,
     substrates = [:M_3pg_c => 1], products = [:M_2pg_c => 1]),
    (id = :R_ENO, enzyme = :ENO, locus = :JCVISYN3A_0213, copies = 998,
     substrates = [:M_2pg_c => 1], products = [:M_pep_c => 1]),
    (id = :R_PYK, enzyme = :PYK, locus = :JCVISYN3A_0221, copies = 551,
     substrates = [:M_adp_c => 1, :M_pep_c => 1],
     products = [:M_atp_c => 1, :M_pyr_c => 1]),
    (id = :R_LDH_L, enzyme = :LDH_L, locus = :JCVISYN3A_0475, copies = 1100,
     substrates = [:M_nadh_c => 1, :M_pyr_c => 1],
     products = [:M_lac__L_c => 1, :M_nad_c => 1]),
]

"""
The three energy currencies this module draws on and supplies but does not
integrate. Nucleotide recycling owns them (spec §11 phase 8); until it is
composed a double holds them.
"""
const GLYCOLYSIS_CURRENCIES = [:M_atp_c, :M_adp_c, :M_pi_c]

"""
The path to the vendored extract. `src/organisms/coreA/data/README.md` records
the upstream file, the commit and the command that regenerates it.
"""
const CENTRAL_GLYCOLYSIS_TABLE =
    joinpath(@__DIR__, "data", "central_glycolysis.tsv")

# ---------------------------------------------------------------------------
# One reaction, compiled
# ---------------------------------------------------------------------------

# Substrate and product counts differ per reaction, so each reaction is its own
# type and the ten are held in a `Tuple`. `map` over a tuple is unrolled, which
# is what keeps the right-hand side type-stable and allocation-free at the
# container decision spec §11 task 1.6 re-took at 32 states.
#
# Every index is resolved once, at construction: `sidx`/`pidx` into the
# concentration vector (the thirteen owned states, then the three currencies),
# `skm`/`pkm`/`kf`/`kr` into the module's 65-value parameter table.
struct GlycolyticRate{NS, NP}
    id::Symbol
    enzyme::Int
    kf::Int
    kr::Int
    sidx::SVector{NS, Int}
    scoef::SVector{NS, Int}
    skm::SVector{NS, Int}
    pidx::SVector{NP, Int}
    pcoef::SVector{NP, Int}
    pkm::SVector{NP, Int}
end

"""
    substrate_terms(r) -> Int
    product_terms(r) -> Int

How many rate-law terms each side of a reaction contributes — one per unit of
stoichiometry, so a coefficient of two counts twice. Their sum over the ten
reactions is the number of Michaelis constants, which is how the model and the
vendored extract check each other.
"""
substrate_terms(r::GlycolyticRate) = sum(r.scoef)
product_terms(r::GlycolyticRate) = sum(r.pcoef)

# ---------------------------------------------------------------------------
# The sub-model
# ---------------------------------------------------------------------------

"""
    CentralGlycolysis(; free = Symbol[], enzyme_conc = nothing,
                        reactions = GLYCOLYTIC_REACTIONS,
                        protein_sources = default_protein_sources(),
                        table_path = CENTRAL_GLYCOLYSIS_TABLE)

The ten reactions from glucose-6-phosphate through lactate.

Owns the eleven `:glycolytic` registry species and the `:redox` pair. Reads
`M_atp_c`, `M_adp_c` and `M_pi_c` through [`inputs`](@ref) and writes back into
them through [`contributions`](@ref), so a composition that owns them sees
glycolysis draw ATP at PFK and supply it at PGK and PYK.

**Every parameter is fixed by default.** Sixty-five free parameters is not an
inference problem anyone wants to start with, and spec §4 D11 puts the first
target set at three to six. `free` names the ones to release, by parameter name
(`:kcatF_R_FBA`, `:km_R_PGI_M_g6p_c`, `:M_g6p_c0`). The module exposes the whole
vector either way, so nothing is hidden — composing it simply does not silently
produce a 65-dimensional posterior.

`enzyme_conc` overrides the ten nominal enzyme concentrations, in reaction
order. `reactions` overrides the stoichiometry, which is what lets a test mutate
one coefficient and show the redox check failing. `protein_sources` names the
ten protein counts this module declares as inputs, so translation (spec §11
phase 11) can supersede the nominal concentrations under whatever names it gives
its proteins.
"""
struct CentralGlycolysis{R <: Tuple} <: AbstractSubModel
    params::Vector{InferParameter}
    values::Vector{Float64}
    free_slot::Vector{Int}
    enzyme_conc::SVector{10, Float64}
    rates::R
    stoich::SMatrix{13, 10, Float64, 130}
    currency_stoich::SMatrix{3, 10, Float64, 30}
    protein_sources::Vector{Symbol}
    notes::Vector{String}
end

"""
    glycolytic_states() -> Vector{Symbol}

The thirteen owned states, in registry order: the eleven glycolytic
intermediates then the redox pair. Taken from the registry rather than typed, so
a registry edit cannot leave this module quietly integrating a different set.
"""
glycolytic_states() = vcat(species_in_group(:glycolytic), species_in_group(:redox))

"""
    default_protein_sources() -> Vector{Symbol}

`:P_<locus>` for each of the ten enzymes, in reaction order. A placeholder for
whatever translation names its proteins; `CentralGlycolysis` takes them as a
keyword so re-pointing the module costs no edit here.
"""
default_protein_sources() = [Symbol(:P_, r.locus) for r in GLYCOLYTIC_REACTIONS]

"""
    nominal_enzyme_concentrations() -> SVector{10, Float64}

The ten published copy numbers as concentrations at the registry's cell, in
reaction order.

**Not** the published `defaultPtnConcentration` of 0.001 mM. That is the value
the source model uses for reactions with *no* gene-protein-reaction rule, which
are decoupled from the boundary by construction; all ten reactions here have
single-gene rules, and taking it would understate every flux by 13× to 67×.

The factor is [`corea_particles_per_mM`](@ref) — 20,180.39 at the registry's
200 nm — derived rather than transcribed, so this agrees with the conversion the
handshake performs rather than with a rounded copy of it. These are nominal
stand-ins for a live count, not imported concentrations, which is why the module
registers them as a reduction declaration.
"""
nominal_enzyme_concentrations() =
    SVector{10, Float64}(r.copies / corea_particles_per_mM() for r in GLYCOLYTIC_REACTIONS)

function CentralGlycolysis(; free::AbstractVector{Symbol} = Symbol[],
                           enzyme_conc = nothing,
                           reactions = GLYCOLYTIC_REACTIONS,
                           protein_sources::AbstractVector{Symbol} = default_protein_sources(),
                           table_path::AbstractString = CENTRAL_GLYCOLYSIS_TABLE)
    length(reactions) == 10 ||
        throw(ArgumentError("Central glycolysis has ten reactions; $(length(reactions)) were given"))
    length(protein_sources) == 10 ||
        throw(ArgumentError("Central glycolysis has ten enzymes; $(length(protein_sources)) protein sources were given"))

    owned = glycolytic_states()
    # The concentration vector a rate law indexes: the thirteen owned states in
    # registry order, then the three currencies in `inputs` order.
    conc_names = vcat(owned, GLYCOLYSIS_CURRENCIES)
    cpos = Dict(s => i for (i, s) in enumerate(conc_names))

    # The extract is read twice: once for `Mode`, which is what the loader
    # imports, and once for `GeometricStd`, which it uses only to classify
    # informedness and does not return. Reading it through the same parser
    # rather than a second ad-hoc one is what keeps the prior's width and the
    # prior's median from ever coming from differently-parsed rows.
    table = read_source_table(table_path; file = CENTRAL_FILE)
    widths = read_source_table(table_path; file = CENTRAL_FILE,
                               value_column = "GeometricStd").values
    tables = [table]

    params = InferParameter[]
    values = Float64[]
    index_of = Dict{String, Int}()

    function import_value!(identifier::AbstractString; name::Symbol, role::Symbol)
        haskey(widths, identifier) || throw(ArgumentError(
            "$table_path has no geometric standard deviation for $identifier"))
        gstd = widths[identifier]
        gstd > 1 || throw(ArgumentError(
            "$identifier has geometric standard deviation $gstd; a log-normal " *
            "prior needs one strictly above 1"))
        p = load_parameter(tables, identifier;
                           name = name, module_id = :CentralGlycolysis,
                           prior = LogNormal(log(table.values[identifier]), log(gstd)),
                           role = role, fixed = !(name in free))
        push!(params, p)
        push!(values, p.value)
        index_of[identifier] = length(values)
        return length(values)
    end

    kf_of = Dict{Symbol, Int}()
    kr_of = Dict{Symbol, Int}()
    for r in reactions
        kf_of[r.id] = import_value!("kcatF_$(r.id)"; name = Symbol("kcatF_$(r.id)"), role = :rate)
    end
    for r in reactions
        kr_of[r.id] = import_value!("kcatR_$(r.id)"; name = Symbol("kcatR_$(r.id)"), role = :rate)
    end
    km_of = Dict{Tuple{Symbol, Symbol}, Int}()
    for r in reactions, (s, _) in vcat(r.substrates, r.products)
        id = "km_$(r.id)_$s"
        km_of[(r.id, s)] = import_value!(id; name = Symbol(id), role = :rate)
    end
    for s in owned
        import_value!("conc_$s"; name = Symbol(s, "0"), role = :initial_condition)
    end

    # Where each free parameter lands in the slice the orchestrator hands
    # `dynamics`. `model_free_params` is the same filter `_build_contexts` uses,
    # so this map and the composed parameter vector cannot disagree.
    free_slot = zeros(Int, length(params))
    slot = 0
    for (k, p) in enumerate(params)
        # The filter `model_free_params` applies, walked in place so a slot is a
        # position and not a search: two parameters could carry identical fields
        # and identity would then pick the wrong one.
        (p.fixed || p.role === :observation) && continue
        slot += 1
        free_slot[k] = slot
    end
    slot == length(model_free_params(params)) || error(
        "free-parameter slots disagree with model_free_params: $slot against " *
        "$(length(model_free_params(params)))")

    rate_of(r) = GlycolyticRate(
        r.id,
        findfirst(x -> x.id === r.id, reactions),
        kf_of[r.id], kr_of[r.id],
        SVector{length(r.substrates), Int}(cpos[s] for (s, _) in r.substrates),
        SVector{length(r.substrates), Int}(n for (_, n) in r.substrates),
        SVector{length(r.substrates), Int}(km_of[(r.id, s)] for (s, _) in r.substrates),
        SVector{length(r.products), Int}(cpos[s] for (s, _) in r.products),
        SVector{length(r.products), Int}(n for (_, n) in r.products),
        SVector{length(r.products), Int}(km_of[(r.id, s)] for (s, _) in r.products))
    rates = Tuple(rate_of(r) for r in reactions)

    # Stoichiometry, split by who integrates the state: `stoich` for the
    # thirteen this module owns, `currency_stoich` for the three it contributes
    # to. Both are built from the same table, so a species cannot appear in a
    # rate law and be missing from a derivative.
    net = zeros(Float64, length(conc_names), 10)
    for (k, r) in enumerate(reactions)
        for (s, n) in r.substrates
            net[cpos[s], k] -= n
        end
        for (s, n) in r.products
            net[cpos[s], k] += n
        end
    end

    ec = enzyme_conc === nothing ? nominal_enzyme_concentrations() :
         SVector{10, Float64}(enzyme_conc)

    return CentralGlycolysis(params, values, free_slot, ec, rates,
                             SMatrix{13, 10, Float64}(net[1:13, :]),
                             SMatrix{3, 10, Float64}(net[14:16, :]),
                             collect(Symbol, protein_sources),
                             _glycolysis_reduction_notes())
end

# ---------------------------------------------------------------------------
# The protocol
# ---------------------------------------------------------------------------

states(::CentralGlycolysis) = glycolytic_states()
parameters(m::CentralGlycolysis) = m.params
formalism(::CentralGlycolysis) = :ode
inference_mode(::CentralGlycolysis) = :differentiable

# The three currencies come first, in `GLYCOLYSIS_CURRENCIES` order, because
# `dynamics` indexes them positionally straight after the owned states. The ten
# protein counts follow: they are not registry species, so the resolver's
# inputs/edges check skips them, and translation supersedes the nominal enzyme
# concentrations through them (spec §11 task 6.4).
inputs(m::CentralGlycolysis) = vcat(GLYCOLYSIS_CURRENCIES, m.protein_sources)

"""
Nine edges, every `peer` unnamed.

The rule that generates them: declare an edge exactly where one of this module's
reactions touches a registry species another Core A′ module also touches, with
the direction mass flows relative to this module. Currency rather than mass for
the three energy species because the published model routes ATP traffic through
a shared pool node rather than module to module, and the edge kind exists to say
so.

`M_pyr_c` is the one that needs stating: this module both produces pyruvate at
PYK and consumes it at LDH_L, so its inbound edge does not describe its own
reactions. It describes the *crossing*, and the crossing is the
phosphotransferase cascade feeding pyruvate in.

No edge names a redox species. NAD⁺ and NADH are touched only by GAPD and
LDH_L, both inside this module, which is exactly why their sum is invariant here
and why check 3 is a per-module test rather than an assembled one.

`peer` is left unnamed throughout: six of the seven counterpart modules do not
exist yet, and the resolver fails an edge whose named peer is absent. The cost
is that a genuinely missing counterpart is not caught at composition, which is
spec §11 phase 13's assertion to make and not this one's.
"""
coupling(::CentralGlycolysis) = CouplingEdge[
    CurrencyEdge(species = :M_atp_c, direction = :in),    # PFK draws
    CurrencyEdge(species = :M_atp_c, direction = :out),   # PGK, PYK supply
    CurrencyEdge(species = :M_adp_c, direction = :in),    # PGK, PYK draw
    CurrencyEdge(species = :M_adp_c, direction = :out),   # PFK supplies
    CurrencyEdge(species = :M_pi_c, direction = :in),     # GAPD draws
    MassEdge(species = :M_g6p_c, direction = :in),        # the PTS cascade supplies
    MassEdge(species = :M_pyr_c, direction = :in),        # the PTS cascade supplies
    MassEdge(species = :M_pep_c, direction = :out),       # the PTS cascade draws
    MassEdge(species = :M_lac__L_c, direction = :out),    # lactate export draws
]

# The four inbound mass edges are on states this module integrates, so they take
# no contribution — the term is already in `dynamics`. Only the three currencies
# cross into another module's derivative.
contributed_states(::CentralGlycolysis) = GLYCOLYSIS_CURRENCIES

reduction_notes(m::CentralGlycolysis) = m.notes

# ---------------------------------------------------------------------------
# The rate law
# ---------------------------------------------------------------------------

# `Rxns.Enzymatic(substrates, products)` of the published simulator
# (`defMetRxns.py:1538`), not the SBtab `KineticLaw` column, which is inert for
# these reactions (spec §4 D1):
#
#   v = E · ( kcatF·Π(Sᵢ/KmSᵢ) − kcatR·Π(Pⱼ/KmPⱼ) )
#           / ( Π(1+Sᵢ/KmSᵢ) + Π(1+Pⱼ/KmPⱼ) − 1 )
#
# One term per unit of stoichiometry, so a coefficient of two appears squared.
@inline function _pval(m::CentralGlycolysis, p, k::Int, ::Type{T}) where {T}
    slot = m.free_slot[k]
    return slot == 0 ? T(m.values[k]) : T(p[slot])
end

@inline function _reaction_rate(r::GlycolyticRate, c, m::CentralGlycolysis, p,
                                ::Type{T}) where {T}
    num_s = one(T)
    den_s = one(T)
    for i in eachindex(r.sidx)
        x = c[r.sidx[i]] / _pval(m, p, r.skm[i], T)
        num_s *= x^r.scoef[i]
        den_s *= (one(T) + x)^r.scoef[i]
    end
    num_p = one(T)
    den_p = one(T)
    for j in eachindex(r.pidx)
        x = c[r.pidx[j]] / _pval(m, p, r.pkm[j], T)
        num_p *= x^r.pcoef[j]
        den_p *= (one(T) + x)^r.pcoef[j]
    end
    forward = _pval(m, p, r.kf, T) * num_s
    reverse = _pval(m, p, r.kr, T) * num_p
    return m.enzyme_conc[r.enzyme] * (forward - reverse) / (den_s + den_p - one(T))
end

"""
    reaction_rates(u, p, t, m::CentralGlycolysis, u_inputs) -> SVector{10}

The ten net fluxes, in `GLYCOLYTIC_REACTIONS` order, in mM/s.

Exposed because spec §4 D8 rules fluxes out of the likelihood and keeps them as
a reported diagnostic, and because the boundary assertions of spec §11 task 6.5
check the mass an outbound currency edge moves against the rate law that moves
it.
"""
function reaction_rates(u, p, t, m::CentralGlycolysis, u_inputs)
    T = promote_type(eltype(u), eltype(u_inputs), eltype(p), Float64)
    c = vcat(SVector{13, T}(u),
             SVector{3, T}(u_inputs[1], u_inputs[2], u_inputs[3]))
    return SVector{10, T}(map(r -> _reaction_rate(r, c, m, p, T), m.rates))
end

dynamics(u, p, t, m::CentralGlycolysis, u_inputs) =
    m.stoich * reaction_rates(u, p, t, m, u_inputs)

# The currencies' terms come from the same ten rates. They are recomputed rather
# than shared: the orchestrator calls `dynamics` and `contributions` as two
# independent functions of (u, p, t, u_inputs), and caching across them would
# mean holding solver state on the module, which is the aliasing hazard the
# out-of-place contract exists to avoid. The cost is one extra rate-law pass per
# right-hand side, which spec §11 task 13.7 is where it gets measured.
contributions(u, p, t, m::CentralGlycolysis, u_inputs) =
    m.currency_stoich * reaction_rates(u, p, t, m, u_inputs)

# ---------------------------------------------------------------------------
# What is ours
# ---------------------------------------------------------------------------

# The eight Michaelis constants where the balanced `Parameter` table and the
# `Quantity` table the simulator reads disagree: identifier, what runs, what is
# imported here. Two of them move by 82× and 227×.
const KM_COLUMN_DISAGREEMENTS = [
    ("km_R_PGI_M_g6p_c", 0.28, 22.9419),
    ("km_R_PGI_M_f6p_c", 0.15, 3.3488),
    ("km_R_PFK_M_atp_c", 0.117, 0.0255),
    ("km_R_FBA_M_fdp_c", 0.005, 0.2447),
    ("km_R_FBA_M_dhap_c", 0.095, 0.0064),
    ("km_R_FBA_M_g3p_c", 1.0, 0.0044),
    ("km_R_GAPD_M_nad_c", 1.3, 5.6969),
    ("km_R_PGM_M_2pg_c", 1.47, 0.0278),
]

function _glycolysis_reduction_notes()
    km_list = join(("$id runs at $ran and is imported at $imported"
                    for (id, ran, imported) in KM_COLUMN_DISAGREEMENTS), "; ")
    return [
        "the terminal oxidase NOX is dropped from the reaction list: LDH_L " *
        "already regenerates NAD+ from NADH so NOX is redundant for redox " *
        "balance, it carries no gene-protein-reaction rule so cutting it " *
        "crosses no boundary, and its catalytic constant was prior-dominated " *
        "so it was never an inference target; this is a departure from the " *
        "published reaction list rather than from the published model's " *
        "behaviour",

        "the 32 Michaelis constants of R_PGI, R_PFK, R_FBA, R_TPI, R_GAPD, " *
        "R_PGK, R_PGM, R_ENO, R_PYK and R_LDH_L are imported from the " *
        "balanced Parameter table rather than the Quantity table the " *
        "published simulator reads, because only the balanced column carries " *
        "an inherited prior; eight of the 32 differ, one by 82 times and one " *
        "by 227 times ($km_list, all in mM)",

        "M_atp_c, M_adp_c and M_pi_c are read and contributed to but not " *
        "integrated here, so any trajectory produced without the nucleotide " *
        "recycling module composed holds them at whatever the double supplies " *
        "and is not a closed energy loop",

        "the ten enzyme concentrations of R_PGI, R_PFK, R_FBA, R_TPI, " *
        "R_GAPD, R_PGK, R_PGM, R_ENO, R_PYK and R_LDH_L are nominal " *
        "stand-ins derived from published copy number at the registry's cell " *
        "volume, not measured concentrations, and translation supersedes them",
    ]
end

export CentralGlycolysis, GLYCOLYTIC_REACTIONS, GLYCOLYSIS_CURRENCIES,
       CENTRAL_GLYCOLYSIS_TABLE, KM_COLUMN_DISAGREEMENTS,
       glycolytic_states, default_protein_sources, nominal_enzyme_concentrations,
       reaction_rates, substrate_terms, product_terms
