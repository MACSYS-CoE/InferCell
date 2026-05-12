struct TranscriptionTranslation <: AbstractSubModel
    params::Vector{InferParameter}
    gene_names::Vector{Symbol}  # [:_default] for single-gene; explicit names for multi-gene
    enzyme_slot::Union{Nothing, Symbol}  # which gene's protein feeds metabolism
end

function TranscriptionTranslation(;
        k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1,
        mRNA0=0.0, protein0=0.0, sigma_obs=0.3)
    params = [
        InferParameter(k_tx,          LogNormal(0, 1),                     false, :k_tx,          :txl, :rate),
        InferParameter(k_tl,          LogNormal(0, 1),                     false, :k_tl,          :txl, :rate),
        InferParameter(gamma_mRNA,    LogNormal(0, 1),                     false, :gamma_mRNA,    :txl, :rate),
        InferParameter(gamma_protein, LogNormal(0, 1),                     false, :gamma_protein, :txl, :rate),
        InferParameter(mRNA0,         Normal(0, 1),                        true,  :mRNA0,         :txl, :initial_condition),
        InferParameter(protein0,      Normal(0, 1),                        true,  :protein0,      :txl, :initial_condition),
        InferParameter(sigma_obs,     truncated(Normal(0, 1); lower=0.0),  false, :sigma_obs,     :txl, :observation),
    ]
    return TranscriptionTranslation(params, [:_default], nothing)
end

# Multi-gene constructor. `genes` lists gene slot names (e.g. [:enzyme, :ribosome, :reporter]).
# `enzyme_slot` (default :enzyme if present in genes) marks which protein feeds metabolism.
# Per-gene rates and ICs default from the kwargs but can be overridden via `overrides`,
# a Dict mapping gene name → NamedTuple of parameter values.
function TranscriptionTranslation(genes::Vector{Symbol};
        k_tx=1.0, k_tl=2.0, gamma_mRNA=0.5, gamma_protein=0.1,
        mRNA0=0.0, protein0=0.0, sigma_obs=0.3,
        enzyme_slot=(:enzyme in genes ? :enzyme : nothing),
        overrides::Dict{Symbol, <:NamedTuple}=Dict{Symbol, NamedTuple}())
    length(genes) >= 1 || error("genes must be non-empty")
    if enzyme_slot !== nothing
        enzyme_slot in genes || error("enzyme_slot :$enzyme_slot not in genes=$genes")
    end

    params = InferParameter[]
    for g in genes
        ov = get(overrides, g, NamedTuple())
        gk_tx   = get(ov, :k_tx, k_tx)
        gk_tl   = get(ov, :k_tl, k_tl)
        gg_m    = get(ov, :gamma_mRNA, gamma_mRNA)
        gg_p    = get(ov, :gamma_protein, gamma_protein)
        gm0     = get(ov, :mRNA0, mRNA0)
        gp0     = get(ov, :protein0, protein0)
        push!(params,
            InferParameter(gk_tx, LogNormal(0, 1), false, Symbol("k_tx_", g),          :txl, :rate),
            InferParameter(gk_tl, LogNormal(0, 1), false, Symbol("k_tl_", g),          :txl, :rate),
            InferParameter(gg_m,  LogNormal(0, 1), false, Symbol("gamma_mRNA_", g),    :txl, :rate),
            InferParameter(gg_p,  LogNormal(0, 1), false, Symbol("gamma_protein_", g), :txl, :rate),
            InferParameter(gm0,   Normal(0, 1),    true,  Symbol(g, "_mRNA0"),         :txl, :initial_condition),
            InferParameter(gp0,   Normal(0, 1),    true,  Symbol(g, "_protein0"),      :txl, :initial_condition),
        )
    end
    push!(params,
        InferParameter(sigma_obs, truncated(Normal(0, 1); lower=0.0), false, :sigma_obs, :txl, :observation))

    return TranscriptionTranslation(params, copy(genes), enzyme_slot)
end

states(m::TranscriptionTranslation) = m.gene_names == [:_default] ?
    [:mRNA, :protein] :
    reduce(vcat, [[Symbol(g, "_mRNA"), Symbol(g, "_protein")] for g in m.gene_names])

parameters(m::TranscriptionTranslation) = m.params
formalism(::TranscriptionTranslation) = :ode
inference_mode(::TranscriptionTranslation) = :differentiable

enzyme_protein_state(m::TranscriptionTranslation) =
    m.enzyme_slot === nothing ? :protein : Symbol(m.enzyme_slot, "_protein")

function dynamics(u, p, t, m::TranscriptionTranslation)
    if m.gene_names == [:_default]
        mRNA, protein = u
        k_tx, k_tl, gamma_mRNA, gamma_protein = p
        return SA[k_tx - gamma_mRNA * mRNA, k_tl * mRNA - gamma_protein * protein]
    else
        n = length(m.gene_names)
        # Per-gene params packed in groups of 4: (k_tx, k_tl, γ_mRNA, γ_protein) for each gene.
        # State packed in pairs of 2: (mRNA, protein) per gene.
        du = Vector{eltype(u)}(undef, 2n)
        for i in 1:n
            mRNA = u[2i - 1]
            protein = u[2i]
            k_tx  = p[4(i - 1) + 1]
            k_tl  = p[4(i - 1) + 2]
            g_m   = p[4(i - 1) + 3]
            g_p   = p[4(i - 1) + 4]
            du[2i - 1] = k_tx - g_m * mRNA
            du[2i]     = k_tl * mRNA - g_p * protein
        end
        return SVector{2n}(du)
    end
end
