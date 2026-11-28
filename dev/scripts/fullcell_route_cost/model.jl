# KIM-195: a mean field of the whole syn3A cell, generic in its number type so that
# ForwardDiff duals pass through it. Included by route_cost_meanfield.jl.
#
# What it is. The CME and the metabolism ODE of the full cell (syn3a-ensemble
# data/model/, the export SpeedyCell runs), as one deterministic ODE on the 2,661 CME
# counts, solved one minute at a time:
#   - first-order reactions (LM type 1): propensity k·x;
#   - second-order reactions (type 2, 363): k·x_a·x_b, the product of means (not exact);
#   - the 175 ODE fluxes (mM/s, from data/model/params/ode_reactions.csv, whose named
#     kinetic constants are tokens K__i), turned into counts/s at the cell volume V;
#   - the hook's cost counters, kept as states and paid into their pools at rate
#     κ·c·p/(p + 1), κ = 1 s⁻¹ (the hook's cadence). The hook pays min(counter, pool)
#     each second (writeback.jl): the whole counter while the pool lasts. This pays
#     the full κ·c while p ≫ 1 count and stops as the pool empties, leaving the unpaid
#     rest on the counter (Tom, 2026-11-27: instant payment with no clamp drove dATP
#     to −2,887 at θ = 0; 2026-11-28: this form, not κ·c·p/(c + p + 1), which pays
#     only about min/2 when counter and pool are comparable);
#     the NMP counters, which the hook pays without a clamp, at κ·c;
#   - at every minute boundary the CME rate constants (rate_rules.csv through
#     nodes.csv), the enzyme concentrations p, the volume V and the surface areas are
#     rebuilt from the current mean state, then held for the minute. SpeedyCell
#     rebuilds the CME rates every minute too; it rebuilds p, V and the ODE state every
#     1 s (the hook), which this model does every 60 s instead.
# The node DAG is evaluated "smoothly" in a run: `round` is the identity, so the one
# rounded node (CellV, read by the replication rules) carries a derivative. In the
# check stage it is evaluated with rounding, against SpeedyCell's own values.
#
# The parameters θ. A block of B log-multipliers, on CME reaction rates (applied after
# each rebuild, so the rebuild cannot overwrite them: placement GE) or on named ODE
# constants K__i (kcatF: placement MET), or both (BOTH). frame.md §1.2.

using SparseArrays
using ForwardDiff

const MODEL = "/fred/oz022/tkimpson/syn3a-ensemble/data/model"

# The hook's constants and lists (SpeedyCell src/Cell/writeback.jl, which copies
# Minimal_Cell CME_ODE/program's in_out).
const NA = 6.022e23
const V_CAP = 6.70e-17
const CELLV_CAP = 670.0
const FRAC_OUT_MEM = 0.513
const AVG_PROT_SA = 28.0
const SA_DICT = [("M_clpn_c", 0.4), ("M_chsterol_c", 0.35), ("M_sm_c", 0.45), ("M_pc_c", 0.55),
                 ("M_pg_c", 0.6), ("M_galfur12dgr_c", 0.6), ("M_12dgr_c", 0.5), ("M_pa_c", 0.5),
                 ("M_cdpdag_c", 0.5)]
const MEM_PROT_IDS = [
    "0005", "0008", "0009", "0010", "0011", "0030", "0034", "0060", "0095", "0113", "0114",
    "0116", "0117", "0132", "0143", "0146", "0164", "0165", "0166", "0167", "0168", "0169",
    "0195", "0196", "0197", "0235", "0239", "0248", "0249", "0296", "0304", "0314", "0317",
    "0326", "0332", "0338", "0345", "0346", "0371", "0372", "0379", "0388", "0398", "0399",
    "0411", "0425", "0426", "0427", "0428", "0439", "0440", "0478", "0481", "0505", "0516",
    "0601", "0639", "0641", "0642", "0643", "0652", "0685", "0686", "0691", "0696", "0706",
    "0707", "0708", "0774", "0777", "0778", "0779", "0787", "0789", "0790", "0791", "0792",
    "0795", "0797", "0822", "0827", "0830", "0835", "0836", "0839", "0852", "0870", "0872",
    "0876", "0878", "0879", "0881", "0908"]
mem_prot_species(id) = id == "0779" ? ["ptsg", "ptsg_P"] : ["M_PTN_JCVISYN3A_$(id)_c"]
const ATP_HYDRO_COUNTERS = ["ATP_translat", "ATP_trsc", "ATP_mRNAdeg", "ATP_DNArep", "ATP_transloc"]
const NTP_COUNTERS = [("ATP_mRNA", "M_atp_c"), ("CTP_mRNA", "M_ctp_c"), ("UTP_mRNA", "M_utp_c"),
                      ("GTP_mRNA", "M_gtp_c"), ("ATP_tRNA", "M_atp_c"), ("CTP_tRNA", "M_ctp_c"),
                      ("UTP_tRNA", "M_utp_c"), ("GTP_tRNA", "M_gtp_c"), ("ATP_rRNA", "M_atp_c"),
                      ("CTP_rRNA", "M_ctp_c"), ("UTP_rRNA", "M_utp_c"), ("GTP_rRNA", "M_gtp_c")]
const NMP_COUNTERS = [("AMP_mRNAdeg", "M_amp_c"), ("UMP_mRNAdeg", "M_ump_c"),
                      ("CMP_mRNAdeg", "M_cmp_c"), ("GMP_mRNAdeg", "M_gmp_c")]
const DNTP_COUNTERS = [("dATP_DNArep", "M_datp_c"), ("dTTP_DNArep", "M_dttp_c"),
                       ("dCTP_DNArep", "M_dctp_c"), ("dGTP_DNArep", "M_dgtp_c")]
const AATRNA_COUNTERS = [("FMET_cost", "M_fmettrna_c", "M_trnamet_c")]

"""Rows of a CSV with a header and no quoting; `limit` keeps a free-text last column whole."""
function csvrows(path; limit = 0)
    lines = readlines(path)
    ncol = length(split(lines[1], ','))
    [split(l, ','; limit = limit > 0 ? limit : ncol) for l in lines[2:end] if !isempty(l)]
end

# ---------------------------------------------------------------------------------------
# The node DAG (nodes.csv), evaluated generically

struct DAG
    op::Vector{Symbol}
    args::Vector{Vector{Int}}
    lit::Vector{Float64}           # const: the literal
    cnt::Vector{Int}               # count: the CME species index
    ids::Vector{Int}               # the nodes this model needs, ascending (topological)
end

function load_dag(names, roots)
    rows = csvrows(joinpath(MODEL, "nodes.csv"); limit = 3)
    n = length(rows)
    sidx = Dict(s => i for (i, s) in enumerate(names))
    op = Vector{Symbol}(undef, n); args = [Int[] for _ in 1:n]
    lit = zeros(n); cnt = zeros(Int, n)
    for (i, r) in enumerate(rows)
        parse(Int, r[1]) == i || error("nodes.csv not in id order")
        o = Symbol(r[2]); op[i] = o
        if o === :const
            lit[i] = parse(Float64, r[3])
        elseif o === :count
            cnt[i] = sidx[r[3]]
        else
            args[i] = parse.(Int, split(r[3], ';'))
        end
    end
    seen = falses(n); stack = collect(roots)
    while !isempty(stack)
        i = pop!(stack); seen[i] && continue
        seen[i] = true; append!(stack, args[i])
    end
    ops = Set(op[findall(seen)])
    ops ⊆ Set([:const, :count, :add, :sub, :mul, :div, :pow, :neg, :abs, :round, :min, :max,
               :select_gt]) || error("unhandled node ops: $(setdiff(ops, Set([:const, :count])))")
    return DAG(op, args, lit, cnt, findall(seen))
end

"""Evaluate the DAG's needed nodes at counts `x` into `v`. `smooth`: `round` is the
identity (a run); otherwise Julia's round half to even, as Python's (the check)."""
function eval_dag!(v::AbstractVector{T}, d::DAG, x, smooth::Bool) where {T}
    @inbounds for i in d.ids
        o = d.op[i]; a = d.args[i]
        v[i] = if o === :const
            T(d.lit[i])
        elseif o === :count
            T(x[d.cnt[i]])
        elseif o === :add
            v[a[1]] + v[a[2]]
        elseif o === :sub
            v[a[1]] - v[a[2]]
        elseif o === :mul
            v[a[1]] * v[a[2]]
        elseif o === :div
            v[a[1]] / v[a[2]]
        elseif o === :pow
            v[a[1]]^v[a[2]]
        elseif o === :neg
            -v[a[1]]
        elseif o === :abs
            abs(v[a[1]])
        elseif o === :round
            smooth ? v[a[1]] : T(round(ForwardDiff.value(v[a[1]])))
        elseif o === :min
            m = v[a[1]]
            for k in 2:length(a)
                v[a[k]] < m && (m = v[a[k]])
            end
            m
        elseif o === :max
            m = v[a[1]]
            for k in 2:length(a)
                v[a[k]] > m && (m = v[a[k]])
            end
            m
        else # select_gt
            v[a[1]] > v[a[2]] ? v[a[3]] : v[a[4]]
        end
    end
    return v
end

# ---------------------------------------------------------------------------------------
# The ODE fluxes, generated from params/ode_reactions.csv as SpeedyCell's ODE/params.jl
# does (bind_params), generic in the number type

struct ODEPart
    ids::Vector{String}            # ODE species, in odecell order
    cme::Vector{Int}               # their CME indices
    enzyme_nodes::Vector{Int}      # p[k] is node enzyme_nodes[k]
    init_nodes::Vector{Int}        # y0's nodes (the check)
    k_names::Vector{String}        # named constants K__i, "reaction.key"
    k_default::Vector{Float64}
    reads::Vector{Vector{Int}}     # flux j reads these ODE species
    S::SparseMatrixCSC{Float64, Int}   # ODE species × reactions
    reaction_ids::Vector{String}
    flux!::Function                # flux!(V, y, p, k)
end

function bind_params(ex, slot, param, enz)
    if ex isa Symbol
        s = String(ex)
        haskey(slot, s) && return :(y[$(slot[s])])
        haskey(param, s) && return param[s]
        m = match(r"^N__(\d+)$", s)
        m === nothing || return :(p[$(enz[parse(Int, m.captures[1])])])
        m = match(r"^K__(\d+)$", s)
        m === nothing || return :(k[$(parse(Int, m.captures[1]))])
        error("unknown identifier $s in a rate")
    elseif ex isa Expr
        ex.head == :call || error("unexpected expression $ex")
        f = ex.args[1]
        f in (:+, :-, :*, :/) || error("unexpected operator $f in a parameterised rate")
        args = map(e -> bind_params(e, slot, param, enz), ex.args[2:end])
        out = Expr(:call, f, args[1], args[2:min(2, end)]...)
        for x in args[3:end]
            out = Expr(:call, f, out, x)      # Python's left-nested a*b*c
        end
        return out
    else
        ex isa Real || error("unexpected literal $ex")
        return ex
    end
end

function species_read(ex, acc = Int[])
    if ex isa Expr
        if ex.head == :ref && ex.args[1] === :y
            push!(acc, ex.args[2])
        else
            foreach(a -> species_read(a, acc), ex.args)
        end
    end
    return acc
end

function load_ode(cme_names)
    cidx = Dict(s => i for (i, s) in enumerate(cme_names))
    sp = sort(csvrows(joinpath(MODEL, "ode", "species.csv"); limit = 5); by = r -> parse(Int, r[1]))
    ids = String[r[2] for r in sp]
    cme = [cidx[s] for s in ids]
    @assert cme == [parse(Int, r[3]) for r in sp] "ode cme_index is not 1-based LM order"
    init_nodes = [parse(Int, r[4]) for r in sp]
    slot = Dict(s => i for (i, s) in enumerate(ids))
    param = Dict(String(r[1]) => parse(Float64, r[2]) for r in csvrows(joinpath(MODEL, "ode", "parameters.csv"); limit = 3))
    prx = sort(csvrows(joinpath(MODEL, "params", "ode_reactions.csv"); limit = 4); by = r -> parse(Int, r[1]))
    enzyme_nodes = sort(unique(parse(Int, m.captures[1]) for r in prx for m in eachmatch(r"N__(\d+)", r[4])))
    enz = Dict(id => k for (k, id) in enumerate(enzyme_nodes))
    exprs = [bind_params(Meta.parse(r[4]), slot, param, enz) for r in prx]
    reads = [sort(unique(species_read(e))) for e in exprs]
    krows = csvrows(joinpath(MODEL, "params", "ode_constants.csv"); limit = 5)
    @assert [r[1] for r in krows] == ["K__$i" for i in 1:length(krows)]
    k_names = ["$(r[2]).$(r[3])" for r in krows]
    k_default = [parse(Float64, r[4]) for r in krows]
    st = csvrows(joinpath(MODEL, "ode", "stoichiometry.csv"))
    nr = length(prx)
    S = sparse([slot[r[1]] for r in st], [parse(Int, r[2]) for r in st], [parse(Float64, r[3]) for r in st],
               length(ids), nr)
    body = [:(V[$j] = $(exprs[j])) for j in 1:nr]
    f = @eval function (V, y, p, k)
        @inbounds begin
            $(body...)
        end
        nothing
    end
    return ODEPart(ids, cme, enzyme_nodes, init_nodes, k_names, k_default, reads, S,
                   String[r[2] for r in prx], f)
end

# ---------------------------------------------------------------------------------------
# The whole model

struct Model
    n::Int                          # states: every CME species
    names::Vector{String}
    x0::Vector{Float64}             # procid 1's initial counts
    nr::Int                         # CME reactions
    lm_type::Vector{Int}
    k_rate_col::Vector{Float64}     # reactions.csv `rate` (minute 1, MinCell_CMEODE.py)
    ra::Vector{Int}; rb::Vector{Int}    # reactants (rb = 0 for type 1)
    ptr::Vector{Int}                # reaction j's net change is entries ptr[j]:ptr[j+1]-1
    row::Vector{Int}; coef::Vector{Float64}     # after the counters are paid
    pay::Vector{Tuple{Int, Int, Vector{Tuple{Int, Float64}}}}   # the hook's payments
    pay_pos::Vector{Vector{Tuple{Int, Int}}}   # payment q, change e → (slot of (row, counter), slot of (row, pool) or 0)
    rule_nodes::Vector{Int}         # reaction j's rate from minute 2 on is node rule_nodes[j]
    dag::DAG
    ode::ODEPart
    # the hook's species
    i_cellsa::Int; i_cellsa_lip::Int; i_cellsa_prot::Int; i_cellv::Int
    lipids::Vector{Tuple{Int, Float64}}
    memprot::Vector{Int}
    # the Jacobian
    proto::SparseMatrixCSC{Float64, Int}
    pos_a::Vector{Int}              # entry e of the stoichiometry → slot of (row, ra)
    pos_b::Vector{Int}              # … slot of (row, rb), or 0
    ode_pos::Vector{Tuple{Int, Int, Int, Float64}}   # (slot, flux j, ODE species read r, S coef)
end

const KAPPA = 1.0                      # s⁻¹: the hook pays every 1 s

"""The hook's payments: (counter, pool or 0 for no clamp, the changes per unit paid)."""
function payments(names)
    idx = Dict(s => i for (i, s) in enumerate(names))
    out = Tuple{Int, Int, Vector{Tuple{Int, Float64}}}[]
    for c in ATP_HYDRO_COUNTERS
        pool, made = occursin("translat", c) ? ("M_gtp_c", "M_gdp_c") : ("M_atp_c", "M_adp_c")
        push!(out, (idx[c], idx[pool], [(idx[c], -1.0), (idx[pool], -1.0), (idx[made], 1.0), (idx["M_pi_c"], 1.0)]))
    end
    for (c, m) in vcat(NTP_COUNTERS, DNTP_COUNTERS)
        push!(out, (idx[c], idx[m], [(idx[c], -1.0), (idx[m], -1.0), (idx["M_ppi_c"], 1.0)]))
    end
    for (c, m) in NMP_COUNTERS
        push!(out, (idx[c], 0, [(idx[c], -1.0), (idx[m], 1.0)]))
    end
    for (c, ch, un) in AATRNA_COUNTERS
        push!(out, (idx[c], idx[ch], [(idx[c], -1.0), (idx[ch], -1.0), (idx[un], 1.0)]))
    end
    return out
end

"""The payment rate and its derivatives in the counter c and the pool p."""
@inline function pay_rate(c, p, clamp::Bool)
    clamp || return KAPPA * c, KAPPA * one(c), zero(c)
    d = p + 1
    return KAPPA * c * p / d, KAPPA * p / d, KAPPA * c / (d * d)
end

function load_model()
    sp = csvrows(joinpath(MODEL, "cme", "species.csv"))
    @assert all(parse(Int, r[1]) == i for (i, r) in enumerate(sp))
    names = String[r[2] for r in sp]
    x0 = Float64[parse(Float64, r[3]) for r in sp]
    n = length(names)
    idx = Dict(s => i for (i, s) in enumerate(names))
    rows = csvrows(joinpath(MODEL, "cme", "reactions.csv"))
    @assert all(parse(Int, r[1]) == j for (j, r) in enumerate(rows))
    nr = length(rows)
    lm_type = [parse(Int, r[2]) for r in rows]
    all(t -> t in (1, 2), lm_type) || error("a reaction of LM type other than 1 or 2")
    k_rate_col = [parse(Float64, r[3]) for r in rows]
    ra = zeros(Int, nr); rb = zeros(Int, nr)
    for (j, r) in enumerate(rows)
        re = parse.(Int, split(r[4], ';'))
        length(re) == lm_type[j] || error("reaction $j: $(length(re)) reactants for type $(lm_type[j])")
        ra[j] = re[1]; lm_type[j] == 2 && (rb[j] = re[2])
    end
    ptr = [1]; row = Int[]; coef = Float64[]
    for r in rows
        for pr in split(r[5], ';')
            sp_, d = split(pr, ':')
            push!(row, parse(Int, sp_)); push!(coef, parse(Float64, d))
        end
        push!(ptr, length(row) + 1)
    end
    pay = payments(names)
    rr = csvrows(joinpath(MODEL, "cme", "rate_rules.csv"))
    @assert all(parse(Int, r[1]) == j for (j, r) in enumerate(rr)) && length(rr) == nr
    rule_nodes = [parse(Int, r[2]) for r in rr]
    ode = load_ode(names)
    dag = load_dag(names, vcat(rule_nodes, ode.enzyme_nodes, ode.init_nodes))
    lipids = [(idx[s], a) for (s, a) in SA_DICT]
    memprot = [idx[s] for id in MEM_PROT_IDS for s in mem_prot_species(id)]

    # The Jacobian's pattern: CME entries, ODE entries, and the diagonal.
    I = Int[]; Jc = Int[]
    for j in 1:nr, e in ptr[j]:(ptr[j + 1] - 1)
        push!(I, row[e]); push!(Jc, ra[j])
        rb[j] > 0 && (push!(I, row[e]); push!(Jc, rb[j]))
    end
    S = ode.S
    for j in 1:size(S, 2), q in nzrange(S, j), r in ode.reads[j]
        push!(I, ode.cme[rowvals(S)[q]]); push!(Jc, ode.cme[r])
    end
    for (c, p, ch) in pay, (r, _) in ch
        push!(I, r); push!(Jc, c)
        p > 0 && (push!(I, r); push!(Jc, p))
    end
    append!(I, 1:n); append!(Jc, 1:n)
    proto = sparse(I, Jc, zeros(length(I)), n, n)
    slotof(i, c) = (rg = nzrange(proto, c); rg[findfirst(==(i), view(rowvals(proto), rg))])
    pos_a = zeros(Int, length(row)); pos_b = zeros(Int, length(row))
    for j in 1:nr, e in ptr[j]:(ptr[j + 1] - 1)
        pos_a[e] = slotof(row[e], ra[j])
        rb[j] > 0 && (pos_b[e] = slotof(row[e], rb[j]))
    end
    ode_pos = Tuple{Int, Int, Int, Float64}[]
    for j in 1:size(S, 2), q in nzrange(S, j), r in ode.reads[j]
        push!(ode_pos, (slotof(ode.cme[rowvals(S)[q]], ode.cme[r]), j, r, nonzeros(S)[q]))
    end
    pay_pos = [[(slotof(r, c), p > 0 ? slotof(r, p) : 0) for (r, _) in ch] for (c, p, ch) in pay]
    return Model(n, names, x0, nr, lm_type, k_rate_col, ra, rb, ptr, row, coef, pay, pay_pos,
                 rule_nodes, dag, ode, idx["CellSA"], idx["CellSA_Lip"], idx["CellSA_Prot"], idx["CellV"],
                 lipids, memprot, proto, pos_a, pos_b, ode_pos)
end

# ---------------------------------------------------------------------------------------
# θ: where each parameter acts

struct Block
    B::Int
    cme_slot::Vector{Int}           # CME reaction j → θ index, or 0
    k_slot::Vector{Int}             # ODE constant K__i → θ index, or 0
end
Block(m::Model) = Block(0, zeros(Int, m.nr), zeros(Int, length(m.ode.k_default)))

# ---------------------------------------------------------------------------------------
# A minute's frozen quantities, rebuilt from the state at its start

"""The volume (L) from the surface area (nm²), capped, as calcCellVolume does; and the
`CellV` value it records. Continuous: no rounding of CellV."""
function cell_volume(sa)
    r = sqrt(sa / 4 / pi) * 1e-9
    V = ((4 / 3) * pi * r^3) * 1000
    V > V_CAP ? (oftype(V, V_CAP), oftype(V, CELLV_CAP)) : (V, V * 1e19)
end

"""calcLipidToSA on the mean state, continuous (no round, no trunc)."""
function surface!(x, m::Model)
    sa = zero(eltype(x))
    for (i, a) in m.lipids
        sa += x[i] * a
    end
    sa_lip = sa * FRAC_OUT_MEM
    np = zero(eltype(x))
    for i in m.memprot
        np += x[i]
    end
    sa_prot = np * AVG_PROT_SA
    x[m.i_cellsa_lip] = sa_lip; x[m.i_cellsa_prot] = sa_prot; x[m.i_cellsa] = sa_lip + sa_prot
    return x
end

struct Minute{T, K}
    k::Vector{T}                    # CME rate constants, multipliers applied
    p::Vector{T}                    # enzyme concentrations, mM
    kode::Vector{K}                 # named ODE constants, multipliers applied
    V::T                            # volume, L
    c::T                            # NA·V/1000: counts per mM
end

"""Rebuild the minute's quantities at state `x` (modified: CellV and, unless `first`,
the surface areas are written as the hook would). `first`: the cell's initial state,
whose surface areas are taken as given."""
function rebuild!(x::AbstractVector{T}, m::Model, blk::Block, θ, kode; first = false, smooth = true) where {T}
    first || surface!(x, m)
    V, cellv = cell_volume(x[m.i_cellsa])
    x[m.i_cellv] = smooth ? cellv : T(round(ForwardDiff.value(cellv)))
    v = Vector{T}(undef, length(m.dag.op))
    eval_dag!(v, m.dag, x, smooth)
    k = Vector{T}(undef, m.nr)
    @inbounds for j in 1:m.nr
        kj = v[m.rule_nodes[j]]
        s = blk.cme_slot[j]
        k[j] = s > 0 ? kj * exp(θ[s]) : kj
    end
    p = T[v[i] for i in m.ode.enzyme_nodes]
    return Minute(k, p, kode, V, V * NA / 1000)
end

"""The named ODE constants with the block's multipliers."""
function ode_constants(m::Model, blk::Block, θ)
    K = eltype(θ) <: ForwardDiff.Dual ? eltype(θ) : Float64
    kode = Vector{K}(undef, length(m.ode.k_default))
    for i in eachindex(kode)
        s = blk.k_slot[i]
        kode[i] = s > 0 ? m.ode.k_default[i] * exp(θ[s]) : K(m.ode.k_default[i])
    end
    return kode
end

# ---------------------------------------------------------------------------------------
# The right-hand side and its Jacobian, for one minute

struct SegRHS{T, K}
    m::Model
    mn::Minute{T, K}
end

function (f::SegRHS)(du, x, _, t)
    m = f.m; mn = f.mn
    fill!(du, zero(eltype(du)))
    @inbounds for j in 1:m.nr
        a = m.rb[j] > 0 ? mn.k[j] * x[m.ra[j]] * x[m.rb[j]] : mn.k[j] * x[m.ra[j]]
        for e in m.ptr[j]:(m.ptr[j + 1] - 1)
            du[m.row[e]] += m.coef[e] * a
        end
    end
    @inbounds for (c, p, ch) in m.pay
        P, _, _ = pay_rate(x[c], p > 0 ? x[p] : zero(eltype(x)), p > 0)
        for (r, d) in ch
            du[r] += d * P
        end
    end
    o = m.ode
    y = [x[i] / mn.c for i in o.cme]
    Vf = similar(y, promote_type(eltype(y), eltype(mn.p), eltype(mn.kode)), size(o.S, 2))
    o.flux!(Vf, y, mn.p, mn.kode)
    S = o.S
    @inbounds for j in 1:size(S, 2), q in nzrange(S, j)
        du[o.cme[rowvals(S)[q]]] += mn.c * nonzeros(S)[q] * Vf[j]
    end
    return nothing
end

struct SegJac{T, K}
    f::SegRHS{T, K}
end

function (jf::SegJac)(J, x, _, t)
    m = jf.f.m; mn = jf.f.mn
    nz = nonzeros(J)
    fill!(nz, zero(eltype(nz)))
    @inbounds for j in 1:m.nr
        kj = mn.k[j]
        for e in m.ptr[j]:(m.ptr[j + 1] - 1)
            c = m.coef[e] * kj
            if m.rb[j] > 0
                nz[m.pos_a[e]] += c * x[m.rb[j]]
                nz[m.pos_b[e]] += c * x[m.ra[j]]
            else
                nz[m.pos_a[e]] += c
            end
        end
    end
    @inbounds for (q, (c, p, ch)) in enumerate(m.pay)
        _, dc, dp = pay_rate(x[c], p > 0 ? x[p] : zero(eltype(x)), p > 0)
        for (e, (r, d)) in enumerate(ch)
            sc, sp_ = m.pay_pos[q][e]
            nz[sc] += d * dc
            sp_ > 0 && (nz[sp_] += d * dp)
        end
    end
    # The ODE part: d(count_i)/dt = c·Σ_j S_ij V_j(x/c), so ∂/∂x_r = Σ_j S_ij ∂V_j/∂y_r.
    o = m.ode
    y = [x[i] / mn.c for i in o.cme]
    fl(yy) = (V = similar(yy, promote_type(eltype(yy), eltype(mn.p), eltype(mn.kode)), size(o.S, 2));
              o.flux!(V, yy, mn.p, mn.kode); V)
    dV = ForwardDiff.jacobian(fl, y)
    @inbounds for (slot, j, r, s) in m.ode_pos
        nz[slot] += s * dV[j, r]
    end
    return nothing
end

seg_tgrad!(dT, u, p, t) = (fill!(dT, zero(eltype(dT))); nothing)
