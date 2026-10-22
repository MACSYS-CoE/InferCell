# Issue #89, B3 close-out: the parameter bias the Jensen gap implies
# (spec/phases/approximate-method.md §12, 2026-10-23). No pass line; AV6 (B7) is
# the test. Linearised: bias ≈ (JᵀWJ + P)⁻¹ JᵀW · gap, in log θ, reported in
# F10's posterior SDs.
#
# Everything is read from the rethink's archive (`RETHINK_ARCHIVE`); nothing is
# simulated.
# - gap: per target, log x̄_SSA (its 200 cells on independent seeds) − (m + δ), and
#   − m uncorrected. δ is the rethink's δ_ref from the source's 800 cells.
# - J: the mean field's Jacobian in log θ, from the archived mean fields. Around
#   m0_low and m0_high, central differences over the ±1 posterior SD targets.
#   Around f10_ref, least squares over the 28 F10 draws (f10_d1–8, f10_p1–20). The
#   prototype mean field rounds counts inside the driver, so ForwardDiff (AV3)
#   does not apply to it.
# - W = (σ_b² I + S)⁻¹, S = Cov(log x̄_50) at the source, by the delta method from
#   its 800 cells. σ_b at KG's 0.10 and at F10's posterior median.
# - P: the targets' prior precisions in log θ (all LogNormal).
# - SE: the gap's sampling error (the target's 200 cells, and δ's 800) propagated
#   linearly.
# Information only: the linear response with F10's posterior covariance in place
# of (JᵀWJ + P)⁻¹, which carries the transcripts' information, at the F10
# targets; and δ taken at m0_v6 (Stage 1's 200 cells), re-read on the F10
# targets' independent seeds.
#
# Usage: julia --project dev/scripts/jensen_bias_89.jl <outdir>
# Run:   sbatch dev/scripts/jensen_bias_89.slurm <outdir outside the worktrees>

using InferCell
using Distributions: LogNormal, params
using Serialization
using Statistics
using LinearAlgebra
using Printf
using Dates

const COMMIT = strip(read(`git rev-parse --short HEAD`, String))
const ARCHIVE = get(ENV, "RETHINK_ARCHIVE", "/fred/oz022/tkimpson/infercell-89-rethink-archive-2026-10-22")
const F10_DRAWS = joinpath(@__DIR__, "f10_16a10_draws.tsv")
const LOGNAMES = [:S_JCVISYN3A_0607, :S_JCVISYN3A_0779, :krnadeg, :kcatF_R_ENO, :kcatF_R_FBA]
const SHORTN = ["S_0607", "S_0779", "krnadeg", "kcatF_ENO", "kcatF_FBA"]
const C = 50          # cells in a bulk sample
const NKG = 200
const KG_SIGMA = 0.10

slurm_job() = get(ENV, "SLURM_JOB_ID", "none")

# F10's draws on the log scale: the targets' columns and σ_b's.
function f10()
    lines = readlines(F10_DRAWS)
    head = split(lines[1], '\t')
    col(n) = findfirst(==(n), head)
    X = [parse(Float64, split(l, '\t')[col(String(n))]) for l in lines[2:end], n in LOGNAMES]
    lσ = [parse(Float64, split(l, '\t')[col("sigma_b")]) for l in lines[2:end]]
    return X, vec(std(X; dims = 1)), cov(X), exp(median(lσ))
end

function load_cells(dir)
    files = sort(filter(f -> occursin(r"^task_\d+\.jls$", f), readdir(dir)))
    isempty(files) && error("no cells in $dir")
    cells = sort(vcat((deserialize(joinpath(dir, f)).cells for f in files)...); by = c -> c.seed)
    allunique(c.seed for c in cells) || error("repeated seeds in $dir")
    return cells
end
ens(pt, s0) = load_cells(joinpath(ARCHIVE, "ens_$(pt)_full_s$(s0)"))
# A target's own seeds: the one full ensemble that is not at the shared 50000.
function own_seed0(pt)
    ds = filter(d -> startswith(d, "ens_$(pt)_full_s") && !endswith(d, "_s50000"), readdir(ARCHIVE))
    length(ds) == 1 || error("expected one independent ensemble for $pt, found $ds")
    return parse(Int, last(split(ds[1], "_s")))
end

mfr(pt) = deserialize(joinpath(ARCHIVE, "meanfield", "$pt.jls"))
mflog(pt) = vec(log.(max.(mfr(pt).counts[:, 2:end], 1.0)))
logs_of(vals) = [log(Dict(vals)[n]) for n in LOGNAMES]

# Cells × (pool, save), saves after t = 0, pools fastest (as vec of panel × saves).
cellmat(cs) = reduce(vcat, [vec(c.counts[:, 2:end])' for c in cs])
logbulk(cs) = log.(max.(vec(mean(cellmat(cs); dims = 1)), 1.0))
# Cov(log x̄_n) by the delta method from per-cell counts.
function logmean_cov(cs, n)
    X = cellmat(cs)
    x̄ = vec(mean(X; dims = 1))
    return Symmetric(cov(X) ./ (x̄ * x̄') ./ n)
end

function prior_sds()
    ps = Dict(p.name => p.prior for m in m0_models() for p in parameters(m))
    return map(LOGNAMES) do n
        d = ps[n]
        d isa LogNormal || error("$n has a $(typeof(d)) prior")
        params(d)[2]
    end
end

# J around m0_low / m0_high: central differences over the ±1 SD targets.
function jac_central(tag, sd)
    m0 = mflog("m0_$tag")
    J = zeros(length(m0), length(LOGNAMES))
    curv = zeros(length(LOGNAMES))
    for (j, n) in enumerate(LOGNAMES)
        mp, mm = mflog("$(tag)_$(n)_+"), mflog("$(tag)_$(n)_-")
        lp, lm = logs_of(mfr("$(tag)_$(n)_+").values)[j], logs_of(mfr("$(tag)_$(n)_-").values)[j]
        J[:, j] = (mp .- mm) ./ (lp - lm)
        curv[j] = norm(mp .+ mm .- 2m0) / max(norm(mp .- mm), eps())
    end
    return J, curv
end

# J around f10_ref: least squares over the 28 F10 draws' mean fields.
function jac_regress()
    m0, l0 = mflog("f10_ref"), logs_of(mfr("f10_ref").values)
    pts = vcat(["f10_d$k" for k in 1:8], ["f10_p$k" for k in 1:20])
    Θ = reduce(vcat, [(logs_of(mfr(p).values) .- l0)' for p in pts])
    Y = reduce(vcat, [(mflog(p) .- m0)' for p in pts])
    Jt = Θ \ Y
    r2 = 1 - sum(abs2, Y .- Θ * Jt) / sum(abs2, Y)
    return Matrix(Jt'), r2, length(pts)
end

kg_total(M, L, σ) = sum(abs2, M .- L) / (2σ^2)

function main()
    out = ARGS[1]
    mkpath(out)
    io = IOBuffer()
    p(args...) = (println(io, args...); println(args...))
    X, sd, ΣF, σF = f10()
    psd = prior_sds()
    P = Diagonal(1 ./ psd .^ 2)
    σs = (KG_SIGMA, σF)
    p("# Issue #89 B3 close-out: the implied parameter bias")
    p()
    p("At $COMMIT (job $(slurm_job())) over `$ARCHIVE`, $(now()). No pass line (§12 2026-10-23).")
    p()
    p("- Bias in F10 posterior SDs (", join([@sprintf("%s %.3g", s, v) for (s, v) in zip(SHORTN, sd)], ", "), ").")
    p("- Prior SDs in log θ: ", join([@sprintf("%s %.3g", s, v) for (s, v) in zip(SHORTN, psd)], ", "), ".")
    p(@sprintf("- σ_b: KG's %.2f, and F10's posterior median %.4f.", KG_SIGMA, σF))
    p("- Each cell is `bias ± SE`, the SE being the gap's sampling error propagated linearly.")
    p()

    Jf, r2, nf = jac_regress()
    Jl, cl = jac_central("low", sd)
    Jh, ch = jac_central("high", sd)
    p("## The Jacobian")
    p()
    p(@sprintf("- f10_ref: least squares over %d F10 draws; R² = %.4f.", nf, r2))
    p("- m0_low, central differences over ±1 SD; curvature ‖m₊ + m₋ − 2m₀‖ / ‖m₊ − m₋‖: ",
      join([@sprintf("%s %.3f", s, v) for (s, v) in zip(SHORTN, cl)], ", "), ".")
    p("- m0_high, the same: ", join([@sprintf("%s %.3f", s, v) for (s, v) in zip(SHORTN, ch)], ", "), ".")
    p()

    groups = (("f10_ref", "f10_d", Jf), ("m0_low", "low_", Jl), ("m0_high", "high_", Jh))
    targets(prefix) = sort(filter(d -> d != "f10_ref" && startswith(d, prefix) &&
                                       !startswith(d, "f10_p"),
                                  [replace(f, ".jls" => "") for f in readdir(joinpath(ARCHIVE, "meanfield"))]))
    rows = Any[]
    for σ in σs
        p(@sprintf("## σ_b = %.4f", σ))
        p()
        for (src, prefix, J) in groups
            rc = ens(src, 40_000)
            length(rc) == 800 || error("$src has $(length(rc)) cells, not 800")
            δ = logbulk(rc) .- mflog(src)
            S = logmean_cov(rc, C)
            W = inv(Symmetric(σ^2 * I + S))
            H = Symmetric(J' * W * J + P)
            A = H \ (J' * W)
            lin_sd = sqrt.(diag(inv(H)))
            Σδ = logmean_cov(rc, length(rc))
            p("### Source $src (800 cells)")
            p()
            p("Linearised posterior SD (bulk + prior only, no transcripts) ÷ F10's: ",
              join([@sprintf("%s %.2f", s, v) for (s, v) in zip(SHORTN, lin_sd ./ sd)], ", "), ".")
            p()
            p("| Target | arm | ", join(SHORTN, " | "), " |")
            p("|---|---|", repeat("---|", length(SHORTN)))
            for t in targets(prefix)
                cs = ens(t, own_seed0(t))[1:NKG]
                L, M = logbulk(cs), mflog(t)
                St = logmean_cov(cs, length(cs))
                for (arm, gap, Σg) in (("uncorrected", L .- M, St), ("δ_ref", L .- (M .+ δ), St + Σδ))
                    b = A * gap ./ sd
                    se = sqrt.(max.(diag(A * Σg * A'), 0.0)) ./ sd
                    push!(rows, (; σ, src, t, arm, b, se))
                    p("| $t | $arm | ", join([@sprintf("%+.2f ± %.2f", x, e) for (x, e) in zip(b, se)], " | "), " |")
                end
            end
            p()
            for arm in ("uncorrected", "δ_ref")
                R = [r for r in rows if r.σ == σ && r.src == src && r.arm == arm]
                B = reduce(hcat, [r.b for r in R]); E = reduce(hcat, [r.se for r in R])
                p("- **$arm:** max |bias| ", join([@sprintf("%s %.2f", s, v) for (s, v) in zip(SHORTN, vec(maximum(abs.(B); dims = 2)))], ", "),
                  "; RMS ", join([@sprintf("%.2f", v) for v in vec(sqrt.(mean(B .^ 2; dims = 2)))], ", "),
                  "; |bias| > 2 SE at ", join([string(count(abs.(B[j, :]) .> 2E[j, :])) for j in axes(B, 1)], ", "),
                  " of $(size(B, 2)).")
            end
            p()
            if src == "f10_ref"
                p("**Information only: F10's posterior covariance in place of (JᵀWJ + P)⁻¹** (it carries the transcripts):")
                p()
                p("| Target | arm | ", join(SHORTN, " | "), " |")
                p("|---|---|", repeat("---|", length(SHORTN)))
                AF = ΣF * J' * W
                for t in targets(prefix)
                    cs = ens(t, own_seed0(t))[1:NKG]
                    L, M = logbulk(cs), mflog(t)
                    for (arm, gap) in (("uncorrected", L .- M), ("δ_ref", L .- (M .+ δ)))
                        p("| $t | $arm | ", join([@sprintf("%+.2f", x) for x in AF * gap ./ sd], " | "), " |")
                    end
                end
                p()
            end
        end
    end

    # δ at m0_v6, re-read on the F10 targets' independent seeds.
    p("## Information only: δ taken at m0_v6, on the F10 targets' independent seeds")
    p()
    v6 = ens("m0_v6", 40_000)
    δv6 = logbulk(v6) .- mflog("m0_v6")
    δref = logbulk(ens("f10_ref", 40_000)) .- mflog("f10_ref")
    rc = ens("f10_ref", 40_000)
    Σv6, Σref = logmean_cov(v6, length(v6)), logmean_cov(rc, length(rc))
    lv6 = logs_of(mfr("m0_v6").values)
    μF = vec(mean(X; dims = 1))
    maha(l) = sqrt(dot(l .- μF, ΣF \ (l .- μF)))
    dmaha = [maha(X[i, :]) for i in axes(X, 1)]
    gtp = findfirst(==(:M_gtp_c), M0_PANEL)
    gtpmean(d) = mean(reshape(d, length(M0_PANEL), :)[gtp, :])
    p("- m0_v6 in F10's posterior: z = ", join([@sprintf("%s %+.2f", s, v) for (s, v) in zip(SHORTN, (lv6 .- μF) ./ sd)], ", "),
      @sprintf("; Mahalanobis %.2f, against the draws' median %.2f and 95%% point %.2f.",
               maha(lv6), quantile(dmaha, 0.5), quantile(dmaha, 0.95)))
    # δ_v6's cells share their seeds with 200 of δ_ref's 800, so the noise in
    # δ_v6 − δ_ref is paired: per-seed influence terms of the two log means.
    Xv, Xr = cellmat(v6), cellmat(rc)
    x̄v, x̄r = vec(mean(Xv; dims = 1)), vec(mean(Xr; dims = 1))
    iv = (Xv .- x̄v') ./ x̄v' ./ size(Xv, 1)
    ir = (Xr .- x̄r') ./ x̄r' ./ size(Xr, 1)
    rseed = Dict(c.seed => k for (k, c) in enumerate(rc))
    shared = [rseed[c.seed] for c in v6 if haskey(rseed, c.seed)]
    length(shared) == length(v6) || error("expected every m0_v6 seed among f10_ref's")
    U = iv .- ir[shared, :]
    rest = setdiff(axes(Xr, 1), shared)
    Σd = U' * U + ir[rest, :]' * ir[rest, :]
    d = δv6 .- δref
    z2 = reshape(d .^ 2 ./ diag(Σd), length(M0_PANEL), :)
    pz = vec(mean(z2; dims = 2))
    p("- δ from m0_v6's $(length(v6)) Stage 1 cells (seeds 40001…, the same seeds as the first $(length(v6)) of δ_ref's $(length(rc))): ",
      @sprintf("‖δ_v6 − δ_ref‖ = %.3f over all pools and saves, against %.3f expected from paired sampling noise (√tr; %.3f if the two were independent). ",
               norm(d), sqrt(tr(Σd)), sqrt(tr(Σv6 + Σref))),
      @sprintf("GTP's mean over saves: δ_v6 %+.3f, δ_ref %+.3f.", gtpmean(δv6), gtpmean(δref)))
    p(@sprintf("- Per pool, the mean over saves of (δ_v6 − δ_ref)² ÷ its paired variance (about 1 for no difference; saves are correlated, so each pool has few degrees of freedom); overall %.2f: ", mean(pz)),
      join([@sprintf("%s %.2f", replace(String(s), "M_" => "", "_c" => ""), v) for (s, v) in zip(M0_PANEL, pz)], ", "), ".")
    p("- KG's noise is the target's split-half statistic, as in `transfer.md`. With δ from 200 cells, a right δ's expected total is about the noise; with δ from 800, about 0.6 of it.")
    p()
    W = inv(Symmetric(KG_SIGMA^2 * I + logmean_cov(rc, C)))
    A = Symmetric(Jf' * W * Jf + P) \ (Jf' * W)
    kg_noise(cs) = (h = length(cs) ÷ 2; sum(abs2, logbulk(cs[1:h]) .- logbulk(cs[(h + 1):(2h)])) / 2 / (2KG_SIGMA^2))
    p("| Target | KG total: uncorrected / δ_ref / δ_v6 (noise) | bias with δ_v6 (σ_b = 0.10): ", join(SHORTN, ", "), " |")
    p("|---|---|---|")
    for t in targets("f10_d")
        cs = ens(t, own_seed0(t))[1:NKG]
        L, M = logbulk(cs), mflog(t)
        tot = [kg_total(M .+ o, L, KG_SIGMA) for o in (0.0, δref, δv6)]
        se = sqrt.(max.(diag(A * (logmean_cov(cs, length(cs)) + Σv6) * A'), 0.0)) ./ sd
        p("| $t | ", @sprintf("%.1f / %.1f / %.1f (%.1f)", tot..., kg_noise(cs)), " | ",
          join([@sprintf("%+.2f ± %.2f", x, e) for (x, e) in zip(A * (L .- (M .+ δv6)) ./ sd, se)], ", "), " |")
    end
    p()
    write(joinpath(out, "bias.md"), take!(io))
    serialize(joinpath(out, "bias.jls"), (; commit = COMMIT, job = slurm_job(), rows, sd, psd, σF, Jf, Jl, Jh, r2, cl, ch))
    println("Wrote ", joinpath(out, "bias.md"))
end

main()
