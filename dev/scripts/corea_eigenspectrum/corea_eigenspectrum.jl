# KIM-181: the Hessian's eigen-spectrum on the reduced model with a wider metabolic
# block (Linear M13). The fit is KIM-178's fitted-likelihood MAP on KIM-92's dataset
# (prereg/kim-178 §2 steps 5 and 6), with the free set widened from ENO and FBA to the
# forward constants of all ten glycolytic reactions. It says how many directions the data
# reach when the metabolic block is wider. Tom's calls of 2026-11-27: the ten kcatF, and
# the MAP refitted in the wider space.
#
# Usage: julia --project dev/scripts/corea_eigenspectrum/corea_eigenspectrum.jl <out>
# Submit it through corea_eigenspectrum.sh, from a worktree at a pushed commit.
#
# One job, in three steps. A step whose .jls exists is skipped, so a killed job resumes:
#   1. map:      the fitted-likelihood MAP over 29 coordinates, warm-started from
#                KIM-178's 21-coordinate MAP (job 21727400) with the 8 new constants at
#                their prior medians, which are their nominal values and the truth   -> map.jls
#   2. hessian:  −H at that MAP, by KIM-92's `hessian_fd`; then three entries checked
#                against second differences of the log density itself                -> hessian.jls
#   3. spectrum: −H whitened by the prior's SDs, its eigenvalues, and the same for the
#                metabolic block's marginal precision                                 -> spectrum.jls, spectrum.md
#
# The dataset. KIM-92's was drawn over COREA_TARGETS, with the 8 other forward constants
# at their nominal values. `mixed_setup` refuses a dataset whose recorded targets differ
# from the fit's, so the record is rebuilt with the wider list: the observations are
# untouched (checked), and the 8 constants' truth is their nominal value. Freeing a
# forward constant derives its reverse through Haldane; at the nominal values every
# derived reverse equals the nominal one to 2.2e-16 relative (checked on 2026-11-27, and
# again below), so the wider model scores the same data at the same truth.
#
# δ and S are KIM-92's (its prep.jls), the reference ensemble at its θ_ref. They are the
# fitted likelihood's, and do not depend on which constants are free.
#
# SMOKE=1 is a plumbing pass, never a result: 2 MAP iterations, and −H's columns only for
# the checked entries' coordinates, the rest of it the prior's precision.

include(joinpath(@__DIR__, "..", "pilot_common.jl"))
using InferCell: derived_ode_values
using ForwardDiff
using Printf
using LinearAlgebra

const KIM92_ARCHIVE = "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17"
const KIM178_MAP = "/fred/oz022/tkimpson/infercell-kim178-2026-11-27/kim92"
const KIM92_SEED = 17_200_001
const MAP_ITERS = SMOKE ? 2 : 200
const FD_STEP = 1e-3          # × the coordinate's prior SD, as KIM-92's metric and KIM-178
const FD2_STEP = 1e-2         # × prior SD, for the second differences of the density
const FD2_TOL = 1e-2          # relative agreement of a checked entry
const ASYM_TOL = 1e-8         # the issue's symmetry check
const ABOVE = 2.0             # whitened −H eigenvalue at which the data's curvature equals the prior's

const REACTIONS = [:PGI, :PFK, :FBA, :TPI, :GAPD, :PGK, :PGM, :ENO, :PYK, :LDH_L]
const FORWARDS = [Symbol(:kcatF_R_, r) for r in REACTIONS]
const REVERSES = [Symbol(:kcatR_R_, r) for r in REACTIONS]
const NEW = setdiff(FORWARDS, COREA_TARGETS)                 # the 8 newly freed
const TARGETS = vcat(COREA_TARGETS, NEW)                     # 28 names, 29 coordinates with σ_b

# d11_models with all ten reactions' forward and reverse constants free. Only the forward
# constants are targets; each reverse is derived through Haldane.
function wide_models()
    ms = corea_models()
    ms[1] = CentralGlycolysis(enzymes = :translated, free = vcat(FORWARDS, REVERSES))
    return ms
end

setup(ds, p) = with_reference(mixed_setup(ds; models = wide_models(), targets = TARGETS), p.delta, p.S)
has(dir, name) = isfile(joinpath(dir, "$name.jls"))

function inputs()
    ds0 = deserialize(joinpath(KIM92_ARCHIVE, "dataset.jls")).result
    p = deserialize(joinpath(KIM92_ARCHIVE, "prep.jls")).result
    ds0.replicate_seed == KIM92_SEED || error("not KIM-92's dataset")
    collect(Symbol, ds0.targets) == COREA_TARGETS || error("KIM-92's dataset is not over COREA_TARGETS")
    ds = merge(ds0, (; targets = TARGETS))
    ds.bulk === ds0.bulk && ds.transcripts === ds0.transcripts && ds.times === ds0.times ||
        error("rebuilding the record changed its observations")
    ms = wide_models()
    nom = nominal_parameter_values(ms)
    der = Dict(derived_ode_values(ms, [f => nom[f] for f in FORWARDS]))
    hal = maximum(abs(der[r] / nom[r] - 1) for r in REVERSES)
    hal < 1e-12 || error(@sprintf("a derived reverse differs from its nominal value by %.3g", hal))
    truth = merge(Dict(ds0.truth.values), Dict(n => nom[n] for n in NEW))
    return (; ds, p, truth, haldane = hal)
end

truth_x(s, truth, ds) = vcat([log(truth[n]) for n in s.names], log(ds.sigma))

# KIM-92's `hessian_fd`, as KIM-178 copied it, with `cols` to compute only some columns
# (a smoke). Columns not computed are left zero.
function hessian_fd(grad, x, h; cols = eachindex(x))
    P = length(x)
    H = zeros(P, P)
    for (k, j) in enumerate(cols)
        e = zeros(P); e[j] = h[j]
        H[:, j] = (grad(x .+ e) .- grad(x .- e)) ./ (2h[j])
        say(@sprintf("hessian: column %d (%d of %d)", j, k, length(cols)))
    end
    return H
end

# A second difference of the log density alone: no gradient in it, so it checks the AD
# gradient and the column differencing together.
function d2(f, x, i, j, hi, hj)
    ei = zeros(length(x)); ei[i] = hi
    ej = zeros(length(x)); ej[j] = hj
    if i == j
        return (f(x .+ ei) - 2f(x) + f(x .- ei)) / hi^2
    end
    return (f(x .+ ei .+ ej) - f(x .+ ei .- ej) - f(x .- ei .+ ej) + f(x .- ei .- ej)) / (4hi * hj)
end

# ---------------------------------------------------------------- the three steps

function step_map(dir, inp)
    has(dir, "map") && return
    s = setup(inp.ds, inp.p)
    names = mixed_parameter_names(s)
    m178 = deserialize(joinpath(KIM178_MAP, "map.jls")).result
    n178 = Symbol.(deserialize(joinpath(KIM178_MAP, "kim92.jls")).result.names)
    x178 = Dict(Symbol(n) => v for (n, v) in zip(n178, m178.fit.x))
    mu = Dict(zip(names, mixed_prior_median(s)))
    x0 = [get(x178, n, mu[n]) for n in names]
    count(n -> haskey(x178, n), names) == length(n178) || error("KIM-178's coordinates are not all in the wider set")
    f(x) = mixed_logdensity(s, x)
    lp0 = f(x0)
    g0 = ForwardDiff.gradient(f, x0)
    t = @elapsed fit = reference_theta(s; x0, maxiters = MAP_ITERS)
    lp = f(fit.x)
    g = ForwardDiff.gradient(f, fit.x)
    sd = mixed_prior_sd(s)
    save(dir, "map", (; fit, names, x0, lp0, lp, grad0 = g0, grad = g, lp_178 = m178.fit.logdensity,
                      haldane = inp.haldane, seconds = t); stage = "map")
    say(@sprintf("map: by %s in %.0f s; log density %.4f at the MAP, %.4f at the warm start (Δ %+.4f); |∇|·sd %.3g at the start, %.3g at the MAP", fit.method, t, lp, lp0, lp - lp0,
                 norm(g0 .* sd), norm(g .* sd)))
end

function step_hessian(dir, inp)
    has(dir, "hessian") && return
    s = setup(inp.ds, inp.p)
    m = load(dir, "map")
    x, names = m.fit.x, m.names
    sd = mixed_prior_sd(s)
    f(z) = mixed_logdensity(s, z)
    grad(z) = ForwardDiff.gradient(f, z)
    idx(n) = findfirst(==(n), names)
    # Three entries: a promoter's diagonal, a new forward constant's diagonal, and a cross
    # term between two forward constants, one of them new.
    checks = [(names[1], names[1]), (:log_kcatF_R_PYK, :log_kcatF_R_PYK),
              (:log_kcatF_R_PFK, :log_kcatF_R_ENO)]
    ij = [(idx(a), idx(b)) for (a, b) in checks]
    any(t -> any(isnothing, t), ij) && error("a checked coordinate is not in the setup")
    cols = SMOKE ? sort!(unique(vcat(first.(ij), last.(ij)))) : eachindex(x)
    t = @elapsed H = hessian_fd(grad, x, FD_STEP .* sd; cols)
    if SMOKE                                  # plumbing only: the prior's precision elsewhere
        for j in setdiff(eachindex(x), cols)
            H[j, j] = -1 / sd[j]^2
        end
        asym = NaN
    else
        asym = maximum(abs, H .- H') / maximum(abs, H)
    end
    Hs = (H .+ H') ./ 2
    fd = map(ij) do (i, j)
        v = d2(f, x, i, j, FD2_STEP * sd[i], FD2_STEP * sd[j])
        (; i, j, a = names[i], b = names[j], H = Hs[i, j], fd = v, rel = abs(Hs[i, j] - v) / abs(v))
    end
    save(dir, "hessian", (; H = Hs, asym, fd, fd_step = FD_STEP, fd2_step = FD2_STEP, seconds = t);
         stage = "hessian")
    say(@sprintf("hessian: %.0f s; asymmetry %.3g", t, asym))
    for c in fd
        say(@sprintf("  check H[%s, %s] = %.6g against %.6g by second differences (rel %.2g)", c.a, c.b, c.H, c.fd, c.rel))
    end
end

# The marginal precision of a block: the inverse of the covariance's block.
marginal_precision(A, b) = inv(Symmetric(inv(Symmetric(A))[b, b]))

function whitened(A, sd)
    D = Diagonal(sd)
    return Symmetric(D * A * D)
end

function step_spectrum(dir, inp)
    s = setup(inp.ds, inp.p)
    m, h = load(dir, "map"), load(dir, "hessian")
    names, x = m.names, m.fit.x
    sd = mixed_prior_sd(s)
    A = Symmetric(-h.H)                                       # the posterior precision at the MAP
    ev_raw = eigvals(A)
    pd = ev_raw[1] > 0
    W = whitened(A, sd)                                       # the prior's precision becomes I
    ev = eigvals(W)
    n_above = count(>(ABOVE), ev)
    met = findall(n -> n in [Symbol(:log_, f) for f in FORWARDS], names)
    ev_met = pd ? eigvals(whitened(marginal_precision(A, met), sd[met])) : fill(NaN, length(met))
    n_met = count(>(ABOVE), ev_met)
    post_sd = pd ? sqrt.(diag(inv(A))) : fill(NaN, length(x))
    xt = truth_x(s, inp.truth, inp.ds)
    lb, ub = mixed_box(s)
    edge = min.(x .- lb, ub .- x) ./ sd
    # KIM-178's run on the narrow block (job 21727400), for the comparison: ENO and FBA.
    l178 = deserialize(joinpath(KIM178_MAP, "laplace.jls")).result
    n178 = Symbol.(deserialize(joinpath(KIM178_MAP, "kim92.jls")).result.names)
    s178 = mixed_setup(merge(inp.ds, (; targets = COREA_TARGETS)); models = d11_models(),
                       targets = COREA_TARGETS)
    sd178 = mixed_prior_sd(s178)
    met178 = findall(n -> n in (:log_kcatF_R_ENO, :log_kcatF_R_FBA), n178)
    A178 = Symmetric(-l178.H)
    ev178 = eigvals(whitened(A178, sd178))
    ev178_met = eigvals(whitened(marginal_precision(A178, met178), sd178[met178]))
    r = (; names, x, truth = xt, prior_sd = sd, post_sd, contraction = post_sd ./ sd, edge,
         ev_raw, ev, n_above, met, ev_met, n_met, positive_definite = pd,
         ev178, n_above178 = count(>(ABOVE), ev178), ev178_met, n_met178 = count(>(ABOVE), ev178_met),
         asym = h.asym, fd = h.fd, lp = m.lp, lp0 = m.lp0, method = m.fit.method,
         seconds = (map = m.seconds, hessian = h.seconds))
    save(dir, "spectrum", r; stage = "spectrum")
    report(dir, r)
end

function report(dir, r)
    io = IOBuffer()
    pr(xs...) = println(io, xs...)
    pr("# KIM-181: the Hessian's eigen-spectrum with the ten glycolytic forward constants free\n")
    pr(@sprintf("Commit `%s`, job %s%s. KIM-92's dataset (seed %d), the fitted likelihood (δ, S from KIM-92's prep).\n",
                COMMIT, JOB, SMOKE ? ", **SMOKE (not a result)**" : "", KIM92_SEED))
    pr(@sprintf("- **MAP:** by %s in %.0f s; log density %.4f, against %.4f at the warm start (Δ %+.4f): %s.",
                r.method, r.seconds.map, r.lp, r.lp0, r.lp - r.lp0, r.lp >= r.lp0 ? "holds" : "FAILS"))
    pr(@sprintf("- **−H:** %.0f s; asymmetry before symmetrising %.3g (check ≤ %.0e: %s); positive definite: %s; nearest the box %.2f prior SD.", r.seconds.hessian, r.asym, ASYM_TOL,
                r.asym <= ASYM_TOL ? "holds" : "FAILS", r.positive_definite, minimum(r.edge)))
    for c in r.fd
        pr(@sprintf("- **Check** H[%s, %s] = %.6g, by second differences of the density %.6g (rel %.2g, check < %.0e: %s).",
                    c.a, c.b, c.H, c.fd, c.rel, FD2_TOL, c.rel < FD2_TOL ? "holds" : "FAILS"))
    end
    pr(@sprintf("- **Directions above the prior's curvature** (whitened −H eigenvalue > %.0f, where the data's curvature exceeds the prior's): **%d of %d** with the wide block; %d of %d with KIM-178's narrow block (ENO, FBA).",
                ABOVE, r.n_above, length(r.ev), r.n_above178, length(r.ev178)))
    pr(@sprintf("- **The metabolic block alone** (marginal precision of the %d forward constants, whitened): **%d of %d** above; the narrow block's ENO and FBA: %d of 2.\n", length(r.met), r.n_met, length(r.met), r.n_met178))
    pr("## Whitened spectrum of −H (prior precision = 1)\n")
    pr("| rank | wide block (29) | narrow block (21) |")
    pr("|---|---|---|")
    e1, e2 = sort(r.ev; rev = true), sort(r.ev178; rev = true)
    for k in eachindex(e1)
        pr(@sprintf("| %d | %.4g | %s |", k, e1[k], k <= length(e2) ? @sprintf("%.4g", e2[k]) : ""))
    end
    pr("\n## The metabolic block's marginal precision, whitened\n")
    pr("| rank | eigenvalue |")
    pr("|---|---|")
    for (k, v) in enumerate(sort(r.ev_met; rev = true))
        pr(@sprintf("| %d | %.4g |", k, v))
    end
    pr("\n## Per coordinate\n")
    pr("| coordinate | truth | MAP | prior SD | Laplace SD | Laplace / prior |")
    pr("|---|---|---|---|---|---|")
    for i in eachindex(r.names)
        pr(@sprintf("| %s | %.4f | %.4f | %.4f | %.4f | %.3f |", r.names[i], r.truth[i], r.x[i],
                    r.prior_sd[i], r.post_sd[i], r.contraction[i]))
    end
    write(joinpath(dir, "spectrum.md"), String(take!(io)))
    print(read(joinpath(dir, "spectrum.md"), String))
end

function main(args)
    length(args) == 1 || error("usage: corea_eigenspectrum.jl <out>")
    dir = abspath(args[1])
    mkpath(dir)
    inp = inputs()
    say(@sprintf("inputs: %d targets (%d new); Haldane at nominal %.2g", length(TARGETS), length(NEW), inp.haldane))
    step_map(dir, inp)
    step_hessian(dir, inp)
    step_spectrum(dir, inp)
end

main(ARGS)
