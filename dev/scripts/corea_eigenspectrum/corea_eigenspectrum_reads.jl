# KIM-181: the numbers in corea_eigenspectrum_result.md that spectrum.md does not print,
# read from the run's saved steps. Login node, a few seconds:
#
#   julia --project=. dev/scripts/corea_eigenspectrum/corea_eigenspectrum_reads.jl <run dir>
#
# <run dir> is the run's outdir (infercell-kim181-2026-11-28/run), or its archive copy.
using InferCell, Serialization, LinearAlgebra, Printf

const REACTIONS = [:PGI, :PFK, :FBA, :TPI, :GAPD, :PGK, :PGM, :ENO, :PYK, :LDH_L]
const ABOVE = 2.0
const ASYM_TOL = 1e-8

load(dir, name) = deserialize(joinpath(dir, "$name.jls"))

function main(dir)
    mf, hf, sf = load(dir, "map"), load(dir, "hessian"), load(dir, "spectrum")
    any(f -> f.smoke, (mf, hf, sf)) && error("$dir holds a smoke run")
    m, h, r = mf.result, hf.result, sf.result
    println("run: commit ", sf.commit, ", job ", sf.job, ", ", sf.ran_on.host)
    names, sd, post = r.names, r.prior_sd, r.post_sd
    A = Symmetric(-h.H)
    W = Symmetric(Diagonal(sd) * A * Diagonal(sd))

    # The counts, and how far the nearest eigenvalues sit from the threshold.
    near(v, t) = sort(v; by = e -> abs(e - t))[1:2]
    @printf("whitened −H: %d of %d above %.0f; nearest %s\n", r.n_above, length(r.ev), ABOVE,
            join((@sprintf("%.3f", e) for e in near(r.ev, ABOVE)), ", "))
    @printf("metabolic marginal: %d of %d above %.0f; nearest %s\n", r.n_met, length(r.ev_met), ABOVE,
            join((@sprintf("%.3f", e) for e in near(r.ev_met, ABOVE)), ", "))
    @printf("likelihood alone: %d of %d; metabolic %d of %d; its smallest %.4f and next %.3f\n",
            r.n_lik, length(r.ev_lik), r.n_lik_met, length(r.ev_lik_met), sort(r.ev_lik_met)[1:2]...)
    @printf("asymmetry %.3g against %.0e (margin %.1f%%)\n", r.asym, ASYM_TOL, 100 * (1 - r.asym / ASYM_TOL))

    # The full spectrum's bottom 10 against the metabolic block's marginal.
    lo = sort(r.ev)[1:length(r.ev_met)]
    @printf("bottom %d of the full whitened spectrum against the metabolic marginal's: max rel diff %.2g\n",
            length(lo), maximum(abs.(lo .- sort(r.ev_met)) ./ sort(r.ev_met)))

    # What each eigenvector of the metabolic block's whitened marginal loads on.
    met = r.met
    Am = inv(Symmetric(inv(A)[met, met]))
    E = eigen(Symmetric(Diagonal(sd[met]) * Am * Diagonal(sd[met])))
    order = sortperm(E.values; rev = true)
    rx = [Symbol(replace(String(names[i]), "log_kcatF_R_" => "")) for i in met]
    println("\nmetabolic block, whitened marginal precision: eigenvalue, then squared loadings ≥ 0.10")
    for k in order
        v = E.vectors[:, k] .^ 2
        top = sort(collect(zip(rx, v)); by = last, rev = true)
        @printf("  %8.3f  %s\n", E.values[k],
                join((@sprintf("%s %.2f", a, b) for (a, b) in top if b >= 0.10), ", "))
    end
    # The subspace the data barely reach (eigenvalue below ABOVE): its weight per reaction.
    weak = order[E.values[order] .<= ABOVE]
    w = vec(sum(E.vectors[:, weak] .^ 2; dims = 2))
    println("\nthe $(length(weak)) directions at or below $ABOVE, weight per reaction (sums to $(length(weak))):")
    println("  ", join((@sprintf("%s %.2f", a, b) for (a, b) in sort(collect(zip(rx, w)); by = last, rev = true)), ", "))

    # Each reaction's own marginal likelihood precision, 1/posterior SD² − 1/prior SD², at the
    # common width: above 1, the data alone pin it tighter than a geometric SD of 2.
    lik = (1 ./ post[met] .^ 2 .- 1 ./ sd[met] .^ 2) .* log(2)^2
    println("\nper reaction, marginal likelihood precision × (ln 2)²:")
    println("  ", join((@sprintf("%s %.3g", a, b) for (a, b) in sort(collect(zip(rx, lik)); by = last, rev = true)), ", "))

    # How close the MAP is: the Newton decrement, ½ gᵀ(−H)⁻¹g, in log-density units.
    g = m.grad
    @printf("\nNewton decrement at the MAP: %.3g (log-density units)\n", 0.5 * dot(g, A \ g))

    # How far the 21 shared coordinates moved from KIM-178's MAP, in this run's Laplace SDs.
    shared = findall(n -> !(n in [Symbol(:log_kcatF_R_, x) for x in REACTIONS]) || n in (:log_kcatF_R_ENO, :log_kcatF_R_FBA), names)
    shift = (m.fit.x[shared] .- m.x0[shared]) ./ post[shared]
    k = argmax(abs.(shift))
    @printf("shared coordinates' shift from KIM-178's MAP: max |Δ| %.2f Laplace SD (%s), median |Δ| %.2f\n",
            abs(shift[k]), names[shared][k], sort(abs.(shift))[cld(length(shift), 2)])
end

main(ARGS[1])
