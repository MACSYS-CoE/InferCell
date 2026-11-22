# KIM-96: the reduced model's uncertainty budget once every parameter is inferred
# (Linear M9, F6). No new simulation: this reads three archives and sets them side by side.
#
#   julia --project dev/scripts/corea_budget_kim96.jl <outdir>
#
# Why there is no split (Decision log 2026-11-19). KIM-42 split Q's parameter variance
# between θ_I (M0's posterior) and θ_H (15 promoters held at their priors), and 91% sat in
# θ_H. On KIM-92's joint posterior every one of the 21 coordinates is reached: the frozen
# threshold is contraction ≥ 0.5 (spec/prereg/kim-92.md §3), and the lowest is 0.919. So the
# prior-dominated set is empty and the held share is 0 by construction. What changes is how
# much parameter variance is left. This script reports it in three states:
#   - prior:     every parameter at its prior, KIM-95's prior arm;
#   - KIM-42:    M0's posterior (V7 replicate 1) with the 15 promoters at their priors,
#                width 1, `per_theta_I` (at width 1 V kept every pair, so `joint` is equal);
#   - posterior: KIM-92's joint posterior, KIM-95's posterior arm.
# Q is the same quantity in all three: `:lactate_out`, the population mean over cells at
# 6,300 s on the hybrid (`predict`'s defaults; KIM-42 passes them explicitly).
#
# Estimator (KIM-95's arms). For each draw d, μ̂_d is the mean over n_d cells and s²_d their
# variance. Var_θ = var(μ̂) − mean(s²/n): the spread across draws less the cell noise, the
# same terms as KIM-95's MC share (mean(s²/n) / var(μ̂), `corea_predict_kim95_score.jl`).
# Its SE is the jackknife over draws. KIM-42's Var(μ) and its SE are its own two-way
# jackknife, read from `budget_rows.tsv`. The three runs are independent, so a ratio's SE
# is the delta method on the log: r √((se_a/a)² + (se_b/b)²).
#
# Every input is checked against its archive's SHA256SUMS before it is read. The archive
# directories can be moved with KIM92_DIR, KIM95_DIR and KIM42_DIR (for the rerun on a copy).

using Printf, SHA, Statistics

const KIM92 = get(ENV, "KIM92_DIR", "/fred/oz022/tkimpson/infercell-kim92-archive-2026-11-17")
const KIM95 = get(ENV, "KIM95_DIR", "/fred/oz022/tkimpson/infercell-kim95-archive-2026-11-19")
const KIM42 = get(ENV, "KIM42_DIR", "/fred/oz022/tkimpson/infercell-m5-kim42-2026-11-05")
const REACH = 0.5          # spec/prereg/kim-92.md §3: "the data reach" at contraction ≥ 0.5
const WIDTH = "1.0"        # KIM-42's default prior width
const TARGET = "Q"

say(x...) = (println(stdout, x...); flush(stdout))

# ------------------------------------------------------------------ inputs, hash-checked

function checked(dir, name)
    path = joinpath(dir, name)
    sums = Dict{String, String}()
    for l in eachline(joinpath(dir, "SHA256SUMS"))
        h, f = split(l; limit = 2)
        sums[strip(lstrip(f, '*'))] = h
    end
    haskey(sums, name) || error("$name is not in $dir/SHA256SUMS")
    h = bytes2hex(open(sha256, path))
    h == sums[name] || error("$path: SHA-256 $h, but SHA256SUMS says $(sums[name])")
    return (; path, name, dir, sha256 = h)
end

function read_tsv(path)
    lines = readlines(path)
    head = split(lines[1], '\t')
    rows = [split(l, '\t') for l in lines[2:end] if !isempty(l)]
    all(r -> length(r) == length(head), rows) || error("$path: ragged rows")
    return Dict(String(h) => [String(r[i]) for r in rows] for (i, h) in enumerate(head)), length(rows)
end

# ------------------------------------------------------------------ the three states

function contraction_table(f)
    t, _ = read_tsv(f.path)
    names = t["coordinate"]
    c = parse.(Float64, t["contraction"])
    return (; names, c, reached = count(>=(REACH), c), n = length(c),
            lo = minimum(c), lo_name = names[argmin(c)])
end

function arm(fs, which)
    mu, D = read_tsv(fs["predictive_$(which)_mu.tsv"].path)
    s2, _ = read_tsv(fs["predictive_$(which)_s2.tsv"].path)
    nn, _ = read_tsv(fs["predictive_$(which)_n.tsv"].path)
    mu["draw"] == s2["draw"] == nn["draw"] || error("$which: the files' draws differ")
    m = parse.(Float64, mu[TARGET]); s = parse.(Float64, s2[TARGET]); n = parse.(Int, nn[TARGET])
    ok = [i for i in 1:D if isfinite(m[i]) && isfinite(s[i]) && n[i] >= 2]
    m, s, n = m[ok], s[ok], n[ok]
    est(idx) = var(m[idx]) - mean(s[idx] ./ n[idx])
    k = length(m)
    v = est(1:k)
    jk = [est([j for j in 1:k if j != i]) for i in 1:k]
    se = sqrt((k - 1) / k * sum(abs2, jk .- mean(jk)))
    mc = mean(s ./ n) / var(m)
    return (; label = which, draws = k, dropped = D - k, mean = mean(m),
            mean_se = std(m) / sqrt(k), var = v, se, mc, cells = extrema(n),
            sigma2 = mean(s))
end

function kim42(f)
    t, _ = read_tsv(f.path)
    pick(q) = only(i for i in eachindex(t["width"]) if t["width"][i] == WIDTH &&
                   t["target"][i] == "q" && t["law"][i] == "per_theta_I" && t["quantity"][i] == q)
    val(q) = (parse(Float64, t["estimate"][pick(q)]), parse(Float64, t["se"][pick(q)]))
    (m, mse), (v, se) = val("mean"), val("total")
    (ib, ibse), (si, sise), (sh, shse) = val("i_between"), val("share_i"), val("share_h")
    return (; label = "KIM-42", mean = m, mean_se = mse, var = v, se, i_between = ib,
            i_between_se = ibse, share_i = si, share_i_se = sise, share_h = sh, share_h_se = shse)
end

ratio(a, b) = (r = a.var / b.var; (r, r * sqrt((a.se / a.var)^2 + (b.se / b.var)^2)))

# ------------------------------------------------------------------ report

g(x) = @sprintf("%.4g", x)
e(x) = @sprintf("%.4e", x)

function report(out, files, ct, pri, k42, post)
    rows = Tuple{String, String, Float64, Float64}[]
    push!(rows, ("contraction", "reached_count", ct.reached, NaN),
          ("contraction", "coordinates", ct.n, NaN), ("contraction", "minimum", ct.lo, NaN))
    for s in (pri, k42, post)
        push!(rows, (s.label, "mean", s.mean, s.mean_se), (s.label, "var_theta", s.var, s.se),
              (s.label, "sd_theta", sqrt(max(s.var, 0)), s.se / (2sqrt(max(s.var, eps())))))
    end
    for s in (pri, post)
        push!(rows, (s.label, "mc_share", s.mc, NaN), (s.label, "draws", s.draws, NaN))
    end
    push!(rows, ("KIM-42", "held_share_I_first", k42.share_i, k42.share_i_se),
          ("KIM-42", "held_share_H_first", k42.share_h, k42.share_h_se),
          ("KIM-42", "i_between", k42.i_between, k42.i_between_se))
    rs = (("prior/posterior", ratio(pri, post)...), ("KIM-42/posterior", ratio(k42, post)...),
          ("prior/KIM-42", ratio(pri, k42)...))
    for (n, r, se) in rs
        push!(rows, ("ratio", n, r, se))
    end
    open(joinpath(out, "budget_rows.tsv"), "w") do io
        println(io, "state\tquantity\testimate\tse")
        for (a, b, x, s) in rows
            println(io, a, '\t', b, '\t', repr(Float64(x)), '\t', repr(s))
        end
    end
    open(joinpath(out, "inputs.tsv"), "w") do io
        println(io, "archive\tfile\tsha256")
        for f in files
            println(io, f.dir, '\t', f.name, '\t', f.sha256)
        end
    end
    open(joinpath(out, "budget.md"), "w") do io
        println(io, "# KIM-96: Q's parameter variance once every parameter is inferred\n")
        println(io, "Q = `lactate_out`, the population mean at 6,300 s (hybrid). Inputs and their SHA-256 are in `inputs.tsv`.\n")
        println(io, "## Contraction (KIM-92's joint posterior, `pilot_coordinates.tsv`)\n")
        println(io, "Reached (contraction ≥ $(REACH), prereg/kim-92 §3): **$(ct.reached) of $(ct.n)**. ",
                "Lowest: $(@sprintf("%.3f", ct.lo)) ($(ct.lo_name)). The prior-dominated set is ",
                ct.reached == ct.n ? "**empty**, so the held share is 0 by construction.\n" :
                "$(ct.n - ct.reached) coordinates.\n")
        println(io, "| coordinate | contraction |\n|---|---|")
        for (n, c) in zip(ct.names, ct.c)
            println(io, "| ", n, " | ", @sprintf("%.3f", c), " |")
        end
        println(io, "\n## Q's parameter variance in three states\n")
        println(io, "| state | source | mean Q ± SE | Var_θ ± SE | √Var_θ | √Var_θ / mean | MC share |")
        println(io, "|---|---|---|---|---|---|---|")
        cells(c) = c[1] == c[2] ? string(c[1]) : "$(c[1])–$(c[2])"
        src = Dict("prior" => "KIM-95 prior arm, $(pri.draws) draws × $(cells(pri.cells)) cells",
                   "KIM-42" => "KIM-42 width 1, `per_theta_I`, two-way jackknife",
                   "posterior" => "KIM-95 posterior arm, $(post.draws) draws × $(cells(post.cells)) cells")
        for s in (pri, k42, post)
            mc = hasproperty(s, :mc) ? @sprintf("%.4f", s.mc) : "—"
            println(io, "| ", s.label, " | ", src[s.label], " | ", e(s.mean), " ± ", g(s.mean_se), " | ",
                    e(s.var), " ± ", g(s.se), " | ", e(sqrt(max(s.var, 0))), " | ",
                    @sprintf("%.4f", sqrt(max(s.var, 0)) / s.mean), " | ", mc, " |")
        end
        println(io, "\nKIM-42's held share (15 promoters at their priors): ",
                @sprintf("%.3f ± %.3f", k42.share_i, k42.share_i_se), " I first, ",
                @sprintf("%.3f ± %.3f", k42.share_h, k42.share_h_se), " H first.")
        println(io, "\n## Ratios of Var_θ (delta method on independent runs)\n")
        println(io, "| ratio | estimate ± SE |\n|---|---|")
        for (n, r, se) in rs
            println(io, "| ", n, " | ", g(r), " ± ", g(se), " |")
        end
        dr = (pri.dropped, post.dropped)
        println(io, "\nDraws dropped for a non-finite μ̂ or s², or n < 2: prior $(dr[1]), posterior $(dr[2]).")
    end
end

function main(out)
    isdir(out) || mkpath(out)
    f92 = checked(KIM92, "pilot_coordinates.tsv")
    f95 = Dict(n => checked(KIM95, n) for n in ("predictive_$(a)_$(k).tsv"
                                                 for a in ("prior", "posterior") for k in ("mu", "s2", "n")))
    f42 = checked(KIM42, "budget_rows.tsv")
    files = vcat([f92], [f95[k] for k in sort(collect(keys(f95)))], [f42])
    say("inputs: all $(length(files)) match their SHA256SUMS")
    ct = contraction_table(f92)
    pri, post = arm(f95, "prior"), arm(f95, "posterior")
    k42 = kim42(f42)
    say("contraction: $(ct.reached) of $(ct.n) reached, lowest $(ct.lo) ($(ct.lo_name))")
    for s in (pri, k42, post)
        say(rpad(s.label, 10), " Var_θ = ", e(s.var), " ± ", g(s.se), ", mean ", e(s.mean))
    end
    say("MC share: prior ", pri.mc, ", posterior ", post.mc)
    say("KIM-42 held share ", k42.share_i, " ± ", k42.share_i_se)
    report(out, files, ct, pri, k42, post)
    say("wrote budget.md, budget_rows.tsv, inputs.tsv to ", out)
end

length(ARGS) == 1 || error("usage: julia --project dev/scripts/corea_budget_kim96.jl <outdir>")
main(ARGS[1])
