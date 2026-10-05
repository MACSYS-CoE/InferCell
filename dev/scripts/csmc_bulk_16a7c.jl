# spec/phases/16-recovery.md task 16a.7c: the path update under bulk metabolites
# (§12 2026-10-05). Per-window update rates and cost per sweep on M0, at
# N ∈ {5, 10, 20, 50}, particle Gibbs.
#
# The population is M0 replicate 16085's 50 cells, observed in bulk at its drawn
# σ, which serves as σ_b. Cells 1 to 5 are updated in turn, each conditioned on
# the other 49 through their summed latent. The rest stay at their true paths,
# which is still a valid conditional update, of these five given the others.
# After each cell's update, the sum is refreshed from the latent its sweep
# returns. Every chain starts from the true paths, an exact posterior draw at
# the true θ, so no burn-in is dropped.
#
# Stages: run <N> <scans>, and merge. Regenerate: csmc_bulk_16a7c.sh, from a
# worktree at a pushed commit.

using InferCell
using Random
using Serialization
using Statistics
using Printf

const DIR = joinpath(@__DIR__, "csmc_bulk_16a7c")
mkpath(DIR)
const REPLICATE = 16085
const NCELLS = 50
const UPDATED = 1:5
const THRESHOLD = 0.10

function setup()
    ms = m0_models()
    ds = m0_dataset(REPLICATE; ncells = NCELLS)
    T = round(Int, M0_HORIZON_S) ÷ 60
    rows = [findfirst(==(s), InferCell._block_names(ms, :ode)) for s in M0_PANEL]
    lrows = [findfirst(==(s), ds.species) for s in M0_PANEL]
    L = [ds.latent[c, lrows, :] for c in 1:NCELLS]            # panel × windows
    cells = map(UPDATED) do c
        seed = ds.seeds[c]
        Random.seed!(seed)
        d = build_m0()
        set_parameters!(d, ms, ds.truth.values)
        base = snapshot(d)
        e = restore(base)
        Random.seed!(seed)
        record_path!(e)
        tm = transcript_map(base, ms)
        tx = zeros(Int, length(tm.states), T)
        for m in 1:T
            for _ in 1:60
                handshake_step!(e)
            end
            tx[:, m] = [e.jump.u[s] for s in tm.states]
            [e.ode.u[i] * e.factor + e.rounding.remainders[i] for i in rows] == L[c][:, m] ||
                error("cell $seed: the recorded run is not the dataset's")
        end
        spec = CSMCSpec(tm, rows; sigma = ds.sigma, cap_floor = bridge_cap(maximum(tx)))
        (; seed, base, spec, data = CellData(tx, ds.bulk), ref = window_events(recorded_path(e), T))
    end
    return (; ds, L, cells, T)
end

function run_scans(N, scans)
    S = setup()
    total = sum(S.L)
    refs = [c.ref for c in S.cells]
    L = copy(S.L)
    changed = [falses(S.T, 0) for _ in UPDATED]
    secs = Float64[]
    # A plain M0 cycle, timed warm, is the cost unit.
    cyc = let b = S.cells[1].base
        run1() = (d = restore(b); Random.seed!(1); for _ in 1:round(Int, M0_HORIZON_S); handshake_step!(d); end)
        run1()
        @elapsed run1()
    end
    @printf("bulk M0, replicate %d: %d cells, σ_b = %.3g; updating cells %s at N = %d; a plain cycle takes %.2f s\n",
            REPLICATE, NCELLS, S.ds.sigma, UPDATED, N, cyc)
    rng = Xoshiro(16_076_000 + N)
    out = joinpath(DIR, "N$(N).jls")
    for s in 1:scans
        for (k, c) in enumerate(UPDATED)
            cell = S.cells[k]
            others = total .- L[c]
            bulk = bulk_window_observations(S.ds.bulk, others, NCELLS)
            t = @elapsed begin
                refs[k], ch, lat = csmc_sweep(rng, cell.base, cell.spec, cell.data, refs[k]; N, bulk)
            end
            total .+= lat .- L[c]
            L[c] = lat
            changed[k] = hcat(changed[k], ch)
            push!(secs, t)
        end
        serialize(out, (; N, cycle = cyc, changed, secs, ncells = NCELLS, updated = collect(UPDATED)))
        @printf("  scan %d: %.1f s; windows changed per cell: %s\n", s,
                sum(secs[end - length(UPDATED) + 1:end]),
                join([count(changed[k][:, end]) for k in eachindex(UPDATED)], " "))
        flush(stdout)
    end
end

function merge()
    files = sort(filter(f -> occursin(r"^N\d+\.jls$", f), readdir(DIR)); by = f -> parse(Int, f[2:end-4]))
    println("| N | Cell-sweeps | Cycles per cell-sweep | Min window rate | Windows below $(THRESHOLD) | Per-window rate, pooled over cells |")
    println("|---|---|---|---|---|---|")
    for f in files
        r = deserialize(joinpath(DIR, f))
        ch = reduce(hcat, r.changed)
        rate = vec(mean(ch; dims = 2))
        n = length(r.updated)
        # The first scan compiles; its sweeps are dropped from the cost.
        cps = mean(r.secs[(n + 1):end]) / r.cycle
        @printf("| %d | %d | %.1f | %.2f (window %d) | %d of %d | %s |\n", r.N, size(ch, 2), cps,
                minimum(rate), argmin(rate), count(<(THRESHOLD), rate), length(rate),
                join([@sprintf("%.2f", x) for x in rate], " "))
    end
end

if ARGS[1] == "run"
    run_scans(parse(Int, ARGS[2]), parse(Int, ARGS[3]))
elseif ARGS[1] == "merge"
    merge()
else
    error("stage must be run or merge")
end
