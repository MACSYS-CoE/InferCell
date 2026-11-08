# Probabilistic trajectories Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `predict` returns time series with uncertainty bands, and a new paper figure shows them on M0, on the reduced model, and with 62 more parameters released to their priors.

**Architecture:** `_ForwardSetup` carries the species and time grid to record, and `_cell` reads them from the `ValidationRun` it already holds, so no simulation changes. A new `src/predict_trajectory.jl` turns recorded cells into bands and a variance split over time. One driver replays the published runs on their seeds, with exact endpoint checks, and runs the new release arm.

**Tech Stack:** Julia 1.10.5, InferCell, Slurm (OzSTAR), Python/matplotlib, LaTeX (TinyTeX).

**Spec:** `docs/superpowers/specs/2026-11-07-probabilistic-trajectories-design.md`

## Global Constraints

- With `trajectory = nothing`, results stay bit for bit identical: `test/test_predict_reference.jl` must pass unchanged, and its fixture is not regenerated.
- Fixed before any run: replicate 1, exact posterior; width 1; the released set is every currently fixed `kcatF_*` and `km_*` of `CentralGlycolysis` and `NucleotideRecycling`, 62 parameters; species `LAC_EXPORT_METER`, `M_lac__L_e`, `M_gtp_c`, `M_pyr_c`.
- Seeds: panel B `Xoshiro(19_542_000)`. Panel A blocks `19_540_000 + 10(r−1) + (b−1)`, truth `19_541_500 + r`, with r = 1. Panel C `19_543_000`, its pre-check `19_543_001`, the smoke pass `19_543_900`.
- A replay mismatch stops the work. Panel C goes ahead whatever its acceptance. Below 10% of pairs kept, the figure drops the filled conditioned band and shows the unconditioned one, labelled.
- Run through Slurm from a detached worktree at a pushed commit, on the nodes `launcher_guard.sh` gives, with `--gcthreads=N,1` on many threads. Parse and lower scripts on the login node first, and wrap top-level script loops in functions.
- Never push to `main`; feature branches and PRs. Paper text: `\pn{}` for every number, "the reduced syn3A model" and not "Core A′", "detectable" and not "resolved", and no calibration claim on the reduced model.

## Review Focus

- A misspelled or jump-block species in `trajectory`: refused before any cell runs, with the name in the message (Task 1).
- A hybrid cell that stalls, in a matrix that records: it keeps `trajectory = nothing`, and bands skip it and count it (Tasks 1, 2).
- Pairs the viability rule rejected: weight 0, so absent from weighted bands, and present under `weights = ones(M, K)` (Task 2).
- Row blocks: `combine_rows` keeps every cell's trajectory, since panel A depends on it (Task 2).
- An `every` that is not a whole number of handshakes: refused before any cell runs, not after hours (Task 1).

---

### Task 1: Record trajectories in `predict`

Start: new worktree `.worktrees/trajectories` from `origin/main`, branch named from a new Linear issue (milestone "M6 Paper 1 complete"). Copy in the design and this plan from `.worktrees/kim-68`, and commit them first.

**Files:**
- Modify: `src/predict_forward.jl` (`CellRun`, `_failed_cell`, `_cell`, `_ForwardSetup`, `_forward_setup`, `_hybrid_cell`, `_mean_field_cell`, `forward_cells`, `intervention_cells`)
- Modify: `src/predict.jl` (`predict` keywords, refusals, docstring)
- Test: `test/test_predict_trajectory.jl` (create; add to `test/runtests.jl` after `test_predict_combine.jl`)

**Interfaces:**
- Produces: `predict(...; trajectory::Union{Nothing, AbstractVector{Symbol}} = nothing, every::Real = 60.0)`, and the same two keywords on `forward_cells` and `intervention_cells`.
- Produces: `CellRun.trajectory::Union{Nothing, NamedTuple{(:t, :species, :values)}}`, with `t::Vector{Float64}`, `species::Vector{Symbol}` and `values::Matrix{Float64}` (one row per time, one column per species). A failed cell holds `nothing`. Keep an 8-argument outer constructor, so existing test call sites stand.
- `every` sets the mean field's `saveat` only when `trajectory` is given. The default path keeps 60 s.

- [ ] **Step 1: Write the failing tests**

```julia
@testset "predict records trajectories" begin
    ms = collect(AbstractSubModel, m0_models())
    pd = _trajectory_test_posterior()   # three draws of M0_TARGETS at nominal values
    M = InferCell.LAC_EXPORT_METER
    S = [M, :M_gtp_c]
    a = predict(ms, pd; horizon = 120.0, n_cells = 2, rng = Xoshiro(7), threads = false)
    b = predict(ms, pd; horizon = 120.0, n_cells = 2, rng = Xoshiro(7), threads = false,
                trajectory = S, every = 60.0)
    @test isequal(a.mu, b.mu) && a.seeds == b.seeds
    @test all(c -> c.trajectory === nothing, vcat(a.cells...))
    c = b.cells[1, 1][1]
    @test c.trajectory.t == [0.0, 60.0, 120.0] && c.trajectory.species == S
    @test size(c.trajectory.values) == (3, 2)
    @test c.trajectory.values[end, 1] - c.trajectory.values[1, 1] ≈ c.q rtol = 1e-12
    @test_throws ArgumentError predict(ms, pd; horizon = 120.0, trajectory = [:M_nope])
    @test_throws ArgumentError predict(ms, pd; horizon = 120.0,
                                       trajectory = [transcript_state(:JCVISYN3A_0607)])
    @test_throws ArgumentError predict(ms, pd; horizon = 120.0, trajectory = S, every = 0.7)
    mf = predict(ms, pd; horizon = 120.0, forward = :mean_field, trajectory = S, threads = false)
    @test mf.cells[1, 1][1].trajectory.t == [0.0, 60.0, 120.0]
end
```

`_trajectory_test_posterior()` builds `PosteriorDraws(M0_TARGETS, θ, nothing, [1, 1, 1], [1, 2, 3], nothing, truth, (kind = :test,))`, with `θ` three rows near the nominal values, as `m0_evidence_kim40.jl` builds the truth's. Match the constructor's actual field meanings when writing it.

Also: a `_hybrid_cell` whose block stalls (reuse the stall set-up in `test/test_predict_forward.jl`) has `trajectory === nothing`.

- [ ] **Step 2: Run them to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test(test_args=["predict_trajectory"])'`, or include the file through a small Slurm job if the suite has no name filter.
Expected: FAIL, unknown keyword `trajectory`.

- [ ] **Step 3: Implement**

Add `trajectory::Vector{Symbol}` and `grid::Vector{Float64}` to `_ForwardSetup`, set by `_forward_setup(...; trajectory = Symbol[], every = 60.0)`; the grid is built as `saveat` is today. `_cell` takes the setup's two fields and reads the rows of `run.t` that match the grid, and the columns of `run.ode_names` that match the species. `predict` refuses before any work: species outside `_block_names(ms, :ode)`, `every <= 0`, and on the hybrid an `every` that is not a whole multiple of `driver.interval`. Thread the setup through `_hybrid_cell` and `_mean_field_cell`.

- [ ] **Step 4: Run the tests, and the frozen reference**

Run: a Slurm job over `test_predict_trajectory.jl`, `test_predict_reference.jl`, `test_predict_forward.jl`, `test_predict.jl`, `test_predict_viability.jl` and `test_predict_intervention.jl`.
Expected: all pass, with the reference unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/predict_forward.jl src/predict.jl test/test_predict_trajectory.jl test/runtests.jl
git commit -m "feat: predict records trajectories of ODE-block species"
```

### Task 2: Bands, the split over time, and the release recipe

**Files:**
- Create: `src/predict_trajectory.jl` (include from `src/InferCell.jl` after `predict_variance.jl`)
- Create: `docs/src/user-guide/prediction-uncertainty.md`
- Modify: `docs/make.jl`, `docs/src/api/predict.md`
- Test: `test/test_predict_trajectory.jl`

**Interfaces:**
- Consumes: `CellRun.trajectory`; `Prediction.weights`, `_variance_split`, `_moments`.
- Produces: `trajectory_bands(p::Prediction, species::Symbol; weights = first(keys(p.weights)), levels = (0.5, 0.9), baseline = false) -> (; t, n_failed, pairs, rows, cells)`. `pairs`, `rows` and `cells` are each `(; mean, lo, hi)`, with `lo`/`hi` a `Dict(level => Vector)`. `weights` is a weight-set name, or an `M × K` matrix (`ones(M, K)` is unconditioned). `baseline = true` subtracts each cell's value at t = 0, so a meter gives Q(t).
- Produces: `trajectory_split(p::Prediction, species::Symbol; weights = first(keys(p.weights)), baseline = false) -> (; t, held_share = (; i_first, h_first), total)`, each a vector over t: `_variance_split` at each t > 0, on μ̂(t) and τ(t) = σ̂²(t)/n.
- Weighted quantile: `_weighted_quantile(x, w, q)`, the type-7 definition generalised by cumulative weight, so that equal weights reproduce `Statistics.quantile`.
- Layers: `pairs` is the weighted quantile of μ̂(t) over pairs. `rows` is the weighted quantile of each row's weighted mean, with row weights `sum(w, dims = 2)`. `cells` pools every cell, each carrying `w[m, k] / n_ok[m, k]`.

- [ ] **Step 1: Write the failing tests**

```julia
@testset "bands and split over time" begin
    x = [1.0, 2.0, 3.0, 10.0]
    @test InferCell._weighted_quantile(x, ones(4), 0.9) ≈ quantile(x, 0.9)
    @test InferCell._weighted_quantile(x, [1.0, 1.0, 1.0, 0.0], 0.5) ≈ 2.0

    ms, pd = d11_models(), _trajectory_test_posterior()
    M = InferCell.LAC_EXPORT_METER
    p = predict(ms, pd; horizon = 120.0, completion = :prior, n_held = 2, n_parameter_draws = 3,
                n_cells = 2, rng = Xoshiro(3), threads = false, trajectory = [M])
    b = trajectory_bands(p, M; baseline = true)
    @test b.t[end] == 120.0 && b.n_failed == 0
    @test b.pairs.mean[end] ≈ p.summary[first(keys(p.summary))].mean rtol = 1e-12
    @test b.pairs.lo[0.9][end] ≈ quantile(vec(p.mu), 0.05) rtol = 1e-12
    @test all(b.pairs.lo[0.9] .<= b.pairs.lo[0.5] .<= b.pairs.hi[0.5] .<= b.pairs.hi[0.9])
    s = trajectory_split(p, M; baseline = true)
    @test s.held_share.i_first[end] ≈ variance_split(p).held_share.i_first rtol = 1e-10

    never = ViabilityRule(; label = "floor", allow_clipped = true,
                          floors = [:M_gtp_c => 1e12])
    r = predict(ms, pd; horizon = 120.0, completion = :viability_conditioned_prior,
                viability = never, n_viability_cells = 0, n_held = 2, n_parameter_draws = 3,
                n_cells = 2, rng = Xoshiro(3), threads = false, trajectory = [M])
    @test_throws ArgumentError trajectory_bands(r, M)          # nothing is weighted
    @test trajectory_bands(r, M; weights = ones(3, 2)).t[end] == 120.0

    e = predict(ms, pd; horizon = 120.0, completion = :prior, n_parameter_draws = 3,
                n_cells = 2, rng = Xoshiro(3), threads = false)
    @test_throws ArgumentError trajectory_bands(e, M)          # none was recorded
    bl = [predict(collect(AbstractSubModel, m0_models()), pd; horizon = 60.0, n_cells = 2,
                  draws = [i], rng = Xoshiro(i), threads = false, trajectory = [M]) for i in 1:2]
    @test combine_rows(bl).cells[2, 1][1].trajectory !== nothing

    # The release recipe: a fixed constant joins the held set, and its reverse is derived.
    rel = d11_models()
    rel[1] = CentralGlycolysis(enzymes = :translated,
        free = [:kcatF_R_ENO, :kcatR_R_ENO, :kcatF_R_FBA, :kcatR_R_FBA,
                :kcatF_R_LDH_L, :kcatR_R_LDH_L])
    h = held_parameters(metered(rel), pd)
    @test :kcatF_R_LDH_L in h.held && !(:kcatR_R_LDH_L in h.held)
    @test any(s -> s.name === :kcatR_R_LDH_L && s.reason === :haldane_derived, h.skipped)
end
```

The test posterior names M0's five targets, which `d11_models()` also frees. Adjust the `ViabilityRule` keywords to the type's own if they differ.

- [ ] **Step 2: Run them to verify they fail**

Run: the Slurm job from Task 1, over this file.
Expected: FAIL, `trajectory_bands` not defined.

- [ ] **Step 3: Implement `src/predict_trajectory.jl`, and export both functions**

- [ ] **Step 4: Write the docs**

`prediction-uncertainty.md` covers: recording a trajectory; the three layers; releasing a parameter with the module's `free =` (a Michaelis constant needs its reaction's `kcatR_*` freed too, and the build says so); that a release is a sensitivity to priors and not a posterior; and that one balanced kcat at geometric SD 1.05 barely moves a band. Add both functions to `api/predict.md`, and the page to `docs/make.jl` under the user guide. Build the docs per `docs/README` or `docs/make.jl`.

- [ ] **Step 5: Run the tests on Slurm, then the full suite (`test/run_tests.slurm`)**

Expected: all pass. Cite the job id.

- [ ] **Step 6: Commit, open the PR to `main`, gate with `/check-PR`, merge on a clean verdict**

```bash
git add src/predict_trajectory.jl src/InferCell.jl docs test/test_predict_trajectory.jl
git commit -m "feat: trajectory bands, the split over time, and the release recipe"
```

### Task 3: The driver, pre-registered, with a smoke run

Branch from `main` at Task 2's merge, in a new worktree.

**Files:**
- Create: `dev/scripts/trajectories_predict.jl`, `dev/scripts/trajectories_predict.sh`
- Modify: `spec/spec.md` (§12 entry)

**Interfaces:**
- Consumes: `predict(...; trajectory, every)`, `trajectory_bands`, `trajectory_split`; `pilot_common.jl` (`save`, `plain`, `say`, `SMOKE`, `COMMIT`); `launcher_guard.sh`.
- Produces, in `<outdir>`: `a_block_<b>.jls`, `a_truth.jls`, `b.jls`, `c_precheck.jls`, `c.jls`, and `tsv/traj_{a,b,c}.tsv` with columns `panel, arm, species, layer, weights, t, mean, lo50, hi50, lo90, hi90`, `tsv/traj_a_truth.tsv` (`species, t, mean, se`), `tsv/traj_a_data.tsv` (`species, t, value`), `tsv/traj_split.tsv` (`t, held_share_i_first, held_share_h_first, total, mc, flagged, mc_undefined, dropped`), `data.jls`, `replay.md`, `bulk_check.md` and `checks.md`.

| Stage | What it runs | Slurm |
|---|---|---|
| `a i`, i = 1…10 | `predict(m0_models(), pd; horizon = 600, n_cells = 887, draws = block_draws(pd, i), rng = Xoshiro(19_540_000 + (i − 1)), threads = false, trajectory = SPECIES)`, with `pd = posterior_draws(V7_DIR/replicate1.jls)` and `block_draws` as in `m0_evidence_kim40.jl` | array, 1 core, 6 GB, 14 h |
| `a 0` | the truth, as `m0_evidence_kim40.jl` builds it: `n_cells = 9_849`, `rng = Xoshiro(19_541_501)`, no intervention | 1 core, 6 GB, 14 h |
| `b` | as in the design's replay call, at 6,300 s | 64 threads, 192 GB, 8 h |
| `c_precheck` | `RELEASED` composition, `:viability_conditioned_prior`, the frozen `RULE`, draws `round.(Int, range(1, 556; length = 4))`, `n_held = 5`, 4 + 4 cells, `rng = Xoshiro(19_543_001)` | 16 threads, 48 GB, 2 h |
| `c` | the same, 40 draws, `n_held = 50`, `rng = Xoshiro(19_543_000)` | 64 threads, 192 GB, 24 h |
| `data` | `m0_dataset(16_200_001; ncells = 50)`, recording how its truth and σ_b differ from replicate 1's | 1 core, 8 GB, 1 h |
| `replay` | every check below, `replay.md`; gates `c` | 1 core, 16 GB, 1 h |
| `report` | the same checks, TSVs, `checks.md` | 1 core, 16 GB, 1 h |

`RELEASED`: `d11_models()`, with `ms[1]` and `ms[3]` rebuilt, freeing every currently fixed `kcatF_*`, `kcatR_*` and `km_*` plus ENO's and FBA's. Stop if `5 + length(p.held)` is not 82, or if any `kcatR_*` is held.

`replay` stops, and writes which value differed, unless all of these hold, the first two exactly (`==`):
- panel B: `p.mu == w1.mu`, and its seeds and θ_H, from `<archive>/runs/budget-kim42/w1.jls`;
- panel A: each block's `mu` equals the row in `<archive>/runs/m0-kim40-kim41/finish/*_r1_v7.jls` `.hybrid.mu`, and the truth's mean equals `truth` for dataset 1, posterior v7, in `score_datasets.tsv`, to its printed digits;
- the regenerated dataset's truth and σ_b equal replicate 1's recorded values, within 1e-12 relative, and `M_lac__L_e`'s meter ratio is constant across t on the truth, within 1% relative (1e-9 was too strong; see §12). Write that ratio R to `replay.md` and `checks.md`; it rescales the data onto Q's axis.

- [ ] **Step 1: Write the §12 entry**

Title it "Probabilistic trajectories: runs and figure". It cites the design, and fixes the seeds, sizes, released set, species, the checks above, the 10% rule, and that these runs add no scored claim.

- [ ] **Step 2: Write the driver and launcher**

Follow `reduced_budget_kim42.{jl,sh}` for stages, `put`, `read_out`, `missing` and `progress`. Under `SMOKE=1`: 630 s on the reduced model, 60 s on M0, 4 draws × 3 held × 2 + 2 cells, B7's `v6` fit, and checks reported without stopping.

- [ ] **Step 3: Parse and lower on the login node**

Run: `julia --project -e 'using InferCell; ex = Meta.parseall(read("dev/scripts/trajectories_predict.jl", String)); for a in ex.args; Meta.lower(Main, a); end'`, with `EBROOTJULIA` set.
Expected: no error and no `:error` or `:incomplete` expression in the result.

- [ ] **Step 4: Push, then smoke from a detached worktree at the pushed commit**

Run: `SMOKE=1 dev/scripts/trajectories_predict.sh <scratch outdir> run`, then `rest`, then read every TSV.
Expected: all stages complete, TSVs have the columns above, and bands are ordered.

- [ ] **Step 5: Commit and PR the driver and §12 entry to `main`, gated with `/check-PR`**

```bash
git add dev/scripts/trajectories_predict.jl dev/scripts/trajectories_predict.sh spec/spec.md
git commit -m "dev: trajectory runs, pre-registered, with exact replay checks"
```

### Task 4: The runs, the result note, and the archive

**Files:**
- Create: `dev/scripts/trajectories_predict_result.md`
- Output: `/fred/oz022/tkimpson/infercell-trajectories-2026-11-<dd>/`, with `SHA256SUMS`

- [ ] **Step 1: Launch `c_precheck`, `a`, `b` and `data` from a detached worktree at the merged commit**

Run: `dev/scripts/trajectories_predict.sh <outdir> run`; poll with `progress`.
Expected: jobs complete with exit 0.

- [ ] **Step 2: Read `replay.md`**

Expected: every check passes. On any mismatch stop, diagnose, and tell the author before more compute. Tell them the pre-check's acceptance and wall time as well, and set `CWALL` from the wall time if 24 h looks short.

- [ ] **Step 3: Launch `c`, then `report`, with `<outdir> rest`**

- [ ] **Step 4: Confirm the data match V7's commit**

Run: `<outdir> bulk <a detached worktree at dfe570a>`, which writes `bulk_check.md`: `m0_dataset(16_200_001; ncells = 50).bulk` against the regenerated `bulk`, exactly. A difference goes in the note and the caption says the data are regenerated at the later commit; it does not stop the work, since the posterior is unchanged.

- [ ] **Step 5: Write the result note**

Commit ids, job ids, the checks, cost against the estimate of about 500 CPU-h, panel C's acceptance, starved fraction and median move (pair means at 6,300 s, with and without the rule, against panel B's), `joint` beside `per_theta_I`, and layer 2's Monte Carlo share at 6,300 s. Write `SHA256SUMS` and copy the archive outside `.worktrees/`.

- [ ] **Step 6: PR the note to `main`; docs only, so ask the author whether `/check-PR` is wanted**

### Task 5: The figure and the paper

Branch from `paper-1-two-column` if #154 is still open, targeting it; otherwise from `main`.

**Files:**
- Create: `paper/figures/data/traj_{a,b,c,a_truth,a_data,split}.tsv` (copied from the run)
- Modify: `paper/figures/make_figures.py` (add `fig_trajectories`)
- Modify: `paper/numbers.tex`, `paper/sections/03-results.tex`, `paper/sections/00-abstract.tex`, `paper/sections/04-discussion.tex`, `paper/sections/02-infercell.tex`, `paper/sections/05-availability.tex`, `paper/STATUS.md`

- [ ] **Step 1: Load the `dataviz` and `publication-figure` skills, then write `fig_trajectories(data_dir)`**

Full width, 90% filled bands. Panel A: export (truth and data × R), GTP (truth, data and an outer band for a new bulk measurement, from σ_b's draws), pyruvate (truth only). Panel B: export and GTP with three nested layers, K_d lines at 2,018 and 2,160, and a sub-panel for the held share against time. Panel C: 20-parameter band against 82-parameter, with the unconditioned outline. "UNSCORED" on B and C, as in `fig_budget`. If the export band is not visible at A's scale, use a deviation from the posterior mean.

- [ ] **Step 2: Run the script, then assert the numbers**

Run: `/fred/oz022/tkimpson/miniconda3/bin/python3 paper/figures/make_figures.py`, then read the PDF as an image.
Expected: panel B's final share prints 0.908, and export's mean ends at the published 9.5765e6; assert both in the script.

- [ ] **Step 3: Add the text**

- **M0 results (the author's decision):** say plainly that the data include a noisy measure of Q through extracellular lactate, so Q's coverage is not a test on an unmeasured quantity, that how far it drives the posterior has not been measured, and that ΔQ is what the data did not see. Check the abstract, the "What this supports" paragraph and the discussion for the same reading, and reword them to match.
- **The figure's paragraph,** with the caption ruling out the anti-goals: A does not validate B or C; no band predicts a real cell; C is a sensitivity to priors, on a posterior fitted with those constants at nominal values.
- **`\pn{}` entries** for every new number, each with its source file.
- **The listing** of the release recipe in `02-infercell.tex`, by the predict API section, checked against the test in Task 2.
- **Availability and `STATUS.md`:** the archive path and the new figure.

- [ ] **Step 4: Build, and confirm the PDF text changed**

Run: `latexmk -pdf main.tex` in `paper/`, from a clean output directory, with TinyTeX on `PATH`.
Expected: no undefined references, and the new caption found by `pdftotext`.

- [ ] **Step 5: Run `/deslop` on the new text and `/check-arxiv-llm-compliance` on the paper, then commit and PR, and gate with `/check-PR`**

```bash
git add paper
git commit -m "paper: probabilistic trajectories, and what M0's Q score tests"
```
