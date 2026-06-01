# Handoff — 2026-05-23

## Goal

Issue #11 — produce a multi-seed empirical calibration curve for the iterative
Level-2 boundary protocol (`iterative_infer` from PR #10) and decide whether it
is honestly calibrated. PR #10's single-seed result hinted at over-confidence;
this issue's job is to confirm or rule out a systematic bias across 20 seeds.

## Status

Branch `issue-11-calibration-sweep` (off `main`). Work complete:

- 20-seed sweep run on Slurm (4 × 5-seed array, milan partition).
- Calibration curve at `examples/figures/step4_calibration_curve.png` (also
  mirrored to `dev/notes/figures/` so it renders in the writeup).
- Interpretation note at `dev/notes/2026-05-22-calibration-curve.md`.
- PR pending: only the writeup + figure + `.gitignore` allowlist tweak go
  into git. The sweep scripts (`examples/calibration_sweep.{jl,slurm}`,
  `examples/step4_calibration_curve.py`) live in the gitignored `examples/`
  scratch area per project convention.

## Result

**Iterative is systematically over-confident** (empirical coverage far below
nominal across all 4 ODE parameters at all 4 nominal levels). The single-pass
baseline is honestly calibrated. Full table and diagnostics in
`dev/notes/2026-05-22-calibration-curve.md`.

Two diagnostics worth highlighting:

- **No seed converged.** All 20 hit `max_iters=4`. Final-step KL has median
  0.568 across seeds (well above `kl_tol=0.05`). The loop is drifting, not
  converging, and the drift tightens posteriors past truth.
- **Iterative costs ~4×** the compute of single-pass (mean 1023 s vs 251 s
  per seed) for a *worse-calibrated* posterior.

## What changed this session

- `examples/calibration_sweep.jl` — per-seed driver. Reads
  `SLURM_ARRAY_TASK_ID`, processes 5 seeds, writes one CSV per shard. Mirrors
  the data-generation pattern in `test/run_integ_iterative.jl`. Writes
  incrementally so timeouts don't lose completed seeds (learned this the hard
  way on the first array submission). Supports `--smoke` for fast local
  validation.
- `examples/calibration_sweep.slurm` — Slurm array wrapper, `--array=0-3`,
  4 h walltime (the first 2 h attempt timed out on shards 0 and 3).
- `examples/step4_calibration_curve.py` — plot script. Globs shards,
  computes empirical coverage per (protocol, param, level), renders 1×2
  subplots with diagonal reference. Mirrors `step4_plot_results.py` style
  (LaTeX, Wong palette, frameless legend).
- `dev/notes/2026-05-22-calibration-curve.md` — writeup with verdict,
  coverage table, mechanism hypothesis (KDE tail truncation in
  cut-posterior loops), follow-up priorities.
- `.gitignore` — allowlist entries for the new note + `dev/notes/figures/`.

## Scope clarification (from user, post-sweep)

`iterative_infer` is a proof-of-concept of cross-module information exchange,
**not** the production protocol. The calibration finding is informational —
the documented follow-ups are *not* queued and should not be pursued unless
the user explicitly changes priority. Single-pass `sequential_infer` is the
honestly calibrated default.

## Next steps

1. PR #21 is open. Merge after review.
2. No iterative-protocol calibration follow-ups planned. If a future
   session is asked to improve `iterative_infer`, check with the user
   first — the protocol may be replaced wholesale rather than fixed
   incrementally.

## Non-obvious context

- **`examples/` is intentionally gitignored** as scratch. The sweep scripts,
  shard CSVs, and PNG live there. The PR ships only the writeup + figure.
  If a future session wants to re-run the sweep, the scripts are on the
  cluster at `/fred/oz022/tkimpson/InferCell/examples/`.
- **The driver writes CSV incrementally** (one row per quantile, flushed
  after each protocol-per-seed). Re-run a single shard with
  `sbatch --array=N examples/calibration_sweep.slurm` if a shard is missing.
- **Sharding scheme**: shard `i` processes seeds `1001 + 5*i` ... `1005 + 5*i`.
  Output goes to `examples/results/step4_calibration_shard_$(i).csv`.
- **`handoff.md` (repo root) is the project's handoff convention.** This file.
  Moved from `docs/handoff.md` in the 2026-05-29 docs restructure so that
  `docs/` holds only published site content. Keep updating it here at root
  (the path the /handoff skill defaults to).
- **`docs/` now holds only published, user-facing site content** (the eventual
  released documentation): `index`, `getting-started`, `positioning`,
  `user-guide/`, `api/`. Everything else lives under `dev/`: `plans/` (the
  design/decision log), `notes/` (research notes incl. the case-for-inference),
  `talks/`, `log/`, `council/`, and `references/`. The dated design docs and
  research notes are deliberately *not* published — they're dev artifacts and may
  be dropped from the released docs. The old `docs/**` + per-file allowlist in
  `.gitignore` is gone; `docs/` is tracked normally, and only build junk plus
  `dev/council/` and `dev/references/` are ignored. Adding a published doc just
  means adding it to the mkdocs `nav`.
