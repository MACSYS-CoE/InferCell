# Handoff — 2026-08-06 16:20

## Goal

Make the ISAB talk's results legible to a general audience by building the deck
on a genuinely three-module minimal cell instead of a rig whose "three modules"
were two.

## Status

Done and pushed. PR #37 open against `main` from `minimal-cell-example`. The
deck compiles clean at 19 pages (13 main line, 5 backup) and every figure comes
from one production run. Nothing is half-finished.

## What changed this session

See `dev/log/2026-08-06-isab-minimal-cell.md` for the narrative, and
`docs/superpowers/specs/2026-08-06-minimal-cell-boundary-figure-design.md` for
the design and the full measured tables.

## Open questions

- **Commit `dev/talks/ISAB/talk.pdf`?** Both previous talk directories track
  their built PDF, but those are finished talks. This one is six days out and
  the committed build came from TinyTeX on the cluster, not the user's local
  TeX. Left out deliberately; trivial to add.
- **Vector figures?** The deck uses the `.png` files because that matches what
  the repo tracks for talk figures. `.pdf` versions of all three exist on disk
  (untracked) if the projector warrants them.

## Blockers / problems

None blocking the talk.

One honest limitation, disclosed on backup slide 16 rather than hidden: the
regulator's promoter rates (`k_on`, `k_off`, `k_tx_burst`) are not individually
identifiable from the current summary statistics, because replicate-mean
trajectories see a telegraph promoter only through
`k_tx_burst·k_on/(k_on+k_off)`. That combination *is* recovered and the boundary
sharpens it 12×. Fixing it properly means adding distributional summaries (Fano
factor, autocorrelation time, full count distribution) to
`src/summary_statistics.jl` — the single most valuable follow-up, and arguably
the sharpest answer to "why carry a stochastic module at all?"

## Next steps

1. Review PR #37 and merge.
2. Rehearse the setup→results pair (slides 10→11). The fix for the old awkward
   transition is that slide 10 now *asks* the question slide 11 answers — say
   the question out loud before advancing, or the repair is wasted.
3. Decide on `talk.pdf` and vector figures (above).
4. After the talk: distributional summary statistics, then the growth/division
   module that was deferred this session.

## Non-obvious context

- **Reproducing the figures:** `sbatch test/run_talk_bursty_boundary.slurm`,
  then run `dev/talks/ISAB/figures/plot_bursty_boundary.py` in its three modes
  (`ssa`, `ode`, `burst`). Packages must be downloaded from the **login node**
  first — compute nodes have no route to the Julia registry; the `.slurm` file
  documents the exact command.
- **Parameter order in the boundary figure is load-bearing.** The three unshared
  parameters are exported first and the three shared ones last, so the
  tightening lands in the bottom-right block and the dashed box sits exactly on
  the coupling. The plot script warns and drops the box rather than drawing a
  wrong one if that ordering ever changes — check the job log after any rerun.
- **Do not trust smoke-setting runs here.** At 100 NUTS samples the
  differentiable chain does not converge, and the boundary faithfully passes the
  resulting overconfidence downstream, which makes the result look far better
  than it is. Both run logs are in the design doc; only the production numbers
  are quotable.
- **Building the deck on the cluster** needs `~/.TinyTeX` on PATH — the system
  TeX has no `beamer`. The user's usual build is on their laptop.
- `examples/` is gitignored, so the exported CSVs behind the figures are not in
  the repo. Rerun the job to regenerate them.
