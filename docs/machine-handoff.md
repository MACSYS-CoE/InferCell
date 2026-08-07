# Machine handoff

_Packed 2026-08-06 on Tom's macOS laptop, for resume on **OzStar**.
To resume: `/switch-machine resume`._

_Lives in `docs/` rather than the repo root because `.gitignore` deliberately
ignores root-level `*.md` — same reason `handoff.md` is in `docs/`._

## Environment

- **Julia**: `Project.toml`, `Manifest.toml` — both committed ✓
  setup: `julia --project -e 'using Pkg; Pkg.instantiate()'`
  On OzStar this must be done from the **login node** — compute nodes have no
  route to the Julia registry. See `test/run_talk_bursty_boundary.slurm`, which
  documents the exact incantation (`JULIA_PKG_PRECOMPILE_AUTO=0`, then
  precompile inside the job; `module load` does not work in non-interactive
  shells and the Julia binary is a wrapper needing `EBROOTJULIA`).

- **Python (figures)**: **DRIFT — no spec file exists.**
  `dev/talks/ISAB/figures/plot_bursty_boundary.py` needs `corner`, `numpy` and
  `matplotlib`, and none of them are declared in any tracked file.
  `docs/requirements.txt` is mkdocs-only and unrelated.
  setup: `pip install corner numpy matplotlib`
  Worth turning into a real spec file at some point; not done here because this
  step is verify-only.

- **LaTeX**: building the deck on OzStar needs `~/.TinyTeX` on PATH — the
  system TeX has no `beamer`.

## Data

- [x] `dev/talks/ISAB/figures/*.png` — small, tracked — travels in git.
      **But all three ISAB corner figures are STALE** (old palette, and the
      blue-only export does not exist yet). They must be regenerated — see
      Regeneration below. This is the session's one blocker.
- [ ] `examples/results/isab_*.csv` — cluster-resident, **gitignored**
      (`.gitignore` line 28: `examples/`) — does **not** travel. Written by job
      15152954 on OzStar. Since you are resuming *on* OzStar, no transfer is
      needed: find them under `examples/results/` in the OzStar checkout, or
      re-run the job.
- [ ] `dev/talks/ISAB/talk.pdf` — untracked, regenerable — left untracked
      deliberately (carried-over open question from the previous handoff).
      Rebuild with `latexmk -pdf talk.tex` in `dev/talks/ISAB/`.

## Cluster access

- host: `ozstar.swin.edu.au`  user: `tkimpson` (already in `~/.ssh/config`)
- No rsync needed for this handoff — you are travelling *toward* the data, not
  away from it. If you later need the CSVs back on the laptop:
  `rsync -avz ozstar.swin.edu.au:<repo>/examples/results/ ./examples/results/`
  (absolute remote repo path not recorded — fill in on arrival.)

## Regeneration

The one thing standing between the deck and presentable. On OzStar, from the
repo root, once the CSVs are in `examples/results/`:

```
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ode
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode alone
python3 dev/talks/ISAB/figures/plot_bursty_boundary.py --mode ssa
```

`alone` and `ssa` **must** come from the same particles, or slides 11 and 12
disagree with each other.

If the CSVs are gone, regenerate them first:
`sbatch test/run_talk_bursty_boundary.slurm` (production settings — do not
quote numbers from a smoke run; at 100 NUTS samples the chain does not converge
and the boundary passes the overconfidence downstream).

Then rebuild: `cd dev/talks/ISAB && latexmk -pdf talk.tex` (needs TinyTeX).

## Secrets

None detected in the working tree. Git access is over SSH; the OzStar checkout
needs its own key already configured for `git@github.com:MACSYS-CoE/InferCell`.

## To resume

Run `/switch-machine resume` on OzStar.
