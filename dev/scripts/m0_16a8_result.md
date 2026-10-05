# 16a.8: M0 and its datasets

Run of record: job 17800819 at `1fec8b8` (`dev/scripts/phase16a_tests.slurm`),
`test/test_m0.jl`, 21 of 21 pass. The same job reran 16a.1 to 16a.4, all passing.

## M0

- **The model.** Transcription, decay and translation run for GAPD
  (`JCVISYN3A_0607`) and ptsG (`JCVISYN3A_0779`) over 600 s. The four ODE modules
  are Core A′'s, with ENO and FBA freed as in `d11_models`. It builds in
  completeness mode.
- **The held proteins.**
  - `HeldEnzymeCounts`, a jump module with no reactions, holds the other 12
    catalytic counts at their proteomics copy numbers.
  - ptsI, ptsH and crr are not produced, so their totals are held.
  - Both holds are M0 reductions, and `reduction_report` names all twelve loci
    and the three carriers.
  - Over a 120 s run the held counts do not move.
- **Size and cost.** 32 ODE states, 34 jump states and 7 jumps. A 600 s cycle
  takes 4.19 s warm, on the login node. A drawn truth had 319 firings.

## M0 datasets (`m0_dataset`)

One seed fixes a replicate:
- a truth over `M0_TARGETS` (purpose `:calibration`), drawn Haldane-consistently;
- σ from `SIGMA_MET_PRIOR`, `LogNormal(log 0.2, 1)`;
- the cell seeds and the noise seed.

The observations are the 15.7 panel's 17 metabolites, lognormal at σ with a
one-particle floor, and M0's two transcripts, exact. t = 0 is dropped.

Tested:
- a replicate carries its truth, σ, seeds and report;
- the transcripts are exact;
- two runs at one seed are identical;
- a different seed draws a different σ and truth.

Two existing functions changed:
- `generate_dataset` gains `names`, passed to `check_truth`, which defaulted to
  D11's six and would refuse an M0 truth.
- `emit_ensemble` takes the transcripts the composition has, rather than
  assuming all 17.
