# Phase 16 compute ledger

The cap is 50k CPU-h for the phase, and 5k for 16a (`spec/phases/16-recovery.md`
D16.6). CPU-h is wall time × cores, as reported by Slurm.

| Job | Task | What | Cores | Wall | CPU-h |
|---|---|---|---|---|---|
| 17754678 | 16a.1 | V2/V2b tests, first run | 1 | 0:03:39 | 0.06 |
| 17754679 | 16a.1 | Full-cycle V2, 3 cells, first run | 1 | 0:08:16 | 0.14 |
| 17755845 | 16a.1 | V2/V2b tests (of record) | 1 | 0:03:40 | 0.06 |
| 17755846 | 16a.1 | Full-cycle V2 and V2b (of record) | 1 | 0:08:25 | 0.14 |
| 17756205 | 16a.2 | 16a tests, identity failed (toy) | 1 | 0:04:21 | 0.07 |
| 17756874 | 16a.2 | 16a tests (of record) | 1 | 0:04:03 | 0.07 |

| 17758130 | 16a.3 | 16a tests, name clash | 1 | 0:04:07 | 0.07 |
| 17758372 | 16a.3/4 | 16a tests, cap bound failed | 1 | 0:04:29 | 0.07 |
| 17758553 | 16a.3/4 | 16a tests (of record) | 1 | 0:05:27 | 0.09 |
| 17800819 | 16a.8 | 16a tests with M0 (of record) | 1 | 0:05:28 | 0.09 |
| 17801360 | 16a.5 | 16a tests with block 1 | 1 | 0:06:43 | 0.11 |
| 17801198 | 16a.2 | V1 on M0, smoke (2 tasks) | 1 | ~0:04 each | 0.13 |
| 17801682 | 16a.2 | V1 on M0, 10⁴ paths (100 tasks) | 1 | ~0:28 each | 34.74 |
| 17801857 | 16a.6 | 16a tests with 16a.6 | 1 | 0:05:34 | 0.09 |
| 17801858, 17802079, 17802347 | 16a.6 | Core A′ clip scans (ENO) | 1 | | 0.29 |
| 17801684–17801688, 17802939–17802941, 17814376 | 16a.5 | V4 grids 1 to 3, merges, 20 chains | 1 | | 29.77 |
| 17814374, 17814375, 17830703 | 16a.6 | Core A′ clip scans (FBA, PGK3) | 1 | | 0.65 |
| 17814555, 17814556, 17814572, 17830702, 17830938 | 16a.5 | V4 grid 4, merge, finals, target check | 1 | | 4.77 |
| 17839207–17839209 | 16a.5 | V4 FBA-fine grid, merge, final (of record) | 1 | | 10.00 |
| 17839362, 17839571, 17840103, 17854107 | 16a.7 | 16a tests with CSMC (corrected after the #81 review: the two earlier rows overlapped) | 1 | | 0.41 |
| v5s_* (two smoke rounds) | 16a.7 | V5 smoke, both cases | 1 | | 0.67 |
| 17854476–17854480 | 16a.7 | V5, M0 (of record) | 1 | | 41.52 |
| 17854481–17854485 | 16a.7 | V5, toy (of record) | 1 | | 5.30 |
| 17877524–17877529 | 16a.7 | Variant smoke, N = 5, three windows | 1 | | 0.42 |
| 17877760, 17877761, 17877922, 17877923 | 16a.7 | Weight-collapse diagnostic (scratch) | 1 | | 0.26 |
| 17878615, 17878616, 17879029, 17879030 | 16a.7 | Collapse attribution, two passes (of record) | 1 | | 0.42 |
| 17879419, 17879420, 17880205, 17880206 | 16a.7 | Rank check, one failed run and the run of record | 1 | | 0.82 |
| 17878617–17878622 | 16a.7 | KP measurement at N = 50 (cell 4 failed under ancestor sampling) | 1 | | 11.00 |
| 17882299 | 16a.7 | 16a tests with the impossible-firing fix (of record) | 1 | 0:06:16 | 0.10 |
| 17882309–17882311 | 16a.7 | KP, cell 4 under ancestor sampling, rerun | 1 | | 2.30 |
| 17884470 | 16a.7 | KP merge (of record) | 1 | 0:01:22 | 0.02 |

| 17904315, 17904316 | 16a.7b | Rank check, five full-scale cells, and merge (of record) | 1 | | 2.56 |
| 17904432 | 16a.7a | 16a tests with annealing (of record) | 1 | 0:13:27 | 0.22 |
| 17905020–17905026 | 16a.7a | Annealed windows on M0, K ∈ {20, 50, 100}, N ∈ {2, 5} (of record) | 1 | | 22.01 |
| 17905037–17905039 | 16a.7a | V5 with annealing at K = 5, M0 (of record) | 1 | | 109.17 |
| 17905042–17905044 | 16a.7a | V5 with annealing at K = 5, toy (of record) | 1 | | 12.54 |
| 17905489 | 16a.7a | Inside one annealed window, K up to 1,000 (of record) | 1 | 1:20:57 | 1.35 |

| 17933661, 17933662 | 16a.7a | Collapse attribution at the cited detection floor (of record) | 1 | | 0.28 |

| 17934962 | 16a.7 | Bulk identifiability from 15.8's run (of record) | 1 | 0:03:31 | 0.06 |

| 17938700–17938702 | 16a.7 | F6 ENO cell: full-path information bound, 40 cells (of record) | 1 | | 5.07 |

| 18049605 | 16a.7c | 16a tests with bulk metabolites (of record) | 1 | 0:12:28 | 0.21 |
| 18052532–18052536 | 16a.7c | Bulk update rates on M0, N ∈ {5, 10, 20, 50} (of record) | 1 | | 35.36 |
| 18052537–18052542 | 16a.7c | V5 bulk half, M0 and toy (of record) | 1 | | (in the line above) |

| 18060223 | 15.7b | Bulk observations for 15.7, two cells resimulated (of record) | 1 | 0:03:22 | 0.06 |

| 18060205–18060213, 18061130–18061135 | 16a.7c | Fresh V5 references, M0 and toy, and the reference checks (of record) | 1 | | 25.01 |

| 18060214–18060217, 18061281, 18061355–18061363 | 16a.7c | V4 on block 1's bulk form: grids 1 to 3, FBA- and ENO-fine grids, chains, check, final (of record) | 1 | | 45.01 |

| 18098390 | 16a.7c | Per-cell update-rate merge (of record) | 1 | | 0.08 |

| 18105666 | 16a.9a | Full suite with the Gibbs chain (Aqua compat failed) | 4 | 0:50:42 | 3.38 |
| 18109288 | 16a.9a | Full suite with the Gibbs chain (of record) | 4 | 0:54:07 | 3.61 |
| 18108335–18108339 | 16a.9a | Gibbs pilot, 4 chains × 200 sweeps on M0, and merge (of record) | 20 | ~11 h each | 896.33 |

| 18157561–18157565 | 16a.9b | V6 stage 1: the pilot's 4 chains, 200 → 700 sweeps, and merge (of record) | 20 | ~28 h each | 2,308.28 |

| 18255035 | 16a.9b | B3 profile, chain 1 at sweep 700 (of record) | 20 | 0:12:10 | 4.06 |
| 18256434–18256437, 18258719, 18258720 | 16a.9b | Thread scaling at 1, 5, 10, 20 threads; the first 1- and 5-thread jobs cancelled by root (of record) | 1–20 | | 14.28 |

| 18271545 | 16a.9b | V6 stage 1 merged with 100 sweeps of burn-in (of record) | 1 | 0:04:09 | 0.07 |

| 18276616, 18276617 | 16a.9b | V7 smoke: 3 replicates, 3 cells, 2 windows, 6 sweeps, and merge | 5 | ~0:06 each | 1.47 |

| 18276511–18276515 | 16a.9b | V6 stage 2: 4 chains, 700 → 1,552 sweeps, and merge | 20 | ~50 h each | 4,038.65 |
| 18420986–18420990 | 16a.9b | V6 top-up: 4 chains, 1,552 → 1,650 sweeps, and merge (of record) | 20 | 5.4–6.9 h each | 462.55 |
| 18421022 | 16a.9b | V6 stage 1's table regenerated from the extended checkpoints | 1 | 0:04:15 | 0.07 |
| 18475428, 18475516 | 16a.9b | V7's patched merge smoked on the smoke checkpoints; the second of record | 1 | 0:02:37, 0:01:23 | 0.07 |
| 18276794 (50 tasks) | 16a.9b | V7: 50 replicates, 34 to 656 sweeps and 16 to the 6-day limit | 5 | 4.5–6.0 d each | 33,296.72 |
| 18628508 (16 tasks) | 16a.9b | V7: the 16 timed-out replicates resumed to 656 sweeps | 5 | 0:29–18:01 each | 626.85 |
| 18689602 | 16a.9b | V7 merge on the complete set (of record) | 1 | 0:01:43 | 0.03 |

| 18510849, 18511245, 18511832, 18512048, 18512335 | 16a.11 | ODE solver benchmark, smoke runs | 2 | 0:04–0:09 each | 0.96 |
| 18512579 | 16a.11 | ODE solver benchmark (of record) | 16 | 0:25:45 | 6.87 |
| 18751515 | 16a.10 | F10 production chain, smoke: 3 cells, 120 s, 2 sweeps | 5 | 0:03:15 | 0.27 |
| 18751635 | 16a.10 | F10 production chain 5 on V6's dataset, 656 sweeps (of record) | 5 | 4-05:30:05 | 507.51 |

**16a total, excluding 16a.9: about 919.2 CPU-h** (the rows' sum, recomputed after
the #81 review, which found the running total had drifted; 403.6 before the ODE
solver benchmark's 7.8, and 411.4 before F10's 507.8).

**16a.9: about 41,656.5 CPU-h**, recorded but not held to 16a's 5k cap (§12
2026-10-06). V7 is 33,923.6 of it, against 25.6k sized and 35.7k projected (§12
2026-10-08, annotated). **Phase 16 total before the v0.3.0 tracks: about 42,575.7 CPU-h** of 50k, so about
7.4k remain.

## v0.3.0 track A (#88), charged to phase 16's remaining cap

| Job | Task | What | Cores | Wall | CPU-h |
|---|---|---|---|---|---|
| 18818909–18818912, 18819650, 18819651, 18819611–18819627 | #88 | A1, A2 and A3 smoke runs and merges (120 s horizon) | 1 | ~0:03 each | 2.20 |
| 18819106, 18819610, 18821310 | #88 | Test suite: Aqua missing from the depot; one failure in the new rebuild test; passing (of record, 28361 of 28361 at `6da2398`) | 4 | 0:00:20, 0:55:06, 0:50:19 | 7.05 |
| 18821440, 18821441 | #88 A1 | M0 ensemble sensitivity, 1470 seeds over 49 tasks, and merge (of record) | 1 | ~0:14 each | 11.10 |
| 18821442, 18821443 | #88 A2 | Core A′ ensemble sensitivity (15.8 rerun), 1470 seeds over 294 tasks, and merge (of record) | 1 | up to 0:37 each | 165.85 |
| 18825531–18825547 | #88 A3 | Reverse-feedback arms: 8 ensembles of 200 cells, 8 path arrays of 40 cells, and merge (of record) | 1 | up to 0:19 each | 86.97 |
| 18836336 | #88 | Test suite after the #92 review (of record, 28362 of 28362 at `2b039a0`) | 4 | 0:50:06 | 3.34 |
| 18836373, 18836375, 18836376 | #88 | A1, A2 and A3 re-merged after the #92 review: every number identical | 1 | ~0:02 each | 0.10 |
| 18836377, 18836378, 18836380–18836396 | #88 | A1 and A3 smoke reruns after the #92 review: A1's task files bitwise equal to the first smoke's | 1 | ~0:03 each | 2.57 |

**#88 total: about 279.2 CPU-h**, against about 320 projected on the issue.
**Phase 16 total with #88: about 42,854.9 CPU-h** of 50k, so about 7.1k remain.

## v0.3.0 track B (#89), charged to phase 16's remaining cap

| Job | Task | What | Cores | Wall | CPU-h |
|---|---|---|---|---|---|
| 18884570 | #89 B1/B2 | Test suite with the transcript factor and the mean field (28627 of 28627) | 4 | 0:53:33 | 3.57 |
| 18892097 | #89 B1/B2 | Test suite after the #93 review (of record, 28630 of 28630 at `c8c9ceb`) | 1 | 1:17:25 | 1.29 |
| 18885264–18885272 | #89 B3 | Jensen-gap smoke: 120 s, 4 cells per M0 truth, mean fields, merge | 1 | 0:01:33–0:08:06 each | 0.86 |
| 18886224–18886232 | #89 B3 | Jensen gap: three M0 ensembles of 200 cells, five mean-field cells, merge (of record) | 1 | 0:01:05–0:07:58 each | 2.20 |
| 18888681–18888686 | #89 B3 | Rounding check: five mean-field cells under fractional carry, merge (of record) | 1 | 0:01:09–0:05:40 each | 0.46 |
| 18892098, 18892099 | #89 B3 | Diagnosis over both archives (of record) | 1 | 0:02:01–0:02:02 each | 0.07 |
| 18892323–18892332 | #89 B3 | Smoke of the #93 review fixes: every stage, through the diagnosis | 1 | 0:00:46–0:04:14 each | 0.70 |
| 19002806 | #89 | Test suite after merging main with #92 (of record, 28640 of 28640 at `ec4c219`) | 4 | 0:50:08 | 3.34 |

| 19012789–19012801 | #89 rethink | Smoke of every stage (4 cells, 120 s) | 1 | 0:00:53–0:03:44 each | 1.11 |
| 19013270–19013346 | #89 rethink | Stage 1 arms, 2A ensembles on shared seeds (superseded), 32 mean fields, score | 1 | up to ~0:10 each | 21.42 |
| 19014973–19014975 | #89 rethink | Seed-range check at f10_ref, m0_low, m0_high | 1 | 0:04:30–0:05:11 each | 1.21 |
| 19015472–19015544 | #89 rethink | 2A on independent seeds, twenty posterior-draw cells; 76 tasks failed on a precompile cache from another CPU (excludes foreign tasks 19015474_1563 and 19015507_579 in the id range) | 1 | 0:00:34–0:15:29 each | 18.50 |
| 19017428, 19017558–19017599 | #89 rethink | Precompile on milan; the missing pieces rerun; transfer (superseded by 19018422), calibrate (of record), score (of record for Stage 1; its 2A tables are the shared-seed run) | 1 | up to ~0:10 each | 5.75 |
| 19018421, 19018422 | #89 rethink | Precompile; transfer rerun with the floor variant (of record) | 1 | 0:01:24, 0:01:05 | 0.04 |
| 19038595 | #89 B3 close-out | Implied-bias report, first run (superseded: m0_v6 section lacked SEs) | 1 | 0:01:42 | 0.03 |
| 19038794 | #89 B3 close-out | Implied-bias report, second run (superseded: unpaired m0_v6 noise) | 1 | 0:00:54 | 0.02 |
| 19039511 | #89 B3 close-out | Implied-bias report (of record, `beeb446`) | 1 | 0:01:21 | 0.02 |
| 19042130 | #89 B4 | Route (ii)'s new tests alone (dev check; FAILED on one test's wrong save time, fixed before `67c9eef`) | 1 | 0:03:15 | 0.05 |
| 19042448 | #89 B4 | Full suite at `67c9eef` (28709 of 28709) | 4 | 0:55:34 | 3.70 |
| 19042449 | #89 B4 | Precompile (out of memory; its stages never ran) | 1 | 0:00:38 | 0.01 |
| 19055277–19055298 | #89 B4 | Precompile, route (ii) and three prototype steps at five truths, merge (`5a56836`; of record except the merge, superseded by 19058585, and corea's 0.25 s prototypes 19055293 and 19055297, which stopped on MaxIters) | 1 | 0:01:17–0:04:20 each | 1.02 |
| 19057735 | #89 B4 | Extrapolated comparison (superseded: corea's 0.25 s missing) | 1 | 0:01:25 | 0.02 |
| 19058579, 19058580, 19058585, 19058586 | #89 B4 | corea's 0.25 s prototypes with the step cap lifted, merge, extrapolation (of record, `e7fb020`) | 1 | 0:00:51–0:04:44 each | 0.18 |
| 19076053 | #89 B5 | Full suite at `4fe4022` (28787 of 28787; superseded by 19084578) | 4 | 0:58:54 | 3.93 |
| 19076054, 19076055 | #89 B5 | Smoke r1: precompile, `theta_ref` (FAILED: the launcher exported NADAPT empty; fixed in `c77f49a`). Its later stages were cancelled before they ran | 1 | 0:02:06, 0:01:05 | 0.05 |
| 19076298, 19076299 | #89 B5 | Smoke r2: precompile, `theta_ref` (cancelled after 2:05:49: the MAP stepped to an extreme θ where one solve does not return; fixed by the support box, `b101190`, §12) | 1 | 0:01:46, 2:05:49 | 2.13 |
| 19081750–19081758 | #89 B5 | Smoke r3 at `b101190`: precompile, `theta_ref`, two ensemble chunks, four chains of 10 + 10, merge (all stages ran; numbers not of record) | 1 | 0:01:26–0:44:57 each | 2.69 |
| 19084578 | #89 B5 | Full suite at `469ca26` (28791 of 28791, of record) | 4 | 0:57:54 | 3.86 |
| 19091859–19091881 | #89 B5 | The full run at `469ca26`: precompile, `theta_ref` (bounded MAP), 16 ensemble chunks (800 cells), four chains of 300 + 300, merge (of record, `dev/scripts/mixed_nuts_89_result.md`) | 1 | 0:01:31–10:27:54 each | 41.48 |
| 19159993, 19160196, 19161734, 19162737, 19162281, 19163602 | #89 B7 | Gradient profile: snapshots at `ddddec5` (the first invalid: it compiled an edited tree), `d39360c` and `6d1320a`, and two compares | 1 | 0:00:45–0:04:13 each | 0.30 |
| 19160188, 19163839, 19163842, 19164153 | #89 B7 | Test files: `test_mixed_likelihood.jl` at `48f1e99` (107 of 107); the mean-field set at `6d1320a` (19163842 FAILED, helpers not included; 19164153 489 of 489) | 1 | 0:00:04–0:11:18 each | 0.40 |
| 19163848–19163870 | #89 B7 | AV6 smoke r1 at `51064a8` (2 replicates, v6, KS, stream arms). The stream arms FAILED (no fit directory; fixed in `ba69134`), so the merge (19163871) never ran | 1 | 0:01:03–0:21:05 each | 4.46 |
| 19171825–19171849, 19177272 | #89 B7 | AV6 smoke r2 at `ba69134`, every stage and the merge; 19177272 reran the merge at `56c573d` with `AV6_COMPARE=1` (numbers not of record) | 1 | 0:01:04–0:22:18 each | 4.86 |
| 19177653–19177825 | #89 B7 | AV6 at `7a102fc`: 50 replicates, v6, three KS refits, two stream arms, merge (167 jobs; running) | 1 | | |
| 11 preps of 19177653–19177825 | #89 B7 | AV6 preps r5, r7, r18, r23, r25, r28, r29, r33, r37, r42, r50: TIMEOUT at 12 h in the bounded MAP (§12 2026-10-25); listed apart as wasted; the row above's CPU-h will include them, so not summed twice | 1 | 12:00:04–12:00:27 each | 132.06 |
| 19200562–19200603 (40 of 42; 19200569, 19200589 are not #89) | #89 B7 | First resume at `7a102fc` (3-day preps). CANCELLED once the MAP trace showed the gradient test would not pass; 11 preps ran 11:16:53 each, the rest never started | 1 | 11:16:53 | 124.10 |
| 19209034, 19209035 | #89 B7 | MAP traces at `af332d7`: r5 (to its 6 h cap) and r46 | 1 | 6:01:28, 0:09:52 | 6.19 |
| 19232720–19232727, 19232740 | #89 B7 | MAP checks of the reltol stop at `6526b45` (r46, r4, r9, r38, r5, r7; 19232727 an empty submission) and the test files (142 of 142) | 1 | 0:00:10–1:07:15 each | 2.03 |
| 19238017–19238056 | #89 B7 | Second resume at `6526b45`: 11 preps and 22 chains, KS low-ENO refits, merge (40 jobs; running) | 1 | | |

**#89 total: about 119.7 CPU-h**: B3's 4.3, six suites' 19.7, the rethink's 48.0,
the close-out's 0.07, B4's 1.3 besides its suite, and B5's 46.4 besides its suites
(the smokes' 4.9 and the full run's 41.5).
The first two suites were missing from #93's rows and are added here from `sacct`.
**Phase 16 total with both tracks: about 42,974.5 CPU-h** of 50k, so about 7.0k
remain.
