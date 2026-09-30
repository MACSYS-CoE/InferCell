# 16a.1: path replay and snapshots

Runs of record: jobs 17755845 (tests) and 17755846 (full cycle), both at
`9bdc3d4`, on dave1 (milan). They supersede the first runs, jobs 17754678 and
17754679 at `1d5d129`, which gave the same V2 results but checked V2b at one
boundary only. Spec: `spec/phases/16-recovery.md` task 16a.1, checks
V2 and V2b.

## V2 at full scale (job 17755846)

These are cells 30001 to 30003 of the 15.7 dataset: truth from `Xoshiro(1507)`, the
published model, 6,300 handshakes each.

| Seed | Firings | 60 s latent equals `generate_dataset` | Replay bitwise at all 6,300 handshakes | Simulate (s) | Replay (s) | Replay / simulate |
|---|---|---|---|---|---|---|
| 30001 | 2,490 | yes | yes | 23.95 | 21.14 | 0.883 |
| 30002 | 2,898 | yes | yes | 23.77 | 21.46 | 0.903 |
| 30003 | 3,738 | yes | yes | 25.54 | 22.71 | 0.889 |

- **"Bitwise" covers every field a handshake can change.** That is both clocks, the
  ODE state and parameters, the jump state and parameters, every rounding carry,
  every deficit and chemostat exchange, the geometry and the census.
- **V2b at full scale.** On cell 30001, a snapshot restored at each of the 104
  window boundaries reproduces the next window bitwise.
- **Recording is transparent on the dataset's own configuration.** The recorded run
  reproduces `generate_dataset`'s 105 saved points exactly.
- **Replay saves 10 to 12% of a cycle.** The stiff solve, not the SSA, is the cost,
  so D16.6's projection stands with at most that discount. Warm simulation here is
  24 to 26 s per cycle, against 13.7's 32.2 s. That is a different truth and a
  different node, so 32.2 s stays the conservative figure.
- **Firing counts are lower than §3 says.** They are 2,490 to 3,738 per cycle at
  this truth, against the parent §3's "about 9,900", which is at nominal
  parameters. That matters for 16a.7's particle proposals.
- **A mid-cycle snapshot is 606 KiB** and takes 3.4 to 4.1 ms to take and 3.2 to
  3.9 ms to restore. At N = 50 particles and 105 windows that is about 0.4 s per
  cell per sweep, which is negligible against 50 × 21 s of replay.

## V2 and V2b in the suite (job 17755845)

`test/test_path_replay.jl`, 16 of 16 pass, over 300 s at a truth from
`Xoshiro(1601)`:
- Recording is transparent: a recorded and an unrecorded run are equal at every
  handshake.
- **V2:** the replay equals the run at every handshake. **Mutation:** removing one
  translation firing changes GTP. A path whose labels differ is refused.
- **V2b:**
  - A snapshot taken at 120 s continues exactly after an intervening particle has
    been changed in every field. That particle ran the SSA for 30 s, set every
    deficit to 7, shifted every carry by 0.25, scaled the conversion factor by 1.01
    and added one to every jump count.
  - At every window boundary (60 to 240 s), a restore reproduces the next window.
  - **Mutations:** a restore carrying the other particle's deficits, or its
    conversion factor, does not continue exactly.
- **The SSA's cached next jump is dead at a boundary.** Corrupting it in a restored
  particle changes nothing when the global stream is seeded alike, and a different
  seed gives a different forward run.

## V2b's SSA clause, amended

V2b asked for a snapshot to carry the SSA's random stream and its pre-drawn next
jump. Neither is possible:
- The SSA draws from JumpProcesses' `DEFAULT_RNG`, which is
  `Random.default_rng()`, the task's global generator. A deep copy does not copy
  it.
- Every handshake calls `reset_aggregated_jumps!`, which redraws the next jump.

The clause is replaced, and per-particle streams move to 16a.7
(`spec/phases/16-recovery.md` §12 2026-09-30, approved).
