# 16a.1: path replay and snapshots

Runs of record: jobs 17754678 (tests) and 17754679 (full cycle), both at
`1d5d129`, on dave1 (milan). Spec: `spec/phases/16-recovery.md` task 16a.1, checks
V2 and V2b.

## V2 at full scale (job 17754679)

These are cells 30001 to 30003 of the 15.7 dataset: truth from `Xoshiro(1507)`, the
published model, 6,300 handshakes each.

| Seed | Firings | 60 s latent equals `generate_dataset` | Replay bitwise at all 6,300 handshakes | Simulate (s) | Replay (s) | Replay / simulate |
|---|---|---|---|---|---|---|
| 30001 | 2,490 | yes | yes | 23.09 | 20.47 | 0.886 |
| 30002 | 2,898 | yes | yes | 24.21 | 21.85 | 0.903 |
| 30003 | 3,738 | yes | yes | 27.15 | 22.46 | 0.828 |

- **"Bitwise" covers every field a handshake can change.** That is both clocks, the
  ODE state and parameters, the jump state and parameters, every rounding carry,
  every deficit and chemostat exchange, the geometry and the census.
- **Recording is transparent on the dataset's own configuration.** The recorded run
  reproduces `generate_dataset`'s 105 saved points exactly.
- **Replay saves 10 to 17% of a cycle.** The stiff solve, not the SSA, is the cost,
  so D16.6's projection stands with at most that discount. Warm simulation here is
  23 to 27 s per cycle, against 13.7's 32.2 s. That is a different truth and a
  different node, so 32.2 s stays the conservative figure.
- **Firing counts are lower than §3 says.** They are 2,490 to 3,738 per cycle at
  this truth, against the parent §3's "about 9,900", which is at nominal
  parameters. That matters for 16a.7's particle proposals.
- **A mid-cycle snapshot is 606 KiB** and takes 3.4 ms to take and 3.1 ms to
  restore. At N = 50 particles and 105 windows that is about 0.35 s per cell per
  sweep, which is negligible against 50 × 21 s of replay.

## V2 and V2b in the suite (job 17754678)

`test/test_path_replay.jl`, 10 of 10 pass, over 300 s at a truth from
`Xoshiro(1601)`:
- Recording is transparent: a recorded and an unrecorded run are equal at every
  handshake.
- **V2:** the replay equals the run at every handshake. **Mutation:** removing one
  translation firing changes GTP. A path whose labels differ is refused.
- **V2b:** a snapshot taken at 120 s continues exactly after an intervening
  particle has been changed in every field. That particle ran the SSA for 30 s, set
  every deficit to 7, shifted every carry by 0.25, scaled the conversion factor by
  1.01 and added one to every jump count. **Mutations:** a restore carrying the
  other particle's deficits, or its conversion factor, does not continue exactly.

## Not met as written: V2b's SSA clause

V2b asks for the intervening particle to differ in its SSA state, and for dropping
the SSA's next-jump state to fail the check. Neither is possible in the driver as
built:
- **The random stream is not the driver's.** The SSA draws from JumpProcesses'
  `DEFAULT_RNG`, which is `Random.default_rng()`, the task's global generator. A
  deep copy does not copy it.
- **The cached next jump is dead at a window boundary.** Every handshake calls
  `reset_aggregated_jumps!`, which redraws it (`handshake_step!`, step 4). So
  leaving it out of a snapshot could not fail.

Replay draws nothing, so V2 and V2b hold for everything a replay uses. Per-particle
random streams are 16a.7's problem, and this needs an amendment. See the handoff.
