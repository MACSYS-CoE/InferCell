# Phase 16 compute ledger

The cap is 50k CPU-h for the phase, and 5k for 16a (`spec/phases/16-recovery.md`
D16.6). CPU-h is wall time × cores, as reported by Slurm.

| Job | Task | What | Cores | Wall | CPU-h |
|---|---|---|---|---|---|
| 17754678 | 16a.1 | V2/V2b tests | 1 | 0:03:39 | 0.06 |
| 17754679 | 16a.1 | Full-cycle V2, 3 cells | 1 | 0:08:16 | 0.14 |

**16a total: 0.20 CPU-h.**
