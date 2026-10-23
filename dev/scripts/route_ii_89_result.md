# Issue #89 B4: route (ii) against the hybrid at a pinned mean-field path

**Built:** `build_mean_field_problem` (`src/mean_field_route.jl`), the reduced syn3A
model with gene expression swapped for its mean field, as one ODE. It executes the
catalytic, drain, rebuild and growth edges (§12 2026-10-23, "B4 as built").
**Suite:** 28709 of 28709 at `67c9eef` (job 19042448). It includes AV3 on route (ii):
ForwardDiff through the solve and its callbacks agrees with central differences to
1e-6 relative on M0's panel.

**The check, as written: FAIL at 0 of 5 truths.** **Amended to a check on
convergence (§12 2026-10-23, after the run): pass at 3 of 5.** corea1 and corea2 fail
the amended max |Δ| on 3PG, from a coarse pair of steps. The finer pair, which
settles it on M0, did not run there.

## Runs

All runs used commit `5a56836` and are archived in
`/fred/oz022/tkimpson/infercell-89-b4-archive-2026-10-23/`.

| What | Job |
|---|---|
| Precompile | 19055277 |
| Route (ii) at m0_v6, m0_low, m0_high, corea1, corea2 | 19055278, 19055282, 19055286, 19055290, 19055294 |
| B3's `meanfield` prototype at 1, 0.5, 0.25 s steps | the three ids after each route job |
| Merge (`merged.md`) | 19055298 |
| Extrapolation (`richardson.md`, at `feb5a75`) | 19057735 |

The first submission (precompile 19042449) ran out of memory on the default
allocation, and its stages never started. Its logs are in
`logs_failed_precompile_19042449/`.

The prototype clipped no drains at any truth. corea1 and corea2's 0.25 s prototypes
(19055293, 19055297) stopped on MaxIters in the hybrid's metabolic block at about
t = 6040 s of 6300.

## The rule as written (job 19055298)

The rule: KG(route, prototype at 1 s) ≤ 0.1 × B3's uncorrected gap total, and
max |Δ log x̄| ≤ 0.02, with σ_b = 10%.

| Truth | KG | bound | max \|Δ\| at 1 s (pool) | 0.5 s | 0.25 s | verdict |
|---|---|---|---|---|---|---|
| m0_v6 | 3.49 | 7.08 | 0.0833 (Pi) | 0.0423 | 0.0210 | FAIL |
| m0_low | 3.88 | 1.24 | 0.0801 (Pi) | 0.0412 | 0.0204 | FAIL |
| m0_high | 2.72 | 2.58 | 0.0632 (Pi) | 0.0315 | 0.0161 | FAIL |
| corea1 | 35.5 | 9.99 | 0.0771 (Pi) | 0.0392 | — | FAIL |
| corea2 | 132 | 20.2 | 0.1130 (GTP) | 0.0534 | — | FAIL |

**Diagnosis.** Every pool's difference shrinks as the prototype's step shrinks, and
the largest ones halve: Pi goes 0.083 → 0.042 → 0.021 at m0_v6. That is a
first-order splitting error in the prototype's 1 s exchange, and route (ii) is its
Δt → 0 limit. The pools it moves most are:

- Pi, mean Δ −0.075;
- 3PG and 2PG, +0.028 and +0.019;
- GDP and GTP.

Per-pool tables are in `merged.md`.

## The amended rule (job 19057735)

Route (ii) is compared with 2·proto(h/2) − proto(h), at the finest pair of steps
that ran. The bounds are the same.

| Truth | pair (s) | KG | bound | max \|Δ\| (pool) | verdict | 3PG's Δ at 1 / 0.5 / 0.25 s |
|---|---|---|---|---|---|---|
| m0_v6 | 0.5, 0.25 | 0.0236 | 7.08 | 0.0085 (Pi) | pass | 0.0316 / 0.0086 / 0.0026 |
| m0_low | 0.5, 0.25 | 0.0138 | 1.24 | 0.0048 (3PG) | pass | 0.0325 / 0.0102 / 0.0031 |
| m0_high | 0.5, 0.25 | 0.0144 | 2.58 | 0.0052 (3PG) | pass | 0.0311 / 0.0091 / 0.0028 |
| corea1 | 1, 0.5 | 3.04 | 9.99 | 0.0211 (3PG) | **FAIL** | 0.0419 / 0.0116 / — |
| corea2 | 1, 0.5 | 6.16 | 20.2 | 0.0347 (3PG) | **FAIL** | 0.0840 / 0.0262 / — |

3PG converges faster than first order: about 3.5× per halving. So extrapolating
from the 1/0.5 s pair overshoots it. On M0 that pair leaves 0.015–0.019 on 3PG,
inside 0.02, and the 0.5/0.25 s pair leaves at most 0.005. The corea truths have
only the coarse pair, so their failure is not settled either way. Settling it needs
the 0.25 s prototype to finish at corea1 and corea2.

## KC's baseline: one route (ii) solve

Route (ii) uses Rodas5P with abstol 1e-10 and reltol 1e-8, and saves every 60 s.
Times are the median of 3 solves.

| Truth | horizon | s/solve | MiB allocated | growth |
|---|---|---|---|---|
| m0_v6 | 600 s | 1.11 | 1487 | 1.0016× |
| m0_low | 600 s | 1.01 | 1562 | 1.0025× |
| m0_high | 600 s | 1.19 | 1630 | 1.0061× |
| corea1 | 6300 s | 9.02 | 17134 | 1.0169× |
| corea2 | 6300 s | 11.55 | 20184 | 1.0119× |

`MeanField.dynamics` is known to allocate. It was not profiled here; speeding it up was out of B4's scope.
