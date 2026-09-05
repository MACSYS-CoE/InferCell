# The 1 s handshake

Source: [`src/handshake.jl`](https://github.com/MACSYS-CoE/InferCell/blob/main/src/handshake.jl).

A mixed ODE/jump composition is advanced by a **split-operator exchange**, not by a single solver call. Every simulated second, counts from the stochastic block become rate-law parameters of the metabolic block, the metabolic block integrates one second, and the costs the stochastic block accrued are debited back against its pools. This mirrors the published whole-cell model's `hookSimulation`, whose communication interval is `delt = 1.0` s.

The two blocks keep **separate state vectors** — the ODE block in mM over `Float64`, the jump block in particles over `Int` — which is what makes the conversion below a real boundary rather than a units comment.

## `HandshakeDriver`

Returned by [`build_problem`](orchestrator.md) for a mixed composition. It holds two live integrators and the exchange records lowered from the declared [coupling edges](corea-interface.md).

```julia
driver = build_problem([metabolism, expression]; tspan = (0.0, 600.0))
record = run_handshake!(driver, 600)
record.census        # (; handshakes, clipped, fraction, deficits)
```

`run_handshake!` refuses to run past the `tspan` the driver was built with — `n_steps * interval` must fit inside it — and returns the handshake times, the ODE state in mM and the jump state in particles at each of them — per handshake rather than per solver step, because the handshake is the only instant at which the two blocks agree on a state. `handshake_step!` runs one exchange.

### Policy keywords

| keyword | default | what it fixes |
|---|---|---|
| `interval` | `1.0` | the exchange period, in simulated seconds |
| `rounding` | `:fractional_carry` | how a continuous pool is written back as whole particles |
| `radius_nm` | `200.0` | the cell radius behind the count↔concentration conversion |
| `ode_solver` | `Rodas5P()` | the stiff integrator for the metabolic block |
| `abstol`, `reltol` | `1e-10`, `1e-8` | pinned tolerances |

## Counts and concentrations

`corea_particles_per_mM(radius_nm)` derives the conversion factor from Avogadro's constant and the volume of a sphere, rather than carrying a transcribed constant: at the published 200 nm radius it returns 20180.39 particles per mM, the 20,180 the scoping note records. `counts_to_mM(n, factor)` is exact division; the other direction needs a policy, because particles are whole.

### Rounding policies

Writing a pool back as whole particles once per handshake, ~6,300 times per cell cycle, gives the three policies three different signatures:

- `:fractional_carry` — the remainder is carried to the next handshake, so the emitted counts track the exact running total and the residual **does not accumulate**. The default.
- `:stochastic` — round up with probability equal to the fraction. Unbiased, so the residual accumulates as the **square root** of the handshake count.
- `:deterministic` — round to nearest, discarding the remainder. Biased, so the residual accumulates **linearly**, reaching up to 0.156 mM per species over a cycle — about 140× the integrator's own bound, and measured at 0.125 mM on the phase-3 toy. Implemented only so that this signature can be demonstrated; it is not a policy this project runs under.

The policy is a field on the driver rather than an argument to the arithmetic, so a trajectory carries the policy that produced it.

## Deferred debits

A [`DeferredCounterEdge`](corea-interface.md) declares that the stochastic block accrues a cost in a named counter and the hook debits it against a pool one step later. Under the published `:clamped_deficit_carried` policy the pool floors at zero and the shortfall is carried to the next handshake — the interface's `max(0, ·)`, and the reason a clamped counter is reported as a gradient obstruction. The two labelled departures also execute: `:unclamped` lets the pool go negative, and `:smoothed` replaces the floor with a softplus of the declared width, which converges to the clamped answer as that width narrows.

One counter may feed several pools — the published charged-tRNA transfer debits the charged pool and credits the uncharged one from a single accrual. The counter is read and cleared once per handshake, so every pool it feeds sees the same accrued value, and a credit matches what the debited pool actually paid rather than what was asked of it.

`clipping_census(driver)` reports how many handshakes ran and at how many of them a counter carried a shortfall material against its cost. The comparison is relative rather than against zero: under `:smoothed` the softplus leaves a strictly positive residue at every step however deep the pool, so an absolute test would report every handshake as clipping. A non-zero count at published parameters means the non-smooth drain is in the operating regime.

## What this layer does not yet do

The 60 s rate-constant rebuild and the volume chain are declared by [`RateConstantEdge`](corea-interface.md) and [`VolumeEdge`](corea-interface.md) and are **not** executed here; they hook into the same loop in later work. Nor does the driver plug into the inference entry points, which build a single SciML problem and call `remake` on it.
