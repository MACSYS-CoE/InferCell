# KIM-99: route (ii)'s log-density and gradient on the reduced model, 6× faster

**Runs:** `dev/scripts/route_ii_speed_kim99.sh`, each from a worktree at a pushed
commit, on milan (AMD EPYC 7543), 1 thread. The setup is KIM-92's pilot profile's:
`mixed_setup` on `d11_models()` over `COREA_TARGETS`, δ = S = 0, default tolerances,
on the pilot's dataset (200 cells, 52 saves to 3120 s, 21 coordinates). The dataset's
recipe is in the script's header.

| what | job | commit | output |
|---|---|---|---|
| before (KIM-92's profile, main's `src/`) | 20454167 | a1be360 | `/fred/oz022/tkimpson/infercell-kim92-pilot-2026-11-10/profile.{md,jls}` |
| after: split, rhs, variants, grad | 20528147–51 | 7b2e209 | `/fred/oz022/tkimpson/infercell-kim99-after-2026-11-10` |
| steps 1–2 only: split | 20529297–300 | 74f2a36 | `/fred/oz022/tkimpson/infercell-kim99-nocache-2026-11-10` |
| smoke profile at main's `src/` | 20526060–63 | 0e618c7 | `/fred/oz022/tkimpson/infercell-kim99-smoke2-2026-11-10` |
| fixture | 20525536 | 95b449f | `test/fixtures/route_ii_density_kim99.jl` |

## Result: the gradient is 6.0× faster at the truth and 6.4× at the prior median

That is inside KIM-99's 5–8× target, above KIM-93's 3× fallback line. Times are
medians of 3.

| x | call | before (20454167) | after (20528151) | speed-up |
|---|---|---|---|---|
| truth | gradient | 277.0 s, 820 GiB | 45.9 s, 58 GiB | 6.0× |
| prior median | gradient | 222.4 s, 679 GiB | 35.0 s, 52 GiB | 6.4× |
| truth | log-density | 79.8 s, 203 GiB | 2.0 s, 2.9 GiB | 40× |
| prior median | log-density | 60.4 s, 144 GiB | 1.6 s, 2.6 GiB | 38× |

**The density is unchanged.** At both x, the log-density and all 21 gradient
components are bitwise equal before and after (every |Δ| is 0.0), and the solver
takes the same steps (2,541 + 883 at the truth, 1,749 + 351 at the prior median).
The RHS calls per log-density fell from 245,920 to 50,263 at the truth.

**The repo test** is `test/test_corea_mixed.jl`'s fixture testset, on
`corea_case(120.0)` at four x. The fixture was frozen at 95b449f, whose `src/` is
main's. The bound is 1e-10 on the fixture's CPU and Julia (EPYC 7543, 1.10.5) and
1e-8 elsewhere: GitHub's runner moved the four values by 1.3e-10 to 3.1e-10 (run
42789293994), about 3e-12 relative.

## The three changes

1. **The Jacobian ran at ForwardDiff chunk 1** (b7da80d): 85 RHS calls per Jacobian.
   The problem is now `FullSpecialize` and the default solver,
   `mean_field_solver(route)`, sets the chunk explicitly (11 here). Each partial is
   computed alone at any chunk, so the Jacobian is the same.
2. **MeanField's kernel was not inferred** (74f2a36): 782 KB and about 110 µs per
   RHS call. Unrolled by an `@generated` function in the same order: 448 B, 4 µs.
3. **The transcript factor redid shared work for each of the 200 cells** (6f70036).
   Each window's distribution is now computed once per (gene, window, starting
   count, cap) and shared; the sum over cells keeps its order.

One dual pass (11 partials) at the truth:

| code | solve | transcript factor | pass |
|---|---|---|---|
| steps 1–2 (74f2a36, job 20529298) | 11.4 s | 50.0 s | 59.1 s |
| steps 1–3 (7b2e209, job 20528148) | 12.2 s | 6.6 s | 19.6 s |

Steps 1–2 alone give about 2× on the gradient, below the fallback line; the cache
is what reaches 6×.

**The Jacobian has to stay exact.** Finite differences (variants stage, job
20528150) move ℓ by 4.07e-10 at the truth and by 4.6e-7 at the prior median, where
they also change the steps (1,785 + 317 against 1,749 + 351). Every ForwardDiff
variant gave Δ = 0.

**A plain `Rodas5P()`** falls back to chunk 1 under a gradient: about 5× the RHS
calls per dual solve, 1.2–1.3× the time (14.78 s against 11.41 s at the truth,
9.19 s against 7.73 s at the prior median).

## What this does not settle

- **The dual solve's split is unmeasured.** It is 12 s per pass. The RHS at
  single-level duals takes 19 µs per call, about 0.8 s over the solve's calls, but
  the Jacobian's calls run on nested duals and were not timed, so how the rest
  divides between them and the dense dual LU is not known. A sparse Jacobian or a
  different linear solver would change the discretisation and break the 1e-10
  criterion, so either needs a Decision-log entry first.
- **The transcript factor's 6.6 s per pass** allocates 28 GiB; in-place
  uniformisation buffers are the next lever.
- **KIM-92's chains at 45.9 s per gradient:** a 200 + 200 chain is about 3.2 days at
  tree depth 4 and 6.6 days at depth 5, close to milan's 7-day limit; depth 6 does
  not fit.

## Tests

At 7b2e209: full suite 29,908/29,908 (job 20528153); test_corea_mixed,
test_mean_field_route and test_mixed_likelihood 454/454 (job 20528152).
