# KIM-184 `library`

- job 21667143 on dave1, AMD EPYC 7543 32-Core Processor; Julia 1.10.5
- threads 1, BLAS threads 1, GC threads 1 (sweep 0)
- window 0–3150.0 s, 52 saves; Rodas5P + KLUFactorization, abstol 1.0e-10, reltol 1.0e-8; chunk 10

One dual solve of 10 partials (B = 20, θ*, partials 1–10) over the whole window, saved on the 52-point grid; no misfit.

| linear solve | seconds | steps (acc+rej) | njacs | nw | s per step |
|---|---|---|---|---|---|
| DualSparseKLU (median of 5) | 0.645 | 144+0 | 144 | 144 | 0.0045 |
| LinearSolve DualLinearCache (1 run) | 554.1 | 144+0 | 144 | 144 | 3.848 |

Ratio 860×. Partials differ by 0.000e+00 of their largest, values by 0.000e+00 relative.

DualSparseKLU runs: 0.668, 0.670, 0.645, 0.627, 0.632

```
Julia Version 1.10.5
Commit 6f3fdf7b362 (2024-08-27 14:19 UTC)
Build Info:
  Official https://julialang.org/ release
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 64 × AMD EPYC 7543 32-Core Processor
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-15.0.7 (ORCJIT, znver3)
Threads: 1 default, 0 interactive, 1 GC (on 64 virtual cores)
Environment:
  LD_LIBRARY_PATH = /apps/modules/software/Julia/1.10.5/lib:/apps/modules/software/Julia/1.10.5/lib/julia:/apps/slurm/latest/lib:/apps/slurm/latest/lib/slurm:/opt/nvidia/latest/usr/lib64
  JULIA_PKG_OFFLINE = 1
```
