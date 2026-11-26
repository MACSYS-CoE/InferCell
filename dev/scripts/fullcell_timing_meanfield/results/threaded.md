# KIM-184 `threaded`

- job 21666253 on dave64, AMD EPYC 7543 32-Core Processor; Julia 1.10.5
- threads 4, BLAS threads 1, GC threads 5 (sweep 1)
- window 0–3150.0 s, 52 saves; Rodas5P + KLUFactorization, abstol 1.0e-10, reltol 1.0e-8; chunk 10

Medians of 5 after a warm-up. A pass is one dual solve of 10 partials.

| B | passes | primal (s) | gradient (s) | gradient / primal | per pass (s) | primal steps (acc+rej) | pass steps | pass njacs | pass nw | s per pass step |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 1 | 0.107 | 0.581 | 5.41 | 0.581 | 174+0 | 144+0 | 144 | 144 | 0.0040 |
| 20 | 2 | 0.101 | 0.620 | 6.12 | 0.603 | 174+0 | 144+0 | 144 | 144 | 0.0042 |
| 40 | 4 | 0.102 | 0.760 | 7.43 | 0.680 | 174+0 | 144+0 | 144 | 144 | 0.0047 |

Every repeat (seconds):

- B = 10: primal 0.107, 0.109, 0.107, 0.108, 0.106; gradient 0.604, 0.624, 0.581, 0.568, 0.579
- B = 20: primal 0.103, 0.102, 0.100, 0.100, 0.101; gradient 0.620, 0.606, 0.622, 0.620, 0.634
- B = 40: primal 0.102, 0.101, 0.103, 0.102, 0.102; gradient 0.961, 0.806, 0.760, 0.678, 0.728

Threaded gradient vs serial (first repeat), max |Δ| per B: 0.000e+00, 0.000e+00, 0.000e+00

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
Threads: 4 default, 0 interactive, 5 GC (on 64 virtual cores)
Environment:
  LD_LIBRARY_PATH = /apps/modules/software/Julia/1.10.5/lib:/apps/modules/software/Julia/1.10.5/lib/julia:/apps/slurm/latest/lib:/apps/slurm/latest/lib/slurm:/opt/nvidia/latest/usr/lib64
  JULIA_PKG_OFFLINE = 1
```
