# KIM-184 `time`

- job 21666252 on dave28, AMD EPYC 7543 32-Core Processor; Julia 1.10.5
- threads 1, BLAS threads 1, GC threads 1 (sweep 0)
- window 0–3150.0 s, 52 saves; Rodas5P + KLUFactorization, abstol 1.0e-10, reltol 1.0e-8; chunk 10

Medians of 5 after a warm-up. A pass is one dual solve of 10 partials.

| B | passes | primal (s) | gradient (s) | gradient / primal | per pass (s) | primal steps (acc+rej) | pass steps | pass njacs | pass nw | s per pass step |
|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 1 | 0.110 | 0.649 | 5.89 | 0.649 | 174+0 | 144+0 | 144 | 144 | 0.0045 |
| 20 | 2 | 0.109 | 1.277 | 11.76 | 0.636 | 174+0 | 144+0 | 144 | 144 | 0.0044 |
| 40 | 4 | 0.109 | 2.551 | 23.39 | 0.638 | 174+0 | 144+0 | 144 | 144 | 0.0044 |

Every repeat (seconds):

- B = 10: primal 0.111, 0.109, 0.110, 0.111, 0.110; gradient 0.650, 0.649, 0.651, 0.641, 0.646
- B = 20: primal 0.121, 0.109, 0.109, 0.108, 0.108; gradient 1.263, 1.286, 1.286, 1.277, 1.273
- B = 40: primal 0.108, 0.108, 0.109, 0.112, 0.110; gradient 2.551, 2.574, 2.544, 2.544, 2.553

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
