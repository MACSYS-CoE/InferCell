# KIM-195 `clamp`

- job 21842560 on gina8, AMD EPYC 7543 32-Core Processor; Julia 1.11.9; LinearSolve 3.87.0
- InferCell 491cb900ed25f6b851af8aa2d12aa7112fab8261; syn3a-ensemble 2075809cbb825140cb31a7bec8bc3529770c8e13 (data/model unmodified)
- threads 1, BLAS threads 1, GC threads 1
- window 0–3150.0 s in minute segments, saves at the whole minutes; Rodas5P, KLUFactorization (primal) and DualSparseKLU (dual, unless a row says the library path); abstol 1.0e-10, reltol 1.0e-8; chunk 10

SpeedyCell: procid 1, 8 seeds, the mean at each minute's end (clamp/speedycell.csv; host gina8 job 21842560 seeds 1:8 minutes 54 wall 199.2 s
syn3a-ensemble 2075809cbb825140cb31a7bec8bc3529770c8e13).
Mean field at θ = 0: steps over the window, κ·c·p/(p + 1) 16487, κ·c·p/(c + p + 1) 6415.

| minute | quantity | SpeedyCell mean | SpeedyCell min at a hook | κ·c·p/(p + 1) | κ·c·p/(c + p + 1) |
|---|---|---|---|---|---|
| 1 | M_datp_c | 7019 | 742 | 6839 | 6840 |
| 1 | M_dttp_c | 8288 | 557 | 8076 | 8077 |
| 1 | M_dctp_c | 8783 | 375 | 8782 | 8782 |
| 1 | M_dgtp_c | 2337 | 318 | 2369 | 2369 |
| 1 | dATP_DNArep | 0 |  | 29.15 | 29.26 |
| 1 | dTTP_DNArep | 0 |  | 29.15 | 29.24 |
| 1 | dCTP_DNArep | 0 |  | 8.546 | 8.554 |
| 1 | dGTP_DNArep | 0 |  | 8.549 | 8.577 |
| 5 | M_datp_c | 2.509e+04 | 14854 | 0.03253 | 93.5 |
| 5 | M_dttp_c | 2.861e+04 | 17431 | 2152 | 2060 |
| 5 | M_dctp_c | 2.374e+04 | 20039 | 1.463e+04 | 1.459e+04 |
| 5 | M_dgtp_c | 1.166e+04 | 8370 | 2455 | 2425 |
| 5 | dATP_DNArep | 0 |  | 2883 | 3074 |
| 5 | dTTP_DNArep | 0 |  | 138.9 | 150.8 |
| 5 | dCTP_DNArep | 0 |  | 46.33 | 47.05 |
| 5 | dGTP_DNArep | 0 |  | 46.34 | 47.84 |
| 10 | M_datp_c | 1.883e+04 | 0 | 0.1685 | 94.49 |
| 10 | M_dttp_c | 2.753e+04 | 669 | 6297 | 117.1 |
| 10 | M_dctp_c | 3.441e+04 | 25253 | 2.735e+04 | 2.37e+04 |
| 10 | M_dgtp_c | 1.339e+04 | 4829 | 8178 | 3675 |
| 10 | dATP_DNArep | 1364 |  | 675.7 | 1.322e+04 |
| 10 | dTTP_DNArep | 551.5 |  | 147.7 | 4426 |
| 10 | dCTP_DNArep | 167.8 |  | 47.82 | 43.7 |
| 10 | dGTP_DNArep | 167.8 |  | 47.83 | 44.15 |
| 20 | M_datp_c | 4165 | 0 | 0.1074 | 79.93 |
| 20 | M_dttp_c | 2.818e+04 | 10997 | 2.443e+04 | 134.8 |
| 20 | M_dctp_c | 5.693e+04 | 47392 | 5.722e+04 | 4.718e+04 |
| 20 | M_dgtp_c | 1.343e+04 | 6798 | 1.529e+04 | 1975 |
| 20 | dATP_DNArep | 1593 |  | 876.9 | 3.752e+04 |
| 20 | dTTP_DNArep | 261 |  | 148 | 4989 |
| 20 | dCTP_DNArep | 76.25 |  | 49.48 | 40.8 |
| 20 | dGTP_DNArep | 76.25 |  | 49.48 | 41.65 |
| 30 | M_datp_c | 1104 | 0 | 3098 | 82.95 |
| 30 | M_dttp_c | 4.773e+04 | 20582 | 5.84e+04 | 3180 |
| 30 | M_dctp_c | 8.747e+04 | 75008 | 9.597e+04 | 7.691e+04 |
| 30 | M_dgtp_c | 1.681e+04 | 8279 | 2.334e+04 | 41.02 |
| 30 | dATP_DNArep | 804 |  | 13.1 | 6.356e+04 |
| 30 | dTTP_DNArep | 75.38 |  | 13.09 | 130.6 |
| 30 | dCTP_DNArep | 29.38 |  | 4.366 | 39.49 |
| 30 | dGTP_DNArep | 29.38 |  | 4.366 | 475.8 |
| 40 | M_datp_c | 1352 | 0 | 3052 | 86.28 |
| 40 | M_dttp_c | 8.108e+04 | 32369 | 9.264e+04 | 1.883e+04 |
| 40 | M_dctp_c | 1.245e+05 | 103722 | 1.373e+05 | 1.134e+05 |
| 40 | M_dgtp_c | 2.117e+04 | 8849 | 2.942e+04 | 41.72 |
| 40 | dATP_DNArep | 1340 |  | 12.12 | 8.745e+04 |
| 40 | dTTP_DNArep | 0 |  | 12.12 | 124.4 |
| 40 | dCTP_DNArep | 0 |  | 3.841 | 39.46 |
| 40 | dGTP_DNArep | 0 |  | 3.842 | 666.2 |
| 52 | M_datp_c | 2532 | 0 | 3996 | 88.65 |
| 52 | M_dttp_c | 1.314e+05 | 63574 | 1.455e+05 | 4.699e+04 |
| 52 | M_dctp_c | 1.775e+05 | 146005 | 1.958e+05 | 1.641e+05 |
| 52 | M_dgtp_c | 3.039e+04 | 12856 | 3.901e+04 | 198 |
| 52 | dATP_DNArep | 1247 |  | 10.79 | 1.167e+05 |
| 52 | dTTP_DNArep | 78.5 |  | 10.78 | 148.2 |
| 52 | dCTP_DNArep | 27.75 |  | 3.525 | 45.07 |
| 52 | dGTP_DNArep | 27.75 |  | 3.525 | 58.09 |
