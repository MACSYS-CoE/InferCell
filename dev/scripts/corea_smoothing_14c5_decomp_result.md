# 14c.5: which departure causes the coupled half's GTP gap

Job 17695360, at commit 2bbb02c, by `dev/scripts/corea_smoothing_14c5_decomp.jl`; analysed by `corea_smoothing_14c5_decomp_merge.jl`. A is the published model (clamp, fractional carry), B the clamp with continuous pools, C the smoothed drain with fractional carry, D the sampler model. Gaps are in particles, the variant minus A. "net" is the gap in the pool minus its carried GTP deficits.

**Continuous pools cause the gap; the smoothed drain does not.** B reproduces D, including the handshake at which the pair decouples, on all four seeds. C never moves GTP by more than one particle and never decouples. The gap builds only during runs of clipped drains, about 3 particles per 10 s, while GDP sits a few particles higher in the continuous-pool model; when the pool refills, GTP comes back higher and GDP lower by about the same amount.

```
=== seed 16184, worst save t=5880, GTP counters [:GTP_mRNA, :GTP_translat]
B: decouples at 6110; GTP gap at t=5880: 21.29 (A pool 654.7); net(pool-deficit) gap 21.29; max|gap| before decoupling 21.29
C: decouples at never; GTP gap at t=5880: -0.13 (A pool 654.7); net(pool-deficit) gap -0.13; max|gap| before decoupling 0.98
D: decouples at 6110; GTP gap at t=5880: 21.29 (A pool 654.7); net(pool-deficit) gap 21.29; max|gap| before decoupling 21.29
  t    A.gtp   D-A gtp  B-A gtp  C-A gtp  D-A net  A.def(sum) D.def(sum) A.clip D.clip  D-A gdp  D-A gmp  D-A atp
 5760     0.11   -0.112   -0.112    0.014    9.374  7.42e+03  7.42e+03   true   true     3.51    -3.40     6.79
 5770    -0.39    0.388    0.388    0.024   12.434  3.02e+03  3.01e+03   true   true     3.27    -3.66     6.76
 5780    -0.34    0.336    0.336   -0.153   15.250       602       587   true   true     2.71    -3.05     7.35
 5790  1911.43    5.810    5.810    0.036    5.810         0         0  false  false    -1.93    -3.88     1.06
 5800   227.67    6.480    6.480    0.023    6.480         0 2.49e-248  false  false    -4.34    -2.14     4.72
 5810   714.73    4.929    4.929   -0.147    4.929         0 2.86e-313  false  false    -0.66    -4.27     3.00
 5820     0.21   -0.212   -0.212   -0.106    3.012  2.26e+03  2.26e+03   true   true     2.26    -2.05     3.42
 5830    -0.06    0.062    0.062    0.504    5.785  3.84e+03  3.84e+03   true   true     2.47    -2.53     4.83
 5840     0.41   -0.412   -0.412   -0.384    8.667  4.05e+03  4.04e+03   true   true     3.47    -3.06     5.58
 5850     0.09   -0.088   -0.088   -0.037   11.487  3.91e+03   3.9e+03   true   true     2.77    -2.68     5.90
 5860    -0.22    0.216    0.216   -0.041   14.630  2.02e+03     2e+03   true   true     2.60    -2.81     6.47
 5870    -0.19    0.186    0.186    0.133   17.908  1.63e+03  1.61e+03   true   true     2.70    -2.89     6.69
 5880   654.68   21.286   21.286   -0.131   21.286         0         0  false  false   -17.93    -3.36     6.92
=== seed 16391, worst save t=3300, GTP counters [:GTP_mRNA, :GTP_translat]
B: decouples at 3315; GTP gap at t=3300: 9.97 (A pool 432.3); net(pool-deficit) gap 9.97; max|gap| before decoupling 9.97
C: decouples at never; GTP gap at t=3300: 0.00 (A pool 432.3); net(pool-deficit) gap 0.00; max|gap| before decoupling 0.05
D: decouples at 3315; GTP gap at t=3300: 9.96 (A pool 432.3); net(pool-deficit) gap 9.96; max|gap| before decoupling 9.96
  t    A.gtp   D-A gtp  B-A gtp  C-A gtp  D-A net  A.def(sum) D.def(sum) A.clip D.clip  D-A gdp  D-A gmp  D-A atp
 3180  2588.21    0.978    0.978    0.000    0.978         0         0  false  false     0.52    -1.50     2.85
 3190  3039.93    0.704    0.704    0.000    0.704         0         0  false  false     0.46    -1.16     2.74
 3200     0.10   -0.102   -0.102    0.000    0.617       348       347   true   true     1.29    -1.18     2.20
 3210     0.33   -0.330   -0.330    0.000    3.031       174       171   true   true     1.74    -1.41     2.98
 3220  3279.73    0.156    0.156    0.000    0.156         0         0  false  false     0.56    -0.72     1.53
 3230  2532.46    0.281    0.281    0.000    0.281         0         0  false  false     0.50    -0.78     1.11
 3240   349.22    0.527    0.527    0.000    0.527         0         0  false  false     0.31    -0.84     1.50
 3250   388.45    1.037    1.037    0.000    1.037         0         0  false  false     0.37    -1.40     1.98
 3260  1912.65    1.089    1.089    0.000    1.089         0         0  false  false     0.20    -1.29     2.53
 3270   111.52    0.756    0.756    0.000    0.756         0         0  false  false     0.85    -1.61     2.16
 3280     0.46   -0.459   -0.459    0.000    1.060  3.57e+03  3.57e+03   true   true     1.83    -1.37     2.41
 3290     0.35   -0.347   -0.347    0.000    4.631  5.86e+03  5.86e+03   true   true     2.35    -2.00     3.63
 3300   432.27    9.955    9.971    0.000    9.955         0         0  false  false    -7.53    -2.42     5.93
=== seed 18227, worst save t=4500, GTP counters [:GTP_mRNA, :GTP_translat]
B: decouples at 4515; GTP gap at t=4500: 18.74 (A pool 242.5); net(pool-deficit) gap 18.74; max|gap| before decoupling 18.74
C: decouples at never; GTP gap at t=4500: 0.13 (A pool 242.5); net(pool-deficit) gap 0.13; max|gap| before decoupling 0.97
D: decouples at 4515; GTP gap at t=4500: 18.74 (A pool 242.5); net(pool-deficit) gap 18.74; max|gap| before decoupling 18.74
  t    A.gtp   D-A gtp  B-A gtp  C-A gtp  D-A net  A.def(sum) D.def(sum) A.clip D.clip  D-A gdp  D-A gmp  D-A atp
 4380   764.58    1.766    1.766    0.018    1.766         0         0  false  false     1.66    -3.43     6.43
 4390  3588.53    1.540    1.540   -0.031    1.540         0         0  false  false     1.65    -3.19     5.87
 4400  1623.83    1.625    1.625    0.008    1.625         0         0  false  false     1.51    -3.14     5.56
 4410  1249.51    1.441    1.441   -0.047    1.441         0         0  false  false     1.80    -3.24     5.50
 4420  2912.25    1.686    1.686    0.010    1.686         0         0  false  false     1.63    -3.32     5.75
 4430  3926.25    1.720    1.720   -0.106    1.720         0         0  false  false     1.22    -2.94     5.55
 4440  2186.73    1.821    1.821    0.054    1.821         0         0  false  false     1.57    -3.39     5.66
 4450  1012.17    1.717    1.717   -0.027    1.717         0         0  false  false     1.44    -3.15     5.67
 4460  1679.39    1.445    1.445   -0.044    1.445         0         0  false  false     1.77    -3.21     5.36
 4470   558.65    1.060    1.060   -0.326    1.060         0         0  false  false     1.96    -3.02     5.24
 4480    -0.43    0.426    0.426    0.188    4.291   5.7e+03  5.69e+03   true   true     2.07    -2.50     6.20
 4490     0.01   -0.012   -0.012    0.020   11.325  3.88e+03  3.87e+03   true   true     3.68    -3.67     8.25
 4500   242.54   18.739   18.739    0.129   18.739         0         0  false  false   -14.38    -4.36    10.13
=== seed 19958, worst save t=5460, GTP counters [:GTP_mRNA, :GTP_translat]
B: decouples at 5477; GTP gap at t=5460: 15.09 (A pool 253.1); net(pool-deficit) gap 15.09; max|gap| before decoupling 15.12
C: decouples at never; GTP gap at t=5460: 0.16 (A pool 253.1); net(pool-deficit) gap 0.16; max|gap| before decoupling 0.95
D: decouples at 5477; GTP gap at t=5460: 15.09 (A pool 253.1); net(pool-deficit) gap 15.09; max|gap| before decoupling 15.12
  t    A.gtp   D-A gtp  B-A gtp  C-A gtp  D-A net  A.def(sum) D.def(sum) A.clip D.clip  D-A gdp  D-A gmp  D-A atp
 5340   544.73    1.531    1.531    0.061    1.531         0         0  false  false     0.35    -1.89     4.08
 5350  5224.47    1.240    1.239    0.070    1.240         0         0  false  false     1.12    -2.36     3.50
 5360   557.66    1.559    1.559    0.008    1.559         0         0  false  false     0.45    -2.01     4.14
 5370  4560.00    1.403    1.403   -0.023    1.403         0         0  false  false     0.94    -2.34     3.61
 5380  2013.35    1.618    1.618   -0.033    1.618         0         0  false  false     0.65    -2.26     4.14
 5390  2654.35    1.925    1.925   -0.004    1.925         0         0  false  false     0.43    -2.36     3.98
 5400     0.28   -0.284   -0.284   -0.735    2.052       691       689   true   true     3.04    -2.76     3.94
 5410  1155.99    3.806    3.806    0.156    3.806         0         0  false  false    -0.05    -3.75     4.20
 5420   961.16    1.591    1.591   -0.035    1.591         0         0  false  false     0.97    -2.56     4.30
 5430    57.43    1.034    1.034   -0.058    1.034         0         0  false  false     1.41    -2.44     4.29
 5440     0.25   -0.248   -0.248   -0.509    3.681   1.6e+03  1.59e+03   true   true     3.23    -2.98     5.13
 5450     0.21   -0.209   -0.209   -0.306    8.192  2.58e+03  2.58e+03   true   true     3.52    -3.32     6.10
 5460   253.06   15.092   15.092    0.164   15.092         0         0  false  false   -11.27    -3.82     8.42
```
