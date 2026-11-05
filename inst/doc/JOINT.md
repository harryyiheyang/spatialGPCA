# Joint U + SPDE decomposition

`svgpc_joint_fit()` fits, at the individual level,

```
X = 1 mu' + U A' + Z B Theta + E,   theta_j ~ N(0, sigma^2 / lambda Q(kappa)^-1)
```

where `U` (N x k) are unpenalized non-spatial factors that replace fixed plink
PCs, `B` is the sparse SPDE projection, `Z` maps individuals to locations and
`lambda` is shared by all SNPs. Only sparse mesh algebra is used: no m x m
dense matrix, no N x p or m x p storage.

```r
model <- svgpc_joint_prepare(mesh, table, kappa = scale$kappa, bed = "merged")
fit <- svgpc_joint_fit(model, k = 10, r = 10, frequencies = "merged.afreq")
fit$individual             # FID, IID, U1..Uk, SPC1..SPCr
fit$individual_orthogonal  # spatial PCs W-orthogonalized against [1, U-bar]
fit$trace                  # one row per BED scan
```

Each iteration is one BED scan. Per block, the scan forms location sums and
`X'[C, Y]` straight from the 2-bit genotype codes. For every SNP, it solves the
individual-level normal equations by block elimination with one sparse Cholesky
of `S = B'WB + lambda Q`:
`t = S^-1 B'W xbar`, `c = G^-1 (C'x - H'B'W xbar)`, `theta = t - H c`, with
`C = [1, U]`, `H = S^-1 B'W Cbar` and `G = C'C - Cbar'WB H`.

**One objective.** U and lambda minimize the same criterion, the individual-level
REML with `C = [1, U]` as the fixed design and sigma^2 profiled (divided by p):

```
V(U, lambda) = (N - q) log(L / (p (N - q))) + log|S| + log|G| - m log lambda - log|Q| - log|C'C|
```

with `L = sum_j PRSS_j` the individual-level penalized residual sum of squares,
`G = C'C - C'ZB S^-1 B'Z'C = C' Sigma^-1 C`. `-log|C'C|` makes V invariant to
`U -> U T + 1 v'`. Every term comes from the sparse Cholesky and the scan; the
tests check V against the dense formula and its U-gradient by finite differences.

Each scan evaluates V exactly at a small lambda grid (one extra sparse solve per
grid point, symbolic factorization reused): the lambda step is exact 1-D
minimization of V.

**U given lambda has a closed form (no backfitting).** Profiling (mu, A, Theta)
integrates the smooth field out: with `Omega^-1 = I - ZB S^-1 B'Z'` (Woodbury,
sparse) and `X_c` the GLS-centred genotypes, the penalized RSS is minimized by
the top k solutions of the generalized eigenproblem

```
X_c X_c' Omega^-1 U = U D,   U' Omega^-1 U = I     (Eckart-Young in the Omega^-1 metric)
```

computed by block subspace iteration with k + `oversample` (6) columns: each scan forms
`X_c' Omega^-1 Y` and `X_c (X_c' Omega^-1 Y)` block by block (two GEMMs per
block, no per-SNP sparse solves) followed by Rayleigh-Ritz. The rate depends on
the eigengap of the whitened data, not on the concurvity between U and the
field. This U ignores the O(1) term `log|G| - log|C'C|` of V, whose effect on U
is O(1/N); V itself, including that term, is what lambda minimizes and what the
trace reports. Columns of U beyond the real structure sit in the noise bulk and
keep rotating slowly.

Block Krylov acceleration (`krylov = 4`): the last four iterates `Y_j` and
their images `M Y_j` (M = X_c X_c' Omega^-1) are kept, and U is the
Rayleigh-Ritz solution over their span. Omega^-1 M is symmetric, so the
projected problem needs only `Omega^-1` applied to the stored iterates (sparse,
no scan). This turns the linear rate of subspace iteration into a Krylov rate:
in the k = 10 simulation the fixed point is reached in 17 scans instead of 91.
Lambda is held while its predicted REML gain is below the tolerance, so the
history stays valid (see Lambda steps below).

The first scan uses the intercept only and a random basis; `init = "table"`
starts the basis from the first k covariate columns (e.g. plink PCs), and
`U_init` from any N x k matrix.
The V subspaces iterate during the main loop, so one final scan (`v_scans = 1`)
at fixed U and lambda computes V, the W-principal components of `F = B Theta`.

Kappa is fixed per fit; compare `fit$reml` over 3-5 candidate kappas.
`fit$concurvity` is the condition number of G (collinearity of U with the field).

## Kernels and speed

**Fused 2-bit kernel** (`kernel = "fused"`, the default for BED input; `cpp_joint_bed_scan`). It
never forms double genotype blocks; x_ij = lut_j[code_ij] (standardized, missing = 0). It has three
OpenMP phases per block of 256 SNPs, all with dynamic scheduling, because hybrid P/E-core CPUs
stall static schedules:
- 1a, over SNPs: location sums and x'x.
- 1b, over chunks of 256 individuals: X'Y, by adding rows of Y into three genotype-class sums with
  thread-private accumulators. Each Y chunk stays in cache across all SNPs of the block.
- 2, over chunks of individuals: X T from per-SNP tables lut[code] * T_j.

Y and the inner sums are single precision, flushed to double every 256 individuals; the relative
error against the double path is about 1e-7. The `"cpp"` and `"blas"` kernels remain for custom
readers.

**Genotype cache (`cache_gb = 8`).** The selected variants' 2-bit rows (N/4 bytes each) are read
into memory once, so later scans do no I/O. A full-chromosome BED does not fit and is disk bound.
For UKB, extract a small BED with exactly the SNPs used for the genetic PCs. `snps` also accepts
a count (first variants) or a vector of variant indices.

**Compiler flags.** Build with `-O3 -march=native`, which roughly halves the kernel time. No
change to the global R configuration is needed:

```
printf 'CXXFLAGS = -O3 -march=native\nCXX17FLAGS = -O3 -march=native\n' > /tmp/Makevars.fast
R_MAKEVARS_USER=/tmp/Makevars.fast R CMD INSTALL spatialGPCA
```

On a cluster, build on the node type the jobs run on. The tall N-row products outside the kernel
use threaded C++ unless R is linked to an optimized BLAS (OpenBLAS/MKL), which is then used.

**Lambda steps and stopping.** Each scan evaluates the REML on a lambda grid and predicts the gain of
the parabola step. Lambda is held when that gain is below tol x |REML|. The REML is flat in lambda
near its optimum (on 7M_1, lambda from 2.3e8 to 3.0e8 gives the same REML), and moving it there
only restarts the Krylov history. The fit stops when dREML and the predicted gain are both below
tol x |REML| (`tol = 1e-7`, relative as in mgcv) and the U subspace change is below 1e-3.

**Warm start (`warm`, default off).** `warm = 10` first fits on every 10th variant. On 7M_1 it did not
reduce the number of full scans, because the subset optimum differs from the full one.

**Local benchmark** (7M_1: N = 487,409, first 20k SNPs, k = 10, i7-14700K on Windows, reference BLAS,
-O3 -march=native):

| version | s / scan | scans | total |
|---|---|---|---|
| first C++ version | 330 | not converged | |
| Krylov + fused kernel (3f4e32f) | 60 | 14 | 905 s |
| lambda hold (baf3429) | 50 | 11 | 639 s |
| cache, threaded tall products, dynamic schedule, float (33dafe3) | 22 | 9 | 225 s |

The 2048-SNP kernel bench (`inst/examples/joint_kernel_bench.R`) takes 13.1 s on 1 thread and
1.94 s on 28 threads. Peak RAM for the 20k run was 7.6 GB.
