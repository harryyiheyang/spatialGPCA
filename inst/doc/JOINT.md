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

Each iteration is one BED scan. Per block the scan decodes and standardizes
genotypes in C++ (OpenMP over SNPs; `threads = 0` uses all cores), forms
location sums, and for every SNP solves the individual-level normal equations
by block elimination with one sparse Cholesky of `S = B'WB + lambda Q`:
`t = S^-1 B'W xbar`, `c = G^-1 (C'x - H'B'W xbar)`, `theta = t - H c`, with
`C = [1, U]`, `H = S^-1 B'W Cbar` and `G = C'C - Cbar'WB H`. The same scan
accumulates `(X - 1 mu' - Z B Theta) A` for the U step.

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

computed by block subspace iteration with k + 10 columns: each scan forms
`X_c' Omega^-1 Y` and `X_c (X_c' Omega^-1 Y)` block by block (two GEMMs per
block, no per-SNP sparse solves) followed by Rayleigh-Ritz. The rate depends on
the eigengap of the whitened data, not on the concurvity between U and the
field. This U ignores the O(1) term `log|G| - log|C'C|` of V, whose effect on U
is O(1/N); V itself, including that term, is what lambda minimizes and what the
trace reports. Columns of U beyond the real structure sit in the noise bulk and
keep rotating slowly; convergence is declared on the change of V and lambda.

Block Krylov acceleration (`krylov = 4`): the last four iterates `Y_j` and
their images `M Y_j` (M = X_c X_c' Omega^-1) are kept, and U is the
Rayleigh-Ritz solution over their span. Omega^-1 M is symmetric, so the
projected problem needs only `Omega^-1` applied to the stored iterates (sparse,
no scan). This turns the linear rate of subspace iteration into a Krylov rate:
in the k = 10 simulation the fixed point is reached in 17 scans instead of 91.
Lambda is held while its REML step is under 1%, so the history stays valid,
and the held step is applied at the end.

The first scan uses the intercept only and a random basis; `init = "table"`
starts the basis from the first k covariate columns (e.g. plink PCs), and
`U_init` from any N x k matrix.
The V subspaces iterate during the main loop, so one final scan (`v_scans = 1`)
at fixed U and lambda computes V, the W-principal components of `F = B Theta`, and U is
rotated by the eigenvectors of `A'A`.

Kappa is fixed per fit; compare `fit$reml` over 3-5 candidate kappas.
Block-coordinate descent converges slowly when U and the smooth field are
nearly collinear (`fit$concurvity` is the condition number of G); in that
regime check the trace and raise `max_iter`.

## Kernels and speed

With a PLINK BED input the default kernel is `"fused"` (`cpp_joint_bed_scan`). It works on the 2-bit
codes directly and never forms an N x block double matrix. For each SNP it adds rows of
`Y = [C, Omega^-1 Y_s]` into three genotype-class sums, which gives `X'Y`, the location sums and
`x'x` in one pass. A second pass over chunks of individuals adds `X T` from per-SNP lookup tables
`lut[code] * T_j`. Both passes use OpenMP: SNP groups x individual chunks for the first, individual
chunks for the second. At N = 200k, 2048 SNPs and 32 columns on 4 cores, this is 2x faster than
decoding plus the C++ products. The `"cpp"` and `"blas"` kernels remain for custom readers.

These loops are plain C++, so the compiler flags matter. Building with `-O3 -march=native` roughly
halves the time again (AVX2/AVX-512). This does not need a change to the global R configuration:

```
printf 'CXXFLAGS = -O3 -march=native\nCXX17FLAGS = -O3 -march=native\n' > /tmp/Makevars.fast
R_MAKEVARS_USER=/tmp/Makevars.fast R CMD INSTALL spatialGPCA
```

On a cluster, build on the same node type the jobs run on.

**Lambda steps.** Each scan predicts the REML gain of the lambda step from the grid parabola. Lambda
is held when that gain is below `tol`. The REML is flat in lambda near its optimum (on 7M_1, lambda
from 2.5e8 to 2.8e8 changes the REML by about 0.01). Without the hold, lambda wanders there and
restarts the Krylov history. Convergence requires dREML and the predicted gain both below tol x |REML| (tol = 1e-7 is relative, as in mgcv), and a U subspace change below 1e-3. On 7M_1 that ends at scan 9 instead of 11; the extra scans only moved the REML in its 8th significant digit.

**Warm start (`warm`, default off).** `warm = 10` first fits on every 10th variant. On 7M_1 it did
not reduce the number of full scans, because the subset optimum of U and lambda differs from the
full one, so it is off by default. `snps` also accepts a vector of variant indices, e.g. an
LD-pruned set, read by the fused kernel with one seek per contiguous run.

**Fused kernel layout.** Phase 1a runs over SNPs and computes the location sums and x'x. Phase 1b runs
over chunks of 256 individuals and computes X'Y with thread-private accumulators, so each Yt chunk
stays in cache across all SNPs of a block. Phase 2 runs over chunks of individuals and computes
X T. Y memory traffic is one pass per block instead of one pass per SNP, which was the bandwidth
limit behind about 4 busy cores out of 28.

**Genotype cache (`cache_gb = 8`).** On 7M_1, scans were bound by disk reads: about 45 MB/s on a
shared machine, so 20k SNPs (2.4 GB of 2-bit rows) took about 50 s per scan with only about 4 busy
cores. When the selected variants fit in `cache_gb`, their 2-bit rows (N/4 bytes each) are read into
memory once, and later scans do no I/O. For UKB, use an LD-pruned set passed as `snps` indices:
200k variants x 487k samples is about 24 GB.

**Precision and width.** In the fused kernel, Y and the inner sums are single precision: twice the
SIMD width and half the cache traffic. Each partial sum over 256 individuals is flushed into double
accumulators; relative error against the double path is about 1e-7. The location sums and x'x stay
double. The kernel uses dynamic OpenMP scheduling, because hybrid CPUs (P and E cores) stall static
schedules. With Krylov acceleration, the oversampling of the U subspace is 6 (`oversample`), down
from 10, which cuts the kernel width from 31 to 27 columns (+16 for the product) with no change in
the number of scans in simulation.

