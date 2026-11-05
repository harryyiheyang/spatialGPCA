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
