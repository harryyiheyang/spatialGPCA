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
accumulates the closed-form U step `(X - 1 mu' - Z B Theta) A (A'A)^-1`.

Lambda uses the location-level pooled REML given U,
`(n - k - 1) log sigma2(lambda) + log|S| + log|G_loc| - m log lambda - log|Q|`,
evaluated exactly in the same scan at a small grid of lambda values (one extra
sparse solve per grid point; the symbolic factorization is reused). The grid
is recentred and narrowed every scan by a parabola step.

Initialization is spatial first: scans with intercept only estimate lambda and
run subspace iteration on `X - 1 mu' - Z B Theta`, whose top k Ritz vectors
start U. `init = "table"` starts from the first k covariate columns of the
table (the two-stage solution) for comparison; compare `fit$reml`.
After convergence, `v_scans` scans at fixed U and lambda compute V, the
W-principal components of `F = B Theta`, by subspace iteration, and U is
rotated by the eigenvectors of `A'A`.

Kappa is fixed per fit; compare `fit$reml` over 3-5 candidate kappas.
Block-coordinate descent converges slowly when U and the smooth field are
nearly collinear (`fit$concurvity` is the condition number of G); in that
regime check the trace and raise `max_iter`.
