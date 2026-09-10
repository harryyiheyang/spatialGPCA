# SPDE reference

Start with the [README workflow](../../README.md#2-construct-the-spde-mesh-and-inspect-its-range).
The `svgpc_spde_` functions return ordinary lists; assign accumulation results
and use `saveRDS/readRDS` for storage.

## Model and units

Use projected planar kilometres, as in `spatial$locations[, c("x", "y")]`.
The model uses linear triangular elements, lumped mass D, stiffness K, and

```
Q = tau^2 (kappa^4 D + 2 kappa^2 K + K D^-1 K).
```

This is the same FEM algebra as mgcvST's `.spde_basis_fem` and
`spde_precision`, evaluated directly in km instead of mgcvST's normalized mesh
coordinates. With normalization width L km, `kappa_normalized=L*kappa_km`.
For identical coefficient precision under that change of units,
`tau_normalized=tau_km/L`; otherwise the common lambda absorbs the overall
precision factor. Keep tau fixed when comparing lambda across runs. There is
no automatic intercept projection in the raw basis, no Q inverse, and no ridge.

`n` is the number of distinct locations, `G` is location by SNP and `counts`
gives location population sizes. `svgpc_spde_locations` prepares locations directly
from the individual table without GP clustering. `prepare_locations` uses the exact
individual-to-location map and averages all supplied numeric covariates. An
intercept is added. SNP QC remains external to the package.

`conditioning="profile"` profiles SNP-specific fixed effects.
`conditioning="qr"` imposes `Xq' sqrt(W) B theta=0` as well as profiling.
They are explicitly different models; qr is the default to match the current
public GP API. Neither method silently substitutes for the other.

## Sparse computation and dense costs

Write `Z=sqrt(W)B`, `H=Z'Xq`, `Y=sqrt(W)G` and `S=B'WB-HH'`.
Each input block computes

```
U = B' W G - H (Xq' Y)
E += ||Y-Xq(Xq'Y)||_F^2
C0 += U U'
```

The first large projection uses the three-nonzero rows of B. We do not first
form `B N`, dense T, or a complete QR complement. Sparse permuted Cholesky
`P Q P'=LL'` whitens the mesh-space S. For qr, a thin QR of `L^-1 P H`
projects the whitened S on both sides in mesh space. Its positive eigenvectors
give `V` satisfying `V'QV=I`, `V'SV=diag(d)` and the chosen constraints.
With `VD=V diag(d^-1/2)`, all-SNP statistics are

```
KJ = VD' C0 VD
a = diag(KJ), b = E-sum(a).
```

In `accumulation="moments"`, this moves the dense eigenvector transformation
outside the SNP loop. Rows of U
for vertices with identically zero observation-basis columns are exactly zero:
only the active rows are accumulated and the matching rows of VD are used in
the final transformation. The full B, full mesh Q, whitening and coefficient
dimension are retained. This does not discard unobserved vertices from the prior.
Accumulation also diagonally equilibrates each active sparse column by its
weighted norm and applies the compensating row scale to VD. The original B
and Q are unchanged. Very weakly identified directions can amplify rounding
when transforming an already accumulated covariance. The default `auto`
uses `directions` when `max(d)/min(d)>1e8`: it computes `J=VD'U` on each
sparse block and then accumulates `J J'`. This costs extra dense active-mesh
by rank multiplication per SNP block but avoids squaring that conditioning
problem in the covariance transformation. The model reports `gram_condition`
and the actual accumulation method. Both algorithms implement the same
full-SNP statistics; neither changes the spectrum, penalty or constraints.
If a is the number of active vertices, covariance accumulation still costs
O(a^2 p), and preparation/whitening,
exact covariance transformation and exact PCA remain dense O(m^3) work with
O(m^2) storage. Active-vertex covariance can still be more costly than
accumulation in estimable rank t when constraints or rank deficiency make t small.
The method is not an arbitrarily large mesh solver or a promise of universal
end-to-end acceleration. It retains O(n q) covariate storage and O(n block)
genotype working memory; no full fitted n by p matrix is needed for PCA.

The REML objective is the existing all-SNP, shared-variance profile criterion:
`sum(log(1+d/lambda)) + r*log((b+sum(a*lambda/(d+lambda)))/E)`, with analytic
boundary treatment and `r=n-rank(X)`. It is not per-SNP REML. GCV remains a
separate residual-space criterion. Raw-F PCA uses `Phi=B VD`,
`R'R=Phi'Phi`, and eigenvectors of `R diag(h) KJ diag(h) R'/p`.
The returned PCs are `B VD R^-1 U`. H2, H and 2H-H2 choose SNP loadings using
the same contract weights, then map them through the same raw fitted F.

## Scale calibration and neighbourhoods

GP geographic cluster centres and SPDE finite-mesh (FM) vertices are two
different knot-generation paths. Neither defines the observation locations:
aggregating individuals at exactly shared coordinates is a separate data layer.
For SPDE, the actual edges of `mesh$tv` define a vertex graph: first ring is
one edge hop; second ring is exactly two hops. This is vertex adjacency,
not triangle-element adjacency. Refining a mesh changes the physical size
of a ring. The audit therefore reports distances in kilometres as well.

For alpha=2 in two dimensions, nu=1 and nominal infinite-plane correlation
is `(kappa*d) K1(kappa*d)`. `svgpc_spde_scale(distance, correlation)` solves
the desired threshold directly. At `d=rho`:

| kappa | nominal correlation |
|---|---:|
| 1.6580762/rho | exp(-1) = 0.367879 |
| 2.5/rho | 0.184727 |
| sqrt(8)/rho | 0.139667 |

These choices have different meanings. 1.658 is only the coefficient matching
the GP kernel's e-folding distance; it is not a fitted optimum or a universal
recommended constant. `sqrt(8)/kappa` is a conventional practical range, not
an exact 0.1 threshold. GP uses B with exponential entries and penalizes with
the knot kernel Q, so its induced prior is `B Q^-1 B'`; matching one nominal
kernel threshold does not equate either entire kernel or the two induced priors.

`svgpc_spde_scale_audit(basis, reference_distance, model)` uses the actual FM
vertices and triangle-edge graph. It reports an exponential reference curve,
nominal Matern and finite-mesh SPDE prior correlation quantiles for first,
second and more distant vertex rings, using small batches of sparse solves.
The reference distance is an explicit physical distance, not a source of mesh
knots. Supplying a qr model also includes its
coefficient constraints. Use these summaries and sensitivity to buffer,
triangle size and kappa to assess a desired one-ring scale. The functions
do not assert that an arbitrary mesh or distance already satisfies it.

B's row has at most three weights: interpolation from a containing triangle.
A column's hat function spans triangles touching one vertex. Marginal prior
correlation and fitted propagation extend beyond those supports. The fitted
smoother also depends on lambda, W and covariates and need not have the same
decay as the unconditioned prior. No choice here imposes a hard one-ring cutoff.

## Mesh and input options

Construct geographic meshes from `svgpc_spde_locations(table)` to retain their
projection. `svgpc_spde_prepare_data(mesh, table, kappa, bed)` preserves the
selected nodes, triangles and coordinate frame after FAM matching. Out-of-mesh
locations cause an error. `mesh$xy` stores vertices; `mesh$tv` stores triangles.

`plot(mesh, kappa = scale$kappa)` draws one map: location-count density, small
mesh vertices and the nominal correlation-0.1 circle. GP uses the same display
with radius `rho * log(10)`. Both circles are centred at the density-grid peak.
Mesh construction retains only aggregated coordinates/counts for the background.
Without kappa, a mesh plot makes no range assumption. Older saved meshes can
supply a location object explicitly with `locations = ...` in the mesh's frame.

The default domain is a buffered rectangle. Custom planar `vertices`, one-based
boundary `segments` and `holes` support nonrectangular domains; coastlines and
barriers are not inferred. RTriangle can add nodes to satisfy area/angle limits
and is needed only for construction. Check the actual vertex count before fitting.

The README uses the bundled coarse Natural Earth US/UK convex envelopes,
projects them into the model's coordinate frame, and offsets the boundary by
50 km for the US or 5 km for the UK, independently of correlation range. `sf` performs that coordinate
operation; `svgpc_spde_mesh()` constructs the triangles. With custom vertices,
the `buffer` argument is metadata: it does not apply another offset. The full
mesh is visible in the plot; coastline details do not force local refinement.

For other genotype sources, `svgpc_spde_accumulate_reader` calls
`reader(start, end)` with consecutive one-based inclusive SNP indices. Each
block must already be standardized and in model location order. Matrix input
uses `svgpc_spde_accumulate_matrix` with the same contract.
