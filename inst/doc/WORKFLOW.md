# spatialGPCA: one-table BED workflow

The first four columns of the individual table are **FID, IID, Lat, Lon**.
All remaining columns are numeric fixed effects, exactly
`X <- as.matrix(table[, -c(1:4), drop = FALSE])`.
GPCs and any user-added income or other columns are treated identically.
An intercept is always included. No missingness, dosage, conditioning,
no-intercept or manually fixed-lambda option is part of this workflow.

```r
library(spatialGPCA)
cls <- svgpc_cluster(table, rho = 300, clusters = 24, bed = "merged")
plot(cls)
cls$rho <- 600
plot(cls)

model <- svgpc_prepare(cls)
svgpc_accumulate(model, "merged.afreq")
parameters <- svgpc_select_lambda(model, method = "REML") # or GCV, once
parameters
result <- svgpc_fit(model, parameters = parameters, components = 10)
plot(result)
PC <- result$outputs$raw_F_PCA$pcs

result$rho <- 300
plot(result)
```

To change the actual clustering, run `svgpc_cluster(table, rho = 300,
clusters = 40, bed = "merged")` again. This function invokes the independent
C++ clustering routine and returns the exact membership and centres used by
preparation. It never reads genotypes or estimates a smoothing penalty.
`clusters` is the initial count. Empty centres are removed during iteration;
the returned count may be smaller, and every individual remains assigned.
The initial count is retained in `cls$clustering$initial_clusters`.
After removing empty centres, set the threshold to floor(mean coordinate-point
count per cluster / 2). Clusters with coordinate-point counts at or below this
threshold are merged into the nearest current cluster. The threshold is computed
once, before small-cluster merging, and is not increased afterwards.
The smallest eligible class is processed first; ties use the existing class order.
The receiver can also be a small class. Each distinct Lat/Lon pair counts once;
the number of individuals sharing coordinates does not affect this rule.
No individual or location is removed.
Centres are updated using the existing distinct-location mean convention.
`cls$clustering` records the nonempty count and the threshold in `merge_at_most`;
the final count is `nrow(cls$centres)`.
The numerical fit consumes these centres without running clustering again.
For geography-only exploration, omit `bed`; provide it when creating the
geographic result that will be used for BED fitting.

`svgpc_fit(model, method = "GCV")` is the shorter equivalent when a separate
parameter-inspection step is unnecessary; it estimates only GCV.
There is no need to supply a numeric lambda. `svgpc_select_lambda()` requires
accumulated genotype statistics and is independent of geographic clustering.

## Inputs and sample alignment

The package matches FID + IID jointly, takes the intersection, and orders it
by FAM. Same-coordinate individuals share an internally generated location.
Counts and means of *all* fixed-effect columns use only matched individuals.
The BED reader uses the complete FAM sample count to decode bytes, skipping
excluded individuals during accumulation. Users do not prepare a location
column, location counts, or separate covariate matrices.

BED/BIM/FAM contain externally prepared high-quality, diploid biallelic
autosomal SNPs. PLINK2 remains the external QC/conversion/frequency tool:

```sh
plink2 --bfile merged --freq --out merged
```

The package reads the resulting ordinary `.afreq` directly and aligns ID,
REF, ALT and ALT_FREQS to the BIM counted allele. Swapped allele direction
uses 1-f. It does not perform QC, drop SNPs, infer strand complements, or
re-estimate frequency. Valid frequencies must be strictly between 0 and 1.
No PLINK executable or Python process is invoked during fitting.

Every missing BED call is filled with the counted-allele mean 2f before
location averaging. Each location SNP mean is centred at 2f and divided by
sqrt(2f(1-f)); its weight remains the complete matched location count.

## Geography and rho

Lat/Lon are input in degrees. The geographic object automatically constructs
a regional spherical azimuthal-equidistant coordinate system centred on the
mean direction of the distinct locations; numerical coordinates are in km.
Its centre and units are stored in `cls$projection`. Pairwise kernel distances
are Euclidean distances **in this projection**, not exact great-circle
separations for all pairs. Do not interpret unprojected degree differences as km.
The plot uses actual geographic coordinates over a map; contiguous-US inputs
show state boundaries, while other regions use a world outline.
The background shows relative individual density, computed by binning location
counts on a fixed display grid and smoothing those counts. It is a visual
summary, not a density per square kilometre. Red dots mark every actual
cluster centre; their size decreases for larger cluster counts. Individual dots,
per-cluster colours and numbered labels are omitted to keep large maps readable.
Display smoothing does not depend on rho or the number of clusters. Exact cluster
membership remains available in `cls$locations$cluster`.
The map shows a circle of radius rho around the cluster centre with the highest
displayed individual density (bilinearly interpolated at the centres). The circle
is constructed in the model's projected km coordinates and mapped back to Lat/Lon.
Its radius and the correlation curve update on rho edits. The map expands to
include the circle; the density calculation and selected centre remain unchanged.
Plots use solid red dots, compact legends, titles and coordinate ticks, without
subtitles, axis titles, captions or explanatory annotations.

rho is supplied by the user and controls the base kernel exp(-distance/rho).
Changing `cls$rho` before preparation changes the next model's kernel.
Changing `result$rho` afterwards changes only the displayed kernel curve;
it does not silently re-estimate the fitted PCs. The original fitted geography
remains in `result$spatial`. Reprepare and accumulate to fit a different rho.
The displayed curve is the base kernel, not an empirical genotype correlation
or the generally nonstationary covariance after the fixed-effect constraint.

## Model

For each SNP the response is its standardized location mean. With W=diag(Ni),
X is the location-mean fixed-effect design plus an intercept, Xq spans sqrt(W)X,
and Z=sqrt(W)B. The response is residualized by fixed-effect regression.
A separate full QR of Z'Xq supplies N spanning ker(Xq'Z), reducing the spatial
coefficient dimension from m to m-rank(Z'Xq). The spatial model uses BN with
precision N'QN. It does not replace this reduction with columnwise residuals
of the original B. The GP coefficients and observational noise remain independent.

REML/GCV select lambda=sigma_e^2/sigma_gp^2 from all accumulated SNPs.
Variances are not reported separately. The independent-noise assumption refers
to the statistical model, not exact orthogonality of realized noise vectors.
See `MATHEMATICS.md` for the weighted equations and zero-boundary limit.

## Installation and US example

```r
install.packages("remotes")
remotes::install_github("harryyiheyang/spatialGPCA")
```

A C++17 compiler and R's BLAS/LAPACK are required (Rtools on Windows).
The example at `inst/examples/us_full_fit.Rmd` simulates 150,000 SNPs and
3,000 US locations, writes a BED fileset and a standard PLINK2-format frequency
table, and fits spatial PCs and SNP-specific geographic effects. These are
synthetic data. The simulation itself does not invoke PLINK2.
