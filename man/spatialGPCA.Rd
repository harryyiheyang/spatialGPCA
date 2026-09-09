\name{spatialGPCA}
\alias{spatialGPCA}
\alias{svgpc_cluster}
\alias{svgpc_prepare}
\alias{svgpc_accumulate}
\alias{svgpc_accumulate_matrix}
\alias{svgpc_select_lambda}
\alias{svgpc_fit}
\alias{svgpc_fitted}
\alias{svgpc_info}
\alias{svgpc_statistics}
\alias{svgpc_save}
\alias{svgpc_load}
\alias{plot.svgpc_cluster}
\alias{plot.svgpc_fit}
\alias{print.svgpc_cluster}
\alias{print.svgpc_fit}
\alias{print.svgpc_parameters}
\title{Spatial Genetic PCA from an Individual Table and BED}
\description{Explore independent geographic clustering, then use prepared BED genotypes to select one smoothing penalty and obtain the spatial PC subspace.}
\usage{
svgpc_cluster(table, rho, clusters = NULL, bed = NULL, seed = 44L)
svgpc_prepare(spatial)
svgpc_accumulate(model, frequencies, block_size = 256L)
svgpc_accumulate_matrix(model, G, block_size = 256L)
svgpc_select_lambda(model, method = c("REML", "GCV"))
svgpc_fit(model, method = c("REML", "GCV"), operators = "raw_F_PCA",
          components = 50L, parameters = NULL)
svgpc_fitted(model, G, parameters, coefficients = FALSE)
svgpc_info(model)
svgpc_statistics(model)
svgpc_save(model, file, compress = FALSE)
svgpc_load(file)
\method{plot}{svgpc_cluster}(x, ...)
\method{plot}{svgpc_fit}(x, ...)
\method{print}{svgpc_cluster}(x, ...)
\method{print}{svgpc_fit}(x, ...)
\method{print}{svgpc_parameters}(x, ...)
}
\arguments{
\item{table}{Individual data frame. First four columns are FID, IID, Lat and Lon, in that order. Every remaining column is a numeric fixed effect.}
\item{rho}{Positive base-kernel distance scale in km; supplied before genotype fitting.}
\item{clusters}{Initial number of geographic clusters; NULL uses min(50, number of distinct locations). Empty centres are removed during iteration, so the returned count may be smaller. All individuals remain assigned. After empty-centre removal, the threshold is floor(mean distinct-coordinate count per class / 2). Classes at or below that coordinate-count threshold are merged into the nearest current class. Multiple individuals sharing one Lat/Lon pair count as one coordinate for this rule. The threshold is held fixed during merging; centres use the existing distinct-location mean convention.}
\item{bed}{Prepared BED/BIM/FAM prefix. If supplied, use the FID/IID intersection in FAM order before clustering. Omit for geography-only exploration.}
\item{seed}{Seed for native geographic clustering.}
\item{spatial}{Result from svgpc_cluster, possibly after editing its rho.}
\item{model}{Prepared native context.}
\item{frequencies}{PLINK2 .afreq filename or frequency data frame. ID/REF/ALT/ALT_FREQS and allele-labelled ID/A1/A2/AF or SNP/A1/A2/MAF tables are accepted.}
\item{block_size}{Number of SNPs per native block.}
\item{G}{Already standardized location-by-SNP matrix for numerical reference or fitted-block evaluation.}
\item{method}{One data-based smoothing selection criterion; only the selected method is evaluated.}
\item{operators}{One or more of raw_F_PCA, H2, H, 2H-H2.}
\item{components}{Number of leading components, at most 50.}
\item{parameters}{Selection result from svgpc_select_lambda for this model. NULL in svgpc_fit performs selection once using method.}
\item{coefficients}{Return original knot-coordinate coefficients instead of fitted location values.}
\item{file}{Native-state cache filename.}
\item{compress}{Passed to saveRDS.}
\item{x}{Geographic, parameter-selection or fitted result.}
\item{...}{Reserved for S3 compatibility.}
}
\value{Clustering returns actual location membership, centres, matched individuals, counts, rho and projection metadata. Preparation returns a native context; accumulation invisibly returns it. Selection returns lambda and criterion diagnostics. Fit returns the requested spatial PC outputs and geographic metadata. Plot draws the map and base-kernel correlation curve and invisibly returns the two ggplot objects.}
\details{The fixed-effect matrix is exactly the columns after the first four; location means use matched individuals, and an intercept is always added. A full QR constructs the coefficient-space null basis N satisfying C'WBN=0, reducing the geographic basis dimension; it is not columnwise residualization of B. The coefficient prior precision is N'QN; the default PCA uses the raw-location metric. Missing BED genotypes are always filled with the counted-allele 2f before aggregation and HWE scaling. Location weights include all matched individuals. There is no PGEN, dosage, missingness or no-intercept switch.

Lat/Lon are transformed automatically to a regional spherical azimuthal-equidistant projection in km. Kernel distances are Euclidean within this projection. The map background shows relative individual density from smoothed location counts on a fixed display grid, not density per square kilometre. Red dots mark every cluster centre, with smaller symbols when there are more centres. A circle with projected radius rho is centred on the cluster centre with the highest displayed individual density. The map expands to include this circle. Display smoothing is independent of rho and clustering. plot reads the current rho each time; changing a fitted result's rho updates its curve without silently refitting the PCs. A new preparation and accumulation is required for a different fitted rho. The curve is the base kernel, not the fixed-effect-constrained GP covariance.}
\section{Computation timings}{After preparation, \code{svgpc_info(model)$seconds} contains wall-clock seconds for the completed native operations. Names beginning with \code{prepare_detail.} distinguish kernel construction, the fixed-effect QR, the constraint crossproduct and full QR, basis/prior transformations, Gram formation, Cholesky factorizations, triangular solves, symmetric eigendecomposition, design construction and the raw-location metric. The numbered detail entries follow execution order. Aggregate stage timings overlap these details; do not sum aggregate and detail entries together. Timings are available after preparation returns and are not live progress messages.}
\examples{
set.seed(9)
table <- data.frame(FID = 1:60, IID = 1:60,
                    Lat = runif(60, 30, 45), Lon = runif(60, -115, -80),
                    GPC1 = rnorm(60), income = rnorm(60))
cls <- svgpc_cluster(table, rho = 300, clusters = 15)
cls
model <- svgpc_prepare(cls)
G <- matrix(rnorm(60 * 100), 60, 100)
svgpc_accumulate_matrix(model, G)
parameters <- svgpc_select_lambda(model, "GCV")
fit <- svgpc_fit(model, parameters = parameters, components = 3)
fit
}
