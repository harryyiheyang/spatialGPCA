\name{svgpc_spde}
\alias{svgpc_spde_correlation}
\alias{svgpc_spde_scale}
\alias{svgpc_spde_mesh}
\alias{svgpc_spde_basis}
\alias{svgpc_spde_prepare}
\alias{svgpc_spde_locations}
\alias{svgpc_spde_prepare_locations}
\alias{svgpc_spde_prepare_data}
\alias{svgpc_spde_select_lambda}
\alias{plot.svgpc_spde_mesh}
\alias{print.svgpc_spde_mesh}
\alias{print.svgpc_spde_parameters}
\alias{print.svgpc_spde_model}
\alias{print.svgpc_spde_locations}
\alias{svgpc_spde_accumulate}
\alias{svgpc_spde_accumulate_matrix}
\alias{svgpc_spde_accumulate_reader}
\alias{svgpc_spde_objective}
\alias{svgpc_spde_fit}
\alias{svgpc_spde_fitted}
\alias{svgpc_spde_scale_audit}
\title{Independent sparse finite-element spatial genetic PCA}
\description{Linear triangular finite elements with lumped mass and the proper
two-dimensional alpha-two Matern SPDE precision. Preserves sparse observation
projection and accumulates shared full-SNP statistics before dense mesh-space
transformations. Does not modify the existing GP path.}
\usage{
svgpc_spde_correlation(distance, kappa)
svgpc_spde_scale(distance, correlation = exp(-1))
svgpc_spde_mesh(locations, max_area, buffer, seed_vertices = NULL,
               vertices = NULL, segments = NULL, holes = NULL, min_angle = 21)
svgpc_spde_basis(mesh, locations, kappa, tau = 1)
svgpc_spde_prepare(basis, pcs = NULL, counts, conditioning = c("qr", "profile"),
                   accumulation = c("auto", "moments", "directions"))
svgpc_spde_locations(table, bed = NULL, projection = NULL)
svgpc_spde_prepare_locations(spatial, mesh, kappa, tau = 1,
                           conditioning = c("qr", "profile"),
                           accumulation = c("auto", "moments", "directions"))
svgpc_spde_prepare_data(mesh, table, kappa, bed = NULL, tau = 1,
                        conditioning = c("qr", "profile"),
                        accumulation = c("auto", "moments", "directions"))
svgpc_spde_accumulate(model, frequencies, block_size = 256L)
svgpc_spde_accumulate_matrix(model, G, block_size = 256L)
svgpc_spde_accumulate_reader(model, reader, snps, block_size = 256L)
svgpc_spde_objective(model, lambda, method = c("REML", "GCV"))
svgpc_spde_select_lambda(model, method = c("REML", "GCV"))
svgpc_spde_fit(model, method = c("REML", "GCV"), lambda = NULL,
               components = 10L, operators = "raw_F_PCA", parameters = NULL)
svgpc_spde_fitted(model, G, lambda, coefficients = FALSE)
svgpc_spde_scale_audit(basis, reference_distance, model = NULL,
                     probes = 16L, block_size = 64L)
\method{plot}{svgpc_spde_mesh}(x, show_vertices = TRUE, vertex_cex = .3,
                               main = NULL, xlab = NULL, ylab = NULL,
                               kappa = NULL, locations = x$locations, ...)
}
\arguments{
\item{x}{An SPDE mesh to plot.}
\item{show_vertices}{Draw dots at mesh vertices.}
\item{vertex_cex}{Mesh-vertex point size in mm; the default matches GP centre size.}
\item{main, xlab, ylab}{Plot title and axis labels.}
\item{...}{Reserved for S3 compatibility.}
\item{distance}{Nonnegative distances in km; a positive scalar for scale calibration.}
\item{kappa}{Positive inverse distance in km to the power minus one. In plot,
supply the scale intended for fitting to draw the nominal correlation-0.1 circle;
NULL plots the mesh and density without assuming a range.}
\item{correlation}{Target nominal Matern correlation strictly between zero and one.}
\item{locations}{For mesh construction, an object from \code{svgpc_spde_locations}
or a location-by-two planar coordinate matrix in km. Geographic objects retain
the projection needed by \code{svgpc_spde_prepare_data}. Basis construction takes
the planar coordinate matrix. In plot, a location object or its locations data
frame supplies coordinates and counts; by default these are retained by the mesh.
For older saved meshes, supply this explicitly in the mesh's projection.}
\item{max_area}{Maximum triangle area in square km.}
\item{buffer}{Positive extension in km for the default rectangular boundary. Required
even with a custom boundary, where it is metadata only and is not applied.}
\item{seed_vertices}{Optional extra mesh seed coordinates in km, used only with
the default rectangular boundary.}
\item{vertices}{Optional custom planar PSLG vertices.}
\item{segments}{Two-column, one-based boundary segment indices for custom vertices.}
\item{holes}{Optional RTriangle hole marker coordinates.}
\item{min_angle}{Triangle minimum angle in degrees, in (0,33].}
\item{mesh}{Mesh from \code{svgpc_spde_mesh}.}
\item{tau}{Positive scale entering the precision matrix as tau squared; fixed
with kappa before lambda selection.}
\item{basis}{Sparse basis from \code{svgpc_spde_basis}.}
\item{pcs}{Location-level fixed-effect matrix; an intercept is added.}
\item{counts}{Positive location population sizes, not numbers of SNPs.}
\item{conditioning}{Explicit coefficient constraint (qr) or fixed-effect profiling (profile).
These are different statistical models.}
\item{accumulation}{Numerical accumulation method, not a statistical model choice.
Moments transforms accumulated sparse-projection covariance once. Directions
transforms each sparse block before covariance to reduce cancellation.
Auto uses directions when the retained Gram condition exceeds 1e8.}
\item{table}{Individual metadata: first four columns FID, IID, Lat, Lon, followed by fixed covariates.}
\item{bed}{Optional prepared BED prefix for exact sample matching.}
\item{projection}{Optional original geographic projection to preserve its coordinate frame.}
\item{spatial}{Result of \code{svgpc_spde_locations}; an existing \code{svgpc_cluster}
object can also supply only observation metadata, without reusing its knots.}
\item{model}{SPDE model; optionally the matching model to audit its constrained prior.}
\item{frequencies}{File path or table with PLINK2 ID/REF/ALT/ALT_FREQS,
PLINK SNP/A1/A2/MAF, or ID/A1/A2/AF columns. Frequencies refer to ALT or A1.}
\item{block_size}{Positive SNP block size, or RHS block size for a scale audit.}
\item{G}{Already standardized location-by-SNP matrix in the model's location
order, finite with no missing entries.}
\item{reader}{Function(start,end), inclusive one-based SNP indices, returning an n by block
standardized matrix in exactly the model's location order. No automatic centering is applied.}
\item{snps}{Total number of SNPs read once, in consecutive blocks.}
\item{lambda}{Nonnegative common smoothing parameter, including infinity; NULL selects it.}
\item{method}{Full-SNP shared-scale REML or residual-space GCV.}
\item{components}{Number of leading nonzero PCA components.}
\item{operators}{Distinct choices: raw_F_PCA is ordinary PCA of raw fitted
location effects. H2, H and 2H-H2 select SNP loadings through the corresponding
smoother operators and return raw and weighted scores; see doc/SPDE.md.}
\item{parameters}{Result of svgpc_spde_select_lambda from this exact model. Reusing it
extracts the requested PCs without rereading BED or reoptimizing lambda.}
\item{coefficients}{Return mesh by SNP coefficients instead of raw location fits.}
\item{reference_distance}{Explicit physical distance in km for the reference exponential curve.}
\item{probes}{Number of deterministically selected mesh-vertex probes with observation support.}
}
\details{The default mesh is an expanded rectangle, not a coastline or barrier mesh.
RTriangle is used only for mesh construction. The audit uses actual stored triangle edges. Matrix and the
package native code handle fitting; sf, fmesher and INLA are not required.
The raw B has at most three nonzeros per row. Sparse Q is factorized without forming
its inverse. Mesh-space Gram matrices, accumulated covariance, whitening eigenvectors
and exact PCA are dense. See installed doc/SPDE.md for formulas and limitations.

The public stages are mesh construction plus plot (no BED required), then
prepare_data(mesh,table,kappa,bed) plus accumulation and select_lambda, then
fit(model,parameters,components). prepare_data uses the selected mesh and its
stored coordinate frame without regenerating nodes or triangles. Construct the
geographic mesh from svgpc_spde_locations(table) to retain this frame. Node count
and triangle count are distinct and are shown by print and plot. Changing PC count
never changes the full spectrum or all-SNP information used for lambda selection.

Mesh plots draw one map with the same density grid, colours and small-node styling
as GP. The heatmap uses log1p-transformed smoothed location counts, not mesh-node density; planar matrix
inputs have unit counts. For supplied kappa, the circle solves
(kappa*d) K1(kappa*d) = 0.1 and is centred at the density-grid peak. It is a nominal
Matern reference, not a fitted influence region or hard cutoff. Geographic maps
show the full buffered mesh as well as the study area.
Plot returns one ggplot invisibly, with a range_circle attribute containing
radius, centre and xy in projected km. Its path contains Lon/Lat in degrees for
geographic input, or x/y in km for planar input.

Accumulation returns a new model list: always assign the returned value. All SNPs
share one lambda and variance scale. A prepared list can be reused for independent
datasets; an accumulated list cannot be accumulated again. Use saveRDS/readRDS to
save these objects. Existing GP cache, plot and fit methods do not dispatch to them.

The scale helper solves z K1(z) = correlation and returns kappa = z/distance.
Matching exp(-1) at the GP distance rho yields approximately 1.6580762/rho;
this is a nominal threshold match, not a data-optimal scale or equality of covariance
operators. GP clustering and SPDE FM vertices are independent knot-generation paths;
the SPDE location preparation function does not run GP clustering. The audit uses
the actual SPDE triangle-edge vertex graph, defining first and second rings by
one and two edge hops, with physical distance summaries and sparse covariance
solves in blocks. Vertex adjacency differs from triangle-element adjacency.
Unconstrained or qr-constrained mesh-vertex prior correlations are not fitted-smoother influence.}
\value{Scale and mesh functions return lists; basis returns sparse B and Q plus FEM
matrices. Preparation and accumulation return an svgpc_spde_model list. Objective
returns a scalar. Fit returns parameter selection and requested PCA outputs. Fitted
returns a dense matrix for the supplied block only. Audit returns aggregate summaries.}
\examples{
svgpc_spde_scale(100)
svgpc_spde_correlation(c(0, 100, 200), 2.5 / 100)
}
