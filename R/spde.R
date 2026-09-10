# Independent sparse finite-element path; the GP interfaces are unchanged.
.spde_positive <- function(x, label) {
    if (length(x) != 1L || !is.finite(x) || x <= 0) stop(label, " must be positive and finite")
    as.numeric(x)
}
.spde_mm <- function(A, B, transA = FALSE, transB = FALSE) {
    if (requireNamespace("CppMatrix", quietly = TRUE))
        return(CppMatrix::matrixMultiply(A, B, transA = transA, transB = transB))
    if (transA) A <- t(A)
    if (transB) B <- t(B)
    A %*% B
}
.spde_eigen <- function(A) {
    A <- (A + t(A)) / 2
    if (requireNamespace("CppMatrix", quietly = TRUE)) {
        E <- CppMatrix::matrixEigen(A)
        ord <- order(E$values, decreasing = TRUE)
        return(list(values = E$values[ord], vectors = E$vectors[, ord, drop = FALSE]))
    }
    eigen(A, symmetric = TRUE)
}
.spde_qr <- function(A) {
    if (!ncol(A)) return(matrix(0, nrow(A), 0L))
    q <- qr(A, LAPACK = TRUE)
    d <- abs(diag(qr.R(q)))
    rank <- sum(d > .Machine$double.eps * max(dim(A)) * d[1L])
    qr.Q(q)[, seq_len(rank), drop = FALSE]
}

svgpc_spde_correlation <- function(distance, kappa) {
    kappa <- .spde_positive(kappa, "kappa")
    if (any(!is.finite(distance)) || any(distance < 0)) stop("distance must be finite and nonnegative")
    z <- distance * kappa
    ans <- rep(1, length(z))
    ans[is.infinite(z)] <- 0
    at <- z > 0 & is.finite(z)
    ans[at] <- exp(log(z[at]) - z[at]) * besselK(z[at], 1, expon.scaled = TRUE)
    ans
}

svgpc_spde_scale <- function(distance, correlation = exp(-1)) {
    distance <- .spde_positive(distance, "distance")
    if (length(correlation) != 1L || !is.finite(correlation) || correlation <= 0 || correlation >= 1)
        stop("correlation must lie strictly between zero and one")
    z <- stats::uniroot(function(z) svgpc_spde_correlation(z, 1) - correlation,
                        c(0, 1000), tol = 1e-12)$root
    list(distance_km = distance, correlation = correlation,
         kappa = z / distance, kappa_units = "km^-1", practical_range_km = sqrt(8) * distance / z)
}

svgpc_spde_mesh <- function(locations, max_area, buffer, seed_vertices = NULL,
                            vertices = NULL, segments = NULL, holes = NULL, min_angle = 21) {
    if (!requireNamespace("RTriangle", quietly = TRUE)) stop("Install RTriangle to construct meshes")
    coordinate_system <- list(name = "Supplied planar coordinates", units = "km", center = c(0, 0), scale = 1)
    obs <- NULL
    if (inherits(locations, "svgpc_spde_locations")) {
        coordinate_system <- locations$projection
        obs <- locations$locations
        locations <- locations$locations[, c("x", "y")]
    }
    custom_boundary <- !is.null(vertices)
    xy <- .matrix_double(locations)
    if (ncol(xy) != 2L || nrow(xy) < 3L || any(!is.finite(xy))) stop("locations must be finite planar km coordinates")
    if (is.null(obs)) obs <- data.frame(x = xy[, 1], y = xy[, 2], count = 1)
    max_area <- .spde_positive(max_area, "max_area (km^2)")
    buffer <- .spde_positive(buffer, "buffer (km)")
    if (length(min_angle) != 1L || !is.finite(min_angle) || min_angle <= 0 || min_angle > 33)
        stop("min_angle must lie in (0, 33]")
    if (is.null(vertices)) {
        # Expanded bounding rectangle: explicit finite-domain approximation.
        x <- range(xy[, 1]) + c(-buffer, buffer)
        y <- range(xy[, 2]) + c(-buffer, buffer)
        vertices <- rbind(c(x[1], y[1]), c(x[2], y[1]), c(x[2], y[2]), c(x[1], y[2]))
        segments <- cbind(1:4, c(2:4, 1))
        if (!is.null(seed_vertices)) vertices <- rbind(vertices, .matrix_double(seed_vertices))
    }
    vertices <- .matrix_double(vertices)
    if (ncol(vertices) != 2L || any(!is.finite(vertices)) || anyDuplicated(as.data.frame(vertices)))
        stop("vertices must be finite distinct planar coordinates")
    if (is.null(segments)) stop("Supply boundary segments with custom vertices")
    segments <- as.matrix(segments)
    if (ncol(segments) != 2L || anyNA(segments) || any(segments < 1 | segments > nrow(vertices)) ||
        any(segments != floor(segments))) stop("Invalid boundary segments")
    p <- RTriangle::pslg(P = vertices, S = segments, H = if (is.null(holes)) NA else holes)
    mesh <- RTriangle::triangulate(p, a = max_area, q = min_angle)
    edges <- rbind(mesh$T[, 1:2], mesh$T[, 2:3], mesh$T[, c(3, 1)])
    edges <- cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
    keys <- paste(edges[, 1], edges[, 2], sep = ":")
    boundary_edges <- edges[!duplicated(keys) & !duplicated(keys, fromLast = TRUE), , drop = FALSE]
    out <- list(xy = mesh$P, tv = mesh$T, locations = obs, units = "km", max_area = max_area,
                buffer = buffer, domain = if (custom_boundary) "custom PSLG" else "buffered rectangle",
                vertex_count = nrow(mesh$P), triangle_count = nrow(mesh$T),
                boundary_edges = boundary_edges, coordinate_system = coordinate_system,
                parameters = list(max_area = max_area, buffer = buffer, min_angle = min_angle,
                                  vertices = vertices, segments = segments, holes = holes,
                                  seed_vertices = seed_vertices, custom_boundary = custom_boundary))
    class(out) <- "svgpc_spde_mesh"
    out
}

svgpc_spde_basis <- function(mesh, locations, kappa, tau = 1) {
    if (!inherits(mesh, "svgpc_spde_mesh")) stop("mesh must come from svgpc_spde_mesh()")
    kappa <- .spde_positive(kappa, "kappa (km^-1)")
    tau <- .spde_positive(tau, "tau")
    xy <- .matrix_double(mesh$xy)
    tv <- as.matrix(mesh$tv)
    loc <- .matrix_double(locations)
    hit <- cpp_spde_project(xy, tv, loc)
    B <- Matrix::drop0(Matrix::sparseMatrix(i = hit$i, j = hit$j, x = hit$x,
                         dims = c(nrow(loc), nrow(xy))))
    a <- xy[tv[, 1], , drop = FALSE]
    b <- xy[tv[, 2], , drop = FALSE]
    c <- xy[tv[, 3], , drop = FALSE]
    area <- abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) -
                (c[, 1] - a[, 1]) * (b[, 2] - a[, 2])) / 2
    if (any(!is.finite(area)) || any(area <= 0)) stop("Degenerate FEM triangles")
    m <- nrow(xy)
    D <- as.numeric(Matrix::sparseMatrix(i = as.vector(tv), j = rep(1L, length(tv)),
                     x = rep(area / 3, 3L), dims = c(m, 1L)))
    if (any(D <= 0)) stop("Mesh contains unused vertices")
    gx <- cbind(b[, 2] - c[, 2], c[, 2] - a[, 2], a[, 2] - b[, 2])
    gy <- cbind(c[, 1] - b[, 1], a[, 1] - c[, 1], b[, 1] - a[, 1])
    nt <- nrow(tv)
    ii <- jj <- integer(nt * 9L)
    xx <- numeric(nt * 9L)
    at <- 0L
    for (j in 1:3) for (k in 1:3) {
        idx <- at + seq_len(nt)
        ii[idx] <- tv[, j]; jj[idx] <- tv[, k]
        xx[idx] <- (gx[, j] * gx[, k] + gy[, j] * gy[, k]) / (4 * area)
        at <- at + nt
    }
    K <- Matrix::forceSymmetric(Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(m, m)))
    Q <- Matrix::forceSymmetric(tau^2 * (kappa^4 * Matrix::Diagonal(m, D) +
         2 * kappa^2 * K + K %*% Matrix::Diagonal(m, 1 / D) %*% K))
    list(B = B, Q = Q, mass = D, stiffness = K, mesh = mesh,
         kappa = kappa, tau = tau, practical_range_km = sqrt(8) / kappa)
}

svgpc_spde_prepare <- function(basis, pcs = NULL, counts, conditioning = c("qr", "profile"),
                               accumulation = c("auto", "moments", "directions")) {
    start <- proc.time()[[3L]]
    conditioning <- match.arg(conditioning)
    accumulation <- match.arg(accumulation)
    B <- basis$B; Q <- basis$Q
    if (!inherits(B, "sparseMatrix") || !inherits(Q, "sparseMatrix")) stop("Use a sparse SPDE basis")
    n <- nrow(B); m <- ncol(B)
    counts <- as.numeric(counts)
    if (length(counts) != n || any(!is.finite(counts)) || any(counts <= 0)) stop("counts must be positive, one per location")
    if (is.null(pcs)) pcs <- matrix(0, n, 0L)
    pcs <- .matrix_double(pcs)
    if (nrow(pcs) != n || any(!is.finite(pcs))) stop("pcs must have one finite row per location")
    Xq <- .spde_qr(cbind(1, pcs) * sqrt(counts))
    r <- n - ncol(Xq)
    if (r < 1) stop("Fixed effects leave no residual degrees of freedom")
    H <- as.matrix(Matrix::crossprod(B, Xq * sqrt(counts)))
    S <- Matrix::crossprod(B, Matrix::Diagonal(n, counts) %*% B)
    projection_scale <- sqrt(as.numeric(Matrix::diag(S)))
    S <- as.matrix(S) - .spde_mm(H, H, transB = TRUE)
    # P Q P' = L L'; sparse factorization and triangular solves, no Q inverse.
    fac <- Matrix::expand(Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE, super = FALSE))
    L <- fac$L; P <- fac$P
    Sw <- as.matrix(Matrix::solve(L, P %*% S %*% Matrix::t(P)))
    Sw <- t(as.matrix(Matrix::solve(L, t(Sw))))
    constraint_rank <- 0L
    if (conditioning == "qr") {
        A <- .spde_qr(as.matrix(Matrix::solve(L, P %*% H)))
        constraint_rank <- ncol(A)
        if (constraint_rank >= m) stop("Covariate constraints remove all mesh directions")
        # Thin projection only in m-space: no full null basis and no dense B N.
        Sw <- Sw - .spde_mm(A, .spde_mm(A, Sw, transA = TRUE))
        Sw <- Sw - .spde_mm(.spde_mm(Sw, A), A, transB = TRUE)
    }
    E <- .spde_eigen(Sw)
    tol <- .Machine$double.eps * max(n, m) * max(1, max(E$values))
    if (min(E$values) < -tol) stop("Material negative whitened Gram eigenvalue")
    keep <- which(E$values > tol)
    if (!length(keep)) stop("No estimable spatial directions")
    d <- E$values[keep]
    V <- as.matrix(Matrix::t(P) %*% Matrix::solve(Matrix::t(L), E$vectors[, keep, drop = FALSE]))
    VD <- sweep(V, 2L, sqrt(d), "/")
    metric <- .spde_mm(VD, as.matrix(Matrix::crossprod(B) %*% VD), transA = TRUE)
    R <- chol((metric + t(metric)) / 2)
    active <- which(Matrix::colSums(abs(B)) > 0)
    gram_condition <- max(d) / min(d)
    if (accumulation == "auto") accumulation <- if (gram_condition > 1e8) "directions" else "moments"
    out <- list(B = B, Q = Q, mesh = basis$mesh, Xq = Xq, H = H, counts = counts, V = V, VD = VD, d = d,
                active = active,
                projection_scale = projection_scale[active],
                accumulation_basis = B[, active, drop = FALSE] %*% Matrix::Diagonal(length(active), 1 / projection_scale[active]),
                accumulation = accumulation, gram_condition = gram_condition,
                accumulation_projection = if (accumulation == "directions") VD[active, , drop = FALSE] * projection_scale[active] else matrix(0, 0L, 0L),
                R = R, n = n, m = m, r = r, conditioning = conditioning,
                constraint_rank = constraint_rank, kappa = basis$kappa, tau = basis$tau,
                practical_range_km = basis$practical_range_km,
                prepare_seconds = proc.time()[[3L]] - start, ready = FALSE)
    class(out) <- "svgpc_spde_model"
    out
}

svgpc_spde_prepare_locations <- function(spatial, mesh, kappa, tau = 1,
                                      conditioning = c("qr", "profile"),
                                      accumulation = c("auto", "moments", "directions")) {
    if (!inherits(spatial, "svgpc_spde_locations") && !inherits(spatial, "svgpc_cluster"))
        stop("Use svgpc_spde_locations() or reuse only the location metadata of svgpc_cluster()")
    X <- as.matrix(spatial$table[, -c(1:4), drop = FALSE])
    n <- nrow(spatial$locations)
    C <- if (ncol(X)) rowsum(X, spatial$sample_location, reorder = TRUE) /
        spatial$locations$count else matrix(0, n, 0L)
    basis <- svgpc_spde_basis(mesh, spatial$locations[, c("x", "y")], kappa, tau)
    out <- svgpc_spde_prepare(basis, C, spatial$locations$count, match.arg(conditioning), match.arg(accumulation))
    out$spatial <- spatial
    out
}

.spde_fresh <- function(model) {
    if (!inherits(model, "svgpc_spde_model")) stop("Expected an SPDE model")
    if (model$ready) stop("Statistics already accumulated; prepare a fresh model to avoid double counting")
}
.spde_finish <- function(model, stats) {
    start <- proc.time()[[3L]]
    if (stats$covariance_space == "directions") {
        model$cov <- stats$cross_covariance
    } else {
        VD <- model$VD[model$active, , drop = FALSE] * model$projection_scale
        model$cov <- .spde_mm(VD, .spde_mm(stats$cross_covariance, VD), transA = TRUE)
    }
    model$cov <- (model$cov + t(model$cov)) / 2
    model$a <- diag(model$cov)
    model$energy <- stats$energy
    model$b <- model$energy - sum(model$a)
    if (model$b < -1e-8 * max(model$energy, 1)) stop("Negative complement energy")
    model$b <- if (length(model$d) == model$r) 0 else max(0, model$b)
    model$p <- stats$p
    model$missing <- stats$missing
    model$timing <- c(stats$timing, transform_covariance = proc.time()[[3L]] - start)
    model$ready <- TRUE
    model
}

svgpc_spde_accumulate_matrix <- function(model, G, block_size = 256L) {
    .spde_fresh(model)
    stats <- cpp_spde_accumulate_matrix(model$accumulation_basis, model$Xq, model$counts, .matrix_double(G),
                                       .positive_integer(block_size, "block_size"), model$accumulation_projection)
    .spde_finish(model, stats)
}
svgpc_spde_accumulate_reader <- function(model, reader, snps, block_size = 256L) {
    .spde_fresh(model)
    if (!is.function(reader)) stop("reader must be a function(start, end) returning standardized location x SNP blocks")
    stats <- cpp_spde_accumulate_reader(model$accumulation_basis, model$Xq, model$counts, reader,
                 .positive_integer(snps, "snps"), .positive_integer(block_size, "block_size"), model$accumulation_projection)
    .spde_finish(model, stats)
}

svgpc_spde_objective <- function(model, lambda, method = c("REML", "GCV")) {
    method <- match.arg(method)
    if (!inherits(model, "svgpc_spde_model") || !model$ready) stop("Accumulate all SNPs before fitting")
    if (length(lambda) != 1L || is.na(lambda) || lambda < 0) stop("lambda must be nonnegative")
    if (model$energy <= 0) stop("No residual variation")
    if (is.infinite(lambda)) return(0)
    d <- model$d; a <- model$a; r <- model$r
    if (length(d) == r) {
        scale <- max(d, lambda)
        v <- d / scale + lambda / scale
        if (method == "REML") return(sum(log(v)) + r * log(sum(a / v) / model$energy))
        return(log(sum(a / v^2) / model$energy) - 2 * log(sum(1 / v) / r))
    }
    z <- lambda / (d + lambda)
    if (method == "REML") {
        if (lambda == 0) return(Inf)
        return(sum(log1p(d / lambda)) + r * log((model$b + sum(a * z)) / model$energy))
    }
    log((model$b + sum(a * z^2)) / model$energy) - 2 * log(1 - sum(1 - z) / r)
}

svgpc_spde_select_lambda <- function(model, method = c("REML", "GCV")) {
    start <- proc.time()[[3L]]
    method <- match.arg(method)
    if (!inherits(model, "svgpc_spde_model") || !model$ready) stop("Accumulate all SNPs before fitting")
    eta <- log(stats::median(model$d)) + seq(-30, 30, by = 0.5)
    vals <- vapply(exp(eta), function(x) svgpc_spde_objective(model, x, method), numeric(1))
    lambda <- Inf; best <- 0; boundary <- "infinity"
    if (method == "GCV" || length(model$d) == model$r) {
        v <- svgpc_spde_objective(model, 0, method)
        if (v < best) { best <- v; lambda <- 0; boundary <- "zero" }
    }
    for (j in 2:(length(eta) - 1L)) if (vals[j] <= vals[j - 1L] && vals[j] <= vals[j + 1L]) {
        opt <- stats::optimize(function(e) svgpc_spde_objective(model, exp(e), method),
                               eta[c(j - 1L, j + 1L)], tol = 1e-8)
        if (opt$objective < best - 1e-12) {
            best <- opt$objective; lambda <- exp(opt$minimum); boundary <- "finite"
        }
    }
    if (min(vals[c(1L, length(vals))]) < best - 1e-8) stop("Optimum reaches search boundary; extend bounds")
    out <- list(lambda = lambda, selector = method, boundary = boundary, objective_relative = best,
                p = model$p, rank = length(model$d), elapsed_seconds = proc.time()[[3L]] - start)
    class(out) <- "svgpc_spde_parameters"
    attr(out, "model") <- model
    out
}

svgpc_spde_fit <- function(model, method = c("REML", "GCV"), lambda = NULL,
                           components = 10L, operators = "raw_F_PCA", parameters = NULL) {
    start <- proc.time()[[3L]]
    supplied_method <- !missing(method)
    method <- match.arg(method)
    if (!inherits(model, "svgpc_spde_model") || !model$ready) stop("Accumulate all SNPs before fitting")
    components <- .positive_integer(components, "components")
    if (!length(operators) || anyNA(operators) || anyDuplicated(operators) ||
        any(!operators %in% c("raw_F_PCA", "H2", "H", "2H-H2"))) stop("Unknown or duplicate PCA operators")
    if (is.null(parameters) && is.null(lambda)) parameters <- svgpc_spde_select_lambda(model, method)
    if (!is.null(parameters)) {
        if (!is.null(lambda)) stop("Supply parameters or lambda, not both")
        if (!inherits(parameters, "svgpc_spde_parameters") || !identical(attr(parameters, "model"), model))
            stop("parameters must be selected from this model by svgpc_spde_select_lambda()")
        if (supplied_method && method != parameters$selector) stop("method disagrees with parameters")
        lambda <- parameters$lambda; method <- parameters$selector; selector <- method
        best <- parameters$objective_relative; boundary <- parameters$boundary
    } else {
        best <- svgpc_spde_objective(model, lambda, method)
        boundary <- "supplied"; selector <- "fixed_lambda"
    }
    h <- model$d / (model$d + lambda)
    outputs <- list()
    for (op in operators) {
        if (!any(h > 0)) {
            outputs[[op]] <- if (op == "raw_F_PCA") list(pcs = matrix(0, model$n, 0L), eigenvalues = numeric()) else
                list(raw_scores = matrix(0, model$n, 0L), weighted_scores = matrix(0, model$n, 0L))
            next
        }
        if (op == "raw_F_PCA") {
            RH <- sweep(model$R, 2L, h, "*")
            E <- .spde_eigen(.spde_mm(.spde_mm(RH, model$cov), RH, transB = TRUE) / model$p)
        } else {
            w <- switch(op, H2 = h^2, H = h, `2H-H2` = 2 * h - h^2)
            E <- .spde_eigen(model$cov * outer(sqrt(w), sqrt(w)))
        }
        idx <- utils::head(which(E$values > .Machine$double.eps * max(dim(model$cov)) * max(1, E$values[1L])), components)
        U <- E$vectors[, idx, drop = FALSE]
        if (op == "raw_F_PCA") {
            map <- .spde_mm(model$VD, backsolve(model$R, U))
            outputs[[op]] <- list(pcs = as.matrix(model$B %*% map), eigenvalues = E$values[idx])
        } else {
            coord <- .spde_mm(model$cov, U * sqrt(w)) * h
            coord <- sweep(coord, 2L, sqrt(E$values[idx]), "/")
            raw <- as.matrix(model$B %*% .spde_mm(model$VD, coord))
            weighted <- raw * sqrt(model$counts)
            weighted <- weighted - .spde_mm(model$Xq, .spde_mm(model$Xq, weighted, transA = TRUE))
            outputs[[op]] <- list(raw_scores = raw, weighted_scores = weighted, eigenvalues = E$values[idx])
        }
    }
    list(lambda = lambda, selector = selector, criterion = method, boundary = boundary, objective_relative = best,
         spatial_df = sum(h), conditioning = model$conditioning, kappa = model$kappa,
         outputs = outputs, p = model$p, elapsed_seconds = proc.time()[[3L]] - start)
}

svgpc_spde_fitted <- function(model, G, lambda, coefficients = FALSE) {
    if (length(lambda) != 1L || is.na(lambda) || lambda < 0) stop("lambda must be nonnegative")
    G <- .matrix_double(G)
    if (nrow(G) != model$n || any(!is.finite(G))) stop("Invalid location genotype block")
    U <- as.matrix(Matrix::crossprod(model$B, G * model$counts)) -
        .spde_mm(model$H, .spde_mm(model$Xq, G * sqrt(model$counts), transA = TRUE))
    J <- .spde_mm(model$VD, U, transA = TRUE)
    theta <- .spde_mm(model$VD, J * (model$d / (model$d + lambda)))
    if (coefficients) return(theta)
    as.matrix(model$B %*% theta)
}
