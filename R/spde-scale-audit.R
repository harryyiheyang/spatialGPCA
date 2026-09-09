svgpc_spde_scale_audit <- function(basis, reference_distance, model = NULL, probes = 16L,
                                  block_size = 64L) {
    xy <- .matrix_double(basis$mesh$xy)
    tv <- basis$mesh$tv
    if (ncol(xy) != 2L || nrow(xy) < 3L || any(!is.finite(xy)) || is.null(tv))
        stop("basis must contain its actual SPDE mesh vertices and triangles")
    reference_distance <- .spde_positive(reference_distance, "reference_distance (km)")
    active <- which(Matrix::colSums(abs(basis$B)) > 0)
    probes <- min(length(active), .positive_integer(probes, "probes"))
    block_size <- .positive_integer(block_size, "block_size")
    edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(3, 1)])
    A <- Matrix::sparseMatrix(i = c(edges[, 1], edges[, 2]), j = c(edges[, 2], edges[, 1]),
                              x = 1, dims = c(nrow(xy), nrow(xy))) != 0
    B <- Matrix::Diagonal(nrow(xy))
    Qf <- Matrix::Cholesky(basis$Q, perm = TRUE, LDL = FALSE)
    correction <- NULL
    if (!is.null(model)) {
        if (!inherits(model, "svgpc_spde_model") || !identical(model$Q, basis$Q))
            stop("model must use this basis precision")
        if (model$conditioning == "qr") {
            H <- .spde_qr(model$H)
            U <- as.matrix(Matrix::solve(Qf, H, system = "A"))
            R <- chol(.spde_mm(H, U, transA = TRUE))
            correction <- t(backsolve(R, t(U), transpose = TRUE))
        }
    }
    # Mesh-vertex variances in small RHS blocks, never a full Q inverse.
    variance <- numeric(nrow(B))
    for (j in seq(1L, nrow(B), by = block_size)) {
        idx <- j:min(nrow(B), j + block_size - 1L)
        V <- as.matrix(Matrix::solve(Qf, Matrix::t(B[idx, , drop = FALSE]), system = "A"))
        variance[idx] <- V[cbind(idx, seq_along(idx))]
    }
    if (!is.null(correction)) variance <- variance - rowSums(correction^2)
    if (any(variance <= 0)) stop("Nonpositive marginal variance at mesh vertices")
    seeds <- active[unique(round(seq(1, length(active), length.out = probes)))]
    V <- as.matrix(Matrix::solve(Qf, Matrix::t(B[seeds, , drop = FALSE]), system = "A"))
    covariance <- V
    if (!is.null(correction)) covariance <- covariance -
        .spde_mm(correction, correction[seeds, , drop = FALSE], transB = TRUE)
    correlation <- covariance / sqrt(outer(variance, variance[seeds]))
    rows <- vector("list", length(seeds))
    for (j in seq_along(seeds)) {
        i <- seeds[j]
        d <- sqrt(rowSums(sweep(xy, 2L, xy[i, ], "-")^2))
        one <- which(A[, i])
        two <- setdiff(which(Matrix::rowSums(A[, one, drop = FALSE]) > 0), c(i, one))
        ring <- rep("beyond_two", nrow(B)); ring[two] <- "second_ring"; ring[one] <- "first_ring"
        rows[[j]] <- data.frame(ring = ring[-i], distance_km = d[-i],
             reference_exponential = exp(-d[-i] / reference_distance),
             nominal_spde = svgpc_spde_correlation(d[-i], basis$kappa),
             finite_mesh_spde = correlation[-i, j])
    }
    pairs <- do.call(rbind, rows)
    summary <- list()
    for (ring in c("first_ring", "second_ring", "beyond_two")) {
        z <- pairs[pairs$ring == ring, -1L, drop = FALSE]
        if (!nrow(z)) next
        for (metric in names(z)) summary[[paste(ring, metric)]] <- data.frame(
            ring = ring, metric = metric, pairs = nrow(z),
            q10 = unname(stats::quantile(z[[metric]], .1)), median = stats::median(z[[metric]]),
            q90 = unname(stats::quantile(z[[metric]], .9)))
    }
    ee <- Matrix::summary(Matrix::triu(A))
    edge_distance <- sqrt(rowSums((xy[ee$i, , drop = FALSE] - xy[ee$j, , drop = FALSE])^2))
    list(summary = do.call(rbind, summary), mesh_vertices = nrow(B), probes = length(seeds),
         reference_distance_km = reference_distance,
         edge_distance_km = stats::quantile(edge_distance, c(.1, .5, .9)), kappa = basis$kappa,
         prior = if (is.null(correction)) "unconstrained" else "qr_constrained",
         adjacency = "Actual SPDE triangle edges: first ring is vertex one-hop, second ring is exactly two-hop",
         probe_selection = "Deterministic vertices with nonzero observation support",
         interpretation = "Mesh-vertex prior correlation, not element adjacency or fitted influence; no hard cutoff")
}
