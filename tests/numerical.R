library(spatialGPCA)
# Fixed-parameter calls are numerical test oracles, not user-facing inputs.
.fit <- function(model, method = c("REML", "GCV"), operators = "raw_F_PCA", components = 50L, lambda = NULL) {
    spatialGPCA:::cpp_fit(model, match.arg(method), operators, as.integer(components),
                    if (is.null(lambda)) NaN else lambda)
}
.fitted <- function(model, G, lambda, coefficients = FALSE) spatialGPCA:::cpp_fitted(model, G, lambda, coefficients)
set.seed(404)
n <- 200L; p <- 500L; m <- 50L
xy <- matrix(runif(2 * n), n, 2)
knots <- xy[sample.int(n, m), ]
C <- matrix(rnorm(n * 15), n, 15)
Ni <- sample(2:30, n, replace = TRUE)
rho <- 0.075
kernel <- function(a, b) exp(-sqrt(outer(a[, 1], b[, 1], "-")^2 +
                                      outer(a[, 2], b[, 2], "-")^2) / rho)
B <- kernel(xy, knots); Q <- kernel(knots, knots)
G <- B %*% matrix(rnorm(m * p), m, p) + C %*% matrix(rnorm(15 * p), 15, p) +
    matrix(rnorm(n * p), n, p) / sqrt(Ni)
errors <- numeric()
check <- function(name, error, tol = 1e-8) {
    errors[name] <<- error
    cat(sprintf("%-45s %.4g (tolerance %.4g)\n", name, error, tol))
    stopifnot(is.finite(error), error < tol)
}
relative <- function(a, b) sqrt(sum((a - b)^2)) / max(sqrt(sum(b^2)), 1)
suberr <- function(a, b) {
    qa <- qr.Q(qr(a)); qb <- qr.Q(qr(b))
    max(abs(1 - svd(crossprod(qa, qb), nu = 0, nv = 0)$d))
}
new_model <- function(mode = "qr", counts = Ni) {
    spatialGPCA:::.prepare_xy(xy, C, counts, knots, rho)
}
for (mode in "qr") {
    model <- new_model(mode)
    s <- spatialGPCA:::cpp_snapshot(model)
    BN <- B %*% s$N; QN <- crossprod(s$N, Q %*% s$N)
    X <- cbind(1, C) * sqrt(Ni); Y <- G * sqrt(Ni)
    Z <- BN * sqrt(Ni); Y0 <- Y - s$Xq %*% crossprod(s$Xq, Y)
    Z0 <- Z - s$Xq %*% crossprod(s$Xq, Z)
    S <- crossprod(Z0); Tcross <- crossprod(s$T, Y0)
    check(paste(mode, "T orthonormal"), max(abs(crossprod(s$T) - diag(ncol(s$T)))))
    if (mode == "qr") check("coefficient constraint", max(abs(crossprod(X, Z))))
    for (lambda in c(0.02, 2, 200)) {
        theta <- s$N %*% solve(S + lambda * QN, crossprod(Z0, Y0))
        fit <- B %*% theta
        check(paste(mode, "A coefficients", lambda), relative(.fitted(model, G, lambda, TRUE), theta))
        check(paste(mode, "A raw fitted", lambda), relative(.fitted(model, G, lambda), fit))
    }
    svgpc_accumulate_matrix(model, G, 37L)
    stat <- svgpc_statistics(model)
    check(paste(mode, "D covariance"), relative(stat$covariance_sum, tcrossprod(Tcross)))
    check(paste(mode, "D energy"), abs(stat$total_energy / sum(Y0^2) - 1))
    whole <- new_model(mode); svgpc_accumulate_matrix(whole, G, p)
    check(paste(mode, "D block vs whole"), relative(stat$covariance_sum, svgpc_statistics(whole)$covariance_sum))
    for (lambda in c(0.2, 20)) {
        V <- diag(n) + Z %*% solve(QN, t(Z)) / lambda
        Vinv <- solve(V)
        P <- Vinv - Vinv %*% X %*% solve(crossprod(X, Vinv %*% X), t(X) %*% Vinv)
        logdet <- function(a) 2 * sum(log(diag(chol(a))))
        direct <- logdet(V) + logdet(crossprod(X, Vinv %*% X)) - logdet(crossprod(X)) +
            (n - ncol(X)) * log(sum(Y * (P %*% Y)) / sum(Y0^2))
        out <- .fit(model, lambda = lambda, components = 5)
        check(paste(mode, "B full REML", lambda), abs(out$objective_relative - direct), 2e-8)
        h <- stat$d / (stat$d + lambda)
        resid <- Y0 - s$T %*% (Tcross * as.vector(h))
        gcv <- log(sum(resid^2) / sum(Y0^2)) - 2 * log(1 - sum(h) / (n - ncol(X)))
        check(paste(mode, "GCV", lambda), abs(.fit(model, "GCV", lambda = lambda)$objective_relative - gcv))
        F <- .fitted(model, G, lambda)
        direct_pc <- svd(F, nu = 5, nv = 0)
        check(paste(mode, "raw F eigenvalues", lambda), relative(out$outputs$raw_F_PCA$eigenvalues,
                                                                  direct_pc$d[1:5]^2 / p))
        check(paste(mode, "raw F subspace", lambda), suberr(out$outputs$raw_F_PCA$pcs, direct_pc$u))
        other <- .fit(model, lambda = lambda, operators = c("H2", "H", "2H-H2"), components = 5)
        for (op in names(other$outputs)) {
            w <- switch(op, H2 = h^2, H = h, `2H-H2` = 2 * h - h^2)
            right <- svd(Tcross * as.vector(sqrt(w)), nu = 0, nv = 5)$v
            check(paste(mode, op, "gene subspace"), suberr(other$outputs[[op]]$raw_scores, F %*% right))
        }
    }
    f0 <- .fitted(model, G, 0)
    check(paste(mode, "zero lambda limit"), relative(.fitted(model, G, 1e-9), f0))
    check(paste(mode, "infinite lambda limit"), max(abs(.fitted(model, G, Inf))))
    lo <- .fit(model, lambda = 0.02, components = 5)
    hi <- .fit(model, lambda = 200, components = 5)
    change <- suberr(lo$outputs$raw_F_PCA$pcs, hi$outputs$raw_F_PCA$pcs)
    cat(mode, "C changed subspace 1-min canonical correlation:", change, "\n")
    stopifnot(change > 0.001)
    reml <- .fit(model, components = 5)
    gcv <- .fit(model, "GCV", components = 5)
    stopifnot(is.finite(reml$lambda), is.finite(gcv$lambda))
    large_w <- new_model(mode, 7 * Ni); svgpc_accumulate_matrix(large_w, G)
    fit7 <- .fit(large_w, components = 5)
    check(paste(mode, "E W scaling lambda"), abs(fit7$lambda / reml$lambda / 7 - 1), 1e-5)
    check(paste(mode, "E W scaling fit"), relative(.fitted(large_w, G, 7 * reml$lambda),
                                                   .fitted(model, G, reml$lambda)))
    other_w <- new_model(mode, rev(Ni)); svgpc_accumulate_matrix(other_w, G)
    stopifnot(relative(.fitted(other_w, G, reml$lambda), .fitted(model, G, reml$lambda)) > 1e-3)
    cache <- tempfile(fileext = ".rds"); svgpc_save(model, cache)
    loaded <- svgpc_load(cache); unlink(cache)
    again <- .fit(loaded, components = 5)
    check(paste(mode, "cache lambda"), abs(again$lambda - reml$lambda))
    check(paste(mode, "cache eigenvalues"), max(abs(again$outputs$raw_F_PCA$eigenvalues - reml$outputs$raw_F_PCA$eigenvalues)))
    stopifnot(inherits(try(svgpc_accumulate_matrix(model, G), silent = TRUE), "try-error"))
    fixed_only <- cbind(1, C) %*% matrix(rnorm(16 * p), 16, p)
    check(paste(mode, "F fixed-only removal"), sqrt(sum(.fitted(model, fixed_only, 1)^2) / sum(fixed_only^2)))
}
model <- new_model("qr")
s <- spatialGPCA:::cpp_snapshot(model)
truth <- B %*% s$N %*% matrix(rnorm(ncol(s$N) * 5), ncol(s$N), 5)
G2 <- truth %*% matrix(rnorm(5 * p), 5, p) + C %*% matrix(rnorm(15 * p), 15, p) +
    matrix(rnorm(n * p, sd = 0.05), n, p) / sqrt(Ni)
svgpc_accumulate_matrix(model, G2)
pc <- .fit(model, components = 5)$outputs$raw_F_PCA$pcs
check("G geographic subspace recovery", suberr(pc, truth), 0.01)
cat("All numerical tests passed; maximum algebraic error:", max(errors[!grepl("recovery", names(errors))]), "\n")

# Every location used as a knot: stable internal zero-boundary selection.
xy <- cbind(seq(0, 1, length.out = 12), rep(0, 12))
model <- spatialGPCA:::.prepare_xy(xy, matrix(numeric(), 12, 0), rep(1, 12), xy, .2)
s <- spatialGPCA:::cpp_snapshot(model)
G <- s$T %*% diag(as.vector(s$d))
svgpc_accumulate_matrix(model, G)
st <- svgpc_statistics(model)
a <- diag(st$covariance_sum)
stopifnot(length(st$d) == st$residual_df, st$b == 0)
for (method in c("REML", "GCV")) {
    expected <- if (method == "REML") sum(log(st$d)) + st$residual_df * log(sum(a / st$d) / st$total_energy) else
        log(sum(a / st$d^2) / st$total_energy) - 2 * log(sum(1 / st$d) / st$residual_df)
    zero <- .fit(model, method, lambda = 0)
    check(paste(method, "full-rank zero objective"), abs(zero$objective_relative - expected))
    near <- .fit(model, method, lambda = 1e-9)
    check(paste(method, "full-rank zero continuity"), abs(near$objective_relative - expected), 1e-6)
    selected <- svgpc_select_lambda(model, method)
    stopifnot(selected$boundary == "zero", selected$lambda == 0)
}
cat("Full-rank REML/GCV boundary tests passed.\n")
