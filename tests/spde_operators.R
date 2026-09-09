library(spatialGPCA)
library(Matrix)

check <- function(label, error, tol = 1e-8) {
    cat(sprintf("%-48s %.4g (tolerance %.4g)\n", label, error, tol))
    stopifnot(length(error) == 1L, is.finite(error), error < tol)
}
relative <- function(A, B) sqrt(sum((A - B)^2)) / max(1, sqrt(sum(B^2)))

# Isolate operator algebra from mesh construction with a synthetic sparse basis.
set.seed(715)
n <- 31L; m <- 12L; p <- 19L; k <- 4L
idx <- t(replicate(n - m, sample.int(m, 3L)))
val <- matrix(runif(3L * (n - m), 0.2, 1), n - m, 3L)
val <- val / rowSums(val)
B <- sparseMatrix(i = c(seq_len(m), rep((m + 1L):n, 3L)),
                  j = c(seq_len(m), as.vector(idx)), x = c(rep(1, m), as.vector(val)),
                  dims = c(n, m))
Q <- bandSparse(m, k = c(-1L, 0L, 1L),
                diagonals = list(rep(-0.3, m - 1L), runif(m, 1, 2), rep(-0.3, m - 1L)))
Q <- forceSymmetric(Q)
basis <- list(B = B, Q = Q, kappa = 1, tau = 1, practical_range_km = sqrt(8))
C <- matrix(rnorm(n * 3L), n, 3L)
W <- runif(n, 1, 9)
G <- matrix(rnorm(n * p), n, p)
X <- cbind(1, C) * sqrt(W)
Xq <- qr.Q(qr(X))
Y <- G * sqrt(W)
Y0 <- Y - Xq %*% crossprod(Xq, Y)
Z <- as.matrix(B) * sqrt(W)
Z0 <- Z - Xq %*% crossprod(Xq, Z)
S <- crossprod(Z0)

for (mode in c("profile", "qr")) {
    model <- svgpc_spde_prepare(basis, C, W, mode)
    model <- svgpc_spde_accumulate_matrix(model, G, block_size = 5L)
    N <- diag(m)
    if (mode == "qr") {
        cls <- qr(crossprod(Z, Xq))
        N <- qr.Q(cls, complete = TRUE)[, (cls$rank + 1L):m, drop = FALSE]
    }
    ZN <- Z0 %*% N
    SN <- crossprod(ZN)
    QN <- crossprod(N, as.matrix(Q) %*% N)
    for (lambda in c(0.2, 2, 50)) {
        tag <- paste(mode, lambda)
        E <- solve(SN + lambda * QN, crossprod(ZN, Y0))
        theta <- N %*% E
        F <- as.matrix(B) %*% theta
        H <- ZN %*% solve(SN + lambda * QN, t(ZN))
        H <- (H + t(H)) / 2
        fit <- svgpc_spde_fit(model, lambda = lambda, components = k,
                             operators = c("raw_F_PCA", "H2", "H", "2H-H2"))
        check(paste(tag, "coefficients"), relative(svgpc_spde_fitted(model, G, lambda, TRUE), theta))
        if (mode == "qr")
            check(paste(tag, "coefficient constraint"), max(abs(crossprod(X, (as.matrix(B) %*% theta) * sqrt(W)))))

        # Direct SNP-space PCA of the raw fitted matrix, in ordinary location L2.
        E <- eigen(crossprod(F) / p, symmetric = TRUE)
        PC <- sweep(F %*% E$vectors[, seq_len(k), drop = FALSE], 2L,
                    sqrt(p * E$values[seq_len(k)]), "/")
        out <- fit$outputs$raw_F_PCA
        check(paste(tag, "raw_F eigenvalues"), relative(out$eigenvalues, E$values[seq_len(k)]))
        check(paste(tag, "raw_F PC subspace"), relative(tcrossprod(out$pcs), tcrossprod(PC)))
        check(paste(tag, "raw_F orthonormality"), max(abs(crossprod(out$pcs) - diag(k))))
        check(paste(tag, "raw_F eigen-equation"),
              relative(F %*% crossprod(F, out$pcs) / p, sweep(out$pcs, 2L, out$eigenvalues, "*")))

        # Extract loadings directly from p x p SNP operators, then apply the same F.
        for (op in c("H2", "H", "2H-H2")) {
            O <- switch(op, H2 = H %*% H, H = H, `2H-H2` = 2 * H - H %*% H)
            E <- eigen(crossprod(Y0, O %*% Y0), symmetric = TRUE)
            raw <- F %*% E$vectors[, seq_len(k), drop = FALSE]
            weighted <- raw * sqrt(W)
            weighted <- weighted - Xq %*% crossprod(Xq, weighted)
            out <- fit$outputs[[op]]
            check(paste(tag, op, "eigenvalues"), relative(out$eigenvalues, E$values[seq_len(k)]))
            check(paste(tag, op, "raw score Gram"), relative(tcrossprod(out$raw_scores), tcrossprod(raw)))
            check(paste(tag, op, "weighted score Gram"), relative(tcrossprod(out$weighted_scores), tcrossprod(weighted)))
            check(paste(tag, op, "weighted fixed effects"), max(abs(crossprod(Xq, out$weighted_scores))))
        }
    }
}

# Full residual rank: very large finite lambda must approach the infinity limit.
n <- 9L
basis <- list(B = as(Diagonal(n), "generalMatrix"), Q = Diagonal(n, seq_len(n)),
              kappa = 1, tau = 1, practical_range_km = sqrt(8))
C <- matrix(rnorm(n * 2L), n, 2L)
W <- runif(n, 1, 4)
G <- matrix(rnorm(n * 11L), n, 11L)
for (mode in c("profile", "qr")) {
    model <- svgpc_spde_prepare(basis, C, W, mode)
    model <- svgpc_spde_accumulate_matrix(model, G, block_size = 4L)
    stopifnot(length(model$d) == model$r)
    for (lambda in c(1e160, 1e200, 1e300, Inf)) {
        check(paste(mode, "full-rank GCV", lambda), abs(svgpc_spde_objective(model, lambda, "GCV")))
        check(paste(mode, "full-rank REML", lambda), abs(svgpc_spde_objective(model, lambda, "REML")))
    }
}
check("correlation overflowing finite product",
      max(abs(svgpc_spde_correlation(c(0, 1e-300, 1e300), 1e300) - c(1, besselK(1, 1), 0))))
cat("All SPDE operator and numerical-boundary checks passed.\n")
