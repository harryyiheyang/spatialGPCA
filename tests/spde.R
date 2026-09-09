library(spatialGPCA)
library(Matrix)
if (requireNamespace("RTriangle", quietly = TRUE)) {
    set.seed(731)
    xy <- cbind(runif(180, 0, 10), runif(180, 0, 8))
    W <- sample(1:20, nrow(xy), TRUE)
    C <- matrix(rnorm(nrow(xy) * 20), nrow(xy), 20)
    G <- matrix(rnorm(nrow(xy) * 47), nrow(xy), 47)
    mesh <- svgpc_spde_mesh(xy, max_area = 4, buffer = 3)
    scale <- svgpc_spde_scale(4)
    basis <- svgpc_spde_basis(mesh, xy, scale$kappa)
    B <- as.matrix(basis$B)
    Q <- as.matrix(basis$Q)
    stopifnot(max(Matrix::rowSums(basis$B != 0)) <= 3,
              max(abs(rowSums(B) - 1)) < 1e-12,
              max(abs(B %*% mesh$xy - xy)) < 1e-8,
              max(abs(basis$stiffness %*% rep(1, ncol(B)))) < 1e-10,
              abs(svgpc_spde_correlation(4, scale$kappa) - exp(-1)) < 1e-10)
    errors <- list()
    for (mode in c("profile", "qr")) {
        R1 <- svgpc_spde_prepare(basis, C, W, mode)
        R2 <- svgpc_spde_accumulate_matrix(R1, G, block_size = 11)
        R3 <- svgpc_spde_accumulate_reader(R1, function(start, end) G[, start:end, drop = FALSE],
                                          ncol(G), block_size = 7)
        stopifnot(max(abs(R2$cov - R3$cov)) / max(abs(R2$cov)) < 1e-10)
        X <- cbind(1, C) * sqrt(W)
        Xq <- qr.Q(qr(X))
        Y <- G * sqrt(W)
        Z <- B * sqrt(W)
        N <- diag(ncol(B))
        if (mode == "qr") {
            A <- qr(crossprod(Z, Xq))
            N <- qr.Q(A, complete = TRUE)[, (A$rank + 1):ncol(B), drop = FALSE]
        }
        Zs <- Z %*% N
        Qs <- crossprod(N, Q %*% N)
        Y0 <- Y - Xq %*% crossprod(Xq, Y)
        Z0 <- Zs - Xq %*% crossprod(Xq, Zs)
        S <- crossprod(Z0)
        for (lam in c(0.01, 1, 100)) {
            theta <- solve(S + lam * Qs, crossprod(Z0, Y0))
            F <- B %*% N %*% theta
            F1 <- svgpc_spde_fitted(R2, G, lam)
            F2 <- B %*% svgpc_spde_fitted(R2, G, lam, coefficients = TRUE)
            V <- diag(nrow(B)) + Zs %*% solve(Qs, t(Zs)) / lam
            Vi <- solve(V)
            H <- crossprod(X, Vi %*% X)
            P <- Vi - Vi %*% X %*% solve(H, t(X) %*% Vi)
            obj <- as.numeric(determinant(V, logarithm = TRUE)$modulus +
                       determinant(H, logarithm = TRUE)$modulus -
                       determinant(crossprod(X), logarithm = TRUE)$modulus) +
                   R2$r * log(sum(Y * (P %*% Y)) / sum(Y0^2))
            obj1 <- svgpc_spde_objective(R2, lam)
            fit <- svgpc_spde_fit(R2, lambda = lam, components = 5,
                                  operators = c("raw_F_PCA", "H2", "H", "2H-H2"))
            E <- eigen(tcrossprod(F) / ncol(G), symmetric = TRUE)
            PC <- fit$outputs$raw_F_PCA$pcs
            err <- c(fitted = max(abs(F - F1)), coefficients = max(abs(F2 - F1)),
                     reml = abs(obj - obj1), eigen = max(abs(E$values[1:5] - fit$outputs$raw_F_PCA$eigenvalues)),
                     pc_space = max(abs(tcrossprod(E$vectors[, 1:5]) - tcrossprod(PC))))
            stopifnot(max(err) < 1e-6)
            errors[[paste(mode, lam)]] <- err
        }
        opt <- svgpc_spde_fit(R2, components = 3)
        stopifnot(is.finite(opt$objective_relative), opt$p == ncol(G))
        stopifnot(svgpc_spde_objective(R2, Inf) == 0,
                  max(abs(svgpc_spde_fitted(R2, G, Inf))) == 0)
        if (mode == "qr") stopifnot(max(abs(crossprod(Xq, Z %*% R2$V))) < 1e-7)
    }
    print(do.call(rbind, errors))
}
