library(spatialGPCA)
library(Matrix)
set.seed(11)
n <- 200; N <- 1500; p <- 300; k <- 2
ll <- data.frame(Lat = 52 + runif(n, -2, 2), Lon = -1.5 + runif(n, -2, 2))
loc <- c(1:n, sample(n, N - n, TRUE))
table <- data.frame(FID = paste0("f", 1:N), IID = paste0("i", 1:N),
                    Lat = ll$Lat[loc], Lon = ll$Lon[loc], PC1 = rnorm(N), PC2 = rnorm(N))
spatial <- svgpc_spde_locations(table)
xy <- as.matrix(spatial$locations[, c("x", "y")])
# Structured triangle mesh covering the projected locations (no RTriangle needed).
gx <- seq(min(xy[, 1]) - 60, max(xy[, 1]) + 60, length.out = 7)
gy <- seq(min(xy[, 2]) - 60, max(xy[, 2]) + 60, length.out = 7)
V <- as.matrix(expand.grid(gx, gy)); id <- matrix(seq_len(nrow(V)), 7, 7)
tv <- do.call(rbind, lapply(1:6, function(i) do.call(rbind, lapply(1:6, function(j)
    rbind(c(id[i, j], id[i + 1, j], id[i + 1, j + 1]), c(id[i, j], id[i + 1, j + 1], id[i, j + 1]))))))
mesh <- structure(list(xy = V, tv = tv, locations = spatial$locations, units = "km",
                       coordinate_system = spatial$projection), class = "svgpc_spde_mesh")
kappa <- svgpc_spde_scale(120)$kappa
model <- svgpc_joint_prepare(mesh, table, kappa)
s <- model$location; w <- model$counts; B <- as.matrix(model$B); Q <- as.matrix(model$Q)
# Simulated standardized genotypes: smooth field + two nonspatial factors + noise.
Lq <- chol(Q + diag(1e-8, nrow(Q)))
Theta <- backsolve(Lq, matrix(rnorm(nrow(Q) * p), nrow(Q))) * 3
Utrue <- cbind(rnorm(N), sample(c(-1, 1), N, TRUE))
X <- (B %*% Theta)[s, ] + Utrue %*% matrix(rnorm(2 * p, sd = .5), 2) + matrix(rnorm(n * p, sd = .3), n)[s, ] + matrix(rnorm(N * p), N)
X <- scale(X, scale = FALSE)
reader <- function(a, b) X[, a:b, drop = FALSE]
sc <- spatialGPCA:::.joint_scan

# 1. U constant within locations reduces to the existing profile model.
Ubar <- matrix(rnorm(n * k), n); Uc <- Ubar[s, ]
lams <- c(.3, 3, 30)
res <- sc(model, reader, p, 64L, Uc, lams, 2L)
xbar <- rowsum(X, s, reorder = TRUE) / w
basis <- list(B = model$B, Q = model$Q, mesh = mesh, kappa = kappa, tau = 1)
old <- svgpc_spde_accumulate_matrix(svgpc_spde_prepare(basis, Ubar, w, "profile"), xbar)
oo <- sapply(lams, function(l) svgpc_spde_objective(old, l))
stopifnot(max(abs(diff(res$reml) - diff(oo))) < 1e-6)
th_old <- svgpc_spde_fitted(old, xbar[, 1:5], lams[2], coefficients = TRUE)
C <- cbind(1, Uc); M <- rbind(cbind(crossprod(C), crossprod(C, B[s, ])),
                             cbind(crossprod(B[s, ], C), crossprod(B[s, ]) + lams[2] * Q))
sol <- solve(M, rbind(crossprod(C, X), crossprod(B[s, ], X)))
stopifnot(max(abs(sol[-(1:(k + 1)), 1:5] - th_old)) < 1e-6 * max(abs(th_old)))

# 2. U varying within locations: objective equals the dense individual-level minimum.
res2 <- sc(model, reader, p, 50L, Utrue, lams[2], 1L)
C <- cbind(1, Utrue); ZB <- B[s, ]
M <- rbind(cbind(crossprod(C), crossprod(C, ZB)), cbind(crossprod(ZB, C), crossprod(ZB) + lams[2] * Q))
rhs <- rbind(crossprod(C, X), crossprod(ZB, X))
dense <- sum(X^2) - sum(rhs * solve(M, rhs))
stopifnot(abs(res2$objective - dense) / dense < 1e-9)
# U step equals the dense closed form.
sol <- solve(M, rhs); A <- t(sol[2:3, ]); R <- X - matrix(1, N, 1) %*% sol[1, , drop = FALSE] - ZB %*% sol[-(1:3), ]
stopifnot(max(abs(res2$Uacc - R %*% A)) < 1e-8 * max(abs(res2$Uacc)))

# 3. Joint fit: monotone objective at fixed lambda, recovers U and finds a finite lambda.
fit <- svgpc_joint_fit(model, k = 2, r = 3, reader = reader, snps = p, verbose = FALSE, max_iter = 30)
tr <- fit$trace
stopifnot(all(is.finite(tr$reml)), fit$lambda > 0, is.finite(fit$lambda))
cc <- cancor(fit$U, Utrue)$cor
stopifnot(min(cc) > .9)
Vi <- fit$V[s, ]
stopifnot(max(abs(crossprod(Vi) - diag(ncol(Vi)))) < 1e-6)
Vo <- fit$V_orthogonal[s, ]
stopifnot(max(abs(crossprod(cbind(1, fit$U), Vo))) < 1e-6 * N)
print(fit); print(round(cc, 4))

# 4. BED decoding matches a direct R decoder, with missing calls set to 2f.
folder <- system.file("extdata", package = "spatialGPCA")
ff <- read.table(file.path(folder, "frequencies.afreq"), header = TRUE, comment.char = "", check.names = FALSE)
fam <- read.table(file.path(folder, "missing.fam"), colClasses = "character")
ns <- nrow(fam); stride <- ceiling(ns / 4); af <- ff$ALT_FREQS
raw <- readBin(file.path(folder, "missing.bed"), "raw", n = 3 + stride * nrow(ff))
rows <- sample(ns, 50)
G <- spatialGPCA:::cpp_joint_bed_block(file.path(folder, "missing.bed"), ns, rows, 4L, 20L, af, 2L)
for (j in 1:20) for (ii in seq_along(rows)) {
    i <- rows[ii]; v <- 4 + j
    code <- bitwAnd(bitwShiftR(as.integer(raw[3 + (v - 1) * stride + (i - 1) %/% 4 + 1]), 2L * ((i - 1L) %% 4L)), 3L)
    g <- c(2, NA, 1, 0)[code + 1L]; if (is.na(g)) g <- 2 * af[v]
    stopifnot(abs(G[ii, j] - (g - 2 * af[v]) / sqrt(2 * af[v] * (1 - af[v]))) < 1e-12)
}
