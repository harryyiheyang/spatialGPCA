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
Theta <- backsolve(Lq, matrix(rnorm(nrow(Q) * p), nrow(Q)))
Theta <- Theta / sd(as.vector(B %*% Theta)) * .6
Utrue <- cbind(rnorm(N), sample(c(-1, 1), N, TRUE))
X <- (B %*% Theta)[s, ] + Utrue %*% matrix(rnorm(2 * p, sd = .5), 2) + matrix(rnorm(n * p, sd = .3), n)[s, ] + matrix(rnorm(N * p), N)
X <- scale(X, scale = FALSE)
reader <- function(a, b) X[, a:b, drop = FALSE]
sc <- spatialGPCA:::.joint_scan

# 1. The scan's REML equals the dense individual-level REML (C = [1, U] fixed, sigma^2 profiled),
#    and U constant within locations reproduces the existing profile model's coefficients.
lams <- c(.3, 3, 30)
ZB <- B[s, ]
dense_reml <- function(U, lam) {
    C <- cbind(1, U); q <- ncol(C)
    Sig <- diag(N) + ZB %*% solve(Q, t(ZB)) / lam
    Si <- solve(Sig); CSC <- crossprod(C, Si %*% C)
    P <- Si - Si %*% C %*% solve(CSC, t(C) %*% Si)
    (N - q) * log(sum(X * (P %*% X)) / (p * (N - q))) + as.numeric(determinant(Sig)$modulus) +
        as.numeric(determinant(CSC)$modulus) - as.numeric(determinant(crossprod(C))$modulus)
}
res <- sc(model, reader, p, 64L, Utrue, lams, 2L)
stopifnot(max(abs(res$reml - sapply(lams, function(l) dense_reml(Utrue, l)))) < 1e-6)
# Invariance to U -> U T + 1 v'
res_t <- sc(model, reader, p, 64L, Utrue %*% matrix(c(2, 1, 0, 3), 2) + 5, lams, 2L)
stopifnot(max(abs(res$reml - res_t$reml)) < 1e-6)
Ubar <- matrix(rnorm(n * k), n); Uc <- Ubar[s, ]
xbar <- rowsum(X, s, reorder = TRUE) / w
basis <- list(B = model$B, Q = model$Q, mesh = mesh, kappa = kappa, tau = 1)
old <- svgpc_spde_accumulate_matrix(svgpc_spde_prepare(basis, Ubar, w, "profile"), xbar)
th_old <- svgpc_spde_fitted(old, xbar[, 1:5], lams[2], coefficients = TRUE)
C <- cbind(1, Uc); M <- rbind(cbind(crossprod(C), crossprod(C, ZB)),
                             cbind(crossprod(ZB, C), crossprod(ZB) + lams[2] * Q))
sol <- solve(M, rbind(crossprod(C, X), crossprod(ZB, X)))
stopifnot(max(abs(sol[-(1:(k + 1)), 1:5] - th_old)) < 1e-6 * max(abs(th_old)))

# 1b. The U proposal is a stationary point condition of the REML: finite-difference check of
#     the gradient at the proposal's defining equation, through the closed form below.
U0 <- Utrue + matrix(rnorm(N * k, sd = .3), N)
r0 <- sc(model, reader, p, 64L, U0, 3, 1L, update_U = TRUE)
grad_fd <- function(U, i, j, h = 1e-5) {
    Up <- U; Up[i, j] <- Up[i, j] + h; Um <- U; Um[i, j] <- Um[i, j] - h
    (dense_reml(Up, 3) - dense_reml(Um, 3)) / (2 * h)
}
# Analytic gradient: (2 / sigma2_p) (U A'A - R0 A) + 2 [(C - ZBH) G^-1 - C (C'C)^-1]_U, sigma2_p = L/df.
df <- N - k - 1; sig <- r0$objective / df
grad_an <- 2 / sig * (U0 %*% r0$AtA - r0$Uacc) + 2 * (r0$Uacc - r0$Uprop %*% r0$AtA) / sig
for (ij in list(c(1, 1), c(17, 2), c(400, 1))) {
    stopifnot(abs(grad_fd(U0, ij[1], ij[2]) - grad_an[ij[1], ij[2]]) < 1e-4 * max(1, abs(grad_an[ij[1], ij[2]])))
}

# 2. U varying within locations: objective equals the dense individual-level minimum.
res2 <- sc(model, reader, p, 50L, Utrue, lams[2], 1L, update_U = TRUE)
C <- cbind(1, Utrue); ZB <- B[s, ]
M <- rbind(cbind(crossprod(C), crossprod(C, ZB)), cbind(crossprod(ZB, C), crossprod(ZB) + lams[2] * Q))
rhs <- rbind(crossprod(C, X), crossprod(ZB, X))
dense <- sum(X^2) - sum(rhs * solve(M, rhs))
stopifnot(abs(res2$objective - dense) / dense < 1e-9)
# U step equals the dense closed form.
sol <- solve(M, rhs); A <- t(sol[2:3, ]); R <- X - matrix(1, N, 1) %*% sol[1, , drop = FALSE] - ZB %*% sol[-(1:3), ]
stopifnot(max(abs(res2$Uacc - R %*% A)) < 1e-8 * max(abs(res2$Uacc)))

# 2b. Subspace iteration converges to the dense generalized eigenproblem X_c X_c' Om^-1 U = U D.
lam <- 3; Si <- solve(crossprod(ZB) + lam * Q)
Om <- diag(N) - ZB %*% Si %*% t(ZB)
o1 <- Om %*% rep(1, N); P1 <- diag(N) - outer(rep(1, N), drop(o1)) / sum(o1)
Xc <- P1 %*% X
eo <- eigen(Om, symmetric = TRUE); Oh <- eo$vectors %*% (sqrt(pmax(eo$values, 0)) * t(eo$vectors))
Ud <- svd(Oh %*% Xc, nu = k)$u
YS <- matrix(rnorm(N * 8), N)
for (it in 1:40) { rr <- sc(model, reader, p, 100L, matrix(0, N, 0), lam, 1L, YS = YS); YS <- rr$sub$next_basis }
Ur <- Oh %*% rr$sub$ritz[, 1:k]
stopifnot(min(svd(crossprod(qr.Q(qr(Ur)), Ud))$d) > 1 - 1e-6)

# 3. Joint fit: monotone objective at fixed lambda, recovers U and finds a finite lambda.
fit <- svgpc_joint_fit(model, k = 2, r = 3, reader = reader, snps = p, verbose = FALSE, max_iter = 60)
tr <- fit$trace
stopifnot(all(is.finite(tr$reml)), fit$lambda > 0, is.finite(fit$lambda))
stopifnot(tail(tr$reml, 1) <= min(tr$reml[-1]) + 1e-4, fit$reml <= min(tr$reml) + 1e-6)
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

# 5. Threaded C++ kernels and BLAS give the same fit.
f1 <- svgpc_joint_fit(model, k = 2, r = 3, reader = reader, snps = p, verbose = FALSE, kernel = "cpp", max_iter = 8)
f2 <- svgpc_joint_fit(model, k = 2, r = 3, reader = reader, snps = p, verbose = FALSE, kernel = "blas", max_iter = 8)
stopifnot(abs(f1$reml - f2$reml) < 1e-8 * abs(f1$reml), abs(log(f1$lambda / f2$lambda)) < 1e-8)

# 6. Fused 2-bit BED kernel: same scan statistics and the same fit as decoding to doubles.
dir <- tempfile("jointbed"); dir.create(dir); prefix <- file.path(dir, "sim")
pb <- 120; afb <- runif(pb, .05, .5)
G <- sapply(afb, function(f) rbinom(N, 2, f))
G[sample(length(G), 200)] <- NA
code <- ifelse(is.na(G), 1L, ifelse(G == 2, 0L, ifelse(G == 1, 2L, 3L)))  # dosage of A1 = COUNTED
strideb <- ceiling(N / 4); bytes <- raw(0)
for (j in seq_len(pb)) {
    cj <- c(code[, j], integer(4 * strideb - N))
    bytes <- c(bytes, as.raw(colSums(matrix(cj, 4) * c(1L, 4L, 16L, 64L))))
}
writeBin(c(as.raw(c(0x6c, 0x1b, 0x01)), bytes), paste0(prefix, ".bed"))
write.table(cbind(table$FID, table$IID, 0, 0, 0, -9), paste0(prefix, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(1, paste0("rs", 1:pb), 0, 1:pb, "A", "G"), paste0(prefix, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(`#CHROM` = 1, ID = paste0("rs", 1:pb), REF = "G", ALT = "A", ALT_FREQS = afb, check.names = FALSE),
            paste0(prefix, ".afreq"), quote = FALSE, row.names = FALSE, sep = "\t")
mb <- svgpc_joint_prepare(mesh, table, kappa, bed = prefix)
bd <- spatialGPCA:::.joint_bed(mb$spatial, NULL)
rd <- function(a, b) spatialGPCA:::cpp_joint_bed_block(bd$path, bd$file_samples, bd$rows, a - 1L, b - a + 1L, bd$af, 2L)
bs <- function(a, b, Yt, q, RtT) spatialGPCA:::cpp_joint_bed_scan(bd$path, bd$file_samples, bd$rows, (a:b) - 1L,
                                                                  bd$af[a:b], mb$location, mb$n, Yt, q, RtT, 2L)
Uq <- matrix(rnorm(N * 2), N); YSq <- matrix(rnorm(N * 5), N); YVq <- qr.Q(qr(matrix(rnorm(mb$n * 4), mb$n)))
r1 <- sc(mb, rd, pb, 50L, Uq, lams, 2L, YS = YSq, YV = YVq, YP = YVq)
r2 <- sc(mb, NULL, pb, 50L, Uq, lams, 2L, YS = YSq, YV = YVq, YP = YVq, bedscan = bs)
stopifnot(max(abs(r1$reml - r2$reml)) < 1e-8, max(abs(r1$sub$next_basis - r2$sub$next_basis)) < 1e-8 * max(abs(r1$sub$next_basis)),
          max(abs(r1$aV$Y - r2$aV$Y)) < 1e-8 * max(abs(r1$aV$Y)), max(abs(r1$aP$Y - r2$aP$Y)) < 1e-8 * max(abs(r1$aP$Y)))
fb1 <- svgpc_joint_fit(mb, k = 2, r = 3, verbose = FALSE, max_iter = 8, threads = 2)
fb2 <- svgpc_joint_fit(mb, k = 2, r = 3, reader = rd, snps = pb, verbose = FALSE, max_iter = 8, kernel = "cpp", threads = 2)
stopifnot(fb1$kernel == "fused", abs(fb1$reml - fb2$reml) < 1e-8 * abs(fb1$reml), abs(log(fb1$lambda / fb2$lambda)) < 1e-8)
cat("fused kernel checks passed\n")

# 7. Variant subsets (index vector) match a reader over the same columns; the warm start reaches the
#    same fixed point as a cold start.
sub <- seq(1, pb, by = 3)
fs1 <- svgpc_joint_fit(mb, k = 2, r = 3, snps = sub, verbose = FALSE, max_iter = 6, threads = 2, warm = 1)
rsub <- function(a, b) rd(1, pb)[, sub[a:b], drop = FALSE]
fs2 <- svgpc_joint_fit(mb, k = 2, r = 3, reader = rsub, snps = length(sub), verbose = FALSE, max_iter = 6, kernel = "cpp", threads = 2)
stopifnot(abs(fs1$reml - fs2$reml) < 1e-8 * abs(fs1$reml))
fw <- svgpc_joint_fit(mb, k = 2, r = 3, verbose = FALSE, max_iter = 40, threads = 2, warm = 2, tol = 1e-6)
fc <- svgpc_joint_fit(mb, k = 2, r = 3, verbose = FALSE, max_iter = 40, threads = 2, warm = 1, tol = 1e-6)
print(c(fw$lambda, fc$lambda, fw$reml, fc$reml)); stopifnot(!is.null(fw$warm_trace), abs(fw$reml - fc$reml) < 1e-2,
          min(cancor(fw$U, fc$U)$cor) > .99)
cat("subset and warm-start checks passed:", nrow(fw$warm_trace), "+", nrow(fw$trace), "vs", nrow(fc$trace), "scans\n")
