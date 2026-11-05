# Joint decomposition X = 1 mu' + U A' + Z B Theta + E (individual level).
# U: k unpenalized individual factors; Theta: SPDE mesh coefficients with
# precision lambda Q(kappa)/sigma^2 shared by all SNPs. Only sparse mesh algebra:
# no m x m dense matrix, no N x p or m x p storage. One genotype scan per iteration.

.chol_logdet <- function(f) {
    d <- tryCatch(Matrix::determinant(f, logarithm = TRUE, sqrt = FALSE),
                  error = function(e) NULL)
    if (!is.null(d)) return(as.numeric(d$modulus))
    2 * sum(log(Matrix::diag(Matrix::expand1(f, "L"))))
}

.joint_bed <- function(spatial, frequencies) {
    prefix <- spatial$bed
    if (is.null(prefix)) stop("Supply bed when preparing the joint model, or a reader when fitting")
    path <- paste0(prefix, ".bed")
    if (!file.exists(path)) stop("Missing genotype file: ", path)
    sm <- .metadata_table(paste0(prefix, ".fam"), FALSE)
    if (ncol(sm) != 6L) stop("FAM must have six columns")
    names(sm)[1:2] <- c("FID", "IID")
    vm <- .metadata_table(paste0(prefix, ".bim"), FALSE)
    if (ncol(vm) != 6L) stop("BIM must have six columns")
    names(vm) <- c("CHROM", "ID", "CM", "POS", "COUNTED", "OTHER")
    if (anyDuplicated(vm$ID) || any(vm$ID %in% c("", "."))) stop("Variant IDs must be unique and nonmissing")
    if (any(grepl(",", vm$COUNTED, fixed = TRUE)) || any(vm$COUNTED == vm$OTHER))
        stop("This version requires biallelic variants")
    idx <- match(.sample_key(spatial$table), .sample_key(sm))
    if (anyNA(idx)) stop("FAM changed after preparation; prepare the joint model again")
    if (is.null(frequencies)) frequencies <- paste0(prefix, ".afreq")
    if (is.character(frequencies) && length(frequencies) == 1L)
        frequencies <- .metadata_table(frequencies)
    if (all(c("ID", "REF", "ALT", "ALT_FREQS") %in% names(frequencies))) {
        ff <- data.frame(ID = frequencies$ID, A1 = frequencies$ALT,
                         A2 = frequencies$REF, AF = frequencies$ALT_FREQS)
    } else if (all(c("SNP", "A1", "A2", "MAF") %in% names(frequencies))) {
        ff <- data.frame(ID = frequencies$SNP, A1 = frequencies$A1,
                         A2 = frequencies$A2, AF = frequencies$MAF)
    } else ff <- frequencies
    .required(ff, c("ID", "A1", "A2", "AF"), "frequencies")
    fidx <- match(vm$ID, ff$ID)
    if (anyNA(fidx)) stop("Missing external allele frequencies for file variants")
    ff <- ff[fidx, , drop = FALSE]
    same <- ff$A1 == vm$COUNTED & ff$A2 == vm$OTHER
    flip <- ff$A2 == vm$COUNTED & ff$A1 == vm$OTHER
    if (any(!same & !flip)) stop("Alleles disagree between frequency and genotype metadata")
    af <- as.numeric(as.character(ff$AF))
    af[flip] <- 1 - af[flip]
    if (any(!is.finite(af)) || any(af <= 0 | af >= 1))
        stop("External allele frequencies must be strictly between zero and one; filter monomorphic SNPs first")
    list(path = normalizePath(path), file_samples = nrow(sm), rows = as.integer(idx),
         af = af, snps = nrow(vm), flipped = sum(flip))
}

svgpc_joint_prepare <- function(mesh, table, kappa, bed = NULL, tau = 1) {
    if (!inherits(mesh, "svgpc_spde_mesh")) stop("Select a mesh from svgpc_spde_mesh() first")
    projection <- mesh$coordinate_system
    if (is.null(projection$origin))
        stop("Build the mesh from svgpc_spde_locations(table) so the projection is retained")
    spatial <- svgpc_spde_locations(table, bed = bed, projection = projection)
    basis <- svgpc_spde_basis(mesh, spatial$locations[, c("x", "y")], kappa, tau)
    w <- as.numeric(spatial$locations$count)
    B <- methods::as(basis$B, "CsparseMatrix")
    BtWB <- Matrix::forceSymmetric(Matrix::crossprod(B, Matrix::Diagonal(length(w), w) %*% B))
    Q <- Matrix::forceSymmetric(basis$Q)
    fq <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)
    out <- list(spatial = spatial, B = B, Q = Q, BtWB = BtWB, counts = w,
                location = as.integer(spatial$sample_location), N = nrow(spatial$table),
                n = length(w), m = ncol(B), kappa = basis$kappa, tau = basis$tau,
                logdetQ = .chol_logdet(fq), mesh = mesh)
    class(out) <- "svgpc_joint_model"
    out
}

print.svgpc_joint_model <- function(x, ...) {
    cat("Joint U + SPDE model:", x$N, "individuals at", x$n, "locations;", x$m, "mesh vertices\n")
    cat("kappa:", x$kappa, "km^-1\n")
    invisible(x)
}

# One pass over all SNPs at fixed U and a small lambda grid.
.joint_scan <- function(model, reader, p, block, U, lambdas, center,
                        update_U = FALSE, YU = NULL, YV = NULL, YP = NULL, fac0 = NULL, threads = 0L,
                        YS = NULL, use_blas = FALSE, bedscan = NULL) {
    # bedscan(start, end, Yt, q, RtT): fused 2-bit kernel returning location sums, x'x and X'Y and
    # adding (X T)' into RtT in place; no N x block double matrix is formed.
    if (!is.null(bedscan) && (update_U || !is.null(YU))) stop("The fused kernel does not support update_U or YU")
    s <- model$location; w <- model$counts; B <- model$B; N <- model$N; n <- model$n; m <- model$m
    C <- cbind(1, U)
    q <- ncol(C); k <- q - 1L
    ZtC <- rowsum(C, s, reorder = TRUE)
    Cbar <- ZtC / w
    BtWC <- as.matrix(Matrix::crossprod(B, ZtC))
    CtC <- crossprod(C)
    CbWCb <- crossprod(Cbar, ZtC)
    g <- length(lambdas)
    facs <- H <- Gi <- vector("list", g)
    logdetS <- logdetG <- numeric(g)
    for (i in seq_len(g)) {
        S <- model$BtWB + lambdas[i] * model$Q
        if (is.null(fac0)) fac0 <- Matrix::Cholesky(S, perm = TRUE, LDL = FALSE)
        facs[[i]] <- Matrix::update(fac0, S)
        H[[i]] <- as.matrix(Matrix::solve(facs[[i]], BtWC))
        Gi[[i]] <- CtC - crossprod(BtWC, H[[i]])
        logdetS[i] <- .chol_logdet(facs[[i]])
        logdetG[i] <- as.numeric(determinant(Gi[[i]], logarithm = TRUE)$modulus)
    }
    Hc <- H[[center]]; Gind <- Gi[[center]]
    if (!is.null(YS)) {
        # Omega^-1 = I - Z B S^-1 B'Z' (Woodbury): the smooth field integrated out at lambda_center.
        omega_inv <- function(Y) Y - (as.matrix(B %*% as.matrix(Matrix::solve(facs[[center]],
            as.matrix(Matrix::crossprod(B, rowsum(Y, s, reorder = TRUE)))))))[s, , drop = FALSE]
        o1 <- drop(omega_inv(matrix(1, N, 1L))); den1 <- sum(o1)
        OY <- omega_inv(YS)
        G0 <- .tcp(YS, OY, threads, use_blas)
        Lg <- chol((G0 + t(G0)) / 2); Li <- backsolve(Lg, diag(ncol(YS)))
        YS <- .tmm(YS, Li, threads, use_blas); OY <- .tmm(OY, Li, threads, use_blas)
        Y1 <- OY - outer(o1, colSums(OY) / den1)
        Racc <- matrix(0, N, ncol(YS)); Tacc <- matrix(0, ncol(YS), ncol(YS))
    }
    prss <- numeric(g); total <- 0; ktime <- 0
    Uacc <- if (k && update_U) matrix(0, N, k) else NULL
    AtA <- matrix(0, k, k)
    acc <- function(Y) if (is.null(Y)) NULL else list(Y = matrix(0, nrow(Y), ncol(Y)), T = matrix(0, ncol(Y), ncol(Y)))
    aU <- acc(YU); aV <- acc(YV); aP <- acc(YP)
    if (!is.null(YP)) Pcoef <- solve(CbWCb, t(ZtC))
    sw <- sqrt(w)
    if (!is.null(bedscan)) {
        Yt <- t(if (is.null(YS)) C else cbind(C, Y1))
        RtT <- if (is.null(YS)) matrix(0, 0L, N) else matrix(0, ncol(YS), N)
    }
    for (st in seq(1L, p, by = block)) {
        en <- min(p, st + block - 1L)
        if (!is.null(bedscan)) {
            tk <- proc.time()[[3L]]
            ag <- bedscan(st, en, Yt, q, RtT)
            ktime <- ktime + proc.time()[[3L]] - tk
            XtW <- ag$XtY
        } else {
            X <- reader(st, en)
            if (nrow(X) != N || ncol(X) != en - st + 1L) stop("Reader returned a block of the wrong size")
            ag <- cpp_joint_aggregate(X, s, n, threads)
            # One threaded pass over the block for C'X and X'(P1' Omega^-1 Y).
            XtW <- if (use_blas) crossprod(X, if (is.null(YS)) C else cbind(C, Y1)) else
                cpp_joint_crossprod(X, if (is.null(YS)) C else cbind(C, Y1), threads)
        }
        zb <- as.matrix(Matrix::crossprod(B, ag$sums))
        CtX <- t(XtW[, seq_len(q), drop = FALSE])
        for (i in seq_len(g)) {
            t_i <- as.matrix(Matrix::solve(facs[[i]], zb))
            e <- CtX - crossprod(H[[i]], zb)
            c_i <- solve(Gi[[i]], e)
            # Individual-level penalized RSS: x'x - z'M^-1 z by block elimination.
            prss[i] <- prss[i] + ag$ss - sum(zb * t_i) - sum(e * c_i)
            if (i == center) { tc <- t_i; cc <- c_i }
        }
        theta <- tc - Hc %*% cc
        F <- as.matrix(B %*% theta)
        if (!is.null(YS)) {
            Tb <- XtW[, -seq_len(q), drop = FALSE]
            if (is.null(bedscan)) Racc <- Racc + if (use_blas) X %*% Tb else cpp_joint_prod(X, Tb, threads)
            Tacc <- Tacc + crossprod(Tb)
        }
        if (k) {
            A <- t(cc[-1L, , drop = FALSE])
            AtA <- AtA + crossprod(A)
            if (update_U) {
                FA <- F %*% A
                Uacc <- Uacc + .spde_mm(X, A) -
                    matrix(drop(crossprod(cc[1L, ], A)), N, k, byrow = TRUE) - FA[s, , drop = FALSE]
            }
        }
        if (!is.null(aU)) {
            M <- crossprod(X, YU) - crossprod(cc, crossprod(C, YU)) -
                crossprod(F, rowsum(YU, s, reorder = TRUE))
            aU$Y <- aU$Y + .spde_mm(X, M) - C %*% (cc %*% M) - (F %*% M)[s, , drop = FALSE]
            aU$T <- aU$T + crossprod(M)
        }
        if (!is.null(aV)) {
            Ft <- F * sw
            M <- crossprod(Ft, YV)
            aV$Y <- aV$Y + Ft %*% M; aV$T <- aV$T + crossprod(M)
        }
        if (!is.null(aP)) {
            Ft <- (F - Cbar %*% (Pcoef %*% F)) * sw
            M <- crossprod(Ft, YP)
            aP$Y <- aP$Y + Ft %*% M; aP$T <- aP$T + crossprod(M)
        }
        total <- total + en - st + 1L
    }
    if (!is.null(bedscan) && !is.null(YS)) Racc <- t(RtT)
    # One objective for U and lambda: individual-level REML with C = [1, U] as the fixed design,
    # sigma^2 profiled, divided by p. -log|C'C| makes it invariant to U -> U T + 1 v'.
    df <- N - q
    logdetCtC <- as.numeric(determinant(CtC, logarithm = TRUE)$modulus)
    reml <- df * log(prss / (total * df)) + logdetS + logdetG - m * log(lambdas) - model$logdetQ - logdetCtC
    objective <- prss[center]
    Uprop <- NULL
    if (k && update_U) {
        # Stationarity of the REML in U: U A'A = R0 A - (L/df) [(C - Z B H) G^-1 - C (C'C)^-1]_U.
        Pc <- C - (as.matrix(B %*% Hc))[s, , drop = FALSE]
        corr <- (objective / df) * (Pc %*% solve(Gind) - C %*% solve(CtC))[, -1L, drop = FALSE]
        Uprop <- (Uacc - corr) %*% solve(AtA)
    }
    sub <- NULL
    if (!is.null(YS)) {
        # Subspace iteration for X_c X_c' Omega^-1 U = U D (U'Omega^-1 U = I), X_c = GLS-centred X.
        e <- eigen((Tacc + t(Tacc)) / 2, symmetric = TRUE)
        sub <- list(ritz = .tmm(YS, e$vectors, threads, use_blas), values = e$values,
                    next_basis = Racc - outer(rep(1, N), colSums(o1 * Racc) / den1),
                    basis = YS, omega_inv = omega_inv)
    }
    list(reml = reml, prss = prss, logdetS = logdetS, logdetG = logdetG, objective = objective, Uacc = Uacc, Uprop = Uprop, AtA = AtA, sub = sub,
         aU = aU, aV = aV, aP = aP, p = total, fac0 = fac0, kernel_seconds = ktime,
         concurvity = if (k) kappa(Gind) else 1)
}

# Tall N x small products. R's reference BLAS (the Windows default) is single-threaded, so these go to
# the threaded C++ kernels unless an optimized BLAS is linked.
.tcp <- function(A, B, threads, blas) if (blas) crossprod(A, B) else cpp_joint_crossprod(A, B, threads)
.tmm <- function(A, T, threads, blas) if (blas) A %*% T else cpp_joint_prod(A, T, threads)

# 1 - smallest cosine of the principal angles between the centred column spans of U1 and U2, from
# small Gram matrices only (no tall QR).
.span_change <- function(U1, U2, threads, blas) {
    U1 <- sweep(U1, 2L, colMeans(U1)); U2 <- sweep(U2, 2L, colMeans(U2))
    r <- tryCatch({
        R1 <- chol(.tcp(U1, U1, threads, blas)); R2 <- chol(.tcp(U2, U2, threads, blas))
        Mx <- backsolve(R1, .tcp(U1, U2, threads, blas), transpose = TRUE) %*% backsolve(R2, diag(ncol(U2)))
        1 - min(svd(Mx, nu = 0, nv = 0)$d)
    }, error = function(e) NA_real_)
    if (is.na(r)) r <- 1 - min(svd(crossprod(qr.Q(qr(U1)), qr.Q(qr(U2))), nu = 0, nv = 0)$d)
    max(0, r)
}

.orth_center <- function(U) {
    U <- sweep(U, 2L, colMeans(U))
    Qr <- qr(U)
    Qm <- qr.Q(Qr)[, seq_len(Qr$rank), drop = FALSE]
    Qm * sqrt(nrow(U))
}

.ritz <- function(Y, a, r) {
    e <- eigen((a$T + t(a$T)) / 2, symmetric = TRUE)
    keep <- seq_len(min(r, sum(e$values > 0)))
    list(vectors = Y %*% e$vectors[, keep, drop = FALSE], values = e$values[keep])
}

svgpc_joint_fit <- function(model, k = 10L, r = 10L, frequencies = NULL, reader = NULL, snps = NULL,
                            init = c("spatial", "table"), U_init = NULL, lambda = NULL,
                            max_iter = 30L, tol = 1e-4, v_scans = 1L, krylov = 4L, warm = 1L, cache_gb = 8,
                            block_size = NULL, threads = 0L, kernel = c("auto", "fused", "cpp", "blas"),
                            seed = 1L, verbose = TRUE) {
    if (!inherits(model, "svgpc_joint_model")) stop("Use svgpc_joint_prepare() first")
    start_time <- proc.time()[[3L]]
    init <- match.arg(init)
    kernel <- match.arg(kernel)
    # Tall-skinny products: optimized BLAS when R is linked to one, else threaded C++ kernels
    # (R's reference BLAS, the Windows default, is single-threaded and unblocked).
    # A PLINK BED input defaults to the fused 2-bit kernel, which never forms double genotype blocks.
    if (kernel == "auto" && is.null(reader)) kernel <- "fused"
    if (kernel == "fused" && !is.null(reader)) stop("kernel = 'fused' reads the BED file directly; drop reader")
    use_blas <- if (kernel == "auto") grepl("openblas|mkl|accelerate|blis|flexiblas|atlas",
        tolower(paste(extSoftVersion()[["BLAS"]], La_library()))) else kernel == "blas"
    N <- model$N; s <- model$location; w <- model$counts
    bedscan <- NULL
    if (is.null(reader)) {
        bed <- .joint_bed(model$spatial, frequencies)
        # snps: NULL (all variants), a count (the first snps variants), or a vector of 1-based variant
        # indices (e.g. an LD-pruned set).
        idx <- if (is.null(snps)) seq_len(bed$snps) else if (length(snps) == 1L)
            seq_len(min(bed$snps, .positive_integer(snps, "snps"))) else sort(unique(as.integer(snps)))
        if (anyNA(idx) || min(idx) < 1L || max(idx) > bed$snps) stop("snps indices out of range")
        snps <- length(idx)
        reader <- function(st, en) cpp_joint_bed_block(bed$path, bed$file_samples, bed$rows,
                                                       st - 1L, en - st + 1L, bed$af, as.integer(threads))
        if (kernel != "fused" && length(snps) && !identical(idx, seq_len(snps)))
            stop("A variant subset needs kernel = 'fused'")
        # The selected variants' 2-bit rows are read into memory once when they fit in cache_gb, so
        # scans after the first do no disk I/O (a busy or slow disk otherwise dominates a scan).
        gcache <- raw(0)
        if (kernel == "fused" && ceiling(bed$file_samples / 4) * snps / 2^30 <= cache_gb) {
            t_read <- proc.time()[[3L]]
            gcache <- cpp_joint_bed_read(bed$path, bed$file_samples, idx - 1L)
            if (verbose) message(sprintf("Cached %d variants (%.2f GB) in %.1fs", snps,
                                         length(gcache) / 2^30, proc.time()[[3L]] - t_read))
        }
        if (kernel == "fused") bedscan <- function(st, en, Yt, q, RtT) {
            ii <- idx[st:en]
            cpp_joint_bed_scan(bed$path, bed$file_samples, bed$rows, ii - 1L, bed$af[ii],
                               s, model$n, Yt, q, RtT, as.integer(threads), gcache, st - 1)
        }
    } else {
        if (!is.function(reader)) stop("reader must be function(start, end) returning N x block standardized genotypes")
        snps <- .positive_integer(snps, "snps")
    }
    p <- snps
    if (is.null(block_size)) block_size <- if (kernel == "fused") 256L else max(16L, min(1024L, as.integer(2e7 / N)))
    block_size <- .positive_integer(block_size, "block_size")
    k <- as.integer(k); r <- as.integer(r)
    if (k < 0L || r < 1L) stop("k must be nonnegative and r positive")
    set.seed(seed)
    # Warm start: fit first on every warm-th variant (1/warm of the cost per scan), then start the
    # full fit from that U and lambda, so the large early moves of lambda happen on cheap scans.
    warm <- as.integer(warm); wf <- NULL
    if (!is.null(bedscan) && warm > 1L && is.null(U_init) && is.null(lambda) && k > 0L &&
        p %/% warm >= 5L * (k + 10L)) {
        if (verbose) message("Warm start on ", length(seq(1L, p, by = warm)), " of ", p, " variants")
        wf <- svgpc_joint_fit(model, k = k, r = r, frequencies = frequencies, snps = idx[seq(1L, p, by = warm)],
                              init = init, max_iter = max_iter, tol = tol, v_scans = 0L, krylov = krylov,
                              warm = 1L, block_size = block_size, threads = threads, kernel = "fused",
                              seed = seed, verbose = verbose)
        U_init <- wf$U; lambda <- wf$lambda
    }
    ell <- function(x) min(x + 10L, model$n - 1L)
    # U given lambda: generalized eigenproblem X_c X_c' Omega^-1 U = U D, solved by block subspace
    # iteration (one scan per step, oversampled by 10 columns). Lambda: exact grid minimization of V.
    ellU <- min(k + 10L, N - 1L)
    U <- matrix(0, N, 0L); YS <- NULL
    start <- NULL
    if (k) {
        start <- if (!is.null(U_init)) .matrix_double(U_init) else if (init == "table") {
            X0 <- as.matrix(model$spatial$table[, -c(1:4), drop = FALSE])
            if (ncol(X0) < k) stop("The table has fewer than k fixed-effect columns for init = 'table'")
            X0[, seq_len(k), drop = FALSE]
        } else NULL
        YS <- cbind(start, matrix(stats::rnorm(N * (ellU - if (is.null(start)) 0L else ncol(start))), N))[, seq_len(ellU), drop = FALSE]
        YS <- qr.Q(qr(sweep(YS, 2L, colMeans(YS))))
        if (!is.null(start)) U <- .orth_center(start)
    }
    eta <- if (is.null(lambda)) log(stats::median(Matrix::diag(model$BtWB)) /
                                    stats::median(Matrix::diag(model$Q))) else log(lambda)
    delta <- if (is.null(lambda)) 8 else 1
    YV <- qr.Q(qr(matrix(stats::rnorm(model$n * ell(r)), model$n, ell(r))))
    # Block Krylov acceleration: the last `krylov` iterates Y_j and their images M Y_j are kept, and U
    # is the Rayleigh-Ritz solution over their span (M = X_c X_c' Omega^-1 is symmetric in the
    # Omega^-1 inner product, so V'Omega^-1 M V needs no extra scan). The history restarts when lambda
    # moves, because M depends on lambda.
    krylov <- max(1L, as.integer(krylov)); histV <- histR <- list(); eta_hist <- NA_real_
    YP <- YV
    trace <- list(); fac0 <- NULL; converged <- FALSE; iter <- 0L; V_last <- NA_real_; Ueig <- numeric()
    repeat {
        iter <- iter + 1L
        t0 <- proc.time()[[3L]]
        offsets <- if (delta > 4) seq(-delta, delta, length.out = 9L) else delta * c(-1, -.5, 0, .5, 1)
        lams <- exp(eta + offsets); ic <- which(offsets == 0)
        res <- .joint_scan(model, reader, p, block_size, U, lams, ic,
                           YV = YV, YP = if (k) YP, YS = YS, fac0 = fac0, threads = as.integer(threads),
                           use_blas = use_blas, bedscan = bedscan)
        fac0 <- res$fac0
        V_now <- res$reml[ic]
        j <- which.min(res$reml); eta_old <- eta
        if (j > 1L && j < length(lams)) {
            h <- offsets[j + 1L] - offsets[j]
            f <- res$reml[(j - 1L):(j + 1L)]
            den <- f[3] - 2 * f[2] + f[1]
            st <- if (den > 0) -h * (f[3] - f[1]) / (2 * den) else 0
            eta <- eta + offsets[j] + max(-h, min(h, st))
            delta <- max(0.05, min(1, abs(eta - eta_old) * 2, delta / 2))
            fmin <- if (den > 0) f[2] - (f[3] - f[1])^2 / (8 * den) else f[2]
        } else {
            eta <- eta + offsets[j]
            delta <- max(1, delta)
            fmin <- res$reml[j]
        }
        # Predicted REML gain of the lambda step. Lambda is held when the gain is below tol: the REML is
        # flat in lambda near its optimum, moving it there only restarts the Krylov history (M depends
        # on lambda) and makes lambda wander. The held step is applied at the end.
        gain <- res$reml[ic] - min(fmin, min(res$reml))
        eta_prop <- eta
        atol <- tol * max(1, abs(res$reml[ic]))  # tol is relative to |REML|
        if (gain < atol) eta <- eta_old
        YV <- qr.Q(qr(res$aV$Y))
        if (k) YP <- qr.Q(qr(res$aP$Y))
        change <- NA_real_
        if (k) {
            if (is.na(eta_hist) || abs(log(lams[ic]) - eta_hist) > 1e-12) histV <- histR <- list()
            eta_hist <- log(lams[ic])
            histV <- c(histV, list(res$sub$basis)); histR <- c(histR, list(res$sub$next_basis))
            if (length(histV) > krylov) { histV <- histV[-1L]; histR <- histR[-1L] }
            Vall <- do.call(cbind, histV); Rall <- do.call(cbind, histR)
            OV <- res$sub$omega_inv(Vall)
            G1 <- .tcp(OV, Vall, threads, use_blas)
            em <- eigen((G1 + t(G1)) / 2, symmetric = TRUE)
            keep <- em$values > 1e-10 * em$values[1L]
            Wm <- em$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(em$values[keep]), sum(keep))
            Gm <- crossprod(Wm, .tcp(OV, Rall, threads, use_blas) %*% Wm)
            eg <- eigen((Gm + t(Gm)) / 2, symmetric = TRUE)
            cf <- Wm %*% eg$vectors
            Unew <- .tmm(Vall, cf[, seq_len(k), drop = FALSE], threads, use_blas)
            Ueig <- eg$values[seq_len(k)] / res$p
            if (ncol(U)) change <- .span_change(U, Unew, threads, use_blas)
            U <- Unew
            YS <- .tmm(Rall, cf[, seq_len(ncol(YS)), drop = FALSE], threads, use_blas)
        }
        dV <- if (is.na(V_last)) NA_real_ else abs(V_last - V_now)
        V_last <- V_now
        trace[[iter]] <- data.frame(iter = iter, lambda = lams[ic], reml = V_now,
                                    objective = res$objective, U_change = change,
                                    seconds = proc.time()[[3L]] - t0, kernel_seconds = res$kernel_seconds)
        if (verbose) message(sprintf("scan %d: lambda=%.4g REML=%.8g dREML=%.2e Uchange=%.2e (%.1fs, kernel %.1fs)",
                                     iter, lams[ic], V_now, dV, change, proc.time()[[3L]] - t0, res$kernel_seconds))
        if (!is.na(dV) && dV < atol && gain < atol && (is.na(change) || change < 1e-3)) { converged <- TRUE; break }
        if (iter >= max_iter) break
    }
    # Final scans at fixed U and lambda: V and its U-orthogonal version by subspace iteration.
    # Final lambda = the one U was solved at; the held step would change the REML by less than tol.
    lam <- exp(eta)
    # v_scans = 0 (used by the warm start): return U and lambda only.
    if (as.integer(v_scans) < 1L) return(list(U = U, lambda = lam, trace = do.call(rbind, trace)))
    ZtC <- rowsum(cbind(1, U), s, reorder = TRUE)
    # Keep the V_orthogonal basis inside the W-complement of [1, U-bar] at the final U.
    Kc <- qr.Q(qr(ZtC / sqrt(w)))
    YP <- qr.Q(qr(YP - Kc %*% crossprod(Kc, YP)))
    v_scans <- as.integer(v_scans)
    for (v in seq_len(v_scans)) {
        t0 <- proc.time()[[3L]]
        res <- .joint_scan(model, reader, p, block_size, U, lam, 1L,
                           update_U = FALSE,
                           YV = YV, YP = YP, fac0 = fac0, threads = as.integer(threads),
                           use_blas = use_blas, bedscan = bedscan)
        if (verbose) message(sprintf("V scan %d: REML=%.6g (%.1fs)", v, res$reml, proc.time()[[3L]] - t0))
        if (v < v_scans) { YV <- qr.Q(qr(res$aV$Y)); YP <- qr.Q(qr(res$aP$Y)) }
    }
    rv <- .ritz(YV, res$aV, r); rp <- .ritz(YP, res$aP, r)
    V <- rv$vectors / sqrt(w); Vp <- rp$vectors / sqrt(w)
    colnames(V) <- paste0("SPC", seq_len(ncol(V))); colnames(Vp) <- paste0("SPC", seq_len(ncol(Vp)))
    if (ncol(U)) {
        U <- sweep(U, 2L, colMeans(U)); U <- sweep(U, 2L, sqrt(colMeans(U^2)), "/")
        colnames(U) <- paste0("U", seq_len(ncol(U)))
    }
    ids <- model$spatial$table[, c("FID", "IID")]
    out <- list(U = U, U_eigenvalues = Ueig,
                V = V, V_eigenvalues = rv$values / res$p,
                V_orthogonal = Vp, V_orthogonal_eigenvalues = rp$values / res$p,
                individual = cbind(ids, U, V[s, , drop = FALSE]),
                individual_orthogonal = cbind(ids, U, Vp[s, , drop = FALSE]),
                locations = model$spatial$locations, lambda = lam, reml = res$reml,
                objective = res$objective, kappa = model$kappa, k = ncol(U), p = res$p,
                concurvity = res$concurvity, converged = converged,
                trace = do.call(rbind, trace), warm_trace = wf$trace, block_size = block_size,
                kernel = if (kernel == "fused") "fused" else if (use_blas) "blas" else "cpp",
                seconds = proc.time()[[3L]] - start_time)
    class(out) <- "svgpc_joint_fit"
    out
}

print.svgpc_joint_fit <- function(x, ...) {
    cat("Joint U + SPDE fit:", x$k, "U factors,", ncol(x$V), "spatial PCs from", x$p, "SNPs\n")
    cat("lambda:", format(x$lambda), "; REML:", format(x$reml), "; kappa:", x$kappa,
        "; converged:", x$converged, "after", nrow(x$trace), "scans;", round(x$seconds, 1), "s\n")
    invisible(x)
}
