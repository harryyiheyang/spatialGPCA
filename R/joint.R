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
                        update_U = TRUE, YU = NULL, YV = NULL, YP = NULL, fac0 = NULL, threads = 0L) {
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
    prss <- numeric(g); total <- 0
    Uacc <- if (k && update_U) matrix(0, N, k) else NULL
    AtA <- matrix(0, k, k)
    acc <- function(Y) if (is.null(Y)) NULL else list(Y = matrix(0, nrow(Y), ncol(Y)), T = matrix(0, ncol(Y), ncol(Y)))
    aU <- acc(YU); aV <- acc(YV); aP <- acc(YP)
    if (!is.null(YP)) Pcoef <- solve(CbWCb, t(ZtC))
    sw <- sqrt(w)
    for (st in seq(1L, p, by = block)) {
        en <- min(p, st + block - 1L)
        X <- reader(st, en)
        if (nrow(X) != N || ncol(X) != en - st + 1L) stop("Reader returned a block of the wrong size")
        ag <- cpp_joint_aggregate(X, s, n, threads)
        zb <- as.matrix(Matrix::crossprod(B, ag$sums))
        CtX <- .spde_mm(C, X, transA = TRUE)
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
        total <- total + ncol(X)
    }
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
    list(reml = reml, prss = prss, objective = objective, Uacc = Uacc, Uprop = Uprop, AtA = AtA,
         aU = aU, aV = aV, aP = aP, p = total, fac0 = fac0,
         concurvity = if (k) kappa(Gind) else 1)
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
                            init = c("spatial", "table"), U = NULL, lambda = NULL,
                            max_iter = 30L, tol = 1e-4, init_scans = 2L, v_scans = 3L,
                            block_size = NULL, threads = 0L, seed = 1L, verbose = TRUE) {
    if (!inherits(model, "svgpc_joint_model")) stop("Use svgpc_joint_prepare() first")
    start_time <- proc.time()[[3L]]
    init <- match.arg(init)
    N <- model$N; s <- model$location; w <- model$counts
    if (is.null(reader)) {
        bed <- .joint_bed(model$spatial, frequencies)
        snps <- bed$snps
        reader <- function(st, en) cpp_joint_bed_block(bed$path, bed$file_samples, bed$rows,
                                                       st - 1L, en - st + 1L, bed$af, as.integer(threads))
    } else {
        if (!is.function(reader)) stop("reader must be function(start, end) returning N x block standardized genotypes")
        snps <- .positive_integer(snps, "snps")
    }
    p <- snps
    if (is.null(block_size)) block_size <- max(16L, min(1024L, as.integer(2e7 / N)))
    block_size <- .positive_integer(block_size, "block_size")
    k <- as.integer(k); r <- as.integer(r)
    if (k < 0L || r < 1L) stop("k must be nonnegative and r positive")
    set.seed(seed)
    ell <- function(x) min(x + 10L, model$n - 1L)
    # Initial U: supplied, table columns (two-stage start), or spatial-first power iteration.
    phase_init <- FALSE
    if (!is.null(U)) {
        U <- .orth_center(.matrix_double(U))
    } else if (k && init == "table") {
        X0 <- as.matrix(model$spatial$table[, -c(1:4), drop = FALSE])
        if (ncol(X0) < k) stop("The table has fewer than k fixed-effect columns for init = 'table'")
        U <- .orth_center(X0[, seq_len(k), drop = FALSE])
    } else if (k) {
        phase_init <- TRUE
        YU <- qr.Q(qr(matrix(stats::rnorm(N * (k + 10L)), N, k + 10L)))
        U <- matrix(0, N, 0L)
    } else U <- matrix(0, N, 0L)
    eta <- if (is.null(lambda)) log(stats::median(Matrix::diag(model$BtWB)) /
                                    stats::median(Matrix::diag(model$Q))) else log(lambda)
    delta <- if (is.null(lambda)) 8 else 1
    YV <- qr.Q(qr(matrix(stats::rnorm(model$n * ell(r)), model$n, ell(r))))
    trace <- list(); fac0 <- NULL; converged <- FALSE
    iter <- 0L; it_init <- 0L
    U_acc <- NULL; Uprop <- NULL; V_acc <- Inf; step <- 1; omega <- 1
    repeat {
        iter <- iter + 1L
        t0 <- proc.time()[[3L]]
        offsets <- if (delta > 4) seq(-delta, delta, length.out = 9L) else delta * c(-1, -.5, 0, .5, 1)
        lams <- exp(eta + offsets); ic <- which(offsets == 0)
        res <- .joint_scan(model, reader, p, block_size, U, lams, ic,
                           update_U = !phase_init && ncol(U) > 0L,
                           YU = if (phase_init) YU else NULL, YV = YV, fac0 = fac0,
                           threads = as.integer(threads))
        fac0 <- res$fac0
        V_now <- res$reml[ic]
        status <- if (phase_init) "init" else "accepted"
        if (!phase_init && ncol(U) && V_now > V_acc + 1e-10 * abs(V_acc)) {
            # The U step increased the REML: backtrack towards the last accepted U.
            status <- "rejected"
            step <- if (step > 1) 1 else step / 2
            omega <- 1
            if (step < 1 / 64) { U <- U_acc; converged <- TRUE; status <- "stalled"; }
            else U <- .orth_center(U_acc + step * (Uprop - U_acc))
        }
        change <- NA_real_
        if (status != "rejected" && status != "stalled") {
            # lambda: exact minimization of the same REML on the grid, then a parabola step.
            j <- which.min(res$reml); eta_old <- eta
            if (j > 1L && j < length(lams)) {
                h <- offsets[j + 1L] - offsets[j]
                f <- res$reml[(j - 1L):(j + 1L)]
                den <- f[3] - 2 * f[2] + f[1]
                st <- if (den > 0) -h * (f[3] - f[1]) / (2 * den) else 0
                eta <- eta + offsets[j] + max(-h, min(h, st))
                delta <- max(0.05, min(1, abs(eta - eta_old) * 2, delta / 2))
            } else {
                eta <- eta + offsets[j]
                delta <- max(1, delta)
            }
            YV <- qr.Q(qr(res$aV$Y))
            if (phase_init) {
                it_init <- it_init + 1L
                if (it_init >= init_scans && (abs(eta - eta_old) < 0.05 || it_init >= init_scans + 4L)) {
                    U <- .orth_center(.ritz(YU, res$aU, k)$vectors)
                    phase_init <- FALSE
                } else YU <- qr.Q(qr(res$aU$Y))
            } else if (ncol(U)) {
                dV <- V_acc - min(res$reml)
                U_acc <- U; V_acc <- min(V_now, min(res$reml))
                # Over-relaxed step along the stationarity update; the next scan's REML
                # check rejects it and backtracks if it fails to decrease the objective.
                omega <- if (step >= 1) min(omega * 1.5, 8) else 1
                step <- omega
                Uprop <- res$Uprop
                Unew <- .orth_center(U + step * (Uprop - U))
                change <- 1 - min(svd(crossprod(U, Unew) / N, nu = 0, nv = 0)$d)
                U <- Unew
                if (is.finite(dV) && dV < tol && abs(eta - eta_old) < 0.01) converged <- TRUE
            } else {
                dV <- V_acc - min(res$reml); V_acc <- min(res$reml)
                if (is.finite(dV) && dV < tol && abs(eta - eta_old) < 0.01) converged <- TRUE
            }
        }
        trace[[iter]] <- data.frame(iter = iter, step = status, lambda = lams[ic], reml = V_now,
                                    objective = res$objective, U_change = change,
                                    seconds = proc.time()[[3L]] - t0)
        if (verbose) message(sprintf("scan %d [%s]: lambda=%.4g REML=%.8g Uchange=%.2e (%.1fs)",
                                     iter, status, lams[ic], V_now, change, proc.time()[[3L]] - t0))
        if (converged) { if (!is.null(U_acc)) U <- U_acc; break }
        if (iter >= max_iter) { if (!is.null(U_acc)) U <- U_acc; break }
    }
    # Final scans at fixed U and lambda: V and its U-orthogonal version by subspace iteration.
    lam <- exp(eta)
    ZtC <- rowsum(cbind(1, U), s, reorder = TRUE)
    YP <- YV
    v_scans <- max(1L, as.integer(v_scans))
    for (v in seq_len(v_scans)) {
        t0 <- proc.time()[[3L]]
        res <- .joint_scan(model, reader, p, block_size, U, lam, 1L, update_U = FALSE,
                           YV = YV, YP = YP, fac0 = fac0, threads = as.integer(threads))
        if (verbose) message(sprintf("V scan %d: REML=%.6g (%.1fs)", v, res$reml, proc.time()[[3L]] - t0))
        if (v < v_scans) { YV <- qr.Q(qr(res$aV$Y)); YP <- qr.Q(qr(res$aP$Y)) }
    }
    rv <- .ritz(YV, res$aV, r); rp <- .ritz(YP, res$aP, r)
    V <- rv$vectors / sqrt(w); Vp <- rp$vectors / sqrt(w)
    colnames(V) <- paste0("SPC", seq_len(ncol(V))); colnames(Vp) <- paste0("SPC", seq_len(ncol(Vp)))
    Ueig <- numeric()
    if (ncol(U)) {
        e <- eigen(res$AtA / res$p, symmetric = TRUE)
        U <- U %*% e$vectors; Ueig <- e$values
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
                trace = do.call(rbind, trace), block_size = block_size,
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
