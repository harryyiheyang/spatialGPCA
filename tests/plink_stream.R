library(spatialGPCA)
folder <- system.file("extdata", package = "spatialGPCA")
D <- as.matrix(read.table(file.path(folder, "coordinates.txt")))
S <- read.table(file.path(folder, "samples.txt"), header = TRUE, colClasses = "character")
ff <- read.table(file.path(folder, "frequencies.afreq"), header = TRUE, comment.char = "", check.names = FALSE)
i <- as.integer(S$location)
table <- data.frame(FID = S$FID, IID = S$IID, Lat = 35 + D[i, 1], Lon = -100 + D[i, 2],
                    GPC1 = D[i, 3], GPC2 = D[i, 4], GPC3 = D[i, 5], income = sin(i))
set.seed(606)
table <- table[-nrow(table), ]
table <- rbind(table, transform(table[1, ], FID = "outside", IID = "outside"))
table <- table[sample(nrow(table)), ]
relative <- function(a, b) sqrt(sum((a - b)^2)) / max(1, sqrt(sum(b^2)))
for (name in c("hard", "missing")) {
    cls <- svgpc_cluster(table, rho = 20, clusters = 24, bed = file.path(folder, name))
    ll <- spatialGPCA:::.geo_inverse(as.matrix(cls$locations[, c("x", "y")]), cls$projection$origin)
    stopifnot(max(abs(as.matrix(ll) - as.matrix(cls$locations[, c("Lat", "Lon")]))) < 1e-10)
    stopifnot(identical(cls$table$IID, S$IID[-nrow(S)]), sum(cls$locations$count) == nrow(S) - 1L)
    raw <- readBin(file.path(folder, paste0(name, ".bed")), "raw", n = file.info(file.path(folder, paste0(name, ".bed")))$size)
    ns <- nrow(S); p <- nrow(ff); stride <- ceiling(ns / 4)
    G <- matrix(NA_real_, ns, p)
    for (j in seq_len(p)) {
        for (i in seq_len(ns)) {
            byte <- as.integer(raw[3 + (j - 1) * stride + (i - 1) %/% 4 + 1])
            code <- bitwAnd(bitwShiftR(byte, 2L * ((i - 1L) %% 4L)), 3L)
            G[i, j] <- c(2, NA, 1, 0)[code + 1L]
        }
        G[is.na(G[, j]), j] <- 2 * ff$ALT_FREQS[j]
    }
    G <- G[cls$file_rows, , drop = FALSE]
    Y <- rowsum(G, cls$sample_location, reorder = TRUE) / cls$locations$count
    Y <- sweep(sweep(Y, 2, 2 * ff$ALT_FREQS, "-"), 2, sqrt(2 * ff$ALT_FREQS * (1 - ff$ALT_FREQS)), "/")
    ref <- svgpc_prepare(cls)
    svgpc_accumulate_matrix(ref, Y, 17L)
    for (flip in c(FALSE, TRUE)) {
        model <- svgpc_prepare(cls)
        frequencies <- if (flip) data.frame(ID = ff$ID, A1 = ff$REF, A2 = ff$ALT, AF = 1 - ff$ALT_FREQS) else file.path(folder, "frequencies.afreq")
        svgpc_accumulate(model, frequencies, block_size = 23L)
        err <- relative(svgpc_statistics(model)$covariance_sum, svgpc_statistics(ref)$covariance_sum)
        cat(name, "intersect BED; frequency flip", flip, "covariance error", err, "\n")
        stopifnot(err < 1e-11, svgpc_info(model)$missing_genotypes == if (name == "missing") 2 else 0)
        s <- spatialGPCA:::cpp_snapshot(model)
        X <- as.matrix(cls$table[, -c(1:4), drop = FALSE])
        C <- rowsum(X, cls$sample_location, reorder = TRUE) / cls$locations$count
        stopifnot(max(abs(s$C - cbind(1, C))) < 1e-12)
        for (method in c("REML", "GCV")) {
            par <- svgpc_select_lambda(model, method)
            fit <- svgpc_fit(model, parameters = par, components = 3)
            stopifnot(fit$selector == method, identical(fit$lambda, par$lambda), svgpc_info(model)$decoded_variants == p)
        }
    }
}
cache <- tempfile(fileext = ".rds")
svgpc_save(model, cache)
loaded <- svgpc_load(cache)
stopifnot(identical(attr(loaded, "spatial")$table, cls$table))
stopifnot(relative(svgpc_statistics(loaded)$covariance_sum, svgpc_statistics(model)$covariance_sum) == 0)
unlink(cache)
pdf(tempfile(fileext = ".pdf"), width = 13, height = 5)
p1 <- plot(cls)
cls$rho <- 40
p2 <- plot(cls)
stopifnot(identical(p1$data, p2$data))
stopifnot(abs(exp(-attr(p2, "range_circle")$radius / 40) - .1) < 1e-12)
fit$rho <- 60
p3 <- plot(fit)
stopifnot(abs(exp(-attr(p3, "range_circle")$radius / 60) - .1) < 1e-12)
dev.off()
stopifnot(!("lambda" %in% names(cls)), !("lambda" %in% names(formals(svgpc_fit))),
          !("conditioning" %in% names(formals(svgpc_prepare))), !("intercept" %in% names(formals(svgpc_prepare))))
cat("All BED intersection, fixed-effect, imputation, selection, cache and editable-rho tests passed.\n")
