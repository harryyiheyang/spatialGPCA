library(spatialGPCA)

relative <- function(A, B) sqrt(sum((A - B)^2)) / max(1, sqrt(sum(B^2)))
folder <- system.file("extdata", package = "spatialGPCA")
D <- as.matrix(read.table(file.path(folder, "coordinates.txt")))
S <- read.table(file.path(folder, "samples.txt"), header = TRUE, colClasses = "character")
ff <- read.table(file.path(folder, "frequencies.afreq"), header = TRUE,
                 comment.char = "", check.names = FALSE)
i <- as.integer(S$location)
table <- data.frame(FID = S$FID, IID = S$IID, Lat = 35 + D[i, 1], Lon = -100 + D[i, 2],
                    GPC1 = D[i, 3], GPC2 = D[i, 4], GPC3 = D[i, 5], income = sin(i))

# Stage 1 has no genotype file or fitted covariates.
cls <- svgpc_cluster(table[, 1:4], rho = 20, clusters = 24)
original <- cls
stopifnot(is.null(cls$bed))
old <- svgpc_prepare(cls)
stopifnot(identical(attr(old, "spatial"), cls),
          !("conditioning" %in% names(formals(svgpc_prepare))))

prefix <- tempfile("gp_stages_")
for (ext in c("bed", "bim", "fam"))
    stopifnot(file.copy(file.path(folder, paste0("hard.", ext)), paste0(prefix, ".", ext)))

# Add spatial signal only to the temporary BED copy so selection retains PCs.
bytes <- readBin(paste0(prefix, ".bed"), "raw", n = file.info(paste0(prefix, ".bed"))$size)
stride <- ceiling(nrow(S) / 4)
for (j in 1:40) {
    coordinate <- D[as.integer(S$location), 1L + (j %% 2L)]
    dose <- as.integer(coordinate > 0.3) + as.integer(coordinate > 0.7)
    code <- c(3L, 2L, 0L)[dose + 1L]
    length(code) <- stride * 4L
    code[is.na(code)] <- 0L
    bytes[3L + (j - 1L) * stride + seq_len(stride)] <-
        as.raw(colSums(matrix(code, 4L) * c(1L, 4L, 16L, 64L)))
}
writeBin(bytes, paste0(prefix, ".bed"))

# A fitting table can change the sample intersection, coordinates and covariates.
table <- table[-c(4L, 9L, nrow(table)), ]
table$Lat[1L] <- table$Lat[1L] + 0.015
table <- rbind(table, transform(table[1L, ], FID = "outside", IID = "outside"))
table <- table[rev(seq_len(nrow(table))), ]
model <- svgpc_prepare(cls, table = table, bed = prefix)
spatial <- attr(model, "spatial")
idx <- match(S$IID, table$IID)
keep <- which(!is.na(idx))
matched <- table[idx[keep], ]
key <- paste(sprintf("%.17g", matched$Lat), sprintf("%.17g", matched$Lon), sep = "\t")
loc <- match(key, unique(key))
counts <- tabulate(loc, length(unique(key)))
stopifnot(identical(cls, original), identical(spatial$centres, original$centres),
          identical(spatial$projection, original$projection),
          identical(spatial$clustering, original$clustering),
          identical(spatial$table$IID, S$IID[keep]), identical(spatial$table$FID, S$FID[keep]),
          identical(spatial$file_rows, keep), identical(spatial$sample_location, loc),
          identical(spatial$locations$count, counts), sum(counts) == nrow(matched))

old_key <- paste(sprintf("%.17g", original$locations$Lat), sprintf("%.17g", original$locations$Lon), sep = "\t")
at <- match(unique(key), old_key)
stopifnot(anyNA(at), identical(spatial$locations$cluster[!is.na(at)], original$locations$cluster[at[!is.na(at)]]))
for (i in which(is.na(at))) {
    d <- (original$centres$x - spatial$locations$x[i])^2 + (original$centres$y - spatial$locations$y[i])^2
    stopifnot(spatial$locations$cluster[i] == which.min(d))
}
ll <- spatialGPCA:::.geo_forward(spatial$locations$Lat, spatial$locations$Lon, original$projection$origin)
state <- spatialGPCA:::cpp_snapshot(model)
C <- rowsum(as.matrix(matched[, -(1:4)]), loc, reorder = TRUE) / counts
stopifnot(relative(state$knots, as.matrix(original$centres[, c("x", "y")])) == 0,
          relative(state$xy, ll) == 0, relative(state$weights, counts) == 0,
          relative(state$C, cbind(1, C)) < 1e-12)

# Each optional preparation input also works independently.
bed_only <- svgpc_prepare(cls, bed = prefix)
stopifnot(identical(attr(bed_only, "spatial")$centres, original$centres),
          identical(attr(bed_only, "spatial")$table$IID, S$IID))
table_only <- svgpc_prepare(cls, table = table)
stopifnot(is.null(attr(table_only, "spatial")$bed),
          nrow(attr(table_only, "spatial")$table) == nrow(table),
          identical(attr(table_only, "spatial")$centres, original$centres))
supplement <- svgpc_prepare(attr(bed_only, "spatial"), table = table)
stopifnot(identical(attr(supplement, "spatial")$table, spatial$table),
          identical(attr(supplement, "spatial")$sample_location, spatial$sample_location))

# Independently decode the small fixture and aggregate the changed location groups.
bytes <- readBin(paste0(prefix, ".bed"), "raw", n = file.info(paste0(prefix, ".bed"))$size)
ns <- nrow(S); p <- nrow(ff); stride <- ceiling(ns / 4)
G <- matrix(NA_real_, ns, p)
for (j in seq_len(p)) {
    for (i in seq_len(ns)) {
        byte <- as.integer(bytes[3L + (j - 1L) * stride + (i - 1L) %/% 4L + 1L])
        code <- bitwAnd(bitwShiftR(byte, 2L * ((i - 1L) %% 4L)), 3L)
        G[i, j] <- c(2, NA, 1, 0)[code + 1L]
    }
    G[is.na(G[, j]), j] <- 2 * ff$ALT_FREQS[j]
}
G <- rowsum(G[keep, , drop = FALSE], loc, reorder = TRUE) / counts
G <- sweep(sweep(G, 2L, 2 * ff$ALT_FREQS, "-"), 2L,
           sqrt(2 * ff$ALT_FREQS * (1 - ff$ALT_FREQS)), "/")
reference <- svgpc_prepare(spatial)
svgpc_accumulate_matrix(reference, G, block_size = 17L)
svgpc_accumulate(model, ff, block_size = 23L)
stopifnot(relative(svgpc_statistics(model)$covariance_sum,
                   svgpc_statistics(reference)$covariance_sum) < 1e-11,
          relative(svgpc_statistics(model)$total_energy,
                   svgpc_statistics(reference)$total_energy) < 1e-11)

# Stage 2 selects once; stage 3 only changes the requested PC count.
parameters <- svgpc_select_lambda(model, "GCV")
stopifnot(is.finite(parameters$lambda))
decoded <- svgpc_info(model)$decoded_variants
stats <- svgpc_statistics(model)
unlink(paste0(prefix, ".bed"))
stopifnot(!file.exists(paste0(prefix, ".bed")))
R2 <- svgpc_fit(model, parameters = parameters, components = 2L)
R3 <- svgpc_fit(model, parameters = parameters, components = 3L)
P2 <- R2$outputs$raw_F_PCA$pcs
P3 <- R3$outputs$raw_F_PCA$pcs
stopifnot(ncol(P2) == 2L, ncol(P3) == 3L,
          identical(R2$lambda, parameters$lambda), identical(R3$lambda, parameters$lambda),
          identical(R2$evaluations, parameters$evaluations), identical(R3$evaluations, parameters$evaluations),
          relative(R2$outputs$raw_F_PCA$eigenvalues, R3$outputs$raw_F_PCA$eigenvalues[1:2]) < 1e-10,
          relative(tcrossprod(P2), tcrossprod(P3[, 1:2, drop = FALSE])) < 1e-9,
          identical(svgpc_statistics(model), stats), svgpc_info(model)$decoded_variants == decoded,
          decoded == p, identical(R3$spatial$centres, original$centres))
unlink(paste0(prefix, c(".bim", ".fam")))
cat("GP three-stage fixed-knots, matching, aggregation and reusable-PC tests passed.\n")
