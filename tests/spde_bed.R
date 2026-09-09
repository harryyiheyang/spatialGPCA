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
    mesh <- svgpc_spde_mesh(svgpc_spde_locations(table), max_area = 300, buffer = 30, min_angle = 21)
    fixed <- serialize(mesh, NULL)
    loc <- svgpc_spde_locations(table, bed = file.path(folder, name))
    stopifnot(identical(loc$sample_location, cls$sample_location),
              identical(loc$table, cls$table),
              max(abs(as.matrix(loc$locations[, c("x", "y", "count")]) -
                      as.matrix(cls$locations[, c("x", "y", "count")]))) < 1e-12)
    ref <- svgpc_spde_prepare_locations(loc, mesh, svgpc_spde_scale(20)$kappa)
    ref <- svgpc_spde_accumulate_matrix(ref, Y, 17L)
    for (flip in c(FALSE, TRUE)) {
        model <- svgpc_spde_prepare_data(mesh, table, svgpc_spde_scale(20)$kappa, bed = file.path(folder, name))
        frequencies <- if (flip) data.frame(ID = ff$ID, A1 = ff$REF, A2 = ff$ALT, AF = 1 - ff$ALT_FREQS) else file.path(folder, "frequencies.afreq")
        model <- svgpc_spde_accumulate(model, frequencies, block_size = 23L)
        err <- relative(model$cov, ref$cov)
        cat(name, "intersect BED; frequency flip", flip, "covariance error", err, "\n")
        stopifnot(err < 1e-11, model$missing == if (name == "missing") 2 else 0)
        stopifnot(identical(fixed, serialize(model$mesh, NULL)))
        model$spatial$bed <- tempfile("deliberately_absent_bed")
        parameters <- svgpc_spde_select_lambda(model)
        R1 <- svgpc_spde_fit(model, parameters = parameters, components = 2)
        R2 <- svgpc_spde_fit(model, parameters = parameters, components = 3)
        stopifnot(identical(R1$lambda, R2$lambda), identical(R1$lambda, parameters$lambda))
    }
}
