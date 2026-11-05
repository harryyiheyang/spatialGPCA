# Local test of the joint U + SPDE fit on a real PLINK BED.
# Usage: Rscript joint_local_test.R <bed prefix> [individual table TSV] [k] [max SNPs]
# The table needs FID, IID, Lat, Lon as its first four columns. Without a table,
# UK coordinates are simulated at 2,000 pseudo-locations: this checks speed and
# memory only, not the spatial result. Frequencies are read from <prefix>.afreq
# (plink2 --bfile <prefix> --freq --out <prefix>).
library(spatialGPCA)
args <- commandArgs(TRUE)
prefix <- sub("\\.bed$", "", path.expand(args[1]))
k <- if (length(args) >= 3) as.integer(args[3]) else 10L
fam <- read.table(paste0(prefix, ".fam"), colClasses = "character")
if (length(args) >= 2 && nzchar(args[2]) && file.exists(args[2])) {
    table <- read.delim(args[2], colClasses = c(FID = "character", IID = "character"))
} else {
    set.seed(1)
    pts <- data.frame(Lat = runif(2000, 50.5, 57.5), Lon = runif(2000, -5, 1))
    i <- sample(nrow(pts), nrow(fam), TRUE)
    table <- data.frame(FID = fam$V1, IID = fam$V2, Lat = pts$Lat[i], Lon = pts$Lon[i])
    message("No table supplied: simulated coordinates (timing test only)")
}
locations <- svgpc_spde_locations(table, bed = prefix)
# Real coordinates: the UK convex envelope. Simulated ones: a rectangle around the locations,
# because simulated points can fall outside the coastline envelope.
mesh <- if (length(args) >= 2 && nzchar(args[2]) && file.exists(args[2])) {
    boundary <- read.delim(system.file("extdata", "uk_boundary.tsv", package = "spatialGPCA"))
    svgpc_spde_mesh(locations, max_area = 1500, buffer = 5, boundary = boundary)
} else svgpc_spde_mesh(locations, max_area = 1500, buffer = 30)
print(mesh)
scale <- svgpc_spde_scale(distance = 100)
t0 <- proc.time()[[3]]
model <- svgpc_joint_prepare(mesh, table, kappa = scale$kappa, bed = prefix)
print(model)
reader <- NULL; snps <- NULL
if (length(args) >= 4) {
    # Limit the number of SNPs for a quick first run.
    snps <- as.integer(args[4])
    bed <- spatialGPCA:::.joint_bed(model$spatial, NULL)
    reader <- function(st, en) spatialGPCA:::cpp_joint_bed_block(bed$path, bed$file_samples, bed$rows,
                                                                st - 1L, en - st + 1L, bed$af, 0L)
}
fit <- svgpc_joint_fit(model, k = k, r = 10, reader = reader, snps = snps, max_iter = 15)
print(fit)
print(fit$trace)
cat("Prepare + fit seconds:", proc.time()[[3]] - t0, "; threads:", spatialGPCA:::cpp_joint_threads(),
    "; kernel:", fit$kernel, "\n")
out <- paste0(prefix, "_joint")
write.table(fit$individual, paste0(out, "_individual.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(fit[c("trace", "lambda", "reml", "U_eigenvalues", "V_eigenvalues", "seconds")], paste0(out, "_summary.rds"))
