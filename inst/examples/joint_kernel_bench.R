# Thread scaling of the fused 2-bit kernel on a real BED, independent of the fit.
# Usage: Rscript joint_kernel_bench.R <bed prefix> [SNPs, default 2048] [threads list, default 1,2,4,8,16,28]
library(spatialGPCA)
args <- commandArgs(TRUE)
prefix <- sub("\\.bed$", "", path.expand(args[1]))
p <- if (length(args) >= 2) as.integer(args[2]) else 2048L
ths <- if (length(args) >= 3) as.integer(strsplit(args[3], ",")[[1]]) else c(1L, 2L, 4L, 8L, 16L, 28L)
ns <- length(readLines(paste0(prefix, ".fam")))
cat("Samples:", ns, "; SNPs:", p, "; omp max threads:", spatialGPCA:::cpp_joint_threads(), "\n")
t0 <- proc.time()[[3]]
cache <- spatialGPCA:::cpp_joint_bed_read(paste0(prefix, ".bed"), ns, 0:(p - 1L))
cat("Cache read:", round(proc.time()[[3]] - t0, 2), "s\n")
set.seed(1)
N <- ns; rows <- seq_len(N); loc <- rep_len(seq_len(2000L), N); af <- rep(0.3, p)
c <- 31L; q <- 11L; Yt <- matrix(rnorm(c * N), c, N)
block <- 256L
for (th in ths) {
    RtT <- matrix(0, c - q, N)
    tm <- system.time(for (st in seq(1L, p, by = block)) {
        en <- min(p, st + block - 1L)
        spatialGPCA:::cpp_joint_bed_scan("", ns, rows, (st:en) - 1L, af[st:en], loc, 2000L, Yt, q, RtT, th, cache, st - 1)
    })
    cat(sprintf("threads=%2d  elapsed=%6.2fs  cpu=%7.2fs  cores busy=%.1f\n", th, tm[["elapsed"]],
                tm[["user.self"]] + tm[["sys.self"]], (tm[["user.self"]] + tm[["sys.self"]]) / tm[["elapsed"]]))
}
