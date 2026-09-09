library(spatialGPCA)
merge_small <- getFromNamespace(".merge_small_clusters", "spatialGPCA")

# Mean 3 gives floor(3/2) = 1: merge the single-coordinate class.
xy <- cbind(c(0, .1, .2, .3, 10, 10.1, 10.2, 10.3, 5), 0)
centres <- rbind(c(.15, 0), c(10.15, 0), c(5, 0))
result <- merge_small(xy, centres, c(rep(1L, 4), rep(2L, 4), 3L))
stopifnot(identical(result$cluster, c(rep(1L, 4), rep(2L, 4), 1L)),
          result$threshold == 1,
          max(abs(result$centres - rbind(c(1.12, 0), c(10.15, 0)))) < 1e-12)

# Mean 4 gives 2: the two-coordinate class is included at the boundary.
xy2 <- cbind(c(0, .1, seq(10, 10.4, .1), seq(30, 30.4, .1)), 0)
centres <- rbind(c(.05, 0), c(10.2, 0), c(30.2, 0))
result <- merge_small(xy2, centres, rep(1:3, c(2, 5, 5)))
stopifnot(identical(result$cluster, c(rep(1L, 7), rep(2L, 5))), result$threshold == 2,
          max(abs(result$centres - rbind(c(7.3, 0), c(30.2, 0)))) < 1e-12)

# A mean below 2 has threshold zero; nonempty classes stay.
xy2 <- cbind(c(0, 1, 10, 11, 12), 0)
centres <- rbind(c(0, 0), c(1, 0), c(11, 0))
result <- merge_small(xy2, centres, c(1L, 2L, 3L, 3L, 3L))
stopifnot(result$threshold == 0, identical(result$centres, centres))

# Replicating people at one coordinate does not change geographic merging.
table <- data.frame(FID = 1:9, IID = 1:9, Lat = 35, Lon = -100 + xy[, 1] / 10, GPC1 = 1:9)
a <- svgpc_cluster(table, rho = 30, clusters = 3)
more <- rbind(table, table[rep(1, 200), ])
more$FID <- more$IID <- seq_len(nrow(more))
b <- svgpc_cluster(more, rho = 30, clusters = 3)
stopifnot(identical(a$centres, b$centres), identical(a$locations$cluster, b$locations$cluster),
          a$clustering$merge_at_most == b$clustering$merge_at_most,
          nrow(b$locations) == 9L, nrow(b$table) == 209L, sum(b$locations$count) == 209L)
