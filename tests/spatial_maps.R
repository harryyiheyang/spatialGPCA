library(spatialGPCA)
set.seed(71)
table <- data.frame(FID = 1:80, IID = 1:80, Lat = runif(80, 30, 45), Lon = runif(80, -115, -80))
table <- rbind(table, table[rep(1:5, 5), ])
table$FID <- table$IID <- seq_len(nrow(table))
cls <- svgpc_cluster(table, rho = 300 / log(10), clusters = 10)
loc <- svgpc_spde_locations(table)
mesh <- svgpc_spde_mesh(loc, max_area = 100000, buffer = 50)
scale <- svgpc_spde_scale(300, correlation = .1)
fixed <- serialize(mesh, NULL)
pdf(tempfile(fileext = ".pdf"), width = 10, height = 6)
p1 <- plot(cls)
p2 <- plot(mesh, kappa = scale$kappa)
stopifnot(inherits(p1, "ggplot"), inherits(p2, "ggplot"), identical(p1$data, p2$data),
          identical(fixed, serialize(mesh, NULL)), sum(mesh$locations$count) == nrow(table))
for (p in list(p1, p2)) {
    circle <- attr(p, "range_circle")
    d <- sqrt(rowSums(sweep(circle$xy, 2, circle$centre, "-")^2))
    stopifnot(max(abs(d - 300)) < 1e-10, nrow(ggplot2::ggplot_build(p)$layout$layout) == 1)
}
stopifnot(abs(exp(-attr(p1, "range_circle")$radius / cls$rho) - .1) < 1e-12,
          abs(svgpc_spde_correlation(attr(p2, "range_circle")$radius, scale$kappa) - .1) < 1e-12)
points1 <- Filter(function(x) inherits(x$geom, "GeomPoint"), p1$layers)
points2 <- Filter(function(x) inherits(x$geom, "GeomPoint"), p2$layers)
stopifnot(length(points1) == 1, length(points2) == 1,
          nrow(points1[[1]]$data) == nrow(cls$centres),
          nrow(points2[[1]]$data) == nrow(mesh$xy),
          identical(points1[[1]]$aes_params$size, points2[[1]]$aes_params$size))
p3 <- plot(mesh, kappa = scale$kappa * 2)
stopifnot(identical(p2$data, p3$data), abs(attr(p3, "range_circle")$radius - 150) < 1e-10,
          identical(fixed, serialize(mesh, NULL)))
mesh0 <- svgpc_spde_mesh(as.matrix(loc$locations[, c("x", "y")]), max_area = 100000, buffer = 50)
p4 <- plot(mesh0, kappa = scale$kappa)
stopifnot(abs(attr(p4, "range_circle")$radius - 300) < 1e-10,
          all(mesh0$locations$count == 1))
dev.off()
