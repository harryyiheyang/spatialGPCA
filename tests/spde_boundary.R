library(spatialGPCA)
xy <- rbind(c(-10, -10), c(10, -10), c(10, 10), c(-10, 10))
ll <- spatialGPCA:::.geo_inverse(xy, c(55, -3))
loc <- svgpc_spde_locations(data.frame(FID = 1:4, IID = 1:4, ll))
boundary <- spatialGPCA:::.geo_inverse(xy, loc$projection$origin)
mesh <- svgpc_spde_mesh(loc, max_area = 20, buffer = 5, boundary = boundary)
v <- mesh$parameters$vertices
stopifnot(max(abs(abs(v) - 15)) < 1e-7)
reversed <- svgpc_spde_mesh(loc, max_area = 20, buffer = 5,
                           boundary = boundary[c(4:1, 4), ])
stopifnot(identical(mesh$tv, reversed$tv), identical(mesh$xy, reversed$xy))
B <- svgpc_spde_basis(mesh, as.matrix(loc$locations[, c("x", "y")]), kappa = .2)$B
stopifnot(max(abs(Matrix::rowSums(B) - 1)) < 1e-10,
          !("sf" %in% loadedNamespaces()))
