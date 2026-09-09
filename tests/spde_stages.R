library(spatialGPCA)
if (requireNamespace("RTriangle", quietly = TRUE)) {
    set.seed(112)
    table <- data.frame(FID = paste0("s", 1:120), IID = paste0("s", 1:120),
                         Lat = runif(120, 35, 36), Lon = runif(120, -100, -99), PC1 = rnorm(120))
    locations <- svgpc_spde_locations(table)
    mesh <- svgpc_spde_mesh(locations, max_area = 100, buffer = 20)
    fixed <- serialize(mesh, NULL)
    plotfile <- tempfile(fileext = ".pdf")
    grDevices::pdf(plotfile)
    plot(mesh)
    grDevices::dev.off()
    stopifnot(file.info(plotfile)$size > 0, mesh$vertex_count == nrow(mesh$xy),
              mesh$triangle_count == nrow(mesh$tv), identical(fixed, serialize(mesh, NULL)))
    model <- svgpc_spde_prepare_data(mesh, table[-(1:5), ], kappa = 0.1)
    stopifnot(identical(model$mesh, mesh),
              identical(model$spatial$projection$origin, locations$projection$origin), model$n == 115)
    G <- matrix(rnorm(model$n * 45), model$n, 45)
    G <- G + outer(model$spatial$locations$x / 10, rnorm(45))
    model <- svgpc_spde_accumulate_matrix(model, G, 8)
    parameters <- svgpc_spde_select_lambda(model)
    R1 <- svgpc_spde_fit(model, parameters = parameters, components = 2)
    R2 <- svgpc_spde_fit(model, parameters = parameters, components = 3)
    stopifnot(identical(R1$lambda, parameters$lambda), identical(R2$lambda, parameters$lambda),
              parameters$p == ncol(G), parameters$rank == length(model$d),
              identical(fixed, serialize(mesh, NULL)), model$p == ncol(G))
    if (is.finite(parameters$lambda)) stopifnot(
        length(R1$outputs$raw_F_PCA$eigenvalues) == 2,
        length(R2$outputs$raw_F_PCA$eigenvalues) == 3,
        max(abs(R1$outputs$raw_F_PCA$eigenvalues - R2$outputs$raw_F_PCA$eigenvalues[1:2])) < 1e-9)
    unlink(plotfile)
    cat("SPDE mesh plot, fixed geometry and separate lambda/PCA stages passed.\n")
}
