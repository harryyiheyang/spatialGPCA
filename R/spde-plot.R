print.svgpc_spde_mesh <- function(x, ...) {
    cat("SPDE finite-element mesh:", nrow(x$xy), "vertices (knots),",
        nrow(x$tv), "triangles\n")
    cat("Coordinates:", x$coordinate_system$name, "; units:", x$units, "\n")
    cat("Boundary:", x$domain, "; maximum area:", x$max_area, "km^2\n")
    invisible(x)
}

print.svgpc_spde_parameters <- function(x, ...) {
    cat("SPDE", x$selector, "lambda:", format(x$lambda), ";", x$boundary, "\n")
    cat("All-SNP count:", x$p, "; full estimable rank:", x$rank, "\n")
    invisible(x)
}

print.svgpc_spde_model <- function(x, ...) {
    cat("SPDE model:", x$n, "locations;", x$m, "mesh vertices;", length(x$d), "estimable directions\n")
    cat("Conditioning:", x$conditioning, "; accumulation:", x$accumulation, "; kappa:", x$kappa, "km^-1\n")
    cat("Statistics:", if (x$ready) paste(x$p, "SNPs accumulated") else "not yet accumulated", "\n")
    invisible(x)
}

print.svgpc_spde_locations <- function(x, ...) {
    cat("SPDE observation metadata:", nrow(x$table), "individuals at", nrow(x$locations), "locations\n")
    cat("Coordinates:", x$projection$name, "; units:", x$projection$units, "\n")
    invisible(x)
}

plot.svgpc_spde_mesh <- function(x, show_vertices = TRUE, vertex_cex = .3,
                                 main = NULL, xlab = NULL, ylab = NULL,
                                 kappa = NULL, locations = x$locations, ...) {
    if (inherits(locations, "svgpc_spde_locations")) locations <- locations$locations
    if (is.null(locations)) stop("Supply locations from svgpc_spde_locations() for this saved mesh")
    vertex_cex <- .spde_positive(vertex_cex, "vertex_cex")
    radius <- if (is.null(kappa)) NULL else svgpc_spde_scale(1, correlation = .1)$kappa /
        .spde_positive(kappa, "kappa (km^-1)")
    if (is.null(main)) main <- paste("SPDE:", nrow(x$xy), "vertices")
    tv <- x$tv
    edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(3, 1)])
    edges <- unique(cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2])))
    p <- .spatial_map(locations, x$xy, x$coordinate_system, radius = radius,
                       edges = edges,
                       node_size = if (show_vertices) vertex_cex else 0,
                       main = main, xlab = xlab, ylab = ylab)
    print(p)
    invisible(p)
}
