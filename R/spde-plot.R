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

plot.svgpc_spde_mesh <- function(x, show_vertices = TRUE, vertex_cex = 0.45,
                                 main = NULL, xlab = "x (km)", ylab = "y (km)", ...) {
    if (is.null(main)) main <- paste("SPDE mesh:", nrow(x$xy), "vertices /", nrow(x$tv), "triangles")
    tv <- x$tv
    edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(3, 1)])
    edges <- unique(cbind(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2])))
    graphics::plot(x$xy, type = "n", asp = 1, xlab = xlab, ylab = ylab, main = main, ...)
    graphics::segments(x$xy[edges[, 1], 1], x$xy[edges[, 1], 2],
                       x$xy[edges[, 2], 1], x$xy[edges[, 2], 2], col = "grey65", lwd = 0.6)
    e <- x$boundary_edges
    graphics::segments(x$xy[e[, 1], 1], x$xy[e[, 1], 2],
                       x$xy[e[, 2], 1], x$xy[e[, 2], 2], col = "#176187", lwd = 1.5)
    if (show_vertices) graphics::points(x$xy, pch = 16, cex = vertex_cex, col = "#17384D")
    invisible(x)
}
