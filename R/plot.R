utils::globalVariables(c("long", "lat", "group", "Lon", "Lat", "density", "lon_end", "lat_end"))

.map_frame <- function(x) {
    ll <- x$locations
    usa <- all(ll$Lon >= -126 & ll$Lon <= -66 & ll$Lat >= 24 & ll$Lat <= 50)
    outline <- maps::map(if (usa) "state" else "world", plot = FALSE, fill = TRUE)
    map <- data.frame(long = outline$x, lat = outline$y,
                       group = cumsum(is.na(outline$x)) + 1L)
    map <- map[is.finite(map$long) & is.finite(map$lat), , drop = FALSE]
    pad <- c(max(diff(range(ll$Lon)) * .08, 1), max(diff(range(ll$Lat)) * .08, 1))
    limits <- list(x = range(ll$Lon) + c(-1, 1) * pad[1],
                   y = range(ll$Lat) + c(-1, 1) * pad[2])
    if (usa) limits <- list(x = c(-125.5, -66), y = c(24, 50))
    list(data = map, limits = limits, usa = usa)
}

.individual_density <- function(loc, limits) {
    # Bin counts before smoothing: work stays bounded as the sample size grows.
    nx <- 240L
    ny <- 120L
    dx <- diff(limits$x) / nx
    dy <- diff(limits$y) / ny
    ix <- pmax(1L, pmin(nx, floor((loc$Lon - limits$x[1]) / dx) + 1L))
    iy <- pmax(1L, pmin(ny, floor((loc$Lat - limits$y[1]) / dy) + 1L))
    bins <- rowsum(loc$count, ix + nx * (iy - 1L), reorder = FALSE)
    H <- matrix(0, nx, ny)
    H[as.integer(rownames(bins))] <- bins[, 1L]
    # Display smoothing is independent of rho and of the cluster count.
    Kx <- exp(-.5 * (outer(seq_len(nx), seq_len(nx), "-") / 2)^2)
    Ky <- exp(-.5 * (outer(seq_len(ny), seq_len(ny), "-") / 2)^2)
    Kx <- sweep(Kx, 2L, colSums(Kx), "/")
    Ky <- sweep(Ky, 2L, colSums(Ky), "/")
    H <- Kx %*% H %*% t(Ky)
    grid <- expand.grid(Lon = limits$x[1] + (seq_len(nx) - .5) * dx,
                        Lat = limits$y[1] + (seq_len(ny) - .5) * dy)
    grid$density <- log1p(as.vector(H)) / log1p(max(H))
    grid
}

.range_circle <- function(density, projection, radius) {
    i <- which.max(density$density)
    centre <- if (is.null(projection$origin)) c(density$Lon[i], density$Lat[i]) else
        as.numeric(.geo_forward(density$Lat[i], density$Lon[i], projection$origin))
    angle <- seq(0, 2 * pi, length.out = 721)
    xy <- sweep(radius * cbind(cos(angle), sin(angle)), 2L, centre, "+")
    path <- if (is.null(projection$origin)) data.frame(Lon = xy[, 1], Lat = xy[, 2]) else
        .geo_inverse(xy, projection$origin)
    list(path = path, xy = xy, centre = centre, radius = radius)
}

.spatial_map <- function(loc, xy, projection, radius = NULL, edges = NULL,
                         node_size = .3, main = NULL,
                         xlab = NULL, ylab = NULL) {
    geographic <- !is.null(projection$origin)
    if (geographic) {
        map <- .map_frame(list(locations = loc))
        nodes <- .geo_inverse(xy, projection$origin)
    } else {
        loc$Lon <- loc$x; loc$Lat <- loc$y
        nodes <- data.frame(Lon = xy[, 1], Lat = xy[, 2])
        pad <- pmax(c(diff(range(loc$x)), diff(range(loc$y))) * .08, .01)
        map <- list(data = NULL, limits = list(x = range(loc$x, xy[, 1]) + c(-1, 1) * pad[1],
                                               y = range(loc$y, xy[, 2]) + c(-1, 1) * pad[2]))
        if (is.null(xlab)) xlab <- "x (km)"
        if (is.null(ylab)) ylab <- "y (km)"
    }
    density <- .individual_density(loc, map$limits)
    circle <- if (is.null(radius)) NULL else .range_circle(density, projection, radius)
    view <- map$limits
    if (!is.null(edges)) {
        view$x <- range(view$x, nodes$Lon)
        view$y <- range(view$y, nodes$Lat)
    }
    if (!is.null(circle)) {
        view$x <- range(view$x, circle$path$Lon)
        view$y <- range(view$y, circle$path$Lat)
        main <- paste0(main, " | 0.1 range: ", format(round(radius, 1), trim = TRUE), " km")
    }
    p <- ggplot2::ggplot(density, ggplot2::aes(Lon, Lat)) +
        ggplot2::geom_raster(ggplot2::aes(fill = density), interpolate = TRUE) +
        ggplot2::scale_fill_gradientn(colours = c("#F7F9FA", "#CDE3ED", "#71AEC7", "#267392", "#12445C"),
                                      limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
                                      name = "Sample density (log scale)")
    if (!is.null(edges)) {
        d <- data.frame(Lon = nodes$Lon[edges[, 1]], Lat = nodes$Lat[edges[, 1]],
                        lon_end = nodes$Lon[edges[, 2]], lat_end = nodes$Lat[edges[, 2]])
        p <- p + ggplot2::geom_segment(data = d,
            ggplot2::aes(xend = lon_end, yend = lat_end), colour = "#687D89",
            linewidth = .12, alpha = .4)
    }
    if (node_size > 0) p <- p + ggplot2::geom_point(data = nodes, colour = "#C73535",
                                                   size = node_size, stroke = 0)
    if (geographic) p <- p + ggplot2::geom_polygon(data = map$data,
        ggplot2::aes(x = long, y = lat, group = group), inherit.aes = FALSE,
        fill = NA, colour = "#253C49", linewidth = .5)
    if (!is.null(circle)) p <- p + ggplot2::geom_path(data = circle$path,
        colour = "#182B35", linewidth = .8)
    p <- p + (if (geographic) ggplot2::coord_quickmap(xlim = view$x, ylim = view$y, expand = FALSE) else
        ggplot2::coord_fixed(xlim = view$x, ylim = view$y, expand = FALSE)) +
        ggplot2::labs(title = main, x = xlab, y = ylab) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom",
                       plot.title = ggplot2::element_text(face = "bold", colour = "#172B3A"),
                       panel.background = ggplot2::element_rect(fill = "#F7F9FA", colour = NA)) +
        ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "top", barwidth = grid::unit(40, "mm"),
                                                       barheight = grid::unit(2, "mm")))
    attr(p, "range_circle") <- circle
    p
}

plot.svgpc_cluster <- function(x, node_size = .3, main = NULL, ...) {
    rho <- .spde_positive(x$rho, "rho (km)")
    node_size <- .spde_positive(node_size, "node_size")
    if (is.null(main)) main <- paste("GP:", nrow(x$centres), "centres")
    p <- .spatial_map(x$locations, as.matrix(x$centres[, c("x", "y")]),
                       x$projection, radius = log(10) * rho, node_size = node_size, main = main)
    print(p)
    invisible(p)
}

plot.svgpc_fit <- function(x, ...) {
    if (is.null(x$spatial)) stop("This fit does not contain a geographic result")
    spatial <- x$spatial
    spatial$rho <- x$rho
    plot(spatial, ...)
}

print.svgpc_parameters <- function(x, ...) {
    cat("spatialGPCA", x$selector, "smoothing selection: lambda =", format(x$lambda), "\n")
    cat("Boundary:", x$boundary, "; spatial degrees of freedom:", format(x$spatial_df), "\n")
    invisible(x)
}

print.svgpc_fit <- function(x, ...) {
    cat("spatialGPCA fitted spatial subspace:", x$selector, "lambda =", format(x$lambda), "\n")
    for (name in names(x$outputs)) {
        output <- x$outputs[[name]]
        score <- if (name == "raw_F_PCA") output$pcs else output$raw_scores
        cat(name, ":", nrow(score), "locations x", ncol(score), "components\n")
    }
    invisible(x)
}
