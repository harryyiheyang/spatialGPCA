utils::globalVariables(c("long", "lat", "group", "Lon", "Lat", "density", "distance", "correlation"))

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
    grid$density <- as.vector(H) / max(H)
    grid
}

.rho_circle <- function(x, density) {
    lon <- sort(unique(density$Lon))
    lat <- sort(unique(density$Lat))
    H <- matrix(density$density, length(lon), length(lat))
    # Interpolate the displayed individual density at each actual centre.
    u <- pmax(1, pmin(length(lon), 1 + (x$centres$Lon - lon[1]) / diff(lon)[1]))
    v <- pmax(1, pmin(length(lat), 1 + (x$centres$Lat - lat[1]) / diff(lat)[1]))
    i <- pmin(floor(u), length(lon) - 1L)
    j <- pmin(floor(v), length(lat) - 1L)
    a <- u - i
    b <- v - j
    value <- (1-a) * (1-b) * H[cbind(i, j)] + a * (1-b) * H[cbind(i+1, j)] +
             (1-a) * b * H[cbind(i, j+1)] + a * b * H[cbind(i+1, j+1)]
    centre <- which.max(value)
    angle <- seq(0, 2*pi, length.out = 721)
    xy <- x$rho * cbind(cos(angle), sin(angle))
    xy <- sweep(xy, 2, as.numeric(x$centres[centre, c("x", "y")]), "+")
    list(path = .geo_inverse(xy, x$projection$origin), centre = centre)
}

plot.svgpc_cluster <- function(x, ...) {
    if (length(x$rho) != 1L || !is.finite(x$rho) || x$rho <= 0)
        stop("rho must be a positive distance in km")
    map <- .map_frame(x)
    loc <- x$locations
    centres <- x$centres
    density <- .individual_density(loc, map$limits)
    centre_size <- max(.8, min(3.5, 3.5 * (24 / nrow(centres))^.25))
    # Distance support is stable when rho is edited, so before/after curves compare directly.
    extent <- sqrt(diff(range(loc$x))^2 + diff(range(loc$y))^2)
    circle <- .rho_circle(x, density)
    view <- list(x = range(map$limits$x, range(circle$path$Lon) + c(-1, 1) * .02 * diff(map$limits$x)),
                 y = range(map$limits$y, range(circle$path$Lat) + c(-1, 1) * .02 * diff(map$limits$y)))
    distance <- seq(0, extent, length.out = 401)
    curve <- data.frame(distance = distance, correlation = exp(-distance / x$rho))
    theme <- ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", colour = "#172B3A"),
        axis.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 8),
        legend.key.size = grid::unit(3, "mm"),
        legend.spacing.x = grid::unit(2, "mm"),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    p1 <- ggplot2::ggplot() +
        ggplot2::geom_polygon(data = map$data, ggplot2::aes(x = long, y = lat, group = group),
                              fill = "#F1F3F5", colour = NA) +
        ggplot2::geom_raster(data = density, ggplot2::aes(x = Lon, y = Lat, fill = density),
                            interpolate = TRUE) +
        ggplot2::geom_polygon(data = map$data, ggplot2::aes(x = long, y = lat, group = group),
                              fill = NA, colour = "#A9B7C1", linewidth = .25) +
        ggplot2::geom_point(data = centres, ggplot2::aes(x = Lon, y = Lat, shape = "Centres"),
                            colour = "#C73535", size = centre_size, stroke = 0) +
        ggplot2::scale_shape_manual(values = c("Centres" = 16), name = NULL) +
        ggplot2::scale_fill_gradientn(colours = c("#F7F9FA", "#CDE3ED", "#71AEC7", "#267392", "#12445C"),
                                      limits = c(0, 1), breaks = c(0, 1), labels = c("Low", "High"),
                                      name = "Density") +
        ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 2)),
                         fill = ggplot2::guide_colourbar(barwidth = grid::unit(18, "mm"),
                                                        barheight = grid::unit(2, "mm"))) +
        ggplot2::geom_path(data = circle$path, ggplot2::aes(x = Lon, y = Lat),
                           colour = "#172B3A", linewidth = .8) +
        ggplot2::coord_quickmap(xlim = view$x, ylim = view$y, expand = FALSE) +
        ggplot2::labs(title = "Geographic clusters", x = NULL, y = NULL) + theme +
        ggplot2::theme(legend.position = "bottom", panel.grid.major = ggplot2::element_blank(),
                         panel.background = ggplot2::element_rect(fill = "#F7F9FA", colour = NA))
    p2 <- ggplot2::ggplot(curve, ggplot2::aes(x = distance, y = correlation)) +
        ggplot2::geom_hline(yintercept = exp(-1), colour = "#BDC7CC", linetype = 3) +
        ggplot2::geom_line(colour = "#176B87", linewidth = 1.2) +
        ggplot2::annotate("segment", x = x$rho, xend = x$rho,
                           y = 0, yend = exp(-1),
                           colour = "#172B3A", linetype = 2, linewidth = .5) +
        ggplot2::annotate("point", x = x$rho, y = exp(-1),
                           colour = "#172B3A", size = 2.5) +
        ggplot2::coord_cartesian(xlim = c(0, extent), ylim = c(0, 1), expand = FALSE) +
        ggplot2::labs(title = paste0("Correlation decay (rho = ", format(x$rho, trim = TRUE), " km)"),
                       x = NULL, y = NULL) + theme
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = c(1.8, 1))))
    print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid::popViewport()
    invisible(list(map = p1, correlation = p2))
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
