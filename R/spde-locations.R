svgpc_spde_locations <- function(table, bed = NULL, projection = NULL) {
    dat <- .individuals(table, bed)
    table <- dat$table
    key <- paste(sprintf("%.17g", table$Lat), sprintf("%.17g", table$Lon), sep = "\t")
    loc <- match(key, unique(key))
    ll <- table[!duplicated(key), c("Lat", "Lon"), drop = FALSE]
    n <- nrow(ll)
    if (n < 3L) stop("At least three distinct locations are needed")
    if (is.null(projection)) {
        rad <- pi / 180
        z <- c(mean(cos(ll$Lat * rad) * cos(ll$Lon * rad)),
               mean(cos(ll$Lat * rad) * sin(ll$Lon * rad)), mean(sin(ll$Lat * rad)))
        origin <- c(atan2(z[3], sqrt(sum(z[1:2]^2))), atan2(z[2], z[1])) / rad
    } else {
        if (!identical(projection$name, "Spherical azimuthal equidistant") ||
            !identical(projection$units, "km") || length(projection$origin) != 2L ||
            any(!is.finite(projection$origin))) stop("Provide the original spherical geographic projection in km")
        origin <- projection$origin
    }
    xy <- .geo_forward(ll$Lat, ll$Lon, origin)
    counts <- tabulate(loc, n)
    dat$locations <- data.frame(ll, xy, count = counts)
    rownames(dat$locations) <- as.character(seq_len(n))
    dat$sample_location <- loc
    dat$projection <- list(name = "Spherical azimuthal equidistant", origin = origin, units = "km")
    class(dat) <- "svgpc_spde_locations"
    dat
}

svgpc_spde_prepare_data <- function(mesh, table, kappa, bed = NULL, tau = 1,
                                   conditioning = c("qr", "profile"),
                                   accumulation = c("auto", "moments", "directions")) {
    if (!inherits(mesh, "svgpc_spde_mesh")) stop("Select a mesh from svgpc_spde_mesh() first")
    projection <- mesh$coordinate_system
    if (is.null(projection$origin))
        stop("For geographic tables, build the mesh from svgpc_spde_locations(table); planar meshes require explicit location-level inputs")
    spatial <- svgpc_spde_locations(table, bed = bed, projection = projection)
    svgpc_spde_prepare_locations(spatial, mesh, kappa, tau, match.arg(conditioning), match.arg(accumulation))
}
