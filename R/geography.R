.geo_forward <- function(lat, lon, origin) {
    rad <- pi / 180
    a <- lat * rad
    b <- (lon - origin[2]) * rad
    a0 <- origin[1] * rad
    z <- pmax(-1, pmin(1, sin(a0) * sin(a) + cos(a0) * cos(a) * cos(b)))
    d <- acos(z)
    if (any(d > pi - 1e-6)) stop("A location is antipodal to the map centre; use a regional dataset")
    k <- rep(1, length(d))
    i <- d > 1e-10
    k[i] <- d[i] / sin(d[i])
    cbind(x = 6371.0088 * k * cos(a) * sin(b),
          y = 6371.0088 * k * (cos(a0) * sin(a) - sin(a0) * cos(a) * cos(b)))
}

.geo_inverse <- function(xy, origin) {
    d <- sqrt(rowSums(xy^2))
    a <- d / 6371.0088
    a0 <- origin[1] * pi / 180
    lat <- rep(a0, nrow(xy))
    lon <- rep(origin[2] * pi / 180, nrow(xy))
    i <- d > 1e-10
    lat[i] <- asin(pmax(-1, pmin(1, cos(a[i]) * sin(a0) +
                     xy[i, 2] * sin(a[i]) * cos(a0) / d[i])))
    lon[i] <- lon[i] + atan2(xy[i, 1] * sin(a[i]),
                             d[i] * cos(a0) * cos(a[i]) - xy[i, 2] * sin(a0) * sin(a[i]))
    data.frame(Lat = lat * 180 / pi, Lon = (lon * 180 / pi + 180) %% 360 - 180)
}

.individuals <- function(table, bed) {
    if (!is.data.frame(table) && !is.matrix(table)) stop("table must be an individual data frame or matrix")
    table <- as.data.frame(table, stringsAsFactors = FALSE)
    if (ncol(table) < 4L) stop("The first four columns must be FID, IID, Lat, Lon")
    names(table)[1:4] <- c("FID", "IID", "Lat", "Lon")
    table$FID <- as.character(table$FID)
    table$IID <- as.character(table$IID)
    key <- .sample_key(table)
    total <- nrow(table)
    file_rows <- seq_len(total)
    file_n <- total
    if (!is.null(bed)) {
        bed <- sub("\\.bed$", "", bed, ignore.case = TRUE)
        fam <- .metadata_table(paste0(bed, ".fam"), FALSE)
        if (ncol(fam) != 6L) stop("FAM must have six columns")
        names(fam)[1:2] <- c("FID", "IID")
        idx <- match(.sample_key(fam), key)
        file_rows <- which(!is.na(idx))
        file_n <- nrow(fam)
        if (!length(file_rows)) stop("No shared FID/IID pairs between table and FAM")
        table <- table[idx[file_rows], , drop = FALSE]
        bed <- normalizePath(bed, mustWork = FALSE, winslash = "/")
    }
    if (!is.numeric(table$Lat) || !is.numeric(table$Lon) ||
        any(!is.finite(table$Lat)) || any(!is.finite(table$Lon)) ||
        any(abs(table$Lat) > 90) || any(abs(table$Lon) > 360))
        stop("Lat and Lon must contain finite coordinates in degrees")
    table$Lon <- (table$Lon + 180) %% 360 - 180
    X <- as.matrix(table[, -c(1:4), drop = FALSE])
    if (ncol(X) && (!is.numeric(X) || any(!is.finite(X))))
        stop("All columns after FID, IID, Lat, Lon must be finite numeric fixed effects")
    rownames(table) <- NULL
    list(table = table, bed = bed, file_rows = file_rows, file_n = file_n,
         matched = c(table = total, fam = file_n, retained = nrow(table)))
}

.merge_small_clusters <- function(xy, centres, cluster) {
    m <- nrow(centres)
    size <- tabulate(cluster, m)
    threshold <- floor(mean(size) / 2)
    sums <- rowsum(xy, cluster, reorder = TRUE)
    parent <- seq_len(m)
    active <- rep(TRUE, m)
    small <- which(active & size <= threshold)
    while (length(small)) {
        from <- small[which.min(size[small])]
        d <- (centres[, 1] - centres[from, 1])^2 + (centres[, 2] - centres[from, 2])^2
        d[!active | seq_len(m) == from] <- Inf
        to <- which.min(d)
        size[to] <- size[to] + size[from]
        sums[to, ] <- sums[to, ] + sums[from, ]
        centres[to, ] <- sums[to, ] / size[to]
        parent[parent == from] <- to
        active[from] <- FALSE
        small <- which(active & size <= threshold)
    }
    keep <- which(active)
    list(centres = centres[keep, , drop = FALSE],
         cluster = match(parent[cluster], keep), threshold = threshold)
}

# The geographic object contains no fitted smoothing penalty.
svgpc_cluster <- function(table, rho, clusters = NULL, bed = NULL, seed = 44L) {
    if (length(rho) != 1L || !is.finite(rho) || rho <= 0) stop("rho must be a positive distance in km")
    dat <- .individuals(table, bed)
    table <- dat$table
    key <- paste(sprintf("%.17g", table$Lat), sprintf("%.17g", table$Lon), sep = "\t")
    loc <- match(key, unique(key))
    ll <- table[!duplicated(key), c("Lat", "Lon"), drop = FALSE]
    n <- nrow(ll)
    if (n < 3L) stop("At least three distinct locations are needed")
    if (is.null(clusters)) clusters <- min(n, 50L)
    clusters <- .positive_integer(clusters, "clusters")
    if (clusters < 2L || clusters > n) stop("clusters must be between 2 and the number of distinct locations")
    rad <- pi / 180
    z <- c(mean(cos(ll$Lat * rad) * cos(ll$Lon * rad)),
           mean(cos(ll$Lat * rad) * sin(ll$Lon * rad)), mean(sin(ll$Lat * rad)))
    origin <- c(atan2(z[3], sqrt(sum(z[1:2]^2))), atan2(z[2], z[1])) / rad
    xy <- .geo_forward(ll$Lat, ll$Lon, origin)
    cls <- cpp_cluster(xy, clusters, as.integer(seed))
    nonempty <- nrow(cls$centres)
    counts <- tabulate(loc, n)
    merged <- .merge_small_clusters(xy, cls$centres, cls$cluster)
    cls$centres <- merged$centres
    cls$cluster <- merged$cluster
    colnames(cls$centres) <- c("x", "y")
    rownames(xy) <- as.character(seq_len(n))
    dat$locations <- data.frame(ll, xy, count = counts, cluster = cls$cluster)
    rownames(dat$locations) <- as.character(seq_len(n))
    dat$sample_location <- loc
    dat$centres <- data.frame(.geo_inverse(cls$centres, origin), cls$centres)
    dat$rho <- as.numeric(rho)
    dat$projection <- list(name = "Spherical azimuthal equidistant", origin = origin, units = "km")
    dat$clustering <- list(initial_clusters = clusters, nonempty_clusters = nonempty,
                           merge_at_most = merged$threshold, seed = seed,
                           iterations = cls$iterations, converged = cls$converged)
    class(dat) <- "svgpc_cluster"
    dat
}

print.svgpc_cluster <- function(x, ...) {
    cat("spatialGPCA geography:", nrow(x$table), "individuals at", nrow(x$locations),
        "locations;", nrow(x$centres), "clusters\n")
    if (!is.null(x$clustering$merge_at_most))
        cat("Initial clusters:", x$clustering$initial_clusters,
            "; merged clusters with <=", x$clustering$merge_at_most, "coordinate points\n")
    cat("rho:", x$rho, "km;", x$projection$name, "coordinates\n")
    cat("Fixed effects:", if (ncol(x$table) > 4L) paste(names(x$table)[-(1:4)], collapse = ", ") else "intercept only", "\n")
    invisible(x)
}
