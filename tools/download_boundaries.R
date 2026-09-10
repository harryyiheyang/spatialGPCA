# Run from the package source directory. Only the small coordinate files are retained.
library(sf)
for (country in c("US", "UK")) {
    resolution <- if (country == "US") "110m" else "10m"
    file <- paste0("ne_", resolution, "_admin_0_countries")
    url <- paste0("https://naturalearth.s3.amazonaws.com/", resolution, "_cultural/", file, ".zip")
    archive <- tempfile(fileext = ".zip")
    download.file(url, archive, mode = "wb")
    folder <- tempfile()
    unzip(archive, exdir = folder)
    world <- st_read(file.path(folder, paste0(file, ".shp")), quiet = TRUE)
    name <- if (country == "US") "United States of America" else "United Kingdom"
    ll <- st_coordinates(st_geometry(world[world$ADMIN == name, ]))[, 1:2]
    if (country == "US") ll <- ll[ll[, 1] > -130 & ll[, 1] < -60 & ll[, 2] > 24 & ll[, 2] < 51, ]
    if (country == "UK") ll <- ll[ll[, 1] > -9 & ll[, 1] < 2 & ll[, 2] > 49 & ll[, 2] < 61, ]
    origin <- if (country == "US") c(38, -97) else c(55, -3)
    crs <- paste0("+proj=aeqd +lat_0=", origin[1], " +lon_0=", origin[2],
                   " +R=6371008.8 +units=km")
    points <- st_as_sf(data.frame(Lon = ll[, 1], Lat = ll[, 2]), coords = c("Lon", "Lat"), crs = 4326)
    xy <- st_coordinates(st_transform(points, crs))
    angle <- 2 * pi * (0:15) / 16
    normals <- cbind(cos(angle), sin(angle))
    h <- apply(xy %*% t(normals), 2, max)
    vertices <- matrix(0, 16, 2)
    for (i in 1:16) {
        j <- i %% 16 + 1L
        vertices[i, ] <- solve(normals[c(i, j), ], h[c(i, j)])
    }
    polygon <- st_sfc(st_polygon(list(rbind(vertices, vertices[1, ]))), crs = crs)
    ll <- st_coordinates(st_transform(polygon, 4326))[1:16, 1:2]
    write.table(unique(data.frame(Lat = round(ll[, 2], 7), Lon = round(ll[, 1], 7))),
                file.path("inst", "extdata", paste0(tolower(country), "_boundary.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
}
