# R owns individual metadata and object lifecycle; numerical work is in src/.
.matrix_double <- function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
}
.model <- function(x) {
    if (!inherits(x, "svgpc_model") || typeof(x) != "externalptr")
        stop("Expected an spatialGPCA model from svgpc_prepare() or svgpc_load()")
    x
}
.positive_integer <- function(x, label) {
    if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x) || x > .Machine$integer.max)
        stop(label, " must be a positive integer")
    as.integer(x)
}

.prepare_xy <- function(coordinates, pcs, counts, knots, rho) {
    coordinates <- .matrix_double(coordinates)
    if (missing(pcs) || is.null(pcs)) pcs <- matrix(numeric(), nrow(coordinates), 0L)
    if (length(counts) != nrow(coordinates)) stop("counts must have one entry per location")
    location_ids <- rownames(coordinates)
    if (is.null(location_ids)) location_ids <- as.character(seq_len(nrow(coordinates)))
    if (anyNA(location_ids) || anyDuplicated(location_ids)) stop("Duplicate or missing location IDs")
    out <- cpp_prepare(coordinates, .matrix_double(pcs), as.numeric(counts),
                       .matrix_double(knots), rho)
    class(out) <- "svgpc_model"
    attr(out, "location_ids") <- location_ids
    out
}

svgpc_prepare <- function(spatial, table = NULL, bed = NULL) {
    if (!inherits(spatial, "svgpc_cluster")) stop("Use a geographic result from svgpc_cluster()")
    if (length(spatial$rho) != 1L || !is.finite(spatial$rho) || spatial$rho <= 0)
        stop("rho must be a positive distance in km")
    if (!is.null(table) || !is.null(bed)) {
        if (is.null(table)) table <- spatial$table
        if (is.null(bed)) bed <- spatial$bed
        dat <- .individuals(table, bed)
        table <- dat$table
        key <- paste(sprintf("%.17g", table$Lat), sprintf("%.17g", table$Lon), sep = "\t")
        loc <- match(key, unique(key))
        ll <- table[!duplicated(key), c("Lat", "Lon"), drop = FALSE]
        xy <- .geo_forward(ll$Lat, ll$Lon, spatial$projection$origin)
        old_key <- paste(sprintf("%.17g", spatial$locations$Lat),
                         sprintf("%.17g", spatial$locations$Lon), sep = "\t")
        old <- match(unique(key), old_key)
        cluster <- spatial$locations$cluster[old]
        for (i in which(is.na(old))) {
            d <- (spatial$centres$x - xy[i, 1L])^2 + (spatial$centres$y - xy[i, 2L])^2
            cluster[i] <- which.min(d)
        }
        spatial[names(dat)] <- dat
        spatial$locations <- data.frame(ll, xy, count = tabulate(loc, nrow(ll)), cluster = cluster)
        rownames(spatial$locations) <- as.character(seq_len(nrow(ll)))
        spatial$sample_location <- loc
    }
    X <- as.matrix(spatial$table[, -c(1:4), drop = FALSE])
    n <- nrow(spatial$locations)
    C <- if (ncol(X)) rowsum(X, spatial$sample_location, reorder = TRUE) /
        spatial$locations$count else matrix(numeric(), n, 0L)
    model <- .prepare_xy(spatial$locations[, c("x", "y")], C,
                         spatial$locations$count, spatial$centres[, c("x", "y")], spatial$rho)
    attr(model, "spatial") <- spatial
    model
}

# G is already standardized, location x SNP.
svgpc_accumulate_matrix <- function(model, G, block_size = 256L) {
    cpp_accumulate_matrix(.model(model), .matrix_double(G),
                          .positive_integer(block_size, "block_size"))
    invisible(model)
}

.metadata_table <- function(file, header = TRUE) {
    if (!file.exists(file)) stop("Missing metadata file: ", file)
    skip <- 0L
    if (header) {
        con <- file(file, "rt")
        on.exit(close(con))
        repeat {
            line <- readLines(con, n = 1L, warn = FALSE)
            if (!length(line)) stop("No metadata header: ", file)
            if (!startsWith(line, "##")) break
            skip <- skip + 1L
        }
    }
    ans <- utils::read.table(file, header = header, skip = skip, sep = "",
                             quote = "", comment.char = "", check.names = FALSE,
                             stringsAsFactors = FALSE, colClasses = "character")
    names(ans) <- sub("^#", "", names(ans))
    ans
}
.required <- function(x, cols, label) {
    if (!all(cols %in% names(x))) stop(label, " requires columns: ", paste(cols, collapse = ", "))
    if (anyNA(x[, cols, drop = FALSE])) stop(label, " has missing metadata")
}
.sample_key <- function(x) {
    .required(x, c("FID", "IID"), "Sample metadata")
    if (any(grepl("[\t\r\n]", x$FID)) || any(grepl("[\t\r\n]", x$IID)))
        stop("Sample IDs must not contain tabs or newlines")
    key <- paste(x$FID, x$IID, sep = "\t")
    if (anyDuplicated(key)) stop("Duplicate FID/IID pairs")
    key
}

svgpc_accumulate <- function(model, frequencies, block_size = 256L) {
    .model(model)
    spatial <- attr(model, "spatial")
    prefix <- spatial$bed
    if (is.null(prefix)) stop("Provide bed to svgpc_prepare(spatial, bed = ...) or svgpc_cluster() before BED fitting")
    path <- paste0(prefix, ".bed")
    if (!file.exists(path)) stop("Missing genotype file: ", path)
    sm <- .metadata_table(paste0(prefix, ".fam"), FALSE)
    if (ncol(sm) != 6L) stop("FAM must have six columns")
    names(sm)[1:2] <- c("FID", "IID")
    vm <- .metadata_table(paste0(prefix, ".bim"), FALSE)
    if (ncol(vm) != 6L) stop("BIM must have six columns")
    names(vm) <- c("CHROM", "ID", "CM", "POS", "COUNTED", "OTHER")
    .required(vm, c("CHROM", "ID", "COUNTED", "OTHER"), "Variant metadata")
    chrom <- sub("^chr", "", vm$CHROM, ignore.case = TRUE)
    if (!all(chrom %in% as.character(1:22)))
        stop("This version requires externally prepared autosomal variants (1..22)")
    if (anyDuplicated(vm$ID) || any(vm$ID %in% c("", "."))) stop("Variant IDs must be unique and nonmissing")
    if (any(grepl(",", vm$COUNTED, fixed = TRUE)) || any(vm$COUNTED == vm$OTHER))
        stop("This version requires biallelic variants")
    samples <- spatial$table
    input_key <- .sample_key(samples)
    file_key <- .sample_key(sm)
    idx <- match(file_key, input_key)
    if (sum(!is.na(idx)) != nrow(samples)) stop("FAM changed after clustering; rebuild the geographic result")
    locations <- integer(length(file_key))
    keep <- which(!is.na(idx))
    locations[keep] <- spatial$sample_location[idx[keep]]
    if (is.character(frequencies) && length(frequencies) == 1L)
        frequencies <- .metadata_table(frequencies)
    if (all(c("ID", "REF", "ALT", "ALT_FREQS") %in% names(frequencies))) {
        ff <- data.frame(ID = frequencies$ID, A1 = frequencies$ALT,
                         A2 = frequencies$REF, AF = frequencies$ALT_FREQS)
    } else if (all(c("SNP", "A1", "A2", "MAF") %in% names(frequencies))) {
        ff <- data.frame(ID = frequencies$SNP, A1 = frequencies$A1,
                         A2 = frequencies$A2, AF = frequencies$MAF)
    } else {
        .required(frequencies, c("ID", "A1", "A2", "AF"), "frequencies")
        ff <- frequencies
    }
    .required(ff, c("ID", "A1", "A2", "AF"), "frequencies")
    if (anyDuplicated(ff$ID)) stop("Duplicate frequency variant IDs")
    fidx <- match(vm$ID, ff$ID)
    if (anyNA(fidx)) stop("Missing external allele frequencies for file variants")
    ff <- ff[fidx, , drop = FALSE]
    same <- ff$A1 == vm$COUNTED & ff$A2 == vm$OTHER
    flip <- ff$A2 == vm$COUNTED & ff$A1 == vm$OTHER
    if (any(!same & !flip)) stop("Alleles disagree between frequency and genotype metadata")
    af <- as.numeric(as.character(ff$AF))
    if (any(!is.finite(af)) || any(af <= 0 | af >= 1))
        stop("External allele frequencies must be strictly between zero and one; no variants are silently dropped")
    af[flip] <- 1 - af[flip]
    cpp_accumulate_file(model, normalizePath(path), as.integer(locations), af,
                        .positive_integer(block_size, "block_size"))
    attr(model, "provenance") <- list(prefix = prefix, missing = "2f",
                                      file_variants = nrow(vm), file_samples = nrow(sm),
                                      matched_samples = length(keep),
                                      supplied_frequency_variants = nrow(frequencies),
                                      flipped_frequency_alleles = sum(flip))
    invisible(model)
}

svgpc_select_lambda <- function(model, method = c("REML", "GCV")) {
    result <- cpp_fit(.model(model), match.arg(method), character(), 1L, NaN)
    result[c("outputs", "sigma2_reml_profile")] <- NULL
    class(result) <- "svgpc_parameters"
    attr(result, "model") <- model
    result
}

svgpc_fit <- function(model, method = c("REML", "GCV"),
                      operators = "raw_F_PCA", components = 50L, parameters = NULL) {
    if (!length(operators) || anyNA(operators) || anyDuplicated(operators)) stop("Provide distinct PCA operators")
    if (is.null(parameters)) parameters <- svgpc_select_lambda(model, match.arg(method))
    if (!inherits(parameters, "svgpc_parameters") || !identical(attr(parameters, "model"), model))
        stop("parameters must be selected from this model by svgpc_select_lambda()")
    if (!missing(method) && match.arg(method) != parameters$selector)
        stop("method disagrees with the supplied parameter-selection result")
    result <- cpp_fit(.model(model), parameters$selector, as.character(operators),
                      .positive_integer(components, "components"),
                      parameters$lambda)
    result$selector <- parameters$selector
    result$boundary <- parameters$boundary
    result$evaluations <- parameters$evaluations
    result$sigma2_reml_profile <- NULL
    result$location_ids <- attr(model, "location_ids")
    result$spatial <- attr(model, "spatial")
    result$rho <- cpp_info(model)$rho
    class(result) <- "svgpc_fit"
    result
}

svgpc_fitted <- function(model, G, parameters, coefficients = FALSE) {
    if (!inherits(parameters, "svgpc_parameters") || !identical(attr(parameters, "model"), model))
        stop("parameters must be selected from this model by svgpc_select_lambda()")
    cpp_fitted(.model(model), .matrix_double(G), parameters$lambda, as.logical(coefficients))
}
svgpc_info <- function(model) {
    out <- cpp_info(.model(model))
    out$provenance <- attr(model, "provenance")
    out
}
svgpc_statistics <- function(model) cpp_statistics(.model(model))
svgpc_save <- function(model, file, compress = FALSE) {
    state <- cpp_snapshot(.model(model))
    state$location_ids <- attr(model, "location_ids")
    state$provenance <- attr(model, "provenance")
    state$spatial <- attr(model, "spatial")
    saveRDS(state, file = file, compress = compress)
    invisible(file)
}
svgpc_load <- function(file) {
    state <- readRDS(file)
    out <- cpp_restore(state)
    class(out) <- "svgpc_model"
    attr(out, "location_ids") <- state$location_ids
    attr(out, "provenance") <- state$provenance
    attr(out, "spatial") <- state$spatial
    out
}
