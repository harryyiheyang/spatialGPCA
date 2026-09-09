svgpc_spde_accumulate <- function(model, frequencies, block_size = 256L) {
    .spde_fresh(model)
    spatial <- model$spatial
    prefix <- spatial$bed
    if (is.null(prefix)) stop("Provide bed to svgpc_spde_prepare_data(mesh, table, kappa, bed = ...) before BED fitting")
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
    if (sum(!is.na(idx)) != nrow(samples)) stop("FAM changed after preparation; reprepare with the same selected mesh and updated metadata")
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
    stats <- cpp_spde_accumulate_bed(model$accumulation_basis, model$Xq, model$counts, normalizePath(path), as.integer(locations), af,
                        .positive_integer(block_size, "block_size"), model$accumulation_projection)
    model$provenance <- list(prefix = prefix, missing = "2f",
                                      file_variants = nrow(vm), file_samples = nrow(sm),
                                      matched_samples = length(keep),
                                      supplied_frequency_variants = nrow(frequencies),
                                      flipped_frequency_alleles = sum(flip))
    .spde_finish(model, stats)
}
