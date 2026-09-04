.object_read_result <- function(original, a, labels, ids, groups, tidy, gds,
  write, verbose, source, chrom = NULL, pos = NULL) {
  .export_flag(tidy, "tidy"); .export_flag(gds, "gds"); .export_flag(write, "write")
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Individual names must be unique and non-missing.")
  markers <- names(labels)
  if (!length(markers) || anyNA(markers) || any(!nzchar(markers)) || anyDuplicated(markers))
    stop("Locus names must be unique and non-missing.")
  if (any(lengths(labels) > 999L)) stop("Tidy GT supports up to 999 allele codes per locus.")
  if (is.null(groups)) groups <- rep("pop", length(ids))
  if (anyNA(groups)) stop("Population assignments contain missing values.")
  if (!tidy && !gds) return(original)
  n <- length(ids); m <- length(markers)
  meta <- tibble::tibble(MARKERS = markers, CHROM = if (is.null(chrom))
    rep("DENOVO", m) else as.character(chrom), LOCUS = markers,
    POS = if (is.null(pos)) seq_len(m) else as.integer(pos),
    REF = vapply(labels, function(z) z[1], ""),
    ALT = vapply(labels, function(z) if (length(z)>1) paste(z[-1],collapse=",") else ".", ""),
    VARIANT_ID = seq_len(m))
  gt <- matrix("000000", n, m)
  odd <- seq.int(1L, ncol(a), 2L)
  aa <- a[, odd, drop = FALSE]; bb <- a[, odd+1L, drop = FALSE]
  aa[is.na(aa)] <- 0L; bb[is.na(bb)] <- 0L
  gt[] <- paste0(sprintf("%03d", aa), sprintf("%03d", bb))
  out <- tibble::tibble(INDIVIDUALS = rep(ids, m), STRATA = rep(as.character(groups), m),
    MARKERS = rep(markers, each = n), GT = as.vector(gt))
  out <- dplyr::left_join(out, meta, by = "MARKERS")
  if (all(lengths(labels) == 2L)) {
    d <- a[, odd, drop = FALSE] + a[, odd+1L, drop = FALSE] - 2L
    out$ALT_DOSAGE <- as.vector(d)
  }
  attr(out, "allele_dictionary") <- labels
  if (write) saveRDS(out, tempfile(pattern = paste0(source, "_"), tmpdir = getwd(),
    fileext = ".rds"))
  if (tidy) return(out)
  if (any(lengths(labels) != 2L))
    stop("GDS conversion here requires biallelic loci; use tidy=TRUE for multiallelic data.")
  geno <- array(NA_integer_, c(2L, n, m))
  for (j in seq_len(m)) geno[, , j] <- t(a[, (2*j-1):(2*j), drop = FALSE]) - 1L
  genometranslator::genome_gds(genotypes = geno, biallelic = TRUE,
    data.source = source, strata = tibble::tibble(INDIVIDUALS = ids,
      STRATA = as.character(groups)), markers.meta = meta, verbose = verbose)
}
