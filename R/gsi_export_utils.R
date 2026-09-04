# Shared diploid GSI export utilities.
.gsi_loci <- function(gds) {
  ids <- SeqArray::seqGetData(gds, "variant.id")
  mm <- genometranslator::extract_markers_metadata(gds)
  definitions <- strsplit(SeqArray::seqGetData(gds, "allele"), ",", fixed = TRUE)
  definitions <- lapply(definitions, function(x) {
    if (length(x) == 2L && x[2] == ".") x <- x[1]
    if (!length(x) || anyNA(x) || any(!nzchar(x) | x == "." |
        grepl("[<>*]", x))) stop("Unsupported missing or symbolic allele definition.")
    x
  })
  loci <- tibble::tibble(EXPORT_ID = paste0("L", seq_along(ids)),
    VARIANT_ID = ids, MARKERS = if ("MARKERS" %in% names(mm))
      as.character(mm$MARKERS[match(ids, mm$VARIANT_ID)]) else as.character(ids),
    CHROM = SeqArray::seqGetData(gds, "chromosome"),
    POS = SeqArray::seqGetData(gds, "position"))
  alleles <- dplyr::bind_rows(lapply(seq_along(ids), function(j)
    tibble::tibble(EXPORT_ID = loci$EXPORT_ID[j],
      CODE = seq_along(definitions[[j]]), ALLELE = definitions[[j]])))
  list(loci = loci, alleles = alleles, counts = lengths(definitions))
}
.gsi_block <- function(gds, ids, counts) {
  SeqArray::seqSetFilter(gds, sample.id = ids, verbose = FALSE)
  gt <- .export_genotypes(gds)
  # Restore requested order: SeqArray returns samples in file order.
  rows <- match(ids, as.character(SeqArray::seqGetData(gds, "sample.id")))
  out <- matrix(NA_integer_, length(ids), 2L * length(counts))
  for (j in seq_along(counts)) {
    a <- t(matrix(gt[, , j], nrow = 2L)) + 1L
    if (any(a > counts[j], na.rm = TRUE)) stop("Allele index exceeds its dictionary.")
    if (any(xor(is.na(a[, 1]), is.na(a[, 2]))))
      stop("Partially missing diploid calls are unsupported; resolve them explicitly.")
    out[, (2L * j - 1L):(2L * j)] <- a[rows, , drop = FALSE]
  }
  out
}
.gsi_metadata <- function(gds, strata) {
  if (is.null(strata)) return(genometranslator::extract_individuals_metadata(gds))
  if (is.character(strata) && length(strata) == 1L && file.exists(strata))
    strata <- readr::read_tsv(strata, col_types = readr::cols(.default = "c"),
      progress = FALSE)
  if (!is.data.frame(strata)) stop("strata must be a metadata table or TSV file.")
  strata
}
