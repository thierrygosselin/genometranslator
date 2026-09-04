.legacy_snapshot <- function(data, strata = NULL, grouped = FALSE,
  biallelic = FALSE, chunk.size = 32L) {
  .export_count(chunk.size, "chunk.size", 1)
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  g <- context$gds
  ids <- as.character(SeqArray::seqGetData(g, "sample.id"))
  info <- .gsi_loci(g)
  if (biallelic && any(info$counts != 2L))
    stop("This format requires biallelic loci; filter explicitly before export.")
  meta <- if (grouped) .export_groups(g, strata, "STRATA", NULL) else NULL
  samples <- if (grouped) meta$samples else tibble::tibble(INDIVIDUALS = ids)
  samples$EXPORT_ID <- paste0("S", seq_along(ids))
  a <- matrix(NA_integer_, length(ids), 2L * nrow(info$loci))
  for (first in seq.int(1L, length(ids), by = chunk.size)) {
    r <- first:min(length(ids), first + chunk.size - 1L)
    a[r, ] <- .gsi_block(g, ids[r], info$counts)
  }
  colnames(a) <- as.vector(rbind(info$loci$EXPORT_ID, paste0(info$loci$EXPORT_ID, ".2")))
  rownames(a) <- ids
  list(a = a, samples = samples, loci = info$loci, alleles = info$alleles,
    counts = info$counts, populations = if (grouped) meta$populations else NULL)
}
.legacy_dosage <- function(x) {
  odd <- seq.int(1L, ncol(x$a), 2L)
  d <- x$a[, odd, drop = FALSE] + x$a[, odd + 1L, drop = FALSE] - 2L
  colnames(d) <- x$loci$EXPORT_ID
  d
}
.legacy_publish <- function(x, filename, folder, suffix, overwrite, writer) {
  .export_flag(overwrite, "overwrite")
  label <- gsub("[^A-Za-z0-9]", "", suffix)
  paths <- .export_paths(filename, label, folder,
    c(suffix, "_samples.tsv", "_loci.tsv", "_alleles.tsv"), overwrite)
  staged <- vapply(paths, function(p) tempfile(), "")
  on.exit(unlink(staged), add = TRUE)
  writer(staged[1])
  readr::write_tsv(x$samples, staged[2])
  readr::write_tsv(x$loci, staged[3])
  readr::write_tsv(x$alleles, staged[4])
  .export_publish(staged, paths, overwrite)
  paths
}
.legacy_object <- function(object, x, filename, folder, suffix, overwrite) {
  attr(object, "loci") <- x$loci
  attr(object, "alleles") <- x$alleles
  attr(object, "samples") <- x$samples
  paths <- character()
  if (!is.null(filename)) paths <- .legacy_publish(x, filename, folder, suffix,
    overwrite, function(p) saveRDS(object, p))
  attr(object, "export.files") <- paths
  object
}
.legacy_start <- function(name, verbose) {
  .export_flag(verbose, "verbose")
  tgbase::startup(package = "genometranslator", f.name = name, verbose = verbose)
}
.legacy_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop(paste("Install the optional package", package, "before using this exporter."))
}
