#' Write BayeScan
#'
#' @description Write codominant population allele counts directly from GDS.
#' @param data Open SeqArray GDS or GDS filename.
#' @param pop.select Optional populations to retain.
#' @param filename Output basename.
#' @param strata Metadata table or TSV containing INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Samples per read block.
#' @param overwrite Replace existing files.
#' @param verbose Display progress messages.
#' @return A list of files, pop.dictionary and markers.dictionary.
#' @details Requires polymorphic loci with called genotypes in every selected
#' population. Violations are errors, never silent filtering. Multiallelic
#' allele counts retain a common dictionary across populations; missing calls
#' contribute no genes. Partial calls are rejected. GDS filters are restored.
#' @references Foll M, Gaggiotti O (2008). A genome-scan method to identify
#' selected loci appropriate for both dominant and codominant markers:
#' a Bayesian perspective. Genetics 180:977-993.
#' @seealso \href{https://github.com/mfoll/BayeScan}{BayeScan format specification}
#' @examples
#' \dontrun{
#' write_bayescan("study.gds", strata = "samples.tsv", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_bayescan <- function(data, pop.select = NULL, filename = NULL,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE,
  verbose = TRUE) {
  start <- .legacy_start("write_bayescan", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(overwrite, "overwrite")
  x <- .legacy_snapshot(data, strata, TRUE, chunk.size = chunk.size)
  if (!is.null(pop.select)) {
    if (!length(pop.select) || anyNA(pop.select) ||
        !all(pop.select %in% x$samples$GROUP)) stop("Invalid pop.select.")
    keep <- x$samples$GROUP %in% pop.select
    x$samples <- x$samples[keep, ]; x$a <- x$a[keep, , drop = FALSE]
  }
  groups <- unique(x$samples$GROUP)
  if (length(groups) < 2L) stop("BayeScan requires at least two populations.")
  counts <- lapply(seq_len(nrow(x$loci)), function(j) {
    z <- vapply(groups, function(p)
      tabulate(as.vector(x$a[x$samples$GROUP == p, (2*j-1):(2*j), drop = FALSE]),
        nbins = x$counts[j]), integer(x$counts[j]))
    if (any(colSums(z) == 0L)) stop("Each locus needs called genotypes in every population.")
    if (sum(rowSums(z)>0) < 2L) stop("BayeScan requires polymorphic loci.")
    z
  })
  pop.dictionary <- tibble::tibble(STRATA = groups, BAYESCAN_POP = seq_along(groups))
  markers.dictionary <- tibble::tibble(MARKERS = x$loci$MARKERS,
    BAYESCAN_MARKERS = seq_len(nrow(x$loci)))
  paths <- .export_paths(filename, "bayescan", path.folder,
    c("_bayescan.txt", "_bayescan_pop_dictionary.tsv",
      "_bayescan_markers_dictionary.tsv", "_bayescan_alleles.tsv"), overwrite)
  staged <- vapply(paths, function(p) tempfile(), "")
  on.exit(unlink(staged), add = TRUE)
  con <- file(staged[1], "wt")
  tryCatch({
    writeLines(c(paste0("[loci]=", length(counts)), "",
      paste0("[populations]=", length(groups)), ""), con)
    for (p in seq_along(groups)) {
      writeLines(paste0("[pop]=", p), con)
      for (j in seq_along(counts))
        writeLines(paste(c(j, sum(counts[[j]][,p]), nrow(counts[[j]]),
          counts[[j]][,p]), collapse = " "), con)
      writeLines("", con)
    }
  }, finally = close(con))
  readr::write_tsv(pop.dictionary, staged[2])
  readr::write_tsv(markers.dictionary, staged[3])
  readr::write_tsv(x$alleles, staged[4])
  .export_publish(staged, paths, overwrite)
  list(files = paths, pop.dictionary = pop.dictionary,
    markers.dictionary = markers.dictionary)
}
