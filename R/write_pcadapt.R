#' Write pcadapt
#'
#' @description Rows are loci, columns are individuals; entries count ALT copies (0/1/2), with missing calls encoded as 9. Use pcadapt::read.pcadapt(..., type='pcadapt').
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param pop.select Optional populations to retain explicitly.
#' @param filename Optional output basename; use path.folder for the directory.
#' @param strata Metadata table or TSV with INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Allow replacement of existing output files.
#' @param verbose Display progress messages.
#' @details GDS selections are restored, including on errors. Partial diploid
#' calls are rejected by allele-based conversions. No implicit biological
#' filtering or imputation is performed. Safe locus/sample IDs have mapping
#' tables. In-memory matrices must fit in RAM; block reads limit temporary memory.
#' @return Export paths and mappings, or the target object as described above.
#' @references Luu, K., Bazin, E., & Blum, M. G. (2017).
#' pcadapt: an R package to perform genome scans for selection based on principal component analysis.
#' Molecular Ecology Resources, 17(1), 67-77.
#' \doi{10.1111/1755-0998.12592}
#' @examples
#' \dontrun{
#' write_pcadapt("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_pcadapt <- function(data, pop.select = NULL, filename = NULL,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE,
  verbose = TRUE) {
  start <- .legacy_start("write_pcadapt", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, strata, grouped = !is.null(pop.select),
    biallelic = TRUE, chunk.size = chunk.size)
  if (!is.null(pop.select)) {
    if (!length(pop.select) || anyNA(pop.select) ||
        !all(pop.select %in% x$samples$GROUP)) stop("Unknown or empty pop.select.")
    keep <- x$samples$GROUP %in% pop.select
    x$samples <- x$samples[keep, ]; x$a <- x$a[keep, , drop = FALSE]
  }
  d <- t(.legacy_dosage(x)); d[is.na(d)] <- 9L
  paths <- .legacy_publish(x, filename, path.folder, "_pcadapt.txt", overwrite,
    function(p) utils::write.table(d, p, row.names = FALSE,
      col.names = FALSE, quote = FALSE))
  list(files = paths, genotype.matrix = d,
    pop.string = if ("GROUP" %in% names(x$samples)) x$samples$GROUP else NULL,
    samples = x$samples, loci = x$loci)
}
