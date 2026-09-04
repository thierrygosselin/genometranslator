#' Write related
#'
#' @description Headerless text: sample ID followed by adjacent allele pairs, missing alleles encoded as zero. Positive allele labels are locus-specific codes, not repeat lengths.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param filename Optional output basename; use path.folder for the directory.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Allow replacement of existing output files.
#' @param verbose Display progress messages.
#' @details GDS selections are restored, including on errors. Partial diploid
#' calls are rejected by allele-based conversions. No implicit biological
#' filtering or imputation is performed. Safe locus/sample IDs have mapping
#' tables. In-memory matrices must fit in RAM; block reads limit temporary memory.
#' @return Export paths and mappings, or the target object as described above.
#' @references Pew J, Muir PH, Wang J, Frasier TR (2015)
#' related: an R package for analysing pairwise relatedness from codominant
#' molecular markers.
#' Molecular Ecology Resources, 15, 557-561.
#' @references Wang, J. 2011.
#' COANCESTRY: A program for simulating, estimating and analysing relatedness
#' and inbreeding coefficients.
#' Molecular Ecology Resources 11(1): 141-145.
#' @examples
#' \dontrun{
#' write_related("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_related <- function(data, filename = NULL, path.folder = getwd(),
  chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  start <- .legacy_start("write_related", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, chunk.size = chunk.size)
  a <- x$a; a[is.na(a)] <- 0L
  paths <- .legacy_publish(x, filename, path.folder, "_related.txt", overwrite,
    function(p) utils::write.table(cbind(x$samples$EXPORT_ID, a), p,
      quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " "))
  invisible(list(files = paths, samples = x$samples, loci = x$loci))
}
