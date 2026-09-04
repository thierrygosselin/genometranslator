#' Write stockR
#'
#' @description Returns a loci-by-individuals ALT dosage matrix (0/1/2, NA missing), with grps and sample.grps attributes. Optional file output uses RDS.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param filename Optional output basename; use path.folder for the directory.
#' @param verbose Display progress messages.
#' @param strata Metadata table or TSV with INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Allow replacement of existing output files.
#' @details GDS selections are restored, including on errors. Partial diploid
#' calls are rejected by allele-based conversions. No implicit biological
#' filtering or imputation is performed. Safe locus/sample IDs have mapping
#' tables. In-memory matrices must fit in RAM; block reads limit temporary memory.
#' @return Export paths and mappings, or the target object as described above.
#' @references Foster SD, Feutry P, Grewe PM, Berry O, Hui FKC, Davies CR.
#' Reliably discriminating stock structure with genetic markers: Mixture models
#' with robust and fast computation. \doi{10.1111/1755-0998.12920}.
#' @examples
#' \dontrun{
#' write_stockr("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_stockr <- function(data, filename = NULL, verbose = TRUE,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE) {
  start <- .legacy_start("write_stockr", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, strata, TRUE, TRUE, chunk.size)
  d <- t(.legacy_dosage(x))
  attr(d, "grps") <- x$samples$GROUP
  attr(d, "sample.grps") <- factor(x$samples$INDIVIDUALS)
  .legacy_object(d, x, filename, path.folder, "_stockr.rds", overwrite)
}
