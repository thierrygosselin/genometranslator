#' Write SNPRelate
#'
#' @description Converts the active SeqArray selection directly using SeqArray::seqGDS2SNP. Requires SNPRelate and biallelic loci. Returns an open read-only SNPRelate GDS; close it with SNPRelate::snpgdsClose.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param filename Optional output basename; use path.folder for the directory.
#' @param verbose Display progress messages.
#' @param path.folder Output directory.
#' @param overwrite Allow replacement of existing output files.
#' @details GDS selections are restored, including on errors. Partial diploid
#' calls are rejected by allele-based conversions. No implicit biological
#' filtering or imputation is performed. Safe locus/sample IDs have mapping
#' tables. In-memory matrices must fit in RAM; block reads limit temporary memory.
#' @return Export paths and mappings, or the target object as described above.
#' @references Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS.
#' A high-performance computing toolset for relatedness and principal component
#' analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
#' doi:10.1093/bioinformatics/bts606
#' @examples
#' \dontrun{
#' write_snprelate("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_snprelate <- function(data, filename = NULL, verbose = TRUE,
  path.folder = getwd(), overwrite = FALSE) {
  start <- .legacy_start("write_snprelate", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .legacy_package("SNPRelate")
  context <- .export_gds(data); on.exit(.export_close(context), add = TRUE)
  info <- .gsi_loci(context$gds)
  if (any(info$counts != 2L)) stop("SNPRelate export requires biallelic loci.")
  ids <- as.character(SeqArray::seqGetData(context$gds, "sample.id"))
  x <- list(samples = tibble::tibble(INDIVIDUALS = ids),
    loci = info$loci, alleles = info$alleles)
  paths <- .legacy_publish(x, filename, path.folder, "_snprelate.gds", overwrite,
    function(p) SeqArray::seqGDS2SNP(context$gds, p, verbose = verbose))
  SNPRelate::snpgdsOpen(paths[1], readonly = TRUE)
}
