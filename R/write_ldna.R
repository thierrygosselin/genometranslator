#' Write LDna
#'
#' @description Computes SNPRelate gametic r-squared and returns the lower triangle (NA diagonal and upper triangle). Memory is quadratic in marker count. max.markers prevents unexpectedly large allocation. Requires SNPRelate.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param filename Optional output basename; use path.folder for the directory.
#' @param parallel.core Number of LD calculation threads.
#' @param path.folder Output directory.
#' @param overwrite Allow replacement of existing output files.
#' @param verbose Display progress messages.
#' @param max.markers Maximum number of loci allowed in the dense LD matrix.
#' @details GDS selections are restored, including on errors. Partial diploid
#' calls are rejected by allele-based conversions. No implicit biological
#' filtering or imputation is performed. Safe locus/sample IDs have mapping
#' tables. In-memory matrices must fit in RAM; block reads limit temporary memory.
#' @return Export paths and mappings, or the target object as described above.
#' @references Kemppainen P, Knight CG, Sarma DK et al. (2015)
#' Linkage disequilibrium network analysis (LDna) gives a global view of
#' chromosomal inversions, local adaptation and geographic structure.
#' Molecular Ecology Resources, 15, 1031-1045.
#' @examples
#' \dontrun{
#' write_ldna("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_ldna <- function(data, filename = NULL, parallel.core = 1L,
  path.folder = getwd(), overwrite = FALSE, verbose = TRUE,
  max.markers = 10000L) {
  start <- .legacy_start("write_ldna", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_count(parallel.core, "parallel.core", 1)
  .export_count(max.markers, "max.markers", 2)
  context <- .export_gds(data); on.exit(.export_close(context), add = TRUE)
  info <- .gsi_loci(context$gds)
  if (nrow(info$loci) > max.markers)
    stop("LD matrix exceeds max.markers; select fewer loci to limit memory.")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  g <- genometranslator::write_snprelate(context$gds, filename = "ld",
    path.folder = folder, verbose = FALSE)
  on.exit(SNPRelate::snpgdsClose(g), add = TRUE, after = FALSE)
  result <- SNPRelate::snpgdsLDMat(g, slide = -1, mat.trim = FALSE,
    method = "r", num.thread = parallel.core, verbose = verbose)
  d <- result$LD^2
  d[upper.tri(d, diag = TRUE)] <- NA_real_
  dimnames(d) <- list(info$loci$EXPORT_ID, info$loci$EXPORT_ID)
  x <- list(samples = tibble::tibble(INDIVIDUALS =
    as.character(SeqArray::seqGetData(context$gds, "sample.id"))),
    loci = info$loci, alleles = info$alleles)
  .legacy_object(d, x, filename, path.folder, "_ldna.rds", overwrite)
}
