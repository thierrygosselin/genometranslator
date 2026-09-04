#' Write fastSTRUCTURE
#'
#' @description Text has no header and two rows per sample, six metadata columns followed by one allele per locus. Missing alleles are -9. Only biallelic loci are accepted. plink.bed=TRUE delegates to write_plink(format='bed').
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param plink.bed Export PLINK BED instead of text.
#' @param pop.levels Population labels in the desired order.
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
#' @references Raj A, Stephens M, Pritchard JK (2014)
#' fastSTRUCTURE: Variational Inference of Population Structure in Large SNP
#' Datasets. Genetics, 197, 573-589.
#' @examples
#' \dontrun{
#' write_faststructure("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_faststructure <- function(data, plink.bed = FALSE, pop.levels = NULL,
  filename = NULL, strata = NULL, path.folder = getwd(), chunk.size = 32L,
  overwrite = FALSE, verbose = TRUE) {
  start <- .legacy_start("write_faststructure", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(plink.bed, "plink.bed")
  if (plink.bed) {
    if (!is.null(pop.levels) || !is.null(strata))
      stop("For BED output, set sample metadata/order in GDS before export.")
    return(genometranslator::write_plink(data, filename = filename,
      format = "bed", path.folder = path.folder, overwrite = overwrite,
      verbose = verbose))
  }
  x <- .legacy_snapshot(data, strata, TRUE, TRUE, chunk.size)
  groups <- x$populations$GROUP
  if (!is.null(pop.levels)) {
    if (anyNA(pop.levels) || anyDuplicated(pop.levels) || !setequal(groups, pop.levels))
      stop("pop.levels must list every population exactly once.")
    groups <- pop.levels
  }
  a <- x$a; a[is.na(a)] <- -9L
  paths <- .legacy_publish(x, filename, path.folder, ".str", overwrite, function(p) {
    con <- file(p, "wt"); on.exit(close(con))
    for (i in seq_len(nrow(a))) for (copy in 1:2)
      writeLines(paste(c(x$samples$EXPORT_ID[i],
        match(x$samples$GROUP[i], groups), 0, 0, 0, 0,
        a[i, seq.int(copy, ncol(a), 2L)]), collapse = "\t"), con)
  })
  invisible(list(files = paths, samples = x$samples, loci = x$loci))
}
