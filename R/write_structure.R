#' Write STRUCTURE
#'
#' @description Two rows per diploid sample, a locus header, sample label and numeric population. Set ONEROWPERIND=0, LABEL=1, POPDATA=1, MARKERNAMES=1 and MISSING=-9 in STRUCTURE. Allele codes are positive locus-specific labels.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
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
#' @references Pritchard JK, Stephens M, Donnelly P. (2000)
#' Inference of population structure using multilocus genotype data.
#' Genetics. Genetics Society of America. 155: 945–959.
#' @examples
#' \dontrun{
#' write_structure("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_structure <- function(data, pop.levels = NULL, filename = NULL,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE,
  verbose = TRUE) {
  start <- .legacy_start("write_structure", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, strata, TRUE, chunk.size = chunk.size)
  groups <- x$populations$GROUP
  if (!is.null(pop.levels)) {
    if (anyNA(pop.levels) || anyDuplicated(pop.levels) || !setequal(groups, pop.levels))
      stop("pop.levels must list every population exactly once.")
    groups <- pop.levels
  }
  x$samples$POP <- match(x$samples$GROUP, groups)
  a <- x$a; a[is.na(a)] <- -9L
  paths <- .legacy_publish(x, filename, path.folder, ".str", overwrite, function(p) {
    con <- file(p, "wt"); on.exit(close(con))
    writeLines(paste(x$loci$EXPORT_ID, collapse = "\t"), con)
    for (i in seq_len(nrow(a))) for (copy in 1:2)
      writeLines(paste(c(x$samples$EXPORT_ID[i], x$samples$POP[i],
        a[i, seq.int(copy, ncol(a), 2L)]), collapse = "\t"), con)
  })
  invisible(list(files = paths, samples = x$samples, loci = x$loci))
}
