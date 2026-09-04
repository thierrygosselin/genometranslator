#' Write HapMap
#'
#' @description Writes the eleven HapMap metadata columns followed by diploid nucleotide calls (NN missing). Biallelic A/C/G/T SNPs with mapped coordinates are required; no coordinates are invented.
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

#' @examples
#' \dontrun{
#' write_hapmap("study.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_hapmap <- function(data, filename = NULL, path.folder = getwd(),
  chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  start <- .legacy_start("write_hapmap", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, biallelic = TRUE, chunk.size = chunk.size)
  if (any(!x$alleles$ALLELE %in% c("A", "C", "G", "T")))
    stop("HapMap export requires biallelic A/C/G/T SNPs.")
  if (anyNA(x$loci$POS) || any(x$loci$POS < 1) ||
      anyNA(x$loci$CHROM) || any(x$loci$CHROM %in% c("", "DENOVO")))
    stop("HapMap requires mapped chromosome and position metadata.")
  paths <- .legacy_publish(x, filename, path.folder, ".hapmap.tsv", overwrite,
    function(p) {
      con <- file(p, "wt"); on.exit(close(con))
      writeLines(paste(c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
        "center", "protLSID", "assayLSID", "panelLSID", "QCcode",
        x$samples$EXPORT_ID), collapse = "\t"), con)
      for (j in seq_len(nrow(x$loci))) {
        dict <- x$alleles$ALLELE[x$alleles$EXPORT_ID == x$loci$EXPORT_ID[j]]
        a <- x$a[, 2*j-1]; b <- x$a[, 2*j]
        calls <- ifelse(is.na(a) | is.na(b), "NN", paste0(dict[a], dict[b]))
        writeLines(paste(c(x$loci$EXPORT_ID[j], paste(dict, collapse = "/"),
          x$loci$CHROM[j], x$loci$POS[j], "+", "NA", "genometranslator",
          rep("NA", 4), calls), collapse = "\t"), con)
      }
    })
  invisible(list(files = paths, samples = x$samples, loci = x$loci))
}
