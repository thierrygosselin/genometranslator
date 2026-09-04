#' Write fineRADstructure
#'
#' @description Write RADpainter input from locus-level nucleotide haplotypes.
#' @param data Open SeqArray GDS or GDS path.
#' @param strata Reserved sample metadata; not needed for RADpainter input.
#' @param filename Output basename.
#' @param path.folder Output directory.
#' @param chunk.size Samples per read block.
#' @param overwrite Replace existing files.
#' @param verbose Display progress messages.
#' @return A list containing files, samples and loci.
#' @details Each GDS variant must already represent one complete locus haplotype
#' with nucleotide sequences as alleles. Multiple SNP records sharing a LOCUS
#' are rejected: concatenating unphased SNPs would invent haplotypes.
#' Equal-length A/C/G/T/N sequences are required within each locus.
#' Two alleles are separated by a slash; missing calls are empty fields.
#' Sample IDs are safe S identifiers with an original-ID mapping.
#' No biological filtering is performed. GDS selections are restored.
#' @references Malinsky M, Trucchi E, Lawson DJ, Falush D (2018).
#' RADpainter and fineRADstructure: Population inference from RADseq data.
#' Molecular Biology and Evolution 35:1284-1290. \doi{10.1093/molbev/msy023}.
#' @seealso \href{https://github.com/millanek/fineRADstructure}{fineRADstructure}
#' @examples
#' \dontrun{
#' write_fineradstructure("haplotypes.gds", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_fineradstructure <- function(data, strata = NULL, filename = NULL,
  path.folder = getwd(), chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  start <- .legacy_start("write_fineradstructure", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  if (!is.null(strata)) stop("Set sample selection in GDS; strata is not used by RADpainter.")
  context <- .export_gds(data); on.exit(.export_close(context), add = TRUE)
  mm <- genometranslator::extract_markers_metadata(context$gds)
  mm <- mm[mm$VARIANT_ID %in% SeqArray::seqGetData(context$gds, "variant.id"), ]
  if ("LOCUS" %in% names(mm) && anyDuplicated(mm$LOCUS))
    stop("Multiple SNPs per LOCUS: supply locus haplotypes; unphased concatenation is unsafe.")
  x <- .legacy_snapshot(context$gds, chunk.size = chunk.size)
  if (any(!grepl("^[ACGTN]+$", x$alleles$ALLELE)))
    stop("RADpainter requires nucleotide haplotype alleles.")
  for (locus in x$loci$EXPORT_ID)
    if (length(unique(nchar(x$alleles$ALLELE[x$alleles$EXPORT_ID==locus]))) != 1L)
      stop("Haplotype sequences must have equal length within each locus.")
  paths <- .legacy_publish(x, filename, path.folder, "_fineradstructure.tsv", overwrite,
    function(p) {
      con <- file(p, "wt"); on.exit(close(con))
      writeLines(paste(x$samples$EXPORT_ID, collapse = "\t"), con)
      for (j in seq_len(nrow(x$loci))) {
        dict <- x$alleles$ALLELE[x$alleles$EXPORT_ID==x$loci$EXPORT_ID[j]]
        a <- x$a[,2*j-1]; b <- x$a[,2*j]
        z <- ifelse(is.na(a) | is.na(b), "", paste(dict[a],dict[b],sep="/"))
        writeLines(paste(z,collapse="\t"), con)
      }
    })
  invisible(list(files=paths, samples=x$samples, loci=x$loci))
}
