# read_genlight ----------------------------------------------------------------
#' @name read_genlight
#' @title Read a genlight object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Reading and writing \code{genlight} objects requires the optional CRAN
#' package \href{https://cran.r-project.org/package=adegenet}{\pkg{adegenet}}.
#' @description Import an \pkg{adegenet} \code{genlight} object, standardize its
#' sample and marker metadata, and optionally produce tidy genomic data or a GDS
#' representation.

#' @param data (path or object) A genlight object in the global environment or
#' path to a genlight file that will be open with \code{readRDS}.


#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds Logical. Generate a genometranslator GDS representation.
#' Default: \code{gds = TRUE}.

#' @param write Logical. Save the tidy table as a uniquely named RDS file.
#' Default: \code{write = FALSE}.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_genlight
#' @return With \code{tidy = TRUE}, a tidy tibble. With
#' \code{tidy = FALSE, gds = TRUE}, the generated GDS result. With both options
#' false, the original \code{genlight} object.
#' tidy=TRUE takes precedence without creating a GDS as a side effect.
#' @param parallel.core Positive integer retained for reader-call compatibility;
#' this conversion is serial.
#' @examples
#' \dontrun{
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   x <- readRDS("genotypes.genlight.rds")
#'   genotypes <- genometranslator::read_genlight(x, gds = FALSE)
#' }
#' }

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @details
#' Missing marker fields are generated when necessary:
#' \enumerate{
#' \item \code{is.null(genlight@pop)}: pop will be integrated
#' in the tidy dataset.
#' \item \code{is.null(data@chromosome)}: \code{DENOVO} is used
#' in the tidy dataset.
#' \item \code{is.null(data@loc.names)}: L1 to \code{ncol(genlight)}
#' will be integrated in the tidy dataset.
#' \item \code{is.null(data@position)}: an integer string of
#' length = \code{ncol(genlight)} will be integrated in the tidy dataset.
#' }
#' Generated fields preserve a unique marker identifier but do not create
#' reference-genome coordinates. When \code{loc.all} is unavailable, symbolic
#' \code{REF = "0"} and \code{ALT = "1"} labels preserve the dosage orientation.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies


read_genlight <- function(data, tidy = TRUE, gds = TRUE, write = FALSE,
  verbose = FALSE, parallel.core = 1L) {
  start <- .legacy_start("read_genlight", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .legacy_package("adegenet")
  .export_count(parallel.core, "parallel.core", 1)
  if (is.character(data) && length(data)==1L) data <- readRDS(data)
  if (!inherits(data, "genlight")) stop("Input is not a genlight object.")
  if (any(adegenet::ploidy(data) != 2L)) stop("Only diploid genlight is supported.")
  d <- as.matrix(data)
  if (any(!is.na(d) & !d %in% 0:2)) stop("Invalid genlight dosage.")
  ids <- adegenet::indNames(data)
  if (is.null(ids)) ids <- paste0("IND-", seq_len(nrow(d)))
  markers <- adegenet::locNames(data)
  if (is.null(markers)) markers <- paste0("L", seq_len(ncol(d)))
  labels <- if (is.null(data@loc.all)) rep(list(c("0", "1")), ncol(d)) else
    strsplit(data@loc.all, "/", fixed = TRUE)
  if (length(labels) != ncol(d) || any(lengths(labels) != 2L))
    stop("genlight requires two allele labels per locus.")
  names(labels) <- markers
  a <- matrix(NA_integer_, nrow(d), 2*ncol(d))
  for (j in seq_len(ncol(d))) {
    a[, 2*j-1] <- ifelse(is.na(d[,j]), NA_integer_, ifelse(d[,j]==2L, 2L, 1L))
    a[, 2*j] <- ifelse(is.na(d[,j]), NA_integer_, ifelse(d[,j]==0L, 1L, 2L))
  }
  .object_read_result(data, a, labels, ids, adegenet::pop(data),
    tidy, gds, write, verbose, "genlight", data@chromosome, data@position)
}

#' Write genlight
#'
#' @description Create the target diploid object directly from GDS allele indices. Multiallelic loci are supported except by genlight, which requires biallelic loci.
#' @param data Open SeqArray GDS or GDS filename.
#' @param write Write the returned object as RDS.
#' @param verbose Display progress messages.
#' @param filename Output basename.
#' @param strata Optional metadata table or TSV with INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Replace existing output files.
#' @return The target object with export.files and locus mapping attributes.
#' @details Uses active GDS whitelists and restores selections on exit.
#' No implicit filtering or imputation. Partial missing calls are rejected.
#' In-memory objects must fit in RAM. Numeric alleles are labels, not repeat sizes.
#' @examples
#' \dontrun{
#' write_genlight("study.gds", strata = "samples.tsv")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_genlight <- function(data, write = FALSE, verbose = FALSE,
  filename = NULL, strata = NULL, path.folder = getwd(), chunk.size = 32L,
  overwrite = FALSE) {
  start <- .legacy_start("write_genlight", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(write, "write"); .legacy_package("adegenet")
  x <- .legacy_snapshot(data, strata, TRUE, TRUE, chunk.size)
  out <- methods::new("genlight", .legacy_dosage(x), ploidy = 2L,
    ind.names = x$samples$INDIVIDUALS, loc.names = x$loci$EXPORT_ID,
    pop = factor(x$samples$GROUP))
  out@chromosome <- factor(x$loci$CHROM)
  out@position <- as.integer(x$loci$POS)
  out@loc.all <- vapply(seq_len(nrow(x$loci)), function(j)
    paste(x$alleles$ALLELE[x$alleles$EXPORT_ID == x$loci$EXPORT_ID[j]],
      collapse = "/"), "")
  .legacy_object(out, x, if (write) if (is.null(filename)) "genlight" else filename
    else NULL, path.folder, "_genlight.rds", overwrite)
}
