# read_genind ------------------------------------------------------------------
#' @name read_genind
#' @title Read a genind object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Reading and writing \code{genind} objects requires the optional CRAN package
#' \href{https://cran.r-project.org/package=adegenet}{\pkg{adegenet}}.
#' @description Import an \pkg{adegenet} \code{genind} object, standardize its
#' sample and locus information, and optionally produce tidy genomic data or a
#' GDS representation. Diploid biallelic and multiallelic objects are supported;
#' GDS generation is currently limited to biallelic data.

#' @param data (path or object) A genind object in the global environment or
#' path to a genind file that will be open with \code{readRDS}.

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds Logical. Generate a genometranslator GDS representation.
#' Currently, for biallelic datasets only.
#' Default: \code{gds = TRUE}.

#' @param write Logical. Save the tidy table as a uniquely named RDS file.
#' Default: \code{write = FALSE}.

#' @note \href{https://github.com/thibautjombart/adegenet}{genind} objects, like
#' genepop, are not optimal genomic format for RADseq datasets,
#' they lack important genotypes and markers metadata: chromosome, locus, snp,
#' position, read depth, allele depth, etc.
#' \href{https://github.com/thibautjombart/adegenet}{genlight} object is a more
#' interesting container and is memory efficient, see \code{\link{read_genlight}}.
#'
#'
#' Alleles are assigned stable codes in the object's allele order. The
#' allele_dictionary attribute retains the original allele labels. Locus names,
#' including periods, are preserved rather than parsed as allele separators.
#'
#'
#' In the unlikely event that the genind object as no stratification/population,
#' \emph{pop} will be added to the strata column.


#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_genind
#' @return With \code{tidy = TRUE}, a tidy tibble. With
#' \code{tidy = FALSE, gds = TRUE}, the generated GDS result. With both options
#' false, the original \code{genind} object. tidy=TRUE takes precedence and does
#' not create a GDS as a side effect. Requesting GDS for multiallelic input is an
#' explicit error; the tidy representation supports multiallelic loci.
#' @examples
#' \dontrun{
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   data("nancycats", package = "adegenet")
#'   cats <- genometranslator::read_genind(nancycats, gds = FALSE)
#' }
#' }
#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies


read_genind <- function(data, tidy = TRUE, gds = TRUE, write = FALSE,
  verbose = FALSE) {
  start <- .legacy_start("read_genind", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .legacy_package("adegenet")
  if (is.character(data) && length(data)==1L) data <- readRDS(data)
  if (!inherits(data, "genind")) stop("Input is not a genind object.")
  if (data@type != "codom" || any(data@ploidy != 2L))
    stop("Only diploid codominant genind objects are supported.")
  labels <- data@all.names
  ids <- adegenet::indNames(data)
  a <- matrix(NA_integer_, nrow(data@tab), 2L * length(labels))
  for (j in seq_along(labels)) {
    cols <- which(as.character(data@loc.fac) == names(labels)[j])
    counts <- data@tab[, cols, drop = FALSE]
    if (ncol(counts) != length(labels[[j]])) stop("Inconsistent genind allele metadata.")
    for (i in seq_len(nrow(counts))) {
      z <- counts[i, ]
      if (all(is.na(z))) next
      if (anyNA(z) || any(z < 0 | z != floor(z)) || sum(z) != 2L)
        stop("Each called genind locus must contain exactly two allele copies.")
      a[i, (2*j-1):(2*j)] <- rep(seq_along(z), times = z)
    }
  }
  .object_read_result(data, a, labels, ids, adegenet::pop(data),
    tidy, gds, write, verbose, "genind")
}


#' Write genind
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
#' write_genind("study.gds", strata = "samples.tsv")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_genind <- function(data, write = FALSE, verbose = FALSE,
  filename = NULL, strata = NULL, path.folder = getwd(), chunk.size = 32L,
  overwrite = FALSE) {
  start <- .legacy_start("write_genind", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(write, "write"); .legacy_package("adegenet")
  x <- .legacy_snapshot(data, strata, TRUE, chunk.size = chunk.size)
  odd <- seq.int(1L, ncol(x$a), 2L)
  calls <- matrix(paste(x$a[, odd], x$a[, odd+1L], sep = "/"), nrow(x$a))
  calls[is.na(x$a[, odd])] <- NA_character_
  colnames(calls) <- x$loci$EXPORT_ID
  out <- adegenet::df2genind(as.data.frame(calls), sep = "/", ploidy = 2L,
    ind.names = x$samples$INDIVIDUALS, pop = factor(x$samples$GROUP),
    type = "codom")
  .legacy_object(out, x, if (write) if (is.null(filename)) "genind" else filename
    else NULL, path.folder, "_genind.rds", overwrite)
}
