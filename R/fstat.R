# read_fstat -------------------------------------------------------------------
#' @name read_fstat
#' @title Read an FSTAT file

#' @description Import and validate diploid data in the FSTAT text format
#' (Goudet, 1995). The header, locus count, unique locus names, and genotype-row
#' widths are checked before genotypes are standardized and returned in long or
#' wide form.

#' @param data Path to an FSTAT file, commonly with extension \code{.dat}, or a
#' one-column data frame containing its lines.

#' @param strata (optional) A tab delimited file with 2 columns. Header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param tidy (optional, logical) With \code{tidy = FALSE},
#' the markers are the variables and the genotypes the observations (wide format).
#' With the default: \code{tidy = TRUE}, markers and genotypes are variables
#' with their own columns (long format).

#' Default: \code{tidy = TRUE}.
#' @param filename (optional) The file name for the tidy data frame
#' written to the working directory.
#' With the default, The tidy data is
#' in the global environment only (i.e. not written in the working directory).

#' Default: \code{filename = NULL}.
#' @return A tibble in long form with \code{STRATA}, \code{INDIVIDUALS},
#' \code{MARKERS}, and six-digit \code{GENOTYPE}, or wide form when
#' \code{tidy = FALSE}. If \code{filename} is supplied, the same table is also
#' written as a tab-separated file.
#'
#' @section Field handling:
#' Population codes become \code{STRATA}; deterministic identifiers of the form
#' \code{IND-1}, \code{IND-2}, ... become \code{INDIVIDUALS}; locus names become
#' \code{MARKERS}; and each allele is padded to three digits. Supplied strata
#' metadata replaces the embedded population code and must contain every
#' generated individual identifier.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_fstat

#' @examples
#' \dontrun{
#' # We will use the fstat dataset provided with adegenet package
#' if (requireNamespace("hierfstat", quietly = TRUE)) {
#'
#' # The simplest form of the function:
#' fstat.file <- genometranslator::read_fstat(
#'     data = system.file(
#'     "extdata/diploid.dat",
#'     package = "hierfstat"
#'     )
#'  )
#'
#' # To output a data frame in wide format, with markers in separate columns:
#' nancycats.wide <- genometranslator::read_fstat(
#'     data = system.file(
#'         "extdata/diploid.dat",
#'         package = "hierfstat"
#'     ),
#' tidy = FALSE
#' )
#' }
#' }


#' @references Goudet J. (1995).
#' FSTAT (Version 1.2): A computer program to calculate F-statistics.
#' Journal of Heredity, 86, 485-486.

#' @seealso \href{https://github.com/jgx65/hierfstat}{hierfstat}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies


read_fstat <- function(data, strata = NULL, tidy = TRUE, filename = NULL,
  verbose = FALSE) {
  start <- .legacy_start("read_fstat", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(tidy, "tidy")
  lines <- if (is.character(data) && length(data) == 1L && file.exists(data))
    readLines(data, warn = FALSE) else if (is.data.frame(data) && ncol(data) == 1L)
      as.character(data[[1]]) else stop("Supply an FSTAT path or one-column table.")
  lines <- trimws(lines)
  if (anyNA(lines) || any(!nzchar(lines))) stop("Empty FSTAT records are not allowed.")
  h <- strsplit(lines[1], "\\s+")[[1]]
  if (length(h) != 4L || any(!grepl("^[0-9]+$", h)))
    stop("FSTAT header requires four positive integers.")
  h <- as.integer(h)
  if (anyNA(h) || any(h < 1L) || !h[4] %in% 1:3 || h[3] >= 10^h[4])
    stop("Invalid FSTAT population/locus count, maximum allele or coding width.")
  if (length(lines) <= h[2] + 1L) stop("FSTAT contains no genotype rows.")
  markers <- lines[seq_len(h[2]) + 1L]
  if (anyDuplicated(markers) || any(markers %in% c("STRATA", "INDIVIDUALS")))
    stop("FSTAT locus names must be unique and not reserved.")
  rows <- strsplit(lines[-seq_len(h[2] + 1L)], "\\s+")
  if (any(lengths(rows) != h[2] + 1L))
    stop("FSTAT rows do not contain one population field and nl genotypes.")
  m <- do.call(rbind, rows)
  if (any(!grepl("^[0-9]+$", m[, 1]))) stop("Invalid FSTAT population code.")
  pop <- as.integer(m[, 1])
  if (anyNA(pop) || any(pop < 1L | pop > h[1])) stop("Population outside FSTAT header range.")
  gt <- m[, -1, drop = FALSE]
  if (any(!grepl("^[0-9]+$", gt)) || any(nchar(gt) != 2L*h[4]))
    stop("Genotype width does not match the FSTAT header.")
  a <- as.integer(substr(gt, 1, h[4])); b <- as.integer(substr(gt, h[4]+1, 2*h[4]))
  if (any(a > h[3] | b > h[3])) stop("Allele exceeds the FSTAT header maximum.")
  gt[] <- paste0(sprintf("%03d", a), sprintf("%03d", b))
  ids <- paste0("IND-", seq_len(nrow(gt)))
  groups <- as.character(pop)
  if (!is.null(strata)) {
    meta <- genometranslator::read_strata(strata, verbose = FALSE)$strata
    groups <- as.character(meta$STRATA[match(ids, meta$INDIVIDUALS)])
    if (anyNA(groups)) stop("Some FSTAT individuals are absent from strata.")
  }
  out <- if (tidy) tibble::tibble(STRATA = rep(groups, h[2]),
    INDIVIDUALS = rep(ids, h[2]), MARKERS = rep(markers, each = length(ids)),
    GENOTYPE = as.vector(gt)) else dplyr::bind_cols(
      tibble::tibble(STRATA = groups, INDIVIDUALS = ids),
      tibble::as_tibble(gt, .name_repair = function(x) markers))
  if (!is.null(filename)) readr::write_tsv(out, filename)
  out
}

#' Write FSTAT
#'
#' @description Write an FSTAT text file and return a hierfstat-style data frame. Allele width follows the dictionary (one to three digits), rather than assuming single-digit alleles.
#' @param data Open SeqArray GDS or GDS filename.
#' @param filename Output basename.
#' @param strata Optional metadata table or TSV with INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Replace existing output files.
#' @param verbose Display progress messages.
#' @return The target object with export.files and locus mapping attributes.
#' @details Uses active GDS whitelists and restores selections on exit.
#' No implicit filtering or imputation. Partial missing calls are rejected.
#' In-memory objects must fit in RAM. Numeric alleles are labels, not repeat sizes.
#' @examples
#' \dontrun{
#' write_hierfstat("study.gds", strata = "samples.tsv")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_hierfstat <- function(data, filename = NULL, strata = NULL,
  path.folder = getwd(), chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  start <- .legacy_start("write_hierfstat", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, strata, TRUE, chunk.size = chunk.size)
  maximum <- max(x$counts)
  if (maximum > 999L) stop("FSTAT supports at most three digits per allele.")
  width <- nchar(as.character(maximum))
  a <- x$a; a[is.na(a)] <- 0L
  odd <- seq.int(1L, ncol(a), 2L)
  calls <- matrix(paste0(sprintf(paste0("%0", width, "d"), a[, odd]),
    sprintf(paste0("%0", width, "d"), a[, odd+1L])), nrow(a))
  pop <- match(x$samples$GROUP, x$populations$GROUP)
  x$samples$FSTAT_POP <- pop
  paths <- .legacy_publish(x, filename, path.folder, "_fstat.dat", overwrite,
    function(p) {
      con <- file(p, "wt"); on.exit(close(con))
      writeLines(c(paste(nrow(x$populations), nrow(x$loci), maximum, width),
        x$loci$EXPORT_ID), con)
      for (i in seq_len(nrow(a)))
        writeLines(paste(c(pop[i], calls[i, ]), collapse = " "), con)
    })
  out <- data.frame(pop = pop, matrix(as.integer(calls), nrow(a)))
  names(out)[-1] <- x$loci$EXPORT_ID
  out[-1] <- lapply(out[-1], function(z) {z[z == 0L] <- NA_integer_; z})
  attr(out, "export.files") <- paths
  attr(out, "loci") <- x$loci
  out
}
