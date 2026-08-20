#' Read a marker whitelist
#'
#' Import and standardise a marker whitelist supplied as an object or a
#' tab-separated file.
#'
#' @param whitelist.markers A data frame or path to a tab-separated whitelist.
#'   It may contain one or more of \code{MARKERS}, \code{CHROM}, \code{LOCUS},
#'   \code{POS}, \code{VARIANT_ID}, and \code{M_SEQ}.
#'   Default: \code{whitelist.markers = NULL}.
#' @param verbose Logical. Display progress messages.
#'   Default: \code{verbose = FALSE}.
#'
#' @return A distinct whitelist with character fields cleaned for genomic
#'   matching, or \code{NULL} when no whitelist is supplied.
#' @seealso \code{\link{clean_markers_names}}
#' @export
#' @template io-dependencies
read_whitelist <- function(whitelist.markers = NULL, verbose = FALSE) {
  if (is.null(whitelist.markers)) return(NULL)

  if (verbose) message("Reading whitelist of markers")

  if (is.character(whitelist.markers) && length(whitelist.markers) == 1L) {
    if (!file.exists(whitelist.markers)) {
      rlang::abort(paste0("whitelist file does not exist: ", whitelist.markers))
    }
    whitelist.markers <- suppressMessages(
      readr::read_tsv(whitelist.markers, col_names = TRUE, show_col_types = FALSE)
    )
  }

  if (!inherits(whitelist.markers, "data.frame")) {
    rlang::abort("whitelist.markers must be a data frame or a file path")
  }

  whitelist.markers <- dplyr::mutate(
    whitelist.markers,
    dplyr::across(tidyselect::everything(), as.character)
  )

  nrow.before <- nrow(whitelist.markers)
  whitelist.markers <- dplyr::distinct(whitelist.markers)
  duplicates <- nrow.before - nrow(whitelist.markers)

  if (duplicates > 0L && verbose) {
    message("Duplicated rows in whitelist of markers: ", duplicates)
    message("    Creating unique whitelist")
    message("    Downstream results may be affected; verify how the whitelist was generated")
  }

  whitelist.markers <- dplyr::mutate(
    whitelist.markers,
    dplyr::across(tidyselect::everything(), clean_markers_names)
  )

  if (tibble::has_name(whitelist.markers, "VARIANT_ID")) {
    whitelist.markers$VARIANT_ID <- as.integer(whitelist.markers$VARIANT_ID)
  }
  if (tibble::has_name(whitelist.markers, "M_SEQ")) {
    whitelist.markers$M_SEQ <- as.integer(whitelist.markers$M_SEQ)
  }

  if (verbose) {
    message("Number of whitelisted markers: ", nrow(whitelist.markers))
  }

  whitelist.markers
}
