#' Common arguments used by genometranslator
#'
#' This documentation-only helper centralizes parameters shared by import,
#' translation, and export functions.
#'
#' @name genometranslator_common_arguments
#' @rdname genometranslator_common_arguments
#' @keywords internal
#' @export
#'
#' @param interactive.filter Logical indicating whether an interactive session
#'   may ask for thresholds.
#'   Default: \code{interactive.filter = TRUE}.
#' @param gds A genome GDS file path or object.
#' Default: \code{gds = NULL}.
#' @param data A supported genomic file, object, or tidy genomic data frame.
#' Default: \code{data = NULL}.
#' @param parallel.core Number of workers available for parallel operations.
#'   Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @param verbose Logical indicating whether progress messages are emitted.
#'   Default: \code{verbose = TRUE}.
#' @param random.seed Optional integer seed for operations involving randomness.
#'   Default: \code{random.seed = NULL}.
#' @param ... Additional arguments passed to lower-level readers, translators,
#'   or writers.
#'
#' @return `NULL`, invisibly. This function exists to share documentation.
genometranslator_common_arguments <- function(
    interactive.filter = TRUE,
    gds = NULL,
    data = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE,
    random.seed = NULL,
    ...
) {
  invisible(NULL)
}
