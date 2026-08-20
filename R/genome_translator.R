# genome_translator -----------------------------------------------------------

#' Translate genomic data between formats
#'
#' A lightweight composition of \code{\link{read_genome}} and
#' \code{\link{write_genome}}. The input format is detected automatically and
#' read with the corresponding format-specific reader. The resulting genome is
#' then sent to the requested writer.
#'
#' Use a format-specific \code{read_*()} or \code{write_*()} function when you
#' need arguments beyond their defaults.
#'
#' @param data A supported genomic file or object.
#' @param strata Optional strata data or filename passed to readers and writers
#'   that support it. Default: \code{strata = NULL}.
#' @param output Character vector naming one or more output formats.
#' @param filename Optional output filename or prefix. Default:
#'   \code{filename = NULL}.
#' @param parallel.core Number of processor cores passed to readers and writers
#'   that support parallel processing. Default:
#'   \code{parallel.core = parallel::detectCores() - 1}.
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#'
#' @return Invisibly returns the result produced by \code{write_genome()}.
#' Writers whose purpose is a file side effect may return \code{NULL}.
#' @export
#'
#' @examples
#' \dontrun{
#' genome_translator(
#'   data = "genomes.vcf",
#'   strata = "strata.tsv",
#'   output = "genepop",
#'   filename = "genomes.gen"
#' )
#' }
genome_translator <- function(
    data,
    strata = NULL,
    output,
    filename = NULL,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE
) {
  if (missing(data)) rlang::abort("`data` is required.")
  if (missing(output) || is.null(output) || !length(output)) {
    rlang::abort("`output` must name at least one genomic output format.")
  }

  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "genome_translator",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  input.type <- detect_genomic_format(data)
  genome <- read_genome(
    data = data,
    strata = strata,
    parallel.core = parallel.core,
    verbose = verbose
  )

  opened.internally <- !identical(input.type, "SeqVarGDSClass") &&
    inherits(genome, "SeqVarGDSClass")
  if (opened.internally) {
    on.exit(try(SeqArray::seqClose(genome), silent = TRUE), add = TRUE)
  }

  result <- write_genome(
    data = genome,
    output = output,
    filename = filename,
    strata = strata,
    parallel.core = parallel.core,
    verbose = verbose
  )

  invisible(result)
}
