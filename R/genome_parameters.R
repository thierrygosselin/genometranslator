# genome_parameters -----------------------------------------------------------
#' Track changes to genomic data
#'
#' Creates, initiates, or updates a tab-separated history of operations that
#' change a genomic dataset.
#'
#' @param generate,initiate,update Logical controls for creating, initiating,
#'   and updating the parameter history.
#' @param parameter.obj Existing parameter-history object.
#' Default: \code{parameter.obj = NULL}.
#' @param data A tidy genomic data frame or supported GDS object.
#' Default: \code{data = NULL}.
#' @param filter.name Name of the operation applied to the data.
#' Default: \code{filter.name = ""}.
#' @param param.name Name of the parameter controlling the operation.
#' Default: \code{param.name = ""}.
#' @param values Parameter value recorded in the history.
#' Default: \code{values = paste(NULL, NULL, sep = " / ")}.
#' @param units Units represented in the before and after summaries.
#' Default: \code{units = "individuals / strata / chrom / locus / markers"}.
#' @param comments Optional comments.
#' Default: \code{comments = ""}.
#' @param path.folder Output directory.
#' Default: \code{path.folder = NULL}.
#' @param file.date Date-time label used in the filename.
#' Default: \code{file.date = NULL}.
#' @param internal,verbose Logical controls for internal use and messages.
#'
#' @return A list containing the genomic summary, parameter row, and history
#'   file path as applicable.
#' @param verbose 
#' Default: \code{verbose = TRUE}.
#' 
#' @param internal 
#' Default: \code{internal = FALSE}.
#' 
#' @param update 
#' Default: \code{update = TRUE}.
#' 
#' @param initiate 
#' Default: \code{initiate = FALSE}.
#' 
#' @param generate 
#' Default: \code{generate = FALSE}.
#' 
#' @export
genome_parameters <- function(
    generate = FALSE,
    initiate = FALSE,
    update = TRUE,
    parameter.obj = NULL,
    data = NULL,
    filter.name = "",
    param.name = "",
    values = paste(NULL, NULL, sep = " / "),
    units = "individuals / strata / chrom / locus / markers",
    comments = "",
    path.folder = NULL,
    file.date = NULL,
    internal = FALSE,
    verbose = TRUE
) {
  if (internal && !verbose) return(NULL)
  res <- list()
  if (is.null(file.date)) file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (is.null(path.folder)) path.folder <- getwd()

  if (!is.null(parameter.obj) && generate && !initiate) {
    generate <- initiate <- update <- FALSE
    res <- parameter.obj
  }
  if (!is.null(parameter.obj) && generate && initiate) generate <- FALSE
  if (is.null(parameter.obj) && update) {
    rlang::abort("`parameter.obj = NULL` is not accepted when `update = TRUE`.")
  }
  if (internal) verbose <- FALSE

  if (generate) {
    parameter.path <- generate_filename(
      name.shortcut = "filters_parameters",
      path.folder = path.folder,
      date = file.date,
      extension = "tsv"
    )$filename

    parameter.obj$filters.parameters.path <- parameter.path
    res$filters.parameters.path <- parameter.path
    res$filters.parameters <- tibble::tibble(
      FILTERS = character(),
      PARAMETERS = character(),
      VALUES = character(),
      BEFORE = character(),
      AFTER = character(),
      BLACKLIST = integer(),
      UNITS = character(),
      COMMENTS = character()
    )
    readr::write_tsv(res$filters.parameters, parameter.path)
    if (verbose) message("Genome parameters file generated: ", basename(parameter.path))
  }

  if (initiate) {
    if (is.null(data)) rlang::abort("A GDS or tidy genomic data object is required.")
    res$info <- parameter.obj$info <- genome_info(data)
    res$filters.parameters.path <- parameter.obj$filters.parameters.path
  }

  if (update) {
    if (is.null(data)) rlang::abort("A GDS or tidy genomic data object is required.")
    info <- parameter.obj$info
    info.new <- genome_info(data)
    before <- unname(unlist(info[c("n.ind", "n.pop", "n.chrom", "n.locus", "n.snp")]))
    after <- unname(unlist(info.new[c("n.ind", "n.pop", "n.chrom", "n.locus", "n.snp")]))

    res$filters.parameters <- tibble::tibble(
      FILTERS = filter.name,
      PARAMETERS = param.name,
      VALUES = if (is.null(values)) "not filtering" else values,
      BEFORE = paste(before, collapse = " / "),
      AFTER = paste(after, collapse = " / "),
      BLACKLIST = paste(before - after, collapse = " / "),
      UNITS = units,
      COMMENTS = comments
    )
    readr::write_tsv(
      res$filters.parameters,
      parameter.obj$filters.parameters.path,
      append = TRUE,
      col_names = FALSE
    )
    if (verbose) message(
      "Genome parameters file updated: ",
      basename(parameter.obj$filters.parameters.path)
    )
    res$info <- info.new
    res$filters.parameters.path <- parameter.obj$filters.parameters.path
  }

  res
}


# genome_info -----------------------------------------------------------------
#' Summarise genomic data dimensions
#'
#' @param x A tidy genomic data frame or supported GDS object.
#' @param verbose Display the summary.
#'
#' Default: \code{verbose = FALSE}.
#' @return A list containing counts of chromosomes, loci, markers, strata, and
#'   individuals.
#' @export
genome_info <- function(x, verbose = FALSE) {
  if (inherits(x, "tbl_df")) {
    res <- list(
      n.pop = if ("POP_ID" %in% names(x)) {
        dplyr::n_distinct(x$POP_ID)
      } else if ("STRATA" %in% names(x)) {
        dplyr::n_distinct(x$STRATA)
      } else NA_integer_,
      n.ind = if ("INDIVIDUALS" %in% names(x)) dplyr::n_distinct(x$INDIVIDUALS) else NA_integer_,
      n.snp = if ("MARKERS" %in% names(x)) dplyr::n_distinct(x$MARKERS) else NA_integer_,
      n.locus = if ("LOCUS" %in% names(x)) dplyr::n_distinct(x$LOCUS) else NA_integer_,
      n.chrom = if ("CHROM" %in% names(x)) dplyr::n_distinct(x$CHROM) else NA_integer_
    )
  } else {
    markers <- extract_markers_metadata(
      gds = x,
      markers.meta.select = c("CHROM", "LOCUS", "MARKERS", "VARIANT_ID"),
      whitelist = TRUE
    )
    marker.id <- if (
      "VARIANT_ID" %in% names(markers) && !all(is.na(markers$VARIANT_ID))
    ) markers$VARIANT_ID else markers$MARKERS
    individuals <- extract_individuals_metadata(
      gds = x,
      ind.field.select = c("STRATA", "INDIVIDUALS"),
      whitelist = TRUE
    )
    res <- list(
      n.pop = dplyr::n_distinct(individuals$STRATA),
      n.ind = dplyr::n_distinct(individuals$INDIVIDUALS),
      n.snp = dplyr::n_distinct(marker.id),
      n.locus = dplyr::n_distinct(markers$LOCUS),
      n.chrom = dplyr::n_distinct(markers$CHROM)
    )
  }

  if (verbose) {
    message("Number of chrom: ", res$n.chrom)
    message("Number of locus: ", res$n.locus)
    message("Number of markers: ", res$n.snp)
    message("Number of strata: ", res$n.pop)
    message("Number of individuals: ", res$n.ind)
  }
  res
}
