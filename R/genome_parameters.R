# genome_parameters -----------------------------------------------------------
#' Track changes to genomic data
#'
#' Creates, initiates, or updates a history of operations that change a genomic
#' dataset. For a genometranslator GDS, the cumulative history is persisted in
#' the GDS metadata and each generated TSV is a snapshot of that history at the
#' beginning of the current operation. The completed operation is then appended
#' to both the TSV and the GDS.
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
#' @return A list containing the genomic summary, current parameter row,
#'   cumulative `filter.history`, and history-file path as applicable.
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
    filter.history <- if (inherits(data, "SeqVarGDSClass")) {
      read_filter_history(data)
    } else {
      empty_filter_history()
    }
    res$filters.parameters <- empty_filter_history()
    res$filter.history <- filter.history
    readr::write_tsv(filter.history, parameter.path)
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

    res$filters.parameters <- normalize_filter_history(tibble::tibble(
      FILTERS = filter.name,
      PARAMETERS = param.name,
      VALUES = if (is.null(values)) "not filtering" else values,
      BEFORE = paste(before, collapse = " / "),
      AFTER = paste(after, collapse = " / "),
      BLACKLIST = paste(before - after, collapse = " / "),
      UNITS = units,
      COMMENTS = comments
    ))
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
    if (inherits(data, "SeqVarGDSClass")) {
      filter.history <- dplyr::bind_rows(
        read_filter_history(data),
        res$filters.parameters
      )
      write_filter_history(data, filter.history)
      res$filter.history <- filter.history
    }
    res$info <- info.new
    res$filters.parameters.path <- parameter.obj$filters.parameters.path
  }

  res
}

#' Import legacy filter-parameter files into a GDS
#'
#' Imports one or more historical `filters_parameters_*.tsv` files into the
#' persistent filtering history stored in a genometranslator GDS. Each element
#' of `paths` may be a parameter file or a filter-results directory containing
#' one or more parameter files. Files are imported in the order supplied.
#'
#' @param gds A writable genome GDS filepath or open `SeqVarGDSClass` object.
#' @param paths Character vector of parameter files or filter-results folders.
#' @param replace Replace the existing GDS history instead of appending to it.
#' @param filename Optional TSV filepath for a snapshot of the resulting
#'   cumulative history.
#' @param verbose Display progress messages.
#'
#' @return Invisibly returns the cumulative filter-history tibble.
#' @export
import_filter_history <- function(
    gds,
    paths,
    replace = FALSE,
    filename = NULL,
    verbose = TRUE
) {
  if (!length(paths) || anyNA(paths)) {
    rlang::abort("`paths` must contain at least one file or directory.")
  }
  files <- unlist(lapply(as.character(paths), function(path) {
    if (dir.exists(path)) {
      list.files(
        path,
        pattern = "^filters_parameters.*[.]tsv$",
        full.names = TRUE
      )
    } else if (file.exists(path)) {
      path
    } else {
      rlang::abort(paste0("Filter-history path does not exist: ", path))
    }
  }), use.names = FALSE)
  if (!length(files)) {
    rlang::abort("No `filters_parameters_*.tsv` files were found.")
  }
  imported <- purrr::map_dfr(files, function(file) {
    normalize_filter_history(readr::read_tsv(
      file,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    ))
  })

  opened.here <- FALSE
  if (is.character(gds) && length(gds) == 1L && file.exists(gds)) {
    gds <- read_genome(gds, verbose = FALSE)
    opened.here <- TRUE
  }
  if (!inherits(gds, "SeqVarGDSClass")) {
    rlang::abort("`gds` must be a GDS filepath or open SeqVarGDSClass object.")
  }
  on.exit({
    if (opened.here) try(SeqArray::seqClose(gds), silent = TRUE)
  }, add = TRUE)

  history <- if (replace) {
    imported
  } else {
    dplyr::bind_rows(read_filter_history(gds), imported)
  }
  write_filter_history(gds, history)
  if (!is.null(filename)) readr::write_tsv(history, filename)
  if (verbose) {
    message(
      "Imported ", nrow(imported), " filtering operation(s) from ",
      length(files), " file(s)."
    )
  }
  invisible(history)
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
