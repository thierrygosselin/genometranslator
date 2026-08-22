# Read package-native genomic data ---------------------------------------------

#' @name read_genome
#' @title Read genomic data
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' This dispatcher uses the dependencies documented by the selected
#' format-specific \code{read_*()} function. Legacy FST/RAD files additionally
#' require the optional CRAN package \pkg{fst}; VCF preparation may require the
#' optional \command{bcftools} executable.
#' @description Detect a supported genomic format and route it to the
#' appropriate format-specific reader. The generic interface uses the defaults
#' of each \code{read_*()} function; call that reader directly when additional
#' control is required. Package-native tabular inputs are normalized with the
#' internal \code{as_tidy_genome()} helper.

#' The function uses \code{\link[arrow]{read_parquet}} or
#' CoreArray Genomic Data Structure (\href{https://github.com/zhengxwen/gdsfmt}{GDS})
#' file system.


#' @param data A file in the working directory ending with .arrow.parquet or .gds,
#' a TSV or legacy RAD/FST file, or an existing wide/tidy genomic table.
#' @param strata Optional strata data or filename passed to readers that support
#' it. Default: \code{strata = NULL}.
#' @param parallel.core Number of processor cores passed to readers that support
#' parallel processing. Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @param columns (optional) For arrow.parquet file.
#' Column names to read.
#' The default is to read all all columns.
#' Default: \code{columns = NULL}.
#' @param allow.dup (optional, logical) To allow the opening of a GDS file with
#' read-only mode when it has been opened in the same R session.
#' Default: \code{allow.dup = FALSE}.
#' @param check (optional, logical) Verify that GDS number of samples and markers
#' match.
#' Default: \code{check = TRUE}.
#' @param import.metadata Logical. Retain columns in addition to the standard
#' genomic columns when normalizing tabular input.
#' Default: \code{import.metadata = TRUE}.
#' @param verbose Logical. Display progress messages. For GDS input, the
#' current number of samples and markers and a summary of active filters are
#' displayed. Default: \code{verbose = TRUE}.

#' @details For GDS file system, \strong{read_genome} will open the GDS connection file
#' set the filters (variants and samples) based on the info found in the file.

#' @return A tidy genomic data frame or
#' GDS object (with read/write permissions) in the global environment.
#' @export
#' @rdname read_genome
#' @seealso
#' \href{https://github.com/apache/arrow/}{arrow}
#' \href{https://github.com/zhengxwen/gdsfmt}{GDS}
#' \code{\link{read_genome}}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @examples
#' \dontrun{
#' shark <- genometranslator::read_genome(data = "data.shark.gds")
#' turtle <- genometranslator::read_genome(data = "data.turtle.arrow.parquet")
#' }

read_genome <- function(
    data,
    strata = NULL,
    columns = NULL,
    allow.dup = FALSE,
    check = TRUE,
    import.metadata = TRUE,
    parallel.core = parallel::detectCores() - 1,
    verbose = TRUE
) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "read_genome",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  # detect format---------------------------------------------------------------
  data.type <- detect_genomic_format(data)

  if (data.type %in% c("genind", "genlight", "gtypes") &&
      is.character(data) && length(data) == 1L) {
    data <- readRDS(data)
  }

  # Format-specific readers ---------------------------------------------------
  # These calls intentionally use each reader's defaults. Users who need
  # format-specific options should call the corresponding reader directly.
  if (identical(data.type, "dart")) {
    return(read_dart(
      data = data,
      strata = strata,
      parallel.core = parallel.core,
      verbose = verbose
    ))
  }
  if (identical(data.type, "fstat.file")) {
    return(read_fstat(data = data, strata = strata, verbose = verbose))
  }
  if (identical(data.type, "genepop.file")) {
    return(read_genepop(data = data, strata = strata, verbose = verbose))
  }
  if (identical(data.type, "genind")) {
    return(read_genind(data = data, verbose = verbose))
  }
  if (identical(data.type, "genlight")) {
    return(read_genlight(
      data = data,
      parallel.core = parallel.core,
      verbose = verbose
    ))
  }
  if (identical(data.type, "gtypes")) {
    return(read_gtypes(data = data, verbose = verbose))
  }
  if (data.type %in% c("plink.tped.file", "plink.bed.file")) {
    return(read_plink(
      data = data,
      parallel.core = parallel.core,
      verbose = verbose
    ))
  }
  if (identical(data.type, "vcf.file")) {
    return(read_vcf(
      data = data,
      strata = strata,
      parallel.core = parallel.core,
      verbose = verbose
    ))
  }

  # An existing GDS connection needs no file import.
  if (identical(data.type, "SeqVarGDSClass")) {
    if (verbose) gds_open_summary(data)
    return(data)
  }

  # arrow parquet file ---------------------------------------------------------
  if ("arrow.parquet" %in% data.type) {
    # for some results, arrow doesn't like having option B only when col_select
    # is NULL from the function call...
    if (is.null(columns)) { #option A
      data <- arrow::read_parquet(file = data)
    } else { #option B
      data <- arrow::read_parquet(file = data, col_select = columns)
    }
    return(as_tidy_genome(data, import.metadata = import.metadata))
  }

  # Legacy FST/RAD file --------------------------------------------------------
  if (identical(data.type, "fst.file")) {
    tgbase::check_package(package = "fst")
    data <- fst::read_fst(path = data, as.data.table = FALSE)
    return(as_tidy_genome(data, import.metadata = import.metadata))
  }

  # TSV or in-memory tabular input ---------------------------------------------
  if (identical(data.type, "tbl_df")) {
    if (is.character(data) && length(data) == 1L) {
      data <- readr::read_tsv(
        file = data,
        col_types = readr::cols(.default = readr::col_character())
      )
    }
    return(as_tidy_genome(data, import.metadata = import.metadata))
  }



  # GDS file -------------------------------------------------------------------
  if ("gds.file" %in% data.type) {
    if (verbose) message("Opening GDS file connection")

    seq_open_temp <- function(data, allow.dup) {
      SeqArray::seqOpen(gds.fn = data, readonly = allow.dup, allow.duplicate = allow.dup)
    }#End seq_open_temp

    safe_seq_open <- purrr::safely(.f = seq_open_temp)

    data.safe <- safe_seq_open(data, allow.dup)

    if (is.null(data.safe$error)) {
      data <- data.safe$result
    } else {
      temp <- gdsfmt::openfn.gds(filename = data, readonly = FALSE, allow.fork = FALSE, allow.duplicate = TRUE)
      if (gdsfmt::get.attr.gdsn(temp$root)$FileFormat == "SNP_ARRAY") {
        gdsfmt::closefn.gds(gdsfile = temp)
        tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)
        temp.name <- stringi::stri_join("radiator_temp_", data)
        message("The input file is a SNP GDS file, conversion using SeqArray::seqSNP2GDS")
        SeqArray::seqSNP2GDS(data, temp.name)
        data <- SeqArray::seqOpen(gds.fn = temp.name, readonly = allow.dup, allow.duplicate = allow.dup)
        message("SeqArray GDS file generated: ", temp.name)
      }
    }

    s <- extract_individuals_metadata(
      gds = data,
      ind.field.select = "INDIVIDUALS",
      whitelist = TRUE
    ) %$%
      INDIVIDUALS
    m <- extract_markers_metadata(
      gds = data,
      markers.meta.select = "VARIANT_ID",
      whitelist = TRUE
    ) %$%
      VARIANT_ID
    SeqArray::seqSetFilter(object = data,
                           variant.id = m,
                           sample.id = as.character(s),
                           verbose = FALSE)

    # Checks--------------------------------------------------------------------
    if (check) {
      check <- SeqArray::seqGetFilter(data)
      if (length(check$sample.sel[check$sample.sel]) != length(s)) {
        rlang::abort("Number of samples don't match, contact author")
      }
      if (length(check$variant.sel[check$variant.sel]) != length(m)) {
        rlang::abort("Number of markers don't match, contact author")
      }
    }
    if (verbose) gds_open_summary(data)
    return(data)
  }#End gds.file

  rlang::abort(paste0(
    "Input format `", data.type, "` is not supported by read_genome()."
  ))
}#End read_genome


# Normalize a wide or tidy genomic table ---------------------------------------
as_tidy_genome <- function(data, import.metadata = FALSE) {
  data %<>% dplyr::rename(STRATA = tidyselect::any_of("POP_ID"))

  if (rlang::has_name(data, "LOCUS") && !rlang::has_name(data, "MARKERS")) {
    data %<>% dplyr::rename(MARKERS = LOCUS)
  }

  if (!"MARKERS" %in% colnames(data) && !"LOCUS" %in% colnames(data)) {
    grouping.cols <- "INDIVIDUALS"
    if (rlang::has_name(data, "STRATA")) {
      grouping.cols <- c("STRATA", "INDIVIDUALS")
    }
    data %<>%
      tgbase::trans_long(
        x = .,
        cols = grouping.cols,
        names_to = "MARKERS",
        values_to = "GT",
        variable_factor = FALSE
      )
  }

  if (!import.metadata) {
    want <- c(
      "STRATA", "INDIVIDUALS", "MARKERS", "CHROM", "LOCUS", "POS",
      "GT", "GT_VCF_NUC", "GT_VCF", "ALT_DOSAGE"
    )
    data %<>% dplyr::select(tidyselect::any_of(want))
  }

  if (rlang::has_name(data, "MARKERS")) {
    data$MARKERS %<>% clean_markers_names(.)
  }
  data$INDIVIDUALS %<>% clean_ind_names(.)
  if (rlang::has_name(data, "STRATA")) {
    data$STRATA %<>% clean_pop_names(.)
  }

  dplyr::ungroup(data)
}

# write_genome------------------------------------------------------------------
#' @name write_genome
#' @title Write genomic data
#' @description Write genomic data using a format-specific \code{write_*()}
#' function. The output format can be supplied explicitly or inferred from an
#' unambiguous filename extension. Call the specialized writer directly when
#' additional control is required. With no \code{output}, the historical
#' Parquet-writing and GDS-closing behaviour is retained.
#'
#'
#' When the object is a CoreArray Genomic Data Structure
#' (\href{https://github.com/zhengxwen/gdsfmt}{GDS}) file system, the function
#' \strong{close the connection with the GDS file}. Before doing so it sets the
#' filters (variants and samples) based on the info found in the file.
#'
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data An object in the global environment: tidy genomic dataset
#' or GDS connection file
#' @param output Optional output format. Supported values include \code{gds},
#' \code{genepop}, \code{fstat}, \code{genind}, \code{genlight}, \code{gtypes},
#' \code{plink}, \code{vcf}, and the other formats having a corresponding
#' package \code{write_*()} function. Default: \code{output = NULL}.
#' @param strata Optional strata data passed to writers that support it.
#' Default: \code{strata = NULL}.
#' @param parallel.core Number of processor cores passed to writers that support
#' it. Default: \code{parallel.core = parallel::detectCores() - 1}.

#' @param filename Name of the Arrow/Parquet file written for tabular data.
#' The argument is not used when closing a GDS connection.

#' @param internal (optional, logical) This is used inside radiator internal code and it stops
#' from writing the file.
#' Default: \code{internal = FALSE}.

#' @param write.message (optional, character) Print a message in the console
#' after writing file.
#' With \code{write.message = NULL}, nothing is printed in the console.
#' Default: \code{write.message = "standard"}. This will print
#' \code{message("File written: ", basename(filename))}.

#' @param verbose (optional, logical) \code{verbose = TRUE} to be chatty
#' during execution.
#' Default: \code{verbose = FALSE}.

#' @return A file written in the working directory or nothing if it's a GDS connection file.
#' @export
#' @rdname write_genome
#'
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @seealso
#' \href{https://arrow.apache.org}{appache arrow}
#'
#' \href{https://github.com/zhengxwen/gdsfmt}{GDS}
#'
#' \code{\link{read_genome}}
#'
#' @examples
#' \dontrun{
#' require(SeqArray)
#' genometranslator::write_genome(data = tidy.data, filename = "data.shark.arrow.parquet")
#' genometranslator::write_genome(data = gds.object)
#' }
#' @template writer-filtering
#' @template io-dependencies


write_genome <- function(
    data,
    output = NULL,
    filename = NULL,
    strata = NULL,
    parallel.core = parallel::detectCores() - 1,
    internal = FALSE,
    write.message = "standard",
    verbose = FALSE
) {

  # Infer formats only when the extension is unambiguous. ---------------------
  if (is.null(output) && !is.null(filename)) {
    lower.filename <- tolower(filename)
    output <- dplyr::case_when(
      grepl("\\.vcf(?:\\.gz)?$", lower.filename) ~ "vcf",
      grepl("\\.gds$", lower.filename) &&
        !inherits(data, "SeqVarGDSClass") ~ "gds",
      grepl("\\.gen$", lower.filename) ~ "genepop",
      grepl("\\.bed$|\\.ped$|\\.tped$", lower.filename) ~ "plink",
      grepl("\\.parquet$", lower.filename) ~ "parquet",
      TRUE ~ NA_character_
    )
    if (is.na(output)) output <- NULL
  }

  # Generic writer routing ----------------------------------------------------
  if (!is.null(output) && !identical(tolower(output), "parquet")) {
    output <- tolower(output)
    aliases <- c(
      fstat = "hierfstat",
      seqarray = "gds"
    )
    output <- unname(ifelse(output %in% names(aliases), aliases[output], output))

    writer.names <- stats::setNames(
      paste0("write_", output),
      output
    )
    missing.writers <- !vapply(
      writer.names,
      exists,
      logical(1),
      mode = "function",
      inherits = TRUE
    )
    if (any(missing.writers)) {
      rlang::abort(paste0(
        "Unsupported output format: ",
        paste(names(writer.names)[missing.writers], collapse = ", "),
        ". Call a format-specific writer if one is available."
      ))
    }

    results <- vector("list", length(writer.names))
    names(results) <- names(writer.names)
    for (i in seq_along(writer.names)) {
      writer <- get(writer.names[[i]], mode = "function", inherits = TRUE)
      writer.args <- list(
        data = data,
        filename = filename,
        strata = strata,
        parallel.core = parallel.core,
        verbose = verbose,
        write = TRUE,
        open = TRUE
      )
      accepted <- names(formals(writer))
      writer.args <- writer.args[names(writer.args) %in% accepted]
      results[[i]] <- do.call(writer, writer.args)
    }

    if (length(results) == 1L) return(results[[1L]])
    return(results)
  }

  if (!internal) {
    # detect format---------------------------------------------------------------
    data.type <- class(data)

    # GDS closing connection and setting filters
    if ("SeqVarGDSClass" %in% data.type) {
      # samples
      s <- extract_individuals_metadata(
        gds = data,
        ind.field.select = "INDIVIDUALS",
        whitelist = TRUE
      ) %$%
        INDIVIDUALS

      # markers
      m <- extract_markers_metadata(
        gds = data,
        markers.meta.select = "VARIANT_ID",
        whitelist = TRUE
      ) %$%
        VARIANT_ID

      if (verbose) message("Setting filters to:")
      if (verbose) message("    number of samples: ", length(s))
      if (verbose) message("    number of markers: ", length(m))
      # ?SeqArray::seqSetFilter
      gds.filename <- data$filename
      SeqArray::seqClose(data)
      data <- SeqArray::seqOpen(gds.fn = gds.filename, readonly = FALSE, allow.duplicate = FALSE)
      SeqArray::seqSetFilter(object = data,
                             variant.id = m,
                             sample.id = as.character(s),
                             verbose = verbose)
      if (verbose) message("Closing connection with GDS file:\n", basename(gds.filename))
      SeqArray::seqClose(data)
      return(gds.filename)
    }

    # using arrow parquet
    if (is.null(filename)) rlang::abort("A filename must be provided")

    file.type <- stringi::stri_extract(
      str = filename,
      regex = "\\.[^\\.]*$"
    )

    # check for .rad file ending
    if (".rad" %in% file.type) {
      message("fst file format deprecated, see function doc")
      message("Changing .rad to .arrow.parquet")
      filename <- stringi::stri_sub_replace(str = filename, from = -4, to = -1, replacement = ".arrow.parquet")
      file.type <- ".parquet"
    }

    if (!".parquet" %in% file.type) {
      rlang::abort("filename with .arrow.parquet must be provided")
    }

    if (verbose) cli::cli_progress_step(msg = "Writing arrow parquet tidy dataset...")

    # writing arrow parquet file
    tibble::as_tibble(data) %>%
      arrow::write_parquet(x = ., sink = filename)

    if (verbose) cli::cli_progress_done()

    if (!is.null(write.message) && verbose) {
      if (write.message == "standard") {
        message("File written: ", basename(filename))
      } else {
        write.message
      }
    }
  }

}#End write_genome
