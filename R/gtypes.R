# read_gtypes ------------------------------------------------------------------

#' @name read_gtypes
#' @title Read a gtypes object into a tidy data frame
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' \code{gtypes} support requires the optional GitHub package
#' \href{https://github.com/EricArcher/strataG}{\pkg{strataG}}. Consult its
#' repository for current installation and troubleshooting information.
#' @description Transform a [strataG gtypes](https://github.com/EricArcher/strataG)
#' object into a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.

#' @param data A gtypes object (>= v.2.0.2) in the global environment.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_gtypes
#' @references Archer FI, Adams PE, Schneiders BB.
#' strataG: An r package for manipulating, summarizing and analysing population
#' genetic data.
#' Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

read_gtypes <- function(data, verbose = FALSE) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "read_gtypes",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("A gtypes object is required.")
  if (!inherits(data, "gtypes")) rlang::abort("Input is not a gtypes object.")

  # import ---------------------------------------------------------------------
  # lots of changes with gtypes so some stuff might be broken...

  input <- tibble::as_tibble(data@data) %>%
    dplyr::rename(
      INDIVIDUALS = id,
      STRATA = stratum,
      MARKERS = tidyselect::any_of("locus"),
      GT = tidyselect::any_of("allele")
      )

  # input <- suppressWarnings(
  #   tibble::as_tibble(data@data) %>%
  #     dplyr::rename(INDIVIDUALS = id, POP_ID = stratum) %>% # ids and strata before
  #     dplyr::mutate(ALLELES = rep(c("A1", "A2"), n() / 2)) %>%
  #     tidyr::pivot_longer(
  #       data = .,
  #       cols = -c("POP_ID", "INDIVIDUALS", "ALLELES"),
  #       names_to = "MARKERS",
  #       values_to = "GT"
  #     )
  # )
  # detect stratg genotype coding ----------------------------------------------
  # For GT = c("A", "C", "G", "T")
  gt.format <- sort(unique(input$GT))

  if (unique(gt.format %in% c("A", "C", "G", "T", NA))) {
    input$GT <- stringi::stri_replace_all_regex(
      str = input$GT,
      pattern = c("A", "C", "G", "T"),
      replacement = c("001", "002", "003", "004"),
      vectorize_all = FALSE
    )
  }

  # For GT = c("1", "2")
  if (unique(gt.format %in% c("1", "2", NA))) {
    input$GT <- stringi::stri_pad_left(str = input$GT, pad = "0", width = 3)
  }

  # For GT coded with only 1 number
  # gtypes.number <- unique(stringi::stri_count_boundaries(str = input$GT))
  # unique(stringi::stri_count_boundaries(str = test))

  # prep tidy ------------------------------------------------------------------
  input %<>%
    dplyr::mutate(
      GT = replace(GT, which(is.na(GT)), "000"),
      GT = stringi::stri_pad_left(str = GT, pad = "0", width = 3),
      STRATA = as.character(STRATA)
    ) %>%
    dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
    dplyr::summarise(GT = stringi::stri_join(GT, collapse = ""), .groups = "drop")

  ## before
  # input %<>%
  #   dplyr::mutate(
  #     GT = replace(GT, which(is.na(GT)), "000"),
  #     POP_ID = as.character(POP_ID)
  #     ) %>%
  #   dplyr::group_by(POP_ID, INDIVIDUALS, MARKERS) %>%
  #   tidyr::pivot_wider(data = ., names_from = "ALLELES", values_from = "GT") %>%
  #   dplyr::ungroup(.) %>%
  #   tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
  #   dplyr::select(POP_ID, INDIVIDUALS, MARKERS, GT)
  return(input)
}#End read_gtypes


# write_gtypes -----------------------------------------------------------------

#' @name write_gtypes
#' @title Write a strataG gtypes object from GDS or tidy data
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' \code{gtypes} support requires the optional GitHub package
#' \href{https://github.com/EricArcher/strataG}{\pkg{strataG}}. Consult its
#' repository for current installation and troubleshooting information.

#' @description Write a
#' \href{https://github.com/EricArcher/strataG}{strataG} object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.

#' @inheritParams genometranslator_common_arguments

#' @param write (logical, optional) To write in the working directory the gtypes
#' object. The file is written with \code{radiator_gtypes_DATE@TIME.RData} if no
#' filename is provided and can be open with load or readRDS.
#' Default: \code{write = FALSE}.
#'
#' @param filename (character, optional) Filename prefix.
#' Default: \code{filename = NULL}.

#' @return An object of the class \href{https://github.com/EricArcher/}{strataG} is returned.



#' @export
#' @rdname write_gtypes
#' @seealso \href{https://github.com/EricArcher/}{strataG}

#' @references Archer FI, Adams PE, Schneiders BB.
#' strataG: An r package for manipulating, summarizing and analysing population
#' genetic data.
#' Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559

#' @examples
#' \dontrun{
#' # require(strataG)
#' # with radiator GDS
#' turtle <- genometranslator::write_gtypes(data = "my.metadata.node.rad")
#'
#' # with tidy data
#' turtle <- genometranslator::write_gtypes(data = "my.radiator.rad")
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering

write_gtypes <- function(data, write = FALSE, filename = NULL) {
  # Loading the namespace registers the S4 `gtypes` class used below.
  if (!requireNamespace("strataG", quietly = TRUE)) {
    rlang::abort(paste0(
      "`write_gtypes()` requires the optional GitHub package strataG. ",
      "Install it with `remotes::install_github(\"EricArcher/strataG\")`."
    ))
  }

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection---------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)

  # Import data ----------------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {

    # Required package -----------------------------------------------------------
    tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)

    if (data.type == "gds.file") data %<>% genometranslator::read_genome(data = .)
    # biallelic <- detect_biallelic_markers(data)# faster with GDS
    markers.meta <- extract_markers_metadata(gds = data, markers.meta.select = "MARKERS", whitelist = TRUE)
    strata <- extract_individuals_metadata(gds = data, ind.field.select = c("INDIVIDUALS", "STRATA"), whitelist = TRUE)

    data <- SeqArray::seqGetData(
      gdsfile = data, var.name = "$dosage_alt") %>%
      magrittr::set_colnames(x = ., value = markers.meta$MARKERS) %>%
      magrittr::set_rownames(x = ., value = strata$INDIVIDUALS) %>%
      tgbase::trans_long(
        x = .,
        cols = "INDIVIDUALS",
        names_to = "MARKERS",
        values_to = "ALT_DOSAGE",
        keep_rownames = "INDIVIDUALS"
        ) %>%
      dplyr::left_join(strata, by = "INDIVIDUALS") %>%
      dplyr::mutate(
        `1` = dplyr::if_else(ALT_DOSAGE == 0L, 1L, ALT_DOSAGE),
        `2` = dplyr::recode(.x = ALT_DOSAGE, `1` = 2L, `0` = 1L),
        ALT_DOSAGE = NULL
      ) %>%
      tgbase::trans_long(
        x = .,
        cols = c("INDIVIDUALS", "STRATA", "MARKERS"),
        names_to = "ALLELES",
        values_to = "GT"
        ) %>%
      tgbase::trans_wide(
        x = .,
        formula = "STRATA + INDIVIDUALS ~ MARKERS + ALLELES",
        values_from = "GT",
        sep = "."
      ) %>%
      dplyr::arrange(STRATA, INDIVIDUALS)

    markers.meta <- strata <- NULL


  } else {#Tidy data
    data %<>% genometranslator::read_genome(data = ., import.metadata = TRUE)

    if (rlang::has_name(data, "ALT_DOSAGE")) {
      data  %<>%
        dplyr::select(MARKERS, STRATA, INDIVIDUALS, ALT_DOSAGE) %>%
        dplyr::mutate(
          `1` = dplyr::if_else(ALT_DOSAGE == 0L, 1L, ALT_DOSAGE),
          `2` = dplyr::recode(.x = ALT_DOSAGE, `1` = 2L, `0` = 1L),
          ALT_DOSAGE = NULL
        ) %>%
        tgbase::trans_long(
          x = .,
          cols = c("INDIVIDUALS", "STRATA", "MARKERS"),
          names_to = "ALLELES",
          values_to = "GT"
        ) %>%
        tgbase::trans_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS + ALLELES",
          values_from = "GT",
          sep = "."
        ) %>%
        dplyr::arrange(STRATA, INDIVIDUALS)
    } else {
      if (!rlang::has_name(data, "GT")) data %<>% calibrate_alleles(data = ., gt = TRUE) %$% input
      data %<>%
        dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT) %>%
        dplyr::arrange(MARKERS, STRATA, INDIVIDUALS) %>%
        dplyr::mutate(
          GT = replace(GT, which(GT == "000000"), NA),
          STRATA = as.character(STRATA),
          `1` = stringi::stri_sub(str = GT, from = 1, to = 3), # most of the time: faster than tidyr::separate
          `2` = stringi::stri_sub(str = GT, from = 4, to = 6),
          GT = NULL
        ) %>%
        tgbase::trans_long(
          x = .,
          cols = c("INDIVIDUALS", "STRATA", "MARKERS"),
          names_to = "ALLELES",
          values_to = "GT"
        ) %>%
        tgbase::trans_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS + ALLELES",
          values_from = "GT",
          sep = "."
        ) %>%
        dplyr::arrange(STRATA, INDIVIDUALS)
    }
  }

  # gtypes----------------------------------------------------------------------
  safe_gtypes <-  purrr::safely(.f = methods::new)

  res <- suppressWarnings(
    safe_gtypes(
      "gtypes",
      gen.data = data[,-c(1,2)],
      ploidy = 2,
      ind.names = data$INDIVIDUALS,
      strata = data$STRATA,
      schemes = NULL,
      sequences = NULL,
      description = NULL,
      other = NULL
    )
  )

  if (is.null(res$error)) {
    res <- res$result
  } else {
    rlang::abort("strataG package must be installed and loaded: library('strataG')")
  }

  if (write) {
    filename.temp <- generate_filename(name.shortcut = filename, extension = "gtypes")
    filename.short <- filename.temp$filename.short
    filename.gtypes <- filename.temp$filename
    saveRDS(object = res, file = filename.gtypes)
    message("File written: ", filename.short)
  }

  return(res)
}# End write_gtypes

# switch_genotypes -------------------------------------------------------------
#' @name switch_genotypes
#' @title switch_genotypes
#' @description todo
#' @rdname switch_genotypes
#' @keywords internal
#' @export
switch_genotypes <- function(x) {
  x <- dplyr::case_when(
    x == 1L ~ 2L,
    x == 2L ~ 2L,
    x == 0L ~ 1L
  )
}
