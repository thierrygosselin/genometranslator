# read_gtypes ------------------------------------------------------------------

#' @name read_gtypes
#' @title Read a gtypes object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' \code{gtypes} support requires the optional GitHub package
#' \href{https://github.com/EricArcher/strataG}{\pkg{strataG}}. Consult its
#' repository for current installation and troubleshooting information.
#' @description Import a diploid \pkg{strataG} \code{gtypes} object and
#' standardize its sample, stratum, locus, and allele fields as a tidy genotype
#' table. Required data-slot columns and diploid allele counts are validated.

#' @param data A gtypes object (>= v.2.0.2), or its saved RDS file.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_gtypes
#' @return A tibble with \code{STRATA}, \code{INDIVIDUALS}, \code{MARKERS}, and
#' six-digit \code{GT}. Nucleotide alleles use a stable A/C/G/T mapping, numeric
#' alleles retain their values, and other labels are encoded within each locus.
#' @examples
#' \dontrun{
#' if (requireNamespace("strataG", quietly = TRUE)) {
#'   x <- readRDS("genotypes.gtypes.rds")
#'   genotypes <- genometranslator::read_gtypes(x)
#' }
#' }
#' @references Archer FI, Adams PE, Schneiders BB.
#' strataG: An r package for manipulating, summarizing and analysing population
#' genetic data.
#' Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies

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
  if (is.character(data) && length(data) == 1L) data <- readRDS(data)
  if (!inherits(data, "gtypes")) rlang::abort("Input is not a gtypes object.")

  # import ---------------------------------------------------------------------
  # lots of changes with gtypes so some stuff might be broken...

  input <- tibble::as_tibble(data@data)
  required.columns <- c("id", "stratum", "locus", "allele")
  missing.columns <- setdiff(required.columns, names(input))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "The gtypes data slot is missing required column(s): ",
      paste(missing.columns, collapse = ", "), "."
    ))
  }
  names(input)[match(required.columns, names(input))] <- c(
    "INDIVIDUALS", "STRATA", "MARKERS", "GT"
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
  gt.format <- sort(unique(as.character(input$GT[!is.na(input$GT)])))
  for (field in c("INDIVIDUALS", "STRATA", "MARKERS")) {
    input[[field]] <- as.character(input[[field]])
    if (anyNA(input[[field]]) || any(!nzchar(input[[field]])))
      rlang::abort("gtypes IDs, populations and loci cannot be missing.")
  }
  raw.ids <- unique(input$INDIVIDUALS)
  if (anyDuplicated(genometranslator::clean_ind_names(raw.ids)))
    rlang::abort("Individual names collide after cleaning.")
  if (length(gt.format) && all(grepl("^[0-9]+$", gt.format)) &&
      any(as.numeric(gt.format) > 999))
    rlang::abort("Numeric allele codes above 999 cannot fit the tidy GT format.")
  if (any(vapply(split(input$GT, input$MARKERS),
      function(z) length(unique(z[!is.na(z)])), integer(1)) > 999L))
    rlang::abort("More than 999 alleles at a locus cannot fit the tidy GT format.")

  # Nucleotide alleles use a stable A/C/G/T integer mapping. Numeric alleles
  # retain their codes; other labels are encoded per locus.
  if (length(gt.format) && all(gt.format %in% c("A", "C", "G", "T"))) {
    input$GT <- stringi::stri_replace_all_regex(
      str = input$GT,
      pattern = c("A", "C", "G", "T"),
      replacement = c("001", "002", "003", "004"),
      vectorize_all = FALSE
    )
  } else if (length(gt.format) && all(grepl("^[0-9]+$", gt.format))) {
    input$GT <- stringi::stri_pad_left(str = input$GT, pad = "0", width = 3)
  } else if (length(gt.format)) {
    input <- input %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::mutate(
        GT = dplyr::if_else(
          is.na(GT),
          NA_character_,
          stringi::stri_pad_left(
            as.character(match(as.character(GT), unique(as.character(GT[!is.na(GT)])))),
            width = 3,
            pad = "0"
          )
        )
      ) %>%
      dplyr::ungroup()
  }

  # For GT coded with only 1 number
  # gtypes.number <- unique(stringi::stri_count_boundaries(str = input$GT))
  # unique(stringi::stri_count_boundaries(str = test))

  # prep tidy ------------------------------------------------------------------
  allele.counts <- input %>%
    dplyr::count(STRATA, INDIVIDUALS, MARKERS, name = "n_alleles")
  if (any(allele.counts$n_alleles != 2L)) {
    rlang::abort("read_gtypes currently requires exactly two allele rows per individual and locus.")
  }

  input %<>%
    dplyr::mutate(
      GT = replace(GT, which(is.na(GT)), "000"),
      GT = stringi::stri_pad_left(str = GT, pad = "0", width = 3),
      INDIVIDUALS = clean_ind_names(INDIVIDUALS),
      STRATA = clean_pop_names(STRATA, factor = FALSE),
      MARKERS = as.character(MARKERS)
    ) %>%
    dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
    dplyr::summarise(GT = stringi::stri_join(GT, collapse = ""), .groups = "drop") %>%
    dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)

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


#' Write gtypes
#'
#' @description Create the target diploid object directly from GDS allele indices. Multiallelic loci are supported except by genlight, which requires biallelic loci.
#' @param data Open SeqArray GDS or GDS filename.
#' @param write Write the returned object as RDS.
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
#' write_gtypes("study.gds", strata = "samples.tsv")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_gtypes <- function(data, write = FALSE, filename = NULL,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE,
  verbose = TRUE) {
  start <- .legacy_start("write_gtypes", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(write, "write"); .legacy_package("strataG")
  x <- .legacy_snapshot(data, strata, TRUE, chunk.size = chunk.size)
  colnames(x$a) <- as.vector(rbind(paste0(x$loci$EXPORT_ID, ".1"),
    paste0(x$loci$EXPORT_ID, ".2")))
  out <- strataG::df2gtypes(data.frame(id = x$samples$INDIVIDUALS,
    strata = x$samples$GROUP, x$a, check.names = FALSE),
    ploidy = 2L, id.col = 1L, strata.col = 2L, loc.col = 3L)
  .legacy_object(out, x, if (write) if (is.null(filename)) "gtypes" else filename
    else NULL, path.folder, "_gtypes.rds", overwrite)
}

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
