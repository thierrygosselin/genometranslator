# write_rubias ---------------------------------------------------------------

#' Write a rubias file
#'
#' @description Convert a supported genomic dataset to the diploid,
#' two-column-per-locus data frame used by \pkg{rubias}. The first four columns
#' are `sample_type`, `repunit`, `collection`, and `indiv`.
#'
#' @param data A GDS file or object, or another genomic object accepted by
#'   [read_genome()].
#' @param strata Optional sample metadata accepted by [read_strata()]. Supply
#'   `SAMPLE_TYPE`, `REPUNIT`, `COLLECTION`, and `INDIVIDUALS` for full control.
#'   Lower-case rubias names are also accepted. When these fields are absent,
#'   `STRATA` or `POP_ID` is used for both `repunit` and `collection`, and all
#'   samples are treated as references. Default: \code{strata = NULL}.
#' @param filename Optional filename prefix. When supplied, the result is written
#'   as `<filename>_rubias.tsv`; a timestamp is added rather than overwriting an
#'   existing file. Default: \code{filename = NULL}.
#' @param parallel.core Number of workers available while importing or
#'   transforming supported inputs. Retained for consistency with other writers.
#'   Default: \code{parallel.core = parallel::detectCores() - 1}.
#' @param verbose Logical. Display progress messages. Default:
#'   \code{verbose = TRUE}.
#'
#' @return A tibble compatible with rubias. When `filename` is supplied, the
#'   same table is also written as a tab-separated file.
#'
#' @details rubias requires one row per diploid individual and four leading
#' character columns:
#' \itemize{
#'   \item `sample_type` must be `"reference"` or `"mixture"`;
#'   \item reference samples require non-missing `repunit` and `collection`;
#'   \item mixture samples require `repunit = NA`, while `collection` identifies
#'     the mixture sample, port, stratum, place, or time group;
#'   \item `indiv` must uniquely identify every individual.
#' }
#' Missing genotypes are written as two `NA` alleles. Partial diploid genotypes
#' are rejected. Locus names cannot contain spaces because rubias uses adjacent
#' columns to identify gene copies.
#'
#' This writer creates an input representation; it does not run the conditional
#' GSI model and does not silently filter individuals or loci. Quality control,
#' marker selection, reference design, and the distinction between reference and
#' mixture samples remain the user's responsibility.
#'
#' @section Dependencies:
#' The output schema follows the documented interface of
#' \href{https://github.com/eriqande/rubias}{rubias}, developed by Eric C.
#' Anderson and Ben Moran. Installing rubias is not required to create the table,
#' but is required to analyse it. This implementation is independently written
#' from the public input specification; no rubias source code is incorporated.
#'
#' @references Moran BM, Anderson EC (2019). Bayesian inference from the
#' conditional genetic stock identification model. Canadian Journal of
#' Fisheries and Aquatic Sciences, 76(4), 551-560.
#' \doi{10.1139/cjfas-2018-0016}.
#'
#' @seealso \href{https://github.com/eriqande/rubias}{rubias}, [write_gsi_sim()]
#' @template writer-filtering
#' @template io-dependencies
#' @export
write_rubias <- function(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
) {
  .start <- tgbase::startup(
    package = "genometranslator", f.name = "write_rubias", verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  if (missing(data)) rlang::abort("`data` is required.")
  if (isTRUE(verbose)) message("Generating rubias data...")

  genome <- genometranslator::read_genome(
    data = data, import.metadata = TRUE, parallel.core = parallel.core
  )
  required <- c("INDIVIDUALS", "MARKERS")
  missing.columns <- setdiff(required, names(genome))
  if (length(missing.columns)) {
    rlang::abort(paste0(
      "The genomic data is missing: ", paste(missing.columns, collapse = ", "), "."
    ))
  }
  if (anyDuplicated(dplyr::select(genome, INDIVIDUALS, MARKERS))) {
    rlang::abort("Each individual-marker combination must occur at most once.")
  }
  if (any(grepl("\\s", as.character(genome$MARKERS)))) {
    rlang::abort("rubias locus names cannot contain spaces.")
  }

  # Sample metadata -----------------------------------------------------------
  if (!is.null(strata)) {
    if (is.character(strata) && length(strata) == 1L) {
      if (!file.exists(strata)) rlang::abort("The `strata` file does not exist.")
      metadata <- vroom::vroom(
        strata, col_types = vroom::cols(.default = "c"),
        show_col_types = FALSE, progress = FALSE
      )
    } else {
      metadata <- tibble::as_tibble(strata)
    }
    rubias.names <- c("sample_type", "repunit", "collection", "indiv")
    rubias.names.upper <- c("SAMPLE_TYPE", "REPUNIT", "COLLECTION", "INDIVIDUALS")
    if (!all(rubias.names %in% names(metadata)) &&
        !all(rubias.names.upper %in% names(metadata))) {
      metadata <- genometranslator::read_strata(
        strata = metadata, keep.two = FALSE, verbose = FALSE
      )$strata
    }
  } else {
    metadata.columns <- intersect(
      c("SAMPLE_TYPE", "REPUNIT", "COLLECTION", "INDIVIDUALS", "STRATA", "POP_ID",
        "sample_type", "repunit", "collection", "indiv"), names(genome)
    )
    metadata <- dplyr::distinct(genome, dplyr::across(dplyr::all_of(metadata.columns)))
  }

  rename.if.present <- function(x, old, new) {
    if (old %in% names(x) && !new %in% names(x)) names(x)[names(x) == old] <- new
    x
  }
  metadata <- rename.if.present(metadata, "SAMPLE_TYPE", "sample_type")
  metadata <- rename.if.present(metadata, "REPUNIT", "repunit")
  metadata <- rename.if.present(metadata, "COLLECTION", "collection")
  metadata <- rename.if.present(metadata, "INDIVIDUALS", "indiv")
  metadata <- rename.if.present(metadata, "POP_ID", "STRATA")

  complete.metadata <- all(
    c("sample_type", "repunit", "collection", "indiv") %in% names(metadata)
  )
  if (!complete.metadata) {
    if (!"indiv" %in% names(metadata)) {
      rlang::abort("Sample metadata must contain `INDIVIDUALS` or `indiv`.")
    }
    if (!"STRATA" %in% names(metadata)) {
      rlang::abort(paste0(
        "Supply the four rubias metadata fields, or provide `STRATA`/`POP_ID` ",
        "for reference samples."
      ))
    }
    metadata <- metadata |>
      dplyr::transmute(
        sample_type = "reference", repunit = as.character(STRATA),
        collection = as.character(STRATA), indiv = as.character(indiv)
      )
  } else {
    metadata <- metadata |>
      dplyr::transmute(
        sample_type = tolower(as.character(sample_type)),
        repunit = as.character(repunit), collection = as.character(collection),
        indiv = as.character(indiv)
      )
  }
  metadata <- dplyr::distinct(metadata)
  if (anyDuplicated(metadata$indiv)) {
    rlang::abort("`indiv`/`INDIVIDUALS` values must be unique in sample metadata.")
  }
  if (length(setdiff(unique(metadata$sample_type), c("reference", "mixture")))) {
    rlang::abort("`sample_type` must contain only `reference` or `mixture`.")
  }
  reference <- metadata$sample_type == "reference"
  mixture <- metadata$sample_type == "mixture"
  if (any(reference & (is.na(metadata$repunit) | !nzchar(metadata$repunit)))) {
    rlang::abort("Every reference sample requires a non-missing `repunit`.")
  }
  if (any(reference & (is.na(metadata$collection) | !nzchar(metadata$collection)))) {
    rlang::abort("Every reference sample requires a non-missing `collection`.")
  }
  if (any(mixture & !is.na(metadata$repunit) & nzchar(metadata$repunit))) {
    rlang::abort("Mixture samples must use `repunit = NA`.")
  }
  if (any(mixture & (is.na(metadata$collection) | !nzchar(metadata$collection)))) {
    rlang::abort("Every mixture sample requires a non-missing `collection`.")
  }

  genome.ids <- unique(as.character(genome$INDIVIDUALS))
  if (!setequal(genome.ids, metadata$indiv)) {
    rlang::abort(paste0(
      "Genomic data and sample metadata contain different individual IDs. ",
      "Missing from metadata: ", paste(setdiff(genome.ids, metadata$indiv), collapse = ", "),
      "; missing from genomic data: ", paste(setdiff(metadata$indiv, genome.ids), collapse = ", "), "."
    ))
  }

  # Diploid alleles -----------------------------------------------------------
  dosage_to_alleles <- function(dosage) {
    dosage <- as.integer(dosage)
    list(
      a1 = dplyr::case_when(dosage %in% 0:2 ~ "0", TRUE ~ NA_character_),
      a2 = dplyr::case_when(
        dosage == 0L ~ "0", dosage %in% c(1L, 2L) ~ "1", TRUE ~ NA_character_
      )
    )
  }
  split_vcf_gt <- function(gt) {
    gt <- as.character(gt)
    parts <- strsplit(gt, "[/|]", perl = TRUE)
    a1 <- vapply(parts, function(x) if (length(x) == 2L) x[1L] else NA_character_, "")
    a2 <- vapply(parts, function(x) if (length(x) == 2L) x[2L] else NA_character_, "")
    missing.gt <- is.na(gt) | gt %in% c("./.", ".|.")
    a1[missing.gt] <- a2[missing.gt] <- NA_character_
    list(a1 = a1, a2 = a2)
  }

  if ("ALT_DOSAGE" %in% names(genome)) {
    alleles <- dosage_to_alleles(genome$ALT_DOSAGE)
  } else if ("GT" %in% names(genome) &&
      all(is.na(genome$GT) | as.character(genome$GT) %in% c("0", "1", "2"))) {
    alleles <- dosage_to_alleles(genome$GT)
  } else if ("GT_VCF_NUC" %in% names(genome)) {
    alleles <- split_vcf_gt(genome$GT_VCF_NUC)
  } else if ("GT_VCF" %in% names(genome)) {
    alleles <- split_vcf_gt(genome$GT_VCF)
  } else if ("GT" %in% names(genome) &&
      all(is.na(genome$GT) | grepl("^[0-9]{6}$", as.character(genome$GT)))) {
    gt <- as.character(genome$GT)
    a1 <- substr(gt, 1L, 3L)
    a2 <- substr(gt, 4L, 6L)
    missing.gt <- is.na(gt) | gt == "000000"
    a1[missing.gt] <- a2[missing.gt] <- NA_character_
    alleles <- list(a1 = a1, a2 = a2)
  } else {
    rlang::abort(paste0(
      "Cannot derive diploid alleles. Supply `ALT_DOSAGE`, dosage `GT`, ",
      "`GT_VCF_NUC`, `GT_VCF`, or six-digit `GT`."
    ))
  }
  if (any(xor(is.na(alleles$a1), is.na(alleles$a2)))) {
    rlang::abort("rubias does not accept partially missing diploid genotypes.")
  }

  allele.data <- tibble::tibble(
    indiv = as.character(genome$INDIVIDUALS),
    marker = as.character(genome$MARKERS), A1 = alleles$a1, A2 = alleles$a2
  )
  markers <- unique(allele.data$marker)
  genotype.wide <- tidyr::pivot_wider(
    allele.data, id_cols = indiv, names_from = marker, values_from = c(A1, A2),
    names_glue = "{marker}.{.value}"
  )
  allele.columns <- as.vector(rbind(paste0(markers, ".A1"), paste0(markers, ".A2")))
  genotype.wide <- dplyr::select(genotype.wide, indiv, dplyr::all_of(allele.columns))

  output <- metadata |>
    dplyr::left_join(genotype.wide, by = "indiv") |>
    dplyr::select(sample_type, repunit, collection, indiv, dplyr::everything())

  if (!is.null(filename)) {
    output.file <- paste0(filename, "_rubias.tsv")
    if (file.exists(output.file)) {
      output.file <- paste0(filename, "_rubias_", format(Sys.time(), "%Y%m%d@%H%M"), ".tsv")
    }
    readr::write_tsv(output, output.file)
    if (isTRUE(verbose)) message("File written: ", basename(output.file))
  }
  output
}
