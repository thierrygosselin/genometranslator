# write_bayescan ---------------------------------------------------------------
# write a bayescan file from a tidy data frame

#' @name write_bayescan
#' @title Write a BayeScan file

#' @description Write a \href{http://cmpg.unibe.ch/software/BayeScan/}{BayeScan}
#' file from a tidy data frame. The data is bi- or multi-allelic.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.

#' @param data A GDS file, open GDS object, tidy data frame object, or a tidy
#' data frame in wide or long format in the working directory. For GDS input,
#' individual metadata stored in the GDS supplies the `STRATA` column.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{genometranslator} \code{\link[genometranslator]{tidy_genome}}.

#' @inheritParams genometranslator::read_strata
#' @inheritParams genometranslator::tidy_genome

#' @param filename (optional) The file name prefix for the bayescan file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{genometranslator_bayescan_}.

#' Default: \code{filename = NULL}.
#' @return A bayescan file is written in the working directory.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = TRUE}.
#' 
#' @param parallel.core Number of workers available for parallel operations.
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#' 
#' @param pop.select (optional, string) Selected list of populations for
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#' and \code{ONT} population samples (out of 20 pops). If \code{pop.labels}
#' argument was used to rename the strata column, use the new names with
#' \code{pop.select}.
#' Default: \code{pop.select = NULL}.
#' 
#' @export
#' @rdname write_bayescan

#' @references Foll, M and OE Gaggiotti (2008) A genome scan method to identify
#' selected loci appropriate
#' for both dominant and codominant markers: A Bayesian perspective.
#' Genetics 180: 977-993

#' @references Foll M, Fischer MC, Heckel G and L Excoffier (2010)
#' Estimating population structure from
#' AFLP amplification intensity. Molecular Ecology 19: 4638-4647

#' @references Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP
#' genome scans detect
#' local adaptation in high-altitude populations of a small rodent (Microtus arvalis).
#' Molecular Ecology 20: 1450-1462

#' @details BayeScan input should contain polymorphic markers represented in
#' every selected stratum. Missingness, minor-allele thresholds, linkage
#' disequilibrium, and other quality-control decisions should be addressed by
#' the user before calling this writer. The function validates common-marker
#' representation and polymorphism, but does not remove failing markers.
#'
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering
#' @template io-dependencies


write_bayescan <- function(
  data,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "write_bayescan",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)

  if (verbose) message("Generating BayeScan file...")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file is missing")

  # Import data ---------------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)
    data <- genometranslator::tidy_genome(
      data = data,
      pop.id = FALSE,
      calibrate.alleles = FALSE,
      parallel.core = parallel.core,
      close.gds = identical(data.type, "gds.file")
    )
  } else if (is.vector(data)) {
    data %<>%
      genometranslator::read_genome(import.metadata = TRUE)
  }

  if (!rlang::has_name(data, "STRATA")) {
    rlang::abort(
      "BayeScan requires population assignments in a `STRATA` column."
    )
  }

  # necessary steps to make sure we work with unique markers and not duplicated LOCUS
  if (rlang::has_name(data, "LOCUS") && !rlang::has_name(data, "MARKERS")) {
    data %<>% dplyr::rename(MARKERS = LOCUS)
  }

  # make sure we use POP_ID and not STRATA here...
  # if (rlang::has_name(data, "STRATA")) data %<>% dplyr::rename(POP_ID = STRATA)

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    message("pop.select: ")
    data %<>% dplyr::filter(STRATA %in% pop.select)
    if (is.factor(data$STRATA)) data$STRATA <- droplevels(data$STRATA)
  }

  # Keep markers represented in every selected stratum -------------------------
  n.strata.total <- dplyr::n_distinct(data$STRATA)
  common.markers <- data %>%
    dplyr::distinct(MARKERS, STRATA) %>%
    dplyr::count(MARKERS, name = "N_STRATA") %>%
    dplyr::filter(N_STRATA == n.strata.total) %>%
    dplyr::pull(MARKERS)
  non.common.markers <- setdiff(unique(data$MARKERS), common.markers)
  if (length(non.common.markers)) {
    rlang::abort(paste0(
      "BayeScan requires markers represented in every selected stratum. ",
      length(non.common.markers), " marker(s) fail this requirement. ",
      "Filter common markers before calling `write_bayescan()`."
    ))
  }

  # detect biallelic markers ---------------------------------------------------
  biallelic <- genometranslator::detect_biallelic_markers(data = data)

  if (!biallelic) {
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "STRATA", "GT_VCF_NUC", "GT")
    data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want)))
    if (rlang::has_name(data, "GT_VCF_NUC")) {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "STRATA", "GT_VCF_NUC")
      data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want))) %>%
        dplyr::rename(GT_HAPLO = GT_VCF_NUC)
    } else {
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "INDIVIDUALS", "STRATA", "GT")
      data <- suppressWarnings(dplyr::select(data, dplyr::one_of(want))) %>%
        dplyr::rename(GT_HAPLO = GT)
    }

    data <- genometranslator::calibrate_alleles(
      data = data,
      biallelic = FALSE,
      parallel.core = parallel.core,
      verbose = TRUE)$input
  }

  # Biallelic and ALT_DOSAGE -------------------------------------------------------
  if (biallelic) {
    data <- genometranslator::calibrate_alleles(
      data = data,
      biallelic = TRUE,
      parallel.core = parallel.core,
      verbose = TRUE
      ) %$%
      input %>%
      dplyr::select(MARKERS, INDIVIDUALS, STRATA, ALT_DOSAGE)
  }

  # Validate polymorphism without filtering -----------------------------------
  if (biallelic) {
    polymorphic.markers <- data %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        CALLED = sum(!is.na(ALT_DOSAGE)),
        ALT_COUNT = sum(ALT_DOSAGE, na.rm = TRUE),
        REF_COUNT = 2L * CALLED - ALT_COUNT,
        .groups = "drop"
      ) %>%
      dplyr::filter(CALLED > 0L, ALT_COUNT > 0L, REF_COUNT > 0L) %>%
      dplyr::pull(MARKERS)
  } else {
    gt.column <- if ("GT_VCF_NUC" %in% names(data)) "GT_VCF_NUC" else "GT_VCF"
    polymorphic.markers <- data %>%
      dplyr::group_by(MARKERS) %>%
      dplyr::summarise(
        POLYMORPHIC = {
          genotypes <- .data[[gt.column]]
          genotypes <- genotypes[!is.na(genotypes) & genotypes != "./."]
          alleles <- unlist(strsplit(genotypes, "[/|]"), use.names = FALSE)
          length(unique(alleles)) > 1L
        },
        .groups = "drop"
      ) %>%
      dplyr::filter(POLYMORPHIC) %>%
      dplyr::pull(MARKERS)
  }
  monomorphic.markers <- setdiff(unique(data$MARKERS), polymorphic.markers)
  if (length(monomorphic.markers)) {
    rlang::abort(paste0(
      "BayeScan requires polymorphic markers. ",
      length(monomorphic.markers), " monomorphic or entirely missing marker(s) ",
      "were detected. Filter them before calling `write_bayescan()`."
    ))
  }

  # prep data wide format ------------------------------------------------------
  n.ind <- dplyr::n_distinct(data$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(data$STRATA)
  n.markers <- dplyr::n_distinct(data$MARKERS)

  data %<>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(
      BAYESCAN_POP = factor(STRATA),
      BAYESCAN_POP = as.integer(BAYESCAN_POP),
      BAYESCAN_MARKERS = factor(MARKERS),
      BAYESCAN_MARKERS = as.integer(BAYESCAN_MARKERS)
    )

  pop.dictionary <- dplyr::distinct(data, STRATA, BAYESCAN_POP)
  markers.dictionary <- dplyr::distinct(data, MARKERS, BAYESCAN_MARKERS) %>%
    dplyr::arrange(BAYESCAN_MARKERS)

  data %<>% dplyr::select(-STRATA, -MARKERS)

  # writing file to directory  ------------------------------------------------
  # Filename: date and time to have unique filenaming
  if (is.null(filename)) {
    filename <- stringi::stri_join("genometranslator_bayescan_", file.date, ".txt")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_bayescan_", file.date, ".txt")
    } else {
      filename <- stringi::stri_join(filename, "_bayescan", ".txt")
    }
  }

  if (biallelic) {
    markers.type <- "biallelic"
  } else {
    markers.type <- "multiallelic"
  }

  message("writing BayeScan file with:
          Number of populations: ", n.pop, "\n    Number of individuals: ", n.ind,
          "\n    Number of ", markers.type, " markers: ", n.markers)

  # Number of markers
  readr::write_file(x = stringi::stri_join("[loci]=", n.markers, "\n\n"), file = filename, append = FALSE)

  # Number of populations
  readr::write_file(x = stringi::stri_join("[populations]=", n.pop, "\n\n"), file = filename, append = TRUE)
  pop.string <- unique(data$BAYESCAN_POP)
  generate_bayescan_biallelic <- function(pop, data) {
    # pop <- "BEA"
    data.pop <- dplyr::filter(data, BAYESCAN_POP %in% pop) %>%
      dplyr::filter(!is.na(ALT_DOSAGE)) %>%
      dplyr::group_by(BAYESCAN_MARKERS) %>%
      dplyr::summarise(
        REF = (length(ALT_DOSAGE[ALT_DOSAGE == 0]) * 2) + (length(ALT_DOSAGE[ALT_DOSAGE == 1])),
        ALT = (length(ALT_DOSAGE[ALT_DOSAGE == 2]) * 2) + (length(ALT_DOSAGE[ALT_DOSAGE == 1]))
      ) %>%
      dplyr::mutate(GENE_N = REF + ALT, ALLELE_N = rep(2, dplyr::n())) %>%
      dplyr::select(BAYESCAN_MARKERS, GENE_N, ALLELE_N, REF, ALT)
    readr::write_file(x = stringi::stri_join("[pop]=", pop, "\n"), file = filename, append = TRUE)
    readr::write_delim(x = data.pop, file = filename, append = TRUE, delim = "  ")
    readr::write_file(x = stringi::stri_join("\n"), file = filename, append = TRUE)
  }
  generate_bayescan_multiallelic <- function(data) {
    pop <- unique(data$BAYESCAN_POP)
    data.pop <- dplyr::select(data, -BAYESCAN_POP)
    readr::write_file(x = stringi::stri_join("[pop]=", pop, "\n"), file = filename, append = TRUE)
    # readr::write_delim(x = data.pop, file = filename, append = TRUE, delim = "  " )
    utils::write.table(x = data.pop, file = filename, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
    readr::write_file(x = stringi::stri_join("\n"), file = filename, append = TRUE)
  }

  if (!biallelic) {
    data.prep <- data %>%
      dplyr::select(GT_VCF, BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::filter(GT_VCF != "./.") %>%
      dplyr::mutate(
        A1 = stringi::stri_sub(str = GT_VCF, from = 1, to = 1),
        A2 = stringi::stri_sub(str = GT_VCF, from = 3, to = 3),
        GT_VCF = NULL
      ) %>%
      tidyr::pivot_longer(
        data = .,
        cols = -c("BAYESCAN_MARKERS", "BAYESCAN_POP"),
        names_to = "ALLELES_GROUP",
        values_to = "ALLELES"
      ) %>%
      dplyr::select(-ALLELES_GROUP)

    allele.count <- data.prep %>%
      dplyr::distinct(BAYESCAN_MARKERS, ALLELES) %>%
      dplyr::group_by(BAYESCAN_MARKERS) %>%
      dplyr::tally(.) %>%
      dplyr::rename(COUNT = n)

    data.prep <- data.prep %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP, ALLELES) %>%
      dplyr::tally(.) %>%
      dplyr::ungroup(.) %>%
      tidyr::complete(data = ., BAYESCAN_POP, tidyr::nesting(BAYESCAN_MARKERS, ALLELES), fill = list(n = 0)) %>%
      dplyr::ungroup(.)

    alleles.markers <- data.prep %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(GENE_N = sum(n)) %>%
      dplyr::ungroup(.) %>%
      dplyr::left_join(allele.count, by = "BAYESCAN_MARKERS") %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(GENE_N = stringi::stri_join(GENE_N, COUNT, sep = " "))

    data.prep <- data.prep %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP, ALLELES) %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(ALLELES = stringi::stri_join(n, collapse = " ")) %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP)

    data <- dplyr::left_join(
      alleles.markers, data.prep, by = c("BAYESCAN_MARKERS", "BAYESCAN_POP")) %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::group_by(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      dplyr::summarise(GT = stringi::stri_join(GENE_N, ALLELES, sep = " ")) %>%
      dplyr::ungroup(.) %>%
      dplyr::arrange(BAYESCAN_MARKERS, BAYESCAN_POP) %>%
      split(x = ., f = .$BAYESCAN_POP)
    data.prep <- alleles.markers <- allele.count <- NULL

    purrr::walk(.x = data, .f = generate_bayescan_multiallelic)
  } else {
    purrr::walk(.x = pop.string, .f = generate_bayescan_biallelic, data = data)
  }


  message("Writting populations dictionary")
  readr::write_tsv(
    x = pop.dictionary,
    file = stringi::stri_replace_all_fixed(
      str = filename, pattern = ".txt",
      replacement = "_pop_dictionary.tsv", vectorize_all = FALSE))
  message("Writting markers dictionary")
  readr::write_tsv(
    x = markers.dictionary,
    file = stringi::stri_replace_all_fixed(
      str = filename, pattern = ".txt",
      replacement = "_markers_dictionary.tsv", vectorize_all = FALSE))

  res <- list(pop.dictionary = pop.dictionary, markers.dictionary = markers.dictionary)
  return(res)
}# End write_bayescan
