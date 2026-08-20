# write a pcadapt file from a tidy data frame

#' @name write_pcadapt
#' @title Write a \href{https://github.com/bcm-uga/pcadapt}{pcadapt}
#' file from a tidy data frame

#' @description Write a
#' \href{https://github.com/bcm-uga/pcadapt}{pcadapt}
#' file from a tidy data frame. The data is biallelic.
#' Used internally in
#' \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.
#'
#' Prepare an appropriately filtered biallelic dataset before export. Consider
#' missingness, minor-allele frequency, linkage disequilibrium, and the sampling
#' design intended for the pcadapt analysis.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{genometranslator} \code{\link{tidy_genome}}.

#' @inheritParams read_strata
#' @inheritParams tidy_genome

#' @param filename (optional) The file name prefix for the pcadapt file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_pcadapt_}.

#' Default: \code{filename = NULL}.
#' @section Details:
#'
#' \strong{Use a filtered dataset}:
#' \enumerate{
#' \item \strong{Control linkage disequilibrium}:
#' Reducing linkage before running genome scan is essential. At least start by
#' removing SNPs on the same RADseq locus (short linkage disequilibrium).
#'
#' \item \strong{Control minor alleles}: Too much Minor Alleles is just noise in the data.
#' Filter using Count, Frequency or Depth.
#'
#' \item \strong{Only use markers that are in common between strata}
#' \item \strong{Only use polymorphic markers}
#'
#' }

#' @return A pcadapt file is written in the working directory a genotype matrix
#' object is also generated in the global environment.

#' @param pop.select (optional, string) Selected list of populations for
#' the analysis. e.g. \code{pop.select = c("QUE", "ONT")} to select \code{QUE}
#' and \code{ONT} population samples (out of 20 pops). If \code{pop.labels}
#' argument was used to rename the strata column, use the new names with
#' \code{pop.select}.
#' Default: \code{pop.select = NULL}.
#' 
#' @export
#' @rdname write_pcadapt
#' @references Luu, K., Bazin, E., & Blum, M. G. (2017).
#' pcadapt: an R package to perform genome scans for selection based on principal component analysis.
#' Molecular Ecology Resources, 17(1), 67-77.

#' @references Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G. (2015).
#' Detecting genomic signatures of natural selection with principal component analysis: application to the 1000 Genomes data.
#' Molecular biology and evolution, msv334.

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering
#' @template io-dependencies


write_pcadapt <- function(
  data,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
) {

  message("Generating pcadapt file...")

  file.date <- format(Sys.time(), "%Y%m%d@%H%M")


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file is missing")

  # Import data ---------------------------------------------------------------
  if (is.vector(data)) data %<>% genometranslator::read_genome(data = ., import.metadata = TRUE)

  # pop.select -----------------------------------------------------------------
  if (!is.null(pop.select)) {
    message("pop.select: ")
    data %<>% dplyr::filter(STRATA %in% pop.select)
    if (is.factor(data$STRATA)) data$STRATA <- droplevels(data$STRATA)
  }


  # detect biallelic markers ---------------------------------------------------
  biallelic <- detect_biallelic_markers(data = data)

  if (!biallelic) rlang::abort("\npcadapt only work with biallelic dataset")


  # Biallelic and ALT_DOSAGE -------------------------------------------------------

  n.ind <- dplyr::n_distinct(data$INDIVIDUALS)
  n.pop <- dplyr::n_distinct(data$STRATA)
  n.markers <- dplyr::n_distinct(data$MARKERS)


  if (!rlang::has_name(data, "ALT_DOSAGE")) {
    data %<>% calibrate_alleles(
      data = ., biallelic = TRUE, alt.dosage = TRUE) %$% input
  }

  data %<>% dplyr::select(MARKERS, INDIVIDUALS, STRATA, ALT_DOSAGE)

  pop.string <- data %>%
    dplyr::distinct(STRATA, INDIVIDUALS) %>%
    dplyr::arrange(STRATA, INDIVIDUALS) %>%
    dplyr::select(STRATA)

  pop.string <- pop.string$STRATA

  data %<>%
    dplyr::select(MARKERS, INDIVIDUALS, STRATA, ALT_DOSAGE) %>%
    dplyr::arrange(STRATA, INDIVIDUALS, MARKERS) %>%
    dplyr::select(-STRATA) %>%
    dplyr::mutate(ALT_DOSAGE = replace(ALT_DOSAGE, which(is.na(ALT_DOSAGE)), 9)) %>%
    tgbase::trans_wide(x = ., formula = "MARKERS ~ INDIVIDUALS", values_from = "ALT_DOSAGE") %>% # could be the other way ...
    dplyr::select(-MARKERS)

  # writing file to directory  ------------------------------------------------

  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_pcadapt_", file.date, ".txt")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_pcadapt_", file.date, ".txt")
    } else {
      filename <- stringi::stri_join(filename, "_pcadapt", ".txt")
    }
  }

  message("writing pcadapt file with:
    Number of populations: ", n.pop, "\n    Number of individuals: ", n.ind,
          "\n    Number of markers: ", n.markers)

  readr::write_delim(x = data, file = filename, col_names = FALSE,
                     append = FALSE, delim = " ")

  data <- as.matrix(data)
  res <- list(genotype.matrix = data, pop.string = pop.string)
  return(res)
}# End write_pcadapt

