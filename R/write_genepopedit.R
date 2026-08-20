# write_genepopedit-----------------------------------------------------------------
#' @name write_genepopedit

#' @title Write a genepopedit flatten object

#' @description Write a genepopedit object from a tidy data frame or GDS file/object.
#'
#' Why not use \code{genepopedit::genepop_flatten}?
#' \itemize{
#' \item genepopedit requires a
#' \href{https://github.com/rystanley/genepopedit#genepopedit}{specific type of genepop},
#' so if you don't want to manipulate your genepop file, radiator is an alternative.
#' \item radiator follows guidelines highlighted here:
#' \href{http://genepop.curtin.edu.au/help_input.html}{genepop format}, but the
#' 3 functions in radiator that reads genepop files:
#' \enumerate{
#' \item \href{https://thierrygosselin.github.io/radiator/reference/genome_translator.html}{genome_translator}
#' \item \code{\link{tidy_genome}}
#' \item the underlying module: \href{https://thierrygosselin.github.io/genometranslator/reference/read_genepop.html}{read_genepop}
#' }
#' imports a \emph{larger variety of genepop alternatives}, similarly to
#' \href{https://github.com/thibautjombart/adegenet}{adegenet} \code{read.genepop}
#'  function, only faster.
#' }

#' @inheritParams genometranslator_common_arguments

#' @return A genepopedit object in the global environment.
#' @seealso \href{https://github.com/rystanley/genepopedit}{genepopedit}
#' @export
#' @rdname write_genepopedit
#' @references Stanley RRE, Jeffery NW, Wringe BF, DiBacco C, Bradbury IR (2017)
#' genepopedit: a simple and flexible tool for manipulating multilocus molecular
#' data in R. Molecular Ecology Resources, 17, 12-18.
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering
#' @template io-dependencies


write_genepopedit <- function(data) {
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") {
      data <- genometranslator::read_genome(data, verbose = FALSE)
    }
    data <- tidy_genome(data = data,
                     pop.id = FALSE,
                     calibrate.alleles = FALSE,
                     parallel.core = parallel::detectCores() - 1) %>%
      dplyr::select(INDIVIDUALS, STRATA, MARKERS, ALT_DOSAGE, REF, ALT)
    data.type <- "tbl_df"
  } else {
    if (is.vector(data)) {
      data <- genometranslator::read_genome(data = data, import.metadata = TRUE)
    }
  }

  if (!rlang::has_name(data, "GT")) data %<>% calibrate_alleles(data = ., gt = TRUE) %$% input

  if (!rlang::has_name(data, "STRATA") && rlang::has_name(data, "POP_ID")) {
    data %<>% dplyr::rename(STRATA = POP_ID)
  }

  # generate genepopedit flatten object
  data %<>%
    dplyr::mutate(SampleNum = INDIVIDUALS) %>%
    dplyr::select(SampleID = INDIVIDUALS, Population = STRATA, SampleNum, MARKERS, GT) %>%
    tgbase::trans_wide(x = ., formula = "SampleID + Population + SampleNum ~ MARKERS", values_from = "GT")
  return(data)
}# End write_genepopedit
