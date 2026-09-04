#' Write HZAR
#'
#' @description Export ALT allele frequencies and called diploid sample sizes.
#' @param data Open SeqArray GDS or GDS path.
#' @param distances Table or TSV containing STRATA and Distance (finite numeric).
#' @param filename Output basename.
#' @param strata Sample metadata with INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Samples per read block.
#' @param overwrite Replace existing files.
#' @param verbose Display progress messages.
#' @return A table with Population, Distance, and paired L*.freq/L*.n columns.
#' Export paths are in the export.files attribute.
#' @details Biallelic diploid data only. Frequencies are calculated from current
#' called genotypes, not historical metrics. Entirely missing population/locus
#' cells have NA frequency and sample size zero. Distances are matched by label,
#' not row order. No missing distances are invented. GDS filters are restored.
#' @seealso \href{https://cran.r-project.org/package=hzar}{HZAR}
#' @examples
#' \dontrun{
#' write_hzar("study.gds", distances = "distances.tsv", strata = "samples.tsv")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_hzar <- function(data, distances, filename = NULL, strata = NULL,
  path.folder = getwd(), chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  start <- .legacy_start("write_hzar", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  x <- .legacy_snapshot(data, strata, TRUE, TRUE, chunk.size)
  if (is.character(distances) && length(distances)==1)
    distances <- readr::read_tsv(distances, show_col_types = FALSE, progress = FALSE)
  if (!is.data.frame(distances) || !all(c("STRATA","Distance") %in% names(distances)) ||
      anyNA(distances$STRATA) || anyDuplicated(distances$STRATA) ||
      !is.numeric(distances$Distance) || any(!is.finite(distances$Distance)))
    stop("Distances require unique STRATA and finite numeric Distance.")
  groups <- x$populations$GROUP
  r <- match(groups, distances$STRATA)
  if (anyNA(r)) stop("Missing distances for exported populations.")
  d <- .legacy_dosage(x)
  out <- tibble::tibble(Population = groups, Distance = distances$Distance[r])
  for (j in seq_len(ncol(d))) {
    n <- vapply(groups, function(p) sum(!is.na(d[x$samples$GROUP == p,j])), integer(1))
    count <- vapply(groups, function(p) sum(d[x$samples$GROUP == p,j],na.rm=TRUE), 0)
    out[[paste0(x$loci$EXPORT_ID[j],".freq")]] <- unname(ifelse(n>0, count/(2*n), NA_real_))
    out[[paste0(x$loci$EXPORT_ID[j],".n")]] <- unname(n)
  }
  paths <- .legacy_publish(x, filename, path.folder, "_hzar.csv", overwrite,
    function(p) readr::write_csv(out, p, na = "NA"))
  attr(out, "export.files") <- paths
  out
}

# Internal nested Function -----------------------------------------------------

#' @title generate_hzar
#' @description function to generate hzar function per groups of markers (to run in parallel)
#' @rdname generate_hzar
#' @keywords internal
#' @export

generate_hzar <- carrier::crate(function(x) {
  `%>%` <- magrittr::`%>%`
  freq.info <- x %>%
    dplyr::mutate(
      A1 = stringi::stri_sub(GT, 1, 3),
      A2 = stringi::stri_sub(GT, 4,6)
    ) %>%
    dplyr::select(MARKERS, POP_ID, INDIVIDUALS, A1, A2) %>%
    tidyr::gather(data = ., key = ALLELES_GROUP, value = ALLELES, -c(MARKERS, INDIVIDUALS, POP_ID)) %>%
    dplyr::select(-ALLELES_GROUP) %>%
    dplyr::group_by(MARKERS, POP_ID, ALLELES) %>%
    dplyr::tally(.) %>%
    dplyr::group_by(MARKERS, POP_ID) %>%
    dplyr::mutate(
      NN = sum(n),
      N = NN / 2) %>%
    dplyr::ungroup(.) %>%
    dplyr::mutate(FREQ = n / NN) %>%
    dplyr::select(-n, -NN)

  sample.n.info <- dplyr::distinct(freq.info, MARKERS, POP_ID, N) %>%
    dplyr::mutate(GROUP = stringi::stri_join(MARKERS, ".N")) %>%
    dplyr::select(MARKERS, POP_ID, GROUP, VALUE = N)

  freq.info <- dplyr::select(freq.info, MARKERS, POP_ID, ALLELES,FREQ) %>%
    dplyr::mutate(GROUP = stringi::stri_join(MARKERS, ALLELES, sep = ".")) %>%
    dplyr::select(MARKERS, POP_ID, GROUP, VALUE = FREQ) %>%
    dplyr::bind_rows(sample.n.info)
  sample.n.info <- NULL
  return(freq.info)
})#End generate_hzar
