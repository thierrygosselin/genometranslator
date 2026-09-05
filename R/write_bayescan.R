#' Write BayeScan
#'
#' @description Convert genotypes stored in GDS into the population allele-count
#' format required by BayeScan. The export preserves multiallelic loci, uses one
#' common allele dictionary across populations, and writes dictionaries that map
#' BayeScan's integer population and locus identifiers back to the GDS metadata.
#' @param data Open SeqArray GDS or GDS filename.
#' @param pop.select Optional populations to retain.
#' @param filename Output basename.
#' @param strata Metadata table or TSV containing INDIVIDUALS and STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Samples per read block.
#' @param overwrite Replace existing files.
#' @param verbose Display progress messages.
#' @return A list of files, pop.dictionary and markers.dictionary.
#' @details
#' \strong{Export safeguards:}
#' \itemize{
#' \item Population membership is obtained from the GDS individual metadata or
#' the supplied \code{strata} table. At least two populations are required.
#' \item Every exported locus must be polymorphic and have called genotypes in
#' every selected population. Violations stop the export; loci are not silently
#' removed.
#' \item Missing diploid genotypes contribute no allele copies. Partial diploid
#' calls are rejected because they cannot be represented unambiguously as
#' population allele counts.
#' \item Alleles are counted as non-negative integers using a single allele
#' dictionary for every population. This prevents population-specific allele
#' ordering from changing the biological meaning of a count column.
#' \item Population and marker dictionaries accompany the input file because
#' BayeScan replaces their original identifiers with consecutive integers.
#' \item Existing GDS selections are restored when the export finishes.
#' }
#'
#' These checks are deliberately stricter than BayeScan 2.1's native input
#' reader. The BayeScan reader assumes that population blocks are complete and
#' that the declared allele classes and their ordering are consistent among
#' populations; it does not comprehensively validate malformed, negative, or
#' internally inconsistent allele counts. Creating the input with this function
#' therefore reduces the risk of a malformed file being accepted or interpreted
#' incorrectly. It does not establish that BayeScan's statistical model is
#' appropriate for the study design.
#'
#' \strong{Interpretation:}
#' BayeScan is a legacy FST-outlier method based on a multinomial-Dirichlet
#' population model. Hierarchical structure, isolation by distance, admixture,
#' bottlenecks, uneven sampling, linked variants, and other departures from its
#' model can affect power or produce false candidates. Dense linked SNPs should
#' not be interpreted as independent confirmation of selection. BayeScan should
#' be used as a complementary sensitivity analysis, with important candidates
#' compared against structure-aware, environmental-association, haplotype-aware,
#' or simulation-based approaches appropriate to the study.
#' @references Foll M, Gaggiotti O (2008). A genome-scan method to identify
#' selected loci appropriate for both dominant and codominant markers:
#' a Bayesian perspective. Genetics 180:977-993.
#' \doi{10.1534/genetics.108.092221}
#'
#' Lotterhos KE, Whitlock MC (2014). Evaluation of demographic history and
#' neutral parameterization on the performance of FST outlier tests. Molecular
#' Ecology 23:2178-2192. \doi{10.1111/mec.12725}
#'
#' de Villemereuil P, Frichot E, Bazin E, Francois O, Gaggiotti OE (2014).
#' Genome scan methods against more complex models: when and how much should we
#' trust them? Molecular Ecology 23:2006-2019. \doi{10.1111/mec.12705}
#' @seealso \href{https://github.com/mfoll/BayeScan}{BayeScan format specification}
#' @examples
#' \dontrun{
#' write_bayescan("study.gds", strata = "samples.tsv", filename = "study")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_bayescan <- function(data, pop.select = NULL, filename = NULL,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE,
  verbose = TRUE) {
  start <- .legacy_start("write_bayescan", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(overwrite, "overwrite")
  x <- .legacy_snapshot(data, strata, TRUE, chunk.size = chunk.size)
  if (!is.null(pop.select)) {
    if (!length(pop.select) || anyNA(pop.select) ||
        !all(pop.select %in% x$samples$GROUP)) stop("Invalid pop.select.")
    keep <- x$samples$GROUP %in% pop.select
    x$samples <- x$samples[keep, ]; x$a <- x$a[keep, , drop = FALSE]
  }
  groups <- unique(x$samples$GROUP)
  if (length(groups) < 2L) stop("BayeScan requires at least two populations.")
  counts <- lapply(seq_len(nrow(x$loci)), function(j) {
    z <- vapply(groups, function(p)
      tabulate(as.vector(x$a[x$samples$GROUP == p, (2*j-1):(2*j), drop = FALSE]),
        nbins = x$counts[j]), integer(x$counts[j]))
    if (any(colSums(z) == 0L)) stop("Each locus needs called genotypes in every population.")
    if (sum(rowSums(z)>0) < 2L) stop("BayeScan requires polymorphic loci.")
    z
  })
  pop.dictionary <- tibble::tibble(STRATA = groups, BAYESCAN_POP = seq_along(groups))
  markers.dictionary <- tibble::tibble(MARKERS = x$loci$MARKERS,
    BAYESCAN_MARKERS = seq_len(nrow(x$loci)))
  paths <- .export_paths(filename, "bayescan", path.folder,
    c("_bayescan.txt", "_bayescan_pop_dictionary.tsv",
      "_bayescan_markers_dictionary.tsv", "_bayescan_alleles.tsv"), overwrite)
  staged <- vapply(paths, function(p) tempfile(), "")
  on.exit(unlink(staged), add = TRUE)
  con <- file(staged[1], "wt")
  tryCatch({
    writeLines(c(paste0("[loci]=", length(counts)), "",
      paste0("[populations]=", length(groups)), ""), con)
    for (p in seq_along(groups)) {
      writeLines(paste0("[pop]=", p), con)
      for (j in seq_along(counts))
        writeLines(paste(c(j, sum(counts[[j]][,p]), nrow(counts[[j]]),
          counts[[j]][,p]), collapse = " "), con)
      writeLines("", con)
    }
  }, finally = close(con))
  readr::write_tsv(pop.dictionary, staged[2])
  readr::write_tsv(markers.dictionary, staged[3])
  readr::write_tsv(x$alleles, staged[4])
  .export_publish(staged, paths, overwrite)
  list(files = paths, pop.dictionary = pop.dictionary,
    markers.dictionary = markers.dictionary)
}
