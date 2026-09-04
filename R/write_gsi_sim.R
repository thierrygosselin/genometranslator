#' Write gsi_sim
#'
#' @description Export active diploid GDS genotypes in the gsi_sim text format.
#' Genotypes are read in sample blocks; no full tidy genotype table is created.
#' @param data An open SeqArray GDS object or a GDS filename.
#' @param pop.levels Character vector listing every population once, in export
#' order. NULL uses first appearance in the active sample order.
#' @param pop.labels Optional non-empty labels corresponding to pop.levels.
#' Repeated labels explicitly merge populations.
#' @param strata Optional sample metadata table or TSV containing INDIVIDUALS
#' and the grouping column. NULL uses the GDS metadata.
#' @param filename Output basename (optionally ending in .txt).
#' @param group.column Population column in metadata. Default: STRATA.
#' @param path.folder Output directory.
#' @param chunk.size Maximum individuals read per genotype block.
#' @param overwrite Allow replacement of existing output files.
#' @param verbose Display progress messages.
#' @return Invisibly, a list containing files and sample, population, locus and
#' allele mapping tables.
#' @details Only diploid GDS input is accepted. Active selections and metadata
#' whitelists are respected and restored, including after errors. No additional
#' biological filtering is performed. Multiallelic loci are supported.
#' Alleles are numbered from one within each GDS locus; zero denotes missing.
#' Partially missing calls are rejected, not silently imputed.
#' Safe S/P/L identifiers are used with companion mapping TSV files.
#' Export baseline and mixture samples separately, using the same marker selection
#' and allele dictionaries. The writer does not run gsi_sim.
#' Phase, coverage and genotype likelihoods are not represented by this format.
#' @references Anderson EC, Waples RS, Kalinowski ST (2008). An improved method
#' for predicting the accuracy of genetic stock identification. Canadian Journal
#' of Fisheries and Aquatic Sciences 65:1475-1486.
#' Anderson EC (2010). Assessing the power of informative subsets of loci for
#' population assignment: standard methods are upwardly biased. Molecular Ecology
#' Resources 10:701-710.
#' @seealso \href{https://github.com/eriqande/gsi_sim}{Official gsi_sim documentation},
#' [write_rubias()]
#' @examples
#' \dontrun{
#' write_gsi_sim("study.gds", strata = "samples.tsv", filename = "baseline")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_gsi_sim <- function(data, pop.levels = NULL, pop.labels = NULL,
  strata = NULL, filename = "gsi_sim", group.column = "STRATA",
  path.folder = getwd(), chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  .export_flag(verbose, "verbose"); .export_flag(overwrite, "overwrite")
  .export_count(chunk.size, "chunk.size", 1)
  start <- tgbase::startup(package = "genometranslator",
    f.name = "write_gsi_sim", verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  gds <- context$gds
  meta <- .export_groups(gds, strata, group.column, pop.levels)
  samples <- meta$samples; populations <- meta$populations
  if (!is.null(pop.labels)) {
    if (!is.character(pop.labels) || length(pop.labels) != nrow(populations) ||
        anyNA(pop.labels) || any(!nzchar(trimws(pop.labels))))
      stop("pop.labels must provide one non-empty label per population.")
    samples$GROUP <- pop.labels[match(samples$GROUP, populations$GROUP)]
    populations <- tibble::tibble(GROUP = unique(pop.labels),
      EXPORT_ID = paste0("P", seq_along(unique(pop.labels))))
  }
  samples$EXPORT_ID <- paste0("S", seq_len(nrow(samples)))
  info <- .gsi_loci(gds)
  filename <- sub("\\.txt$", "", filename, ignore.case = TRUE)
  paths <- .export_paths(filename, "gsi_sim", path.folder,
    c(".txt", "_samples.tsv", "_populations.tsv", "_loci.tsv", "_alleles.tsv"),
    overwrite)
  staged <- vapply(paths, function(x) tempfile(), "")
  on.exit(unlink(staged), add = TRUE)
  con <- file(staged[1], "wt")
  on.exit(if (!is.null(con)) close(con), add = TRUE)
  writeLines(c(paste(nrow(samples), nrow(info$loci)), info$loci$EXPORT_ID), con)
  for (p in seq_len(nrow(populations))) {
    rows <- which(samples$GROUP == populations$GROUP[p])
    writeLines(paste("POP", populations$EXPORT_ID[p]), con)
    for (first in seq.int(1L, length(rows), by = chunk.size)) {
      r <- rows[first:min(length(rows), first + chunk.size - 1L)]
      block <- .gsi_block(gds, samples$INDIVIDUALS[r], info$counts)
      block[is.na(block)] <- 0L
      for (k in seq_along(r))
        writeLines(paste(c(samples$EXPORT_ID[r[k]], block[k, ]), collapse = " "), con)
    }
  }
  close(con)
  # Remove the connection cleanup once closed; use a sentinel instead.
  con <- NULL
  readr::write_tsv(samples, staged[2])
  readr::write_tsv(populations, staged[3])
  readr::write_tsv(info$loci, staged[4])
  readr::write_tsv(info$alleles, staged[5])
  .export_publish(staged, paths, overwrite)
  if (verbose) .export_message("gsi_sim export: ", nrow(samples), " samples; ",
    nrow(info$loci), " loci. File: ", paths[1])
  invisible(list(files = paths, samples = samples, populations = populations,
    loci = info$loci, alleles = info$alleles))
}
