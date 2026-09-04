#' Write an Arlequin file
#'
#' @description Export active diploid GDS genotypes as an Arlequin project,
#' with sample, population, and allele-code mapping tables.
#' @param data A GDS filepath or open `SeqVarGDSClass` connection.
#' @param pop.levels Complete population order; never a population filter.
#' @param filename Output basename, optionally ending in `.arp`.
#' @param strata Sample metadata data frame or TSV file. NULL uses GDS metadata.
#' Must cover every active sample, with unique `INDIVIDUALS` identifiers.
#' @param group.column Population column in sample metadata.
#' @param path.folder Output directory.
#' @param chunk.size Number of individuals read per block.
#' @param overwrite Replace existing output files.
#' @param verbose Display messages, wrapped to 80 columns.
#' @details Writes diploid `STANDARD` data with unknown gametic phase.
#' Individuals have frequency one and two allele rows. Missing alleles are `?`,
#' matching the profile declaration. Each biallelic or multiallelic record is
#' one locus, with REF coded 1 and successive ALTs coded 2, 3, etc. Original
#' sequences are recorded in the allele map. Phase is not exported. Haploid,
#' polyploid, and presence/absence data are not supported.
#'
#' Exports the intersection of active selections and metadata whitelists.
#' Temporary selections are restored on success or error. Filepaths are opened
#' read-only and closed; caller-owned connections remain open. Safe identifiers
#' replace sample/population names in the project; mapping tables retain the
#' original names. The Structure section puts all populations in one group,
#' without inferring a biological hierarchy.
#' @return Invisibly, a list of output paths, sample/population maps, and
#' exported variant IDs. The input GDS is not modified.
#' @references Excoffier, L., Laval, G., and Schneider, S. (2005).
#' Arlequin ver. 3.0: An integrated software package for population genetics
#' data analysis. Evolutionary Bioinformatics Online, 1, 47-50.
#'
#' Official Arlequin format specification (sections 4 and 5):
#' \url{https://cmpg.unibe.ch/software/arlequin3522/man/Arlequin35.pdf}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' export <- genometranslator::write_arlequin(
#'   "study.gds", strata = "samples.tsv", filename = "study",
#'   path.folder = "exports"
#' )
#' export$files
#' }
write_arlequin <- function(
    data, pop.levels = NULL, filename = NULL, strata = NULL,
    group.column = "STRATA", path.folder = getwd(), chunk.size = 32L,
    overwrite = FALSE, verbose = TRUE
) {
  force(data)
  .export_count(chunk.size, "chunk.size", 1)
  .export_flag(verbose, "verbose")
  .export_flag(overwrite, "overwrite")
  start <- tgbase::startup(package = "genometranslator",
    f.name = "write_arlequin", verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  gds <- context$gds
  metadata <- .export_groups(gds, strata, group.column, pop.levels)
  ids <- SeqArray::seqGetData(gds, "variant.id")
  chrom <- SeqArray::seqGetData(gds, "chromosome")
  position <- SeqArray::seqGetData(gds, "position")
  alleles <- strsplit(SeqArray::seqGetData(gds, "allele"), ",", fixed = TRUE)
  sample.map <- metadata$samples
  sample.map$EXPORT_ID <- paste0("S", seq_len(nrow(sample.map)))
  paths <- .export_paths(filename, "arlequin", path.folder,
    c(".arp", "_samples.tsv", "_populations.tsv", "_alleles.tsv"), overwrite)
  stage <- tempfile("arlequin-export-")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  staged <- file.path(stage, basename(paths))
  con <- file(staged[[1]], "wt")
  on.exit(if (!is.null(con)) close(con), add = TRUE)
  writeLines(c("[Profile]", 'Title="genometranslator Arlequin export"',
    paste0("NbSamples=", nrow(metadata$populations)), "GenotypicData=1",
    "LocusSeparator=WHITESPACE", "GameticPhase=0", "MissingData='?'",
    "DataType=STANDARD", "[Data]", "[[Samples]]"), con)
  for (group in metadata$populations$GROUP) {
    pop <- metadata$populations$EXPORT_ID[metadata$populations$GROUP == group]
    samples <- sample.map[sample.map$GROUP == group, , drop = FALSE]
    writeLines(c(paste0('SampleName="', pop, '"'),
      paste0("SampleSize=", nrow(samples)), "SampleData={"), con)
    blocks <- split(seq_len(nrow(samples)),
      ceiling(seq_len(nrow(samples)) / chunk.size))
    for (block in blocks) {
      selected <- samples$INDIVIDUALS[block]
      SeqArray::seqSetFilter(gds, sample.id = selected, variant.id = ids,
        verbose = FALSE)
      gt <- .export_genotypes(gds)
      actual <- as.character(SeqArray::seqGetData(gds, "sample.id"))
      for (j in seq_along(selected)) {
        a <- gt[, match(selected[[j]], actual), , drop = FALSE]
        dim(a) <- c(2L, length(ids))
        for (copy in 1:2) {
          value <- ifelse(is.na(a[copy, ]), "?", as.character(a[copy, ] + 1L))
          prefix <- if (copy == 1L) {
            paste(samples$EXPORT_ID[block[[j]]], "1")
          } else ""
          writeLines(paste(prefix, paste(value, collapse = " ")), con)
        }
      }
    }
    writeLines("}", con)
  }
  writeLines(c("[[Structure]]", 'StructureName="One group"', "NbGroups=1",
    "Group={", paste0('"', metadata$populations$EXPORT_ID, '"'), "}"), con)
  close(con)
  con <- NULL
  allele.map <- dplyr::bind_rows(lapply(seq_along(ids), function(i) {
    tibble::tibble(LOCUS_INDEX = i, VARIANT_ID = ids[[i]],
      CHROM = chrom[[i]], POS = position[[i]],
      ALLELE_CODE = seq_along(alleles[[i]]), ALLELE = alleles[[i]])
  }))
  readr::write_tsv(sample.map, staged[[2]])
  readr::write_tsv(metadata$populations, staged[[3]])
  readr::write_tsv(allele.map, staged[[4]])
  .export_publish(staged, paths, overwrite)
  if (verbose) .export_message("Arlequin export: ", nrow(sample.map),
    " samples, ", length(ids), " loci in ", nrow(metadata$populations),
    " populations.\nFile: ", paths[[1]])
  invisible(list(files = paths, samples = sample.map,
    populations = metadata$populations, variant.id = ids))
}
