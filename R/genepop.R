#' Read a Genepop file
#'
#' @description Validate and read diploid Genepop text into a long or wide
#' tibble. Two- and three-digit allele codes may differ between loci, but must
#' be consistent within each locus. Continued individual records are supported.
#' @param data Filepath or one-column data frame containing all lines,
#' including the title.
#' @param strata Optional sample metadata accepted by [read_strata()]. Must
#' cover every individual; its STRATA column replaces numbered population blocks.
#' @param tidy TRUE returns long STRATA/INDIVIDUALS/MARKERS/GT data;
#' FALSE returns one genotype column per locus.
#' @param filename Optional TSV output path.
#' @param verbose Display progress messages.
#' @details Zero denotes a missing allele. Partial missing calls are retained;
#' each allele is padded to three digits in the returned GT string. Numeric
#' allele labels are not nucleotide identities or inferred repeat lengths.
#' Blank/duplicate cleaned individual names receive generated identifiers and
#' an id_conversion attribute. The original title is retained as genepop_title.
#' Population blocks are numbered in file order unless strata is supplied.
#' Locus names are trimmed; embedded spaces are removed, and collisions rejected.
#'
#' Genepop itself also supports haploid loci, but this reader deliberately
#' rejects two- or three-character haploid genotypes rather than inventing
#' diploid calls. Genotypes must have four or six numeric characters.
#' Empty populations, internal blank lines, malformed rows, duplicate loci and
#' inconsistent coding are errors. Continuations cannot cross a Pop separator.
#' Genepop does not retain read depth, likelihoods, phase or genomic coordinates;
#' prefer GDS or VCF when these are needed.
#' @return A long or wide tibble, optionally also saved as TSV.
#' @references Raymond M. and Rousset F. (1995). GENEPOP (version 1.2):
#' population genetics software for exact tests and ecumenicism.
#' Journal of Heredity 86:248-249.
#' Rousset F. (2008). genepop'007: a complete re-implementation of the genepop
#' software for Windows and Linux. Molecular Ecology Resources 8:103-106.
#' \doi{10.1111/j.1471-8286.2007.01931.x}
#' Official input specification:
#' \url{https://f-rousset.r-universe.dev/genepop/doc/the-input-file.html}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' x <- read_genepop("study.gen")
#' wide <- read_genepop("study.gen", tidy = FALSE)
#' }
read_genepop <- function(data, strata = NULL, tidy = TRUE, filename = NULL,
                         verbose = FALSE) {
  .export_flag(tidy, "tidy")
  .export_flag(verbose, "verbose")
  start <- tgbase::startup(package = "genometranslator",
    f.name = "read_genepop", verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  if (is.character(data) && length(data) == 1L && file.exists(data)) {
    lines <- readLines(data, warn = FALSE)
  } else if (is.data.frame(data) && ncol(data) == 1L) {
    lines <- as.character(data[[1L]])
  } else stop("data must be a Genepop filepath or one-column table.")
  if (length(lines) < 3L || anyNA(lines)) stop("Incomplete Genepop file.")
  title <- lines[[1L]]
  lines <- trimws(lines[-1L])
  while (length(lines) && !nzchar(utils::tail(lines, 1L))) lines <- utils::head(lines, -1L)
  if (any(!nzchar(lines))) stop("Internal blank lines are not valid Genepop records.")
  separators <- which(tolower(lines) == "pop")
  if (!length(separators)) stop("No Genepop population separator was found.")
  if (separators[1] == 1L) stop("No locus names precede the first Pop.")
  markers <- trimws(unlist(strsplit(lines[seq_len(separators[1] - 1L)], ",",
    fixed = TRUE)))
  markers <- gsub("[[:space:]]", "", markers)
  if (!length(markers) || any(!nzchar(markers)) || anyDuplicated(markers) ||
      any(markers %in% c("INDIVIDUALS", "STRATA")))
    stop("Locus names must be nonempty, unique, and not reserved metadata names.")
  ids <- character(); records <- character(); groups <- integer()
  boundaries <- c(separators, length(lines) + 1L)
  for (p in seq_along(separators)) {
    first <- boundaries[p] + 1L; last <- boundaries[p + 1L] - 1L
    if (first > last) stop("Empty Genepop population block.")
    current <- 0L
    for (line in lines[first:last]) {
      comma <- regexpr(",", line, fixed = TRUE)[1]
      if (comma > 0L) {
        if (grepl(",", substring(line, comma + 1L), fixed = TRUE))
          stop("Only one individual/genotype separator comma is allowed.")
        ids <- c(ids, trimws(substr(line, 1L, comma - 1L)))
        records <- c(records, trimws(substring(line, comma + 1L)))
        groups <- c(groups, p)
        current <- length(records)
      } else {
        if (!current) stop("A continuation cannot start a population block.")
        records[current] <- paste(records[current], line)
      }
    }
  }
  rows <- strsplit(trimws(records), "[[:space:]]+")
  if (any(lengths(rows) != length(markers)))
    stop("Genepop genotype count does not match the number of loci.")
  gt <- do.call(rbind, rows)
  for (j in seq_along(markers)) {
    gt[, j] <- .genepop_codes(gt[, j])
  }
  cleaned <- genometranslator::clean_ind_names(ids)
  conversion <- NULL
  if (anyNA(cleaned) || any(!nzchar(cleaned)) || anyDuplicated(cleaned)) {
    cleaned <- paste0("genometranslator-individual-", seq_along(ids))
    conversion <- tibble::tibble(INDIVIDUALS = cleaned, BAD_ID = ids)
    warning("Blank or duplicate sample IDs replaced; see id_conversion.", call. = FALSE)
  }
  pop <- as.character(groups)
  if (!is.null(strata)) {
    meta <- genometranslator::read_strata(strata = strata, verbose = FALSE)$strata
    pop <- as.character(meta$STRATA[match(cleaned, meta$INDIVIDUALS)])
    if (anyNA(pop)) stop("Some Genepop individuals are absent from supplied strata.")
  }
  if (tidy) {
    out <- tibble::tibble(STRATA = rep(pop, times = length(markers)),
      INDIVIDUALS = rep(cleaned, times = length(markers)),
      MARKERS = rep(markers, each = length(ids)), GT = as.vector(gt))
  } else {
    colnames(gt) <- markers
    out <- dplyr::bind_cols(tibble::tibble(STRATA = pop, INDIVIDUALS = cleaned),
      tibble::as_tibble(gt))
  }
  attr(out, "genepop_title") <- title
  if (!is.null(conversion)) attr(out, "id_conversion") <- conversion
  if (!is.null(filename)) readr::write_tsv(out, filename)
  if (verbose) .export_message(length(ids), " individuals; ", length(markers),
    " loci; ", length(unique(pop)), " populations.")
  out
}

.genepop_codes <- function(x) {
  x <- as.character(x)
  widths <- unique(nchar(x))
  if (anyNA(x) || any(!grepl("^[0-9]+$", x)))
    stop("Genepop genotypes must contain numeric allele codes.")
  if (length(widths) != 1L || !widths %in% c(4L, 6L))
    stop("Use consistent four- or six-digit diploid codes within each locus.")
  width <- widths / 2L
  paste0(sprintf("%03d", as.integer(substr(x, 1L, width))),
    sprintf("%03d", as.integer(substring(x, width + 1L))))
}

#' Write a Genepop file
#'
#' @description Export diploid GDS calls in blocks, or validated long tidy
#' genotype data, with explicit population order and reversible ID mappings.
#' @param data Open GDS, GDS filepath, or long data frame containing STRATA,
#' INDIVIDUALS, MARKERS and GT (four- or six-digit diploid codes).
#' @param pop.levels Complete population order, not a population filter.
#' @param genepop.header Optional single-line ASCII title.
#' @param markers.line TRUE writes comma-separated loci on one line; FALSE
#' writes one locus per line.
#' @param filename Basename, optionally ending in .gen.
#' @param strata Optional sample metadata for GDS input; NULL uses GDS metadata.
#' @param group.column Population column for GDS metadata.
#' @param path.folder Output directory.
#' @param chunk.size Number of GDS individuals read per block.
#' @param overwrite Replace existing outputs.
#' @param verbose Display messages wrapped to 80 columns.
#' @details GDS selections and metadata whitelists are respected and restored.
#' Missing alleles are 000, including partially missing genotypes. Each GDS
#' allele receives a locus-specific code: REF=001, first ALT=002, etc., up to
#' 999 alleles. These codes are labels, not nucleotide identities or repeat
#' lengths. An allele mapping table retains original strings. Phase, coverage,
#' likelihoods and coordinates are not carried by Genepop; coordinate metadata
#' are supplied separately in the locus map. Invariant records are retained.
#' Haploid/polyploid and SilicoDArT presence/absence data are not supported.
#'
#' Export identifiers are S1, S2, etc. for samples and L1, L2, etc. for loci,
#' avoiding delimiters and reserved words. Sidecars preserve original names.
#' Population blocks follow pop.levels or first appearance. Tidy input must
#' contain exactly one call per sample/locus and consistent sample populations.
#' No global R options are changed. Caller-owned GDS connections remain open.
#' @return Invisibly, a list with files, sample/population/locus mappings.
#' @references Raymond M. and Rousset F. (1995). GENEPOP (version 1.2):
#' population genetics software for exact tests and ecumenicism.
#' Journal of Heredity 86:248-249.
#' Rousset F. (2008). genepop'007: a complete re-implementation of the genepop
#' software for Windows and Linux. Molecular Ecology Resources 8:103-106.
#' \doi{10.1111/j.1471-8286.2007.01931.x}
#' \url{https://f-rousset.r-universe.dev/genepop/doc/the-input-file.html}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' write_genepop("study.gds", strata = "samples.tsv", filename = "study",
#'   path.folder = "exports", pop.levels = c("North", "South"))
#' }
write_genepop <- function(data, pop.levels = NULL, genepop.header = NULL,
    markers.line = TRUE, filename = NULL, strata = NULL, group.column = "STRATA",
    path.folder = getwd(), chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  .export_flag(markers.line, "markers.line")
  .export_flag(overwrite, "overwrite"); .export_flag(verbose, "verbose")
  .export_count(chunk.size, "chunk.size", 1)
  start <- tgbase::startup(package = "genometranslator",
    f.name = "write_genepop", verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  is.gds <- inherits(data, "SeqVarGDSClass") || is.character(data)
  if (is.gds) {
    context <- .export_gds(data)
    on.exit(.export_close(context), add = TRUE)
    gds <- context$gds
    metadata <- .export_groups(gds, strata, group.column, pop.levels)
    samples <- metadata$samples; populations <- metadata$populations
    ids <- SeqArray::seqGetData(gds, "variant.id")
    mm <- genometranslator::extract_markers_metadata(gds)
    loci <- tibble::tibble(EXPORT_ID = paste0("L", seq_along(ids)),
      VARIANT_ID = ids, MARKERS = if ("MARKERS" %in% names(mm))
        as.character(mm$MARKERS[match(ids, mm$VARIANT_ID)]) else as.character(ids),
      CHROM = SeqArray::seqGetData(gds, "chromosome"),
      POS = SeqArray::seqGetData(gds, "position"))
    definitions <- strsplit(SeqArray::seqGetData(gds, "allele"), ",", fixed = TRUE)
    definitions <- lapply(definitions, function(a) a[!a %in% c("", ".")])
    if (any(lengths(definitions) < 1L | lengths(definitions) > 999L))
      stop("Genepop supports between 1 and 999 allele codes per locus.")
    alleles <- dplyr::bind_rows(lapply(seq_along(ids), function(j)
      tibble::tibble(EXPORT_ID = loci$EXPORT_ID[j],
        CODE = sprintf("%03d", seq_along(definitions[[j]])),
        ALLELE = definitions[[j]])))
  } else {
    if (!is.data.frame(data) ||
        !all(c("STRATA", "INDIVIDUALS", "MARKERS", "GT") %in% names(data)))
      stop("Tidy input requires STRATA, INDIVIDUALS, MARKERS and GT.")
    if (!is.null(strata)) stop("Supply population assignments in tidy STRATA.")
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    if (!nrow(data)) stop("Tidy input contains no genotypes.")
    data$GT <- as.character(data$GT)
    for (column in c("STRATA", "INDIVIDUALS", "MARKERS")) {
      data[[column]] <- as.character(data[[column]])
      if (anyNA(data[[column]]) || any(!nzchar(data[[column]])))
        stop("Tidy IDs and population assignments cannot be missing or empty.")
    }
    samples <- unique(data[c("INDIVIDUALS", "STRATA")])
    names(samples)[2] <- "GROUP"
    if (anyDuplicated(samples$INDIVIDUALS))
      stop("Each individual must belong to exactly one population.")
    if (is.null(pop.levels)) pop.levels <- unique(samples$GROUP)
    if (anyNA(pop.levels) || anyDuplicated(pop.levels) ||
        !setequal(pop.levels, samples$GROUP)) stop("pop.levels must include every population.")
    populations <- tibble::tibble(GROUP = pop.levels,
      EXPORT_ID = paste0("P", seq_along(pop.levels)))
    marker.ids <- unique(data$MARKERS)
    if (anyDuplicated(data[c("INDIVIDUALS", "MARKERS")]) ||
        nrow(data) != nrow(samples) * length(marker.ids))
      stop("Tidy input requires exactly one call per individual and locus.")
    loci <- tibble::tibble(EXPORT_ID = paste0("L", seq_along(marker.ids)),
      MARKERS = marker.ids)
    for (marker in marker.ids) {
      rows <- data$MARKERS == marker
      data$GT[rows] <- .genepop_codes(data$GT[rows])
    }
    alleles <- tibble::tibble(EXPORT_ID = character(), CODE = character(),
      ALLELE = character())
  }
  samples$EXPORT_ID <- paste0("S", seq_len(nrow(samples)))
  if (is.null(genepop.header)) genepop.header <- paste("genometranslator Genepop",
    format(Sys.Date(), "%Y-%m-%d"))
  if (length(genepop.header) != 1L || is.na(genepop.header) ||
      !nzchar(genepop.header) || grepl("[\r\n]", genepop.header) ||
      is.na(iconv(genepop.header, to = "ASCII"))) stop("Use a single-line ASCII title.")
  if (!is.null(filename)) filename <- sub("\\.gen$", "", filename, ignore.case = TRUE)
  paths <- .export_paths(filename, "genepop", path.folder,
    c(".gen", "_samples.tsv", "_populations.tsv", "_loci.tsv", "_alleles.tsv"),
    overwrite)
  stage <- tempfile("genepop-"); dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  staged <- file.path(stage, basename(paths))
  con <- file(staged[1], "wt")
  on.exit(if (!is.null(con)) close(con), add = TRUE)
  writeLines(genepop.header, con)
  writeLines(if (markers.line) paste(loci$EXPORT_ID, collapse = ", ") else
    loci$EXPORT_ID, con)
  for (group in populations$GROUP) {
    writeLines("Pop", con)
    rows <- which(samples$GROUP == group)
    blocks <- split(rows, ceiling(seq_along(rows) / chunk.size))
    for (block in blocks) {
      selected <- samples$INDIVIDUALS[block]
      if (is.gds) {
        SeqArray::seqSetFilter(gds, sample.id = selected, variant.id = ids,
          verbose = FALSE)
        gt <- .export_genotypes(gds)
        actual <- SeqArray::seqGetData(gds, "sample.id")
      }
      for (k in seq_along(block)) {
        if (is.gds) {
          calls <- gt[, match(selected[k], actual), , drop = FALSE]
          dim(calls) <- c(2L, length(ids))
          if (any(!is.na(calls) & (calls < 0 |
              calls >= rep(lengths(definitions), each = 2L))))
            stop("GDS genotype allele index exceeds the locus allele dictionary.")
          codes <- matrix(sprintf("%03d", ifelse(is.na(calls), 0L, calls + 1L)),
            nrow = 2L)
          value <- paste0(codes[1, ], codes[2, ])
        } else {
          d <- data[data$INDIVIDUALS == selected[k], ]
          value <- d$GT[match(loci$MARKERS, d$MARKERS)]
        }
        writeLines(paste0(samples$EXPORT_ID[block[k]], ", ",
          paste(value, collapse = " ")), con)
      }
    }
  }
  close(con); con <- NULL
  readr::write_tsv(samples, staged[2]); readr::write_tsv(populations, staged[3])
  readr::write_tsv(loci, staged[4]); readr::write_tsv(alleles, staged[5])
  .export_publish(staged, paths, overwrite)
  if (verbose) .export_message(nrow(samples), " individuals; ", nrow(loci),
    " loci exported. File: ", paths[1])
  invisible(list(files = paths, samples = samples, populations = populations,
    loci = loci))
}
