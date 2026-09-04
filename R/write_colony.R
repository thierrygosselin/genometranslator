#' Write COLONY
#'
#' @description Export diploid allele-index calls from a SeqArray GDS to a
#' COLONY input file, with sample, locus and allele correspondence tables.
#' This function exports data; it does not run COLONY or infer relationships.
#' @param data Open SeqArray GDS or an existing GDS path. Active selections and
#' genometranslator whitelists are respected and restored after export.
#' @param strata Optional sample metadata table or TSV file with unique
#' INDIVIDUALS. Must cover every active sample; matching is by ID, not row order.
#' @param sample.markers Optional number of loci to sample without replacement.
#' @param pop.select Optional STRATA values to retain before locus sampling.
#' @param allele.freq NULL to let COLONY estimate frequencies, "overall" to
#' estimate them from all exported samples, or STRATA values defining the
#' frequency reference subset. Every allele observed in the export must occur
#' in that subset. Empirical frequencies are not rounded to two decimals.
#' @param inbreeding 0 for no inbreeding, 1 to allow inbreeding.
#' @param mating.sys.males,mating.sys.females 0 for polygamy, 1 for monogamy.
#' These describe mating over the sampled cohorts, not only one breeding season.
#' @param clone 0 to disable clone inference, 1 to enable it.
#' @param run.length 1, 2, 3 or 4 for short, medium, long or very long runs.
#' @param analysis 0 for pairwise-likelihood score, 1 for full likelihood,
#' or 2 for their combined method.
#' @param allelic.dropout,error.rate Scalar allelic-dropout and other typing
#' error probabilities from zero (inclusive) to one (exclusive).
#' Applied to every exported locus.
#' @param random.seed Positive integer used for locus sampling and written to
#' COLONY. Default 1234. Sampling restores the caller's R random-number state.
#' @param role.column Optional metadata column containing offspring, father or
#' mother (case-insensitive; whitespace is trimmed). Each sample has one role.
#' NULL treats all samples as offspring. At least one offspring is required.
#' @param probability.father,probability.mother Candidate-parent inclusion
#' probabilities from zero to one, inclusive. Required when candidates are present;
#' otherwise zero. They are not estimated from the number of candidates.
#' @param num.runs Number of replicate COLONY runs.
#' @param filename Output basename, optionally ending in .dat.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Allow replacement of existing export files.
#' @param verbose Display messages wrapped to 80 characters.
#' @return Invisibly, named paths to the COLONY .dat file and sample, locus
#' and allele mapping TSV files. COLONY uses S1... sample IDs and L1... locus
#' IDs; interpret results through these tables. Run COLONY in a separate folder
#' for each dataset: its output prefix is the fixed safe name colony_result.
#' @details
#' \itemize{
#' \item Exports codominant, diploid SNP or multiallelic calls. Allele codes
#' are positive integers local to each locus. Fully missing calls become 0 0;
#' partially missing calls, symbolic alleles and presence/absence data are rejected.
#' \item Uses dioecious diploid settings, full-sibship scaling, no sibship
#' prior, no known/excluded relationships, and very-high likelihood precision.
#' Haplodiploid, dominant-marker and overlapping sample-role exports are not
#' supported. Candidate parents are optional, never guessed from sample names.
#' \item Rejects loci without at least two observed alleles across exported
#' samples. No imputation, LD pruning or biological null-allele calling occurs.
#' Assess population structure, batch effects, linkage and typing errors before
#' inference. Multiple SNPs from one RAD locus are not independent evidence.
#' \item Depth, genotype likelihoods and phase are not represented in this
#' hard-call format. Preserve the original GDS. Error-rate parameters are
#' locus-level model settings, not per-genotype GL or PL values.
#' \item Reads the GDS in sample blocks, but retains the exported allele matrix
#' in RAM. This is not an out-of-core whole-genome inference workflow.
#' \item Files are staged before publication; existing files are protected
#' unless overwrite is TRUE. The GDS is never modified.
#' }
#' @references
#' Jones, O. R., & Wang, J. (2010). COLONY: a program for parentage
#' and sibship inference from multilocus genotype data. Molecular Ecology
#' Resources, 10, 551–555. \doi{10.1111/j.1755-0998.2009.02787.x}.
#'
#' Wang, J. (2012). Computationally efficient sibship and parentage assignment
#' from multilocus marker data. Genetics, 191, 183–194.
#' \doi{10.1534/genetics.111.138149}.
#' @seealso \url{https://www.zsl.org/about-zsl/resources/software/colony}
#' for the official software and user guide, section 4.1 (input specification).
#' @examples
#' \dontrun{
#' files <- genometranslator::write_colony("study.gds",
#'   filename = "study", path.folder = "colony_export", random.seed = 1234)
#' # For candidate parents, metadata must assign each sample one explicit role:
#' files <- genometranslator::write_colony("study.gds", strata = "samples.tsv",
#'   role.column = "ROLE", probability.father = 0.5,
#'   probability.mother = 0.5, filename = "parentage")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_colony <- function(data, strata = NULL, sample.markers = NULL,
  pop.select = NULL, allele.freq = NULL, inbreeding = 0,
  mating.sys.males = 0, mating.sys.females = 0, clone = 0,
  run.length = 2, analysis = 1, allelic.dropout = 0, error.rate = 0.02,
  random.seed = 1234L, role.column = NULL, probability.father = NULL,
  probability.mother = NULL, num.runs = 1L, filename = NULL,
  path.folder = getwd(), chunk.size = 32L, overwrite = FALSE, verbose = TRUE) {
  .export_flag(verbose, "verbose")
  .export_flag(overwrite, "overwrite")
  start <- proc.time()[[3]]
  if (verbose) .export_message("genometranslator::write_colony")
  settings <- list(inbreeding = inbreeding, mating.sys.males = mating.sys.males,
    mating.sys.females = mating.sys.females, clone = clone,
    run.length = run.length, analysis = analysis, num.runs = num.runs,
    random.seed = random.seed)
  for (nm in names(settings)) .export_count(settings[[nm]], nm,
    if (nm %in% c("run.length", "num.runs", "random.seed")) 1 else 0)
  if (any(unlist(settings[1:4]) > 1) || run.length > 4 || analysis > 2 ||
      num.runs > .Machine$integer.max || random.seed > .Machine$integer.max)
    rlang::abort("COLONY integer settings are outside their supported ranges.")
  rate <- function(x, name, inclusive = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 ||
        if (inclusive) x > 1 else x >= 1)
      rlang::abort(paste0(name, " must be a scalar probability in ",
        if (inclusive) "[0, 1]." else "[0, 1)."))
    x
  }
  rate(allelic.dropout, "allelic.dropout"); rate(error.rate, "error.rate")
  if (!is.null(allele.freq) && (!is.character(allele.freq) ||
      !length(allele.freq) || anyNA(allele.freq) || any(!nzchar(allele.freq))))
    rlang::abort("allele.freq must be NULL, overall, or STRATA values.")
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  g <- context$gds
  ids <- as.character(SeqArray::seqGetData(g, "sample.id"))
  meta <- .gsi_metadata(g, strata)
  if (anyDuplicated(names(meta)) || !"INDIVIDUALS" %in% names(meta) ||
      anyNA(meta$INDIVIDUALS) || anyDuplicated(meta$INDIVIDUALS) ||
      any(!nzchar(trimws(as.character(meta$INDIVIDUALS)))))
    rlang::abort("Metadata must have unique, non-missing INDIVIDUALS.")
  rows <- match(ids, as.character(meta$INDIVIDUALS))
  if (anyNA(rows)) rlang::abort("Metadata must cover every active sample.")
  meta <- meta[rows, , drop = FALSE]
  if (!is.null(pop.select)) {
    if (!is.character(pop.select) || !length(pop.select) || anyNA(pop.select) ||
        !"STRATA" %in% names(meta) || !all(pop.select %in% meta$STRATA))
      rlang::abort("pop.select must name observed STRATA values.")
    keep <- meta$STRATA %in% pop.select
    meta <- meta[keep, , drop = FALSE]; ids <- ids[keep]
    SeqArray::seqSetFilter(g, sample.id = ids, verbose = FALSE)
  }
  if (!is.null(sample.markers)) {
    .export_count(sample.markers, "sample.markers", 1)
    variants <- SeqArray::seqGetData(g, "variant.id")
    if (sample.markers > length(variants))
      rlang::abort("sample.markers exceeds the number of active loci.")
    old.seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    restore <- function() {
      if (is.null(old.seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          rm(".Random.seed", envir = .GlobalEnv)
      } else assign(".Random.seed", old.seed, envir = .GlobalEnv)
    }
    chosen <- tryCatch({set.seed(random.seed)
      variants[sample.int(length(variants), sample.markers)]
    }, finally = restore())
    SeqArray::seqSetFilter(g, variant.id = chosen, verbose = FALSE)
  }
  roles <- rep("offspring", length(ids))
  if (!is.null(role.column)) {
    if (!is.character(role.column) || length(role.column) != 1L ||
        is.na(role.column) || !role.column %in% names(meta))
      rlang::abort("role.column must name a metadata column.")
    roles <- tolower(trimws(as.character(meta[[role.column]])))
    if (anyNA(roles) || any(!roles %in% c("offspring", "father", "mother")))
      rlang::abort("Every role must be offspring, father or mother.")
  }
  if (!any(roles == "offspring")) rlang::abort("At least one offspring is required.")
  inclusion <- function(p, role) {
    if (any(roles == role)) rate(p, paste0("probability.", role), TRUE)
    else {
      if (!is.null(p) && rate(p, paste0("probability.", role), TRUE) != 0)
        rlang::abort(paste0("No ", role, " candidates: inclusion must be zero."))
      0
    }
  }
  probs <- c(inclusion(probability.father, "father"),
    inclusion(probability.mother, "mother"))
  x <- .legacy_snapshot(g, chunk.size = chunk.size)
  x$samples <- tibble::tibble(INDIVIDUALS = ids,
    EXPORT_ID = paste0("S", seq_along(ids)), ROLE = roles)
  observed <- lapply(seq_along(x$counts), function(j)
    sort(unique(as.vector(x$a[, (2*j-1):(2*j), drop = FALSE]))))
  if (any(lengths(observed) < 2L))
    rlang::abort("Exported loci must have at least two observed alleles.")
  freqs <- NULL
  if (!is.null(allele.freq)) {
    ref <- rep(TRUE, length(ids))
    if (!identical(allele.freq, "overall")) {
      if (!"STRATA" %in% names(meta) || !all(allele.freq %in% meta$STRATA))
        rlang::abort("Frequency reference values must occur in STRATA.")
      ref <- meta$STRATA %in% allele.freq
    }
    freqs <- lapply(seq_along(x$counts), function(j) {
      a <- as.vector(x$a[ref, (2*j-1):(2*j), drop = FALSE])
      n <- vapply(observed[[j]], function(k) sum(a == k, na.rm = TRUE), 0)
      if (any(n == 0)) rlang::abort(paste(
        "Frequency reference lacks an allele observed in exported locus",
        x$loci$EXPORT_ID[j]))
      n/sum(n)
    })
  }
  x$a[is.na(x$a)] <- 0L
  if (!is.null(filename)) filename <- sub("\\.dat$", "", filename,
    ignore.case = TRUE)
  paths <- .legacy_publish(x, filename, path.folder, ".dat", overwrite,
    function(path) {
      collapse <- function(x) paste(x, collapse = " ")
      lines <- c("colony_project", "colony_result", sum(roles == "offspring"),
        nrow(x$loci), random.seed, 0, 2, inbreeding, 0,
        collapse(c(mating.sys.males, mating.sys.females)), clone, 1, 0,
        as.integer(!is.null(freqs)))
      if (!is.null(freqs)) {
        lines <- c(lines, collapse(lengths(observed)))
        for (j in seq_along(freqs)) lines <- c(lines, collapse(observed[[j]]),
          collapse(sprintf("%.17g", freqs[[j]])))
      }
      lines <- c(lines, num.runs, run.length, 0, 10000, 0, analysis, 3,
        collapse(x$loci$EXPORT_ID), "0@",
        paste0(sprintf("%.17g", allelic.dropout), "@"),
        paste0(sprintf("%.17g", error.rate), "@"))
      calls <- function(role) vapply(which(roles == role), function(i)
        collapse(c(x$samples$EXPORT_ID[i], x$a[i, ])), "")
      lines <- c(lines, calls("offspring"), collapse(probs),
        collapse(c(sum(roles == "father"), sum(roles == "mother"))),
        calls("father"), calls("mother"), "0 0", "0 0", rep("0", 6))
      writeLines(lines, path, useBytes = TRUE)
    })
  names(paths) <- c("data", "samples", "loci", "alleles")
  if (verbose) {
    .export_message("Exported ", length(ids), " samples and ", nrow(x$loci),
      " loci. COLONY input: ", paths[1])
    .export_message("Completed write_colony in ",
      round(proc.time()[[3]] - start, 2), " seconds.")
  }
  invisible(paths)
}
