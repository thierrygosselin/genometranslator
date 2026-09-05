#' Write OutFLANK
#'
#' @description Calculate population FST components from diploid biallelic GDS
#' genotypes and write the per-locus table consumed by OutFLANK. Genotypes are
#' read in locus blocks; no full sample-by-locus matrix is constructed.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param strata Metadata table or TSV containing INDIVIDUALS and STRATA.
#' Defaults to the individual metadata stored in GDS.
#' @param pop.select Optional population names to retain.
#' @param filename Output basename, default "genometranslator". All files
#' receive a shared YYYYMMDD@HHMM timestamp before the extension, for example
#' genometranslator_OutFLANK_YYYYMMDD@HHMM.tsv.
#' @param path.folder Output directory.
#' @param chunk.size Maximum loci per genotype read block.
#' @param min.called Minimum fully called individuals per population and locus,
#' at least two. This is an estimability safeguard, not a recommended study size.
#' @param invalid.loci Either "error" (default) or "exclude". Exclusion of
#' unestimable or monomorphic loci must be requested explicitly and is audited.
#' @param verbose Display progress messages.
#' @details
#' \itemize{
#' \item Requires at least two populations, matched by sample ID rather than
#' metadata row order. GDS filters are restored on success and error.
#' \item Existing files are never replaced. A same-minute filename collision
#' stops the export; choose another basename or rerun in a later minute.
#' \item Requires diploid biallelic hard calls. Presence/absence, multiallelic
#' and symbolic alleles, and partial diploid calls are rejected. GL/PL and
#' read depth are not used; no genotype imputation is performed.
#' \item Computes the diploid biallelic Weir--Cockerham (1984) estimator
#' theta and the corresponding no-sampling-correction statistic required by
#' OutFLANK. This is a local implementation of the equations, not a call to
#' MakeDiploidFSTMat.
#' \item Equal population allele frequencies retain valid He and allele
#' frequency. Negative corrected FST and zero uncorrected FST are not clipped.
#' meanAlleleFreq refers to REF, the first allele in the GDS dictionary,
#' matching OutFLANK's first-homozygote convention; REF is not necessarily
#' ancestral. Allele-label reversal leaves FST and He unchanged.
#' \item Each locus must meet min.called in every selected population and be
#' polymorphic across those populations. Requested exclusions are recorded in
#' the locus audit. Unsupported genotype encodings always stop the export.
#' \item LocusName preserves MARKERS; missing or duplicate marker identifiers
#' are rejected. Dictionaries retain GDS variant IDs and original locus order.
#' Genotype memory scales with samples times chunk.size; output tables remain
#' in memory. The population audit contains per-population call counts.
#' }
#' \strong{Which FST is exported?}
#' \itemize{
#' \item \code{FST} is the per-locus Weir--Cockerham (1984) theta estimate:
#' \eqn{a/(a+b+c)}, where a, b and c are the among-population,
#' among-individual-within-population and within-individual components.
#' \code{T1} is a and \code{T2} is a+b+c. The calculation accounts for
#' finite and unequal population sample sizes using fully called individuals
#' at each locus. Negative estimates are possible and are retained.
#' \item \code{FSTNoCorr} is OutFLANK's corresponding statistic without the
#' finite-sample correction, with numerator \code{T1NoCorr} and denominator
#' \code{T2NoCorr}. It is not interchangeable with corrected WC84 theta.
#' OutFLANK fits its trimmed null distribution to this uncorrected statistic;
#' both sets of components are supplied in its expected input schema.
#' \item These estimators are chosen for compatibility with the OutFLANK
#' method, not because WC84 is universally preferable to Hudson's, Nei's or
#' other estimators. Substituting another estimator changes the statistic and
#' requires reconsidering the null-model calibration. WC84 is appropriate
#' here for individual diploid genotype data; this implementation is not a
#' pooled-read or genotype-likelihood estimator.
#' \item A multilocus summary, if needed, is the ratio of summed numerators
#' to summed denominators over the intended locus set, not the arithmetic
#' mean of per-locus FST values. The writer itself exports per-locus values.
#' }
#' This function prepares statistics only: it does not fit the null model,
#' calculate p-values, prune LD or establish selection. Assess population
#' structure, admixture, batch effects, missingness and linkage before fitting.
#' Sample selection and complete-case requirements can change ascertainment.
#' OutFLANK's downstream handling of zero FST, missing flags and numerical
#' fitting must be checked for the installed version; this writer does not
#' repair those routines or depend on proposed upstream patches. In particular,
#' do not silently delete zero-FST loci merely to make a scan run.
#'
#' \strong{Important for users of dartR's gl.outflank():}
#' The reviewed dartR.popgen 1.2.2 source (dated 23 March 2026) has the
#' following limitations; these should not be attributed to every release.
#' \itemize{
#' \item dartR uses embedded OutFLANK routines. Updating the separate
#' OutFLANK package does not automatically update those calculations.
#' \item Reproduced problems include incorrect equal-frequency statistics,
#' failure when calls remain in only one population, and likelihood overflow.
#' \item This writer safeguards input preparation, but cannot repair
#' downstream fitting or establish that the scan is appropriate for the data.
#' }
#' See \code{vignette("converter_comparison", package = "genometranslator")},
#' section \emph{dartR.popgen: OutFLANK wrapper}, for evidence, reproducible
#' examples, identifier and indexing cautions, and the limits of this review.
#' @return A list with files (named TSV paths), FstDataFrame (exported statistics),
#' samples, populations, loci (all selected loci with STATUS and REASON), and
#' population.audit (called and missing counts per locus and population).
#' @references Weir BS, Cockerham CC (1984). Estimating F-statistics for the
#' analysis of population structure. Evolution 38:1358-1370.
#' \doi{10.1111/j.1558-5646.1984.tb05657.x}
#'
#' Whitlock MC, Lotterhos KE (2015). Reliable detection of loci responsible
#' for local adaptation: inference of a null model through trimming the
#' distribution of FST. The American Naturalist 186(S1):S24-S36.
#' \doi{10.1086/682949}
#' @seealso \href{https://github.com/whitlock/OutFLANK}{OutFLANK documentation},
#' \code{\link{write_pcadapt}}, \code{\link{write_bayescan}}
#' @examples
#' \dontrun{
#' exported <- genometranslator::write_outflank("study.gds",
#'   strata = "samples.tsv", filename = "study")
#' # Review exported$loci and exported$population.audit before analysis.
#' # OutFLANK is optional and is not required by the writer.
#' fit <- OutFLANK::OutFLANK(exported$FstDataFrame,
#'   NumberOfSamples = nrow(exported$populations))
#' # Equivalently, reload the statistics TSV without changing its row order:
#' fst <- read.delim(exported$files[["statistics"]], check.names = FALSE,
#'   colClasses = c(LocusName = "character"))
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_outflank <- function(data, strata = NULL, pop.select = NULL,
  filename = NULL, path.folder = getwd(), chunk.size = 1000L,
  min.called = 2L, invalid.loci = c("error", "exclude"),
  verbose = TRUE) {
  start <- .legacy_start("write_outflank", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  .export_count(chunk.size, "chunk.size", 1)
  .export_count(min.called, "min.called", 2)
  invalid.loci <- match.arg(invalid.loci)
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  g <- context$gds
  meta <- .export_groups(g, strata, "STRATA", NULL)
  if (!is.null(pop.select)) {
    if (!is.character(pop.select) || !length(pop.select) || anyNA(pop.select) ||
        anyDuplicated(pop.select) || !all(pop.select %in% meta$samples$GROUP))
      stop("Invalid pop.select.", call. = FALSE)
    meta$samples <- meta$samples[meta$samples$GROUP %in% pop.select, ]
    meta$populations <- meta$populations[meta$populations$GROUP %in% pop.select, ]
  }
  groups <- meta$populations$GROUP
  if (length(groups) < 2L) stop("OutFLANK requires at least two populations.")
  if (any(meta$populations$N_INDIVIDUALS < min.called))
    stop("Each selected population must have at least min.called individuals.")
  SeqArray::seqSetFilter(g, sample.id = meta$samples$INDIVIDUALS, verbose = FALSE)
  info <- .gsi_loci(g)
  if (any(info$counts != 2L)) stop("OutFLANK export requires biallelic loci.")
  loci <- info$loci
  if (anyNA(loci$MARKERS) || any(!nzchar(trimws(loci$MARKERS))) ||
      anyDuplicated(loci$MARKERS)) stop("Marker IDs must be unique and non-missing.")
  ids <- loci$VARIANT_ID
  results <- vector("list", length(ids))
  audit <- vector("list", length(ids))
  reasons <- rep("", length(ids))
  for (first in seq.int(1L, length(ids), by = chunk.size)) {
    rows <- first:min(length(ids), first + chunk.size - 1L)
    SeqArray::seqSetFilter(g, variant.id = ids[rows], verbose = FALSE)
    gt <- .export_genotypes(g)
    if (any(gt > 1L, na.rm = TRUE)) stop("Invalid biallelic allele index.")
    for (k in seq_along(rows)) {
      j <- rows[k]
      a <- matrix(gt[, , k], nrow = 2L)
      if (any(xor(is.na(a[1, ]), is.na(a[2, ]))))
        stop("Partial diploid calls are unsupported; resolve them explicitly.")
      dose <- colSums(a)
      counts <- t(vapply(groups, function(p) {
        z <- dose[meta$samples$GROUP == p]
        tabulate(z[!is.na(z)] + 1L, nbins = 3L)
      }, integer(3)))
      n <- unname(rowSums(counts))
      audit[[j]] <- tibble::tibble(VARIANT_ID = ids[j], MARKERS = loci$MARKERS[j],
        STRATA = groups, N_CALLED = n,
        N_MISSING = meta$populations$N_INDIVIDUALS - n)
      if (any(n < min.called)) reasons[j] <- "insufficient_called_per_population"
      else if (sum(counts[, 2] + 2 * counts[, 3]) %in% c(0, 2 * sum(n)))
        reasons[j] <- "monomorphic"
      else {
        results[[j]] <- .outflank_components(counts)
        if (any(!is.finite(results[[j]]))) reasons[j] <- "undefined_components"
      }
    }
    if (verbose) .export_message("OutFLANK statistics: ", max(rows),
      "/", length(ids), " loci.")
  }
  keep <- reasons == ""
  if (any(!keep) && invalid.loci == "error")
    stop(sum(!keep), " loci failed estimability checks. ",
      "Filter explicitly or use invalid.loci = 'exclude' for an audit.",
      call. = FALSE)
  if (!any(keep)) stop("No estimable polymorphic loci remain; no files written.")
  fst <- as.data.frame(do.call(rbind, results[keep]))
  fst <- cbind(LocusName = loci$MARKERS[keep], fst)
  loci$STATUS <- ifelse(keep, "exported", "excluded")
  loci$REASON <- reasons
  audit <- dplyr::bind_rows(audit)
  suffixes <- paste0(c("_OutFLANK_", "_OutFLANK_samples_",
    "_OutFLANK_populations_", "_OutFLANK_loci_", "_OutFLANK_population_audit_"),
    file.date, ".tsv")
  # Resolve names only; publication below always prohibits replacement.
  paths <- .export_paths(if (is.null(filename)) "genometranslator" else filename,
    "outflank", path.folder, suffixes, TRUE)
  if (any(file.exists(paths)))
    stop("Output exists for this timestamp. Choose another basename or minute.",
      call. = FALSE)
  names(paths) <- c("statistics", "samples", "populations", "loci", "population.audit")
  staged <- vapply(paths, function(p) tempfile(), "")
  on.exit(unlink(staged), add = TRUE)
  tables <- list(fst, meta$samples, meta$populations, loci, audit)
  for (i in seq_along(staged)) readr::write_tsv(tables[[i]], staged[i])
  .export_publish(staged, paths, FALSE)
  if (verbose) .export_message("Exported ", sum(keep), " loci; excluded ", sum(!keep), ".")
  list(files = paths, FstDataFrame = fst, samples = meta$samples,
    populations = meta$populations, loci = loci, population.audit = audit)
}

# Counts are individuals in REF/REF, REF/ALT and ALT/ALT classes.
# WC84 corrected components and OutFLANK's no-sampling-correction convention.
.outflank_components <- function(counts) {
  n <- rowSums(counts)
  r <- length(n)
  nbar <- mean(n)
  nc <- (sum(n) - sum(n^2)/sum(n))/(r - 1)
  p <- (counts[, 1] + counts[, 2]/2)/n
  pbar <- sum(n * p)/sum(n)
  s2 <- sum(n * (p - pbar)^2)/((r - 1) * nbar)
  hbar <- sum(counts[, 2])/sum(n)
  a <- nbar/nc * (s2 - (pbar*(1-pbar) - (r-1)/r*s2 - hbar/4)/(nbar-1))
  b <- nbar/(nbar-1) * (pbar*(1-pbar) - (r-1)/r*s2 - (2*nbar-1)/(4*nbar)*hbar)
  c <- hbar/2
  an <- nbar/nc*s2
  bn <- pbar*(1-pbar) - (r-1)/r*s2 - hbar/2
  c(He = 2*pbar*(1-pbar), FST = a/(a+b+c), T1 = a, T2 = a+b+c,
    FSTNoCorr = an/(an+bn+c), T1NoCorr = an, T2NoCorr = an+bn+c,
    meanAlleleFreq = pbar)
}
