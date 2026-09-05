#' Write pcadapt
#'
#' @description Prepare diploid biallelic genotypes for a pcadapt genome scan,
#' preserving sample and marker order and writing identifier dictionaries.
#' The supported pcadapt text format contains one locus per row and one sample
#' per column, with ALT-copy counts of 0, 1, or 2 and 9 for missing genotypes.
#' Import this file with \code{pcadapt::read.pcadapt(type = "pcadapt",
#' type.out = "bed")} to create the compact BED input used by pcadapt 4.
#' @param data Open SeqArray GDS or GDS path. Active whitelists are respected.
#' @param pop.select Optional populations to retain explicitly.
#' @param filename Optional output basename; use path.folder for the directory.
#' @param strata Metadata table or TSV with INDIVIDUALS and STRATA, used for
#' population selection and the sample dictionary. Populations are not required
#' by pcadapt itself.
#' @param path.folder Output directory.
#' @param chunk.size Number of samples per GDS read block.
#' @param overwrite Allow replacement of existing output files.
#' @param verbose Display progress messages.
#' @details
#' \itemize{
#' \item Requires diploid, biallelic hard genotype calls. Multiallelic and
#' symbolic allele definitions and partial diploid calls are rejected. GL/PL,
#' phase and allele depth are not represented in this analysis format.
#' \item After population selection, every locus must have called genotypes
#' and both alleles observed; every sample must have at least one called locus.
#' Failures stop the export. No implicit filtering or imputation is performed.
#' \item Complete missing calls become 9 in the text file. Matrix input to
#' pcadapt instead requires NA; do not pass the returned matrix unchanged.
#' \item Active GDS selections are respected and restored, including on errors.
#' Samples and loci in the dictionaries follow the exported matrix order.
#' \item This writer uses the supported text-to-BED import route, not pcadapt's
#' deprecated VCF/PED converters. For large datasets, direct BED export with
#' \code{\link{write_plink}} is another route. This writer assembles matrices
#' in RAM; chunk.size limits temporary reads, not total memory use.
#' }
#' Keep the object returned by read.pcadapt: its generated BED stores dimensions
#' in that object and does not require companion BIM/FAM files for this route.
#' Use a fresh basename for each changed dataset: pcadapt can reuse an existing
#' text-file-plus-.bed output without rebuilding it. Overwriting the text alone
#' does not refresh that derived BED file.
#'
#' Choose K from PCA diagnostics and biological context. Check associations of
#' PCs with missingness, batch, relatedness and large linked regions. pcadapt
#' detects unusual associations with PCs; these are candidates for selection,
#' not proof of adaptation. Examine LD and sensitivity to K and min.maf.
#' Internal LD clumping uses marker-index windows, not physical distances.
#' The supplied pcadapt 4.4.1 source has a scaling defect in the optional
#' componentwise test with multiple PCs; the examples use the default
#' Mahalanobis test, which uses a separate calculation.
#' @return A list containing \code{files} (genotype text, sample, locus and
#' allele dictionary paths, in that order), \code{genotype.matrix} (loci by
#' samples, missing = 9), \code{pop.string} (assignments when available, otherwise
#' NULL), \code{samples}, and \code{loci}.
#' @references Luu, K., Bazin, E., & Blum, M. G. (2017).
#' pcadapt: an R package to perform genome scans for selection based on principal component analysis.
#' Molecular Ecology Resources, 17(1), 67-77.
#' \doi{10.1111/1755-0998.12592}
#'
#' Prive F, Luu K, Vilhjalmsson BJ, Blum MGB (2020). Performing Highly Efficient
#' Genome Scans for Local Adaptation with R Package pcadapt Version 4.
#' Molecular Biology and Evolution 37:2153-2154. \doi{10.1093/molbev/msaa053}
#'
#' Duforet-Frebourg N, Luu K, Laval G, Bazin E, Blum MGB (2016). Detecting Genomic
#' Signatures of Natural Selection with Principal Component Analysis:
#' Application to the 1000 Genomes Data. Molecular Biology and Evolution
#' 33:1082-1093. \doi{10.1093/molbev/msv334}
#' @seealso \href{https://bcm-uga.github.io/pcadapt/reference/read.pcadapt.html}{pcadapt input documentation},
#' \code{\link{write_plink}}
#' @examples
#' \dontrun{
#' exported <- genometranslator::write_pcadapt(
#'   "study.gds", filename = "study", path.folder = "pcadapt_input")
#' input <- pcadapt::read.pcadapt(exported$files[1], type = "pcadapt",
#'   type.out = "bed")
#' # Inspect the scree plot; choose K for the study, not automatically as 2.
#' exploratory <- pcadapt::pcadapt(input, K = 5, pca.only = TRUE)
#' plot(exploratory$singular.values, type = "b", xlab = "PC",
#'   ylab = "Relative singular value")
#' scan <- pcadapt::pcadapt(input, K = 2, method = "mahalanobis", min.maf = 0.05)
#' results <- exported$loci
#' results$pvalue <- scan$pvalues
#' tested <- is.finite(results$pvalue)
#' results$BH_FDR <- NA_real_
#' results$BH_FDR[tested] <- stats::p.adjust(results$pvalue[tested], method = "BH")
#' candidates <- results[tested & !is.na(results$BH_FDR) &
#'   results$BH_FDR <= 0.05, ]
#'
#' # Optional population subset; polymorphism is checked again after selection.
#' subset <- genometranslator::write_pcadapt("study.gds",
#'   strata = "samples.tsv", pop.select = c("North", "South"),
#'   filename = "north_south", path.folder = "pcadapt_input")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_pcadapt <- function(data, pop.select = NULL, filename = NULL,
  strata = NULL, path.folder = getwd(), chunk.size = 32L, overwrite = FALSE,
  verbose = TRUE) {
  start <- .legacy_start("write_pcadapt", verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  .export_flag(overwrite, "overwrite")
  x <- .legacy_snapshot(data, strata,
    grouped = !is.null(pop.select) || !is.null(strata),
    biallelic = TRUE, chunk.size = chunk.size)
  if (!is.null(pop.select)) {
    if (!length(pop.select) || anyNA(pop.select) ||
        !all(pop.select %in% x$samples$GROUP)) stop("Unknown or empty pop.select.")
    keep <- x$samples$GROUP %in% pop.select
    x$samples <- x$samples[keep, ]; x$a <- x$a[keep, , drop = FALSE]
  }
  d <- t(.legacy_dosage(x))
  called <- rowSums(!is.na(d))
  if (any(called == 0L)) stop("pcadapt requires called genotypes at every locus.")
  alt <- rowSums(d, na.rm = TRUE)
  if (any(alt == 0 | alt == 2 * called))
    stop("pcadapt requires polymorphic loci after sample selection.")
  if (any(colSums(!is.na(d)) == 0L))
    stop("pcadapt cannot use samples with all genotypes missing.")
  d[is.na(d)] <- 9L
  paths <- .legacy_publish(x, filename, path.folder, "_pcadapt.txt", overwrite,
    function(p) utils::write.table(d, p, row.names = FALSE,
      col.names = FALSE, quote = FALSE))
  list(files = paths, genotype.matrix = d,
    pop.string = if ("GROUP" %in% names(x$samples)) x$samples$GROUP else NULL,
    samples = x$samples, loci = x$loci)
}
