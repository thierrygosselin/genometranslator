#' Write a dadi spectrum
#'
#' @description Build a joint site-frequency spectrum from active diploid GDS
#' calls and write a native dadi `.fs` file. Replaces the legacy tidy-data SNP
#' table writer; the old FASTA and sumstats arguments are no longer used.
#' @inheritParams write_fastsimcoal
#' @param filename Output basename without the `.fs` suffix.
#' @param mask.corners Mask the all-ancestral and all-derived cells (default
#' TRUE), as normally required for dadi inference. Counts remain in the file
#' and audit even when masked. Folded spectra additionally mask cells above
#' half the pooled chromosome count. Set FALSE only for a model that explicitly
#' handles the monomorphic cells.
#' @details Uses the same site selection and global folding as
#' [write_fastsimcoal()]. Complete calls are required in every selected sample;
#' no imputation, projection, or likelihood-based estimation is performed.
#' Eligible invariant records are counted automatically, including in an
#' all-sites GDS. Missing genotypes are never replaced with reference calls.
#' For an unfolded spectrum, supply ancestral alleles for invariant records
#' as well as SNPs. REF is not assumed ancestral.
#'
#' Population IDs in the file are P1, P2, etc.; the population sidecar maps
#' these IDs to the original names in `pop.levels` order. Dimensions are
#' chromosome counts plus one; the last population varies fastest.
#' Load in Python with `dadi.Spectrum.from_file("study.fs")`.
#'
#' Counts of reported calls are not independent evidence of callability.
#' Review missingness, caller thresholds, ascertainment and upstream filtering
#' before inference. MAC/MAF filtering changes the SFS. Invariant sites alone
#' do not establish an unbiased callable-site denominator.
#' Input data and selections are preserved; existing files require overwrite.
#' @section Future work:
#' Blocked gVCF support, population-specific projection, and genotype-likelihood
#' SFS estimation are not implemented. Block-aware import must preserve
#' boundaries and quality summaries without inventing per-base information.
#' @return Invisibly, output paths, population mapping, flattened spectrum,
#' dimensions, mask, per-record audit, and export parameters.
#' @references dadi frequency-spectrum file specification:
#' \url{https://dadi.readthedocs.io/en/latest/user-guide/importing-data/}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' result <- genometranslator::write_dadi(
#'   "study.gds", strata = "samples.tsv", pop.levels = c("North", "South"),
#'   filename = "study", path.folder = "exports"
#' )
#' result$parameters
#' }
write_dadi <- function(
    data, strata = NULL, group.column = "STRATA", pop.levels = NULL,
    folded = TRUE, ancestral.alleles = NULL, invariant.sites = 0,
    filename = NULL, path.folder = getwd(), chunk.size = 1000L,
    max.cells = 1e7, overwrite = FALSE, verbose = TRUE, mask.corners = TRUE
) {
  .export_flag(mask.corners, "mask.corners")
  .write_sfs(data, strata, group.column, pop.levels, folded, ancestral.alleles,
    invariant.sites, filename, path.folder, chunk.size, max.cells, overwrite,
    verbose, format = "dadi", mask.corners = mask.corners)
}
