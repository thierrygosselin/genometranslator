#' Write VCF or BCF
#'
#' @description Export selected GDS genotypes as VCF text, BGZF-compressed VCF,
#' or binary BCF. The filename extension selects the encoding. Legacy tidy
#' genotype tables remain supported through the existing GT-only translator.
#' @param data An open SeqVarGDSClass, GDS filepath, or legacy tidy genotype table.
#' @param strata Legacy tidy-table option to include population information.
#' Not supported for GDS input; population metadata are not VCF genotype fields.
#' @param filename Output path ending in `.vcf`, `.vcf.gz`, or `.bcf`.
#' A basename without an extension receives `.vcf`. Default is timestamped.
#' @param source Legacy tidy-table source label. For GDS, existing annotations
#' are used; leave NULL.
#' @param empty Write a legacy empty VCF template; default FALSE.
#' @param index NULL (default) creates a CSI index for compressed VCF and BCF.
#' FALSE disables indexing. Plain VCF cannot be indexed.
#' @param overwrite Replace existing output; default FALSE. Existing sidecars
#' must also be handled safely: export refuses to leave an obsolete index.
#' @param bcftools.path Path or command name of bcftools, required for compressed
#' VCF and BCF. These outputs are coordinate-sorted before indexing.
#' @param verbose Print progress messages.
#' @details GDS export uses SeqArray::seqGDS2VCF and preserves available VCF
#' INFO/FORMAT annotations, genotype allele indices and missing calls, including
#' invariant and multiallelic records. Active selections and metadata whitelists
#' are respected and restored. No genotype filtering or imputation is performed.
#' Annotations not retained on import cannot be recovered. Original site-level
#' INFO summaries (such as AF, AC, AN, NS or DP) are not recalculated after sample
#' selection; recalculate these with downstream tools before interpreting them.
#' Package-specific GDS metadata are not automatically VCF header annotations.
#'
#' BCF is a binary encoding, not a different genotype model or blocked gVCF.
#' Temporary output is created before publication; a failed conversion does not
#' replace the destination. Plain VCF preserves GDS record order; compressed
#' outputs are sorted with bcftools. The legacy tidy path only exports fields
#' constructed by the old translator, not GDS depth or likelihood information.
#' @return Invisibly, a named character vector with the output file and, when
#' requested, its CSI index.
#' @references Danecek P, Auton A, Abecasis G et al. (2011). The variant call
#' format and VCFtools. Bioinformatics 27:2156-2158.
#' \url{https://samtools.github.io/hts-specs/VCFv4.3.pdf}
#' \url{https://samtools.github.io/bcftools/bcftools}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' write_vcf(genome, filename = "study.vcf")
#' write_vcf(genome, filename = "study.vcf.gz")
#' write_vcf(genome, filename = "study.bcf")
#' }
write_vcf <- function(data, strata = FALSE, filename = NULL, source = NULL,
                      empty = FALSE, index = NULL, overwrite = FALSE,
                      bcftools.path = "bcftools", verbose = TRUE) {
  .export_flag(empty, "empty")
  .export_flag(overwrite, "overwrite")
  .export_flag(verbose, "verbose")
  if (is.null(filename)) filename <- paste0("genometranslator_",
    format(Sys.time(), "%Y%m%d_%H%M%S"), ".vcf")
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) ||
      !nzchar(filename)) stop("filename must be one output path.")
  compressed <- grepl("\\.(bcf|vcf\\.gz)$", filename, ignore.case = TRUE)
  bcf <- grepl("\\.bcf$", filename, ignore.case = TRUE)
  if (!compressed && !grepl("\\.vcf$", filename, ignore.case = TRUE))
    filename <- paste0(filename, ".vcf")
  if (is.null(index)) index <- compressed
  .export_flag(index, "index")
  if (index && !compressed) stop("Indexing requires .vcf.gz or .bcf output.")
  folder <- dirname(filename)
  if (!dir.exists(folder)) stop("Output directory does not exist.")
  filename <- file.path(normalizePath(folder), basename(filename))
  sidecars <- paste0(filename, c(".csi", ".tbi"))
  if (!overwrite && any(file.exists(c(filename, sidecars))))
    stop("Output exists; choose another filename or set overwrite = TRUE.")
  # Do not silently leave an old alternate index when replacing data.
  if (any(file.exists(sidecars[if (index) 2L else 1:2])))
    stop("An existing index would become stale; choose a new output filename.")
  exe <- NULL
  if (compressed) {
    exe <- Sys.which(bcftools.path)
    if (!nzchar(exe) && file.exists(bcftools.path)) exe <- bcftools.path
    if (!nzchar(exe)) stop("bcftools is required; supply bcftools.path.")
  }
  start <- tgbase::startup(package = "genometranslator", f.name = "write_vcf",
    verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  stage <- tempfile("vcf-export-"); dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  raw <- file.path(stage, "calls.vcf")
  gds.input <- !empty && (inherits(data, "SeqVarGDSClass") ||
    (is.character(data) && length(data) == 1L &&
      grepl("\\.gds(\\.rad)?$", data, ignore.case = TRUE)))
  if (gds.input) {
    if (!identical(strata, FALSE) || !is.null(source))
      stop("strata and source overrides apply only to legacy tidy input.")
    context <- .export_gds(data)
    on.exit(.export_close(context), add = TRUE)
    SeqArray::seqGDS2VCF(context$gds, vcf.fn = raw, verbose = FALSE)
  } else {
    .write_vcf_tidy(data = if (empty) NULL else data, strata = strata,
      filename = file.path(stage, "calls"), source = source, empty = empty)
  }
  staged <- raw
  if (compressed) {
    staged <- file.path(stage, if (bcf) "calls.bcf" else "calls.vcf.gz")
    genometranslator::bcftools_exec(exe,
      c("sort", if (bcf) "-Ob" else "-Oz", "-o", staged,
        "-T", file.path(stage, "sort"), raw), verbose = FALSE)
  }
  paths <- c(file = filename)
  if (index) {
    genometranslator::bcftools_exec(exe, c("index", "--csi", staged),
      verbose = FALSE)
    staged <- c(staged, paste0(staged, ".csi"))
    paths <- c(paths, index = paste0(filename, ".csi"))
  }
  .export_publish(staged, paths, overwrite)
  if (verbose) .export_message("Output written: ", filename)
  invisible(paths)
}
