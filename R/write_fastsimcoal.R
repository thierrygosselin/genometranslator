#' Write a fastsimcoal SFS
#'
#' @description Build a joint site-frequency spectrum from active diploid GDS
#' genotypes and write the fastsimcoal multidimensional observed-SFS format.
#' @param data A GDS filepath or open `SeqVarGDSClass` connection.
#' @param strata Sample metadata data frame or TSV; NULL uses GDS metadata.
#' @param group.column Population column in sample metadata.
#' @param pop.levels Complete population order, matching the fastsimcoal model.
#' Default is first appearance among active GDS samples.
#' @param folded Write a minor-allele spectrum (`_MSFS.obs`). FALSE writes a
#' derived-allele spectrum (`_DSFS.obs`) and requires ancestral alleles.
#' @param ancestral.alleles Data frame with unique `VARIANT_ID` and `ANCESTRAL`
#' columns. Ancestral bases must match REF or ALT; unknown/mismatched states
#' are excluded and audited. REF is never assumed ancestral automatically.
#' @param invariant.sites Number of additional callable invariant sites absent
#' from the input records, added to the all-zero bin. For unfolded spectra,
#' these must be ancestral in all exported samples. Never include sites already
#' represented in the GDS, filtered SNPs, or uncallable bases in this count.
#' @param filename Output basename without the `_MSFS.obs`/`_DSFS.obs` suffix.
#' @param path.folder Output directory.
#' @param chunk.size Number of variants read per genotype block.
#' @param max.cells Maximum number of joint-SFS cells, guarding memory use.
#' @param overwrite Replace existing output files.
#' @param verbose Display messages wrapped to 80 columns.
#' @details Exports one full joint spectrum, including for one or two groups,
#' for use with fastsimcoal's `--multiSFS` option and `-m` (folded) or `-d`
#' (unfolded). Population indices are zero-based; sample sizes are chromosome
#' counts, not individual counts. The last population varies fastest in the
#' flattened spectrum, as required by the format specification.
#'
#' Biallelic A/C/G/T SNPs and invariant A/C/G/T records with complete calls in every
#' selected sample contribute. No imputation, likelihood-based SFS estimation,
#' or projection is performed. Excluded sites and reasons are written to an
#' audit table. Complete-case exclusion can bias an SFS when missingness is
#' nonrandom; inspect the audit before demographic inference. MAC/MAF filtering,
#' discovery panels, and RAD locus ascertainment also alter the SFS.
#'
#' Folding uses the allele count pooled over all selected populations, not
#' independent per-population folding. A pooled frequency tie contributes 0.5
#' to each complementary configuration. Observed invariant records remain in
#' the spectrum. Invariant counts are never inferred from chromosome lengths.
#' Invariant records with ALT missing are counted automatically; their calls
#' must be reference homozygotes. Unfolded spectra require an ancestral base
#' for these records too. An ancestral base different from REF places an
#' invariant record in the fixed-derived bin. Depth/quality thresholds are not
#' applied here: a called genotype is not independent proof of callability.
#' Blocked gVCF is not supported (future work). Duplicate mapped positions
#' are excluded, not counted repeatedly. The audit distinguishes input record
#' type and whether retained calls are invariant among the selected samples.
#' A SNP-only spectrum without a callable-site denominator requires an
#' appropriate polymorphic-site-conditioned inference strategy; it is not a
#' complete sequence SFS for absolute demographic scaling.
#'
#' The active selection and metadata whitelists are respected. Selections are
#' restored even on failure; filepath-owned connections are closed. GDS data
#' and metadata are not changed. Diploid genotype calls are required, not
#' SilicoDArT presence/absence values.
#' @return Invisibly, a list with output files, population mapping, a flattened
#' spectrum, dimensions, per-record audit, mask, and export parameters.
#' @references Official fastsimcoal manual, sections on minor-allele and
#' multidimensional site-frequency spectra:
#' \url{https://cmpg.unibe.ch/software/fastsimcoal2/man/fastsimcoal28.pdf}
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
#' @examples
#' \dontrun{
#' sfs <- genometranslator::write_fastsimcoal(
#'   "study.gds", strata = "samples.tsv", pop.levels = c("North", "South"),
#'   filename = "study", path.folder = "exports", folded = TRUE
#' )
#' sfs$populations
#' # Supply a validated VARIANT_ID/ANCESTRAL table for an unfolded spectrum.
#' }
write_fastsimcoal <- function(
    data, strata = NULL, group.column = "STRATA", pop.levels = NULL,
    folded = TRUE, ancestral.alleles = NULL, invariant.sites = 0,
    filename = NULL, path.folder = getwd(), chunk.size = 1000L,
    max.cells = 1e7, overwrite = FALSE, verbose = TRUE
) {
  .write_sfs(data, strata, group.column, pop.levels, folded, ancestral.alleles,
    invariant.sites, filename, path.folder, chunk.size, max.cells, overwrite,
    verbose, format = "fastsimcoal")
}

# Common non-mutating SFS engine. Format-specific serialization is below.
.write_sfs <- function(
    data, strata, group.column, pop.levels, folded, ancestral.alleles,
    invariant.sites, filename, path.folder, chunk.size, max.cells, overwrite,
    verbose, format, mask.corners = TRUE
) {
  force(data)
  .export_flag(folded, "folded")
  .export_flag(overwrite, "overwrite")
  .export_flag(verbose, "verbose")
  .export_count(chunk.size, "chunk.size", 1)
  .export_count(max.cells, "max.cells", 1)
  .export_count(invariant.sites, "invariant.sites", 0)
  start <- tgbase::startup(package = "genometranslator",
    f.name = paste0("write_", format), verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  if (!folded && (!is.data.frame(ancestral.alleles) ||
      !all(c("VARIANT_ID", "ANCESTRAL") %in% names(ancestral.alleles))))
    rlang::abort("Unfolded spectra require a VARIANT_ID/ANCESTRAL table.")
  if (!folded && (anyNA(ancestral.alleles$VARIANT_ID) ||
      anyDuplicated(ancestral.alleles$VARIANT_ID)))
    rlang::abort("Ancestral VARIANT_ID values must be unique and non-missing.")
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  gds <- context$gds
  metadata <- .export_groups(gds, strata, group.column, pop.levels)
  populations <- metadata$populations
  populations$N_CHROMOSOMES <- 2L * populations$N_INDIVIDUALS
  sizes <- populations$N_CHROMOSOMES
  dims <- sizes + 1
  if (prod(dims) > max.cells)
    rlang::abort("Joint SFS exceeds max.cells; reduce the sampled data first.")
  spectrum <- numeric(prod(dims))
  weights <- vapply(seq_along(dims), function(i) {
    if (i == length(dims)) 1 else prod(dims[(i + 1L):length(dims)])
  }, numeric(1))
  ids <- SeqArray::seqGetData(gds, "variant.id")
  definitions <- strsplit(SeqArray::seqGetData(gds, "allele"), ",", fixed = TRUE)
  types <- .vcf_site_types(SeqArray::seqGetData(gds, "allele"))
  if (any(types == "GVCF_REFERENCE"))
    rlang::abort("Blocked gVCF is not supported; use site-level records.")
  end.node <- gdsfmt::index.gdsn(gds, "annotation/info/END", silent = TRUE)
  if (!is.null(end.node)) {
    ends <- SeqArray::seqGetData(gds, "annotation/info/END")
    if (any(types == "INVARIANT_RECORD" &
        ends > SeqArray::seqGetData(gds, "position"), na.rm = TRUE))
      rlang::abort("Blocked gVCF is not supported; use site-level records.")
  }
  invariant <- types == "INVARIANT_RECORD"
  valid <- vapply(definitions, function(x)
    length(x) == 2L && all(x %in% c("A", "C", "G", "T")) && x[1] != x[2],
    logical(1))
  valid <- valid | invariant
  chrom <- SeqArray::seqGetData(gds, "chromosome")
  pos <- SeqArray::seqGetData(gds, "position")
  key <- paste(chrom, pos, sep = ":")
  duplicate <- (duplicated(key) | duplicated(key, fromLast = TRUE)) &
    !is.na(chrom) & !chrom %in% c("", "DENOVO", "CHROM_1") &
    !is.na(pos) & pos > 0
  ancestral <- if (!folded) as.character(ancestral.alleles$ANCESTRAL[
    match(ids, ancestral.alleles$VARIANT_ID)]) else rep(NA_character_, length(ids))
  status <- ifelse(valid, "KEPT", "NOT_BIALLELIC_SNP")
  status[duplicate] <- "DUPLICATE_POSITION"
  valid[duplicate] <- FALSE
  observed.invariant <- rep(NA, length(ids))
  samples <- metadata$samples
  groups <- lapply(populations$GROUP, function(x) which(samples$GROUP == x))
  blocks <- split(seq_along(ids), ceiling(seq_along(ids) / chunk.size))
  for (block in blocks) {
    SeqArray::seqSetFilter(gds, sample.id = samples$INDIVIDUALS,
      variant.id = ids[block], verbose = FALSE)
    gt <- .export_genotypes(gds)
    a <- matrix(gt[1L, , ], nrow = nrow(samples))
    b <- matrix(gt[2L, , ], nrow = nrow(samples))
    for (j in seq_along(block)) {
      i <- block[[j]]
      if (!valid[[i]]) next
      if (anyNA(a[, j]) || anyNA(b[, j])) {
        status[[i]] <- "INCOMPLETE_CALLS"
        next
      }
      allowed <- if (invariant[[i]]) 0L else 0:1
      if (any(!a[, j] %in% allowed) || any(!b[, j] %in% allowed)) {
        status[[i]] <- "INVALID_ALLELE_INDEX"
        next
      }
      counts <- vapply(groups, function(rows)
        sum(a[rows, j] + b[rows, j]), numeric(1))
      observed.invariant[[i]] <- sum(counts) %in% c(0, sum(sizes))
      if (!folded) {
        allowed.ancestor <- if (invariant[[i]]) c("A", "C", "G", "T") else definitions[[i]]
        if (is.na(ancestral[[i]]) || !ancestral[[i]] %in% allowed.ancestor) {
          status[[i]] <- "UNKNOWN_ANCESTRAL"
          next
        }
        if (ancestral[[i]] != definitions[[i]][[1L]]) counts <- sizes - counts
      } else if (sum(counts) > sum(sizes) / 2) {
        counts <- sizes - counts
      }
      index <- 1 + sum(counts * weights)
      if (folded && sum(counts) == sum(sizes) / 2) {
        other <- 1 + sum((sizes - counts) * weights)
        spectrum[[index]] <- spectrum[[index]] + 0.5
        spectrum[[other]] <- spectrum[[other]] + 0.5
      } else spectrum[[index]] <- spectrum[[index]] + 1
    }
  }
  if (!any(status == "KEPT")) rlang::abort("No eligible complete site records.")
  spectrum[[1L]] <- spectrum[[1L]] + invariant.sites
  suffix <- if (folded) "_MSFS.obs" else "_DSFS.obs"
  if (format == "dadi") suffix <- ".fs"
  paths <- .export_paths(filename, format, path.folder,
    c(suffix, "_sfs_populations.tsv", "_sfs_audit.tsv", "_sfs_parameters.tsv"),
    overwrite)
  stage <- tempfile("sfs-export-"); dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  staged <- file.path(stage, basename(paths))
  mask <- rep(FALSE, length(spectrum))
  if (format == "dadi") {
    if (mask.corners) mask[c(1L, length(mask))] <- TRUE
    if (folded) {
      totals <- numeric(length(spectrum))
      for (k in seq_along(dims)) totals <- totals +
        ((seq_along(spectrum) - 1) %/% weights[k]) %% dims[k]
      mask <- mask | totals > sum(sizes) / 2
    }
    # Safe IDs avoid quoting ambiguity; original population names are in TSV.
    writeLines(c(paste(c(dims, if (folded) "folded" else "unfolded",
      paste0('"', populations$EXPORT_ID, '"')), collapse = " "),
      paste(format(spectrum, digits = 17, scientific = FALSE, trim = TRUE),
        collapse = " "), paste(as.integer(mask), collapse = " ")), staged[[1L]])
  } else writeLines(c("1 observations. No. of demes and sample sizes are on next line",
    paste(c(length(sizes), sizes), collapse = " "),
    paste(format(spectrum, digits = 17, scientific = FALSE, trim = TRUE),
      collapse = " ")),
    staged[[1L]])
  audit <- tibble::tibble(VARIANT_ID = ids, SITE_TYPE = types,
    OBSERVED_INVARIANT = observed.invariant, STATUS = status)
  parameters <- tibble::tibble(
    FOLDED = folded, ADDITIONAL_INVARIANT_SITES = invariant.sites,
    INPUT_VARIANTS = length(ids), KEPT_VARIANTS = sum(status == "KEPT"),
    INVARIANT_RECORDS_KEPT = sum(invariant & status == "KEPT"),
    OBSERVED_INVARIANT_KEPT = sum(observed.invariant & status == "KEPT", na.rm = TRUE),
    FORMAT = format, MASK_CORNERS = format == "dadi" && mask.corners,
    EXCLUDED_VARIANTS = sum(status != "KEPT"), CHUNK_SIZE = chunk.size,
    MISSING_POLICY = "complete diploid calls in all selected samples",
    POPULATION_ORDER = paste(populations$GROUP, collapse = ";"))
  readr::write_tsv(populations, staged[[2L]])
  readr::write_tsv(audit, staged[[3L]])
  readr::write_tsv(parameters, staged[[4L]])
  .export_publish(staged, paths, overwrite)
  if (verbose) {
    .export_message(format, " SFS: ", sum(status == "KEPT"),
      " records kept; ", sum(status != "KEPT"), " excluded.\nFile: ", paths[[1L]])
    .export_message("Invariant input records counted: ",
      sum(invariant & status == "KEPT"), ". Additional supplied sites: ",
      invariant.sites, ".")
    if (invariant.sites == 0 && !any(invariant & status == "KEPT")) .export_message(
      "No additional invariant sites supplied. A SNP-only input is not a full ",
      "sequence SFS; account for ascertainment and callable sites in inference.")
  }
  invisible(list(files = paths, populations = populations,
    spectrum = spectrum, dimensions = dims, audit = audit, mask = mask,
    parameters = parameters))
}
