#' Write rubias
#'
#' @description Create the two-column-per-locus rubias tibble directly from active
#' diploid GDS genotypes, optionally writing it and mapping tables to disk.
#' @param data An open SeqArray GDS object or a GDS filename.
#' @param strata Metadata table or TSV. Use INDIVIDUALS and STRATA for an
#' all-reference dataset, or INDIVIDUALS, SAMPLE_TYPE, REPUNIT and COLLECTION.
#' Lowercase rubias names (indiv, sample_type, repunit, collection) also work.
#' NULL uses GDS metadata. Extra metadata samples are ignored; every active sample
#' must be covered exactly once.
#' @param filename Optional output basename; NULL returns the table without files.
#' @param verbose Display progress messages.
#' @param path.folder Output directory.
#' @param chunk.size Maximum individuals read per genotype block.
#' @param overwrite Allow replacement of existing output files.
#' @return A rubias tibble. Attributes loci, alleles and files contain the mapping
#' tables and any output paths. Genetic columns start at column 5.
#' @details The first columns are sample_type, repunit, collection and indiv.
#' Reference collections must belong to exactly one reporting unit.
#' Mixture repunit values must be NA, and collection identifies the mixture group.
#' Incomplete explicit rubias metadata is rejected rather than replaced by
#' all-reference defaults. GDS selections and whitelists are respected and
#' selections restored after success or failure.
#'
#' Adjacent allele columns use safe L1/L1.2 names with an original-marker map.
#' Allele labels start at one within each locus, including multiallelic loci.
#' Missing diploid calls have two NA values. Partial calls are rejected.
#' This exporter supports diploids only, not rubias's special haploid convention.
#' Block reads reduce temporary memory, but the returned wide tibble must fit in
#' memory. Export reference and mixture together before separating their rows to
#' ensure identical locus columns and allele coding.
#' The export does not perform stock identification, marker selection or quality
#' filtering, and does not preserve phase, depth or genotype likelihoods.
#' The legacy unused parallel.core argument has been removed.
#' @references Moran BM, Anderson EC (2019). Bayesian inference from the
#' conditional genetic stock identification model. Canadian Journal of Fisheries
#' and Aquatic Sciences 76:551-560. \doi{10.1139/cjfas-2018-0016}.
#' @seealso \href{https://eriqande.github.io/rubias/articles/rubias-overview.html}{Official rubias input documentation},
#' [write_gsi_sim()]
#' @examples
#' \dontrun{
#' x <- write_rubias("study.gds", strata = "samples.tsv", filename = "study")
#' # rubias functions use gen_start_col = 5.
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_rubias <- function(data, strata = NULL, filename = NULL, verbose = TRUE,
  path.folder = getwd(), chunk.size = 32L, overwrite = FALSE) {
  .export_flag(verbose, "verbose"); .export_flag(overwrite, "overwrite")
  .export_count(chunk.size, "chunk.size", 1)
  start <- tgbase::startup(package = "genometranslator",
    f.name = "write_rubias", verbose = verbose)
  on.exit(tgbase::teardown(start), add = TRUE)
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  gds <- context$gds
  ids <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  meta <- .gsi_metadata(gds, strata)
  aliases <- c(INDIVIDUALS = "indiv", SAMPLE_TYPE = "sample_type",
    REPUNIT = "repunit", COLLECTION = "collection")
  for (old in names(aliases)) {
    new <- aliases[[old]]
    if (old %in% names(meta) && new %in% names(meta))
      stop("Use only one spelling of each rubias metadata field.")
    if (old %in% names(meta)) names(meta)[names(meta) == old] <- new
  }
  if (!"indiv" %in% names(meta)) stop("Metadata requires INDIVIDUALS or indiv.")
  meta$indiv <- as.character(meta$indiv)
  if (anyNA(meta$indiv) || any(!nzchar(trimws(meta$indiv))) ||
      anyDuplicated(meta$indiv)) stop("Metadata individual IDs must be unique and non-missing.")
  rows <- match(ids, meta$indiv)
  if (anyNA(rows)) stop("Metadata must cover every active sample.")
  meta <- meta[rows, , drop = FALSE]
  fields <- c("sample_type", "repunit", "collection")
  if (!any(fields %in% names(meta))) {
    if (!"STRATA" %in% names(meta)) stop("Reference metadata requires STRATA.")
    meta$sample_type <- "reference"
    meta$repunit <- meta$collection <- as.character(meta$STRATA)
  } else if (!all(fields %in% names(meta))) {
    stop("Explicit rubias metadata requires sample_type, repunit and collection.")
  }
  output <- tibble::as_tibble(lapply(meta[c(fields, "indiv")], as.character))
  if (anyNA(output$sample_type) ||
      any(!output$sample_type %in% c("reference", "mixture")))
    stop("sample_type must be reference or mixture.")
  if (anyNA(output$collection) || any(!nzchar(trimws(output$collection))))
    stop("Every sample requires a collection.")
  ref <- output$sample_type == "reference"
  if (any(is.na(output$repunit[ref]) | !nzchar(trimws(output$repunit[ref]))))
    stop("Reference samples require a repunit.")
  if (any(!is.na(output$repunit[!ref]))) stop("Mixture repunit must be NA.")
  mapping <- unique(output[ref, c("collection", "repunit")])
  if (anyDuplicated(mapping$collection))
    stop("Each reference collection must belong to one repunit.")
  if (length(intersect(output$collection[ref], output$collection[!ref])))
    stop("Reference and mixture collection names must be distinct.")
  info <- .gsi_loci(gds)
  columns <- as.vector(rbind(info$loci$EXPORT_ID, paste0(info$loci$EXPORT_ID, ".2")))
  genetics <- matrix(NA_integer_, length(ids), length(columns),
    dimnames = list(NULL, columns))
  for (first in seq.int(1L, length(ids), by = chunk.size)) {
    r <- first:min(length(ids), first + chunk.size - 1L)
    genetics[r, ] <- .gsi_block(gds, ids[r], info$counts)
  }
  output <- dplyr::bind_cols(output, tibble::as_tibble(genetics))
  paths <- character()
  if (!is.null(filename)) {
    paths <- .export_paths(filename, "rubias", path.folder,
      c("_rubias.tsv", "_loci.tsv", "_alleles.tsv"), overwrite)
    staged <- vapply(paths, function(x) tempfile(), "")
    on.exit(unlink(staged), add = TRUE)
    readr::write_tsv(output, staged[1], na = "NA")
    readr::write_tsv(info$loci, staged[2])
    readr::write_tsv(info$alleles, staged[3])
    .export_publish(staged, paths, overwrite)
  }
  attr(output, "loci") <- info$loci
  attr(output, "alleles") <- info$alleles
  attr(output, "files") <- paths
  if (verbose) .export_message("rubias table: ", length(ids), " samples; ",
    nrow(info$loci), " loci. Genetic columns start at column 5.")
  output
}
