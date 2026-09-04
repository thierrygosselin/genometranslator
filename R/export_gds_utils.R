# Shared helpers for non-mutating GDS exporters.
.export_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x))
    rlang::abort(paste0("`", name, "` must be TRUE or FALSE."))
}
.export_count <- function(x, name, minimum = 0) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != floor(x))
    rlang::abort(paste0("`", name, "` must be a whole number >= ", minimum, "."))
}
.export_message <- function(...) {
  lines <- unlist(lapply(strsplit(paste0(...), "\n", fixed = TRUE)[[1]],
    function(x) strwrap(x, width = 80)))
  lines <- unlist(lapply(lines, function(x) {
    starts <- seq.int(1L, max(1L, nchar(x)), by = 80L)
    substring(x, starts, starts + 79L)
  }))
  message(paste(lines, collapse = "\n"))
}
.export_gds <- function(data) {
  owned <- !inherits(data, "SeqVarGDSClass")
  if (owned && (!is.character(data) || length(data) != 1L ||
      is.na(data) || !file.exists(data)))
    rlang::abort("`data` must be an open SeqArray GDS or an existing GDS path.")
  gds <- if (owned) SeqArray::seqOpen(data, readonly = TRUE) else data
  SeqArray::seqFilterPush(gds)
  success <- FALSE
  on.exit(if (!success) {
    SeqArray::seqFilterPop(gds)
    if (owned) SeqArray::seqClose(gds)
  }, add = TRUE)
  source <- genometranslator::extract_data_source(gds)
  if (any(grepl("silico", source, ignore.case = TRUE)))
    rlang::abort("Presence/absence observations are not diploid genotypes.")
  samples <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  variants <- SeqArray::seqGetData(gds, "variant.id")
  im <- genometranslator::extract_individuals_metadata(gds, whitelist = FALSE)
  mm <- genometranslator::extract_markers_metadata(gds, whitelist = FALSE)
  samples <- samples[samples %in% im$INDIVIDUALS[
    is.na(im$FILTERS) | im$FILTERS == "whitelist"]]
  variants <- variants[variants %in% mm$VARIANT_ID[
    is.na(mm$FILTERS) | mm$FILTERS == "whitelist"]]
  if (!length(samples) || !length(variants))
    rlang::abort("No active whitelisted samples or variants to export.")
  SeqArray::seqSetFilter(gds, sample.id = samples, variant.id = variants,
    verbose = FALSE)
  success <- TRUE
  list(gds = gds, owned = owned)
}
.export_close <- function(context) {
  SeqArray::seqFilterPop(context$gds)
  if (context$owned) SeqArray::seqClose(context$gds)
}
.export_groups <- function(gds, strata, column, order) {
  if (!is.character(column) || length(column) != 1L || is.na(column))
    rlang::abort("`group.column` must name one metadata column.")
  if (is.null(strata)) {
    strata <- genometranslator::extract_individuals_metadata(gds)
  } else if (is.character(strata) && length(strata) == 1L && file.exists(strata)) {
    strata <- readr::read_tsv(strata, col_types = readr::cols(.default = "c"),
      progress = FALSE)
  }
  if (!is.data.frame(strata) || !all(c("INDIVIDUALS", column) %in% names(strata)))
    rlang::abort("Metadata must contain INDIVIDUALS and the population column.")
  ids <- as.character(strata$INDIVIDUALS)
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids))
    rlang::abort("Metadata sample IDs must be unique and non-missing.")
  selected <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  rows <- match(selected, ids)
  if (anyNA(rows)) rlang::abort("Metadata must cover every active sample.")
  groups <- as.character(strata[[column]][rows])
  if (anyNA(groups) || any(!nzchar(trimws(groups))))
    rlang::abort("Every exported sample must have a population assignment.")
  if (is.null(order)) order <- unique(groups)
  if (!is.character(order) || anyNA(order) || anyDuplicated(order) ||
      !setequal(order, groups))
    rlang::abort("`pop.levels` must list every observed population exactly once.")
  list(samples = tibble::tibble(INDIVIDUALS = selected, GROUP = groups),
    populations = tibble::tibble(GROUP = order,
      EXPORT_ID = paste0("P", seq_along(order)),
      POP_INDEX = seq_along(order) - 1L,
      N_INDIVIDUALS = vapply(order, function(x) sum(groups == x), integer(1))))
}
.export_genotypes <- function(gds) {
  gt <- SeqArray::seqGetData(gds, "genotype")
  if (length(dim(gt)) != 3L || dim(gt)[[1L]] != 2L)
    rlang::abort("Export requires diploid allele-index genotypes.")
  gt[!is.na(gt) & gt < 0] <- NA_integer_
  gt
}
.export_paths <- function(filename, label, folder, suffixes, overwrite) {
  if (is.null(filename)) filename <- paste0("genometranslator_", label, "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"))
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) ||
      !nzchar(filename) || basename(filename) != filename)
    rlang::abort("`filename` must be a basename; use `path.folder` for its path.")
  filename <- sub("\\.arp$", "", filename, ignore.case = TRUE)
  if (!dir.exists(folder)) dir.create(folder, recursive = TRUE)
  paths <- file.path(normalizePath(folder, mustWork = TRUE),
    paste0(filename, suffixes))
  if (!overwrite && any(file.exists(paths)))
    rlang::abort("Output exists. Choose a new filename or set overwrite = TRUE.")
  paths
}
.export_publish <- function(staged, paths, overwrite) {
  if (!all(file.copy(staged, paths, overwrite = overwrite)))
    rlang::abort("Could not publish all export files; inspect the output folder.")
}
