#' @title Clean marker's names for genomic workflows
#' @description Function to clean marker's name
#' of weird separators
#' that interfere with some packages
#' or codes. \code{/}, \code{:}, \code{-} and \code{.} are changed to an underscore
#' \code{_}.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.
#' @param x (character string) Markers character string.
#' @rdname clean_markers_names
#' @keywords internal
#' @export
clean_markers_names <- function(x) {
  x <- stringi::stri_replace_all_fixed(
    str = as.character(x),
    pattern = c("/", ":", "-", "."),
    replacement = "_",
    vectorize_all = FALSE)
}#End clean_markers_names

# Safe chromosome component used inside MARKERS. CHROM itself remains the
# authoritative, unmodified reference/assembly sequence identifier.
safe_marker_chrom <- function(chrom) {
  original <- as.character(chrom)
  de.novo <- is.na(original) | !nzchar(original) |
    original %in% c("CHROM", "DENOVO")
  safe <- clean_markers_names(original)
  safe[de.novo] <- "DENOVO"
  safe <- gsub("[^A-Za-z0-9_]", "_", safe)

  map <- unique(data.frame(original = original, safe = safe, stringsAsFactors = FALSE))
  map <- map[!is.na(map$original) & nzchar(map$original), , drop = FALSE]
  collisions <- map$safe[duplicated(map$safe) | duplicated(map$safe, fromLast = TRUE)]
  if (length(collisions)) {
    details <- vapply(unique(collisions), function(id) {
      paste(map$original[map$safe == id], collapse = ", ")
    }, character(1))
    rlang::abort(paste0(
      "Distinct CHROM values produce the same safe MARKERS component: ",
      paste(paste0(unique(collisions), " <- ", details), collapse = "; "),
      ". Rename the conflicting reference sequences explicitly."
    ))
  }
  safe
}

make_marker_id <- function(chrom, locus, position) {
  stringi::stri_join(
    safe_marker_chrom(chrom),
    clean_markers_names(locus),
    clean_markers_names(position),
    sep = "__"
  )
}

#' @title Clean individual's names for genomic workflows
#' @description Function to clean individual's name
#' that interfere with some packages
#' or codes. \code{_} and \code{:} are changed to a dash \code{-}.
#' Whitespaces are removed.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.
#' @param x (character string) Individuals character string.
#' @rdname clean_ind_names
#' @export
clean_ind_names <- function(x) {
  x <- stringi::stri_replace_all_fixed(
    str = as.character(x),
    pattern = c("_", ":", " ", ","),
    replacement = c("-", "-", "", ""),
    vectorize_all = FALSE) %>%
    stringi::stri_replace_all_regex(
      str = .,
      pattern = "\\s+",
      replacement = "",
      vectorize_all = FALSE
    )
}#End clean_ind_names

#' @title Clean population's names for genomic workflows
#' @description Function to clean pop's name
#' that interfere with some packages
#' or codes. Space is changed to an underscore \code{_}.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.
#' @param x (character string) Population character string.
#' @param factor (logical) Default: \code{factor = TRUE}. Will also keep or transform the factor levels.
#' @rdname clean_pop_names
#' @export
clean_pop_names <- function(x, factor = TRUE) {
  clean_pop <- function(x) {
    stringi::stri_replace_all_fixed(
      str = as.character(x),
      pattern = " ",
      replacement = "_",
      vectorize_all = FALSE)
  }

  if (is.factor(x)) {
    pop.levels <- clean_pop(as.character(levels(x)))
  } else {
    pop.levels <- clean_pop(as.character(unique(x)))
  }
  x <- clean_pop(as.character(x))
  if (factor) x <- factor(x, levels = pop.levels)
  return(x)
}#End clean_pop_names


# Genomic output filename policy ----------------------------------------------
generate_filename <- function(
    name.shortcut = NULL,
    path.folder = NULL,
    date = TRUE,
    extension = c(
      "tsv", "gds.rad", "rad", "gds", "gen", "dat", "genind",
      "genlight", "gtypes", "vcf", "colony", "bayescan", "gsisim",
      "hierfstat", "hzar", "ldna", "pcadapt", "plink", "related",
      "stockr", "structure", "arlequin", "arrow.parquet"
    )
) {
  extension <- match.arg(extension)
  stem <- if (is.null(name.shortcut)) "genometranslator" else name.shortcut
  timestamp <- if (is.character(date)) date else if (isTRUE(date)) format(Sys.time(), "%Y%m%d@%H%M") else NULL
  dated <- function(x) if (is.null(timestamp)) x else paste0(x, "_", timestamp)

  spec <- switch(
    extension,
    gen = list(stem = paste0(stem, "_genepop"), extension = "gen"),
    dat = list(stem = paste0(stem, "_fstat"), extension = "dat"),
    hierfstat = list(stem = paste0(stem, "_hierfstat"), extension = "dat"),
    structure = list(stem = paste0(stem, "_structure"), extension = "str"),
    genind = list(stem = paste0(stem, "_genind"), extension = "RData"),
    genlight = list(stem = paste0(stem, "_genlight"), extension = "RData"),
    gtypes = list(stem = paste0(stem, "_gtypes"), extension = "RData"),
    stockr = list(stem = paste0(stem, "_stockr"), extension = "RData"),
    bayescan = list(stem = paste0(stem, "_bayescan"), extension = "txt"),
    pcadapt = list(stem = paste0(stem, "_pcadapt"), extension = "txt"),
    related = list(stem = paste0(stem, "_related"), extension = "txt"),
    hzar = list(stem = paste0(stem, "_hzar"), extension = "csv"),
    arlequin = list(stem = paste0(stem, "_arlequin"), extension = "csv"),
    list(stem = stem, extension = extension)
  )

  tgbase::build_filename(
    stem = dated(spec$stem),
    extension = spec$extension,
    path.folder = path.folder,
    date = FALSE,
    create = TRUE
  )
}

# Read the end of a text file without relying on the Unix `tail` executable.
read_last_lines <- function(path, n = 20L, chunk.size = 65536L) {
  n <- max(1L, as.integer(n)[1L])
  con <- file(path, open = "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = 0L, origin = "end")
  size <- seek(con)
  position <- size
  raw <- raw(0)

  while (position > 0L) {
    take <- min(as.double(chunk.size), position)
    position <- position - take
    seek(con, where = position, origin = "start")
    raw <- c(readBin(con, what = "raw", n = take), raw)
    if (sum(raw == as.raw(10L)) > n) break
  }

  lines <- strsplit(rawToChar(raw), "\n", fixed = TRUE)[[1L]]
  lines <- sub("\r$", "", lines)
  utils::tail(lines[nzchar(lines)], n)
}

# Split slash-delimited diploid genotypes into adjacent allele columns.
split_genotype_columns <- function(x) {
  pieces <- lapply(x, function(column) {
    values <- stringi::stri_split_fixed(
      as.character(column),
      pattern = "/",
      simplify = TRUE
    )
    tibble::as_tibble(values, .name_repair = "minimal")
  })
  dplyr::bind_cols(pieces)
}

# Join sample and marker metadata back to a tidy genotype table.
join_genome_metadata <- function(x, s, m, g = NULL, env.arg = NULL) {
  if (!is.null(g)) {
    x <- dplyr::left_join(x, g, by = intersect(names(x), names(g)))
  }
  x <- dplyr::left_join(x, s, by = intersect(names(x), names(s)))
  x <- dplyr::select(x, -tidyselect::any_of(c("ID_SEQ", "STRATA_SEQ")))
  x <- dplyr::left_join(x, m, by = intersect(names(x), names(m)))
  x <- dplyr::select(x, -tidyselect::any_of("M_SEQ"))
  wanted <- c(
    "VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT",
    "INDIVIDUALS", "STRATA", "POP_ID", "GT_VCF", "GT_VCF_NUC", "GT",
    "ALT_DOSAGE"
  )
  dplyr::select(x, tidyselect::any_of(wanted), tidyselect::everything())
}

genome_results_message <- function(
    rad.message = NULL,
    filters.parameters,
    internal = FALSE,
    verbose = TRUE
) {
  if (internal) return(invisible(NULL))
  if (!is.null(rad.message)) message(rad.message)
  values <- filters.parameters$filters.parameters
  message("Number of individuals / strata / chrom / locus / SNP:")
  if (verbose) message("    Before: ", values$BEFORE)
  message("    Blacklisted: ", values$BLACKLIST)
  if (verbose) message("    After: ", values$AFTER)
  invisible(NULL)
}

dots_keepers_core <- function() {
  unique(c(
    "path.folder", "filename", "parameters", "internal", "random.seed",
    "blacklist.genotypes", "blacklist.id", "whitelist.markers",
    "gt", "alt.dosage", "gt.vcf", "gt.vcf.nuc", "calibrate.alleles",
    "keep.allele.names", "keep.gds", "markers.info", "vcf.metadata",
    "vcf.stats", "wide", "write.tidy", "tidy.check", "tidy.vcf",
    "tidy.dart", "pop.levels", "pop.labels", "pop.select"
  ))
}

dots_deprecated_core <- function() {
  c(
    "maf.thresholds", "common.markers", "max.marker", "monomorphic.out",
    "snp.ld", "filter.call.rate", "filter.markers.coverage",
    "filter.markers.missing", "number.snp.reads", "mixed.genomes.analysis",
    "duplicate.genomes.analysis", "ref.calibration"
  )
}

# Determine whether input coordinates come from a reference-guided assembly.
detect_ref_genome <- function(data = NULL, verbose = TRUE) {
  if (is.null(data)) return(FALSE)
  data.type <- detect_genomic_format(data)

  if (identical(data.type, "vcf.file")) {
    header <- check_header_source_vcf(vcf = data)$check.header
    result <- !is.null(header$reference) ||
      (!is.null(header$contig) && nrow(header$contig) > 0L)
  } else if (identical(data.type, "dart")) {
    result <- isTRUE(detect_dart_format(data, verbose = FALSE)$ref.genome)
  } else if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    opened <- identical(data.type, "gds.file")
    if (opened) data <- read_genome(data, verbose = FALSE)
    if (opened) on.exit(SeqArray::seqClose(data), add = TRUE)

    root <- genome_metadata_path(data)
    node <- gdsfmt::index.gdsn(
      data,
      paste0(root, "/reference.genome"),
      silent = TRUE
    )
    stored <- if (is.null(node)) NULL else tryCatch(
      gdsfmt::read.gdsn(node),
      error = function(e) NULL
    )
    if (is.logical(stored) && length(stored) == 1L) {
      result <- stored
    } else {
      markers <- extract_markers_metadata(
        data,
        markers.meta.select = c("CHROM", "POS"),
        whitelist = TRUE,
        verbose = FALSE
      )
      result <- "CHROM" %in% names(markers) &&
        any(!is.na(markers$CHROM) & nzchar(as.character(markers$CHROM)))
    }
  } else {
    rlang::abort("Reference-genome detection requires VCF, DArT, or GDS input.")
  }

  if (verbose) {
    message("Reads assembly: ", if (isTRUE(result)) "reference-assisted" else "de novo")
  }
  isTRUE(result)
}
