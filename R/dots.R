# Package-specific optional-argument policy ----------------------------------
genometranslator_dots_defaults <- function(keepers) {
  defaults <- stats::setNames(vector("list", length(keepers)), keepers)
  true.names <- intersect(keepers, c(
    "keep.gds", "vcf.stats", "vcf.metadata", "filter.common.markers",
    "filter.monomorphic", "ld.figures", "dart.sequence", "force.stats"
  ))
  false.names <- intersect(keepers, c(
    "keep.allele.names", "calibrate.alleles", "long.ld.missing",
    "detect.mixed.genomes", "detect.duplicate.genomes", "dp", "internal",
    "heatmap.fst", "wide", "filter.hwe", "gt", "alt.dosage", "gt.vcf",
    "gt.vcf.nuc", "filter.haplotype.format"
  ))
  defaults[true.names] <- rep(list(TRUE), length(true.names))
  defaults[false.names] <- rep(list(FALSE), length(false.names))
  if ("filter.strands" %in% keepers) defaults[["filter.strands"]] <- "blacklist"
  if ("ld.method" %in% keepers) defaults[["ld.method"]] <- "r2"
  if ("iteration.subsample" %in% keepers) defaults[["iteration.subsample"]] <- 1L
  defaults
}

genometranslator_dots <- function(
  func.name = as.list(sys.call())[[1]],
  fd = NULL,
  args.list = NULL,
  dotslist = list(),
  keepers = character(),
  deprecated = character(),
  env = parent.frame(),
  assign = TRUE,
  verbose = TRUE
) {
  if (is.null(fd)) fd <- character()
  if (is.null(args.list)) args.list <- list()

  formal.names <- intersect(names(args.list), fd)
  formal.summary <- tibble::tibble(
    ARGUMENTS = formal.names,
    VALUES = vapply(args.list[formal.names], tgbase::format_argument, character(1)),
    GROUPS = "formal"
  )

  dots.summary <- tgbase::resolve_dots(
    dots = dotslist,
    allowed = keepers,
    deprecated = deprecated,
    defaults = genometranslator_dots_defaults(keepers),
    env = env,
    assign = assign,
    verbose = verbose
  )

  if (isTRUE(verbose)) message("Resolved optional arguments for ", func.name)
  dplyr::bind_rows(formal.summary, dots.summary)
}
