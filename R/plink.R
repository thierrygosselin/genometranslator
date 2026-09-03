# PLINK helpers ----------------------------------------------------------------

.plink2_path <- function(plink.path = NULL) {
  candidates <- if (is.null(plink.path)) c("plink2", "plink") else plink.path
  paths <- Sys.which(candidates)
  paths <- unname(paths[nzchar(paths)])
  if (!length(paths)) {
    rlang::abort(c(
      "PLINK 2 is required for this conversion.",
      "i" = "Install PLINK 2 and place `plink2` on PATH, or supply `plink.path`."
    ))
  }
  version <- suppressWarnings(system2(
    paths[[1]], "--version", stdout = TRUE, stderr = TRUE
  ))
  if (!any(grepl("PLINK v?2|PLINK 2", version, ignore.case = TRUE))) {
    rlang::abort(c(
      "The detected executable is not PLINK 2.",
      "i" = paste(version, collapse = " ")
    ))
  }
  paths[[1]]
}

.plink_prefix <- function(data, format) {
  extension <- switch(
    format,
    "plink.pgen.file" = "\\.pgen$",
    "plink.bed.file" = "\\.bed$",
    "plink.ped.file" = "\\.ped$",
    "plink.tped.file" = "\\.tped$",
    rlang::abort(paste0("Unsupported PLINK format: ", format))
  )
  sub(extension, "", data, ignore.case = TRUE)
}

.plink2_input_args <- function(data, format) {
  prefix <- .plink_prefix(data, format)
  switch(
    format,
    "plink.pgen.file" = c(
      "--pfile", prefix,
      if (!file.exists(paste0(prefix, ".pvar")) &&
          file.exists(paste0(prefix, ".pvar.zst"))) "vzs"
    ),
    "plink.bed.file" = c("--bfile", prefix),
    "plink.ped.file" = c("--file", prefix),
    "plink.tped.file" = c("--tfile", prefix)
  )
}

.run_plink2 <- function(plink.path, args, verbose = TRUE) {
  output <- system2(
    command = plink.path,
    args = vapply(args, shQuote, character(1)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    rlang::abort(c(
      "PLINK 2 conversion failed.",
      "i" = paste(utils::tail(output, 12L), collapse = "\n")
    ))
  }
  if (verbose && length(output)) {
    message(paste(utils::tail(output, 3L), collapse = "\n"))
  }
  invisible(output)
}

.plink_to_vcf <- function(data, format, plink.path = NULL, verbose = TRUE) {
  executable <- .plink2_path(plink.path)
  temp.prefix <- tempfile(pattern = "genometranslator_plink_")
  args <- c(
    .plink2_input_args(data, format),
    "--allow-extra-chr",
    "--export", "vcf",
    "--out", temp.prefix
  )
  .run_plink2(executable, args, verbose = verbose)
  vcf <- paste0(temp.prefix, ".vcf")
  if (!file.exists(vcf)) {
    rlang::abort("PLINK 2 completed without creating the expected VCF file.")
  }
  list(vcf = vcf, prefix = temp.prefix)
}

# Read PLINK -------------------------------------------------------------------
#' @name read_plink
#' @title Read PLINK
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Binary BED input uses the declared Bioconductor dependencies \pkg{SeqArray}
#' and \pkg{gdsfmt}. PLINK 2 is required to normalize PGEN/PVAR/PSAM,
#' PED/MAP, and TPED/TFAM input.
#' @description Read PLINK 1 BED/BIM/FAM, PED/MAP, or TPED/TFAM files and
#' PLINK 2 PGEN/PVAR/PSAM files. Inputs are converted to
#' a connection SeqArray \href{https://github.com/zhengxwen/SeqArray}{SeqArray}
#' GDS object/file of class \code{SeqVarGDSClass} (Zheng et al. 2017).
#' The Genomic Data Structure (GDS) file format is detailed in
#' \href{https://github.com/zhengxwen/gdsfmt}{gdsfmt}.
#'
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users.

#' @param data Path to the primary PLINK genotype file: \code{.pgen},
#' \code{.bed}, \code{.ped}, or \code{.tped}. Companion files with the same
#' prefix must be present.
#' \itemize{
#' \item PLINK 2: \code{pgen/pvar/psam} (\code{pvar.zst} is accepted).
#' \item PLINK 1 binary: \code{bed/bim/fam}.
#' \item PLINK 1 text: \code{ped/map} or \code{tped/tfam}.
#' }
#' @param strata Optional sample metadata passed to the VCF/GDS importer for
#' formats normalized with PLINK 2.
#' @param plink.path Optional path to the PLINK 2 executable. With
#' \code{NULL}, \code{plink2} and then \code{plink} are searched on PATH.

#' @param filename (optional) The file name of the Genomic Data Structure (GDS) file.
#' genometranslator will append \code{.gds} to the filename.
#' If the filename chosen exists in the working directory,
#' the default \code{radiator_datetime.gds} is chosen.
#' Default: \code{filename = NULL}.

#' @inheritParams tidy_genome
#' @inheritParams genometranslator_common_arguments


#' @details BED is imported directly with SeqArray. PGEN, PED, and TPED are
#' converted to a temporary VCF with PLINK 2 and then imported to GDS. The
#' temporary conversion files are removed after a successful import.


#### To do ....

# @section Advance mode:
#
# \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
# \enumerate{
# \item \code{path.folder}: to write ouput in a specific path
# (used internally in radiator). Default: \code{path.folder = getwd()}.
# If the supplied directory doesn't exist, it's created.
# \item \code{random.seed}: (integer, optional) For reproducibility, set an integer
# that will be used inside codes that uses randomness. With default,
# a random number is generated, printed and written in the directory.
# Default: \code{random.seed = NULL}.
# \item \code{subsample.markers.stats}: By default, when no filters are
# requested and that the number of markers is > 200K,
# 0.2 of markers are randomly selected to generate the
# statistics (individuals and markers). This is an all-around and
# reliable number.
# In doubt, overwrite this value by using 1 (all markers selected) and
# expect a small computational cost.
# }

#' @return
#' A writable \code{SeqVarGDSClass} connection for every supported flavour.

#' @export
#' @rdname read_plink

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795


#' @examples
#' \dontrun{
#' modern <- genometranslator::read_plink("study.pgen")
#' admixture <- genometranslator::read_plink("study.bed")
#' legacy <- genometranslator::read_plink("study.tped")
#' }

#' @seealso
#' \href{https://www.cog-genomics.org/plink/2.0/formats}{PLINK 2 formats}


read_plink <- function(
  data,
  filename = NULL,
  strata = NULL,
  plink.path = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
) {

  # #Test
  # filename <- NULL
  # parallel.core <- parallel::detectCores() - 1
  # verbose <- TRUE

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "read_plink",
    verbose = verbose
  )
  file.date <- .start$file.date
  on.exit(tgbase::teardown(.start), add = TRUE)
  res <- list()

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("PLINK file missing")

  # Warning concerning plink files----------------------------------------------
  plink.format <- genometranslator::detect_genomic_format(data = data)

  # PLINK 2 and PLINK 1 text formats are normalized with the official PLINK 2
  # engine. This avoids maintaining independent PED/TPED/PGEN parsers and makes
  # all readers return the same persistent GDS representation.
  if (plink.format %in% c(
    "plink.pgen.file", "plink.ped.file", "plink.tped.file"
  )) {
    if (verbose) message("Normalizing ", plink.format, " with PLINK 2...")
    converted <- .plink_to_vcf(
      data = data,
      format = plink.format,
      plink.path = plink.path,
      verbose = verbose
    )
    on.exit(unlink(
      Sys.glob(paste0(converted$prefix, ".*")),
      recursive = FALSE,
      force = TRUE
    ), add = TRUE)
    return(genometranslator::read_vcf(
      data = converted$vcf,
      strata = strata,
      filename = filename,
      parallel.core = parallel.core,
      verbose = verbose
    ))
  }

  # PLINK TPED -----------------------------------------------------------------
  if (FALSE && plink.format == "plink.tped.file") {
    message("Reading PLINK tped...")

    # Get file size
    plink.size <- file.size(data)
    if (!plink.size <= 100000000) {
      message("\n\nNote: huge PLINK tped files are no longer recommended inside radiator")
      message("Convert to PLINK BED with: --tfile [] --make-bed --allow-no-sex --allow-extra-chr --chr-set [] with PLINK command")
      rlang::abort("Read the documentation of genometranslator::read_plink")
    }

    # Importing file -----------------------------------------------------------
    fam.file <- stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".tped",
      replacement = ".tfam",
      vectorize_all = FALSE
    )
    res$strata <- readr::read_delim(
      file = fam.file,
      delim = " ",
      col_names = c("STRATA", "INDIVIDUALS"),
      col_types = "cc____"
    ) %>%
      genometranslator::read_strata(strata = .) %$%
      strata
    fam.file <- NULL

    # preparing header for tped file
    tped.header.prep <- res$strata %>%
      dplyr::select(INDIVIDUALS) %>%
      dplyr::mutate(
        NUMBER = seq(1, n()),
        ALLELE1 = rep("A1", n()), ALLELE2 = rep("A2", n())
      ) %>%
      tgbase::trans_long(
        x = .,
        cols = c("INDIVIDUALS", "NUMBER"),
        names_to = "ALLELES_GROUP",
        values_to = "ALLELES"
      ) %>%
      dplyr::arrange(NUMBER) %>%
      dplyr::select(-ALLELES_GROUP) %>%
      tidyr::unite(INDIVIDUALS_ALLELES, c(INDIVIDUALS, ALLELES), sep = "_", remove = FALSE) %>%
      dplyr::arrange(NUMBER) %>%
      dplyr::mutate(NUMBER = seq(from = (1 + 4), to = n() + 4)) %>%
      dplyr::select(-ALLELES)
    tped.header.names <- c("CHROM", "LOCUS", "POS", tped.header.prep$INDIVIDUALS_ALLELES)
    tped.header.integer <- c(1, 2, 4, tped.header.prep$NUMBER)
    tped.header.prep <- NULL

    # import PLINK
    res$data <- data.table::fread(
      input = data,
      sep = " ",
      header = FALSE,
      stringsAsFactors = FALSE,
      verbose = FALSE,
      select = tped.header.integer,
      col.names = tped.header.names,
      showProgress = TRUE,
      data.table = FALSE) %>%
      tibble::as_tibble(.) %>%
      dplyr::mutate(
        CHROM = as.character(CHROM),
        LOCUS = as.character(LOCUS),
        POS = as.character(POS),
        MARKERS = make_marker_id(CHROM, LOCUS, POS)
      )

    # Unused objects
    tped.header.integer <- tped.header.names <- NULL
    return(res)
  }

  # PLINK BED ------------------------------------------------------------------
  if (plink.format == "plink.bed.file") {
    message("Reading PLINK bed file...")

    # Required package -----------------------------------------------------------
    tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)

    # Function call and dotslist -------------------------------------------------
    rad.dots <- genometranslator_dots(
      func.name = as.list(sys.call())[[1]],
      fd = rlang::fn_fmls_names(),
      args.list = as.list(environment()),
      dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
      keepers = c("internal", "path.folder", "parameters"),
      verbose = FALSE
    )

    # Generate folders and filenames ---------------------------------------------
    wf <- path.folder <- tgbase::generate_folder(
      folder = "read_plink",
      path.folder = path.folder,
      prefix.int = FALSE,
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    radiator.folder <- tgbase::generate_folder(
      folder = "import_gds",
      path.folder = path.folder,
      prefix.int = TRUE,
      internal = internal,
      file.date = file.date,
      verbose = verbose)

    # write the dots file
    tgbase::write_tgbase_tsv(
      data = rad.dots,
      path.folder = path.folder,
      filename = "radiator_read_plink_args",
      date = TRUE,
      internal = internal,
      write.message = "Function call and arguments stored in: ",
      verbose = verbose
    )

    filename <- generate_filename(
      name.shortcut = filename,
      path.folder = radiator.folder,
      extension = "gds")

    filename.short <- filename$filename.short
    filename <- filename$filename

    timing.plink <- proc.time()

    # read_plink -----------------------------------------------------------------
    bim.file <- stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".bed",
      replacement = ".bim",
      vectorize_all = FALSE
    )
    fam.file <- stringi::stri_replace_all_fixed(
      str = data,
      pattern = ".bed",
      replacement = ".fam",
      vectorize_all = FALSE
    )

    if (verbose) message("Generating GDS...")
    gds <- SeqArray::seqBED2GDS(
      bed.fn = data,
      fam.fn = fam.file,
      bim.fn = bim.file,
      out.gdsfn = filename,
      compress.geno = "ZIP_RA",
      compress.annotation = "ZIP_RA",
      verbose = verbose
    ) %>%
      SeqArray::seqOpen(gds.fn = ., readonly = FALSE)
    if (verbose) message("done! timing: ", round((proc.time() - timing.plink)[[3]]), " sec\n")

    # PLINK: Summary ----------------------------------------------------------------
    summary_gds(gds = gds, verbose = TRUE)
    if (verbose) message("\nFile written: ", filename.short)

    # Generate radiator skeleton -------------------------------------------------
    metadata.node <- genome_gds_skeleton(gds)

    # data.source ----------------------------------------------------------------
    update_genome_gds(gds = gds, node.name = "data.source", value = "plink")

    # bi- or multi-alllelic--------------------------------------------------
    biallelic <- detect_biallelic_markers(data = gds, verbose = verbose)

    # Clean sample id---------------------------------------------------------
    individuals.plink <- tibble::tibble(
      INDIVIDUALS_PLINK = SeqArray::seqGetData(gds, "sample.id")) %>%
      dplyr::mutate(INDIVIDUALS_CLEAN = clean_ind_names(INDIVIDUALS_PLINK))

    if (!identical(individuals.plink$INDIVIDUALS_PLINK, individuals.plink$INDIVIDUALS_CLEAN)) {
      if (verbose) message("Cleaning PLINK's sample names")
      clean.id.filename <- stringi::stri_join("cleaned.plink.id.info_", file.date, ".tsv")
      readr::write_tsv(x = individuals.plink,
                       file = file.path(radiator.folder, clean.id.filename))

      update_genome_gds(gds = gds, node.name = "id.clean", value = individuals.plink)
    }
    # replace id
    update_genome_gds(
      gds = gds,
      metadata.node = FALSE,
      node.name = "sample.id",
      value = individuals.plink$INDIVIDUALS_CLEAN,
      replace = TRUE)

    individuals <- dplyr::select(individuals.plink, INDIVIDUALS = INDIVIDUALS_CLEAN)
    individuals.plink <- NULL

    # sync id with STRATA---------------------------------------------------------
    if (verbose) message("Using .fam file for strata...")
    strata <- readr::read_delim(
      file = fam.file,
      delim = " ",
      col_names = c("STRATA", "INDIVIDUALS"),
      col_types = "cc____"
    ) %>%
      genometranslator::read_strata(strata = .) %$%
      strata

    id.levels <- individuals$INDIVIDUALS
    individuals %<>%
      dplyr::left_join(
        join_strata(individuals, strata, verbose = verbose) %>%
          dplyr::mutate(FILTERS = "whitelist")
        , by = "INDIVIDUALS"
      ) %>%
      dplyr::mutate(FILTERS = tidyr::replace_na(data = FILTERS, replace = "filter.stata"))

    strata <- generate_strata(data = dplyr::filter(individuals, FILTERS == "whitelist"), pop.id = FALSE)

    #Update GDS node
    update_genome_gds(gds = gds, node.name = "individuals.meta", value = individuals, sync = TRUE)

    # PLINK: Markers metadata  ------------------------------------------------------
    markers.meta <- extract_markers_metadata(gds = gds)

    # PLINK: reference genome or de novo -------------------------------------------
    ref.genome <- detect_ref_genome(data = gds, verbose = verbose)

    # Generate MARKERS column and fix types --------------------------------------
    markers.meta %<>%
      dplyr::mutate(
        dplyr::across(
          .cols = c(LOCUS, POS),
          .fns = clean_markers_names
        )
      ) %>%
      dplyr::mutate(
        MARKERS = make_marker_id(CHROM, LOCUS, POS),
        REF = SeqArray::seqGetData(gdsfile = gds, var.name = "$ref"),
        ALT = SeqArray::seqGetData(gdsfile = gds, var.name = "$alt")
      )

    # PLINK file with duplicate markers... sometimes tagged ISOFORMS...
    dup.markers <- length(markers.meta$MARKERS) - length(unique(markers.meta$MARKERS))
    if (dup.markers > 0) {
      message("\nNumber of duplicate MARKERS id: ", dup.markers)
      message("Adding integer to differentiate...")
      markers.meta %<>%
        dplyr::arrange(MARKERS) %>%
        dplyr::mutate(MARKERS_NEW = MARKERS) %>%
        dplyr::group_by(MARKERS_NEW) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_join(MARKERS, seq(1, n(), by = 1), sep = "_")
        ) %>%
        dplyr::ungroup(.) %>%
        dplyr::select(-MARKERS_NEW)
    }

    # ADD MARKERS META to GDS
    update_genome_gds(gds = gds, node.name = "markers.meta", value = markers.meta)

    # replace chromosome info in GDS
    # Why ? well snp ld e.g. will otherwise be performed by chromosome and with de novo data = by locus...
    update_genome_gds(
      gds = gds,
      metadata.node = FALSE,
      node.name = "chromosome",
      value = markers.meta$CHROM,
      replace = TRUE
    )
    # # radiator_parameters: generate --------------------------------------------
    filters.parameters <- genome_parameters(
      generate = TRUE,
      initiate = FALSE,
      update = FALSE,
      parameter.obj = parameters,
      path.folder = radiator.folder,
      file.date = file.date,
      verbose = verbose,
      internal = internal)

    # radiator_parameters: initiate --------------------------------------------
    # with original PLINK values
    filters.parameters <- genome_parameters(
      generate = FALSE,
      initiate = TRUE,
      update = TRUE,
      parameter.obj = filters.parameters,
      data = gds,
      filter.name = "plink",
      param.name = "original values in PLINK + strata",
      values = "",
      path.folder = path.folder,
      file.date = file.date,
      internal = internal,
      verbose = verbose
    )
    return(gds)
  } # end plink's BED format
} # End read_plink


# tidy_plink -------------------------------------------------------------------
#' @name tidy_plink
#' @title Tidy PLINK tped and bed files

#' @description Transform bi-allelic PLINK files in \code{.tped} or \code{.bed} formats
#' into a tidy dataset.
#'
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.
#' Ensure marker type, missingness, allele coding, and any linkage pruning are
#' appropriate for the intended PLINK analysis before export.

#' @param data The PLINK file.
#' \itemize{
#' \item {bi-allelic data only}. For haplotypes use VCF.
#' \item \code{tped} file format: the corresponding \code{tfam} file must be in the directory.
#' \item \code{bed} file format: IS THE PREFERRED format, the corresponding
#' \code{fam} and \code{bim} files must be in the directory.
#' }

#' @inheritParams genometranslator_common_arguments

#' @section Advance mode:
#'
#' \emph{dots-dots-dots ...} allows to pass several arguments for fine-tuning the function:
#' \enumerate{
#' \item \code{calibrate.alleles}: logical. For \code{tped} files, if
#' \code{calibrate.alleles = FALSE} the function runs faster
#' but REF/ALT alleles may not be calibrated. The default assumes the users or
#' sotware producing the PLINK file calibrated the alleles.
#' Default: \code{calibrate.alleles = FALSE}.
#' }


#' @references Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS,
#' Laurie C, Levine D (2017). SeqArray -- A storage-efficient high-performance
#' data format for WGS variant calls.
#' Bioinformatics.
#'
#' @references
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795

#' @examples
#' \dontrun{
#' data <- genometranslator::tidy_plink(data = "my_plink_file.bed", verbose = TRUE)
#'
#'
#' # when conversion is required from TPED to BED, in Terminal:
#' # plink --tfile my_plink_file --make-bed --allow-no-sex --allow-extra-chr --chr-set 95
#' }

#' @seealso
#' \href{https://www.cog-genomics.org/plink/1.9/}{PLINK}
#'
#' \code{\link{read_plink}}

#' @return
#' A tidy tibble of the PLINK file.
#'


#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname tidy_plink
#' @keywords internal
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}

tidy_plink <- function(
  data,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
) {
  # Test
  # verbose = TRUE
  # strata = NULL
  # calibrate.alleles = TRUE
  # parallel.core = parallel::detectCores() - 1
  # parameters <- NULL
  # filename=NULL
  # internal = FALSE
  # path.folder = getwd()

  # Cleanup---------------------------------------------------------------------
  tgbase::function_header(f.name = "tidy_plink", verbose = verbose)
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")
  if (verbose) message("Execution date@time: ", file.date)
  old.dir <- getwd()
  opt.change <- getOption("width")
  options(future.globals.maxSize = Inf)
  options(width = 70)
  timing <- tgbase::tic()
  #back to the original directory and options
  on.exit(setwd(old.dir), add = TRUE)
  on.exit(options(width = opt.change), add = TRUE)
  on.exit(tgbase::toc(timing), add = TRUE)
  on.exit(tgbase::function_header(f.name = "tidy_plink", start = FALSE, verbose = verbose), add = TRUE)

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("PLINK file missing")

  # Function call and dotslist -------------------------------------------------
  path.folder <- filename <- calibrate.alleles <- NULL
  rad.dots <- genometranslator_dots(
    func.name = as.list(sys.call())[[1]],
    fd = rlang::fn_fmls_names(),
    args.list = as.list(environment()),
    dotslist = rlang::dots_list(..., .homonyms = "error", .check_assign = TRUE),
    deprecated = c("blacklist.id", "pop.select", "pop.levels", "pop.labels"),
    keepers = c("calibrate.alleles", "filename", "internal", "path.folder", "parameters"),
    verbose = FALSE
  )

  # Warning concerning plink files----------------------------------------------
  plink.format <- genometranslator::detect_genomic_format(data = data)

  data <- genometranslator::read_plink(
    data = data,
    filename = filename,
    parallel.core = parallel.core,
    verbose = FALSE,
    internal = TRUE,
    path.folder = path.folder,
    parameters = parameters
  )

  if (FALSE && plink.format == "plink.tped.file") {
    # Make tidy ------------------------------------------------------------------
    if (verbose) message("Tidying the PLINK tped file ...")
    # Filling GT and new separating INDIVIDUALS from ALLELES
    # combining alleles
    strata <- data$strata
    data <- tgbase::trans_long(
      x = data$data,
      cols = c("MARKERS", "CHROM", "LOCUS", "POS"),
      names_to = "INDIVIDUALS_ALLELES",
      values_to = "GT"
    )

    # detect GT coding
    detect.gt.coding <- unique(sample(x = data$GT, size = 100, replace = FALSE))
    gt.letters <- c("A", "C", "G", "T")

    if (TRUE %in% unique(gt.letters %in% detect.gt.coding)) {
      if (verbose) message("Genotypes coded with letters")
      gt.letters.df <- tibble::tibble(
        GT = c("A", "C", "G", "T", "0"),
        NEW_GT = c("001", "002", "003", "004", "000")
      )
      data  %<>%
        dplyr::left_join(
          gt.letters.df, by = "GT") %>%
        dplyr::select(-GT) %>%
        dplyr::rename(GT = NEW_GT)
      gt.letters.df <- NULL
    } else {
      if (verbose) message("Genotypes coded with integers")
      data  %<>%
        dplyr::mutate(GT = stringi::stri_pad_left(str = GT, width = 3, pad = "0"))
    }
    detect.gt.coding <- gt.letters <- NULL

    data %<>%
      tidyr::separate(
        data = .,
        col = INDIVIDUALS_ALLELES,
        into = c("INDIVIDUALS", "ALLELES"),
        sep = "_") %>%
      tgbase::trans_wide(
        x = .,
        formula = "MARKERS + CHROM + LOCUS + POS + INDIVIDUALS ~ ALLELES",
        values_from = "GT"
      ) %>%
      dplyr::ungroup(.) %>%
      tidyr::unite(data = ., col = GT, A1, A2, sep = "") %>%
      dplyr::select(MARKERS, CHROM, LOCUS, POS, INDIVIDUALS, GT)

    # population levels and strata
    if (verbose) message("Integrating the tfam/strata file...")

    data %<>% dplyr::left_join(strata, by = "INDIVIDUALS")
    strata <- NULL

    # removing untyped markers across all-pop
    remove.missing.gt <- data %>%
      dplyr::select(LOCUS, GT) %>%
      dplyr::filter(GT != "000000")

    untyped.markers <- dplyr::n_distinct(data$LOCUS) - dplyr::n_distinct(remove.missing.gt$LOCUS)
    if (untyped.markers > 0) {
      if (verbose) message("Number of marker with 100 % missing genotypes: ", untyped.markers)
      data <- suppressWarnings(
        dplyr::semi_join(data,
                         remove.missing.gt %>%
                           dplyr::distinct(LOCUS, .keep_all = TRUE),
                         by = "LOCUS")
      )
    }

    # Unused objects
    remove.missing.gt <- NULL

    # detect if biallelic give vcf style genotypes
    # biallelic <- detect_biallelic_markers(input)
    # filename <- internal <- parameters <- path.folder <- calibrate.alleles <- NULL
    # rm(filename, internal, parameters, path.folder, calibrate.alleles)

    if (calibrate.alleles) {
      data %<>% calibrate_alleles(data = ., verbose = verbose)
      return(res = list(input = data$input, biallelic = data$biallelic))
    } else {
      return(res = list(input = data, biallelic = detect_biallelic_markers(data)))
    }
  } #End tidy tped

  if (plink.format %in% c(
    "plink.pgen.file", "plink.bed.file",
    "plink.ped.file", "plink.tped.file"
  )) {
    # Get info markers and individuals -----------------------------------------
    gds.info <- genometranslator::summary_gds(gds = data, verbose = FALSE)
    n.markers <- gds.info$n.markers
    n.individuals <- gds.info$n.ind

    cat("\n\n################################## IMPORTANT ###################################\n")
    message("Tidying PLINK file with ", n.markers, " SNPs is not optimal:")
    message("    1. a computer with lots of RAM is required")
    message("    2. it's very slow to generate")
    message("    3. it's very slow to run codes after")
    message("    4. for most non model species this number of markers is not realistic...")
    message("\nRecommendation:")
    message("    1. filter your dataset. e.g. with explore_genomes")
    message("\nIdeally target a maximum of ~ 10 000 - 20 0000 unlinked SNPs\n")

    message("\nGenotypes formats generated with ", n.markers, " SNPs: ")
    message("    ALT_DOSAGE (the dosage of ALT allele: 0, 1, 2 NA)")

    biallelic <- detect_biallelic_markers(data)
    tidy.data <- tidy_genome(
      data = data,
      markers.meta = NULL,
      pop.id = FALSE,
      calibrate.alleles = calibrate.alleles,
      close.gds = TRUE
    )
    return(res = list(input = tidy.data, biallelic = biallelic))
  }#End tidy bed
} # End tidy_plink


# write_plink ------------------------------------------------------------------

#' @name write_plink
#' @title Write PLINK

#' @description Export genomic data to PLINK 2 PGEN/PVAR/PSAM or PLINK 1
#' BED/BIM/FAM, PED/MAP, or TPED/TFAM. The active data are first represented as
#' VCF and the official PLINK 2 engine creates the requested fileset.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @param data A tidy data frame object in the global environment or
#' a tidy data frame in wide or long format in the working directory.
#' \emph{How to get a tidy data frame ?}
#' Look into \pkg{genometranslator} \code{\link{tidy_genome}}.

#' @param filename Output prefix without a PLINK extension. With
#' \code{NULL}, a timestamped prefix is created in \code{path.folder}.
#' @param format Output flavour: \code{"pgen"} (PGEN/PVAR/PSAM),
#' \code{"bed"} (BED/BIM/FAM), \code{"ped"} (PED/MAP), or \code{"tped"}
#' (TPED/TFAM). Default: \code{"pgen"}.
#' @param path.folder Output directory. Default: current working directory.
#' @param plink.path Optional path to PLINK 2. With \code{NULL}, executable
#' names \code{plink2} and then \code{plink} are searched on PATH.
#' @param plink.args Additional PLINK 2 command-line arguments. These are
#' appended before \code{--out}; use only documented PLINK 2 flags.
#' @param overwrite Logical. Permit replacement of an existing fileset.
#' Default: \code{FALSE}.
#' @param verbose Logical. Print PLINK 2 progress. Default: \code{TRUE}.
#' @details PGEN is the recommended general PLINK output. BED is useful for
#' ADMIXTURE and other software that expects PLINK 1 binary input. PED and TPED
#' are uncompressed legacy exchange formats. BED, PED, and TPED workflows should
#' normally use filtered biallelic markers; this function does not silently
#' remove unsuitable variants.
#' @return Invisibly, the paths of every file in the generated PLINK fileset.
#' @examples
#' \dontrun{
#' # Modern PLINK 2 files
#' genometranslator::write_plink(genome, "study", format = "pgen")
#'
#' # PLINK 1 binary files for ADMIXTURE and legacy pipelines
#' genometranslator::write_plink(genome, "study_admixture", format = "bed")
#'
#' # Text-based legacy exports
#' genometranslator::write_plink(genome, "study", format = "ped")
#' genometranslator::write_plink(genome, "study", format = "tped")
#' }
#' @export
#' @rdname write_plink
#' @references Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR,
#' Bender D, et al.
#' PLINK: a tool set for whole-genome association and population-based linkage
#' analyses.
#' American Journal of Human Genetics. 2007: 81: 559–575. doi:10.1086/519795


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering
#' @template io-dependencies

write_plink <- function(
    data,
    filename = NULL,
    format = c("pgen", "bed", "ped", "tped"),
    path.folder = getwd(),
    plink.path = NULL,
    plink.args = character(),
    overwrite = FALSE,
    verbose = TRUE
) {
  format <- match.arg(format)
  if (!dir.exists(path.folder)) {
    dir.create(path.folder, recursive = TRUE, showWarnings = FALSE)
  }
  path.folder <- normalizePath(path.folder, mustWork = TRUE)

  if (is.null(filename)) {
    filename <- paste0(
      "genometranslator_plink_",
      base::format(Sys.time(), "%Y%m%d@%H%M%S")
    )
  }
  filename <- basename(sub(
    "\\.(pgen|pvar|psam|bed|bim|fam|ped|map|tped|tfam)$",
    "", filename, ignore.case = TRUE
  ))
  output.prefix <- file.path(path.folder, filename)
  extensions <- switch(
    format,
    pgen = c(".pgen", ".pvar", ".psam"),
    bed = c(".bed", ".bim", ".fam"),
    ped = c(".ped", ".map"),
    tped = c(".tped", ".tfam")
  )
  output.files <- paste0(output.prefix, extensions)
  existing <- output.files[file.exists(output.files)]
  if (length(existing) && !overwrite) {
    rlang::abort(c(
      "The requested PLINK fileset already exists.",
      "i" = paste(basename(existing), collapse = ", "),
      "i" = "Use a different `filename` or set `overwrite = TRUE`."
    ))
  }

  executable <- .plink2_path(plink.path)
  temp.prefix <- tempfile(pattern = "genometranslator_write_plink_")
  vcf <- paste0(temp.prefix, ".vcf")
  opened.gds <- NULL

  if (inherits(data, "SeqVarGDSClass")) {
    SeqArray::seqGDS2VCF(data, vcf.fn = vcf, verbose = verbose)
  } else if (is.character(data) && length(data) == 1L &&
             grepl("\\.gds$", data, ignore.case = TRUE)) {
    opened.gds <- genometranslator::read_genome(data, verbose = FALSE)
    on.exit(SeqArray::seqClose(opened.gds), add = TRUE)
    SeqArray::seqGDS2VCF(opened.gds, vcf.fn = vcf, verbose = verbose)
  } else if (is.character(data) && length(data) == 1L &&
             grepl("\\.vcf(?:\\.gz)?$", data, ignore.case = TRUE)) {
    vcf <- normalizePath(data, mustWork = TRUE)
  } else {
    genometranslator::write_vcf(data = data, filename = temp.prefix)
  }

  if (!file.exists(vcf)) {
    rlang::abort("Unable to create the intermediate VCF for PLINK export.")
  }
  if (startsWith(vcf, tempdir())) {
    on.exit(unlink(vcf, force = TRUE), add = TRUE)
  }

  action <- switch(
    format,
    pgen = c("--make-pgen"),
    bed = c("--make-bed"),
    ped = c("--export", "ped"),
    tped = c("--export", "tped")
  )
  args <- c(
    "--vcf", vcf,
    "--allow-extra-chr",
    action,
    as.character(plink.args),
    "--out", output.prefix
  )
  args <- args[nzchar(args)]
  .run_plink2(executable, args, verbose = verbose)

  generated <- output.files[file.exists(output.files)]
  if (length(generated) != length(output.files)) {
    rlang::abort(c(
      "PLINK 2 did not create the complete requested fileset.",
      "i" = paste("Expected:", paste(basename(output.files), collapse = ", "))
    ))
  }
  if (verbose) {
    message("PLINK ", toupper(format), " files written:")
    message(paste0("  ", generated, collapse = "\n"))
  }
  invisible(generated)
} # end write_plink
