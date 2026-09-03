# read_genlight ----------------------------------------------------------------
#' @name read_genlight
#' @title Read a genlight object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Reading and writing \code{genlight} objects requires the optional CRAN
#' package \href{https://cran.r-project.org/package=adegenet}{\pkg{adegenet}}.
#' @description Import an \pkg{adegenet} \code{genlight} object, standardize its
#' sample and marker metadata, and optionally produce tidy genomic data or a GDS
#' representation.

#' @param data (path or object) A genlight object in the global environment or
#' path to a genlight file that will be open with \code{readRDS}.

#' @inheritParams genometranslator_common_arguments

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds Logical. Generate a genometranslator GDS representation.
#' Default: \code{gds = TRUE}.

#' @param write Logical. Write the tidy table as an Arrow/Parquet file.
#' Default: \code{write = FALSE}.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_genlight
#' @return With \code{tidy = TRUE}, a tidy tibble. With
#' \code{tidy = FALSE, gds = TRUE}, the generated GDS result. With both options
#' false, the original \code{genlight} object.
#' @examples
#' \dontrun{
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   x <- readRDS("genotypes.genlight.rds")
#'   genotypes <- genometranslator::read_genlight(x, gds = FALSE)
#' }
#' }

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @details
#' Missing marker fields are generated when necessary:
#' \enumerate{
#' \item \code{is.null(genlight@pop)}: pop will be integrated
#' in the tidy dataset.
#' \item \code{is.null(data@chromosome)}: \code{DENOVO} is used
#' in the tidy dataset.
#' \item \code{is.null(data@loc.names)}: LOCUS1 to \code{ncol(genlight)}
#' will be integrated in the tidy dataset.
#' \item \code{is.null(data@position)}: an integer string of
#' length = \code{ncol(genlight)} will be integrated in the tidy dataset.
#' }
#' Generated fields preserve a unique marker identifier but do not create
#' reference-genome coordinates. When \code{loc.all} is unavailable, symbolic
#' \code{REF = "0"} and \code{ALT = "1"} labels preserve the dosage orientation.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies


read_genlight <- function(
    data,
    tidy = TRUE,
    gds = TRUE,
    write = FALSE,
    verbose = FALSE,
    parallel.core = parallel::detectCores() - 1
) {
  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "read_genlight",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  # Test
  # data = "radiator_genlight_20191211@1836.RData"
  # tidy = TRUE
  # gds = TRUE
  # write = FALSE
  # verbose = TRUE
  # parallel.core = 12L


  # Package requirement --------------------------------------------------------
  tgbase::check_package(package = "adegenet")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("A genlight object or RDS file is required.")

  # Import data ---------------------------------------------------------------
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) rlang::abort("The genlight RDS file does not exist.")
    data <- readRDS(data)
  }
  if (!inherits(data, "genlight")) rlang::abort("Input is not a genlight object.")
  original.data <- data
  requested.tidy <- isTRUE(tidy)

  if (verbose) message("genlight info:")
  # strata ?
  strata <- tibble::tibble(INDIVIDUALS = data@ind.names)
  if (is.null(data@pop)) {
    if (verbose) message("    strata: no")
    if (verbose) message("    'pop' will be added")
    strata %<>% dplyr::mutate(STRATA = "pop")
  } else {
    if (verbose) message("    strata: yes")
    strata$STRATA = data@pop
  }

  n.markers <- dim(data)[2]
  n.ind <- nrow(strata)

  # Chromosome ?
  if (is.null(data@chromosome)) {
    if (verbose) message("    Chromosome/contig/scaffold: no")
    data@chromosome <- factor(rep("DENOVO", n.markers))
    chrom.info <- FALSE
  } else {
    if (verbose) message("    Chromosome/contig/scaffold: yes")
    chrom.info <- TRUE
  }

  # Locus ?
  if (is.null(data@loc.names)) {
    if (verbose) message("    Locus: no")
    locus.info <- FALSE
    data@loc.names <- stringi::stri_join("LOCUS", seq(from = 1, to = n.markers, by = 1))
  } else {
    if (verbose) message("    Locus: yes")
    locus.info <- TRUE
  }

  # POS ?
  if (is.null(data@position)) {
    if (verbose) message("    POS: no")
    pos.info <- FALSE
    # data@position <- rlang::as_integer(seq(from = 1, to = n.markers, by = 1))
    data@position <- vctrs::vec_cast(x = seq(from = 1, to = n.markers, by = 1), to = integer())
  } else {
    if (verbose) message("    POS: yes")
    pos.info <- TRUE
  }


  # markers
  markers <- tibble::tibble(
    CHROM = data@chromosome,#adegenet::chromosome(data),
    LOCUS = data@loc.names,#adegenet::locNames(data),
    POS = data@position#adegenet::position(data)
  ) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), .fns = as.character)) %>%
    dplyr::mutate(dplyr::across(c(LOCUS, POS), .fns = clean_markers_names)) %>%
    dplyr::mutate(MARKERS = make_marker_id(CHROM, LOCUS, POS)) %>%
    dplyr::select(MARKERS, CHROM, LOCUS, POS)

  # Nuc info
  nuc.data <- data@loc.all
  if (!is.null(nuc.data)) {
    markers %<>%
      dplyr::mutate(
        REF = stringi::stri_sub(str = nuc.data, from = 1, to = 1),
        ALT = stringi::stri_sub(str = nuc.data, from = 3, to = 3)
      )
  } else {
    # genlight dosage is defined relative to its second allele even when the
    # nucleotide labels were not retained. Symbolic labels preserve that
    # orientation for tidy and GDS conversion.
    markers %<>% dplyr::mutate(REF = "0", ALT = "1")
  }

  if (gds) tidy <- TRUE

  if (tidy) {
    if (write) {
      filename.temp <- generate_filename(extension = "arrow.parquet")
      filename.short <- filename.temp$filename.short
      filename.genlight <- filename.temp$filename
    }

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT","STRATA", "INDIVIDUALS",
              "GT_VCF", "ALT_DOSAGE", "GT")

    if (verbose) message("Generating tidy data...")
    tidy.data <- data.frame(data) %>%
      magrittr::set_colnames(x = ., value = markers$MARKERS) %>%
      tibble::add_column(.data = ., INDIVIDUALS = rownames(.), .before = 1) %>%
      tgbase::trans_long(
        x = .,
        cols = "INDIVIDUALS",
        names_to = "MARKERS",
        values_to = "ALT_DOSAGE"
      ) %>%
      dplyr::full_join(markers, by = "MARKERS") %>%
      dplyr::full_join(strata, by =  "INDIVIDUALS") %>%
      dplyr::mutate(
        INDIVIDUALS = clean_ind_names(INDIVIDUALS),
        STRATA = clean_pop_names(STRATA)
      ) %>%
      dplyr::arrange(MARKERS, STRATA, INDIVIDUALS) %>%
      dplyr::select(tidyselect::any_of(want)) %>%
      calibrate_alleles(
        data = .,
        biallelic = TRUE,
        verbose = verbose
      ) %$%
      input

    if (write) genometranslator::write_genome(data = tidy.data, filename = filename.genlight, verbose = verbose)

  }#End tidy genlight

  if (gds) {
    # generate the GDS --------------------------------------------------------------
    # markers %<>% dplyr::mutate(VARIANT_ID = as.integer(factor(MARKERS)))

    # generate genotypes format for easy reading into GDS
    tidy.data %<>%
      dplyr::mutate(
        GDS_A1 = dplyr::case_when(
          ALT_DOSAGE == 0 ~ 0,
          ALT_DOSAGE == 1 ~ 0,
          ALT_DOSAGE == 2 ~ 1,
          is.na(ALT_DOSAGE) ~ NA_integer_),
        GDS_A2 = dplyr::case_when(
          ALT_DOSAGE == 0 ~ 0,
          ALT_DOSAGE == 1 ~ 1,
          ALT_DOSAGE == 2 ~ 1,
          is.na(ALT_DOSAGE) ~ NA_integer_)
      )

    gds.filename <- genome_gds(
      data.source = "genlight",
      genotypes = gt2array(
        genotypes = tidy.data,
        n.ind = n.ind,
        n.snp = n.markers
      ),
      strata = strata,
      biallelic = TRUE,
      markers.meta = markers,
      filename = NULL,
      verbose = verbose
    )
    # if (verbose) message("Written: GDS filename: ", gds.filename)
  }# End gds genlight

  if (requested.tidy) return(tidy.data)
  if (gds) return(gds.filename)
  return(original.data)
} # End read_genlight

# write_genlight ----------------------------------------------------------------
#' @name write_genlight
#' @title Write a genlight object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Reading and writing \code{genlight} objects requires the optional CRAN
#' package \href{https://cran.r-project.org/package=adegenet}{\pkg{adegenet}}.
#' @description Write a \code{genlight} object from a tidy data frame or GDS file or object.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and might be of interest for users. \code{genlight} is a formal (S4) class
#' for storing genotypes of binary SNPs in a compact way, using a bit-level
#' coding scheme. This storage is most efficient with haploid data,
#' where the memory taken to represent data can be reduced more than 50 times.
#' However, \code{genlight} can be used for any level of ploidy,
#' and still remain an efficient storage mode.

#' @inheritParams genometranslator_common_arguments

#' @param write (logical, optional) To write in the working directory the genlight
#' object. The file is written with \code{radiator_genlight_DATE@TIME.RData} and
#' can be open with load or readRDS.
#' Default: \code{write = FALSE}.

#' @param dartr (logical, optional) For non-dartR users who wants to have a genlight
#' object ready for the dartR package. This option transfer or generates:
#' \code{CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG,
#' ONE_RATIO_REF, ONE_RATIO_SNP}. These markers metadata are stored into
#' the genlight slot: \code{genlight.obj@other$loc.metrics}.
#' \strong{Use the radiator generated GDS data for best result}.
#' Default: \code{dartr = FALSE}.


#' @param biallelic (logical, optional) If you already know that the data is
#' biallelic use this argument to speed up the function.
#' Default: \code{biallelic = TRUE}.

#' @param parallel.core Number of workers available for parallel operations.
#' Default: \code{parallel.core = parallel::detectCores() - 2}.
#' 
#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname write_genlight

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.

#' @examples
#' \dontrun{
#' # With defaults:
#' turtle <- genometranslator::write_genlight(data = "my.metadata.node.rad")
#'
#' # Write gl object in directory:
#' turtle <- genometranslator::write_genlight(data = "my.metadata.node.rad", write = TRUE)
#'
#' # Generate a dartR ready genlight and verbose = TRUE:
#' turtle <- genometranslator::write_genlight(
#'     data = "my.metadata.node.rad",
#'     write = TRUE,
#'     dartr = TRUE,
#'     verbose = TRUE
#'  )
#' }

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering


write_genlight <- function(
    data,
    write = FALSE,
    dartr = FALSE,
    verbose = FALSE,
    parallel.core = parallel::detectCores() - 2,
    biallelic = TRUE
) {


  # NOTE: Make sure it's the same levels.... id, markers etc...
  # TEST
  # biallelic = TRUE
  # write = TRUE
  # dartr = TRUE
  # verbose = TRUE
  # parallel.core = parallel::detectCores() - 1

  # Checking for missing and/or default arguments ------------------------------
  tgbase::check_package(package = "adegenet")
  if (missing(data)) rlang::abort("Input file missing")
  if (verbose) message("Generating genlight...")

  # File type detection---------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------
  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    # Check that SeqArray is installed (it requires automatically: SeqArray and gdsfmt)
    tgbase::check_package(package = "SeqArray", cran = FALSE, bioc = TRUE)

    if (data.type == "gds.file") {
      data <- genometranslator::read_genome(data, verbose = verbose)
    }
    biallelic <- detect_biallelic_markers(data)# faster with GDS
    if (!biallelic) rlang::abort("genlight object requires biallelic genotypes")

    data.source <- genometranslator::extract_data_source(data)
    # markers.meta <- extract_markers_metadata(gds = data, whitelist = TRUE)
    # data.bk <- data
    # data.bk -> data


    # no longer needed
    # if (dartr && !any(unique(c("1row", "2rows") %in% data.source))) {
      # counts data and data with read depth alleles depth info...

      # if (!"counts" %in% data.source) {
      #   genotypes.meta <- genometranslator::extract_genotypes_metadata(gds = data, whitelist = TRUE)
      # } else {
      #   genotypes.meta <- NULL
      # }

      data <- tidy_genome(
        data = data,
        wide = FALSE,
        pop.id = FALSE,
        parallel.core = parallel.core,
        calibrate.alleles = FALSE
      )
      markers.levels <- unique(data$MARKERS)
      ind.levels <- unique(data$INDIVIDUALS)

      # genotypes.meta <- NULL
      want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT",
                "CALL_RATE",
                "AVG_COUNT_REF",
                "AVG_COUNT_SNP",
                "REP_AVG",
                "ONE_RATIO_REF",
                "ONE_RATIO_SNP"
                )
      markers.meta <- data %>%
        dplyr::select(tidyselect::any_of(want)) %>%
        dplyr::distinct(MARKERS, .keep_all = TRUE)

      data %<>%
        tidyr::pivot_wider(
          id_cols = c(INDIVIDUALS, STRATA),
          names_from = MARKERS,
          values_from = ALT_DOSAGE
        )
    # }

    # no longer needed
    # else {
    #   markers.meta <- genometranslator::extract_markers_metadata(
    #     gds = data
    #   )
    #   data.check <- genometranslator::extract_individuals_metadata(
    #     gds = data, ind.field.select = c("INDIVIDUALS", "STRATA"), whitelist = TRUE) %>%
    #     dplyr::bind_cols(
    #       SeqArray::seqGetData(
    #         gdsfile = data, var.name = "$dosage_alt") %>%
    #         magrittr::set_colnames(x = ., value = markers.meta$MARKERS) %>%
    #         tibble::as_tibble(.)
    #     ) %>%
    #     dplyr::rename(POP_ID = STRATA)
    # }

    # dartR-----------------------------------------------------------------------
    if (dartr) {
      # That bits of code below generate what's necessary for dartR
      if (verbose) message("Calculating read depth for each SNP\n")
      if ("dart" %in% data.source) {
        # if (any(unique(c("1row", "2rows") %in% data.source))) {
        markers.meta %<>%
          dplyr::mutate(
            N_IND = SeqArray::seqApply(
              gdsfile = data,
              var.name = "$dosage_alt",
              FUN = function(g) length(g[!is.na(g)]),
              margin = "by.variant", as.is = "integer",
              parallel = parallel.core),
            rdepth = round(
              ((N_IND * ONE_RATIO_REF * AVG_COUNT_REF) +
                 (N_IND * ONE_RATIO_SNP * AVG_COUNT_SNP)) / N_IND
            )
          ) %>%
          dplyr::ungroup(.)
        # data.bk <- NULL
      } else {
        if (rlang::has_name(data, "READ_DEPTH")) {
          # dart 2 rows counts...
          if (rlang::has_name(markers.meta, "AVG_COUNT_REF")) {

            #NOTE TO MYSELF: will have to keep REP_AVG from count somehow... to do
            # for now, I expect people will use this with filtered data so not important
            # to fill with 1
            if ("counts" %in% data.source) rep.avg <- markers.meta$REP_AVG
            not.wanted <- c("CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP",
                            "REP_AVG", "ONE_RATIO_REF", "ONE_RATIO_SNP")
            markers.meta %<>% dplyr::select(-tidyselect::any_of(not.wanted))
          }

          suppressWarnings(
            markers.meta %<>%
              dplyr::left_join(
                data %>%
                  dplyr::group_by(MARKERS) %>%
                  dplyr::summarise(
                    rdepth = mean(READ_DEPTH, na.rm = TRUE),
                    AVG_COUNT_REF = mean(ALLELE_REF_DEPTH, na.rm = TRUE),
                    AVG_COUNT_SNP = mean(ALLELE_ALT_DEPTH, na.rm = TRUE),
                    CALL_RATE = length(INDIVIDUALS[!is.na(ALT_DOSAGE)]) / length(INDIVIDUALS),
                    ONE_RATIO_REF = length(INDIVIDUALS[ALT_DOSAGE == 0]) + length(INDIVIDUALS[ALT_DOSAGE == 1]) / length(INDIVIDUALS),
                    ONE_RATIO_SNP = length(INDIVIDUALS[ALT_DOSAGE == 2]) + length(INDIVIDUALS[ALT_DOSAGE == 1]) / length(INDIVIDUALS)
                  ) %>%
                  dplyr::mutate(REP_AVG = 1L)
                , by = "MARKERS"
              ) %>%
              dplyr::ungroup(.)
          )
          genotypes.meta <- NULL
          #Note to myself: will have to remove duplicated info in code between genotypes.meta and data...
          if ("counts" %in% data.source) {
            markers.meta$REP_AVG <- rep.avg
          }
        }
      }

      # ~25 times slower
      # Calculate Read Depth (from Arthur Georges)
      # gl.obj@other$loc.metrics$rdepth <- array(NA, adegenet::nLoc(gl.obj))
      # for (i in 1:adegenet::nLoc(gl.obj)){
      #   called.ind <- round(adegenet::nInd(gl.obj) * gl.obj@other$loc.metrics$CallRate[i],0)
      #   ref.count <- called.ind * gl.obj@other$loc.metrics$OneRatioRef[i]
      #   alt.count <- called.ind * gl.obj@other$loc.metrics$OneRatioSnp[i]
      #   sum.count.ref <- ref.count * gl.obj@other$loc.metrics$AvgCountRef[i]
      #   sum.count.alt <- ref.count * gl.obj@other$loc.metrics$AvgCountSnp[i]
      #   gl.obj@other$loc.metrics$rdepth[i] <- round((sum.count.alt + sum.count.ref) / called.ind, 1)
      # }

      # gl.obj@other$loc.metrics$rdepth
      data %<>%
        dplyr::select(MARKERS, STRATA, INDIVIDUALS, ALT_DOSAGE) %>%
        dplyr::mutate(ALT_DOSAGE = as.integer(ALT_DOSAGE)) %>%
        tgbase::trans_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS",
          values_from = "ALT_DOSAGE"
        )
    }#End dartr

    data.type <- "tbl_df"
  } else {
    # Tidy data
    # if (rlang::has_name(data, "STRATA")) data %<>% dplyr::rename(POP_ID = STRATA)
    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "STRATA", "INDIVIDUALS",
              "REF", "ALT", "GT_VCF", "ALT_DOSAGE",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "ONE_RATIO_REF", "ONE_RATIO_SNP")
    data %<>%
      genometranslator::read_genome(data = ., import.metadata = TRUE) %>%
      dplyr::select(tidyselect::any_of(want)) %>%
      dplyr::arrange(STRATA, INDIVIDUALS)

    # Detect if biallelic data ---------------------------------------------------
    if (is.null(biallelic)) biallelic <- detect_biallelic_markers(data)
    if (!biallelic) rlang::abort("genlight object requires biallelic genotypes")

    want <- c("MARKERS", "CHROM", "LOCUS", "POS", "REF", "ALT",
              "CALL_RATE", "AVG_COUNT_REF", "AVG_COUNT_SNP", "REP_AVG",
              "ONE_RATIO_REF", "ONE_RATIO_SNP")
    data %<>% dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)
    markers.meta <- dplyr::select(data, tidyselect::any_of(want)) %>%
      dplyr::distinct(MARKERS, .keep_all = TRUE) %>%
      separate_markers(
        data = .,
        sep = "__",
        markers.meta.all.only = TRUE,
        biallelic = TRUE,
        verbose = verbose)

    if (!rlang::has_name(data, "ALT_DOSAGE") && rlang::has_name(data, "GT_VCF")) {
      data$ALT_DOSAGE <- stringi::stri_replace_all_fixed(
        str = data$GT_VCF,
        pattern = c("0/0", "1/1", "0/1", "1/0", "./."),
        replacement = c("0", "2", "1", "1", NA),
        vectorize_all = FALSE
      )
    }

    data %<>%
      dplyr::select(MARKERS, STRATA, INDIVIDUALS, ALT_DOSAGE) %>%
      dplyr::mutate(ALT_DOSAGE = as.integer(ALT_DOSAGE)) %>%
      tgbase::trans_wide(
        x = .,
        formula = "STRATA + INDIVIDUALS ~ MARKERS",
        values_from = "ALT_DOSAGE"
      )
  }# End tidy data

  # Generate genlight
  parallel.core.temp <- FALSE
  if (length(markers.meta$MARKERS) > 10000) parallel.core.temp <- parallel.core

  # longer
  # gl.obj2 <- methods::new(
  #   "genlight",
  #   gen = data[,-(1:2)],
  #   ploidy = 2,
  #   ind.names = data$INDIVIDUALS,
  #   chromosome = markers.meta$CHROM,
  #   loc.names = markers.meta$LOCUS,
  #   position = markers.meta$POS,
  #   pop = data$POP_ID,
  #   parallel = parallel.core.temp)
  # tictoc::toc()


  gl.obj <- methods::new(
    "genlight",
    data[,-(1:2)],
    parallel = parallel.core.temp
  )
  adegenet::indNames(gl.obj)   <- data$INDIVIDUALS
  adegenet::pop(gl.obj)        <- data$STRATA
  adegenet::chromosome(gl.obj) <- markers.meta$CHROM
  adegenet::locNames(gl.obj)   <- markers.meta$LOCUS
  adegenet::position(gl.obj)   <- markers.meta$POS

  if (dartr) {
    gl.obj@other$loc.metrics <- markers.meta %>%
      dplyr::select(
        OneRatioRef = ONE_RATIO_REF,
        OneRatioSnp = ONE_RATIO_SNP,
        AvgCountRef = AVG_COUNT_REF,
        AvgCountSnp = AVG_COUNT_SNP,
        CallRate = CALL_RATE,
        RepAvg = REP_AVG,
        rdepth
      )
  }#End dartr


  # Check
  # gl.obj@n.loc
  # gl.obj@ind.names
  # gl.obj@chromosome
  # gl.obj@position
  # length(gl.obj@position)
  # gl.obj@loc.names
  # length(gl.obj@loc.names)
  # gl.obj@pop
  # gl.obj@strata
  # adegenet::nLoc(gl.obj)
  # adegenet::popNames(gl.obj)
  # adegenet::indNames(gl.obj)
  # adegenet::nPop(gl.obj)
  # adegenet::NA.posi(gl.obj)
  # names(gl.obj@other$loc.metrics)


  if (write) {
    filename.temp <- generate_filename(extension = "genlight")
    filename.short <- filename.temp$filename.short
    filename.genlight <- filename.temp$filename
    saveRDS(object = gl.obj, file = filename.genlight)
    if (verbose) message("File written: ", filename.short)
  }
  return(gl.obj)
} # End write_genlight
