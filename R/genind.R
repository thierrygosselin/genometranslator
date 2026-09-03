# read_genind ------------------------------------------------------------------
#' @name read_genind
#' @title Read a genind object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Reading and writing \code{genind} objects requires the optional CRAN package
#' \href{https://cran.r-project.org/package=adegenet}{\pkg{adegenet}}.
#' @description Import an \pkg{adegenet} \code{genind} object, standardize its
#' sample and locus information, and optionally produce tidy genomic data or a
#' GDS representation. Diploid biallelic and multiallelic objects are supported;
#' GDS generation is currently limited to biallelic data.

#' @param data (path or object) A genind object in the global environment or
#' path to a genind file that will be open with \code{readRDS}.

#' @param tidy (logical) Generate a tidy dataset.
#' Default: \code{tidy = TRUE}.

#' @param gds Logical. Generate a genometranslator GDS representation.
#' Currently, for biallelic datasets only.
#' Default: \code{gds = TRUE}.

#' @param write Logical. Write the tidy table as an Arrow/Parquet file.
#' Default: \code{write = FALSE}.

#' @note \href{https://github.com/thibautjombart/adegenet}{genind} objects, like
#' genepop, are not optimal genomic format for RADseq datasets,
#' they lack important genotypes and markers metadata: chromosome, locus, snp,
#' position, read depth, allele depth, etc.
#' \href{https://github.com/thibautjombart/adegenet}{genlight} object is a more
#' interesting container and is memory efficient, see \code{\link{read_genlight}}.
#'
#'
#' By default allele names will be kept for the tidy dataset,
#' if the alleles is numeric and length < 3.
#'
#'
#' In the unlikely event that the genind object as no stratification/population,
#' \emph{pop} will be added to the strata column.


#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_genind
#' @return With \code{tidy = TRUE}, a tidy tibble. With
#' \code{tidy = FALSE, gds = TRUE}, the generated GDS result. With both options
#' false, the original \code{genind} object. Multiallelic input disables GDS
#' generation and returns tidy data when requested.
#' @examples
#' \dontrun{
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'   data("nancycats", package = "adegenet")
#'   cats <- genometranslator::read_genind(nancycats, gds = FALSE)
#' }
#' }
#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies


read_genind <- function(
    data,
    tidy = TRUE,
    gds = TRUE,
    write = FALSE,
    verbose = FALSE
) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "read_genind",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  # TEST
  # tidy = TRUE
  # gds = TRUE
  # write = FALSE
  # verbose = TRUE
  if (verbose) cli::cli_progress_step("Reading genind")

  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("A genind object or RDS file is required.")
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) rlang::abort("The genind RDS file does not exist.")
    data <- readRDS(data)
  }
  if (!inherits(data, "genind")) rlang::abort("Input is not a genind object.")
  original.data <- data
  requested.tidy <- isTRUE(tidy)


  # Working on individuals and strata ------------------------------------------
  strata <- tibble::tibble(INDIVIDUALS = rownames(data@tab))
  if (is.null(data@pop)) {
    strata %<>% dplyr::mutate(STRATA = "pop")
  } else {
    strata$STRATA = data@pop
  }

  n.ind <- length(strata$INDIVIDUALS)
  n.snp <- length(unique(data@loc.fac))

  biallelic <- max(unique(data@loc.n.all)) == 2
  if (!biallelic) gds <- FALSE
  if (gds) tidy <- TRUE

  if (write) {
    filename.temp <- generate_filename(extension = "arrow.parquet")
    filename.short <- filename.temp$filename.short
    filename.genind <- filename.temp$filename
  }

  if (tidy) {
    # detect A2... facilitate the conversion if made by radiator...
    A2 <- FALSE
    codom <- FALSE
    colnames.coding <- sample(colnames(data@tab), min(ncol(data@tab), 100))
    A2 <- TRUE %in% unique(stringi::stri_detect_fixed(str = colnames.coding, pattern = ".A2"))
    if (data@type == "codom") codom <- TRUE

    if (biallelic) {
      if (verbose) cli::cli_progress_step("Preparing biallelic tidy data")
      markers.meta <- colnames(data@tab)
      n.snp <- length(markers.meta)/2

      alt.alleles <- tibble::tibble(
        MARKERS_ALLELES = markers.meta,
        COUNT = colSums(x = data@tab, na.rm = TRUE)
      ) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_extract_first_regex(str = markers.meta, pattern = "^[^.]+"),
          ALLELES = stringi::stri_extract_last_regex(str = markers.meta, pattern = "(?<=\\.).*"),
          ALLELES = stringi::stri_replace_all_fixed(str = ALLELES, pattern = c("A1__", "A2__"), replacement = c("", ""), vectorize_all = FALSE),
        ) %>%
        dplyr::arrange(MARKERS, COUNT, ALLELES) %>%
        dplyr::mutate(REF = rep(c("ALT", "REF"), n() / 2))

      # check that REF/ALT are A, C, G, T before going further .....
      discard.alleles <- FALSE
      unique.alleles <- unique(alt.alleles$ALLELES)
      if (length(unique.alleles) > 4) discard.alleles <- TRUE
      if (any(stringi::stri_detect_regex(str = unique.alleles, pattern = "[0-9]"))) discard.alleles <- TRUE

      if (!discard.alleles) {
        ref.alt <- dplyr::select(alt.alleles, MARKERS, REF, ALLELES) %>%
          tgbase::trans_wide(x = ., formula = "MARKERS ~ REF", values_from = "ALLELES")
      }

      alt.alleles %<>% dplyr::filter(REF == "ALT") %$% MARKERS_ALLELES
      if (length(alt.alleles) != n.snp) rlang::abort("Contact author problem with tidying the genind")

      data <- tibble::as_tibble(t(data@tab), rownames = "MARKERS") %>%
        dplyr::filter(MARKERS %in% alt.alleles) %>%
        dplyr::mutate(
          MARKERS = stringi::stri_extract_first_regex(str = MARKERS, pattern = "^[^.]+"),
          MARKERS = stringi::stri_replace_all_fixed(
            str = MARKERS,
            pattern = c("__A1", "__A2"),
            replacement = c("", ""),
            vectorize_all = FALSE
          ),
          VARIANT_ID = as.integer(factor(MARKERS))) %>%
        dplyr::arrange(VARIANT_ID)

      alt.alleles <- NULL
      markers.meta <- dplyr::distinct(data, VARIANT_ID, MARKERS)

      if (!discard.alleles) {
        markers.meta %<>%
          dplyr::left_join(ref.alt, by = "MARKERS")
      }
      ref.alt <- NULL

      # generate markers metadata
      markers.meta %<>% separate_markers(data = ., generate.ref.alt = TRUE, biallelic = TRUE)

      data <- tgbase::trans_long(
        x = data,
        cols = c("MARKERS", "VARIANT_ID"),
        names_to = "INDIVIDUALS",
        values_to = "ALT_DOSAGE"
      ) %>%
        dplyr::left_join(strata, by = "INDIVIDUALS") %>%
        dplyr::left_join(markers.meta, by = c("MARKERS", "VARIANT_ID")) %>%
        dplyr::mutate(
          INDIVIDUALS = clean_ind_names(INDIVIDUALS),
          STRATA = clean_pop_names(STRATA)
        ) %>%
        dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)

    } else {
      if (verbose) cli::cli_progress_step("Preparing multi-allelic tidy data")
      data <- tgbase::trans_long(
          x = tibble::as_tibble(data@tab, rownames = "INDIVIDUALS"),
          cols = "INDIVIDUALS",
          names_to = c("MARKERS", "ALLELES"),
          values_to = "COUNT",
          names_sep = "\\.",
          tidy = TRUE
        ) %>%
        dplyr::left_join(strata, by = "INDIVIDUALS")

      numeric.alleles <- stringi::stri_detect_regex(
        str = unique(data$ALLELES),
        pattern = "^[0-9]+$"
      )
      if (!all(numeric.alleles)) {
        data %<>% dplyr::mutate(ALLELES = as.numeric(factor(ALLELES)))
      }
      data %<>%
        dplyr::mutate(ALLELES = stringi::stri_pad_left(str = ALLELES, pad = "0", width = 3)) %>%
        dplyr::filter(is.na(COUNT) | COUNT != 0)

      multi_genind <- function(x) {
        count.type <- unique(x$COUNT)
        if (is.na(count.type)) {
          x %<>%
            dplyr::distinct(STRATA, INDIVIDUALS, MARKERS) %>%
            dplyr::mutate(GT = rep("000000", dplyr::n())) %>%
            dplyr::ungroup(.)
        } else if (count.type == 2L) {
          x %<>%
            dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
            dplyr::summarise(GT = stringi::stri_join(ALLELES, ALLELES, sep = ""), .groups = "drop")
        } else {
          x %<>%
            dplyr::group_by(STRATA, INDIVIDUALS, MARKERS) %>%
            dplyr::summarise(GT = stringi::stri_join(ALLELES, collapse = ""), .groups = "drop")
        }
      }

      if (n.ind * n.snp >= 5000000) {
        data <- tgbase::parallel_map(
          .x = data,
          .f = multi_genind,
          flat.future = "dfr",
          split.with = "COUNT"
        )
      } else {
        data %<>%
          dplyr::group_split(COUNT) %>%
          purrr::map_dfr(.x = ., .f = multi_genind)
      }
    }#End for multi-allelic or weird genind

    if (write) genometranslator::write_genome(data = data, filename = filename.genind, verbose = verbose)

  }# End tidy genind
  if (gds) {
    if (verbose) cli::cli_progress_step("Working on the GDS dataset")
    if (!rlang::has_name(x = data, "ALT_DOSAGE")) rlang::abort("Missing ALT_DOSAGE format to generate GDS, contact author")
    n.snp <- dplyr::n_distinct(data$MARKERS)
    n.ind <- dplyr::n_distinct(data$INDIVIDUALS)

    # generate genotypes format for easy reading into GDS
    data %<>%
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
      data.source = "genind",
      genotypes = gt2array(
        genotypes = data,
        n.ind = n.ind,
        n.snp = n.snp
      ),
      strata = strata,
      biallelic = TRUE,
      markers.meta = markers.meta,
      filename = NULL,
      verbose = verbose
    )
  }# End gds genind

  if (requested.tidy) return(data)
  if (gds) return(gds.filename)
  return(original.data)
} # End read_genind


# write_genind ------------------------------------------------------------------

#' @name write_genind
#' @title Write a genind object
#' @section Dependencies:
#' Required package dependencies are declared in DESCRIPTION and installed
#' with genometranslator. Run `genometranslator_dependencies()` to inspect
#' core packages, optional packages, and external executables.
#'
#' Reading and writing \code{genind} objects requires the optional CRAN package
#' \href{https://cran.r-project.org/package=adegenet}{\pkg{adegenet}}.

#' @description Write a genind object from a tidy data frame.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @inheritParams genometranslator_common_arguments

#' @param write (logical, optional) To write in the working directory the genind
#' object. The file is written with \code{radiator_genind_DATE@TIME.RData} and
#' can be open with load or readRDS.
#' Default: \code{write = FALSE}.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname write_genind

#' @references Jombart T (2008) adegenet: a R package for the multivariate
#' analysis of genetic markers. Bioinformatics, 24, 1403-1405.
#' @references Jombart T, Ahmed I (2011) adegenet 1.3-1:
#' new tools for the analysis of genome-wide SNP data.
#' Bioinformatics, 27, 3070-3071.


#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering


write_genind <- function(data, write = FALSE, verbose = FALSE) {


  # TEST
  # write = TRUE
  # verbose = TRUE


  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)


  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") data %<>% genometranslator::read_genome(data = ., verbose = verbose)
    data <- tidy_genome(data = data, parallel.core = parallel::detectCores() - 1, pop.id = FALSE)
    data.type <- "tbl_df"
  } else {
    want <- c("MARKERS", "STRATA", "INDIVIDUALS", "REF", "ALT", "GT", "ALT_DOSAGE")
    data %<>%
      genometranslator::read_genome(data = ., import.metadata = TRUE) %>%
      dplyr::select(tidyselect::any_of(want))
  }

  #
  # Check strata and pop.levels
  pop.levels <- unique(data$STRATA)
  if (rlang::has_name(data, "STRATA") && is.factor(data$STRATA)) {
    if (dplyr::n_distinct(data$STRATA) != nlevels(data$STRATA)) {
      data$STRATA <- droplevels(data$STRATA)
      pop.levels <- levels(data$STRATA)
    }
    pop.levels <- levels(data$STRATA)
  }

  # Make sure that STRATA and INDIVIDUALS are character
  data$INDIVIDUALS <- as.character(data$INDIVIDUALS)
  data$STRATA <- as.character(data$STRATA)

  # Isolate the strata
  # we convert STRATA to factor because adegenet does it automatically...
  pop.num <- unique(stringi::stri_detect_regex(str = pop.levels, pattern = "^[0-9]+$"))
  if (length(pop.num) == 1 && pop.num) pop.levels <- as.character(sort(as.numeric(pop.levels)))

  # When ALT_DOSAGE available
  if (rlang::has_name(data, "ALT_DOSAGE")) {
    if (rlang::has_name(data, "REF")) {
      data <- dplyr::bind_rows(
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, REF, n = ALT_DOSAGE) %>%
          dplyr::mutate(
            REF = stringi::stri_join("A1", REF, sep = "__"),
            MARKERS_ALLELES = stringi::stri_join(MARKERS, REF, sep = "."),
            MARKERS = NULL,
            REF = NULL,
            n = as.integer(abs(n - 2))
          ),
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, ALT, n = ALT_DOSAGE) %>%
          dplyr::mutate(
            ALT = stringi::stri_join("A2", ALT, sep = "__"),
            MARKERS_ALLELES = stringi::stri_join(MARKERS, ALT, sep = "."),
            MARKERS = NULL,
            ALT = NULL
          )
      ) %>%
        tgbase::trans_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS_ALLELES",
          values_from = "n"
        ) %>%
        dplyr::mutate(STRATA = factor(as.character(STRATA), levels = pop.levels)) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(STRATA, INDIVIDUALS)
    } else {
      data <- dplyr::bind_rows(
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, n = ALT_DOSAGE) %>%
          dplyr::mutate(
            MARKERS_ALLELES = stringi::stri_join(MARKERS, "A1", sep = "."),
            MARKERS = NULL,
            REF = NULL,
            n = as.integer(abs(n - 2))
          ),
        dplyr::select(data, MARKERS, STRATA, INDIVIDUALS, n = ALT_DOSAGE) %>%
          dplyr::mutate(
            MARKERS_ALLELES = stringi::stri_join(MARKERS, "A2", sep = "."),
            MARKERS = NULL,
            ALT = NULL
          )
      ) %>%
        tgbase::trans_wide(
          x = .,
          formula = "STRATA + INDIVIDUALS ~ MARKERS_ALLELES",
          values_from = "n"
        ) %>%
        dplyr::mutate(STRATA = factor(as.character(STRATA), levels = pop.levels)) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(STRATA, INDIVIDUALS)
    }
  } else {
    missing.geno <- dplyr::ungroup(data) %>%
      dplyr::select(MARKERS, INDIVIDUALS, GT) %>%
      dplyr::filter(GT == "000000") %>%
      dplyr::select(MARKERS, INDIVIDUALS) #%>%  dplyr::mutate(MISSING = rep("blacklist", n()))

    data <- suppressWarnings(
      dplyr::ungroup(data) %>%
        dplyr::select(MARKERS, INDIVIDUALS, GT, STRATA) %>%
        dplyr::filter(GT != "000000") %>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6),
          GT = NULL
        ) %>%
        tgbase::trans_long(
          x = .,
          cols = c("INDIVIDUALS", "STRATA", "MARKERS"),
          names_to = "ALLELES",
          values_to = "GT"
        ) %>%
        dplyr::arrange(MARKERS, STRATA, INDIVIDUALS, GT) %>%
        dplyr::count(x = ., STRATA, INDIVIDUALS, MARKERS, GT) %>%
        dplyr::ungroup(.) %>%
        tidyr::complete(data = ., tidyr::nesting(INDIVIDUALS, STRATA), tidyr::nesting(MARKERS, GT), fill = list(n = 0)) %>%
        dplyr::mutate(MARKERS_ALLELES = stringi::stri_join(MARKERS, GT, sep = ".")) %>%
        dplyr::anti_join(missing.geno, by = c("MARKERS", "INDIVIDUALS")) %>%
        dplyr::select(-MARKERS, -GT) %>%
        dplyr::mutate(STRATA = factor(as.character(STRATA), levels = pop.levels)) %>%# xvalDapc doesn't accept pop as ordered factor
        dplyr::arrange(MARKERS_ALLELES, INDIVIDUALS) %>%
        tgbase::trans_wide(x = ., formula = "STRATA + INDIVIDUALS ~ MARKERS_ALLELES", values_from = "n") %>%
        dplyr::arrange(STRATA, INDIVIDUALS))
  }

  strata.genind <- dplyr::distinct(.data = data, INDIVIDUALS, STRATA) %>%
    dplyr::mutate(INDIVIDUALS = factor(INDIVIDUALS, levels = unique(data$INDIVIDUALS)))

  # genind arguments common to all data.type

  ind <- data$INDIVIDUALS
  pop <- data$STRATA

  data <- dplyr::ungroup(data) %>%
    dplyr::select(-c(INDIVIDUALS, STRATA))

  suppressWarnings(rownames(data) <- ind)

  # Remove individuals and loci entirely missing ---------------------------------
  data.matrix <- as.matrix(data)

  all.na.ind <- rowSums(is.na(data.matrix)) == ncol(data.matrix)
  all.na.loc <- colSums(is.na(data.matrix)) == nrow(data.matrix)

  if (any(all.na.ind)) {
    message(
      "number of individuals missing at all loci and removed: ",
      sum(all.na.ind)
    )
    data.matrix <- data.matrix[!all.na.ind, , drop = FALSE]
    pop <- pop[!all.na.ind]
    ind <- ind[!all.na.ind]
    strata.genind <- strata.genind %>%
      dplyr::filter(INDIVIDUALS %in% ind)

    strata.genind$INDIVIDUALS <- factor(strata.genind$INDIVIDUALS, levels = ind)

    strata.genind <- dplyr::arrange(strata.genind, INDIVIDUALS)

    }

  if (any(all.na.loc)) {
    message(
      "number of loci missing in all individuals and removed: ",
      sum(all.na.loc)
    )
    data.matrix <- data.matrix[, !all.na.loc, drop = FALSE]
  }

  data <- tibble::as_tibble(data.matrix)
  suppressWarnings(rownames(data) <- ind)

   if (nrow(data) == 0L) rlang::abort("No individuals left after removing all-NA individuals")
  if (ncol(data) == 0L) rlang::abort("No loci left after removing all-NA loci")

  # genind constructor
  prevcall <- match.call()
  res <- adegenet::genind(
    tab = data,
    pop = pop,
    prevcall = prevcall,
    ploidy = 2,
    type = "codom",
    strata = strata.genind,
    hierarchy = NULL
  )
  data <- strata.genind <- NULL

  if (write) {
    filename.temp <- generate_filename(extension = "genind")
    filename.short <- filename.temp$filename.short
    filename.genind <- filename.temp$filename
    saveRDS(object = res, file = filename.genind)
    if (verbose) message("File written: ", filename.short)
  }

  return(res)
} # End write_genind
