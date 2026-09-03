# read_genepop------------------------------------------------------------------
#' @name read_genepop

#' @title Read a Genepop file

#' @description Import and validate diploid data in the Genepop text format.
#' Both comma-separated locus names and one-locus-per-line headers are accepted,
#' as are one-, two-, and three-digit allele codes. Population blocks, unique
#' locus names, genotype-row widths, and genotype coding are checked before a
#' long or wide tibble is returned.

#' @param data Path to a Genepop file, commonly with extension \code{.gen}, or a
#' one-column data frame containing its lines.

#' @param strata (optional) A tab delimited file with 2 columns. Header:
#' \code{INDIVIDUALS} and \code{STRATA}.
#' The \code{STRATA} column can be any hierarchical grouping.
#' To create a strata file see \code{\link{individuals2strata}}.
#' Default: \code{strata = NULL}.

#' @param tidy (optional, logical) With \code{tidy = FALSE},
#' the markers are the variables and the genotypes the observations (wide format).
#' With the default: \code{tidy = TRUE}, markers and genotypes are variables
#' with their own columns (long format).

#' Default: \code{tidy = TRUE}.
#' @param filename (optional) The file name for the tidy data frame
#' written to the working directory.
#' With the default, The tidy data is
#' in the global environment only (i.e. not written in the working directory).

#' Default: \code{filename = NULL}.
#' @return A tibble in long form with \code{STRATA}, \code{INDIVIDUALS},
#' \code{MARKERS}, and six-digit \code{GT}, or wide form when
#' \code{tidy = FALSE}. If duplicate individual names are replaced, the
#' original-to-generated mapping is stored in the \code{id_conversion}
#' attribute. If \code{filename} is supplied, the table is also written as a
#' tab-separated file.
#'
#' @section Field handling:
#' Genepop population blocks become \code{STRATA}; cleaned sample names become
#' \code{INDIVIDUALS}; locus names become \code{MARKERS}; and each allele is
#' padded to three digits. Supplied strata metadata replaces embedded population
#' blocks and must contain every Genepop individual.

#' @details \href{https://genepop.curtin.edu.au/help_input.html}{genepop format}
#' \enumerate{
#' \item \strong{First line:} This line is used to store information about
#' your data, any characters are allowed.
#' This line is not kept inside \code{read_genepop}.
#' \item \strong{Second line:} 2 options i) wide: all the locus are stored on this same
#' line with comma (or a comma+space) separator;
#' ii) long: the name of the first locus
#' (the remaining locus are on separate and subsequent rows with the long format:
#' not recommended with genomic datasets with thousands of markers...).
#'
#' The remaining lines are blocks of population and genotypes,
#' \href{https://genepop.curtin.edu.au/help_input.html}{for genepop format examples}.
#' \item \strong{population identifier:} The population block are separated by
#' the word: \code{POP}, or \code{Pop} or \code{pop}. Flavors of the genepop
#' software uses the first or the last identifier of every sub-population,
#' in all output files, to name populations
#' \href{https://genepop.curtin.edu.au/help_input.html}{(more info)}.
#' This is not very convenient for population naming and prone to errors,
#' this is where the \code{strata} argument inside \code{read_genepop}
#' (described above) becomes handy.
#' \item \strong{individual identifier:} After the population identifier
#' the individuals that belong to the same population are found subsequently on
#' separate lines. The genepop format specify you can use any character,
#' including a blank space or tab. Spaces are allowed in the identifier names.
#' You may leave it blank except for a comma if you wish. The comma between
#' the individual identifier and the list of genotypes is required.
#' For good naming habit however the function \code{read_genepop} will
#' replace \code{"_", ":"} with \code{"-"}, the comma \code{","} along any white
#' space characters defined as \code{"\t", "\n", "\f", "\r", "\p{Z}"} found in the individual
#' name will be trimmed.
#' \item \strong{genotypes:} For each locus, genotypes are separated by one or
#' more blank spaces or tab. 0101 indicates that this individual is homozygous
#' for the 01 allele at the first locus.
#' An alternative input format exists, where each allele is coded by three digits
#' (instead of two as described above). However, the total number of
#' different alleles, for each locus, should not be higher than 99.
#' Missing data is coded with zeros: \code{0000 or 000000}.
#' }

#' @note
#' \href{https://genepop.curtin.edu.au/help_input.html}{genepop format notes:}
#' \itemize{
#' \item No constraint on blanks separating the various fields.
#' \item tabs or spaces allowed.
#' \item Loci names can appear on separate lines (long), or on one line (wide)
#' if separated by commas.
#' \item Individual identifier may have blanks but must end with a comma.
#' \item Alleles are numbered from 01 to 99 (or 001 to 999).
#' Consecutive numbers to designate alleles are not required.
#' \item Populations are defined by the position of the "Pop" separator.
#' To group various populations, just remove relevant "Pop" separators.
#' \item Missing data should be indicated as 00 (or 000) rather than blanks.
#' There are three possibilities for missing data :
#' no information (0000) or (000000), partial information
#' for first allele (1000) or (010000), partial information
#' for second allele (0010) or (000010).
#' \item The number of locus names should correspond to the number of genotypes
#' in each row.
#' If you remove one or several loci from your input file,
#' you should remove both their names and the corresponding genotypes.
#' \item No empty lines should be found within the file.
#' \item No more than one empty line should be present at the end of file.
#' }
#' \strong{not an ideal genomic format: } The nice thing about RADseq dataset is
#' that you have several important genotypes and markers metadata
#' (chromosome, locus, snp, position, read depth, allele depth, etc.) available,
#' these are all lacking in the genepop format. This format is kept for archival
#' reasons and data exchange. VCF or GDS should be preferred when sequencing
#' likelihoods, depth, alleles, coordinates, and other metadata must be retained.

#' @param verbose Logical indicating whether progress messages are emitted.
#' Default: \code{verbose = FALSE}.
#' 
#' @export
#' @rdname read_genepop
#' @examples
#' \dontrun{
#' # We will use the genepop dataset provided with adegenet package
#' if (requireNamespace("adegenet", quietly = TRUE)) {
#'
#' # The simplest form of the function:
#' nancycats.tidy <- genometranslator::read_genepop(
#'     data = system.file(
#'         "files/nancycats.gen",
#'         package = "adegenet"
#'         )
#'     )
#'
#' # To output a data frame in wide format, with markers in separate columns:
#' nancycats.wide <- genometranslator::read_genepop(
#'     data = system.file(
#'     "files/nancycats.gen",
#'     package="adegenet"
#' ),
#'     tidy = FALSE
#' )
#' }
#' }


#' @references Raymond M. & Rousset F, (1995).
#' GENEPOP (version 1.2): population genetics software for exact tests
#' and ecumenicism.
#' J. Heredity, 86:248-249
#' @references Rousset F.
#' genepop'007: a complete re-implementation of the genepop software
#' for Windows and Linux.
#' Molecular Ecology Resources.
#' 2008, 8: 103-106.
#' doi:10.1111/j.1471-8286.2007.01931.x

#' @seealso \href{https://genepop.curtin.edu.au}{genepop}

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template io-dependencies


read_genepop <- function(
    data,
    strata = NULL,
    tidy = TRUE,
    filename = NULL,
    verbose = FALSE
) {

  # Common startup -------------------------------------------------------------
  .start <- tgbase::startup(
    package = "genometranslator",
    f.name = "read_genepop",
    verbose = verbose
  )
  on.exit(tgbase::teardown(.start), add = TRUE)

  # Checking for missing and/or default arguments-------------------------------
  if (missing(data)) rlang::abort("A Genepop file or table is required.")

  # Import data ------------------------------------------------------------------
  if (is.character(data) && length(data) == 1L) {
    if (!file.exists(data)) rlang::abort("The Genepop file does not exist.")
    data <- readr::read_delim(
      file = data,
      delim = "?",
      skip = 1,
      trim_ws = TRUE,
      col_names = "data",
      col_types = "c")
  } else if (inherits(data, "data.frame")) {
    if (ncol(data) != 1L) {
      rlang::abort("An in-memory Genepop table must contain exactly one column.")
    }
    names(data) <- "data"
    data <- data %>%
      dplyr::slice(-1) %>% # removes genepop header
      tibble::as_tibble()
  } else {
    rlang::abort("`data` must be one Genepop file path or a one-column table.")
  }

  # Replace white space with only 1 space
  data$data %<>%
    stringi::stri_replace_all_regex(
      str = .,
      pattern = "\\s+",
      replacement = " ",
      vectorize_all = FALSE
    ) %>%
    # Remove unnecessary spaces
    stringi::stri_trim_right(str = ., pattern = "\\P{Wspace}")

  # Pop indices ----------------------------------------------------------------
  pop.indices <- which(data$data %in% c("Pop", "pop", "POP"))
  npop <- length(pop.indices)
  if (!npop) rlang::abort("No Genepop population separator (`Pop`) was found.")
  if (pop.indices[[1]] == 1L) rlang::abort("No locus names were found before the first `Pop` separator.")

  # Markers --------------------------------------------------------------------
  # With this function, it doesn't matter if the markers are on one row or one per row.
  markers <- dplyr::slice(.data = data, 1:(pop.indices[1] - 1)) %>% purrr::flatten_chr(.x = .)
  markers <- unlist(stringi::stri_split_fixed(str = markers, pattern = ","))
  markers <- stringi::stri_replace_all_fixed(
    str = markers,
    pattern = " ",
    replacement = "",
    vectorize_all = FALSE
  )
  markers <- markers[nzchar(markers)]
  if (!length(markers)) rlang::abort("The Genepop file contains no locus names.")
  if (anyDuplicated(markers)) rlang::abort("Genepop locus names must be unique.")

  # Remove markers from the dataset --------------------------------------------
  data %<>% dplyr::slice(-(1:(pop.indices[1] - 1)))

  # New pop indices and pop string ---------------------------------------------
  pop.indices <- which(data$data %in% c("Pop", "pop", "POP"))
  pop.indices.lagged <- diff(c(pop.indices, (nrow(data) + 1))) - 1
  pop <- factor(rep.int(1:npop, times = pop.indices.lagged))

  # remove pop indices from dataset---------------------------------------------
  data %<>% dplyr::slice(-pop.indices)

  # Scan for genotypes split on 2 lines ----------------------------------------
  # looks for lines without a comma
  problem <- which(!seq_along(data$data) %in% grep(",", data$data, fixed = TRUE))
  if (length(problem) > 0) {
    for (i in sort(problem, decreasing = TRUE)) {
      if (i == 1L) {
        rlang::abort("A continued genotype row appears before any individual record.")
      }
      data$data[i - 1] <- paste(data$data[i - 1], data$data[i], sep = " ")
    }
    data <- dplyr::slice(.data = data, -problem)
    pop <- pop[-problem]
  }

  # preparing the dataset ------------------------------------------------------
  # separate the individuals based on the comma (",")
  data %<>%
    tidyr::separate(
      data = .,
      col = data,
      into = c("INDIVIDUALS", "GT"),
      sep = ","
    )

  # Individuals
  individuals <- dplyr::select(.data = data, INDIVIDUALS) %>%
    dplyr::mutate(INDIVIDUALS = clean_ind_names(x = INDIVIDUALS))

  # Check for duplicate individual names
  # some genepop format don't provide individuals
  id.conversion <- NULL
  if (anyDuplicated(individuals$INDIVIDUALS)) {
    rlang::warn(c(
      "Genepop individual identifiers are duplicated.",
      "i" = "Unique genometranslator identifiers were generated; inspect the `id_conversion` attribute."
    ))
    bad.id <- dplyr::select(individuals, BAD_ID = INDIVIDUALS) %>%
      dplyr::mutate(INDIVIDUALS = stringi::stri_join("genometranslator-individual-", seq_len(nrow(.)))) %>%
      dplyr::select(INDIVIDUALS, BAD_ID)

    individuals <- dplyr::select(bad.id, INDIVIDUALS)
    id.conversion <- bad.id
  }

  # isolate the genotypes
  data %<>%
    dplyr::select(GT) %>%
    dplyr::mutate(
      GT = stringi::stri_replace_all_fixed(
        str = GT,
        pattern = c("\t", ","), # remove potential comma and tab
        replacement = c(" ", ""),
        vectorize_all = FALSE
      ),
      # replace white space character: [\t\n\f\r\p{Z}]
      GT = stringi::stri_replace_all_regex(
        str = GT,
        pattern = "\\s+",
        replacement = " ",
        vectorize_all = FALSE
      ),
      # trim unnecessary whitespaces at start and end of string
      GT = stringi::stri_trim_both(
        str = GT,
        pattern = "\\P{Wspace}"
      )
    )

  # create a data frame --------------------------------------------------------
  # separate the dataset by space
  genotype.rows <- stringi::stri_split_fixed(str = data$GT, pattern = " ")
  genotype.counts <- lengths(genotype.rows)
  if (any(genotype.counts != length(markers))) {
    bad.rows <- which(genotype.counts != length(markers))
    rlang::abort(paste0(
      "Genepop genotype count does not match the number of loci for individual row(s): ",
      paste(utils::head(bad.rows, 10L), collapse = ", "), "."
    ))
  }
  data <- tibble::as_tibble(
    do.call(rbind, genotype.rows),
    .name_repair = "minimal"
  ) %>%
    magrittr::set_colnames(x = ., markers)
  markers <- NULL

  # Population info ------------------------------------------------------------
  # Strata
  if (!is.null(strata)) {
    #join strata and data
    individuals %<>%
      dplyr::left_join(
        genometranslator::read_strata(strata = strata, verbose = FALSE) %$% strata,
        by = "INDIVIDUALS"
      )
    if (anyNA(individuals$STRATA)) {
      rlang::abort("Some Genepop individuals are absent from the supplied strata metadata.")
    }

  } else {
    # add pop based on internal genepop: integer and reorder the columns
    individuals %<>% dplyr::mutate(STRATA = pop)
  }

  # combine the individuals back to the dataset
  data <- dplyr::bind_cols(individuals, data)
  individuals <- NULL

  # Scan for genotype coding and tidy ------------------------------------------
  gt.coding <- dplyr::select(.data = data, -INDIVIDUALS, -STRATA) %>%
    purrr::flatten_chr(.) %>%
    unique(.) %>%
    nchar(.) %>%
    unique(.)

  if (length(gt.coding) != 1) {
    rlang::abort("Mixed genotype codings are not supported:
  use 1, 2 or 3 characters/numbers for alleles")
  } else {
    if (!gt.coding %in% c(2L, 4L, 6L)) {
      rlang::abort("Genepop genotypes must use one-, two-, or three-digit alleles.")
    }
    gt.sep <- gt.coding / 2L
    data <- tgbase::trans_long(
      x = data,
      cols = c("STRATA", "INDIVIDUALS"),
      names_to = "MARKERS",
      values_to = "GT"
    ) %>%
      tidyr::separate(
        data = ., col = GT, into = c("A1", "A2"),
        sep = gt.sep, remove = TRUE, extra = "drop"
      ) %>%
      dplyr::mutate(
        A1 = stringi::stri_pad_left(str = A1, pad = "0", width = 3),
        A2 = stringi::stri_pad_left(str = A2, pad = "0", width = 3),
        GT = stringi::stri_join(A1, A2),
        A1 = NULL, A2 = NULL
      ) %>%
      dplyr::arrange(MARKERS, STRATA, INDIVIDUALS)

    if (!tidy) {
      data <- tgbase::trans_wide(
        x = data,
        formula = "STRATA + INDIVIDUALS ~ MARKERS",
        values_from = "GT"
      ) %>%
        dplyr::arrange(STRATA, INDIVIDUALS)
    }
  }

  # writing to a file  ---------------------------------------------------------
  if (!is.null(filename)) readr::write_tsv(x = data, file = filename, col_names = TRUE)

  if (!is.null(id.conversion)) attr(data, "id_conversion") <- id.conversion
  return(data)
} # end read_genepop


# write_genepop-----------------------------------------------------------------
#' @name write_genepop

#' @title Write a Genepop file

#' @description Write a genepop file from a tidy data frame or GDS file/object.
#' Used internally in \href{https://github.com/thierrygosselin/genometranslator}{genometranslator}
#' and \href{https://github.com/thierrygosselin/assigner}{assigner}
#' and might be of interest for users.

#' @inheritParams genometranslator_common_arguments

#' @param pop.levels (optional, string) A character string with your populations ordered.
#' Default: \code{pop.levels = NULL}. Described in \code{\link{read_strata}}.

#' @param genepop.header The first line of the Genepop file.
#' With the default, Will use "radiator genepop with date".

#' Default: \code{genepop.header = NULL}.
#' @param markers.line (optional, logical) In the genepop and structure
#' file, you can write the markers on a single line separated by
#' commas \code{markers.line = TRUE},
#' or have markers on a separate line, i.e. in one column, for the genepop file
#' (not very useful with thousands of markers) and not printed at all for the
#' structure file.
#' Default: \code{markers.line = TRUE}.

#' @param filename (optional) The file name prefix for the genepop file
#' written to the working directory. With default: \code{filename = NULL},
#' the date and time is appended to \code{radiator_genepop_}.

#' Default: \code{filename = NULL}.
#' @param ... other parameters passed to the function.

#' @return A genepop file is saved to the working directory.

#' @export
#' @rdname write_genepop
#' @references Raymond M. & Rousset F, (1995).
#' GENEPOP (version 1.2): population genetics software for exact tests
#' and ecumenicism.
#' J. Heredity, 86:248-249
#' @references Rousset F.
#' genepop'007: a complete re-implementation of the genepop software
#' for Windows and Linux.
#' Molecular Ecology Resources.
#' 2008, 8: 103-106.
#' doi:10.1111/j.1471-8286.2007.01931.x

#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @template writer-filtering
#' @template io-dependencies


write_genepop <- function(
    data,
    pop.levels = NULL,
    genepop.header = NULL,
    markers.line = TRUE,
    filename = NULL,
    ...
) {


  # # For testing
  # pop.levels = NULL
  # genepop.header = NULL
  # markers.line = TRUE
  # filename = NULL

  options(stringsAsFactors = FALSE)
  cli::cli_progress_step("Reading data")
  # Checking for missing and/or default arguments ------------------------------
  if (missing(data)) rlang::abort("Input file missing")

  # File type detection----------------------------------------------------------
  data.type <- genometranslator::detect_genomic_format(data)

  # Import data ---------------------------------------------------------------

  if (data.type %in% c("SeqVarGDSClass", "gds.file")) {
    if (data.type == "gds.file") data %<>% genometranslator::read_genome(data = .)
    data <- tidy_genome(data = data, pop.id = FALSE, parallel.core = parallel::detectCores() - 1)
    data.type <- "tbl_df"
  } else {
    if (is.vector(data)) data %<>% genometranslator::read_genome(data = ., )
  }

  if (!rlang::has_name(data, "GT")) {
    cli::cli_progress_step("Recoding genotypes...")
    data <- gt_recoding(x = data, gt = TRUE, alt.dosage = FALSE, gt.vcf = FALSE, gt.vcf.nuc = FALSE)
  }

  cli::cli_progress_step("Preparing data")

  data %<>% dplyr::select(STRATA, INDIVIDUALS, MARKERS, GT)

  # pop.levels -----------------------------------------------------------------
  if (!is.null(pop.levels)) {
    data %<>%
      dplyr::mutate(
        STRATA = factor(STRATA, levels = pop.levels, ordered = TRUE),
        STRATA = droplevels(STRATA)
      ) %>%
      dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  } else {
    data %<>%
      dplyr::mutate(STRATA = factor(STRATA)) %>%
      dplyr::arrange(STRATA, INDIVIDUALS, MARKERS)
  }

  # Create a marker vector  ------------------------------------------------
  markers <- dplyr::distinct(data, MARKERS) %>% dplyr::arrange(MARKERS) %$% MARKERS

  # Wide format ----------------------------------------------------------------
  data  %<>%
    dplyr::arrange(MARKERS) %>%
    tgbase::trans_wide(x = ., formula = "STRATA + INDIVIDUALS ~ MARKERS", values_from = "GT") %>%
    dplyr::arrange(STRATA, INDIVIDUALS) %>%
    dplyr::mutate(INDIVIDUALS = paste(INDIVIDUALS, ",", sep = ""))

  # Write the file in genepop format -------------------------------------------
  # Date and time
  file.date <- format(Sys.time(), "%Y%m%d@%H%M")

  # Filename -------------------------------------------------------------------
  if (is.null(filename)) {
    filename <- stringi::stri_join("radiator_genepop_", file.date, ".gen")
  } else {
    filename.problem <- file.exists(filename)
    if (filename.problem) {
      filename <- stringi::stri_join(filename, "_genepop_", file.date, ".gen")
    } else {
      filename <- stringi::stri_join(filename, "_genepop", ".gen")
    }
  }

  # genepop header  ------------------------------------------------------------
  if (is.null(genepop.header)) {
    genepop.header <- stringi::stri_join("radiator genepop ", file.date)
  }

  # genepop construction
  # could probably use purrr here...
  cli::cli_progress_step("Writing genepop")

  filename.connection <- file(filename, "w") # open the connection to the file
  writeLines(text = genepop.header, con = filename.connection, sep = "\n") # write the genepop header
  if (markers.line) { # write the markers on a single line
    writeLines(text = stringi::stri_join(markers, sep = ",", collapse = ", "), con = filename.connection, sep = "\n")
  } else {# write the markers on a single column (separate lines)
    writeLines(text = stringi::stri_join(markers, sep = "\n"), con = filename.connection, sep = "\n")
  }
  close(filename.connection) # close the connection

  write_gen <- function(x, filename) {
    readr::write_delim(x = as.data.frame("pop"), file = filename, delim = "\n", append = TRUE, col_names = FALSE)
    readr::write_delim(x = x, file = filename, delim = " ", append = TRUE, col_names = FALSE)
  }

  purrr::walk(
    .x = dplyr::group_split(.tbl = data, STRATA, .keep = FALSE),
    .f = write_gen,
    filename = filename
  )
  cli::cli_progress_step(stringi::stri_join("Genepop file: ", filename))

  invisible(filename)
}# End write_genepop
