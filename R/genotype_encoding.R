# separate_gt ------------------------------------------------------------------
#' @title separate_gt
#' @description Separate genotype field
#' @rdname separate_gt
#' @keywords internal
#' @param parallel.core 
#' Default: \code{parallel.core = parallel::detectCores() - 1}.
#' 
#' @param split.chunks 
#' Default: \code{split.chunks = 3}.
#' 
#' @param filter.missing 
#' Default: \code{filter.missing = FALSE}.
#' 
#' @param remove 
#' Default: \code{remove = TRUE}.
#' 
#' @param alleles.naming 
#' Default: \code{alleles.naming = c("A1", "A2")}.
#' 
#' @param exclude 
#' Default: \code{exclude = c("LOCUS", "INDIVIDUALS", "POP_ID")}.
#' 
#' @param haplotypes 
#' Default: \code{haplotypes = FALSE}.
#' 
#' @param gather 
#' Default: \code{gather = TRUE}.
#' 
#' @param gt 
#' Default: \code{gt = "GT_VCF_NUC"}.
#' 
#' @export
separate_gt <- function(
    x,
    gt = "GT_VCF_NUC",
    gather = TRUE,
    haplotypes = FALSE,
    exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"),
    alleles.naming = c("A1", "A2"),
    remove = TRUE,
    filter.missing = FALSE,
    split.chunks = 3,
    parallel.core = parallel::detectCores() - 1
) {

  ## TEST
  # gather = TRUE
  # haplotypes = FALSE
  # exclude = c("LOCUS", "INDIVIDUALS", "POP_ID")
  # alleles.naming = c("A1", "A2")
  # remove = TRUE
  # filter.missing = FALSE
  # split.chunks = 3


  separate_genotype <-  carrier::crate(function(x, gt, alleles.naming, remove, filter.missing, gather, haplotypes, exclude){
    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`

    # remove SPLIT_VEC from exclude if detected
    if (rlang::has_name(x, "SPLIT_VEC")) {
      exclude <- c(exclude, "SPLIT_VEC")
    }



    # discard the other gt format
    gt.format <- c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC")
    not.wanted <- setdiff(gt.format, gt)
    x %<>% dplyr::select(-tidyselect::any_of(not.wanted))

    if (gt == "GT_VCF_NUC") {
      if (filter.missing) x  %<>% dplyr::filter(GT_VCF_NUC != "./.")
      x %<>%
        dplyr::bind_cols(
          stringi::stri_split_fixed(str = x$GT_VCF_NUC, pattern = "/", simplify = TRUE) %>%
            magrittr::set_colnames(x = ., value = alleles.naming) %>%
            tibble::as_tibble()
        )
    }
    if (gt == "GT_VCF") {
      if (filter.missing) x  %<>% dplyr::filter(GT_VCF != "./.")
      x %<>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT_VCF, from = 1, to = 1),
          A2 = stringi::stri_sub(str = GT_VCF, from = 3, to = 3)
        )
    }
    if (gt == "ALT_DOSAGE") {
      if (filter.missing) x  %<>% dplyr::filter(!is.na(ALT_DOSAGE))
      x %<>%
        dplyr::mutate(
          A1 = dplyr::if_else(ALT_DOSAGE == 0L, 1L, ALT_DOSAGE),
          A2 = dplyr::if_else(ALT_DOSAGE != 2L, ALT_DOSAGE + 1L, ALT_DOSAGE)
        )
    }
    if (gt == "GT") {
      if (filter.missing) x  %<>% dplyr::filter(GT != "000000")
      x %<>%
        dplyr::mutate(
          A1 = stringi::stri_sub(str = GT, from = 1, to = 3),
          A2 = stringi::stri_sub(str = GT, from = 4, to = 6)
        )
    }

    if (remove) x %<>% dplyr::select(-tidyselect::any_of(gt.format))

    if (gather) {
      x %<>%
        tgbase::trans_long(
          x = .,
          cols = exclude,
          names_to = "ALLELES_GROUP",
          values_to = "ALLELES",
          variable_factor = FALSE
        )
      if (haplotypes) x %<>% dplyr::rename(HAPLOTYPES = ALLELES)
    }

    return(x)
  })#End separate_genotype


  if (split.chunks > 1) {
    x %<>%
      tgbase::parallel_map(
        .x = .,
        .f = separate_genotype,
        flat.future = "dfr",
        split.vec = TRUE,
        split.with = NULL,
        split.chunks = split.chunks,
        parallel.core = split.chunks,
        forking = TRUE,
        gt = gt,
        alleles.naming = alleles.naming,
        remove = remove,
        filter.missing = filter.missing,
        gather = gather,
        haplotypes = haplotypes,
        exclude = exclude
      )
  } else {
    x %<>%
      separate_genotype(
        x = .,
        gt = gt,
        alleles.naming = alleles.naming,
        remove = remove,
        filter.missing = filter.missing,
        gather = gather,
        haplotypes = haplotypes,
        exclude = exclude)
  }
  return(x)
}#End separate_gt

# detect_gt ---------------------------------------------------------------------
#' @rdname detect_gt
#' @title detect_gt
#' @description Detect the genotype format used in the data set.
#' @param x The data
#' @param gt.format (character)
#' Default: \code{gt.format = c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC")}.
#' @param keep.one (logical) Will return only one format if \code{keep.one = TRUE}.
#' Default: \code{keep.one = TRUE}.
#' @param favorite If more than one format is present and \code{keep.one = TRUE},
#' the favorite will be returned, if present. Otherwise, the first format in
#' `gt.format` is returned.
#' Default: \code{favorite = "ALT_DOSAGE"}.
#' @keywords internal
#' @export
detect_gt <- function(x, gt.format = c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC"), keep.one = TRUE, favorite = "ALT_DOSAGE") {

  detect.gt <- intersect(gt.format, colnames(x))

  if (length(detect.gt) == 0L) detect.gt <- NULL

  if (keep.one) {
    if (length(detect.gt) > 1) {
      if (favorite %in% detect.gt) {
        detect.gt <- favorite
      } else {
        detect.gt <- detect.gt[[1L]]
      }
    }
  }
  return(detect.gt)
}#End detect_gt



# gt_recoding ---------------------------------------------------------------------
#' @rdname gt_recoding
#' @title gt_recoding
#' @description Translate among supported genotype encodings.
#' @keywords internal
#' @param arrange 
#' Default: \code{arrange = TRUE}.
#' 
#' @param gt.vcf.nuc 
#' Default: \code{gt.vcf.nuc = TRUE}.
#' 
#' @param gt.vcf 
#' Default: \code{gt.vcf = TRUE}.
#' 
#' @param alt.dosage 
#' Default: \code{alt.dosage = TRUE}.
#' 
#' @param gt 
#' Default: \code{gt = TRUE}.
#' 
#' @export
gt_recoding <- function(x, gt = TRUE, alt.dosage = TRUE, gt.vcf = TRUE, gt.vcf.nuc = TRUE, arrange = TRUE) {
  # what genotype format we have
  detect.gt <- detect_gt(x) #utils

  # if it's just GT no need to calibrate or generate other genotypes formats
  if (all(length(detect.gt) == 1, detect.gt == "GT", !any(alt.dosage, gt.vcf, gt.vcf.nuc))) return(x)

  # conditions and checks
  if (arrange) x %<>% dplyr::mutate(IDTEMP = seq_len(dplyr::n()))

  # if (gt.vcf.nuc && !rlang::has_name(x, "REF")) gt.vcf.nuc <- FALSE
  # e.g. if we have GT and we want ALT_DOSAGE or any other format that requires nucleotide info
  # we could issue a message...

  if (all(detect.gt == "GT") && !rlang::has_name(x, "REF") && any(alt.dosage, gt.vcf, gt.vcf.nuc)) {
    message("\n\nGenerating automatically 3 new genotypes formats, see doc...")
    message("Format(s) chosen requires nucleotide information not found in the data")
    message("001 converted to A")
    message("002 converted to C")
    message("003 converted to G")
    message("004 converted to T")

    # generate all the other format from GT
    if (all(!rlang::has_name(x, "A1"), rlang::has_name(x, "REF"))) {
      x  %<>%
        dplyr::mutate(
          A1 = dplyr::recode(REF, "A" = "001", "C" = "002", "G" = "003", "T" = "004"),
          A2 = dplyr::recode(ALT, "A" = "001", "C" = "002", "G" = "003", "T" = "004")
        )
      remove.extra <- TRUE
    } else {# if no REF
      x1 <- dplyr::select(x, -tidyselect::any_of(c("POP_ID", "IDTEMP")))
      gt.col <- purrr::discard(.x = colnames(x1), .p = colnames(x1) %in% c("GT", "GT_VCF", "ALT_DOSAGE", "GT_VCF_NUC"))

      # need to do it in 2 steps to get the het genotypes
      x1 %<>%
        separate_gt(x = ., gt = "GT", gather = FALSE, exclude = gt.col) %>%
        dplyr::mutate(
          HET = dplyr::case_when(
            A1 == "000" ~ NA_integer_,
            A1 == A2 ~ 0L,
            A1 != A2 ~ 1L
          ),
          SPLIT_VEC = NULL
        ) %>%
        # gather the results at this stage
        tgbase::trans_long(
          x = .,
          cols = c(gt.col, "HET"),
          names_to = "ALLELES_GROUP",
          values_to = "ALLELES",
          variable_factor = FALSE
        )


      # checking what type of GT format... inspired by the codes in gtypes.R
      gt.wanted <- sort(unique(x1$ALLELES))
      gt.wanted <- purrr::keep(.x = gt.wanted, .p = gt.wanted %in% c("001", "002", "003", "004"))
      if (length(gt.wanted) != 4L) rlang::abort("Contact author problem with genotype format used")

      missing <- x1 %>%
        dplyr::filter(ALLELES == "000") %>%
        dplyr::distinct(MARKERS, ALLELES) %>%
        dplyr::mutate(DOS_ALT = NA_integer_)

      x2 <- x1 %>%
        dplyr::filter(ALLELES != "000") %>%
        dplyr::count(MARKERS, ALLELES) %>%
        dplyr::group_by(MARKERS) %>%
        dplyr::mutate(ALLELES_GROUP = dplyr::if_else(n == max(n, na.rm = TRUE), "REF", "ALT", "EQUAL")) %>%
        dplyr::group_by(MARKERS, ALLELES_GROUP) %>%
        dplyr::mutate(EQUAL_COUNTS = dplyr::n()) %>%
        dplyr::ungroup(.)


      equal.ref <- x2 %>%
        dplyr::filter(EQUAL_COUNTS == 2L)

      x2 %<>%
        dplyr::filter(EQUAL_COUNTS == 1L) %>%
        dplyr::select(MARKERS, ALLELES, ALLELES_GROUP)

      if (nrow(equal.ref) > 0) {
        equal.ref %<>%
          dplyr::group_by(MARKERS) %>%
          dplyr::mutate(ALLELES_GROUP = c("REF", "ALT")) %>%
          dplyr::select(MARKERS, ALLELES, ALLELES_GROUP)

        x2 %<>% dplyr::bind_rows(equal.ref)
        equal.ref <- NULL
      }

      x2 %<>%
        dplyr::mutate(
          NUC = dplyr::case_when(
            ALLELES == "001" ~ "A",
            ALLELES == "002" ~ "C",
            ALLELES == "003" ~ "G",
            ALLELES == "004" ~ "T"),
          GROUP = dplyr::case_when(
            ALLELES_GROUP == "REF" ~ "A1",
            ALLELES_GROUP == "ALT" ~ "A2"),
          DOS_ALT = dplyr::case_when(
            ALLELES_GROUP == "REF" ~ 0L,
            ALLELES_GROUP == "ALT" ~ 1L)
        )

      missing %<>%
        dplyr::bind_rows(
          x2 %>%
            dplyr::select(MARKERS, ALLELES, DOS_ALT)
        )

      x2 %<>%
        dplyr::select(-DOS_ALT) %>%
        tgbase::trans_wide(x = ., formula = MARKERS ~ ALLELES_GROUP + GROUP, values_from = c("NUC", "ALLELES")) %>%
        dplyr::rename(REF = NUC_REF_A1, ALT = NUC_ALT_A2, A1 = ALLELES_REF_A1, A2 = ALLELES_ALT_A2)

      x1 %<>%
        dplyr::select(-ALLELES_GROUP) %>%
        dplyr::left_join(x2, by = "MARKERS") %>%
        dplyr::left_join(missing, by = c("MARKERS", "ALLELES")) %>%
        dplyr::select(-ALLELES, -A1, -A2)

      x %<>%
        dplyr::left_join(dplyr::distinct(x2, MARKERS, REF, ALT), by = "MARKERS")

      x2 <- missing <- NULL

      het <- x1 %>%
        dplyr::filter(HET == 1L) %>%
        dplyr::distinct(MARKERS, INDIVIDUALS, ALT, REF) %>%
        dplyr::mutate(
          GT_VCF = "0/1",
          ALT_DOSAGE = 1L,
          GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/"),
          ALT = NULL, REF = NULL
        )

      missing <- x1 %>%
        dplyr::filter(is.na(HET)) %>%
        dplyr::distinct(MARKERS, INDIVIDUALS) %>%
        dplyr::mutate(
          GT_VCF = "./.",
          ALT_DOSAGE = NA_integer_,
          GT_VCF_NUC = "./."
        )

      hom <- x1 %>%
        dplyr::filter(HET != 1L & !is.na(HET)) %>%
        dplyr::select(INDIVIDUALS, MARKERS, ALT, REF, DOS_ALT) %>%
        dplyr::group_by(MARKERS, INDIVIDUALS, ALT, REF) %>%
        dplyr::summarise(
          ALT_DOSAGE = sum(DOS_ALT, na.rm = TRUE), .groups = "drop"
        )

      hom <- dplyr::bind_rows(
        hom %>%
          dplyr::filter(ALT_DOSAGE == 0L) %>%
          dplyr::mutate(
            GT_VCF = "0/0",
            GT_VCF_NUC = stringi::stri_join(REF, REF, sep = "/"),
            ALT = NULL, REF = NULL
          ),
        hom %>%
          dplyr::filter(ALT_DOSAGE == 2L) %>%
          dplyr::mutate(
            GT_VCF = "1/1",
            GT_VCF_NUC = stringi::stri_join(ALT, ALT, sep = "/"),
            ALT = NULL, REF = NULL
          )
      )
      x1 <- NULL
      hom %<>% dplyr::bind_rows(het, missing)
      het <- missing <- NULL

      x %<>%
        dplyr::left_join(hom, by = c("MARKERS", "INDIVIDUALS"))
      hom <- NULL
      gt <- FALSE
      remove.extra <- FALSE

      # complete tidy dataset with calibrated alleles...
    }
  } else {
    remove.extra <- FALSE

    if (gt) {
      if (all(!rlang::has_name(x, "A1"), rlang::has_name(x, "REF"))) {
        x  %<>%
          dplyr::mutate(
            A1 = dplyr::recode(REF, "A" = "001", "C" = "002", "G" = "003", "T" = "004"),
            A2 = dplyr::recode(ALT, "A" = "001", "C" = "002", "G" = "003", "T" = "004")
          )
      }

      if (all(!rlang::has_name(x, "A1"), !rlang::has_name(x, "REF"))) {
        x  %<>% dplyr::mutate(A1 =  "001", A2 = "002")
      }
      remove.extra <- TRUE
    }


    gt_map <- function(
    x,
    gt.format = c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC"),
    gt = TRUE,
    alt.dosage = TRUE,
    gt.vcf = TRUE,
    gt.vcf.nuc = TRUE
    ) {

      gt.format <- match.arg(
        arg = gt.format,
        choices = c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC"),
        several.ok = FALSE
      )

      #start with missing genotypes
      if (gt.format == "ALT_DOSAGE") {
        dosage.value <- unique(x$ALT_DOSAGE)
        if (is.na(dosage.value)) {
          x %<>%
            {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "./.") else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = "./.") else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = "000000") else .}
        } else {
          if (dosage.value == 0L) {
            x %<>%
              {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "0/0") else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, REF, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A1)) else .}
          }
          if (dosage.value == 1L) {
            x %<>%
              {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "0/1") else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A2)) else .}
          }
          if (dosage.value == 2L) {
            x %<>%
              {if (gt.vcf) dplyr::mutate(.data = ., GT_VCF = "1/1") else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(ALT, ALT, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A2, A2)) else .}
          }
        }
      }#End ALT_DOSAGE

      if (gt.format == "GT_VCF") {
        gt.vcf <- unique(x$GT_VCF)
        if (gt.vcf == "./.") {
          x %<>%
            {if (alt.dosage) dplyr::mutate(.data = ., ALT_DOSAGE = NA_integer_) else .} %>%
            {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = "./.") else .} %>%
            {if (gt) dplyr::mutate(.data = ., GT = "000000") else .}
        } else {
          if (gt.vcf == "0/0") {
            x %<>%
              {if (alt.dosage) dplyr::mutate(.data = ., ALT_DOSAGE = 0L) else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, REF, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A1)) else .}
          }
          if (gt.vcf == "1/0") {
            x %<>%
              {if (alt.dosage) dplyr::mutate(.data = ., ALT_DOSAGE = 1L) else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A2)) else .}
          }
          if (gt.vcf == "0/1") {
            x %<>%
              {if (alt.dosage) dplyr::mutate(.data = ., ALT_DOSAGE = 1L) else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(REF, ALT, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A1, A2)) else .}
          }
          if (gt.vcf == "1/1") {
            x %<>%
              {if (alt.dosage) dplyr::mutate(.data = ., ALT_DOSAGE = 2L) else .} %>%
              {if (gt.vcf.nuc) dplyr::mutate(.data = ., GT_VCF_NUC = stringi::stri_join(ALT, ALT, sep = "/")) else .} %>%
              {if (gt) dplyr::mutate(.data = ., GT = stringi::stri_join(A2, A2)) else .}
          }
        }
      }#End GT_VCF

      if (gt.format == "GT_VCF_NUC") {
        message("Not implemented, yet...")
      }
      return(x)
    }#End gt_map


    # split the data
    x %<>%
      dplyr::group_split(.data[[detect.gt]]) %>%
      purrr::map_dfr(.x = ., .f = gt_map, gt.format = detect.gt, gt = gt, alt.dosage = alt.dosage, gt.vcf = gt.vcf, gt.vcf.nuc = gt.vcf.nuc)

  }

  if (remove.extra) x %<>% dplyr::select(-c(A1, A2))
  if (arrange) x %<>% dplyr::arrange(IDTEMP) %>% dplyr::select(-IDTEMP)
  return(x)
}#End gt_recoding
