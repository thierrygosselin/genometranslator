#' Split genotype strings into allele columns
#'
#' Splits a vector of genotype strings into one row per genotype and one column
#' per allele. The implementation is vectorized and does not require parallel
#' processing.
#'
#' @param x Character vector of genotype strings.
#' @param separator Fixed separator between alleles.
#' Default: \code{separator = "/"}.
#' @param alleles.naming Optional column names. When `NULL`, names are generated
#'   as `A1`, `A2`, and so on.
#'
#' Default: \code{alleles.naming = NULL}.
#' @return A tibble with one row per genotype and one allele per column.
#' @export
split_genotypes <- function(x, separator = "/", alleles.naming = NULL) {
  if (!is.character(x)) x <- as.character(x)
  if (length(separator) != 1L || is.na(separator) || !nzchar(separator)) {
    rlang::abort("`separator` must be one non-empty string.")
  }

  alleles <- stringi::stri_split_fixed(
    str = x,
    pattern = separator,
    simplify = TRUE
  )

  n.alleles <- ncol(alleles)
  if (is.null(alleles.naming)) {
    alleles.naming <- paste0("A", seq_len(n.alleles))
  }
  if (length(alleles.naming) != n.alleles || anyDuplicated(alleles.naming)) {
    rlang::abort("`alleles.naming` must contain one unique name per allele column.")
  }

  colnames(alleles) <- alleles.naming
  tibble::as_tibble(alleles)
}
