#' Check genometranslator dependencies
#'
#' Reports required and optional R packages and checks whether the optional
#' \command{bcftools} and \command{plink2} executables are available. This function is intentionally
#' diagnostic: it does not modify the R library or a Conda environment.
#'
#' @param verbose Logical. Print installation guidance for missing components.
#' Default: \code{verbose = TRUE}.
#'
#' @return A tibble with component, source, requirement level, and availability.
#' @export
genometranslator_dependencies <- function(verbose = TRUE) {
  required.cran <- c(
    "arrow", "carrier", "cli", "data.table", "dplyr", "magrittr",
    "matrixStats", "purrr", "readr", "rlang", "stringi", "sys", "tibble",
    "tidyr", "tidyselect", "vctrs", "vroom", "withr"
  )
  required.bioc <- c("gdsfmt", "Rsamtools", "SeqArray")
  required.github <- "tgbase"
  optional.cran <- c("adegenet", "fst", "knitr", "rmarkdown")
  optional.bioc <- "SNPRelate"
  optional.github <- "strataG"

  packages <- tibble::tibble(
    component = c(
      required.cran, required.bioc, required.github,
      optional.cran, optional.bioc, optional.github
    ),
    source = c(
      rep("CRAN", length(required.cran)),
      rep("Bioconductor", length(required.bioc)),
      "GitHub: thierrygosselin/tgbase",
      rep("CRAN", length(optional.cran)),
      "Bioconductor",
      "GitHub: EricArcher/strataG"
    ),
    required = c(
      rep(TRUE, length(required.cran) + length(required.bioc) + 1L),
      rep(FALSE, length(optional.cran) + 2L)
    )
  )
  packages$available <- vapply(
    packages$component,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )

  executables <- tibble::tibble(
    component = c("bcftools", "plink2"),
    source = "Conda, Homebrew, or system PATH",
    required = FALSE,
    available = c(
      nzchar(Sys.which("bcftools")),
      nzchar(Sys.which("plink2"))
    )
  )
  result <- dplyr::bind_rows(packages, executables)

  if (verbose && any(!result$available)) {
    message("Missing genometranslator components:")
    apply(result[!result$available, , drop = FALSE], 1L, function(x) {
      message("  ", x[["component"]], " [", x[["source"]], "]")
    })
    message("See the Installation section in the genometranslator README.")
  }

  result
}
