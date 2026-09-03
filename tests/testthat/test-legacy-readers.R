testthat::test_that("read_fstat validates and standardizes FSTAT data", {
  path <- tempfile(fileext = ".dat")
  writeLines(c(
    "2 2 2 2",
    "L1",
    "L2",
    "1 0101 0102",
    "2 0202 0000"
  ), path)

  long <- genometranslator::read_fstat(path)
  wide <- genometranslator::read_fstat(path, tidy = FALSE)

  testthat::expect_named(long, c("STRATA", "INDIVIDUALS", "MARKERS", "GENOTYPE"))
  testthat::expect_equal(nrow(long), 4L)
  testthat::expect_true(all(nchar(long$GENOTYPE) == 6L))
  testthat::expect_equal(nrow(wide), 2L)
  testthat::expect_true(all(c("STRATA", "INDIVIDUALS", "L1", "L2") %in% names(wide)))
})

testthat::test_that("read_fstat rejects inconsistent rows", {
  path <- tempfile(fileext = ".dat")
  writeLines(c("1 2 2 2", "L1", "L2", "1 0101"), path)
  testthat::expect_error(
    genometranslator::read_fstat(path),
    "do not contain one population field"
  )
})

testthat::test_that("read_genepop handles long and wide output", {
  path <- tempfile(fileext = ".gen")
  writeLines(c(
    "Example Genepop data",
    "L1, L2",
    "Pop",
    "sample-1, 0101 0102",
    "Pop",
    "sample-2, 0202 0000"
  ), path)

  long <- genometranslator::read_genepop(path)
  wide <- genometranslator::read_genepop(path, tidy = FALSE)

  testthat::expect_named(long, c("STRATA", "INDIVIDUALS", "MARKERS", "GT"))
  testthat::expect_equal(nrow(long), 4L)
  testthat::expect_true(all(nchar(long$GT) == 6L))
  testthat::expect_equal(nrow(wide), 2L)
  testthat::expect_true(all(c("STRATA", "INDIVIDUALS", "L1", "L2") %in% names(wide)))
})

testthat::test_that("read_genepop validates separators and genotype counts", {
  no.pop <- tempfile(fileext = ".gen")
  writeLines(c("Example", "L1", "sample-1, 0101"), no.pop)
  testthat::expect_error(genometranslator::read_genepop(no.pop), "separator")

  bad.count <- tempfile(fileext = ".gen")
  writeLines(c("Example", "L1,L2", "Pop", "sample-1, 0101"), bad.count)
  testthat::expect_error(genometranslator::read_genepop(bad.count), "count does not match")
})

testthat::test_that("read_genind modern return paths work", {
  testthat::skip_if_not_installed("adegenet")
  x <- adegenet::df2genind(
    data.frame(L1 = c("A/A", "A/B"), L2 = c("C/C", "C/T")),
    sep = "/",
    ploidy = 2,
    pop = factor(c("one", "two"))
  )

  tidy <- genometranslator::read_genind(x, gds = FALSE)
  original <- genometranslator::read_genind(x, tidy = FALSE, gds = FALSE)

  testthat::expect_s3_class(tidy, "tbl_df")
  testthat::expect_true(all(c("STRATA", "INDIVIDUALS", "MARKERS") %in% names(tidy)))
  testthat::expect_s4_class(original, "genind")
})

testthat::test_that("read_genlight generates explicit de novo metadata", {
  testthat::skip_if_not_installed("adegenet")
  x <- adegenet::as.genlight(matrix(c(0, 1, 2, NA), nrow = 2L))
  adegenet::indNames(x) <- c("sample-1", "sample-2")
  adegenet::pop(x) <- factor(c("one", "two"))

  tidy <- genometranslator::read_genlight(x, gds = FALSE)
  original <- genometranslator::read_genlight(x, tidy = FALSE, gds = FALSE)

  testthat::expect_s3_class(tidy, "tbl_df")
  testthat::expect_true(all(tidy$CHROM == "DENOVO"))
  testthat::expect_s4_class(original, "genlight")
})

testthat::test_that("read_gtypes rejects unsupported input clearly", {
  testthat::expect_error(
    genometranslator::read_gtypes(data.frame()),
    "not a gtypes object"
  )
})

testthat::test_that("read_gtypes standardizes diploid allele rows", {
  if (!methods::isClass("gtypes")) {
    methods::setClass("gtypes", slots = c(data = "data.frame"))
  }
  x <- methods::new(
    "gtypes",
    data = data.frame(
      id = rep(c("sample_1", "sample_2"), each = 2L),
      stratum = rep(c("South East", "North"), each = 2L),
      locus = "locus-1",
      allele = c("A", "C", "G", NA_character_),
      stringsAsFactors = FALSE
    )
  )

  result <- genometranslator::read_gtypes(x)
  by.individual <- dplyr::arrange(result, INDIVIDUALS)

  testthat::expect_named(result, c("STRATA", "INDIVIDUALS", "MARKERS", "GT"))
  testthat::expect_equal(by.individual$INDIVIDUALS, c("sample-1", "sample-2"))
  testthat::expect_equal(by.individual$STRATA, c("South_East", "North"))
  testthat::expect_equal(by.individual$GT, c("001002", "003000"))
})
