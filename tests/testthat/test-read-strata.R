testthat::test_that("read_strata prepares sample metadata", {
  metadata <- data.frame(
    INDIVIDUALS = c("sample_1", "sample 2"),
    POP_ID = c("South East", "South East"),
    SEX = c("F", "M"),
    stringsAsFactors = FALSE
  )

  result <- genometranslator::read_strata(metadata)

  testthat::expect_named(
    result,
    c("strata", "pop.levels", "pop.labels", "pop.select", "blacklist.id")
  )
  testthat::expect_equal(
    as.character(result$strata$INDIVIDUALS),
    c("sample-1", "sample2")
  )
  testthat::expect_equal(as.character(result$strata$STRATA), rep("South_East", 2))
  testthat::expect_equal(result$strata$SEX, c("F", "M"))
  testthat::expect_true(is.factor(result$strata$STRATA))
  testthat::expect_false("POP_ID" %in% names(result$strata))
  testthat::expect_false("pop.id" %in% names(formals(genometranslator::read_strata)))
})

testthat::test_that("read_strata validates its metadata contract", {
  testthat::expect_error(
    genometranslator::read_strata(data.frame(INDIVIDUALS = "sample-1")),
    "missing required column.*STRATA"
  )
  testthat::expect_error(
    genometranslator::read_strata(data.frame(
      INDIVIDUALS = "",
      STRATA = "population"
    )),
    "INDIVIDUALS.*missing or empty"
  )
  testthat::expect_error(
    genometranslator::read_strata(data.frame(
      INDIVIDUALS = "sample-1",
      STRATA = NA_character_
    )),
    "STRATA.*missing or empty"
  )
})

testthat::test_that("read_strata checks population labelling arguments", {
  metadata <- data.frame(
    INDIVIDUALS = c("sample-1", "sample-2"),
    STRATA = c("A", "B")
  )

  testthat::expect_error(
    genometranslator::read_strata(metadata, pop.labels = c("one", "two")),
    "pop.levels is required"
  )
  testthat::expect_error(
    genometranslator::read_strata(
      metadata,
      pop.levels = c("A", "B"),
      pop.labels = "one"
    ),
    "different length"
  )
})
