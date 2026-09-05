test_that("pcadapt import preserves genotype orientation and missingness", {
  f <- export_fixture(incomplete = TRUE)
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  SeqArray::seqSetFilter(f$gds, variant.id = 1:2, verbose = FALSE)
  before <- SeqArray::seqGetFilter(f$gds)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  x <- write_pcadapt(f$gds, strata = f$strata, filename = "pc",
    path.folder = folder, verbose = FALSE)
  expect_equal(x$pop.string, f$strata$STRATA)
  expect_equal(x$samples$INDIVIDUALS, paste0("s", 1:4))
  expect_equal(x$loci$VARIANT_ID, 1:2)
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  skip_if_not_installed("pcadapt")
  expect_warning(bed <- pcadapt::read.pcadapt(x$files[1], type = "pcadapt",
    type.out = "bed"), "more individuals than SNPs")
  expect_true(file.exists(bed))
  expect_equal(attr(bed, "n"), 4)
  expect_equal(attr(bed, "p"), 2)
  observed <- t(pcadapt::bed2matrix(bed))
  expected <- x$genotype.matrix
  expected[expected == 9L] <- NA_integer_
  expect_equal(unname(observed), unname(expected))
})

test_that("pcadapt checks polymorphism after population selection", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  SeqArray::seqSetFilter(f$gds, variant.id = 1:2, verbose = FALSE)
  before <- SeqArray::seqGetFilter(f$gds)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  expect_error(write_pcadapt(f$gds, strata = f$strata,
    pop.select = "Population A", path.folder = folder, verbose = FALSE),
    "polymorphic")
  expect_length(list.files(folder), 0)
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
})
