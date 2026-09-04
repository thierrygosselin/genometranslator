test_that("Genepop reads mixed locus widths and guarded continuations", {
  make <- function(lines) data.frame(lines = lines)
  x <- read_genepop(make(c("Title", "L1,L2", " pOp ",
    "s1, 0100", "002003", "s2, 0002 000000")))
  expect_equal(x$GT, c("001000", "000002", "002003", "000000"))
  expect_equal(attr(x, "genepop_title"), "Title")
  expect_error(read_genepop(make(c("t", "L1", "Pop", "s, 01"))), "diploid")
  expect_error(read_genepop(make(c("t", "L1", "Pop", "s, 0x00"))), "numeric")
  expect_error(read_genepop(make(c("t", "L1", "Pop", "s, 0101", "Pop", "0101"))),
    "continuation")
  expect_error(read_genepop(make(c("t", "L1", "Pop", "Pop", "s, 0101"))), "Empty")
  expect_error(read_genepop(make(c("t", "L1", "Pop", "s, 0101", "b, 001001"))),
    "consistent")
})

test_that("Genepop GDS writer roundtrips multi-alleles and missing copies", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  before <- SeqArray::seqGetFilter(f$gds)
  x <- write_genepop(f$gds, strata = f$strata, filename = "test.gen",
    path.folder = folder, markers.line = FALSE, chunk.size = 1, verbose = FALSE,
    pop.levels = rev(unique(f$strata$STRATA)))
  expect_true(all(file.exists(x$files)))
  expect_equal(basename(x$files[1]), "test.gen")
  y <- read_genepop(x$files[1])
  expect_equal(nrow(y), 12)
  expect_equal(y$GT[y$INDIVIDUALS == "S1" & y$MARKERS == "L3"], "002003")
  expect_equal(y$GT[y$INDIVIDUALS == "S2" & y$MARKERS == "L3"], "000000")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_genepop(f$gds, strata = f$strata, filename = "test",
    path.folder = folder, verbose = FALSE), "Output exists")
  expect_error(write_genepop(f$gds, strata = f$strata, pop.levels = "absent",
    path.folder = folder, verbose = FALSE), "every observed")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
})

test_that("Genepop tidy writer validates completeness and preserves partial calls", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  d <- data.frame(STRATA = "p", INDIVIDUALS = c("a", "b"),
    MARKERS = "m", GT = c("001000", "000002"))
  x <- write_genepop(d, path.folder = folder, filename = "tidy", verbose = FALSE)
  expect_equal(read_genepop(x$files[1])$GT, d$GT)
  expect_error(write_genepop(rbind(d, d[1, ]), path.folder = folder), "exactly one")
})
