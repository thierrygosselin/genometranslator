
test_that("Arlequin exports correct diploid rows and restores selections", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  before <- SeqArray::seqGetFilter(f$gds)
  x <- write_arlequin(f$gds, strata = f$strata, filename = "test",
    path.folder = folder, chunk.size = 1, verbose = FALSE)
  lines <- readLines(x$files[[1]])
  expect_true("S1 1 1 1 2" %in% lines)
  expect_true(" 1 1 3" %in% lines)
  expect_true("S2 1 1 1 ?" %in% lines)
  expect_true("MissingData='?'" %in% lines)
  expect_true('SampleName="P1"' %in% lines)
  expect_true(all(file.exists(x$files)))
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_arlequin(f$gds, strata = f$strata,
    pop.levels = "Population A", path.folder = folder), "every observed")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_arlequin(f$gds, strata = f$strata, filename = "test",
    path.folder = folder), "Output exists")
  SeqArray::seqSetFilter(f$gds, sample.id = c("s1", "s3"),
    variant.id = 2L, verbose = FALSE)
  selected <- SeqArray::seqGetFilter(f$gds)
  y <- write_arlequin(f$gds, strata = f$strata, filename = "subset",
    path.folder = folder, verbose = FALSE)
  expect_equal(y$variant.id, 2L)
  expect_equal(nrow(y$samples), 2)
  expect_identical(SeqArray::seqGetFilter(f$gds), selected)
})

test_that("fastsimcoal folding uses global counts and splits ties", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  before <- SeqArray::seqGetFilter(f$gds)
  x <- write_fastsimcoal(f$gds, strata = f$strata, filename = "test",
    path.folder = folder, invariant.sites = 10, chunk.size = 1, verbose = FALSE)
  expect_equal(unname(x$dimensions), c(5, 5))
  expect_equal(sum(x$spectrum), 12)
  expect_equal(x$spectrum[c(1, 2, 9, 17)], c(10, 1, 0.5, 0.5))
  expect_equal(x$audit$STATUS, c("KEPT", "KEPT", "NOT_BIALLELIC_SNP"))
  lines <- readLines(x$files[[1]])
  expect_equal(lines[[2]], "2 4 4")
  values <- scan(text = lines[[3]], quiet = TRUE)
  expect_equal(values, x$spectrum)
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  y <- write_fastsimcoal(f$gds, strata = f$strata, filename = "reverse",
    pop.levels = rev(unique(f$strata$STRATA)), path.folder = folder,
    chunk.size = 100, verbose = FALSE)
  expect_equal(y$spectrum[[6]], 1)
  expect_equal(y$populations$GROUP, rev(unique(f$strata$STRATA)))
  expect_error(write_fastsimcoal(f$gds, strata = f$strata,
    folded = FALSE), "ANCESTRAL")
  expect_error(write_fastsimcoal(f$gds, strata = f$strata,
    max.cells = 10), "max.cells")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
})

test_that("unfolded SFS uses explicit ancestral states", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  x <- write_fastsimcoal(f$gds, strata = f$strata, filename = "derived",
    path.folder = folder, folded = FALSE,
    ancestral.alleles = data.frame(VARIANT_ID = 1:2, ANCESTRAL = c("G", "C")),
    verbose = FALSE)
  expect_equal(x$spectrum[c(2, 17)], c(1, 1))
  expect_equal(sum(x$spectrum), 2)
  expect_match(x$files[[1]], "_DSFS.obs$")
  y <- write_fastsimcoal(f$gds, strata = f$strata, filename = "unknown",
    path.folder = folder, folded = FALSE,
    ancestral.alleles = data.frame(VARIANT_ID = 1:2, ANCESTRAL = c("N", "C")),
    verbose = FALSE)
  expect_equal(y$audit$STATUS[[1]], "UNKNOWN_ANCESTRAL")
  expect_equal(sum(y$spectrum), 1)
})

test_that("export messages have a hard 80-column limit", {
  text <- capture.output(.export_message(paste(rep("x", 205), collapse = "")),
    type = "message")
  expect_true(all(nchar(text) <= 80))
})

test_that("SFS audits incomplete records and filepath connections close", {
  f <- export_fixture(incomplete = TRUE)
  SeqArray::seqClose(f$gds)
  on.exit(unlink(f$path), add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  x <- write_fastsimcoal(f$path, strata = f$strata,
    path.folder = folder, verbose = FALSE)
  expect_equal(x$audit$STATUS[[2]], "INCOMPLETE_CALLS")
  expect_equal(sum(x$spectrum), 1)
  g <- SeqArray::seqOpen(f$path)
  expect_equal(length(SeqArray::seqGetData(g, "variant.id")), 3)
  SeqArray::seqClose(g)
})

test_that("metadata blacklists are respected without mutating the selection", {
  f <- export_fixture()
  SeqArray::seqClose(f$gds)
  g <- SeqArray::seqOpen(f$path, readonly = FALSE)
  on.exit({SeqArray::seqClose(g); unlink(f$path)}, add = TRUE)
  mm <- extract_markers_metadata(g)
  mm$FILTERS[[3]] <- "test_blacklist"
  update_genome_gds(g, node.name = "markers.meta", value = mm, sync = FALSE)
  before <- SeqArray::seqGetFilter(g)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  x <- write_arlequin(g, strata = f$strata, path.folder = folder,
    verbose = FALSE)
  expect_equal(x$variant.id, 1:2)
  expect_identical(SeqArray::seqGetFilter(g), before)
})
