sfs_all_fixture <- function() {
  vcf <- tempfile(fileext = ".vcf"); path <- tempfile(fileext = ".gds")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
    "1\t10\ti1\tA\t.\t.\tPASS\t.\tGT\t0/0\t0/0",
    "1\t20\ti2\tC\t.\t.\tPASS\t.\tGT\t0/0\t./.",
    "1\t30\ts1\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/1",
    "1\t40\ts2\tA\tG\t.\tPASS\t.\tGT\t1/1\t1/1"), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  unlink(vcf)
  list(path = path, gds = SeqArray::seqOpen(path, readonly = FALSE),
    strata = data.frame(INDIVIDUALS = c("a", "b"), STRATA = c("A", "B")))
}

test_that("all-sites spectra count invariant records and audit missing calls", {
  f <- sfs_all_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  before <- SeqArray::seqGetFilter(f$gds)
  x <- write_fastsimcoal(f$gds, strata = f$strata, path.folder = folder,
    filename = "fsc", verbose = FALSE, chunk.size = 1)
  expect_equal(sum(x$spectrum), 3)
  expect_equal(x$spectrum[1], 2)
  expect_equal(x$spectrum[5], 1)
  expect_equal(x$parameters$INVARIANT_RECORDS_KEPT, 1)
  expect_equal(x$parameters$OBSERVED_INVARIANT_KEPT, 2)
  expect_equal(x$audit$STATUS, c("KEPT", "INCOMPLETE_CALLS", "KEPT", "KEPT"))
  y <- write_dadi(f$gds, strata = f$strata, path.folder = folder,
    filename = "dadi", verbose = FALSE)
  expect_equal(y$spectrum, x$spectrum)
  expect_equal(which(y$mask), c(1L, 6L, 8L, 9L))
  lines <- readLines(y$files[1])
  expect_equal(lines[1], '3 3 folded "P1" "P2"')
  expect_equal(scan(text = lines[2], quiet = TRUE), y$spectrum)
  expect_equal(scan(text = lines[3], quiet = TRUE), as.integer(y$mask))
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_dadi(f$gds, strata = f$strata, path.folder = folder,
    filename = "dadi", verbose = FALSE), "Output exists")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  u <- write_dadi(f$gds, strata = f$strata, path.folder = folder,
    filename = "unfolded", verbose = FALSE, folded = FALSE,
    ancestral.alleles = data.frame(VARIANT_ID = 1:4,
      ANCESTRAL = c("G", "C", "C", "A")))
  expect_equal(u$spectrum[9], 2) # invariant REF differs from ancestral + fixed ALT
  expect_equal(u$spectrum[5], 1)
  expect_equal(which(u$mask), c(1L, 9L))
  expect_error(write_dadi(f$gds, strata = f$strata, folded = FALSE), "ancestral|ANCESTRAL")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  gdsfmt::add.gdsn(f$gds, "position", c(10L, 20L, 30L, 30L), replace = TRUE)
  SeqArray::seqClose(f$gds)
  f$gds <- SeqArray::seqOpen(f$path, readonly = FALSE)
  d <- write_fastsimcoal(f$gds, strata = f$strata, path.folder = folder,
    filename = "duplicates", verbose = FALSE)
  expect_equal(d$audit$STATUS[3:4], rep("DUPLICATE_POSITION", 2))
  expect_equal(sum(d$spectrum), 1)
  gdsfmt::add.gdsn(f$gds, "allele",
    c("A,<NON_REF>", "C", "C,T", "A,G"), replace = TRUE)
  expect_error(write_dadi(f$gds, strata = f$strata, verbose = FALSE), "gVCF")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
})
