test_that("COLONY exports GDS calls and requested settings faithfully", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  out <- write_colony(f$gds, filename = "test.dat", path.folder = folder,
    inbreeding = 1, mating.sys.males = 1, mating.sys.females = 1,
    clone = 1, run.length = 4, analysis = 2, allelic.dropout = .03,
    error.rate = .15, random.seed = 987, num.runs = 3, verbose = FALSE,
    chunk.size = 1)
  z <- readLines(out[["data"]])
  expect_equal(z[3:14], c("4", "3", "987", "0", "2", "1", "0",
    "1 1", "1", "1", "0", "0"))
  expect_equal(z[15:23], c("3", "4", "0", "10000", "0", "2", "3",
    "L1 L2 L3", "0@"))
  expect_equal(as.numeric(sub("@", "", z[24:25])), c(.03, .15))
  expect_equal(z[26:29], c("S1 1 1 1 1 2 3", "S2 1 2 1 1 0 0",
    "S3 2 2 1 2 1 2", "S4 1 2 1 1 1 1"))
  expect_equal(z[30:39], c("0 0", "0 0", "0 0", "0 0", rep("0", 6)))
  expect_true(all(file.exists(out)))
  expect_equal(readr::read_tsv(out[["samples"]], show_col_types = FALSE)$INDIVIDUALS,
    paste0("s", 1:4))
  expect_equal(SeqArray::seqGetData(f$gds, "sample.id"), paste0("s",1:4))
  expect_equal(SeqArray::seqGetData(f$gds, "variant.id"), 1:3)
  expect_error(write_colony(f$gds, filename = "test", path.folder = folder,
    verbose = FALSE), "Output exists")
})

test_that("COLONY handles metadata roles and complete frequency records", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  meta <- f$strata
  meta$ROLE <- c(" Offspring ", "offspring", "Father", "mother")
  out <- write_colony(f$gds, strata = meta[4:1, ], role.column = "ROLE",
    probability.father = .6, probability.mother = .4, allele.freq = "overall",
    filename = "parents", path.folder = folder, verbose = FALSE)
  z <- readLines(out[1])
  expect_equal(z[3], "2")
  expect_equal(z[14:15], c("1", "2 2 3"))
  expect_equal(z[16], "1 2")
  expect_equal(scan(text = z[17], quiet = TRUE), c(.5,.5))
  expect_equal(scan(text = z[19], quiet = TRUE), c(.875,.125))
  expect_equal(scan(text = z[21], quiet = TRUE), c(.5,1/3,1/6))
  expect_equal(z[35:38], c("0.6 0.4", "1 1",
    "S3 2 2 1 2 1 2", "S4 1 2 1 1 1 1"))
  expect_error(write_colony(f$gds, strata = meta, role.column = "ROLE",
    verbose = FALSE), "probability.father")
  expect_error(write_colony(f$gds, strata = meta, allele.freq = "Population A",
    verbose = FALSE), "lacks an allele")
  expect_error(write_colony(f$gds, strata = meta[-1, ], verbose = FALSE),
    "cover every")
  expect_equal(SeqArray::seqGetData(f$gds, "sample.id"), paste0("s",1:4))
})

test_that("COLONY sampling is reproducible and does not change RNG or filters", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  set.seed(8); before <- .Random.seed
  a <- write_colony(f$gds, sample.markers = 2, random.seed = 7,
    filename = "a", path.folder = folder, verbose = FALSE)
  expect_identical(.Random.seed, before)
  set.seed(9)
  copied <- tempfile(fileext = ".gds")
  file.copy(f$path, copied)
  on.exit(unlink(copied), add = TRUE)
  b <- write_colony(copied, sample.markers = 2, random.seed = 7,
    filename = "b", path.folder = folder, verbose = FALSE)
  expect_identical(readLines(a[1]), readLines(b[1]))
  expect_identical(readLines(a[3]), readLines(b[3]))
  expect_equal(SeqArray::seqGetData(f$gds, "variant.id"), 1:3)
  expect_error(write_colony(f$gds, sample.markers = 4, verbose = FALSE),
    "exceeds")
  expect_error(write_colony(f$gds, error.rate = 1, verbose = FALSE), "probability")
  expect_error(write_colony(f$gds, inbreeding = 2, verbose = FALSE), "ranges")
  expect_error(write_colony(f$gds, random.seed = NA, verbose = FALSE), "whole")
  expect_error(write_colony(f$strata, verbose = FALSE), "GDS")
  SeqArray::seqSetFilter(f$gds, sample.id = "s1", variant.id = 1,
    verbose = FALSE)
  expect_error(write_colony(f$gds, verbose = FALSE), "two observed alleles")
  expect_equal(SeqArray::seqGetData(f$gds, "sample.id"), "s1")
})

test_that("COLONY retains rare allele frequencies and rejects partial calls", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  vcf <- file.path(folder, "rare.vcf"); path <- file.path(folder, "rare.gds")
  header <- c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", paste0("sample", 1:200)), collapse = "\t"))
  record <- function(calls) paste(c("1", "1", "m1", "A", "G", ".",
    "PASS", ".", "GT", calls), collapse = "\t")
  writeLines(c(header, record(c("0/1", rep("0/0", 199)))), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  out <- write_colony(path, filename = "rare", path.folder = folder,
    allele.freq = "overall", verbose = FALSE)
  expect_equal(scan(text = readLines(out[1])[17], quiet = TRUE),
    c(399/400, 1/400))
  # A failed export must not publish any files.
  partial <- file.path(folder, "partial.gds")
  writeLines(c(header, record(c("./1", "0/1", rep("0/0", 198)))), vcf)
  SeqArray::seqVCF2GDS(vcf, partial, verbose = FALSE)
  expect_error(write_colony(partial, filename = "bad", path.folder = folder,
    verbose = FALSE), "Partially missing")
  expect_false(file.exists(file.path(folder, "bad.dat")))
})
