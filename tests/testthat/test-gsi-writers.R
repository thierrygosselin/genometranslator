test_that("GSI writers preserve diploid alleles and active selections", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  before <- SeqArray::seqGetFilter(f$gds)
  r <- write_rubias(f$gds, strata = f$strata, chunk.size = 1,
    filename = "r", path.folder = folder, verbose = FALSE)
  expect_equal(names(r)[1:4], c("sample_type", "repunit", "collection", "indiv"))
  expect_equal(r$L1, c(1L, 1L, 2L, 1L))
  expect_equal(r$L1.2, c(1L, 2L, 2L, 2L))
  expect_equal(r$L3, c(2L, NA, 1L, 1L))
  expect_equal(r$L3.2, c(3L, NA, 2L, 1L))
  expect_true(all(file.exists(attr(r, "files"))))
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  x <- write_gsi_sim(f$gds, strata = f$strata, filename = "g.txt",
    path.folder = folder, chunk.size = 1, verbose = FALSE,
    pop.levels = rev(unique(f$strata$STRATA)))
  lines <- readLines(x$files[1])
  expect_equal(lines[1:4], c("4 3", "L1", "L2", "L3"))
  expect_true("S3 2 2 1 2 1 2" %in% lines)
  expect_true("S2 1 2 1 1 0 0" %in% lines)
  expect_true(all(file.exists(x$files)))
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_gsi_sim(f$gds, strata = f$strata, filename = "g",
    path.folder = folder, verbose = FALSE), "Output exists")
  SeqArray::seqSetFilter(f$gds, sample.id = c("s2", "s4"),
    variant.id = 1L, verbose = FALSE)
  selected <- SeqArray::seqGetFilter(f$gds)
  r <- write_rubias(f$gds, strata = f$strata, verbose = FALSE)
  expect_equal(r$indiv, c("s2", "s4"))
  expect_equal(ncol(r), 6)
  expect_identical(SeqArray::seqGetFilter(f$gds), selected)
})

test_that("rubias validates reference and mixture design", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  m <- data.frame(indiv = paste0("s", 1:4),
    sample_type = c("reference", "reference", "mixture", "mixture"),
    repunit = c("R", "R", NA, NA), collection = c("A", "A", "M", "M"))
  r <- write_rubias(f$gds, strata = m, verbose = FALSE)
  expect_equal(r$sample_type, m$sample_type)
  expect_true(all(is.na(r$repunit[3:4])))
  bad <- m; bad$repunit[2] <- "Q"
  expect_error(write_rubias(f$gds, strata = bad, verbose = FALSE), "one repunit")
  bad <- m; bad$sample_type[1] <- NA
  expect_error(write_rubias(f$gds, strata = bad, verbose = FALSE), "sample_type")
  bad <- m; bad$repunit[3] <- "Q"
  expect_error(write_rubias(f$gds, strata = bad, verbose = FALSE), "Mixture")
  expect_error(write_rubias(f$gds, strata = m[1:3, ], verbose = FALSE), "cover")
  expect_error(write_rubias(f$gds, strata = m[, -3], verbose = FALSE), "Explicit")
  expect_error(write_rubias(data.frame(), verbose = FALSE), "GDS")
})

test_that("partial calls fail without changing selections or publishing files", {
  vcf <- tempfile(fileext = ".vcf"); path <- tempfile(fileext = ".gds")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=100>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "1\t1\tm\tA\tG\t.\tPASS\t.\tGT\t0/."), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  g <- SeqArray::seqOpen(path)
  folder <- tempfile(); dir.create(folder)
  on.exit({SeqArray::seqClose(g); unlink(c(vcf, path));
    unlink(folder, recursive = TRUE)}, add = TRUE)
  m <- data.frame(INDIVIDUALS = "s1", STRATA = "A")
  before <- SeqArray::seqGetFilter(g)
  expect_error(write_rubias(g, strata = m, verbose = FALSE), "Partially")
  expect_identical(SeqArray::seqGetFilter(g), before)
  expect_error(write_gsi_sim(g, strata = m, path.folder = folder,
    verbose = FALSE), "Partially")
  expect_length(list.files(folder), 0)
  expect_identical(SeqArray::seqGetFilter(g), before)
})
