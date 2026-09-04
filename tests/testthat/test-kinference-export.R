test_that("kinference export preserves calls, IDs and metadata", {
  skip_if_not_installed("gbasics", "1.2")
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  SeqArray::seqSetFilter(f$gds, variant.id = 1:2, verbose = FALSE)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  metadata <- f$strata[4:1, ]; metadata$plate <- 4:1
  x <- write_kinference(f$gds, strata = metadata, filename = "test",
    path.folder = folder, chunk.size = 1L, verbose = FALSE)
  expect_s3_class(x, "snpgeno")
  expect_equal(unname(as.character(x)), matrix(c("AAO","AB","BBO","AB",
    "AAO","AAO","AB","AAO"), 4, 2))
  expect_identical(x$info$INDIVIDUALS, paste0("s", 1:4))
  expect_identical(x$info$plate, 1:4)
  expect_identical(x$locinfo$A, c("A", "C"))
  expect_identical(x$locinfo$B, c("G", "T"))
  expect_identical(gbasics::rowid_field(x), "INDIVIDUALS")
  y <- readRDS(attr(x, "export.files")[1])
  expect_identical(as.character(y), as.character(x))
  expect_identical(y$info, x$info)
  expect_equal(SeqArray::seqGetData(f$gds, "variant.id"), 1:2)
  expect_equal(SeqArray::seqGetData(f$gds, "sample.id"), paste0("s",1:4))
  expect_error(write_kinference(f$gds, path.folder = folder,
    filename = "test", verbose = FALSE), "Output exists")
  expect_error(write_kinference(f$gds, strata = metadata[-1,],
    path.folder = folder, verbose = FALSE), "cover every")
  if (requireNamespace("kinference", quietly = TRUE)) {
    z <- kinference::est_ALF_nonulls(x)
    expect_equal(unname(z$locinfo$pbonzer[,"A"]), c(.5, .875))
    expect_true(all(z$locinfo$pbonzer[,"O"] == 0))
  }
})

test_that("kinference rejects missing, multiallelic and invariant selections", {
  skip_if_not_installed("gbasics", "1.2")
  f <- export_fixture(incomplete = TRUE)
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add = TRUE)
  expect_error(write_kinference(f$gds, verbose = FALSE), "biallelic")
  SeqArray::seqSetFilter(f$gds, variant.id = 1:2, verbose = FALSE)
  expect_error(write_kinference(f$gds, verbose = FALSE), "Missing calls")
  expect_equal(SeqArray::seqGetData(f$gds,"variant.id"), 1:2)
  SeqArray::seqSetFilter(f$gds, sample.id = "s1", variant.id = 1,
    verbose = FALSE)
  expect_error(write_kinference(f$gds, verbose = FALSE), "invariant")
})

test_that("exported objects reach kinference power and pair screening", {
  skip_if_not_installed("gbasics", "1.2")
  skip_if_not_installed("kinference")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  vcf <- file.path(folder,"calls.vcf"); path <- file.path(folder,"calls.gds")
  calls <- withr::with_seed(17, matrix(sample(c("0/0","0/1","1/1"),
    8000, replace=TRUE, prob=c(.25,.5,.25)), 80, 100))
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",
      paste0("s",1:80)),collapse="\t"),
    vapply(1:100,function(j) paste(c("1",j,paste0("m",j),"A","G",".",
      "PASS",".","GT",calls[,j]),collapse="\t"),"")),vcf)
  SeqArray::seqVCF2GDS(vcf,path,verbose=FALSE)
  x <- write_kinference(path,path.folder=folder,verbose=FALSE)
  x <- kinference::est_ALF_nonulls(x)
  x <- kinference::kin_power(x,k=.5)
  expect_true(all(is.finite(x$locinfo$E_HSP)))
  pairs <- suppressMessages(kinference::find_HSPs(x,keep_thresh=0,
    limit_pairs=100,ij_numeric=TRUE))
  expect_true(all(c("i","j","PLOD") %in% names(pairs)))
  expect_true(all(pairs$i %in% 1:80 & pairs$j %in% 1:80))
})
