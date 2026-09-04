test_that("VCF encodings preserve GDS calls and annotations", {
  exe <- Sys.which("bcftools")
  skip_if(!nzchar(exe), "bcftools unavailable")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  vcf <- file.path(folder, "input.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
    "1\t10\ti1\tA\t.\t.\tPASS\t.\tGT\t0/0\t./.",
    "1\t20\ts1\tC\tT\t40\tPASS\t.\tGT:DP:AD:PL\t0/0:10:10,0:0,30,60\t1|0:10:5,5:30,0,30"), vcf)
  path <- file.path(folder, "input.gds")
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  gds <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(gds), add = TRUE)
  SeqArray::seqSetFilter(gds, sample.id = "b", verbose = FALSE)
  before <- SeqArray::seqGetFilter(gds)
  queries <- list()
  for (ext in c("vcf", "vcf.gz", "bcf")) {
    file <- file.path(folder, paste0("output.", ext))
    result <- write_vcf(gds, filename = file, verbose = FALSE)
    expect_true(all(file.exists(result)))
    expect_false(file.exists(paste0(file, ".vcf")))
    queries[[ext]] <- bcftools_exec(exe, c("query", "-f",
      "%POS\t%REF\t%ALT[\t%GT\t%DP\t%AD\t%PL]\n", file), verbose = FALSE)$stdout
    expect_identical(SeqArray::seqGetFilter(gds), before)
    expect_error(write_vcf(gds, filename = file, verbose = FALSE), "Output exists")
    if (ext != "vcf") {
      expect_equal(bcftools_exec(exe, c("index", "-n", file), verbose = FALSE)$stdout, "2\n")
    }
  }
  expect_equal(queries$vcf, queries$bcf)
  expect_equal(queries$vcf, queries$vcf.gz)
  expect_match(queries$bcf, "1\\|0\\t10\\t5,5\\t30,0,30")
  expect_match(queries$bcf, "10\\tA\\t\\.\\t\\./\\.")
  expect_error(write_vcf(gds, filename = file.path(folder, "bad.vcf"), index = TRUE),
    "Indexing requires")
  expect_identical(SeqArray::seqGetFilter(gds), before)
  legacy <- data.frame(MARKERS = "m1", CHROM = "1", LOCUS = "l1", POS = 10L,
    REF = "001", ALT = "003", INDIVIDUALS = "a", GT_VCF = "0/1")
  legacy.file <- file.path(folder, "legacy.bcf")
  write_vcf(legacy, filename = legacy.file, verbose = FALSE)
  expect_equal(bcftools_exec(exe, c("query", "-f", "%REF\t%ALT\n", legacy.file),
    verbose = FALSE)$stdout, "A\tG\n")
})
