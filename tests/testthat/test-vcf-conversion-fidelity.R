test_that("GT validation checks later records without inventing calls", {
  path <- tempfile(fileext = ".vcf")
  on.exit(unlink(path), add = TRUE)
  header <- c("##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts")
  good <- "chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT:PL\t0/1:."
  writeLines(c(header, good), path)
  expect_true(.validate_vcf_gt_records(path, chunk.size = 1L))
  writeLines(c(header, good,
    "chr1\t2\t.\tA\tG\t.\tPASS\t.\tDP\t10"), path)
  expect_error(.validate_vcf_gt_records(path, chunk.size = 1L),
    "record 2.*chr1:2.*lacks a leading GT")
  folder <- tempfile()
  expect_error(read_vcf(path, path.folder = folder, verbose = FALSE),
    "lacks a leading GT")
  expect_false(dir.exists(folder))
  compressed <- tempfile(fileext = ".vcf.gz")
  on.exit(unlink(compressed), add = TRUE)
  con <- gzfile(compressed, "wt")
  writeLines(c(header, good), con); close(con)
  expect_true(.validate_vcf_gt_records(compressed))
})

test_that("read_vcf preserves chromosome names, GT and FORMAT cardinality", {
  skip_if(!nzchar(Sys.which("bcftools")), "bcftools unavailable")
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  path <- file.path(folder, "input.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb",
    "chr1\t10\tm1\tA\tG\t50\tPASS\t.\tGT:DP:AD:PL\t0/0:10:10,0:0,30,60\t0/1:12:6,6:30,0,30",
    "chr1\t20\tm2\tA\tG\t50\tPASS\t.\tGT:PL\t0/0:.\t0/1:."), path)
  g <- read_vcf(path, filename = "test", path.folder = file.path(folder, "import"),
    parallel.core = 1, verbose = FALSE, vcf.metadata = TRUE,
    filter.haplotype.format = FALSE)
  on.exit(SeqArray::seqClose(g), add = TRUE, after = FALSE)
  expect_equal(SeqArray::seqGetData(g, "chromosome"), rep("chr1", 2))
  calls <- SeqArray::seqGetData(g, "genotype")
  expect_equal(as.integer(calls[, , 2]), c(0L, 0L, 0L, 1L))
  out <- file.path(folder, "output.vcf")
  write_vcf(g, filename = out, verbose = FALSE)
  text <- readLines(out)
  expect_true(any(grepl("ID=AD,Number=R", text, fixed = TRUE)))
  expect_true(any(grepl("ID=PL,Number=G", text, fixed = TRUE)))
  expect_true(any(grepl("chr1\t20\tm2\tA\tG\t50\tPASS", text, fixed = TRUE)))
  expect_true(any(grepl("0/0:.\t0/1:.", text, fixed = TRUE)))
})
