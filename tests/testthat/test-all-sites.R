test_that("site classification distinguishes invariant records and symbolic alleles", {
  expect_equal(.vcf_site_types(c("A", "A,", "A,.", "C,T", "A,C,G",
    "AT,A", "A,<NON_REF>", "A,*", "N")),
    c(rep("INVARIANT_RECORD", 3), "SNP", "SNP", "OTHER",
      "GVCF_REFERENCE", "SYMBOLIC", "OTHER"))
})

test_that("SeqArray retains invariant GT calls and missingness", {
  vcf <- tempfile(fileext = ".vcf")
  path <- tempfile(fileext = ".gds")
  on.exit(unlink(c(vcf, path)), add = TRUE)
  writeLines(c("##fileformat=VCFv4.2", "##source=Stacks v2.68",
    "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t10\t1:0:+\tA\t.\t.\tPASS\t.\tGT\t0/0\t./.",
    "chr1\t11\t1:1:+\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0"), vcf)
  SeqArray::seqVCF2GDS(vcf, path, verbose = FALSE)
  gds <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(gds), add = TRUE)
  expect_equal(.vcf_site_types(SeqArray::seqGetData(gds, "allele")),
    c("INVARIANT_RECORD", "SNP"))
  gt <- SeqArray::seqGetData(gds, "genotype")
  expect_equal(as.integer(gt[, 1, 1]), c(0L, 0L))
  expect_true(all(is.na(gt[, 2, 1])))
  folder <- tempfile(); dir.create(folder)
  imported <- read_vcf(vcf, all.sites = TRUE, parallel.core = 1,
    path.folder = folder, verbose = FALSE)
  on.exit(SeqArray::seqClose(imported), add = TRUE)
  expect_equal(length(SeqArray::seqGetData(imported, "variant.id")), 2L)
  expect_equal(extract_markers_metadata(imported)$SITE_TYPE,
    c("INVARIANT_RECORD", "SNP"))
  content <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(imported,
    genome_metadata_path(imported, "vcf.site.content")))
  expect_true(as.logical(content$ALL_SITES))
  expect_equal(as.numeric(content$INVARIANT_RECORDS), 1)
  expect_true(all(is.na(SeqArray::seqGetData(imported, "genotype")[, 2, 1])))
})
