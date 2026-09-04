export_fixture <- function(incomplete = FALSE) {
  vcf <- tempfile(fileext = ".vcf")
  gds <- tempfile(fileext = ".gds")
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "1\t10\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    if (incomplete) {
      "1\t20\tm2\tC\tT\t.\tPASS\t.\tGT\t./.\t0/0\t0/1\t0/0"
    } else "1\t20\tm2\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/0",
    "1\t30\tm3\tA\tC,G\t.\tPASS\t.\tGT\t1/2\t./.\t0/1\t0/0"
  ), vcf)
  SeqArray::seqVCF2GDS(vcf, gds, verbose = FALSE)
  unlink(vcf)
  list(path = gds, gds = SeqArray::seqOpen(gds),
    strata = data.frame(INDIVIDUALS = paste0("s", 1:4),
      STRATA = c("Population A", "Population A", "Population B", "Population B")))
}
