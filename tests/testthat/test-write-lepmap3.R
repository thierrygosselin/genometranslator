lepmap_fixture <- function(folder) {
  vcf <- file.path(folder, "input.vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Likelihoods">',
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Log10 likelihoods">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tdad\tmom\tc1\tc2",
    paste("chr1", 10, "s1", "A", "G", ".", "PASS", ".", "GT:PL:GL",
      "0/1:30,0,30:-3,0,-3", "0/0:0,30,60:0,-3,-6",
      "0/1:30,0,30:-3,0,-3", "./.:.:.", sep = "\t"),
    paste("chr1", 20, "s2", "C", "T", ".", "PASS", ".", "GT:PL:GL",
      "0/0:0,30,60:0,-3,-6", "0/1:30,0,30:-3,0,-3",
      "0/0:0,30,60:0,-3,-6", "0/1:30,0,30:-3,0,-3", sep = "\t")), vcf)
  path <- file.path(folder, "input.gds")
  SeqArray::seqVCF2GDS(vcf, path, ignore.chr.prefix = "", verbose = FALSE)
  path
}
lepmap_ped <- function() data.frame(FAMILY = "F1",
  INDIVIDUALS = c("c2", "mom", "dad", "c1"),
  FATHER = c("dad", "0", "0", "dad"),
  MOTHER = c("mom", "0", "0", "mom"), SEX = c(0, 2, 1, 0))

test_that("Lep-MAP3 bundle matches IDs, retains likelihoods and restores filters", {
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  gds <- SeqArray::seqOpen(lepmap_fixture(folder))
  on.exit(SeqArray::seqClose(gds), add = TRUE)
  before <- SeqArray::seqGetFilter(gds)
  for (field in c("auto", "PL", "GL")) {
    paths <- write_lepmap3(gds, lepmap_ped(), field, filename = field,
      path.folder = folder, verbose = FALSE)
    expect_true(all(file.exists(paths)))
    expect_identical(SeqArray::seqGetFilter(gds), before)
    ped <- readLines(paths[["pedigree"]])
    expect_length(ped, 6)
    expect_equal(ped[2], "CHR\tPOS\tdad\tmom\tc1\tc2")
    selected <- if (field == "auto") "PL" else field
    lines <- readLines(paths[["vcf"]])
    records <- lines[!startsWith(lines, "#")]
    expect_true(all(grepl(paste0("GT:", selected), records, fixed = TRUE)))
    expect_false(any(grepl(paste0("ID=", if (selected == "PL") "GL" else "PL"), lines)))
    expect_match(records[1], "\\./\\.:\\.$")
    expect_match(records[1], "^chr1\\t10")
    expect_error(write_lepmap3(gds, lepmap_ped(), filename = field,
      path.folder = folder, verbose = FALSE), "Output exists")
    java <- Sys.getenv("LEPMAP3_JAVA")
    bin <- Sys.getenv("LEPMAP3_BIN")
    if (nzchar(java) && nzchar(bin)) {
      err <- tempfile(); out <- tempfile()
      status <- system2(java, c("-cp", shQuote(bin), "ParentCall2",
        shQuote(paste0("data=", paths[["pedigree"]])),
        shQuote(paste0("vcfFile=", paths[["vcf"]]))), stdout = out, stderr = err)
      expect_equal(status, 0)
      expect_false(any(grepl("Error|Exception|strange", readLines(err))))
      expect_gte(length(readLines(out)), 8)
      unlink(c(err, out))
    }
  }
  SeqArray::seqSetFilter(gds, sample.id = c("dad", "mom", "c1"),
    variant.id = 2L, verbose = FALSE)
  selected <- SeqArray::seqGetFilter(gds)
  expect_error(write_lepmap3(gds, lepmap_ped(), path.folder = folder,
    verbose = FALSE), "match active samples")
  expect_identical(SeqArray::seqGetFilter(gds), selected)
  p <- lepmap_ped(); p <- p[p$INDIVIDUALS != "c2", ]
  paths <- write_lepmap3(gds, p, filename = "subset", path.folder = folder,
    verbose = FALSE)
  expect_length(readLines(paths[["vcf"]])[!startsWith(readLines(paths[["vcf"]]), "#")], 1)
  expect_identical(SeqArray::seqGetFilter(gds), selected)
})

test_that("Lep-MAP3 rejects inconsistent pedigrees", {
  p <- lepmap_ped(); ids <- c("dad", "mom", "c1", "c2")
  expect_equal(.lepmap3_pedigree(p, ids)$INDIVIDUALS, ids)
  bad <- p; bad$SEX[3] <- 2
  expect_error(.lepmap3_pedigree(bad, ids), "matching SEX")
  bad <- p; bad$FATHER[1] <- "absent"
  expect_error(.lepmap3_pedigree(bad, ids), "Every named parent")
  bad <- p; bad$FATHER[3] <- "dad"
  expect_error(.lepmap3_pedigree(bad, ids), "own parent")
  bad <- p; bad$FAMILY[3] <- "F2"
  expect_error(.lepmap3_pedigree(bad, ids), "share")
  bad <- p; bad$SEX[4] <- 1; bad$FATHER[3] <- "c1"
  expect_error(.lepmap3_pedigree(bad, ids), "cycle")
  bad <- p; bad$FATHER <- "0"; bad$MOTHER <- "0"
  expect_error(.lepmap3_pedigree(bad, ids), "relationship")
})

test_that("Lep-MAP3 validates actual vectors without a GT fallback", {
  source <- tempfile(); dest <- tempfile()
  on.exit(unlink(c(source, dest)), add = TRUE)
  header <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta"
  check <- function(value, field = "PL") {
    writeLines(c(header, paste0("1\t1\tx\tA\tC\t.\tPASS\t.\tGT:", field,
      "\t", value)), source)
    .lepmap3_vcf(source, dest, "a", field)
  }
  expect_error(check("0/1:0,.,3"), "Invalid PL")
  expect_error(check("0/1:0,3"), "Invalid PL")
  expect_error(check("0/1:-1,0,3"), "Invalid PL")
  expect_error(check("0/1:."), "no observed PL")
  expect_error(check("1:0,3,5"), "Diploid")
  expect_error(check("0/1:0,Inf,-3", "GL"), "Invalid GL")
})
