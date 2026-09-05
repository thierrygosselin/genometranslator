test_that("OutFLANK components preserve equal frequencies and allele symmetry", {
  expect_equal(.outflank_components(matrix(2, 2, 3)),
    c(He=.5, FST=-2/15, T1=-1/30, T2=.25, FSTNoCorr=0,
      T1NoCorr=0, T2NoCorr=.25, meanAlleleFreq=.5))
  z <- rbind(c(3,2,1), c(1,2,3))
  expect_equal(.outflank_components(z),
    c(He=.5, FST=.1, T1=1/36, T2=5/18, FSTNoCorr=.2,
      T1NoCorr=1/18, T2NoCorr=5/18, meanAlleleFreq=.5))
  z <- rbind(c(5,2,1), c(1,3,2), c(1,1,3))
  x <- .outflank_components(z); y <- .outflank_components(z[,3:1])
  expect_equal(x[1:7], y[1:7])
  expect_equal(x[8] + y[8], c(meanAlleleFreq=1))
})

test_that("OutFLANK export is chunk invariant and matches samples by ID", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add=TRUE)
  SeqArray::seqSetFilter(f$gds, variant.id=1:2, verbose=FALSE)
  before <- SeqArray::seqGetFilter(f$gds)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive=TRUE), add=TRUE)
  x <- write_outflank(f$gds, strata=f$strata[4:1,], filename="one",
    path.folder=folder, chunk.size=1, verbose=FALSE)
  y <- write_outflank(f$gds, strata=f$strata, filename="two",
    path.folder=folder, chunk.size=20, verbose=FALSE)
  expect_equal(x$FstDataFrame, y$FstDataFrame)
  expect_equal(x$population.audit, y$population.audit)
  expect_equal(x$samples$INDIVIDUALS, paste0("s",1:4))
  expect_equal(x$loci$VARIANT_ID, 1:2)
  expect_equal(unname(unlist(x$FstDataFrame[1,-1])),
    unname(.outflank_components(rbind(c(1,1,0),c(0,1,1)))))
  expect_equal(read.delim(x$files[["statistics"]], check.names=FALSE,
    colClasses=c(LocusName="character")),
    x$FstDataFrame)
  expect_true(all(file.exists(x$files)))
  expect_true(all(grepl("_[0-9]{8}@[0-9]{4}\\.tsv$", x$files)))
  expect_false("overwrite" %in% names(formals(write_outflank)))
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_outflank(f$gds, strata=f$strata, filename="one",
    path.folder=folder, verbose=FALSE), "Output exists")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
})

test_that("OutFLANK exclusions are explicit and audited", {
  f <- export_fixture(incomplete=TRUE)
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add=TRUE)
  SeqArray::seqSetFilter(f$gds, variant.id=1:2, verbose=FALSE)
  before <- SeqArray::seqGetFilter(f$gds)
  folder <- tempfile(); dir.create(folder)
  on.exit(unlink(folder, recursive=TRUE), add=TRUE)
  expect_error(write_outflank(f$gds, strata=f$strata,
    path.folder=folder, verbose=FALSE), "estimability")
  expect_length(list.files(folder), 0)
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  x <- write_outflank(f$gds, strata=f$strata, invalid.loci="exclude",
    path.folder=folder, verbose=FALSE)
  expect_equal(nrow(x$FstDataFrame), 1)
  expect_equal(x$loci$STATUS, c("exported", "excluded"))
  expect_equal(x$loci$REASON[2], "insufficient_called_per_population")
  expect_equal(x$population.audit$N_CALLED, c(2,2,1,2))
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
})

test_that("OutFLANK rejects unsupported inputs and insufficient groups", {
  f <- export_fixture()
  on.exit({SeqArray::seqClose(f$gds); unlink(f$path)}, add=TRUE)
  before <- SeqArray::seqGetFilter(f$gds)
  expect_error(write_outflank(f$gds, strata=f$strata, verbose=FALSE), "biallelic")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_outflank(f$gds, strata=f$strata, min.called=1,
    verbose=FALSE), "min.called")
  expect_error(write_outflank(f$gds, strata=f$strata, pop.select="Population A",
    verbose=FALSE), "two populations")
  expect_error(write_outflank(f$gds, strata=f$strata[-1,],
    verbose=FALSE), "every active sample")
  expect_identical(SeqArray::seqGetFilter(f$gds), before)
  expect_error(write_outflank(matrix(0,2,2), verbose=FALSE), "GDS")
})

test_that("OutFLANK path input handles zero FST and audits unusable loci", {
  vcf <- tempfile(fileext=".vcf"); path <- tempfile(fileext=".gds")
  folder <- tempfile(); dir.create(folder)
  on.exit({unlink(c(vcf,path)); unlink(folder,recursive=TRUE)}, add=TRUE)
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "1\t10\ta\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\t0/0\t1/1",
    "1\t20\tb\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0",
    "1\t30\tc\tA\tG\t.\tPASS\t.\tGT\t./.\t./.\t0/0\t1/1",
    "1\t40\td\tA\tG\t.\tPASS\t.\tGT\t./.\t./.\t./.\t./."), vcf)
  SeqArray::seqVCF2GDS(vcf,path,verbose=FALSE)
  meta <- data.frame(INDIVIDUALS=paste0("s",1:4), STRATA=rep(c("A","B"),each=2))
  x <- write_outflank(path, strata=meta, invalid.loci="exclude",
    filename="test",path.folder=folder,verbose=FALSE)
  expect_equal(x$FstDataFrame$He,.5)
  expect_equal(x$FstDataFrame$FSTNoCorr,0)
  expect_equal(x$loci$REASON,c("", "monomorphic",
    "insufficient_called_per_population", "insufficient_called_per_population"))
  g <- SeqArray::seqOpen(path)
  on.exit(SeqArray::seqClose(g),add=TRUE,after=FALSE)
  SeqArray::seqSetFilter(g,variant.id=2:4,verbose=FALSE)
  before <- SeqArray::seqGetFilter(g)
  expect_error(write_outflank(g,strata=meta,invalid.loci="exclude",
    filename="none",path.folder=folder,verbose=FALSE), "No estimable")
  expect_identical(SeqArray::seqGetFilter(g),before)
  expect_false(any(grepl("none",list.files(folder))))
})

test_that("OutFLANK refuses partial calls without publishing files", {
  vcf <- tempfile(fileext=".vcf"); path <- tempfile(fileext=".gds")
  folder <- tempfile(); dir.create(folder)
  on.exit({unlink(c(vcf,path)); unlink(folder,recursive=TRUE)},add=TRUE)
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=1,length=1000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "1\t10\ta\tA\tG\t.\tPASS\t.\tGT\t0/.\t1/1\t0/0\t1/1"),vcf)
  SeqArray::seqVCF2GDS(vcf,path,verbose=FALSE)
  meta <- data.frame(INDIVIDUALS=paste0("s",1:4),STRATA=rep(c("A","B"),each=2))
  expect_error(write_outflank(path,strata=meta,path.folder=folder,
    verbose=FALSE), "Partial")
  expect_length(list.files(folder),0)
})
