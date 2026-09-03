make_plink_files <- function(extensions) {
  prefix <- tempfile("plink_files_")
  paths <- paste0(prefix, extensions)
  file.create(paths)
  list(prefix = prefix, paths = paths)
}

make_fake_plink2 <- function() {
  script <- tempfile("fake_plink2_")
  writeLines(c(
    "#!/bin/sh",
    "if [ \"$1\" = \"--version\" ]; then",
    "  echo 'PLINK v2.00a test'",
    "  exit 0",
    "fi",
    "out=''",
    "previous=''",
    "action=''",
    "for argument in \"$@\"; do",
    "  if [ \"$previous\" = \"--out\" ]; then out=\"$argument\"; fi",
    "  if [ \"$previous\" = \"--export\" ]; then action=\"$argument\"; fi",
    "  case \"$argument\" in",
    "    --make-pgen) action='pgen' ;;",
    "    --make-bed) action='bed' ;;",
    "  esac",
    "  previous=\"$argument\"",
    "done",
    "case \"$action\" in",
    "  pgen) touch \"${out}.pgen\" \"${out}.pvar\" \"${out}.psam\" ;;",
    "  bed) touch \"${out}.bed\" \"${out}.bim\" \"${out}.fam\" ;;",
    "  ped) touch \"${out}.ped\" \"${out}.map\" ;;",
    "  tped) touch \"${out}.tped\" \"${out}.tfam\" ;;",
    "  vcf) printf '##fileformat=VCFv4.3\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' > \"${out}.vcf\" ;;",
    "  *) echo 'unknown action' >&2; exit 2 ;;",
    "esac"
  ), script)
  Sys.chmod(script, mode = "0755")
  script
}

test_that("PLINK 1 and 2 primary files are detected", {
  pgen <- make_plink_files(c(".pgen", ".pvar", ".psam"))
  bed <- make_plink_files(c(".bed", ".bim", ".fam"))
  ped <- make_plink_files(c(".ped", ".map"))
  tped <- make_plink_files(c(".tped", ".tfam"))

  expect_identical(detect_genomic_format(pgen$paths[[1]]), "plink.pgen.file")
  expect_identical(detect_genomic_format(bed$paths[[1]]), "plink.bed.file")
  expect_identical(detect_genomic_format(ped$paths[[1]]), "plink.ped.file")
  expect_identical(detect_genomic_format(tped$paths[[1]]), "plink.tped.file")
})

test_that("compressed PVAR companion is accepted", {
  files <- make_plink_files(c(".pgen", ".pvar.zst", ".psam"))
  expect_identical(detect_genomic_format(files$paths[[1]]), "plink.pgen.file")
  args <- genometranslator:::.plink2_input_args(
    files$paths[[1]], "plink.pgen.file"
  )
  expect_true("vzs" %in% args)
})

test_that("missing PLINK companion files are reported", {
  pgen <- make_plink_files(".pgen")
  expect_error(
    detect_genomic_format(pgen$paths[[1]]),
    "Missing PLINK 2 companion file"
  )
})

test_that("PLINK 2 routes each requested export flavour", {
  fake <- make_fake_plink2()
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.3",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    "1\t1\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1"
  ), vcf)
  output <- tempfile("plink_output_")
  dir.create(output)
  output <- normalizePath(output)

  expected <- list(
    pgen = c(".pgen", ".pvar", ".psam"),
    bed = c(".bed", ".bim", ".fam"),
    ped = c(".ped", ".map"),
    tped = c(".tped", ".tfam")
  )
  for (flavour in names(expected)) {
    files <- write_plink(
      data = vcf,
      filename = flavour,
      format = flavour,
      path.folder = output,
      plink.path = fake,
      verbose = FALSE
    )
    expect_setequal(files, file.path(output, paste0(flavour, expected[[flavour]])))
  }
})

test_that("PLINK normalization builds a temporary VCF", {
  fake <- make_fake_plink2()
  files <- make_plink_files(c(".pgen", ".pvar", ".psam"))
  converted <- genometranslator:::.plink_to_vcf(
    files$paths[[1]], "plink.pgen.file", plink.path = fake, verbose = FALSE
  )
  expect_true(file.exists(converted$vcf))
  expect_match(readLines(converted$vcf, n = 1L), "VCFv4.3")
})
