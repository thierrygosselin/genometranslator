# SeqArray's genotype storage needs GT on every sample-bearing record.
# Scan in bounded chunks before creating outputs; never infer calls from DP/PL.
.validate_vcf_gt_records <- function(path, chunk.size = 1000L) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path))
    stop("data must name an existing VCF file.", call. = FALSE)
  con <- gzfile(path, "rt")
  on.exit(close(con), add = TRUE)
  record <- 0L
  repeat {
    lines <- readLines(con, n = chunk.size, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    for (line in lines) {
      record <- record + 1L
      fields <- stringi::stri_split_fixed(line, "\t", n = 10L)[[1]]
      if (length(fields) > 9L &&
          !startsWith(paste0(fields[9], ":"), "GT:")) {
        stop(paste0("VCF record ", record, " (", fields[1], ":", fields[2],
          ") lacks a leading GT field.\n",
          "This GDS importer requires GT on every sample-bearing record.\n",
          "Retain the original VCF. In a separate copy, add explicitly missing\n",
          "GT calls using known sample ploidy, or perform upstream genotyping.\n",
          "Do not infer genotype calls from depth alone."), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}
