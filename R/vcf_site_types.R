# Record classification is independent of allele frequencies in selected samples.
# An ALT-bearing record may become monomorphic after sample selection.
.vcf_site_types <- function(alleles) {
  vapply(strsplit(alleles, ",", fixed = TRUE), function(a) {
    if (any(a %in% c("<NON_REF>", "<*>"))) return("GVCF_REFERENCE")
    if (any(grepl("<|>|\\*|\\[|\\]", a))) return("SYMBOLIC")
    if (length(a) >= 1L && a[1L] %in% c("A", "C", "G", "T") &&
        (length(a) == 1L || all(a[-1L] %in% c("", "."))))
      return("INVARIANT_RECORD")
    if (length(a) >= 2L && all(a %in% c("A", "C", "G", "T")))
      return("SNP")
    "OTHER"
  }, character(1), USE.NAMES = FALSE)
}
