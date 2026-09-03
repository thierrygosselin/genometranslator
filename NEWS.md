# genometranslator 0.0.0.9000

* Removed the obsolete `tidy_vcf()` wrapper. Use `read_genome()` for automatic
  format detection, `read_vcf()` for VCF-specific import controls, and
  `tidy_genome()` only when an in-memory table is required. Genomic filtering
  remains the responsibility of `radr`; the historical wrapper remains in
  `radiator` for legacy workflows.

* First commit: initial version of genometranslator




