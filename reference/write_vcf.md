# Write a vcf file from a tidy data frame

Write a vcf file (file format version 4.3, see details below) from a
tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users. This writer translates the supplied
calls; it does not apply VCF quality filters. Resolve marker, sample,
genotype-quality, and missingness filtering before export according to
the downstream use of the VCF.

## Usage

``` r
write_vcf(data, strata = FALSE, filename = NULL, source = NULL, empty = FALSE)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- strata:

  (optional, logical) Should the strata information be included in the
  FORMAT field (along the GT info for each samples ?). To make the VCF
  population-ready use `strata = TRUE`. The strata information must be
  included in the `STRATA` column of the tidy dataset. Default:
  `strata = FALSE`. Experimental.

- filename:

  (optional) The file name prefix for the vcf file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `radiator_vcf_file_`. Default: `filename = NULL`.

- source:

  source of vcf Default: `source = NULL`.

- empty:

  generate an empty vcf. Default: `empty = FALSE`.

## Details

**VCF file format version:**

If you need a different file format version than the current one, just
change the version inside the newly created VCF, that should do the
trick. [For more information on Variant Call Format
specifications](https://vcftools.github.io/specs.html).

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## References

Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and
VCFtools. Bioinformatics, 27, 2156-2158.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
