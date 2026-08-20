# Write a faststructure file from a tidy data frame

Write a faststructure file from a tidy data frame. For biallelic dataset
only. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users. Prepare quality-controlled biallelic markers and
consider missingness, minor-allele frequency, and linkage disequilibrium
before export.

## Usage

``` r
write_faststructure(
  data,
  plink.bed = FALSE,
  pop.levels = NULL,
  filename = NULL,
  ...
)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- plink.bed:

  To export in plink BED format. This will write 3 files: `.bed`,
  `.bim`, `.fam`. For this option to work, the argument data must be a
  GDS file or object. Default: `plink.bed = FALSE`.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. White spaces in population
  names are replaced by underscore. Default: `pop.levels = NULL`.

- filename:

  (optional) The file name prefix for the faststructure file written to
  the working directory. With default: `filename = NULL`, the date and
  time is appended to `radiator_faststructure_`. Default:
  `filename = NULL`.

- ...:

  other parameters passed to the function.

## Value

A faststructure file is saved to the working directory.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

The standard text output has no additional R-package dependency. Setting
`plink.bed = TRUE` requires the optional Bioconductor package
[SNPRelate](https://bioconductor.org/packages/SNPRelate). This writer
prepares input files; it does not install or run fastStructure.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## References

Raj A, Stephens M, Pritchard JK (2014) fastSTRUCTURE: Variational
Inference of Population Structure in Large SNP Datasets. Genetics, 197,
573-589.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
