# Write a genind object from a tidy data frame or GDS file or object.

Write a genind object from a tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users.

## Usage

``` r
write_genind(data, write = FALSE, verbose = FALSE)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- write:

  (logical, optional) To write in the working directory the genind
  object. The file is written with `radiator_genind_DATE@TIME.RData` and
  can be open with load or readRDS. Default: `write = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

Reading and writing `genind` objects requires the optional CRAN package
[adegenet](https://cran.r-project.org/package=adegenet).

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## References

Jombart T (2008) adegenet: a R package for the multivariate analysis of
genetic markers. Bioinformatics, 24, 1403-1405.

Jombart T, Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of
genome-wide SNP data. Bioinformatics, 27, 3070-3071.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
