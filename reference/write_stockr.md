# Write a stockR dataset from a tidy data frame or GDS file or object.

Write a stockR dataset (Fost et al. submitted). Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_stockr(data, filename = NULL, verbose = TRUE)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- filename:

  (optional) The stockr object is written in the working directory. The
  file is written with `radiator_stockr_DATE@TIME.RData` and can be open
  with readRDS. Default: `filename = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

## Value

The object generated is a matrix with dimension: MARKERS x INDIVIDUALS.
The genotypes are coded like PLINK: 0, 1 or 2 alternate allele. 0:
homozygote for the reference allele, 1: heterozygote, 2: homozygote for
the alternate allele. Missing genotypes have NA. The object also as 2
attributes. `attributes(data)$grps` with `STRATA/POP_ID` of the
individuals and `attributes(data)$sample.grps` filled with
`INDIVIDUALS`. Both attributes can be used inside *stockR*.

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

Foster et al. submitted

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
