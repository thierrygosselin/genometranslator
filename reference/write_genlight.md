# Write a `genlight` object from: a tidy data frame, GDS file or object.

Write a `genlight` object from a tidy data frame or GDS file or object.
Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users. `genlight` is a formal (S4) class
for storing genotypes of binary SNPs in a compact way, using a bit-level
coding scheme. This storage is most efficient with haploid data, where
the memory taken to represent data can be reduced more than 50 times.
However, `genlight` can be used for any level of ploidy, and still
remain an efficient storage mode.

## Usage

``` r
write_genlight(
  data,
  write = FALSE,
  dartr = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 2,
  biallelic = TRUE
)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- write:

  (logical, optional) To write in the working directory the genlight
  object. The file is written with `radiator_genlight_DATE@TIME.RData`
  and can be open with load or readRDS. Default: `write = FALSE`.

- dartr:

  (logical, optional) For non-dartR users who wants to have a genlight
  object ready for the dartR package. This option transfer or generates:
  `CALL_RATE, AVG_COUNT_REF, AVG_COUNT_SNP, REP_AVG, ONE_RATIO_REF, ONE_RATIO_SNP`.
  These markers metadata are stored into the genlight slot:
  `genlight.obj@other$loc.metrics`. **Use the radiator generated GDS
  data for best result**. Default: `dartr = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 2`.

- biallelic:

  (logical, optional) If you already know that the data is biallelic use
  this argument to speed up the function. Default: `biallelic = TRUE`.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

Reading and writing `genlight` objects requires the optional CRAN
package [adegenet](https://cran.r-project.org/package=adegenet).

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

## Examples

``` r
if (FALSE) { # \dontrun{
# With defaults:
turtle <- genometranslator::write_genlight(data = "my.metadata.node.rad")

# Write gl object in directory:
turtle <- genometranslator::write_genlight(data = "my.metadata.node.rad", write = TRUE)

# Generate a dartR ready genlight and verbose = TRUE:
turtle <- genometranslator::write_genlight(
    data = "my.metadata.node.rad",
    write = TRUE,
    dartr = TRUE,
    verbose = TRUE
 )
} # }
```
