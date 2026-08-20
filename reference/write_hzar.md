# Write a HZAR file from a tidy data frame.

Write a HZAR file (Derryberry et al. 2013), from a tidy data frame. Used
internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_hzar(
  data,
  distances = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- distances:

  (optional) A file with 2 columns, `POP_ID` and the distance
  information per populations. With default: `distances = NULL`, the
  column is left empty. Default: `distances = NULL`.

- filename:

  (optional) The file name prefix for the HZAR file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `radiator_hzar_`. Default: `filename = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

## Value

A HZAR file is written in the working directory.

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

Derryberry EP, Derryberry GE, Maley JM, Brumfield RT. hzar: hybrid zone
analysis using an R software package. Molecular Ecology Resources.
2013;14: 652-663. doi:10.1111/1755-0998.12209

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# The simplest form of the function:
hzar.data <- genometranslator::write_hzar(data = tidydata)


# Using genepop dataset, nancycats, from adegenet
# require(adegenet)
nancycats <- system.file("files/nancycats.gen", package = "adegenet")


# using genome_translator:
nanycats.hzar <- genometranslator::genome_translator(data = nancycats, output = "hzar")


# using the separate modules:
# tidy the genepop file then pipe the result in write_hzar
nanycats.hzar <- genometranslator::read_genepop(data = nancycats) %>%
    genometranslator::write_hzar(data = .)
} # }
```
