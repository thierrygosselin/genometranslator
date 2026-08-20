# Write a related file from a tidy data frame

Write a related file from a tidy data frame. This output file format
enables to run the data in the
[related](https://github.com/timothyfrasier/related) R package (Pew et
al. 2015), which is essantially the R version of
[COANCESTRY](https://www.zsl.org/science/software/coancestry) fortran
program developed by Jinliang Wang. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_related(
  data,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- filename:

  (optional) The file name prefix for the related file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `radiator_related_`. Default: `filename = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

- ...:

  other parameters passed to the function.

## Value

A related file is saved to the working directory.

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

Pew J, Muir PH, Wang J, Frasier TR (2015) related: an R package for
analysing pairwise relatedness from codominant molecular markers.
Molecular Ecology Resources, 15, 557-561.

Wang, J. 2011. COANCESTRY: A program for simulating, estimating and
analysing relatedness and inbreeding coefficients. Molecular Ecology
Resources 11(1): 141-145.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
