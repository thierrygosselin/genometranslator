# Write an arlequin file from a tidy data frame

Write a arlequin file from a tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users.

## Usage

``` r
write_arlequin(data, pop.levels = NULL, filename = NULL, ...)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- pop.levels:

  (optional, string) A character string with your populations ordered.
  Default: `pop.levels = NULL`.

- filename:

  (optional) The file name prefix for the arlequin file written to the
  working directory. With default: `filename = NULL`, the filename
  generated follow this `radiator_arlequin_DATE@TIME.csv`. Default:
  `filename = NULL`.

- ...:

  other parameters passed to the function.

## Value

An arlequin file is saved to the working directory.

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

Excoffier, L.G. Laval, and S. Schneider (2005) Arlequin ver. 3.0: An
integrated software package for population genetics data analysis.
Evolutionary Bioinformatics Online 1:47-50.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
