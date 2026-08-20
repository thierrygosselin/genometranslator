# Write a genepop file

Write a genepop file from a tidy data frame or GDS file/object. Used
internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users.

## Usage

``` r
write_genepop(
  data,
  pop.levels = NULL,
  genepop.header = NULL,
  markers.line = TRUE,
  filename = NULL,
  ...
)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- pop.levels:

  (optional, string) A character string with your populations ordered.
  Default: `pop.levels = NULL`. Described in
  [`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md).

- genepop.header:

  The first line of the Genepop file. With the default, Will use
  "radiator genepop with date". Default: `genepop.header = NULL`.

- markers.line:

  (optional, logical) In the genepop and structure file, you can write
  the markers on a single line separated by commas
  `markers.line = TRUE`, or have markers on a separate line, i.e. in one
  column, for the genepop file (not very useful with thousands of
  markers) and not printed at all for the structure file. Default:
  `markers.line = TRUE`.

- filename:

  (optional) The file name prefix for the genepop file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `radiator_genepop_`. Default: `filename = NULL`.

- ...:

  other parameters passed to the function.

## Value

A genepop file is saved to the working directory.

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

Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population
genetics software for exact tests and ecumenicism. J. Heredity,
86:248-249

Rousset F. genepop'007: a complete re-implementation of the genepop
software for Windows and Linux. Molecular Ecology Resources. 2008, 8:
103-106. doi:10.1111/j.1471-8286.2007.01931.x

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
