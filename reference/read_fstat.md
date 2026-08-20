# Read an FSTAT file into a tidy or wide data frame

Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users. The function `read_fstat` reads a
file in the [fstat](http://www2.unil.ch/popgen/softwares/fstat.htm) file
(Goudet, 1995) into a wide or long/tidy data frame

To manipulate and prune the dataset prior to tidying, use the functions
[`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md)
and
[`genome_translator`](https://thierrygosselin.github.io/genometranslator/reference/genome_translator.md),
that uses blacklist and whitelist along several other filtering options.

## Usage

``` r
read_fstat(data, strata = NULL, tidy = TRUE, filename = NULL, verbose = FALSE)
```

## Arguments

- data:

  A [fstat](http://www2.unil.ch/popgen/softwares/fstat.htm) filename
  with extension `.dat`.

- strata:

  (optional) A tab delimited file with 2 columns. Header: `INDIVIDUALS`
  and `STRATA`. The `STRATA` column can be any hierarchical grouping. To
  create a strata file see
  [`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md).
  Default: `strata = NULL`.

- tidy:

  (optional, logical) With `tidy = FALSE`, the markers are the variables
  and the genotypes the observations (wide format). With the default:
  `tidy = TRUE`, markers and genotypes are variables with their own
  columns (long format). Default: `tidy = TRUE`.

- filename:

  (optional) The file name for the tidy data frame written to the
  working directory. With the default, The tidy data is in the global
  environment only (i.e. not written in the working directory). Default:
  `filename = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

## Value

The output in your global environment is a wide or long/tidy data frame.
If `filename` is provided, the wide or long/tidy data frame is also
written to the working directory.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## References

Goudet J. (1995). FSTAT (Version 1.2): A computer program to calculate
F-statistics. Journal of Heredity 86:485-486

## See also

[hierfstat](https://github.com/jgx65/hierfstat)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# We will use the fstat dataset provided with adegenet package
require("hierfstat")

# The simplest form of the function:
fstat.file <- genometranslator::read_fstat(
    data = system.file(
    "extdata/diploid.dat",
    package = "hierfstat"
    )
 )

# To output a data frame in wide format, with markers in separate columns:
nancycats.wide <- genometranslator::read_fstat(
    data = system.file(
        "extdata/diploid.dat",
        package = "hierfstat"
    ),
tidy = FALSE
)
} # }
```
