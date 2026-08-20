# Write genomic data

Write genomic data using a format-specific `write_*()` function. The
output format can be supplied explicitly or inferred from an unambiguous
filename extension. Call the specialized writer directly when additional
control is required. With no `output`, the historical Parquet-writing
and GDS-closing behaviour is retained.

When the object is a CoreArray Genomic Data Structure
([GDS](https://github.com/zhengxwen/gdsfmt)) file system, the function
**close the connection with the GDS file**. Before doing so it sets the
filters (variants and samples) based on the info found in the file.

Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users.

## Usage

``` r
write_genome(
  data,
  output = NULL,
  filename = NULL,
  strata = NULL,
  parallel.core = parallel::detectCores() - 1,
  internal = FALSE,
  write.message = "standard",
  verbose = FALSE
)
```

## Arguments

- data:

  An object in the global environment: tidy genomic dataset or GDS
  connection file

- output:

  Optional output format. Supported values include `gds`, `genepop`,
  `fstat`, `genind`, `genlight`, `gtypes`, `plink`, `vcf`, and the other
  formats having a corresponding package `write_*()` function. Default:
  `output = NULL`.

- filename:

  Name of the Arrow/Parquet file written for tabular data. The argument
  is not used when closing a GDS connection.

- strata:

  Optional strata data passed to writers that support it. Default:
  `strata = NULL`.

- parallel.core:

  Number of processor cores passed to writers that support it. Default:
  `parallel.core = parallel::detectCores() - 1`.

- internal:

  (optional, logical) This is used inside radiator internal code and it
  stops from writing the file. Default: `internal = FALSE`.

- write.message:

  (optional, character) Print a message in the console after writing
  file. With `write.message = NULL`, nothing is printed in the console.
  Default: `write.message = "standard"`. This will print
  `message("File written: ", basename(filename))`.

- verbose:

  (optional, logical) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = FALSE`.

## Value

A file written in the working directory or nothing if it's a GDS
connection file.

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

## See also

[appache arrow](https://arrow.apache.org)

[GDS](https://github.com/zhengxwen/gdsfmt)

[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
require(SeqArray)
genometranslator::write_genome(data = tidy.data, filename = "data.shark.arrow.parquet")
genometranslator::write_genome(data = gds.object)
} # }
```
