# Read genomic data

Detect a supported genomic format and route it to the appropriate
format-specific reader. The generic interface uses the defaults of each
`read_*()` function; call that reader directly when additional control
is required. Package-native tabular inputs are normalized with the
internal `as_tidy_genome()` helper. The function uses
[`read_parquet`](https://arrow.apache.org/docs/r/reference/read_parquet.html)
or CoreArray Genomic Data Structure
([GDS](https://github.com/zhengxwen/gdsfmt)) file system.

## Usage

``` r
read_genome(
  data,
  strata = NULL,
  columns = NULL,
  allow.dup = FALSE,
  check = TRUE,
  import.metadata = TRUE,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE
)
```

## Arguments

- data:

  A file in the working directory ending with .arrow.parquet or .gds, a
  TSV or legacy RAD/FST file, or an existing wide/tidy genomic table.

- strata:

  Optional strata data or filename passed to readers that support it.
  Default: `strata = NULL`.

- columns:

  (optional) For arrow.parquet file. Column names to read. The default
  is to read all all columns. Default: `columns = NULL`.

- allow.dup:

  (optional, logical) To allow the opening of a GDS file with read-only
  mode when it has been opened in the same R session. Default:
  `allow.dup = FALSE`.

- check:

  (optional, logical) Verify that GDS number of samples and markers
  match. Default: `check = TRUE`.

- import.metadata:

  Logical. Retain columns in addition to the standard genomic columns
  when normalizing tabular input. Default: `import.metadata = TRUE`.

- parallel.core:

  Number of processor cores passed to readers that support parallel
  processing. Default: `parallel.core = parallel::detectCores() - 1`.

- verbose:

  (optional, logical) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = FALSE`.

## Value

A tidy genomic data frame or GDS object (with read/write permissions) in
the global environment.

## Details

For GDS file system, **read_genome** will open the GDS connection file
set the filters (variants and samples) based on the info found in the
file.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

This dispatcher uses the dependencies documented by the selected
format-specific `read_*()` function. Legacy FST/RAD files additionally
require the optional CRAN package fst; VCF preparation may require the
optional `bcftools` executable.

## See also

[arrow](https://github.com/apache/arrow/)
[GDS](https://github.com/zhengxwen/gdsfmt) `read_genome`

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
shark <- genometranslator::read_genome(data = "data.shark.gds")
turtle <- genometranslator::read_genome(data = "data.turtle.arrow.parquet")
} # }
```
