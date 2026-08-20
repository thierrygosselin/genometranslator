# Read a marker whitelist

Import and standardise a marker whitelist supplied as an object or a
tab-separated file.

## Usage

``` r
read_whitelist(whitelist.markers = NULL, verbose = FALSE)
```

## Arguments

- whitelist.markers:

  A data frame or path to a tab-separated whitelist. It may contain one
  or more of `MARKERS`, `CHROM`, `LOCUS`, `POS`, `VARIANT_ID`, and
  `M_SEQ`. Default: `whitelist.markers = NULL`.

- verbose:

  Logical. Display progress messages. Default: `verbose = FALSE`.

## Value

A distinct whitelist with character fields cleaned for genomic matching,
or `NULL` when no whitelist is supplied.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## See also

[`clean_markers_names`](https://thierrygosselin.github.io/genometranslator/reference/clean_markers_names.md)
