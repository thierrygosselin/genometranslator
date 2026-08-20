# Separate markers column into chrom, locus and pos

Genometranslator uses unique marker names by combining `CHROM`, `LOCUS`,
`POS` columns, with double underscore separators, into
`MARKERS = CHROM__LOCUS__POS`.

Used internally in `genometranslator` and might be of interest for users
who need to get back to the original metadata the function provides an
easy way to do it.

## Usage

``` r
separate_markers(
  data,
  sep = "__",
  markers.meta.lists.only = FALSE,
  markers.meta.all.only = FALSE,
  generate.markers.metadata = TRUE,
  generate.ref.alt = FALSE,
  biallelic = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- data:

  An object with a column named `MARKERS`. If `CHROM`, `LOCUS`, `POS`
  are already present, the function returns the dataset untouched. The
  data can be whitelists and blacklists of markers or tidy datasets or
  genome GDS object.

- sep:

  (optional, character) Separator used to identify the different field
  in the `MARKERS` column.

  When the `MARKERS` column doesn't have separator and the function is
  used to generate markers metadata column:
  `"CHROM", "LOCUS", "POS", "REF", "ALT"`, use `sep = NULL`. Default:
  `sep = "__"`.

- markers.meta.lists.only:

  (logical, optional) Allows to keep only the markers metadata:
  `"VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS"`, useful for
  whitelist or blacklist. Default: `markers.meta.lists.only = FALSE`.

- markers.meta.all.only:

  (logica, optionall) Allows to keep all available markers metadata:
  `"VARIANT_ID", "MARKERS", "CHROM", "LOCUS", "POS", "COL", "REF", "ALT"`,
  useful inside genomic translation workflows. Default:
  `markers.meta.all.only = FALSE`.

- generate.markers.metadata:

  (logical, optional) Generate missing markers metadata when missing.
  `"CHROM", "LOCUS", "POS"`.

  Default: `generate.markers.metadata = TRUE`.

- generate.ref.alt:

  (logical, optional) Generate missing REF/ALT alleles with: REF = A and
  ALT = C (for biallelic datasets, only). It is turned off automatically
  when argument `markers.meta.lists.only = TRUE` and on automatically
  when argument `markers.meta.all.only = TRUE` Default:
  `generate.ref.alt = FALSE`.

- biallelic:

  (logical) Speed up the function execution by entering if the dataset
  is biallelic or not. Used internally for verification, before
  generating REF/ALT info. The argument is required when
  `generate.ref.alt = TRUE` and REF/ALT information is absent. Biallelic
  detection is analysis policy and is not performed automatically by
  this package. Default: `biallelic = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

## Value

The same data in the global environment, with 3 new columns: `CHROM`,
`LOCUS`, `POS`. Additionnal columns may be genrated, see arguments
documentation.

## See also

[`generate_markers_metadata`](https://thierrygosselin.github.io/genometranslator/reference/generate_markers_metadata.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
whitelist <- genometranslator::separate_markers(data = whitelist.markers)
tidy.data <- genometranslator::separate_markers(data = bluefintuna.data)
} # }
```
