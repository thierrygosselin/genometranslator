# Generate markers metadata

Generate markers metadata: `CHROM, LOCUS, POS, REF, ALT` when missing
from tidy datasets.

## Usage

``` r
generate_markers_metadata(
  data,
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

- generate.markers.metadata:

  (logical, optional) Generate missing markers metadata when missing.
  `"CHROM", "LOCUS", "POS"`. Default:
  `generate.markers.metadata = TRUE`.

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

Depending on argument's value, the same data is returned in the global
environment, with potential these additional columns:
`CHROM, LOCUS, POS, REF, ALT`.

## See also

[`separate_markers`](https://thierrygosselin.github.io/genometranslator/reference/separate_markers.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
tidy.data <- genometranslator::generate_markers_metadata(data = bluefintuna.data)
} # }
```
