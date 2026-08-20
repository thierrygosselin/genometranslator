# Filter markers based on the VCF FILTER column

Filter markers based on the VCF FILTER column

## Usage

``` r
filter_vcf_filter_column(
  gds,
  markers.meta,
  filters.parameters,
  path.folder,
  file.date,
  verbose = TRUE
)
```

## Arguments

- gds:

  A SeqArray GDS object containing the imported VCF.

- markers.meta:

  The current markers metadata tibble (from radiator).

- filters.parameters:

  The filters.parameters object to update.

- path.folder:

  Path to the filtering results folder.

- file.date:

  Character string used for naming result files.

- verbose:

  (logical) Verbosity.

  Default: `verbose = TRUE`.

## Value

A list with:

- `markers.meta`: updated tibble;

- `filters.parameters`: updated parameters object.
