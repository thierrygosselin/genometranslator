# extract_markers_metadata

Import markers metadata from a genometranslator or SeqArray GDS.

## Usage

``` r
extract_markers_metadata(
  gds,
  markers.meta.select = NULL,
  metadata.node = TRUE,
  whitelist = FALSE,
  blacklist = FALSE,
  verbose = FALSE
)
```

## Arguments

- gds:

  The GDS object (genometranslator, legacy radiator, or SeqArray).

- markers.meta.select:

  (optional, character) Names of metadata fields to import. For package
  GDS files, these are the column names in the metadata node. For plain
  SeqArray GDS, standardized names `VARIANT_ID`, `CHROM`, `LOCUS`, `POS`
  are mapped to `variant.id`, `chromosome`, `annotation/id`, `position`.
  Default: `markers.meta.select = NULL`.

- metadata.node:

  (logical) Whether to prefer the package metadata node if present. Both
  `genometranslator` and legacy `radiator` nodes are supported. If the
  node is missing or empty, the function falls back to SeqArray nodes.
  Default: `metadata.node = TRUE`.

- whitelist:

  (logical) If `TRUE` and a `FILTERS` column is present, only rows with
  `FILTERS == "whitelist"` are returned. Default: `whitelist = FALSE`.

- blacklist:

  (logical) If `TRUE` and a `FILTERS` column is present, only rows with
  `FILTERS != "whitelist"` are returned. Default: `blacklist = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.
