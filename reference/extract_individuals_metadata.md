# extract_individuals_metadata

Import individuals metadata from a package or SeqArray node

## Usage

``` r
extract_individuals_metadata(
  gds,
  ind.field.select = NULL,
  metadata.node = TRUE,
  whitelist = FALSE,
  blacklist = FALSE,
  verbose = FALSE
)
```

## Arguments

- gds:

  The gds object.

- ind.field.select:

  (optional, character) Default:`ind.field.select = NULL`. Default:
  `ind.field.select = NULL`.

- metadata.node:

  (optional, logical) Default:`metadata.node = TRUE`. Default:
  `metadata.node = TRUE`.

- whitelist:

  (optional, logical) Default:`whitelist = FALSE`. Default:
  `whitelist = FALSE`.

- blacklist:

  (optional, logical) Default:`blacklist = FALSE`. Default:
  `blacklist = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.
