# extract_genotypes_metadata

Import gds or radiator genotypes meta node

## Usage

``` r
extract_genotypes_metadata(
  gds,
  genotypes.meta.select = NULL,
  genotypes = FALSE,
  metadata.node = TRUE,
  index.only = FALSE,
  sync.markers.individuals = TRUE,
  whitelist = FALSE,
  blacklist = FALSE,
  verbose = FALSE
)
```

## Arguments

- gds:

  The gds object.

- genotypes.meta.select:

  (optional, character) Default:`genotypes.meta.select = NULL`. Default:
  `genotypes.meta.select = NULL`.

- genotypes:

  (optional, character) Default: `genotypes = FALSE`.

- metadata.node:

  (optional, logical) Default: `metadata.node = TRUE`.

- index.only:

  (optional, logical) Default: `index.only = FALSE`.

- sync.markers.individuals:

  (optional, logical) Default: `sync.markers.individuals = TRUE`.

- whitelist:

  (optional, logical) Default: `whitelist = FALSE`.

- blacklist:

  (optional, logical) Default: `blacklist = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.
