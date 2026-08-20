# dart2gds

Transform dart to GDS format

## Usage

``` r
dart2gds(
  genotypes,
  strata = NULL,
  markers.meta,
  filename.gds,
  dart.check,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- strata:

  Default: `strata = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Default: `verbose = TRUE`.
