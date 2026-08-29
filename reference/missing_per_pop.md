# missing_per_pop

Generate missingness per markers per pop helper table using the variants
and samples currently active in the GDS. The incoming GDS filter is
restored before the function returns, including after an error.

## Usage

``` r
missing_per_pop(gds, strata, parallel.core = parallel::detectCores() - 1)
```

## Arguments

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.
