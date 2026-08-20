# Generate strata object from the data

Generate a strata object from the data. The function uses the `POP_ID`
or `STRATA` columns along the `INDIVIDUALS`.

## Usage

``` r
generate_strata(data, pop.id = FALSE)
```

## Arguments

- data:

  A tidy dataset object. Documented in
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- pop.id:

  (logical) When `pop.id = TRUE`, the strata returns the stratification
  colname `POP_ID`. With the default, Returns `STRATA`. Default:
  `pop.id = FALSE`.

## See also

[`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md),
[`summary_strata`](https://thierrygosselin.github.io/genometranslator/reference/summary_strata.md),
[`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md),
[`change_pop_names`](https://thierrygosselin.github.io/genometranslator/reference/change_pop_names.md),
[`join_strata`](https://thierrygosselin.github.io/genometranslator/reference/join_strata.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
