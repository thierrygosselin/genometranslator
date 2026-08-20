# Join the strata with the data

The function first filters individuals in data then include the strata.

## Usage

``` r
join_strata(data, strata = NULL, pop.id = FALSE, verbose = TRUE)
```

## Arguments

- data:

  A tidy dataset object. Documented in
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- strata:

  (path or object) The strata file or object. Additional documentation
  is available in
  [`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md).
  Use that function to whitelist/blacklist populations/individuals.
  Option to set `pop.levels/pop.labels` is also available. Default:
  `strata = NULL`.

- pop.id:

  (logical) When `pop.id = TRUE`, the strata returns the stratification
  colname `POP_ID`. With the default, Returns `STRATA`. Default:
  `pop.id = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

## Value

The data filtered by the strata by individuals.

## See also

[`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md),
[`summary_strata`](https://thierrygosselin.github.io/genometranslator/reference/summary_strata.md),
[`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md),
[`change_pop_names`](https://thierrygosselin.github.io/genometranslator/reference/change_pop_names.md),
[`generate_strata`](https://thierrygosselin.github.io/genometranslator/reference/generate_strata.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data <- genometranslator::join_strata(
    data = my_tidy_dataset_object,
    strata = my_strata_object)
} # }
```
