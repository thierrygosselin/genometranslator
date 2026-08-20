# Transform into a factor the STRATA column, change names and reorder the levels

Transform into a factor the STRATA column, change names and reorder the
levels.

## Usage

``` r
change_pop_names(data, pop.levels = NULL, pop.labels = NULL)
```

## Arguments

- data:

  A GDS filename or open `SeqVarGDSClass` object.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. White spaces in population
  names are replaced by underscore. Default: `pop.levels = NULL`.

- pop.labels:

  (optional, string) Use this argument to rename/relabel your pop or
  combine your pop. e.g. To combine `"QUE"` and `"ONT"` into a new pop
  called `"NEW"`: (1) First, define the levels for your pop with
  `pop.levels` argument: `pop.levels = c("QUE", "ONT", "ALB")`. (2)
  then, use `pop.labels` argument:
  `pop.labels = c("NEW", "NEW", "ALB")`. To rename `"QUE"` to `"TAS"`:
  `pop.labels = c("TAS", "ONT", "ALB")`. Default: `pop.labels = NULL`.
  White spaces in population names are replaced by underscore.

## See also

[`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md),
[`summary_strata`](https://thierrygosselin.github.io/genometranslator/reference/summary_strata.md),
[`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md),
[`join_strata`](https://thierrygosselin.github.io/genometranslator/reference/join_strata.md),
[`generate_strata`](https://thierrygosselin.github.io/genometranslator/reference/generate_strata.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
