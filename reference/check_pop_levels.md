# Check the use of pop.levels, pop.labels and pop.select arguments.

Check that pop.levels and pop.labels and pop.select arguments are used
correctly and that the values are cleaned for spaces.

## Usage

``` r
check_pop_levels(pop.levels = NULL, pop.labels = NULL, pop.select = NULL)
```

## Arguments

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

- pop.select:

  (optional, string) Selected list of populations for the analysis. e.g.
  `pop.select = c("QUE", "ONT")` to select `QUE` and `ONT` population
  samples (out of 20 pops). If `pop.labels` argument was used to rename
  the strata column, use the new names with `pop.select`. Default:
  `pop.select = NULL`.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
