# Write a gsi_sim file from a data frame (wide or long/tidy).

Write a gsi_sim file from a data frame (wide or long/tidy). Used
internally in [assigner](https://github.com/thierrygosselin/assigner)
and
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_gsi_sim(
  data,
  pop.levels = NULL,
  pop.labels = NULL,
  strata = NULL,
  filename = "gsi_sim.unname.txt"
)
```

## Arguments

- data:

  A tidy genomic data set in the working directory tidy formats. *How to
  get a tidy data frame ?* Look for genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- pop.levels:

  (option, string) This refers to the levels in a factor. In this case,
  the id of the pop. Use this argument to have the pop ordered your way
  instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. Default: `pop.levels = NULL`.
  If you find this too complicated, there is also the `strata` argument
  that can do the same thing, see below.

- pop.labels:

  (optional, string) Use this argument to rename/relabel your pop or
  combine your pop. e.g. To combine `"QUE"` and `"ONT"` into a new pop
  called `"NEW"`: (1) First, define the levels for your pop with
  `pop.levels` argument: `pop.levels = c("QUE", "ONT", "ALB")`. (2)
  then, use `pop.labels` argument:
  `pop.levels = c("NEW", "NEW", "ALB")`.#' To rename `"QUE"` to `"TAS"`:
  `pop.labels = c("TAS", "ONT", "ALB")`. Default: `pop.labels = NULL`.
  If you find this too complicated, there is also the `strata` argument
  that can do the same thing, see below.

- strata:

  (optional) A tab delimited file with 2 columns with header:
  `INDIVIDUALS` and `STRATA`. Default: `strata = NULL`. Use this
  argument to rename or change the populations id with the new `STRATA`
  column. The `STRATA` column can be any hierarchical grouping.

- filename:

  The name of the file written to the working directory. Default:
  `filename = "gsi_sim.unname.txt"`.

## Value

A gsi_sim input file is saved to the working directory.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## References

Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008) An
improved method for predicting the accuracy of genetic stock
identification. Canadian Journal of Fisheries and Aquatic Sciences 65,
7:1475-1486.

Anderson, E. C. (2010) Assessing the power of informative subsets of
loci for population assignment: standard methods are upwardly biased.
Molecular ecology resources 10, 4:701-710.

## See also

[gsi_sim](https://github.com/eriqande/gsi_sim) and
[rubias](https://github.com/eriqande/rubias): genetic stock
identification (GSI) in the tidyverse.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
