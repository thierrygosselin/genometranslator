# Create a strata file from a list of individuals

If your individuals have a consistent naming scheme (e.g.
SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020), use this
function to rapidly create a strata file. Several functions in
genometranslator and assigner requires a `strata` argument, i.e. a data
frame with the individuals and associated groupings. If you have already
run [stacks](http://catchenlab.life.illinois.edu/stacks/) on your data,
the strata file is similar to a stacks `population map file`, make sure
you have the required column names (`INDIVIDUALS` and `STRATA`).

## Usage

``` r
individuals2strata(data, strata.start, strata.end, filename = NULL)
```

## Arguments

- data:

  A file or data frame object with individuals in a column. The column
  name is `INDIVIDUALS`.

- strata.start:

  (integer) The start of your strata id. See details for more info.

- strata.end:

  (integer) The end of your strata id. See details for more info.

- filename:

  (optional) The file name for the strata object if you want to save it
  in the working directory. With the default, The starta object is in
  the global environment only (i.e. not written in the working
  directory). Default: `filename = NULL`.

## Value

a strata object and file, if requested. The file is tab delimited with 2
columns named: `INDIVIDUALS` and `STRATA`. The `STRATA` column can be
any hierarchical grouping.

## Details

`strata.start` and `strata.end` The info must be found within the name
of your individual sample. If not, you'll have to create a strata file
by hand, the old fashion way. e.g. if your individuals are identified in
this form : SPECIES-POPULATION-MATURITY-YEAR-ID = CHI-QUE-ADU-2014-020,
then, to have the population id in the `STRATA` column,
`strata.start = 5` and `strata.end = 7`. The `STRATA` column can be any
hierarchical grouping.

## See also

[`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md),
[`summary_strata`](https://thierrygosselin.github.io/genometranslator/reference/summary_strata.md),
[`change_pop_names`](https://thierrygosselin.github.io/genometranslator/reference/change_pop_names.md),
[`join_strata`](https://thierrygosselin.github.io/genometranslator/reference/join_strata.md),
[`generate_strata`](https://thierrygosselin.github.io/genometranslator/reference/generate_strata.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
strata.abalone <- individuals2strata(
data = "individuals.abalone.tsv",
strata.start = 5,
strata.end = 7,
filename = "strata.abalone.tsv"
)
} # }
```
