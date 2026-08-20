# Summary of strata

Summarise the information of a strata file or object. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
summary_strata(strata)
```

## Arguments

- strata:

  (path or object) The strata file or object.

## Value

1.  Number of strata/populations

2.  Number of individuals

3.  Number of individuals per populations

4.  Number of duplicate ids.

## See also

[`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md),
[`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md),
[`change_pop_names`](https://thierrygosselin.github.io/genometranslator/reference/change_pop_names.md),
[`join_strata`](https://thierrygosselin.github.io/genometranslator/reference/join_strata.md),
[`generate_strata`](https://thierrygosselin.github.io/genometranslator/reference/generate_strata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
genometranslator::summary_strata(strata)
} # }
```
