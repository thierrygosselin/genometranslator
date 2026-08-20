# Clean individual's names for genomic workflows

Function to clean individual's name that interfere with some packages or
codes. `_` and `:` are changed to a dash `-`. Whitespaces are removed.
Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
clean_ind_names(x)
```

## Arguments

- x:

  (character string) Individuals character string.
