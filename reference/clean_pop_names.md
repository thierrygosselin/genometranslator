# Clean population's names for genomic workflows

Function to clean pop's name that interfere with some packages or codes.
Space is changed to an underscore `_`. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
clean_pop_names(x, factor = TRUE)
```

## Arguments

- x:

  (character string) Population character string.

- factor:

  (logical) Default: `factor = TRUE`. Will also keep or transform the
  factor levels.
