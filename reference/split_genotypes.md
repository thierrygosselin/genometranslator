# Split genotype strings into allele columns

Splits a vector of genotype strings into one row per genotype and one
column per allele. The implementation is vectorized and does not require
parallel processing.

## Usage

``` r
split_genotypes(x, separator = "/", alleles.naming = NULL)
```

## Arguments

- x:

  Character vector of genotype strings.

- separator:

  Fixed separator between alleles. Default: `separator = "/"`.

- alleles.naming:

  Optional column names. When `NULL`, names are generated as `A1`, `A2`,
  and so on.

  Default: `alleles.naming = NULL`.

## Value

A tibble with one row per genotype and one allele per column.
