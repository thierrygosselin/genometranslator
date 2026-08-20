# id_qc_helper

Generate helpers tables and blacklists for individual's QC

## Usage

``` r
id_qc_helper(
  x,
  x.sum,
  stats.id = c("missing", "heterozygosity", "coverage"),
  path.folder = NULL
)
```

## Arguments

- stats.id:

  Default: `stats.id = c("missing", "heterozygosity", "coverage")`.

- path.folder:

  Default: `path.folder = NULL`.
