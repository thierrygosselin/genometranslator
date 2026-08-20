# separate_gt

Separate genotype field

## Usage

``` r
separate_gt(
  x,
  gt = "GT_VCF_NUC",
  gather = TRUE,
  haplotypes = FALSE,
  exclude = c("LOCUS", "INDIVIDUALS", "POP_ID"),
  alleles.naming = c("A1", "A2"),
  remove = TRUE,
  filter.missing = FALSE,
  split.chunks = 3,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- gt:

  Default: `gt = "GT_VCF_NUC"`.

- gather:

  Default: `gather = TRUE`.

- haplotypes:

  Default: `haplotypes = FALSE`.

- exclude:

  Default: `exclude = c("LOCUS", "INDIVIDUALS", "POP_ID")`.

- alleles.naming:

  Default: `alleles.naming = c("A1", "A2")`.

- remove:

  Default: `remove = TRUE`.

- filter.missing:

  Default: `filter.missing = FALSE`.

- split.chunks:

  Default: `split.chunks = 3`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.
