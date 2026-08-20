# generate_stats

Generate individuals and markers statistics

## Usage

``` r
generate_stats(
  gds,
  individuals = TRUE,
  markers = TRUE,
  missing = TRUE,
  heterozygosity = TRUE,
  coverage = TRUE,
  allele.coverage = TRUE,
  mac = TRUE,
  snp.position.read = TRUE,
  snp.per.locus = TRUE,
  subsample = NULL,
  exhaustive = TRUE,
  force.stats = TRUE,
  path.folder = NULL,
  plot = TRUE,
  digits = 6,
  file.date = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- individuals:

  Default: `individuals = TRUE`.

- markers:

  Default: `markers = TRUE`.

- missing:

  Default: `missing = TRUE`.

- heterozygosity:

  Default: `heterozygosity = TRUE`.

- coverage:

  Default: `coverage = TRUE`.

- allele.coverage:

  Default: `allele.coverage = TRUE`.

- mac:

  Default: `mac = TRUE`.

- snp.position.read:

  Default: `snp.position.read = TRUE`.

- snp.per.locus:

  Default: `snp.per.locus = TRUE`.

- subsample:

  Default: `subsample = NULL`.

- exhaustive:

  Default: `exhaustive = TRUE`.

- force.stats:

  Default: `force.stats = TRUE`.

- path.folder:

  Default: `path.folder = NULL`.

- plot:

  Default: `plot = TRUE`.

- digits:

  Default: `digits = 6`.

- file.date:

  Default: `file.date = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Default: `verbose = TRUE`.
