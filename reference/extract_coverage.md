# check_coverage

Check that the coverage info is in the GDS. By default, it will look for
the DP and AD info in the FORMAT field.

Extract coverage information from a GDS file

## Usage

``` r
check_coverage(
  gds,
  genotypes.metadata.check = FALSE,
  stacks.haplo.check = FALSE,
  dart.check = FALSE
)

extract_coverage(
  gds,
  individuals = TRUE,
  markers = TRUE,
  coverage = TRUE,
  allele.coverage = TRUE,
  coverage.stats = c("sum", "mean", "median", "iqr"),
  subsample.info = 1,
  verbose = TRUE,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- gds:

  The gds object.

- genotypes.metadata.check:

  (optional, logical) Look for already extracted coverage information in
  the radiator genotypes_metadata field of the GDS. Default:
  `genotypes.metadata.check = FALSE`.

- stacks.haplo.check:

  (optional, logical) stacks haplotypes VCF header is baddly generated.
  It will say you have Read and allele Depth info, but you don't.
  Default: `stacks.haplo.check = FALSE`.

- dart.check:

  (optional, logical) DArT have different reporting for coverage
  information. Will no longer report the average coverage stats from 1
  and 2-rows DArT format. Default: `dart.check = FALSE`.

- individuals:

  (optional, logical) Default: `individuals = TRUE`.

- markers:

  (optional, logical) Default: `markers = TRUE`.

- coverage:

  (optional, logical) Default: `coverage = TRUE`.

- allele.coverage:

  (optional, logical) Default: `ad = TRUE`. Default:
  `allele.coverage = TRUE`.

- coverage.stats:

  (optional, character string). Choice of stats to use with coverage.
  Default: `coverage.stats = c("sum", "mean", "median", "iqr")`.

- subsample.info:

  (optional, double) Default: `subsample.info = 1`. The subsample
  proportion used (e.g. 0.3 or none the default).

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.
