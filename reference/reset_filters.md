# Reset filters (individuals and markers) in radiator GDS object.

List current active filters in radiator GDS object. Reset specific
filters or all filters at once.

## Usage

``` r
reset_filters(
  gds,
  list.filters = TRUE,
  reset.all = FALSE,
  filter.individuals.missing = FALSE,
  filter.individuals.heterozygosity = FALSE,
  filter.individuals.coverage.total = FALSE,
  detect.mixed.genomes = FALSE,
  detect.duplicate.genomes = FALSE,
  filter.reproducibility = FALSE,
  filter.monomorphic = FALSE,
  filter.common.markers = FALSE,
  filter.ma = FALSE,
  filter.mean.coverage = FALSE,
  filter.genotyping = FALSE,
  filter.snp.position.read = FALSE,
  filter.snp.number = FALSE,
  filter.short.ld = FALSE,
  filter.long.ld = FALSE,
  filter.hwe = FALSE,
  filter.whitelist = FALSE
)
```

## Arguments

- gds:

  A genome GDS file path or object. Default: `gds = NULL`.

- list.filters:

  filters (logical, optional) List current active filters for
  individuals and markers. Default: `list.filters = TRUE`.

- reset.all:

  (logical, optional) Reset all individuals and markers filters.
  Default: `reset.all = FALSE`.

- filter.individuals.missing:

  (logical, optional) Default: `filter.individuals.missing = FALSE`.

- filter.individuals.heterozygosity:

  (logical, optional) Default:
  `filter.individuals.heterozygosity = FALSE`.

- filter.individuals.coverage.total:

  (logical, optional) Default:
  `filter.individuals.coverage.total = FALSE`.

- detect.mixed.genomes:

  (logical, optional) Default: `detect.mixed.genomes = FALSE`.

- detect.duplicate.genomes:

  (logical, optional) Default: `detect.duplicate.genomes = FALSE`.

- filter.reproducibility:

  (logical, optional) Default: `filter.reproducibility = FALSE`.

- filter.monomorphic:

  (logical, optional) Default: `filter.monomorphic = FALSE`.

- filter.common.markers:

  (logical, optional) Default: `filter.common.markers = FALSE`.

- filter.ma:

  (logical, optional) Default: `filter.ma = FALSE`.

- filter.mean.coverage:

  (logical, optional) Default: `filter.mean.coverage = FALSE`.

- filter.genotyping:

  (logical, optional) Default: `filter.genotyping = FALSE`.

- filter.snp.position.read:

  (logical, optional) Default: `filter.snp.position.read = FALSE`.

- filter.snp.number:

  (logical, optional) Default: `filter.snp.number = FALSE`.

- filter.short.ld:

  (logical, optional) Default: `filter.short.ld = FALSE`.

- filter.long.ld:

  (logical, optional) Default: `filter.long.ld = FALSE`.

- filter.hwe:

  (logical, optional) Default: `filter.hwe = FALSE`.

- filter.whitelist:

  (logical, optional) Default: `filter.whitelist = FALSE`.

## See also

[`sync_gds`](https://thierrygosselin.github.io/genometranslator/reference/sync_gds.md),
[`list_filters`](https://thierrygosselin.github.io/genometranslator/reference/list_filters.md).

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# List active filters for individuals and markers
reset_filters(gds)

# You changed your decision concerning the genotyping threshold or
# entered the wrong one:
reset_filters(gds, filter.genotyping = TRUE)
} # }
```
