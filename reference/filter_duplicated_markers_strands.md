# Filter duplicated SNPs on different strands (+/-)

Filter duplicated SNPs on different strands (+/-)

## Usage

``` r
filter_duplicated_markers_strands(
  gds,
  markers.meta,
  filters.parameters,
  path.folder,
  file.date,
  filter.strands = c("keep.both", "blacklist", "best.stats"),
  detect.strand = FALSE,
  parallel.core = 1L,
  internal = TRUE,
  verbose = TRUE
)
```

## Arguments

- gds:

  A SeqArray GDS object.

- markers.meta:

  Tibble of markers metadata (must include `VARIANT_ID`, `MARKERS`,
  `CHROM`, `LOCUS`, `POS`, `FILTERS`).

- filters.parameters:

  The filters.parameters tracking object.

- path.folder:

  Path to the main filtering folder.

- file.date:

  Character string used in filenames.

- filter.strands:

  (character) Strategy to handle duplicated strands:

  - `"keep.both"` - detect but do not filter;

  - `"blacklist"` - blacklist all duplicated markers;

  - `"best.stats"` - select markers based on missingness and MAC (see
    details in code) and blacklist them (current behaviour).

  Default: `filter.strands = c("keep.both", "blacklist", "best.stats")`.

- detect.strand:

  (logical) Whether strand information was detected in `LOCUS` /
  `STRANDS`. If `FALSE`, the function returns without modification.
  Default: `detect.strand = FALSE`.

- parallel.core:

  (integer) Number of cores to use for
  [`SeqArray::seqAlleleCount()`](https://rdrr.io/pkg/SeqArray/man/seqAlleleFreq.html)
  and
  [`SeqArray::seqMissing()`](https://rdrr.io/pkg/SeqArray/man/seqMissing.html)
  when `filter.strands = "best.stats"`. Default: `parallel.core = 1L`.

- internal:

  (logical) Passed to the internal results-message helper and
  [`genome_parameters()`](https://thierrygosselin.github.io/genometranslator/reference/genome_parameters.md).
  Default: `internal = TRUE`.

- verbose:

  (logical) Verbosity of messages.

  Default: `verbose = TRUE`.

## Value

A list with:

- `markers.meta`: updated tibble;

- `filters.parameters`: updated parameters object.
