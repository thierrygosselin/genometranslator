# Upgrade a legacy radiator GDS file

Creates a new GDS file and migrates its package metadata namespace from
`/radiator` to `/genometranslator`. The source file is never modified.
Standard SeqArray nodes and all metadata values are preserved.

## Usage

``` r
upgrade_genome_gds(
  gds,
  filename = NULL,
  overwrite = FALSE,
  open = FALSE,
  verbose = TRUE
)
```

## Arguments

- gds:

  Path to a GDS file or an open GDS connection.

- filename:

  Output filename. When `NULL`, `_genometranslator.gds` is appended to
  the source stem. Default: `filename = NULL`.

- overwrite:

  Whether an existing output file may be replaced. Default:
  `overwrite = FALSE`.

- open:

  Whether to return an open, writable SeqArray connection instead of the
  output path. Default: `open = FALSE`.

- verbose:

  Whether to report migration progress.

  Default: `verbose = TRUE`.

## Value

Invisibly returns the output path, or an open SeqArray connection when
`open = TRUE`.
