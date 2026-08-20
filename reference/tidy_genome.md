# Convert a GDS genome to a tidy table

Materialize selected individuals and markers from a GDS file or open GDS
object as a tidy genomic table.

## Usage

``` r
tidy_genome(
  data,
  markers.meta = NULL,
  markers.meta.select = NULL,
  wide = FALSE,
  individuals = NULL,
  pop.id = FALSE,
  calibrate.alleles = TRUE,
  strip.rad = FALSE,
  parallel.core = parallel::detectCores() - 1,
  close.gds = FALSE,
  ...
)
```

## Arguments

- data:

  A GDS filename or open `SeqVarGDSClass` object.

- markers.meta:

  Default: `markers.meta = NULL`.

- markers.meta.select:

  Default: `markers.meta.select = NULL`.

- wide:

  Default: `wide = FALSE`.

- individuals:

  Default: `individuals = NULL`.

- pop.id:

  Logical. When `pop.id = TRUE`, the strata column is named `POP_ID`;
  otherwise it is named `STRATA`. Default: `pop.id = FALSE`.

- calibrate.alleles:

  Default: `calibrate.alleles = TRUE`.

- strip.rad:

  Default: `strip.rad = FALSE`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

- close.gds:

  Logical. Close an open GDS object supplied by the caller. A GDS opened
  internally from a filename is always closed. Default:
  `close.gds = FALSE`.

## Value

A tidy genomic data frame.
