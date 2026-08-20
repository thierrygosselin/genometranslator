# Common arguments used by genometranslator

This documentation-only helper centralizes parameters shared by import,
translation, and export functions.

## Usage

``` r
genometranslator_common_arguments(
  interactive.filter = TRUE,
  gds = NULL,
  data = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  random.seed = NULL,
  ...
)
```

## Arguments

- interactive.filter:

  Logical indicating whether an interactive session may ask for
  thresholds. Default: `interactive.filter = TRUE`.

- gds:

  A genome GDS file path or object. Default: `gds = NULL`.

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- random.seed:

  Optional integer seed for operations involving randomness. Default:
  `random.seed = NULL`.

- ...:

  Additional arguments passed to lower-level readers, translators, or
  writers.

## Value

`NULL`, invisibly. This function exists to share documentation.
