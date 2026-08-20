# sync_gds

Synchronize gds with samples and markers. If left NULL, the info is
first search in the radiator node, if not found, it goes in the level
above. An argument also allows to reset the filters.

## Usage

``` r
sync_gds(
  gds,
  samples = NULL,
  variant.id = NULL,
  reset.gds = FALSE,
  reset.filters.m = FALSE,
  reset.filters.i = FALSE,
  verbose = FALSE
)
```

## Arguments

- gds:

  The gds object.

- samples:

  (optional, character string). Will sync the gds object/file with these
  samples. With default, uses the individuals in the radiator node. If
  not found, goes a level above and uses the individuals in the main
  GDS. Default: `samples = NULL`.

- variant.id:

  (optional, integer string). Will sync the gds object/file with these
  variant.id With default, uses the variant.id in the radiator node. If
  not found, goes a level above and uses the variant.id in the main GDS.
  Default: `variant.id = NULL`.

- reset.gds:

  (optional, logical) Default: `reset.gds = FALSE`.

- reset.filters.m:

  (optional, logical) To reset only markers/variant. Default:
  `reset.filters.m = FALSE`.

- reset.filters.i:

  (optional, logical) To reset only individuals. Default:
  `reset.filters.i = FALSE`.

- verbose:

  (optional, logical) Default: `verbose = FALSE`.

## See also

`sync_gds`,
[`list_filters`](https://thierrygosselin.github.io/genometranslator/reference/list_filters.md).
