# Track changes to genomic data

Creates, initiates, or updates a tab-separated history of operations
that change a genomic dataset.

## Usage

``` r
genome_parameters(
  generate = FALSE,
  initiate = FALSE,
  update = TRUE,
  parameter.obj = NULL,
  data = NULL,
  filter.name = "",
  param.name = "",
  values = paste(NULL, NULL, sep = " / "),
  units = "individuals / strata / chrom / locus / markers",
  comments = "",
  path.folder = NULL,
  file.date = NULL,
  internal = FALSE,
  verbose = TRUE
)
```

## Arguments

- generate, initiate, update:

  Logical controls for creating, initiating, and updating the parameter
  history.

- initiate:

  Default: `initiate = FALSE`.

- update:

  Default: `update = TRUE`.

- parameter.obj:

  Existing parameter-history object. Default: `parameter.obj = NULL`.

- data:

  A tidy genomic data frame or supported GDS object. Default:
  `data = NULL`.

- filter.name:

  Name of the operation applied to the data. Default:
  `filter.name = ""`.

- param.name:

  Name of the parameter controlling the operation. Default:
  `param.name = ""`.

- values:

  Parameter value recorded in the history. Default:
  `values = paste(NULL, NULL, sep = " / ")`.

- units:

  Units represented in the before and after summaries. Default:
  `units = "individuals / strata / chrom / locus / markers"`.

- comments:

  Optional comments. Default: `comments = ""`.

- path.folder:

  Output directory. Default: `path.folder = NULL`.

- file.date:

  Date-time label used in the filename. Default: `file.date = NULL`.

- internal, verbose:

  Logical controls for internal use and messages.

- verbose:

  Default: `verbose = TRUE`.

- internal:

  Default: `internal = FALSE`.

- generate:

  Default: `generate = FALSE`.

## Value

A list containing the genomic summary, parameter row, and history file
path as applicable.
