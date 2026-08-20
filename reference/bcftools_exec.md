# Run a bcftools command and log stderr

Internal helper around
[`sys::exec_internal()`](https://jeroen.r-universe.dev/sys/reference/exec.html)
that:

- runs `bcftools.path` with the supplied `args`,

- appends stderr to a log file (if provided),

- optionally aborts on non-zero exit status.

## Usage

``` r
bcftools_exec(
  bcftools.path,
  args,
  log.file = NULL,
  label = NULL,
  verbose = TRUE,
  fail_on_status = TRUE
)
```

## Arguments

- bcftools.path:

  (character) Path or name of the `bcftools` executable.

- args:

  (character) Vector of arguments passed to bcftools.

- log.file:

  (optional, character) Path to a log file. If not `NULL`, stderr is
  appended to this file with a small header. Default: `log.file = NULL`.

- label:

  (optional, character) Short label used in the log header. Default:
  `label = NULL`.

- verbose:

  (logical) If `TRUE`, the constructed command line is printed. Default:
  `verbose = TRUE`.

- fail_on_status:

  (logical) If `TRUE` (default), a non-zero exit status triggers an
  error. If `FALSE`, the status is returned invisibly.

  Default: `fail_on_status = TRUE`.

## Value

Invisibly returns a list with:

- `status` – exit status (integer),

- `stdout` – character scalar with STDOUT,

- `stderr` – character scalar with STDERR.
