# Translate genomic data between formats

A lightweight composition of
[`read_genome`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.md)
and
[`write_genome`](https://thierrygosselin.github.io/genometranslator/reference/write_genome.md).
The input format is detected automatically and read with the
corresponding format-specific reader. The resulting genome is then sent
to the requested writer.

## Usage

``` r
genome_translator(
  data,
  strata = NULL,
  output,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- data:

  A supported genomic file or object.

- strata:

  Optional strata data or filename passed to readers and writers that
  support it. Default: `strata = NULL`.

- output:

  Character vector naming one or more output formats.

- filename:

  Optional output filename or prefix. Default: `filename = NULL`.

- parallel.core:

  Number of processor cores passed to readers and writers that support
  parallel processing. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Value

Invisibly returns the result produced by
[`write_genome()`](https://thierrygosselin.github.io/genometranslator/reference/write_genome.md).
Writers whose purpose is a file side effect may return `NULL`.

## Details

Use a format-specific `read_*()` or `write_*()` function when you need
arguments beyond their defaults.

## Examples

``` r
if (FALSE) { # \dontrun{
genome_translator(
  data = "genomes.vcf",
  strata = "strata.tsv",
  output = "genepop",
  filename = "genomes.gen"
)
} # }
```
