# Write data in rubias format

Convert a supported genomic dataset to the diploid, two-column-per-locus
data frame used by rubias. The first four columns are `sample_type`,
`repunit`, `collection`, and `indiv`.

## Usage

``` r
write_rubias(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE
)
```

## Arguments

- data:

  A GDS file or object, or another genomic object accepted by
  [`read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.md).

- strata:

  Optional sample metadata accepted by
  [`read_strata()`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md).
  Supply `SAMPLE_TYPE`, `REPUNIT`, `COLLECTION`, and `INDIVIDUALS` for
  full control. Lower-case rubias names are also accepted. When these
  fields are absent, `STRATA` or `POP_ID` is used for both `repunit` and
  `collection`, and all samples are treated as references. Default:
  `strata = NULL`.

- filename:

  Optional filename prefix. When supplied, the result is written as
  `<filename>_rubias.tsv`; a timestamp is added rather than overwriting
  an existing file. Default: `filename = NULL`.

- parallel.core:

  Number of workers available while importing or transforming supported
  inputs. Retained for consistency with other writers. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical. Display progress messages. Default: `verbose = TRUE`.

## Value

A tibble compatible with rubias. When `filename` is supplied, the same
table is also written as a tab-separated file.

## Details

rubias requires one row per diploid individual and four leading
character columns:

- `sample_type` must be `"reference"` or `"mixture"`;

- reference samples require non-missing `repunit` and `collection`;

- mixture samples require `repunit = NA`, while `collection` identifies
  the mixture sample, port, stratum, place, or time group;

- `indiv` must uniquely identify every individual.

Missing genotypes are written as two `NA` alleles. Partial diploid
genotypes are rejected. Locus names cannot contain spaces because rubias
uses adjacent columns to identify gene copies.

This writer creates an input representation; it does not run the
conditional GSI model and does not silently filter individuals or loci.
Quality control, marker selection, reference design, and the distinction
between reference and mixture samples remain the user's responsibility.

## Dependencies

The output schema follows the documented interface of
[rubias](https://github.com/eriqande/rubias), developed by Eric C.
Anderson and Ben Moran. Installing rubias is not required to create the
table, but is required to analyse it. This implementation is
independently written from the public input specification; no rubias
source code is incorporated.

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## References

Moran BM, Anderson EC (2019). Bayesian inference from the conditional
genetic stock identification model. Canadian Journal of Fisheries and
Aquatic Sciences, 76(4), 551-560.
[doi:10.1139/cjfas-2018-0016](https://doi.org/10.1139/cjfas-2018-0016) .

## See also

[rubias](https://github.com/eriqande/rubias),
[`write_gsi_sim()`](https://thierrygosselin.github.io/genometranslator/reference/write_gsi_sim.md)
