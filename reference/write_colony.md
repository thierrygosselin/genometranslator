# Write a `COLONY` input file

Write a `COLONY` input file.

## Usage

``` r
write_colony(
  data,
  strata = NULL,
  sample.markers = NULL,
  pop.select = NULL,
  allele.freq = NULL,
  inbreeding = 0,
  mating.sys.males = 0,
  mating.sys.females = 0,
  clone = 0,
  run.length = 2,
  analysis = 1,
  allelic.dropout = 0,
  error.rate = 0.02,
  print.all.colony.opt = FALSE,
  random.seed = NULL,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1,
  filename = NULL,
  ...
)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- strata:

  (path or object) The strata file or object. Additional documentation
  is available in
  [`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md).
  Use that function to whitelist/blacklist populations/individuals.
  Option to set `pop.levels/pop.labels` is also available. Default:
  `strata = NULL`.

- sample.markers:

  (number) `COLONY` can take a long time to run, use a random subsample
  of your markers to speed test `COLONY` e.g. `sample.markers = 500` to
  use only 500 randomly chosen markers. With the default, Will use all
  markers. Default: `sample.markers = NULL`.

- pop.select:

  (optional, string) Selected list of populations for the analysis. e.g.
  `pop.select = c("QUE", "ONT")` to select `QUE` and `ONT` population
  samples (out of 20 pops). If `pop.labels` argument was used to rename
  the strata column, use the new names with `pop.select`. Default:
  `pop.select = NULL`.

- allele.freq:

  (optional, string) Allele frequency can be computed from a select
  group. e.g. `allele.freq = "QUE"` or `allele.freq = c("QUE", "ONT")`.
  Using `allele.freq = "overall"` will use all the samples to compute
  the allele frequency. With the default, Will not compute allele
  frequency. Default: `allele.freq = NULL`.

- inbreeding:

  (boolean) 0/1 no inbreeding/inbreeding. Default: `inbreeding = 0`.

- mating.sys.males:

  (boolean) Mating system in males. 0/1 polygyny/monogyny. Default:
  `mating.sys.males = 0`.

- mating.sys.females:

  (boolean) Mating system in females. 0/1 polygyny/monogyny. Default:
  `mating.sys.females = 0`.

- clone:

  (boolean) Should clones and duplicated individuals be inferred. 0/1,
  yes/no. Default: `clone = 0`.

- run.length:

  (integer) Length of run. 1 (short), 2 (medium), 3 (long), 4 (very
  long). Start with short or medium run and consider longer run if your
  estimates probability are not stable or really good. Default:
  `run.length = 2`.

- analysis:

  (integer) Analysis method. 0 (Pairwise-Likelihood Score), 1 (Full
  Likelihood), 2 (combined Pairwise-Likelihood Score and Full
  Likelihood). Default: `analysis = 1`.

- allelic.dropout:

  Locus allelic dropout rate. Default : `allelic.dropout = 0`. Default:
  `allelic.dropout = 0`.

- error.rate:

  Locus error rate. Default: `error.rate = 0.02`.

- print.all.colony.opt:

  (logical) Should all `COLONY` options be printed in the file.

  **This require manual curation, for the file to work directly with
  `COLONY`**. Default = `print.all.colony.opt = FALSE`. Default:
  `print.all.colony.opt = FALSE`.

- random.seed:

  (integer, optional) For reproducibility, set an integer that will be
  used inside the function that requires randomness. With default, a
  random number is generated and printed in the appropriate output.
  Default: `random.seed = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- filename:

  Name of the acronym for filenaming in the working directory. Default:
  `filename = NULL`.

- ...:

  Additional arguments passed to lower-level readers, translators, or
  writers.

## Value

A `COLONY` file in your working directory (2 if you selected imputations
arguments...)

## Details

**It is highly recommended to read (twice!) the user guide distributed
with `COLONY` to find out the details for input and output of the
software.**

Not all options are provided here.

But to ease the process, all the required options to properly run
`COLONY` will be printed in the file written in your working directory.
Change the values accordingly and wisely.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## References

Jones OR, Wang J (2010) COLONY: a program for parentage and sibship
inference from multilocus genotype data. Molecular Ecology Resources,
10, 551–555.

Wang J (2012) Computationally Efficient Sibship and Parentage Assignment
from Multilocus Marker Data. Genetics, 191, 183–194.

## See also

`COLONY` is available on Jinliang Wang web site
<https://www.zsl.org/science/software/colony>

[colony installation
instructions](https://thierrygosselin.github.io/radiator/articles/rad_genomics_computer_setup.html#colony)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simplest way to run the function with a tidy dataset:
colony.file <- genometranslator::write_colony(data = "turtle.data.rad")
} # }
```
