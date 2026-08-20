# Write a [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) file from a tidy data frame

Write a [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) file from a
tidy data frame. The data is bi- or multi-allelic. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_bayescan(
  data,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  A GDS file, open GDS object, tidy data frame object, or a tidy data
  frame in wide or long format in the working directory. For GDS input,
  individual metadata stored in the GDS supplies the codeSTRATA column.
  *How to get a tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- pop.select:

  (optional, string) Selected list of populations for the analysis. e.g.
  `pop.select = c("QUE", "ONT")` to select `QUE` and `ONT` population
  samples (out of 20 pops). If `pop.labels` argument was used to rename
  the strata column, use the new names with `pop.select`. Default:
  `pop.select = NULL`.

- filename:

  (optional) The file name prefix for the bayescan file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `genometranslator_bayescan_`. Default:
  `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

## Value

A bayescan file is written in the working directory.

## Details

BayeScan input should contain polymorphic markers represented in every
selected stratum. Missingness, minor-allele thresholds, linkage
disequilibrium, and other quality-control decisions should be addressed
by the user before calling this writer. The function validates
common-marker representation and polymorphism, but does not remove
failing markers.

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

Foll, M and OE Gaggiotti (2008) A genome scan method to identify
selected loci appropriate for both dominant and codominant markers: A
Bayesian perspective. Genetics 180: 977-993

Foll M, Fischer MC, Heckel G and L Excoffier (2010) Estimating
population structure from AFLP amplification intensity. Molecular
Ecology 19: 4638-4647

Fischer MC, Foll M, Excoffier L and G Heckel (2011) Enhanced AFLP genome
scans detect local adaptation in high-altitude populations of a small
rodent (Microtus arvalis). Molecular Ecology 20: 1450-1462

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
