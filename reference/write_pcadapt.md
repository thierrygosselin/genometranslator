# Write a [pcadapt](https://github.com/bcm-uga/pcadapt) file from a tidy data frame

Write a [pcadapt](https://github.com/bcm-uga/pcadapt) file from a tidy
data frame. The data is biallelic. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

Prepare an appropriately filtered biallelic dataset before export.
Consider missingness, minor-allele frequency, linkage disequilibrium,
and the sampling design intended for the pcadapt analysis.

## Usage

``` r
write_pcadapt(
  data,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- pop.select:

  (optional, string) Selected list of populations for the analysis. e.g.
  `pop.select = c("QUE", "ONT")` to select `QUE` and `ONT` population
  samples (out of 20 pops). If `pop.labels` argument was used to rename
  the strata column, use the new names with `pop.select`. Default:
  `pop.select = NULL`.

- filename:

  (optional) The file name prefix for the pcadapt file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `radiator_pcadapt_`. Default: `filename = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

## Value

A pcadapt file is written in the working directory a genotype matrix
object is also generated in the global environment.

## Details

**Use a filtered dataset**:

1.  **Control linkage disequilibrium**: Reducing linkage before running
    genome scan is essential. At least start by removing SNPs on the
    same RADseq locus (short linkage disequilibrium).

2.  **Control minor alleles**: Too much Minor Alleles is just noise in
    the data. Filter using Count, Frequency or Depth.

3.  **Only use markers that are in common between strata**

4.  **Only use polymorphic markers**

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

Luu, K., Bazin, E., & Blum, M. G. (2017). pcadapt: an R package to
perform genome scans for selection based on principal component
analysis. Molecular Ecology Resources, 17(1), 67-77.

Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G.
(2015). Detecting genomic signatures of natural selection with principal
component analysis: application to the 1000 Genomes data. Molecular
biology and evolution, msv334.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
